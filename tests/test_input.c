/* Synthetic integration tests; no GRMHD/scientific computation or HDF/GSL.
 * From the core repository root:
 *   cc -std=c99 -Wall -Wextra -Wno-unknown-pragmas -fopenmp \
 *      -Isrc -Imodel/bhac tests/test_input.c src/axion.c -lm -o /tmp/test_input
 *   /tmp/test_input
 * Variants: add -DTEST_POL=0, -DTEST_METRIC=1 (BL), -DTEST_METRIC=3 (KS),
 * or -DTEST_METRIC=8 (unsupported CKS). Optional: -fsanitize=address,undefined.
 * The minimal definitions below replace only build-time physics scaffolding;
 * parser, startup and axion implementation are the actual production sources.
 */
#define _POSIX_C_SOURCE 200809L
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define DEFINITIONS_H
#define NDIM 4
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef TEST_POL
#define TEST_POL 1
#endif
#ifndef TEST_METRIC
#define TEST_METRIC 6
#endif
#define POL TEST_POL
#define BL 1
#define KS 3
#define MKSBHAC 6
#define MKSHARM 5
#define CKS 8
#define metric TEST_METRIC
#define num_frequencies 1
#define num_pixels_1d 10
#define tot_pixels 100
#define max_steps 1
#define nspec 4
#define FREQLOG 1
#define FREQFILE 0
#define FREQS FREQLOG
#define SMR 1
#define AMR 0
#define USERSPEC 0
#define UNIF 0
#define MSUN 1.989e33
#define KPCTOCM 3.086e21
struct Camera {
    double IQUV[tot_pixels][num_frequencies][4];
    double tau[tot_pixels][num_frequencies];
    double tauF[tot_pixels][num_frequencies];
    double alpha[tot_pixels], beta[tot_pixels], dx[2];
};

#include "../src/core.c"

int num_blocks, tot_blocks, tot_blocks_alloc, max_level;
double BLOCK_SIZE_X, BLOCK_SIZE_Y;
double R_GRAV, C_CONST, R_ISCO, CUTOFF_INNER, JANSKY_FACTOR;
double a, hslope;
static int stage, init_calls, free_calls, parallel_queries;
static char table_path[256];
static int expect_enabled;

static int counted_init(AxionModel *model, const char *path,
                         double alpha, double omega, double norm, double phase,
                         double spin, double slope, AxionCoordinates coords) {
#ifdef _OPENMP
    assert(!omp_in_parallel());
#endif
    assert(stage == 2 && ++init_calls == 1);
    assert(model == &axion_model);
    assert(spin == 0.5 && slope == 0.3);
    if (norm != 0) {
        assert(alpha == 0.3 && omega == 0.271 && norm == 2 && phase == 3);
#if (metric == BL)
        assert(coords == AXION_BL);
#elif (metric == KS)
        assert(coords == AXION_KS);
#elif (metric == MKSBHAC)
        assert(coords == AXION_MKSBHAC);
#else
        assert(0 && "unsupported enabled metric must not initialize");
#endif
    }
    int status = axion_model_init(model, path, alpha, omega, norm, phase,
                                  spin, slope, coords);
    if (status == AXION_OK) stage = 3;
    return status;
}
static void counted_free(AxionModel *model) {
#ifdef _OPENMP
    assert(!omp_in_parallel());
#endif
    assert(stage == 6 && ++free_calls == 1);
    axion_model_free(model);
}
#define main raptor_main
#define axion_model_init counted_init
#define axion_model_free counted_free
#include "../src/main.c"
#undef main
#undef axion_model_init
#undef axion_model_free

void init_model(void) {
    assert(stage == 0);
    a = 0.5;
    hslope = 0.3;
    stage = 1;
}
void set_constants(void) {
    assert(stage == 1);
    stage = 2;
}
void init_camera(struct Camera **camera) {
    assert(stage == 3 && init_calls == 1);
    assert(MBH == 4.14e6 * MSUN && SOURCE_DIST == 8.127 * KPCTOCM);
    assert(strcmp(GRMHD_FILE, "synthetic-grmhd") == 0 && TIME_INIT == 7.5);
    if (expect_enabled) {
        size_t na, nr;
        double radial;
        assert(axion_model.table != NULL);
        assert(axion_table_dimensions(axion_model.table, &na, &nr) == AXION_OK);
        assert(na == 2 && nr == 3);
        assert(axion_table_interpolate(axion_model.table, 0.3, 3, &radial) == AXION_OK);
        assert(fabs(radial - 3) < 1e-12);
        /* Queries still succeed after deletion: no lazy/per-ray table load. */
        assert(unlink(table_path) == 0);
    } else {
        assert(axion_model.table == NULL && axion_model.norm == 0);
    }
    *camera = calloc(1, sizeof(**camera));
    assert(*camera != NULL);
    (*camera)->dx[0] = (*camera)->dx[1] = 1;
    tot_blocks = 1;
    stage = 4;
}
void integrate_geodesic(double alpha, double beta, double *path,
                        int *steps, double cutoff) {
    (void)alpha; (void)beta; (void)path; (void)cutoff;
    *steps = 0;
}
void radiative_transfer_polarized(double *path, int steps, double frequency,
                                  double *fx, double *fy, double *p, int print,
                                  double *iquv, double *tau, double *tauf,
                                  const AxionModel *axion, int block, int pixel) {
    (void)path; (void)steps; (void)frequency; (void)fx; (void)fy; (void)p;
    (void)print; (void)iquv; (void)tau; (void)tauf; (void)block; (void)pixel;
    assert(stage >= 4 && init_calls == 1 && axion == &axion_model);
    double X[4] = {0, 3, 0.7, 0.2}, value = NAN;
#if (metric == MKSBHAC)
    X[1] = log(3);
#endif
    assert(axion_field_value(axion, X, &value) == AXION_OK);
    assert(isfinite(value));
    if (!expect_enabled) assert(value == 0);
#ifdef _OPENMP
    assert(omp_in_parallel());
#pragma omp atomic
#endif
    parallel_queries++;
}
double radiative_transfer_unpolarized(double *path, int steps, double *freq,
                                      double iquv[num_frequencies][4],
                                      double tau[num_frequencies]) {
    (void)path; (void)steps; (void)freq; (void)iquv; (void)tau;
    assert(!expect_enabled && init_calls == 1);
    return 0;
}
void prerun_refine(struct Camera **camera) {
    assert(stage == 4 && init_calls == 1);
    double frequencies[num_frequencies] = {230e9};
    calculate_image_block(*camera, frequencies, 0);
    stage = 5;
}
void output_files(struct Camera *camera, double spec[num_frequencies][nspec],
                   double frequencies[num_frequencies]) {
    (void)camera; (void)spec; (void)frequencies;
    assert(stage == 5 && init_calls == 1);
#if POL
    assert(parallel_queries == 2 * tot_pixels);
#endif
    stage = 6;
}

static const char *triples[] = {
    "MBH (Msun) 4.14e6", "DISTANCE (kpc) 8.127", "M_UNIT (g) 1e18",
    "R_LOW (-) 1", "R_HIGH (-) 1", "INCLINATION (deg) 50",
    "IMG_WIDTH (pixels) 200", "IMG_HEIGHT (pixels) 200",
    "CAM_SIZE_X (Rg) 40", "CAM_SIZE_Y (Rg) 40", "FREQS_PER_DEC (-) 1",
    "FREQ_MIN (Hz) 230e9", "STEPSIZE (-) 0.01", "AXION_ALPHA (-) 0.3",
    "AXION_OMEGA (-) 0.271", "AXION_NORM (-) 0", "AXION_PHASE (-) 3",
    NULL, "MAX_LEVEL (-) 1"
};
static char input_path[256];
static unsigned tests;

/* replace=-1 for none; omit=-1 for none; reverse tests arbitrary ordering. */
static void make_input(int omit, int replace, const char *replacement,
                        int compact, int reverse, int enabled,
                        const char *extra) {
    FILE *f = fopen(input_path, "w");
    assert(f != NULL);
    for (int k = 0; k < 19; ++k) {
        int i = reverse ? 18 - k : k;
        if (i == omit) continue;
        if (i == replace) fprintf(f, "%s", replacement);
        else if (i == 17) fprintf(f, "AXION_TABLE (-) %s", table_path);
        else if (i == 15 && enabled) fprintf(f, "AXION_NORM (-) 2");
        else fprintf(f, "%s", triples[i]);
        fputs(compact ? " " : "\n", f);
    }
    if (extra != NULL) fputs(extra, f);
    assert(fclose(f) == 0);
}
#if POL && (metric == BL || metric == KS || metric == MKSBHAC)
static void make_table(void) {
    FILE *f = fopen(table_path, "w");
    assert(f != NULL);
    fputs("0.2,2,1\n0.2,5,4\n0.2,30,2\n"
          "0.4,2,3\n0.4,5,6\n0.4,30,4\n", f);
    assert(fclose(f) == 0);
}
#endif

/* Fork isolates exit-on-error parser and global startup state. Capture stderr
 * to verify useful field/migration diagnostics, not merely nonzero status. */
static void run_case(const char *name, int success, const char *diagnostic,
                      int startup, int enabled, int argc_override,
                      const char *time_arg, const char *grmhd_arg) {
    int fd[2];
    assert(pipe(fd) == 0);
    fflush(NULL);
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0) {
        close(fd[0]);
        assert(dup2(fd[1], STDERR_FILENO) >= 0);
        close(fd[1]);
        FILE *null = freopen("/dev/null", "w", stdout);
        assert(null != NULL);
        char *argv[] = {"RAPTOR", input_path,
                        (char *)(grmhd_arg ? grmhd_arg : "synthetic-grmhd"),
                        (char *)(time_arg ? time_arg : "7.5"), NULL};
        expect_enabled = enabled;
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(2);
#endif
        if (startup) {
            int result = raptor_main(argc_override < 0 ? 4 : argc_override, argv);
            if (result == 0) {
                assert(stage == 6 && init_calls == 1 && free_calls == 1);
                assert(axion_model.table == NULL);
            }
            exit(result);
        }
        read_model(argv);
        assert(stage == 0 && init_calls == 0 && axion_model.table == NULL);
        assert(axion_alpha == 0.3 || axion_alpha == 0);
        assert(axion_omega == 0.271 || axion_omega == 0);
        assert(phase_index >= 0 && phase_index <= 15);
        assert(strcmp(axion_table, table_path) == 0);
        exit(0);
    }
    close(fd[1]);
    char errors[16384];
    size_t n = 0;
    ssize_t got;
    while ((got = read(fd[0], errors + n, sizeof(errors) - 1 - n)) > 0)
        n += (size_t)got;
    errors[n] = '\0';
    close(fd[0]);
    int status;
    assert(waitpid(pid, &status, 0) == pid);
    if (!WIFEXITED(status) || (WEXITSTATUS(status) == 0) != success ||
        (diagnostic && !strstr(errors, diagnostic))) {
        fprintf(stderr, "FAIL %s: status=%d expected success=%d diagnostic=%s\n%s\n",
                name, status, success, diagnostic ? diagnostic : "(none)", errors);
        exit(1);
    }
    ++tests;
}
#define PARSE(name, ok, msg) run_case(name, ok, msg, 0, 0, -1, NULL, NULL)

int main(void) {
    char dir[] = "/tmp/raptor-input-XXXXXX";
    assert(mkdtemp(dir) != NULL);
    assert(snprintf(input_path, sizeof(input_path), "%s/model.in", dir) > 0);
    assert(snprintf(table_path, sizeof(table_path), "%s/table.csv", dir) > 0);
    make_input(-1, -1, NULL, 0, 0, 0, NULL);
    PARSE("full baseline without table", 1, NULL);
    make_input(-1, -1, NULL, 1, 0, 0, NULL);
    PARSE("no newlines", 1, NULL);
    make_input(-1, -1, NULL, 0, 1, 0, NULL);
    PARSE("reverse named order", 1, NULL);

    for (int i = 0; i < 19; ++i) {
        char label[64];
        if (i == 17) strcpy(label, "AXION_TABLE");
        else {
            size_t len = strcspn(triples[i], " ");
            memcpy(label, triples[i], len);
            label[len] = '\0';
        }
        make_input(i, -1, NULL, 0, 0, 0, NULL);
        PARSE("each field required", 0, label);
    }
    /* Every strict prefix of a complete file is rejected, including truncation
     * between label/unit/value and after complete preceding triples. */
    make_input(-1, -1, NULL, 1, 0, 0, NULL);
    FILE *full = fopen(input_path, "r");
    assert(full != NULL);
    char contents[4096];
    size_t bytes = fread(contents, 1, sizeof(contents) - 1, full);
    assert(!ferror(full) && fclose(full) == 0 && bytes < sizeof(contents) - 1);
    contents[bytes] = '\0';
    for (size_t cut = 0; cut < bytes; ++cut) {
        if (contents[cut] != ' ') continue;
        /* The final separator is a valid EOF, not a truncated input. */
        if (cut == bytes - 1) continue;
        FILE *shortfile = fopen(input_path, "w");
        assert(shortfile != NULL);
        assert(fwrite(contents, 1, cut, shortfile) == cut);
        assert(fclose(shortfile) == 0);
        PARSE("each truncated token-boundary prefix", 0, "Input error");
    }
    make_input(13, -1, NULL, 1, 0, 0, NULL);
    PARSE("legacy omega-only compact file", 0, "Legacy omega-only");
    FILE *f = fopen(input_path, "w");
    assert(f != NULL);
    fputs("MBH (Msun) 4.14e6 DISTANCE (kpc) 8.127 M_UNIT (g) 1e18 "
          "R_LOW (-) 1 R_HIGH (-) 1 INCLINATION (deg) 50 "
          "IMG_WIDTH (pixels) 200 IMG_HEIGHT (pixels) 200 "
          "CAM_SIZE_X (Rg) 40 CAM_SIZE_Y (Rg) 40 FREQS_PER_DEC (-) 1 "
          "FREQ_MIN (Hz) 230e9 STEPSIZE (-) 0.01 AXION_NORM (-) 1 "
          "AXION_OMEGA (-) 0.2 AXION_PHASE (-) 0 MAX_LEVEL (-) 1", f);
    assert(fclose(f) == 0);
    PARSE("complete original 17-field old example", 0, "Run.py needs redesign");
    f = fopen(input_path, "w");
    assert(f != NULL && fclose(f) == 0);
    PARSE("empty file", 0, "AXION_ALPHA");
    make_input(-1, 18, "MAX_LEVEL (-)", 0, 0, 0, NULL);
    PARSE("short final triple", 0, "incomplete");
    make_input(-1, -1, NULL, 0, 0, 0, "AXION_ALPHA (-) 0.3");
    PARSE("duplicate field", 0, "duplicate");
    make_input(-1, -1, NULL, 0, 0, 0, "GARBAGE");
    PARSE("unknown extra field", 0, "unknown field");
    struct BadInput { int field; const char *text, *diagnostic; } bad[] = {
        {13, "axion_alpha (-) 0.3", "unknown field"},
        {13, "AXION_ALPHA (Hz) 0.3", "incorrect unit"},
        {14, "OMEGA (-) 0.271", "unknown field"},
        {0, "MBH (Msun) nan", "finite number"},
        {0, "MBH (Msun) inf", "finite number"},
        {0, "MBH (Msun) 1e309", "finite number"},
        {0, "MBH (Msun) 1e300", "cgs conversion"},
        {0, "MBH (Msun) 4junk", "finite number"},
        {0, "MBH (Msun) 0", "must be positive"},
        {2, "M_UNIT (g) -1", "allowed range"},
        {5, "INCLINATION (deg) 181", "allowed range"},
        {6, "IMG_WIDTH (pixels) 201", "image grid"},
        {6, "IMG_WIDTH (pixels) 100", "image grid"},
        {6, "IMG_WIDTH (pixels) 0", "allowed range"},
        {6, "IMG_WIDTH (pixels) 2147483648", "in-range integer"},
        {6, "IMG_WIDTH (pixels) 200.0", "in-range integer"},
        {10, "FREQS_PER_DEC (-) 0", "must be positive"},
        {11, "FREQ_MIN (Hz) -1", "allowed range"},
        {12, "STEPSIZE (-) nan", "finite number"},
        {13, "AXION_ALPHA (-) -0.1", "allowed range"},
        {14, "AXION_OMEGA (-) nan", "finite number"},
        {15, "AXION_NORM (-) -1", "allowed range"},
        {16, "AXION_PHASE (-) 16", "allowed range"},
        {16, "AXION_PHASE (-) -1", "allowed range"},
        {16, "AXION_PHASE (-) 1.5", "in-range integer"},
        {18, "MAX_LEVEL (-) 0", "allowed range"},
        {17, "AXION_TABLE (-) filename withspaces", "unknown field"}
    };
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        make_input(-1, bad[i].field, bad[i].text, 0, 0, 0, NULL);
        PARSE(bad[i].text, 0, bad[i].diagnostic);
    }
    char long_token[400], long_triple[450];
    memset(long_token, 'x', sizeof(long_token) - 1);
    long_token[sizeof(long_token) - 1] = '\0';
    snprintf(long_triple, sizeof(long_triple), "AXION_TABLE (-) %s", long_token);
    make_input(-1, 17, long_triple, 0, 0, 0, NULL);
    PARSE("bounded table filename", 0, "255 bytes");
    make_input(-1, -1, NULL, 0, 0, 0, NULL);
    FILE *binary = fopen(input_path, "a");
    assert(binary != NULL);
    assert(fwrite("\0ignored", 1, 8, binary) == 8);
    assert(fclose(binary) == 0);
    PARSE("embedded NUL rejected", 0, "read error");
    make_input(-1, -1, NULL, 0, 0, 0, NULL);
    run_case("bounded GRMHD filename", 0, "255 bytes", 0, 0, -1, NULL, long_token);
    run_case("TIME_INIT finite", 0, "TIME_INIT", 0, 0, -1, "nan", NULL);
    run_case("TIME_INIT trailing text", 0, "TIME_INIT", 0, 0, -1, "7.5junk", NULL);
    run_case("argc too short", 0, "Usage:", 1, 0, 1, NULL, NULL);
    run_case("argc too long", 0, "Usage:", 1, 0, 5, NULL, NULL);
    make_input(-1, -1, NULL, 0, 0, 0,
               NULL);
    run_case("baseline single initialization no missing-table load", 1, NULL,
             1, 0, -1, NULL, NULL);
    make_input(-1, 13, "AXION_ALPHA (-) 0", 0, 0, 0, NULL);
    FILE *zero = fopen(input_path, "r+");
    assert(zero != NULL);
    bytes = fread(contents, 1, sizeof(contents) - 1, zero);
    contents[bytes] = '\0';
    char *omega = strstr(contents, "AXION_OMEGA (-) 0.271");
    assert(omega != NULL);
    /* Equal-width replacement to retain the rest of the input. */
    memcpy(omega, "AXION_OMEGA (-) 0    ", strlen("AXION_OMEGA (-) 0.271"));
    rewind(zero);
    assert(fwrite(contents, 1, bytes, zero) == bytes && fclose(zero) == 0);
    run_case("all-zero alpha/omega baseline startup", 1, NULL,
             1, 0, -1, NULL, NULL);
    make_input(-1, 13, "AXION_ALPHA (-) 0", 0, 0, 0, NULL);
    PARSE("zero baseline alpha", 1, NULL);
    make_input(-1, 14, "AXION_OMEGA (-) 0", 0, 0, 0, NULL);
    PARSE("zero baseline omega", 1, NULL);
    make_input(-1, 16, "AXION_PHASE (-) 15", 0, 0, 0, NULL);
    PARSE("phase upper endpoint", 1, NULL);
    make_input(-1, 16, "AXION_PHASE (-) 0", 0, 0, 0, NULL);
    PARSE("phase lower endpoint", 1, NULL);
    make_input(-1, 13, "AXION_ALPHA (-) 0", 0, 0, 1, NULL);
    PARSE("enabled alpha must be positive", 0, "both must be positive");
    make_input(-1, 14, "AXION_OMEGA (-) 0", 0, 0, 1, NULL);
    PARSE("enabled omega must be positive", 0, "both must be positive");
    make_input(-1, -1, NULL, 0, 0, 1, NULL);
#if !POL
    run_case("POL0 enabled rejected", 0, "requires POL=1", 1, 1, -1, NULL, NULL);
#elif metric != BL && metric != KS && metric != MKSBHAC
    run_case("unsupported enabled coordinates rejected", 0, "unsupported", 1, 1, -1, NULL, NULL);
#else
    run_case("missing enabled table startup fails", 0, "Axion initialization failed", 1, 1, -1, NULL, NULL);
    make_table();
    run_case("enabled single initialization before parallel queries", 1, NULL, 1, 1, -1, NULL, NULL);
#endif
    unlink(table_path);
    assert(unlink(input_path) == 0 && rmdir(dir) == 0);
    printf("PASS %u synthetic input/startup tests (POL=%d metric=%d)\n", tests, POL, metric);
    return 0;
}
