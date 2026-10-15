/* Standalone tests, no HDF or full RAPTOR build. Run only newly built binary:
 * cc -std=c99 -O2 -Wall -Wextra -Wpedantic -Werror -fopenmp \
 *    -Isrc src/axion.c tests/test_axion.c -lm -o /tmp/test-axion
 * scratch=$(mktemp -d /tmp/axion-tests.XXXXXX)
 * /tmp/test-axion "$scratch" [optional existing table paths ...]
 * rm -r "$scratch"
 * The scratch directory must exist. Optional tables are opened read-only.
 */
#include "axion.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.141592653589793238462643383279502884
#define CHECK(expr) do { if (!(expr)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr); exit(1); \
} } while (0)
#define OK(expr) CHECK((expr) == AXION_OK)
#define FAIL(expr) CHECK((expr) != AXION_OK)

static char path[4096];
static const char rectangular[] =
    "0.1,3,1\n0.1,5,3\n0.1,25,5\n"
    "0.3,3,5\n0.3,5,7\n0.3,25,9\n";

static void close_to(double actual, double expected)
{
    CHECK(isfinite(actual));
    CHECK(fabs(actual - expected) <= 3e-12 * fmax(1.0, fabs(expected)));
}

static void fixture_bytes(const char *text, size_t size)
{
    FILE *fp = fopen(path, "wb");
    CHECK(fp != NULL);
    CHECK(fwrite(text, 1, size, fp) == size);
    CHECK(fclose(fp) == 0);
}

static void fixture(const char *text)
{
    fixture_bytes(text, strlen(text));
}

static void test_tables(void)
{
    AxionTable *table = NULL;
    size_t na = 0, nr = 0;
    double v = -123;
    int i;
    const double queries[][3] = {
        {0.1,3,1}, {0.1,25,5}, {0.3,3,5}, {0.3,25,9},
        {0.2,4,4}, {0.1,5,3}, {0.3,5,7}, {0.2,15,6}
    };
    fixture(rectangular);
    OK(axion_table_load(path, &table));
    CHECK(table != NULL);
    OK(axion_table_dimensions(table, &na, &nr));
    CHECK(na == 2 && nr == 3);
    for (i = 0; i < (int)(sizeof(queries)/sizeof(queries[0])); ++i) {
        OK(axion_table_interpolate(table, queries[i][0], queries[i][1], &v));
        close_to(v, queries[i][2]);
    }
    v = -123;
    CHECK(axion_table_interpolate(table, nextafter(0.1, -INFINITY), 3, &v) == AXION_ERR_DOMAIN);
    CHECK(v == -123);
    FAIL(axion_table_interpolate(table, nextafter(0.3, INFINITY), 3, &v));
    FAIL(axion_table_interpolate(table, 0.1, nextafter(3, -INFINITY), &v));
    FAIL(axion_table_interpolate(table, 0.1, nextafter(25, INFINITY), &v));
    FAIL(axion_table_interpolate(table, NAN, 3, &v));
    FAIL(axion_table_interpolate(table, 0.1, INFINITY, &v));
    FAIL(axion_table_interpolate(table, 0.1, NAN, &v));
    FAIL(axion_table_interpolate(table, 0.1, 3, NULL));
    FAIL(axion_table_dimensions(table, NULL, &nr));
    axion_table_free(table);

    fixture(" \t0.2 , 3 , 2 \r\n0.2, 5, 6"); /* no final newline */
    OK(axion_table_load(path, &table));
    OK(axion_table_dimensions(table, &na, &nr));
    CHECK(na == 1 && nr == 2);
    OK(axion_table_interpolate(table, 0.2, 4, &v)); close_to(v, 4);
    OK(axion_table_interpolate(table, 0.2, 3, &v)); close_to(v, 2);
    OK(axion_table_interpolate(table, 0.2, 5, &v)); close_to(v, 6);
    FAIL(axion_table_interpolate(table, nextafter(0.2, 0), 4, &v));
    FAIL(axion_table_interpolate(table, nextafter(0.2, 1), 4, &v));
    axion_table_free(table);
    fixture(".2,3,1e-320\n.2,5,2e-320\n");
    OK(axion_table_load(path, &table)); /* finite representable subnormal R */
    OK(axion_table_interpolate(table, .2, 3, &v));
    CHECK(v > 0 && v < DBL_MIN);
    axion_table_free(table);
    axion_table_free(NULL);
    FAIL(axion_table_interpolate(NULL, 0.1, 3, &v));
    FAIL(axion_table_dimensions(NULL, &na, &nr));
    FAIL(axion_table_load(path, NULL));
    FAIL(axion_table_load(NULL, &table)); CHECK(table == NULL);
    FAIL(axion_table_load("", &table)); CHECK(table == NULL);
    FAIL(axion_table_load("/nonexistent/axion/test.csv", &table)); CHECK(table == NULL);
    puts("PASS table interpolation, inclusive endpoints, single alpha, errors");
}

static void test_bad_tables(void)
{
    const char *bad[] = {
        "", "\n", "alpha,r,R\n", "0.1,3,1\n", /* empty/header/short */
        "0.1,3\n0.1,4,2\n", "0.1,3,1,9\n0.1,4,2\n",
        "0.1,3,1garbage\n0.1,4,2\n", "0.1;3;1\n0.1,4,2\n",
        "0.1,,1\n0.1,4,2\n", "0.1,3,1\n\n0.1,4,2\n",
        "nan,3,1\nnan,4,2\n", "0.1,inf,1\n0.1,4,2\n",
        "0.1,3,nan\n0.1,4,2\n", "0.1,3,inf\n0.1,4,2\n",
        "0.1,3,-1\n0.1,4,2\n", "0.1,3,1e9999\n0.1,4,2\n",
        "0.1,3,1e-9999\n0.1,4,2\n", /* no silent textual underflow */
        "0.1,3,1\n0.1,3,2\n", "0.1,4,1\n0.1,3,2\n",
        "0.3,3,1\n0.3,4,2\n0.1,3,3\n0.1,4,4\n", /* alpha decreases */
        "0.1,3,1\n0.1,4,2\n0.3,3,3\n", /* incomplete last slice */
        "0.1,3,1\n0.1,4,2\n0.3,3,3\n0.3,5,4\n", /* different r */
        "0.1,3,1\n0.1,4,2\n0.3,4,3\n0.3,3,4\n", /* unsorted slice */
        "0.1,3,1\n0.3,3,2\n", /* nr=1 */
        "0.1,3,1\n0.1,4,2\n0.3,3,3\n0.3,4,4\n0.1,3,5\n0.1,4,6\n",
        "0.1,3,1\n0.1,4,2\n0.3,3,3\n0.4,4,4\n" /* nonconstant alpha slice */
    };
    size_t i;
    for (i = 0; i < sizeof(bad)/sizeof(bad[0]); ++i) {
        AxionTable *table = NULL;
        fixture(bad[i]);
        FAIL(axion_table_load(path, &table));
        CHECK(table == NULL);
    }
    {
        const char binary[] = "0.1,3,1\0junk\n0.1,4,2\n";
        AxionTable *table = NULL;
        fixture_bytes(binary, sizeof(binary) - 1);
        FAIL(axion_table_load(path, &table)); CHECK(table == NULL);
    }
    puts("PASS malformed/short/nonrectangular/nonmonotonic/nonfinite tables");
}

static void test_dynamic_table(void)
{
    AxionTable *table = NULL;
    size_t na, nr;
    int i, j;
    double v;
    FILE *fp = fopen(path, "w");
    CHECK(fp != NULL);
    /* Exceeds BOTH old hardcoded dimensions. Includes a very long CSV line. */
    for (i = 0; i < 37; ++i)
        for (j = 0; j < 3101; ++j) {
            if (i == 0 && j == 0) {
                int k;
                for (k = 0; k < 12000; ++k) CHECK(fputc(' ', fp) != EOF);
            }
            CHECK(fprintf(fp, "%.17g,%d,%d\n", 0.1 + 0.01 * i, 3+j, 2*i+j) > 0);
        }
    CHECK(fclose(fp) == 0);
    OK(axion_table_load(path, &table));
    OK(axion_table_dimensions(table, &na, &nr)); CHECK(na == 37 && nr == 3101);
    OK(axion_table_interpolate(table, 0.1, 3, &v)); close_to(v, 0);
    OK(axion_table_interpolate(table, 0.1 + 0.01 * 36, 3103, &v)); close_to(v, 3172);
    OK(axion_table_interpolate(table, 0.225, 125.5, &v)); close_to(v, 147.5);
    axion_table_free(table);
    /* Extreme finite axis differences must not cause division by infinity. */
    fp = fopen(path, "w"); CHECK(fp != NULL);
    CHECK(fprintf(fp, "%.17g,%.17g,0\n%.17g,%.17g,2\n"
                      "%.17g,%.17g,2\n%.17g,%.17g,4\n",
                  -DBL_MAX, -DBL_MAX, -DBL_MAX, DBL_MAX,
                  DBL_MAX, -DBL_MAX, DBL_MAX, DBL_MAX) > 0);
    CHECK(fclose(fp) == 0);
    OK(axion_table_load(path, &table));
    OK(axion_table_interpolate(table, 0, 0, &v)); close_to(v, 2);
    axion_table_free(table);
    puts("PASS dynamic dimensions, long CSV lines, extreme finite axes");
}

static void init_model(AxionModel *m, double alpha, double omega, double norm,
                       double phase, double spin, double hslope, AxionCoordinates coords)
{
    OK(axion_model_init(m, path, alpha, omega, norm, phase, spin, hslope, coords));
}

static double field(const AxionModel *m, const double X[4])
{
    double result = -999;
    OK(axion_field_value(m, X, &result));
    return result;
}

static void test_baseline_and_parameters(void)
{
    AxionModel m;
    double v = 3;
    double X[4] = {NAN, NAN, NAN, NAN};
    int coords;
    for (coords = AXION_BL; coords <= AXION_CKS + 1; ++coords) {
        OK(axion_model_init(&m, "/missing/never-read.csv", NAN, NAN, 0, NAN,
                            NAN, NAN, (AxionCoordinates)coords));
        CHECK(m.norm == 0 && m.table == NULL);
        OK(axion_field_value(&m, NULL, &v)); CHECK(v == 0);
        OK(axion_field_value(&m, X, &v)); CHECK(v == 0);
        /* Explicitly prove zero fastpath precedes all parameter checks. */
        m.spin = NAN; m.omega = NAN; m.alpha = NAN; m.phase = NAN;
        m.hslope = NAN; m.coords = AXION_CKS;
        OK(axion_field_value(&m, NULL, &v)); CHECK(v == 0);
        axion_model_free(&m); axion_model_free(&m);
    }
    fixture(rectangular);
#define BAD_INIT(a,o,n,p,s,h,c) do { \
    FAIL(axion_model_init(&m, path, a,o,n,p,s,h,c)); \
    CHECK(m.table == NULL && m.norm == 0); axion_model_free(&m); \
} while (0)
    BAD_INIT(0, .2, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(-.1, .2, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(NAN, .2, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, 0, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, -.1, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, INFINITY, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, -1, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, NAN, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, INFINITY, 0, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, -1, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, 16, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, .5, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, NAN, .5, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, 0, 1, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, 0, -1, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, 0, NAN, 0, AXION_BL);
    BAD_INIT(.1, .2, 1, 0, .5, NAN, AXION_BL);
    BAD_INIT(.1, .2, 1, 0, .5, 0, AXION_MKS);
    BAD_INIT(.1, .2, 1, 0, .5, 0, AXION_MKSHARM);
    BAD_INIT(.1, .2, 1, 0, .5, 0, AXION_CKS);
    BAD_INIT(.1, .2, 1, 0, .5, 0, (AxionCoordinates)99);
    BAD_INIT(.09, .2, 1, 0, .5, 0, AXION_BL);
    BAD_INIT(.31, .2, 1, 0, .5, 0, AXION_BL);
#undef BAD_INIT
    FAIL(axion_model_init(NULL, path, .1, .2, 1, 0, .5, 0, AXION_BL));
    FAIL(axion_model_init(&m, "/missing/axion.csv", .1, .2, 1, 0, .5, 0, AXION_BL));
    CHECK(m.table == NULL);
    axion_model_free(NULL);
    FAIL(axion_field_value(NULL, NULL, &v));
    FAIL(axion_field_value(&m, NULL, NULL));
    CHECK(strlen(axion_error(AXION_ERR_DOMAIN)) > 0);
    CHECK(strlen(axion_error(-1)) > 0);
    puts("PASS zero baseline bypass, validation, unsupported coordinate failfast");
}

static void test_field_formula(void)
{
    AxionModel m, other;
    double X[4] = {-1022, 5, PI/2, 0};
    double v, w = exp(-pow(5.0/25.0,15.0));
    int phase;
    fixture(rectangular);
    for (phase = 0; phase < 16; ++phase) {
        init_model(&m, .2, .017, 2*PI, phase, .9375, 0, AXION_BL);
        close_to(field(&m, X), 5 * cos(phase*PI/8) * w);
        axion_model_free(&m);
    }
    init_model(&m, .2, .017, 2*PI, 0, .9375, 0, AXION_BL);
    X[0] += 17;
    close_to(field(&m, X), 5*cos(17*.017)*w);
    X[0] += 2*PI/.017;
    close_to(field(&m, X), 5*cos(17*.017)*w); /* omega sets period, NOT alpha */
    init_model(&other, .2, .031, 2*PI, 0, .9375, 0, AXION_BL);
    CHECK(fabs(field(&m, X) - field(&other, X)) > .01);
    axion_model_free(&other);
    X[0] = -1022;
    init_model(&other, .1, .017, 2*PI, 0, .9375, 0, AXION_BL);
    close_to(field(&other, X), 3*w); /* alpha changes R, NOT frequency */
    axion_model_free(&other);
    X[1] = 25;
    close_to(field(&m, X), 7*exp(-1.0)); /* exact W25,15 retained */
    X[1] = 3;
    close_to(field(&m, X), 3*exp(-pow(3.0/25.0,15)));
    X[2] = 0;
    close_to(field(&m, X), 0);
    X[2] = PI/2; X[3] = PI;
    CHECK(field(&m, X) < 0); /* negative FIELD is valid, R itself is nonnegative */
    X[3] = 0; X[1] = 26; v = 987;
    CHECK(exp(-pow(26.0/25.0,15)) > 0);
    CHECK(axion_field_value(&m, X, &v) == AXION_ERR_DOMAIN); CHECK(v == 987);
    X[1] = 35;
    CHECK(exp(-pow(35.0/25.0,15)) > 0);
    FAIL(axion_field_value(&m, X, &v)); /* tiny nonzero window is NOT a cutoff */
    X[1] = 100;
    CHECK(exp(-pow(100.0/25.0,15)) == 0);
    close_to(field(&m, X), 0); /* exact zero window beyond table: no query */
    X[1] = DBL_MAX;
    close_to(field(&m, X), 0);
    X[1] = 1 + sqrt((1-.9375)*(1+.9375));
    FAIL(axion_field_value(&m, X, &v));
    X[1] = 1; FAIL(axion_field_value(&m, X, &v));
    X[1] = 5; X[0] = NAN; FAIL(axion_field_value(&m, X, &v));
    X[0] = DBL_MAX; /* finite t+omega need not overflow with small omega */
    CHECK(isfinite(field(&m, X)));
    FAIL(axion_field_value(&m, NULL, &v));
    FAIL(axion_field_value(&m, X, NULL));
    axion_model_free(&m);
    fixture(".2,3,1\n.2,5,1\n");
    init_model(&m, .2, .17, 1, 0, 0, 0, AXION_BL);
    axion_model_free(&m);
    FAIL(axion_model_init(&m, path, nextafter(.2,1), .17,1,0,0,0,AXION_BL));
    puts("PASS field alpha/omega independence, 16 phases, fixed window/t_ref, domains");
}

static void test_window_profile(void)
{
    AxionModel m;
    const double radii[] = {3, 15, 20, 24, 25, 26, 28, 30, 35};
    size_t i;
    fixture(".2,3,1\n.2,35,1\n");
    init_model(&m, .2, .017, 2*PI, 0, .9375, 0, AXION_BL);
    for (i = 0; i < sizeof(radii)/sizeof(radii[0]); ++i) {
        double X[4] = {-1022, radii[i], PI/2, 0};
        double expected = exp(-pow(radii[i]/25,15));
        double actual = field(&m, X);
        CHECK(expected > 0 && actual > 0);
        CHECK(fabs(actual/expected-1) < 1e-14);
    }
    axion_model_free(&m);
    puts("PASS unchanged W25,15 including tiny nonzero far-window values");
}

static void ks_shifts(double r, double spin, double *f, double *g)
{
    /* Independent direct expression, as in existing GRmath.c. */
    double s = sqrt(1-spin*spin);
    double ratio = log((r-(1+s))/(r-(1-s)));
    *f = log(r*r-2*r+spin*spin) + ratio/s;
    *g = spin/(2*s)*ratio;
}

static void test_covariance(void)
{
    const double spins[] = {0, .5, -.7, .9375, .99, .999999};
    const double radii[] = {3, 4, 5, 10, 25};
    size_t i, j;
    fixture(rectangular);
    for (i = 0; i < sizeof(spins)/sizeof(spins[0]); ++i) {
        AxionModel bl, ks, bhac;
        init_model(&bl, .2, .231, .7, 3, spins[i], .4, AXION_BL);
        init_model(&ks, .2, .231, .7, 3, spins[i], .4, AXION_KS);
        init_model(&bhac, .2, .231, .7, 3, spins[i], .4, AXION_MKSBHAC);
        for (j = 0; j < sizeof(radii)/sizeof(radii[0]); ++j) {
            double r = radii[j], f, g, x2 = .61;
            double theta = x2+.5*.4*sin(2*x2);
            double Xbl[4] = {-1007, r, theta, .8};
            double Xks[4], Xbhac[4];
            ks_shifts(r, spins[i], &f, &g);
            Xks[0] = Xbl[0]+f; Xks[1] = r;
            Xks[2] = theta; Xks[3] = Xbl[3]+g;
            memcpy(Xbhac, Xks, sizeof(Xks));
            Xbhac[1] = log(r); Xbhac[2] = x2;
            close_to(field(&bl, Xbl), field(&ks, Xks));
            close_to(field(&bl, Xbl), field(&bhac, Xbhac));
        }
        axion_model_free(&bl); axion_model_free(&ks); axion_model_free(&bhac);
    }
    /* A representable point immediately outside Schwarzschild horizon. */
    fixture(".2,2,1\n.2,3,1\n");
    {
        AxionModel ks;
        double X[4] = {0, nextafter(2, INFINITY), .5, .3};
        init_model(&ks, .2, .231, 1, 3, 0, 0, AXION_KS);
        CHECK(isfinite(field(&ks, X)));
        X[1] = 2;
        { double v; FAIL(axion_field_value(&ks, X, &v)); }
        axion_model_free(&ks);
    }
    /* Near-extremal spin, at the closest representable exterior radius. */
    fixture(".2,1,1\n.2,3,1\n");
    {
        const double near_extremal[] = {nextafter(1,0), nextafter(-1,0), .9375};
        for (i = 0; i < sizeof(near_extremal)/sizeof(near_extremal[0]); ++i) {
            AxionModel ks;
            double spin = near_extremal[i];
            double rp = 1+sqrt((1-spin)*(1+spin));
            double X[4] = {0, nextafter(rp, INFINITY), .5, .3};
            init_model(&ks, .2, .231, 1, 3, spin, 0, AXION_KS);
            CHECK(isfinite(field(&ks, X)));
            axion_model_free(&ks);
        }
    }
    puts("PASS BL/KS/MKSBHAC same-event covariance, background spin, horizon limit");
}

static void test_parallel_reads(void)
{
    AxionModel m;
    unsigned char before[sizeof(m)];
    double X[4] = {-1015, 5, .7, .2};
    double reference;
    int i, failures = 0;
    fixture(rectangular);
    init_model(&m, .2, .027, .5, 9, .9375, .2, AXION_KS);
    memcpy(before, &m, sizeof(m));
    reference = field(&m, X);
    /* Prove queries do not reopen the CSV. */
    CHECK(remove(path) == 0);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:failures)
#endif
    for (i = 0; i < 100000; ++i) {
        double v = -1, rvalue = -1, invalid = 99;
        size_t na, nr;
        if (axion_field_value(&m, X, &v) || v != reference ||
            axion_table_interpolate(m.table, .2, 4, &rvalue) || fabs(rvalue-4) > 1e-14 ||
            axion_table_dimensions(m.table, &na, &nr) || na != 2 || nr != 3 ||
            axion_table_interpolate(m.table, .4, 4, &invalid) != AXION_ERR_DOMAIN ||
            invalid != 99 || !axion_error(AXION_ERR_DOMAIN)) ++failures;
    }
    CHECK(failures == 0);
    CHECK(memcmp(before, &m, sizeof(m)) == 0);
    axion_model_free(&m);
#ifdef _OPENMP
    printf("PASS immutable parallel reads (%d OpenMP max threads)\n", omp_get_max_threads());
#else
    puts("PASS immutable serial reads (compile with -fopenmp for concurrent test)");
#endif
}

static void test_existing_table(const char *name)
{
    AxionTable *table;
    size_t na, nr;
    OK(axion_table_load(name, &table));
    OK(axion_table_dimensions(table, &na, &nr));
    printf("PASS read-only existing CSV %s: %lu x %lu\n", name,
           (unsigned long)na, (unsigned long)nr);
    axion_table_free(table);
}

int main(int argc, char **argv)
{
    int i, length;
    if (argc < 2) { fprintf(stderr, "usage: %s EXISTING_SCRATCH_DIR [read-only CSV ...]\n", argv[0]); return 2; }
    length = snprintf(path, sizeof(path), "%s/test_axion.csv", argv[1]);
    CHECK(length > 0 && (size_t)length < sizeof(path));
    test_tables();
    test_bad_tables();
    test_dynamic_table();
    test_baseline_and_parameters();
    test_field_formula();
    test_window_profile();
    test_covariance();
    test_parallel_reads();
    for (i = 2; i < argc; ++i) test_existing_table(argv[i]);
    puts("All axion tests passed.");
    return 0;
}
