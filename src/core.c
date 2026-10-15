/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include "time.h"

// GLOBAL VARS
//////////////

char GRMHD_FILE[256];

double MBH, M_UNIT, TIME_INIT, INCLINATION;
double axion_alpha, axion_omega, axion_norm;
int phase_index;
char axion_table[256];
AxionModel axion_model;
double R_HIGH, R_LOW;
double FREQS_PER_DEC, FREQ_MIN, FREQ_MAX;

double SOURCE_DIST; // Distance to M87 (cm); for Sgr A* use (2.47e22)

int IMG_WIDTH, IMG_HEIGHT;
double CAM_SIZE_X, CAM_SIZE_Y;
double STEPSIZE;

// FUNCTIONS
////////////

/* Whitespace-delimited named triples also accept old scripts' newline-free
 * serialization, but never their ambiguous omega-only axion semantics. */
static void input_error(FILE *input, const char *path, const char *field,
                        const char *reason) {
    fprintf(stderr, "Input error in %s: %s: %s\n", path, field, reason);
    if (input != NULL)
        fclose(input);
    exit(EXIT_FAILURE);
}

/* Unlike a width-limited scanf, reject rather than split an overlong token. */
static int input_token(FILE *input, char token[256]) {
    int c;
    size_t n = 0;
    do {
        c = fgetc(input);
    } while (c != EOF && isspace((unsigned char)c));
    if (c == EOF)
        return ferror(input) ? -1 : 0;
    do {
        if (n == 255 || c == '\0')
            return -1;
        token[n++] = (char)c;
        c = fgetc(input);
    } while (c != EOF && !isspace((unsigned char)c));
    token[n] = '\0';
    return ferror(input) ? -1 : 1;
}

static double input_number(FILE *input, const char *path, const char *field,
                           const char *token) {
    char *end;
    errno = 0;
    double value = strtod(token, &end);
    if (errno != 0 || end == token || *end != '\0' || !isfinite(value))
        input_error(input, path, field, "expected a finite number (no trailing text)");
    return value;
}

static int input_integer(FILE *input, const char *path, const char *field,
                         const char *token) {
    char *end;
    errno = 0;
    long value = strtol(token, &end, 10);
    if (errno != 0 || end == token || *end != '\0' ||
        value < INT_MIN || value > INT_MAX)
        input_error(input, path, field, "expected an in-range integer");
    return (int)value;
}

void read_model(char *argv[]) {
    enum { REAL, INTEGER, PATH };
    struct InputField {
        const char *name, *unit;
        int type;
        void *value;
        double min, max;
        int seen;
    } fields[] = {
        {"MBH", "(Msun)", REAL, &MBH, 0, HUGE_VAL, 0},
        {"DISTANCE", "(kpc)", REAL, &SOURCE_DIST, 0, HUGE_VAL, 0},
        {"M_UNIT", "(g)", REAL, &M_UNIT, 0, HUGE_VAL, 0},
        {"R_LOW", "(-)", REAL, &R_LOW, 0, HUGE_VAL, 0},
        {"R_HIGH", "(-)", REAL, &R_HIGH, 0, HUGE_VAL, 0},
        {"INCLINATION", "(deg)", REAL, &INCLINATION, 0, 180, 0},
        {"IMG_WIDTH", "(pixels)", INTEGER, &IMG_WIDTH, 1, INT_MAX, 0},
        {"IMG_HEIGHT", "(pixels)", INTEGER, &IMG_HEIGHT, 1, INT_MAX, 0},
        {"CAM_SIZE_X", "(Rg)", REAL, &CAM_SIZE_X, 0, HUGE_VAL, 0},
        {"CAM_SIZE_Y", "(Rg)", REAL, &CAM_SIZE_Y, 0, HUGE_VAL, 0},
        {"FREQS_PER_DEC", "(-)", REAL, &FREQS_PER_DEC, 0, HUGE_VAL, 0},
        {"FREQ_MIN", "(Hz)", REAL, &FREQ_MIN, 0, HUGE_VAL, 0},
        {"STEPSIZE", "(-)", REAL, &STEPSIZE, 0, HUGE_VAL, 0},
        {"AXION_ALPHA", "(-)", REAL, &axion_alpha, 0, HUGE_VAL, 0},
        {"AXION_OMEGA", "(-)", REAL, &axion_omega, 0, HUGE_VAL, 0},
        {"AXION_NORM", "(-)", REAL, &axion_norm, 0, HUGE_VAL, 0},
        {"AXION_PHASE", "(-)", INTEGER, &phase_index, 0, 15, 0},
        {"AXION_TABLE", "(-)", PATH, axion_table, 0, 0, 0},
        {"MAX_LEVEL", "(-)", INTEGER, &max_level, 1, INT_MAX, 0}
    };
    const size_t count = sizeof(fields) / sizeof(fields[0]);
    char label[256], unit[256], token[256];
    if (argv == NULL || argv[1] == NULL || argv[2] == NULL || argv[3] == NULL)
        input_error(NULL, "command line", "arguments",
                    "usage: RAPTOR model.in GRMHD_FILE TIME_INIT");
    const char *path = argv[1];
    fprintf(stdout, "\nUsing model parameter file %s\n", path);
    FILE *input = fopen(path, "r");
    if (input == NULL)
        input_error(NULL, path, "file", "cannot open parameter file");

    int status;
    while ((status = input_token(input, label)) == 1) {
        size_t i;
        for (i = 0; i < count && strcmp(label, fields[i].name) != 0; ++i) {}
        if (i == count)
            input_error(input, path, label,
                        "unknown field; use exact names, including AXION_ALPHA, "
                        "AXION_OMEGA, AXION_NORM, AXION_PHASE and AXION_TABLE. "
                        "Legacy omega-only inputs require migration to independent "
                        "alpha and true omega (see docs/axion-input.md)");
        struct InputField *field = &fields[i];
        if (field->seen)
            input_error(input, path, label, "duplicate field");
        if (input_token(input, unit) != 1 || input_token(input, token) != 1)
            input_error(input, path, label, "incomplete LABEL UNIT VALUE triple or token exceeds 255 bytes");
        if (strcmp(unit, field->unit) != 0)
            input_error(input, path, label, "incorrect unit token (see run/model.in)");
        if (field->type == PATH) {
            memcpy(field->value, token, strlen(token) + 1);
        } else {
            double value = field->type == INTEGER
                ? input_integer(input, path, label, token)
                : input_number(input, path, label, token);
            if (value < field->min || value > field->max)
                input_error(input, path, label, "value outside allowed range (see docs/axion-input.md)");
            if (field->type == INTEGER)
                *(int *)field->value = (int)value;
            else
                *(double *)field->value = value;
        }
        field->seen = 1;
    }
    if (status < 0)
        input_error(input, path, "token", "read error or token exceeds 255 bytes");
    /* Prioritize the migration diagnostic over other missing fields. */
    if (!fields[13].seen)
        input_error(input, path, "AXION_ALPHA",
                    "required mass/coupling parameter missing. Legacy omega-only "
                    "inputs are rejected: provide independent AXION_ALPHA and true "
                    "AXION_OMEGA plus AXION_TABLE; do not copy the old mixed-use "
                    "OMEGA blindly. See docs/axion-input.md; old Run.py needs redesign.");
    for (size_t i = 0; i < count; ++i)
        if (!fields[i].seen)
            input_error(input, path, fields[i].name, "required field missing (also required for AXION_NORM=0)");

    if (MBH <= 0 || SOURCE_DIST <= 0 || M_UNIT <= 0 || R_LOW <= 0 ||
        R_HIGH <= 0 || CAM_SIZE_X <= 0 || CAM_SIZE_Y <= 0 ||
        FREQS_PER_DEC <= 0 || FREQ_MIN <= 0 || STEPSIZE <= 0)
        input_error(input, path, "physical scales", "mass, distance, temperature ratios, camera sizes, frequency and stepsize must be positive");
    if (axion_norm > 0 && (axion_alpha <= 0 || axion_omega <= 0))
        input_error(input, path, "AXION_ALPHA/AXION_OMEGA", "both must be positive for AXION_NORM>0; zero is allowed only for a disabled baseline");
#if (!POL)
    if (axion_norm != 0)
        input_error(input, path, "AXION_NORM", "nonzero axion rotation requires POL=1; refusing to silently ignore it");
#endif
    if (IMG_WIDTH != IMG_HEIGHT || IMG_WIDTH % num_pixels_1d != 0 ||
        IMG_HEIGHT % num_pixels_1d != 0 ||
        (double)(IMG_WIDTH / num_pixels_1d) * (IMG_HEIGHT / num_pixels_1d) > INT_MAX)
        input_error(input, path, "image grid", "must be square, positive multiples of num_pixels_1d, with block count <= INT_MAX");
    if (!isfinite(MBH * MSUN) || !isfinite(SOURCE_DIST * KPCTOCM))
        input_error(input, path, "physical scales", "cgs conversion overflows");
    if (strlen(argv[2]) >= sizeof(GRMHD_FILE))
        input_error(input, path, "GRMHD_FILE", "path exceeds 255 bytes");
    memcpy(GRMHD_FILE, argv[2], strlen(argv[2]) + 1);
    TIME_INIT = input_number(input, path, "TIME_INIT", argv[3]);

    fprintf(stderr, "\nModel parameters:\n\n");
    fprintf(stderr, "MBH \t\t= %g Msun\n", MBH);
    fprintf(stderr, "DISTANCE \t= %g kpc\n", SOURCE_DIST);
    fprintf(stderr, "M_UNIT \t\t= %g grams\n", M_UNIT);
    fprintf(stderr, "R_LOW \t\t= %g \n", R_LOW);
    fprintf(stderr, "R_HIGH \t\t= %g \n", R_HIGH);
    fprintf(stderr, "INCLINATION \t= %g deg\n", INCLINATION);

    fprintf(stderr, "METRIC \t\t= ");
#if (metric == MKSBHAC)
    fprintf(stderr, "MKS BHAC\n");
#elif (metric == BL)
    fprintf(stderr, "Boyer-Lindquist\n");
#elif (metric == KS)
    fprintf(stderr, "Kerr-Schild\n");
#elif (metric == CKS)
    fprintf(stderr, "CKS BHAC\n");
#elif (metric == MKSHARM)
    fprintf(stderr, "MKS HARM3D\n");
#endif

    fprintf(stderr, "\nObserver parameters:\n\n");
    fprintf(stderr, "IMG_WIDTH \t= %d \n", IMG_WIDTH);
    fprintf(stderr, "IMG_HEIGHT \t= %d \n", IMG_HEIGHT);
    fprintf(stderr, "CAM_SIZE_X \t= %g GM/c2\n", CAM_SIZE_X);
    fprintf(stderr, "CAM_SIZE_Y \t= %g GM/c2\n", CAM_SIZE_Y);
    fprintf(stderr, "FREQS_PER_DEC \t= %lf \n", FREQS_PER_DEC);
    fprintf(stderr, "FREQ_MIN \t= %g Hz\n", FREQ_MIN);
    fprintf(stderr, "STEPSIZE \t= %g \n", STEPSIZE);
    fprintf(stderr, "AXION \talpha=%g, omega=%g, norm=%g, phase=%d, table=%s\n",
            axion_alpha, axion_omega, axion_norm, phase_index, axion_table);
    // to cgs units
    MBH *= MSUN;
    SOURCE_DIST *= KPCTOCM;

    fclose(input);
}

// For a single block this function will iterate over the pixels and call
// geodesic integrations as well as radiation transfer
void calculate_image_block(struct Camera *intensityfield,
                           double frequencies[num_frequencies],int block) {
#if (!POL)
    (void)block;
#endif

#pragma omp parallel for shared(frequencies, intensityfield, p)                \
    schedule(static, 1)
    for (int pixel = 0; pixel < tot_pixels; pixel++) {
        int steps = 0;

        double *lightpath2 = malloc(9 * max_steps * sizeof(double));

#if (POL)
        double f_x = 0.;
        double f_y = 0.;
        double p = 0.;
#endif
        // INTEGRATE THIS PIXEL'S GEODESIC
        integrate_geodesic((*intensityfield).alpha[pixel],
                           (*intensityfield).beta[pixel], lightpath2, &steps,
                           CUTOFF_INNER);
        // PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
#if (POL)
        for (int f = 0; f < num_frequencies; f++) {

            radiative_transfer_polarized(lightpath2, steps, frequencies[f],
                                         &f_x, &f_y, &p, 0,
                                         (*intensityfield).IQUV[pixel][f],
                                         &(*intensityfield).tau[pixel][f],
                                         &(*intensityfield).tauF[pixel][f],
                                         &axion_model, block, pixel);
        }

#else
        radiative_transfer_unpolarized(lightpath2, steps, frequencies,
                                       (*intensityfield).IQUV[pixel],
                                       (*intensityfield).tau[pixel]);
        for (int f = 0; f < num_frequencies; f++) {
            (*intensityfield).IQUV[pixel][f][0] *= pow(frequencies[f], 3.);
        }
#endif
        free(lightpath2);
    }
#pragma omp barrier
}

// Functions that computes a spectrum at every frequency
// by integrating over the image struct
void compute_spec(struct Camera *intensityfield,
                  double energy_spectrum[num_frequencies][nspec]) {
    double dA;
#if (POL)
    double S_I, S_Q, S_U, S_V;
#endif

    for (int block = 0; block < tot_blocks; block++) {
        dA = (intensityfield)[block].dx[0] * (intensityfield)[block].dx[1];
        for (int pixel = 0; pixel < tot_pixels; pixel++) {
            for (int freq = 0; freq < num_frequencies; freq++) {
#if (POL)

                S_I = (intensityfield)[block].IQUV[pixel][freq][0];
                S_Q = (intensityfield)[block].IQUV[pixel][freq][1];
                S_U = (intensityfield)[block].IQUV[pixel][freq][2];
                S_V = (intensityfield)[block].IQUV[pixel][freq][3];

                // Stokes I
                energy_spectrum[freq][0] += S_I * dA;

                // Stokes Q
                energy_spectrum[freq][1] += S_Q * dA;

                // Stokes U
                energy_spectrum[freq][2] += S_U * dA;

                // stokes V
                energy_spectrum[freq][3] += S_V * dA;

#else
                energy_spectrum[freq][0] +=
                    (intensityfield)[block].IQUV[pixel][freq][0] * dA;
#endif
            }
        }
    }
}
