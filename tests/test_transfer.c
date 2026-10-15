/* Synthetic integration test of the REAL transfer -> plasma -> tetrad -> axion
 * chain. Only upstream metric/GRmath/emission/model dependencies are replaced.
 * No transfer, tetrad, field, or rotation algorithm is copied into this file.
 */
#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_global_vars.h"
#include "axion_transport.h"

#include <stdint.h>

static const char *test_name = "setup";
#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL [%s] line %d: %s\n", test_name, __LINE__, #c); \
    exit(1); } } while (0)

/* Units chosen so dl_scaled*C_CONST=dl_geometric at frequency=1. */
double L_unit = 1., B_unit = 1.;
double C_CONST = PLANCK_CONSTANT /
    (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

enum FluidMode { NO_FLUID_UNTOUCHED, NO_FLUID_GARBAGE, NE_ZERO,
                 SIGMA_CUT_ZERO, PITCH_ALIGNED, ALL_PLASMA };
static enum FluidMode fluid_mode;
static int fluid_calls, vacuum_calls, pitch_calls, freq_calls, coeff_calls;
static int tetrad_calls, field_calls, connection_calls, normalize_calls;
static int mutate_lookup;
static double absorption;

static void reset(enum FluidMode mode)
{
    fluid_mode = mode;
    fluid_calls = vacuum_calls = pitch_calls = freq_calls = coeff_calls = 0;
    tetrad_calls = field_calls = connection_calls = normalize_calls = 0;
}

void metric_dd(double X[4], double g[4][4])
{
    (void)X;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            g[i][j] = i == j ? (i == 0 ? -1. : 1.) : 0.;
}
void metric_uu(double X[4], double g[4][4]) { metric_dd(X, g); }
void lower_index(double X[4], double u[4], double d[4])
{
    (void)X;
    for (int i = 0; i < 4; ++i) d[i] = (i == 0 ? -1. : 1.) * u[i];
}
void raise_index(double X[4], double d[4], double u[4])
{
    lower_index(X, d, u);
}
double inner_product(double *X, double *u, double *v)
{
    (void)X;
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3];
}
double four_velocity_norm(double X[4], double u[4])
{
    return inner_product(X, u, u);
}
void normalize_null(double X[4], double u[4])
{
    (void)X;
    ++normalize_calls;
    u[0] = sqrt(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
}
double get_r(double X[4]) { return X[1]; }
void connection_udd(double X[4], double g[4][4][4])
{
    (void)X;
    ++connection_calls;
    memset(g, 0, 4*4*4*sizeof(double));
}
void connection_num_udd(double X[4], double g[4][4][4])
{
    connection_udd(X, g);
}

/* Fluid output is complete and finite whenever success is returned. Sigma cut
 * models the upstream contract n_e=0, but deliberately returns 1 to verify the
 * transfer density gate independently of get_fluid_params' return value. */
int get_fluid_params(double X[4], struct GRMHD *m)
{
    const int source = fabs(X[1] - 8.) < 1e-12;
    ++fluid_calls;
    if (!source && fluid_mode != ALL_PLASMA) ++vacuum_calls;
    if (!source && fluid_mode == NO_FLUID_UNTOUCHED) return 0;
    if (!source && fluid_mode == NO_FLUID_GARBAGE) {
        for (int i = 0; i < 4; ++i)
            m->U_u[i] = m->U_d[i] = m->B_u[i] = m->B_d[i] = NAN;
        m->n_e = m->B = m->theta_e = m->sigma = NAN;
        m->sigma_min = m->dx_local = m->beta = NAN;
        m->igrid_c = -999;
        if (mutate_lookup) { X[0] += 3.; X[1] = 1.; X[3] += .7; }
        return 0;
    }
    memset(m, 0, sizeof(*m));
    m->U_u[0] = 1.; m->U_d[0] = -1.;
    m->B_u[3] = m->B_d[3] = -1.;
    m->n_e = m->B = m->theta_e = 1.;
    m->sigma = .1; m->sigma_min = 1.;
    m->dx_local = .25; m->beta = 1.; m->igrid_c = 0;
    if (!source && (fluid_mode == NE_ZERO || fluid_mode == SIGMA_CUT_ZERO)) {
        m->n_e = 0.;
        if (fluid_mode == SIGMA_CUT_ZERO) m->sigma = 2.*SIGMA_CUT;
    }
    if (!source && fluid_mode == PITCH_ALIGNED) {
        m->B_u[1] = m->B_d[1] = -1.;
        m->B_u[3] = m->B_d[3] = 0.;
    }
    return 1;
}
double pitch_angle(double *X, double *k, double *B, double *U)
{
    (void)X;
    ++pitch_calls;
    CHECK(isfinite(inner_product(X, U, U)));
    CHECK(fabs(inner_product(X, U, U) + 1.) < 1e-12);
    CHECK(isfinite(k[0]) && k[0] > 0.);
    return B[1] == -1. ? 0. : M_PI/2.;
}
double freq_in_plasma_frame(double U[4], double k[4])
{
    ++freq_calls;
    CHECK(isfinite(U[0]) && isfinite(k[0]) && k[0] < 0.);
    return 1.;
}
static void coefficient_args(double te, double ne, double nu, double B,
                             double pitch)
{
    ++coeff_calls;
    CHECK(isfinite(te) && te > 0. && isfinite(ne) && ne > 0.);
    CHECK(isfinite(nu) && nu == 1. && isfinite(B) && B > 0.);
    CHECK(isfinite(pitch) && pitch > 0.);
}
#define EMISSION(name, value) \
    double name(double te, double ne, double nu, double B, double pitch) { \
        coefficient_args(te, ne, nu, B, pitch); return (value); }
EMISSION(j_I, 2.)
EMISSION(j_Q, .6)
EMISSION(j_V, .2)
EMISSION(rho_Q, absorption*.4)
EMISSION(rho_V, absorption*.7)
#define ABSORPTION(name, value) \
    double name(double te, double ne, double nu, double B, double pitch, double j) { \
        coefficient_args(te, ne, nu, B, pitch); CHECK(isfinite(j)); return (value); }
ABSORPTION(a_I, absorption)
ABSORPTION(a_Q, absorption*.2)
ABSORPTION(a_V, absorption*.1)

/* Linker wrappers observe calls, then delegate to the real implementation.
 * create_observer_tetrad is in a separate production object, so --wrap sees
 * both plasma and camera calls. Field evaluation likewise remains real. */
void __real_create_observer_tetrad(double *, double *, double *, double *, double [][4]);
void __wrap_create_observer_tetrad(double *X, double *k, double *U, double *B,
                                  double t[][4])
{
    ++tetrad_calls;
    __real_create_observer_tetrad(X, k, U, B, t);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) CHECK(isfinite(t[i][j]));
}
int __real_axion_field_value(const AxionModel *, const double [4], double *);
int __wrap_axion_field_value(const AxionModel *m, const double X[4], double *v)
{
    ++field_calls;
    return __real_axion_field_value(m, X, v);
}

static double *make_ray(int n)
{
    /* Exactly n records; ASan catches any access to record n or -1. */
    double *ray = malloc((size_t)n * 9 * sizeof(*ray));
    CHECK(ray != NULL);
    const double dl = 4./(n-1);
    for (int p = 0; p < n; ++p) {
        double *x = ray + 9*p;
        x[0] = -p*dl; x[1] = 4. + p*dl;
        x[2] = M_PI/2.; x[3] = .23;
        x[4] = 1.; x[5] = -1.; x[6] = x[7] = 0.;
        x[8] = -dl; /* Production must take fabs of backward ray dl. */
    }
    /* There is no segment beyond the source. Reading this dl is a bug. */
    ray[9*(n-1)+8] = NAN;
    return ray;
}
static void transfer(double *ray, int n, const AxionModel *m, double out[4])
{
    double tau = 0., tauF = 0.;
    radiative_transfer_polarized(ray, n, 1., NULL, NULL, NULL, 0,
                                out, &tau, &tauF, m, 7, 11);
    for (int i = 0; i < 4; ++i) CHECK(isfinite(out[i]));
    CHECK(isfinite(tau) && isfinite(tauF));
    if (absorption == 0.) CHECK(tau == 0. && tauF == 0.);
}
static void close_real(double actual, double expected, double scale)
{
    CHECK(isfinite(actual) && isfinite(expected));
    if (fabs(actual - expected) > 3e-10 * scale) {
        fprintf(stderr, "actual %.17g expected %.17g scale %.17g\n",
                actual, expected, scale);
        CHECK(0);
    }
}
static void counters(int n, enum FluidMode mode, int active)
{
    const int plasma = mode == ALL_PLASMA ? n-1 : 1;
    CHECK(fluid_calls == n-1);
    CHECK(vacuum_calls == (mode == ALL_PLASMA ? 0 : n-2));
    CHECK(pitch_calls == (mode == PITCH_ALIGNED ? n-1 : plasma));
    CHECK(freq_calls == plasma);
    CHECK(coeff_calls == 8*plasma);
    CHECK(tetrad_calls == plasma+1); /* plasma plus final camera */
    CHECK(field_calls == (active ? 2*(n-1) : 0));
    CHECK(connection_calls == 4*(n-1)); /* real GR RK4 */
    CHECK(normalize_calls == 0);
}
static double endpoint_delta(const AxionModel *m, double *ray, int n)
{
    double observer, source;
    CHECK(axion_field_value(m, ray, &observer) == AXION_OK);
    CHECK(axion_field_value(m, ray + 9*(n-1), &source) == AXION_OK);
    return observer-source;
}
static void check_rotation(const double base[4], const double turned[4], double delta)
{
    const double complex P0 = base[1] + I*base[2];
    const double complex P = turned[1] + I*turned[2];
    CHECK(cabs(P0) > 0.);
    const double complex ratio = P/P0;
    const double complex expected = cexp(2.*I*delta);
    close_real(creal(ratio), creal(expected), 1.);
    close_real(cimag(ratio), cimag(expected), 1.);
    close_real(turned[0], base[0], fabs(base[0]));
    close_real(turned[3], base[3], fabs(base[0]));
    close_real(turned[1]*turned[1]+turned[2]*turned[2],
               base[1]*base[1]+base[2]*base[2], base[0]*base[0]);
}
static void empty_rays(const AxionModel *active)
{
    test_name = "steps 0/1 must not dereference lightpath or call dependencies";
    for (int n = 0; n <= 1; ++n) {
        double out[4] = {9.,9.,9.,9.};
        reset(NO_FLUID_GARBAGE);
        /* Intentionally inaccessible even for steps=1: no segment exists. */
        transfer(n ? (double *)(uintptr_t)1 : NULL, n, active, out);
        for (int i = 0; i < 4; ++i) CHECK(out[i] == 0.);
        CHECK(fluid_calls == 0 && pitch_calls == 0 && coeff_calls == 0);
        CHECK(tetrad_calls == 0 && field_calls == 0 && connection_calls == 0);
    }
}
static void one_source(int n, enum FluidMode mode, const AxionModel *active,
                       const AxionModel *disabled)
{
    double *ray = make_ray(n), base[4], off[4], turned[4];
    double delta = endpoint_delta(active, ray, n);
    CHECK(fabs(delta) > .01 && fabs(sin(2.*delta)) > .01);
    reset(mode); transfer(ray, n, NULL, base); counters(n, mode, 0);
    /* Independent constant-emission expectation, not a duplicate solver. */
    close_real(base[0], 2.*4./(n-1), base[0]);
    close_real(hypot(base[1],base[2]), .6*4./(n-1), base[0]);
    close_real(base[3], .2*4./(n-1), base[0]);
    reset(mode); transfer(ray, n, disabled, off); counters(n, mode, 0);
    for (int i = 0; i < 4; ++i) CHECK(off[i] == base[i]);
    reset(mode); transfer(ray, n, active, turned); counters(n, mode, 1);
    check_rotation(base, turned, delta);
    free(ray);
}
static void plasma_baseline(const AxionModel *disabled)
{
    test_name = "norm=0 leaves nonzero upstream plasma absorption/Faraday unchanged";
    double *ray = make_ray(9), base[4], off[4];
    absorption = .1;
    reset(ALL_PLASMA); transfer(ray, 9, NULL, base); counters(9, ALL_PLASMA, 0);
    reset(ALL_PLASMA); transfer(ray, 9, disabled, off); counters(9, ALL_PLASMA, 0);
    CHECK(base[0] > 0. && base[0] < 8.);
    for (int i = 0; i < 4; ++i) CHECK(off[i] == base[i]);
    absorption = 0.;
    free(ray);
}
int main(int argc, char **argv)
{
    AxionModel active, disabled;
    CHECK(argc == 2);
    /* CSV fixture is synthetic, generated in the runner's private /tmp dir.
     * alpha is intentionally NOT omega; no inference from mass is allowed. */
    CHECK(axion_model_init(&active, argv[1], .2, .137, 3., 5., 0., 0.,
                           AXION_BL) == AXION_OK);
    CHECK(active.alpha != active.omega);
    /* Baseline needs no table; malicious-path parameter validation belongs to
     * test_axion.c. Here we test its interaction with the production chain. */
    CHECK(axion_model_init(&disabled, NULL, .2, .137, 0., 5., 0., 0.,
                           AXION_BL) == AXION_OK);
    CHECK(disabled.table == NULL);
    empty_rays(&active);
    test_name = "2 exactly allocated records, one forward segment";
    one_source(2, NO_FLUID_UNTOUCHED, &active, &disabled);
    test_name = "source then no-fluid return 0 without resetting output";
    one_source(7, NO_FLUID_UNTOUCHED, &active, &disabled);
    test_name = "no-fluid garbage and mutated lookup coordinates cannot leak";
    mutate_lookup = 1;
    one_source(7, NO_FLUID_GARBAGE, &active, &disabled);
    mutate_lookup = 0;
    test_name = "fluid success but n_e=0 is vacuum";
    one_source(7, NE_ZERO, &active, &disabled);
    test_name = "sigma-cut n_e=0 remains vacuum even with fluid success";
    one_source(7, SIGMA_CUT_ZERO, &active, &disabled);
    test_name = "pitch-aligned early plasma return still rotates every segment";
    one_source(7, PITCH_ALIGNED, &active, &disabled);
    test_name = "max_steps exactly allocated, no p+1 or p-1 overrun";
    one_source((int)max_steps, NO_FLUID_GARBAGE, &active, &disabled);
    plasma_baseline(&disabled);
    axion_model_free(&active);
    axion_model_free(&disabled);
    puts("PASS: real polarized transfer/tetrad/axion synthetic chain");
    return 0;
}
