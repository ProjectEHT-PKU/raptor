/* Actual MKSBHAC metric + tetrad, synthetic null rays only; no fluid/HDF IO. */
#include "definitions.h"
#include "functions.h"
#include "axion_transport.h"

/* Only the globals needed by the retained metric/tetrad/link sections. */
double a = .9375, hslope = 0., R0 = 0., B_unit = 1.;

static void require(int condition, const char *what) {
    if (!condition) { fprintf(stderr, "FAIL curved screen: %s\n", what); exit(1); }
}

static double complex project(double g[4][4], double basis[4][4], int axis,
                               double complex f[4]) {
    double complex result = 0;
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
        result += g[i][j] * basis[i][axis] * f[j];
    return result;
}

int main(void) {
    const double radii[] = {1.4, 2., 3., 5., 10., 50., 240., 10000.};
    const double thetas[] = {.4, 1.2, 2.5};
    int cases = 0;
    for (int spin = 0; spin < 3; spin++) {
        a = spin == 0 ? 0. : spin == 1 ? .9375 : -.9375;
        for (int slope = 0; slope < 2; slope++) {
            hslope = slope ? .2 : 0.;
            for (size_t r = 0; r < sizeof(radii)/sizeof(*radii); r++) {
                if (radii[r] <= 1. + sqrt(1. - a*a)) continue;
                for (size_t th = 0; th < sizeof(thetas)/sizeof(*thetas); th++) {
                for (int direction_index = 0; direction_index < 3; direction_index++) {
                    double X[4] = {0., log(radii[r]), thetas[th], .7};
                    double gu[4][4], gd[4][4], U[4], Ud[4];
                    metric_uu(X, gu); metric_dd(X, gd);
                    for (int i = 0; i < 4; i++) U[i] = -gu[i][0]/sqrt(-gu[0][0]);
                    lower_index(X, U, Ud);
                    /* Radial and mixed spatial vectors in the observer's rest space. */
                    double trial[4] = {0., direction_index == 2 ? -1. : 1.,
                                       direction_index == 0 ? 0. : .13,
                                       direction_index == 0 ? 0. : -.21};
                    double beta = 0.;
                    for (int i = 0; i < 4; i++) beta += Ud[i]*trial[i];
                    double n[4], norm = 0.;
                    for (int i = 0; i < 4; i++) n[i] = trial[i] + beta*U[i];
                    norm = sqrt(inner_product(X, n, n));
                    double k[4];
                    for (int i = 0; i < 4; i++) k[i] = U[i] + n[i]/norm;
                    require(fabs(inner_product(X, k, k)) < 3e-11, "null ray");
                    double direction[4] = {0., 0., 1., 0.}, tetrad[4][4];
                    create_observer_tetrad(X, k, U, direction, tetrad);
                    double complex f[4], gauged[4];
                    /* Elliptic Jones state exercises nonzero Q, U and V. */
                    const double complex e1 = .8, e2 = .36 + .48*I;
                    for (int i = 0; i < 4; i++) {
                        f[i] = tetrad[i][1]*e1 + tetrad[i][2]*e2;
                        gauged[i] = f[i] + (.3 + .4*I)*k[i];
                    }
                    const double delta = .237;
                    require(axion_rotate_polarization(X,k,delta,f), "rotation status");
                    require(axion_rotate_polarization(X,k,delta,gauged), "gauge rotation status");
                    double complex f1 = project(gd,tetrad,1,f), f2 = project(gd,tetrad,2,f);
                    double complex q1 = project(gd,tetrad,1,gauged), q2 = project(gd,tetrad,2,gauged);
                    double Q0 = creal(e1*conj(e1)-e2*conj(e2));
                    double U0 = 2.*creal(conj(e1)*e2), V0 = -2.*cimag(conj(e1)*e2);
                    double Q = creal(f1*conj(f1)-f2*conj(f2));
                    double Uq = 2.*creal(conj(f1)*f2), V = -2.*cimag(conj(f1)*f2);
                    require(cabs((Q+I*Uq)-(Q0+I*U0)*cexp(2.*I*delta)) < 3e-10, "spin-2 sign/basis");
                    require(fabs(V-V0) < 3e-10, "V invariant");
                    require(fabs(creal(f1*conj(f1)+f2*conj(f2))-1.) < 3e-10, "norm invariant");
                    require(cabs(f1-q1)+cabs(f2-q2) < 3e-10, "screen gauge invariant");
                    require(axion_rotate_polarization(X,k,-delta,f), "reverse status");
                    require(cabs(project(gd,tetrad,1,f)-e1)+cabs(project(gd,tetrad,2,f)-e2) < 3e-10, "reverse rotation");
                    cases++;
                }
                }
            }
        }
    }
    printf("PASS curved MKSBHAC screen: %d spin/slope/radius/theta cases\n", cases);
    return 0;
}
