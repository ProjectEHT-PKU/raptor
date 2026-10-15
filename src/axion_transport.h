/* Axion-only transport in a deterministic local screen, independent of fluid. */
#ifndef RAPTOR_AXION_TRANSPORT_H
#define RAPTOR_AXION_TRANSPORT_H

#include <complex.h>

/* delta_chi is a forward EVPA increment, not twice that angle. */
int axion_rotate_polarization(const double X[4], const double k[4],
                              double delta_chi, double complex f[4]);

#endif
