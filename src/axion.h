#ifndef RAPTOR_AXION_H
#define RAPTOR_AXION_H

#include <stddef.h>

/* Standalone C99, M=1. No RAPTOR, HDF, RNG, or lazy initialization dependency. */
typedef struct AxionTable AxionTable;

typedef enum {
    AXION_BL = 0,       /* X = (t_BL, r, theta, phi_BL), linear r */
    AXION_MKSBHAC,      /* X = (t_KS, log(r), x2, phi_KS), R0=0;
                          theta=x2+hslope*sin(2*x2)/2 */
    AXION_MKS,         /* reserved: unsupported until R0 contract is specified */
    AXION_MKSHARM,     /* reserved: unsupported until R0/phi contract is specified */
    AXION_KS,          /* X = (t_KS, r, theta, phi_KS), linear r */
    AXION_CKS          /* unsupported */
} AxionCoordinates;

typedef enum {
    AXION_OK = 0,
    AXION_ERR_ARGUMENT,
    AXION_ERR_IO,
    AXION_ERR_MEMORY,
    AXION_ERR_FORMAT,
    AXION_ERR_GRID,
    AXION_ERR_DOMAIN,
    AXION_ERR_COORDINATES,
    AXION_ERR_NUMERIC
} AxionStatus;

/* Owns table. Treat all members as immutable after successful initialization.
 * Do not copy an owning model or mutate it while readers are active. */
typedef struct {
    const AxionTable *table;
    double alpha;       /* mass used ONLY for radial-table lookup */
    double omega;       /* independent positive Re(M omega), for time phase */
    double norm;        /* nonnegative N, field normalization N/(2*pi) */
    double phase;       /* integer index 0..15, phase radians = phase*pi/8 */
    double spin;        /* background metric a, not inferred from the CSV */
    double hslope;
    AxionCoordinates coords;
} AxionModel;

/* All functions return AXION_OK on success (except free/error).
 * Error text depends only on status; no shared or per-model writable error state.
 * Queries write only the caller's output, and leave it unchanged on failure. */
const char *axion_error(int status);

/* Load headerless alpha,r,R CSV, sorted by alpha then r. Whitespace/CRLF are
 * allowed, but blank rows, extra columns and trailing garbage are rejected.
 * Require nalpha>=1, nr>=2, finite entries, R>=0, strictly increasing alpha,
 * and an identical strictly increasing r grid in every slice. No renormalizing.
 * On failure *out=NULL. Call load/free only outside concurrent query regions.
 * Numeric parsing uses the caller's C numeric locale; never change locale while
 * loading. The table contains no spin metadata; provenance is caller-owned. */
int axion_table_load(const char *path, AxionTable **out);
void axion_table_free(const AxionTable *table);
int axion_table_dimensions(const AxionTable *table, size_t *nalpha, size_t *nr);
int axion_table_interpolate(const AxionTable *table, double alpha, double r,
                            double *value);

/* Initialize single-threaded BEFORE ray/OpenMP regions. An uninitialized or
 * freed destination is allowed; free a live model before reinitializing it.
 * Any failed init leaves a zeroed, safely freeable model.
 * norm=0 loads nothing and ignores path/alpha/omega/phase/spin/hslope/coords.
 * Otherwise alpha,omega>0; integer phase in [0,15]; |spin|<1; finite hslope;
 * BL/KS/MKSBHAC only. alpha must be in the table domain (single slice: exact).
 * Fixed historical choices t_ref=1022, window radius=25 and power=15 remain.
 * All nonzero-norm evaluations require r > 1+sqrt(1-spin^2).
 * Outside the table: error unless exp(-(r/25)^15) evaluates to EXACT zero.
 * No extrapolation or approximate small-window cutoff is used. */
int axion_model_init(AxionModel *model, const char *table_path,
                     double alpha, double omega, double norm, double phase,
                     double spin, double hslope, AxionCoordinates coords);
void axion_model_free(AxionModel *model);
int axion_field_value(const AxionModel *model, const double X[4], double *value);

#endif
