#include "axion.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AXION_PI 3.141592653589793238462643383279502884

struct AxionTable {
    size_t nalpha, nr;
    double *alpha, *r, *values;
};

typedef struct { double alpha, r, value; } AxionRow;

const char *axion_error(int status)
{
    switch (status) {
    case AXION_OK: return "success";
    case AXION_ERR_ARGUMENT: return "invalid axion argument or model parameter";
    case AXION_ERR_IO: return "cannot open, read, or close axion table";
    case AXION_ERR_MEMORY: return "axion allocation failed or size overflow";
    case AXION_ERR_FORMAT: return "expected exactly three finite CSV numbers with R>=0";
    case AXION_ERR_GRID: return "axion table must be rectangular with increasing alpha and shared increasing r, nalpha>=1, nr>=2";
    case AXION_ERR_DOMAIN: return "axion query outside table or exterior-horizon domain";
    case AXION_ERR_COORDINATES: return "unsupported axion coordinates (only BL, KS, MKSBHAC supported)";
    case AXION_ERR_NUMERIC: return "nonfinite axion coordinate transformation or field result";
    default: return "unknown axion status";
    }
}

/* Grow buffers without losing the original allocation on realloc failure. */
static void *reserve(void *buffer, size_t *capacity, size_t need, size_t item)
{
    size_t next;
    void *p;
    if (need <= *capacity) return buffer;
    if (need > SIZE_MAX / item) return NULL;
    next = *capacity ? *capacity : 64;
    while (next < need) {
        if (next > SIZE_MAX / 2) { next = need; break; }
        next *= 2;
    }
    if (next > SIZE_MAX / item) next = need;
    p = realloc(buffer, next * item);
    if (!p) return NULL;
    *capacity = next;
    return p;
}

/* A dynamic C99 line reader: no fixed-size truncation, no POSIX getline. */
static int read_line(FILE *fp, char **line, size_t *capacity, int *done)
{
    size_t length = 0;
    int c;
    *done = 0;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '\0') return AXION_ERR_FORMAT;
        if (length > SIZE_MAX - 2) return AXION_ERR_MEMORY;
        {
            char *next = reserve(*line, capacity, length + 2, sizeof(char));
            if (!next) return AXION_ERR_MEMORY;
            *line = next;
        }
        (*line)[length++] = (char)c;
        if (c == '\n') break;
    }
    if (ferror(fp)) return AXION_ERR_IO;
    if (!length) { *done = 1; return AXION_OK; }
    (*line)[length] = '\0';
    return AXION_OK;
}

static int parse_row(const char *line, AxionRow *row)
{
    double numbers[3];
    const char *p = line;
    char *end;
    size_t i;
    for (i = 0; i < 3; ++i) {
        while (isspace((unsigned char)*p)) ++p;
        errno = 0;
        numbers[i] = strtod(p, &end);
        /* Accept representable subnormals (strtod may set ERANGE), but never
         * silently turn a nonzero textual number into zero by underflow. */
        if (end == p || (errno == ERANGE && numbers[i] == 0) ||
            !isfinite(numbers[i]))
            return AXION_ERR_FORMAT;
        p = end;
        while (isspace((unsigned char)*p)) ++p;
        if (i < 2) {
            if (*p != ',') return AXION_ERR_FORMAT;
            ++p;
        } else if (*p != '\0') return AXION_ERR_FORMAT;
    }
    if (numbers[2] < 0) return AXION_ERR_FORMAT;
    row->alpha = numbers[0]; row->r = numbers[1]; row->value = numbers[2];
    return AXION_OK;
}

void axion_table_free(const AxionTable *table)
{
    if (!table) return;
    free(table->alpha);
    free(table->r);
    free(table->values);
    free((void *)table);
}

int axion_table_load(const char *path, AxionTable **out)
{
    FILE *fp;
    AxionRow *rows = NULL;
    AxionTable *table = NULL;
    char *line = NULL;
    size_t line_capacity = 0, capacity = 0, count = 0, nr, na, i, j;
    int status = AXION_OK, done;
    if (!out) return AXION_ERR_ARGUMENT;
    *out = NULL;
    if (!path || !*path) return AXION_ERR_ARGUMENT;
    fp = fopen(path, "r");
    if (!fp) return AXION_ERR_IO;
    for (;;) {
        AxionRow row;
        status = read_line(fp, &line, &line_capacity, &done);
        if (status || done) break;
        status = parse_row(line, &row);
        if (status) break;
        if (count == SIZE_MAX) { status = AXION_ERR_MEMORY; break; }
        {
            AxionRow *next = reserve(rows, &capacity, count + 1, sizeof(*rows));
            if (!next) { status = AXION_ERR_MEMORY; break; }
            rows = next;
        }
        rows[count++] = row;
    }
    if (fclose(fp) != 0 && !status) status = AXION_ERR_IO;
    free(line);
    if (status) goto cleanup;
    if (count < 2) { status = AXION_ERR_GRID; goto cleanup; }
    for (nr = 1; nr < count && rows[nr].alpha == rows[0].alpha; ++nr) {}
    if (nr < 2 || count % nr) { status = AXION_ERR_GRID; goto cleanup; }
    na = count / nr;
    for (i = 0; i < na; ++i) {
        size_t base = i * nr;
        if (i && !(rows[base].alpha > rows[base - nr].alpha)) {
            status = AXION_ERR_GRID; goto cleanup;
        }
        for (j = 0; j < nr; ++j) {
            if (rows[base + j].alpha != rows[base].alpha ||
                rows[base + j].r != rows[j].r ||
                (j && !(rows[base + j].r > rows[base + j - 1].r))) {
                status = AXION_ERR_GRID; goto cleanup;
            }
        }
    }
    /* count is already bounded by SIZE_MAX/sizeof(AxionRow). */
    table = calloc(1, sizeof(*table));
    if (!table) { status = AXION_ERR_MEMORY; goto cleanup; }
    table->nalpha = na; table->nr = nr;
    table->alpha = malloc(na * sizeof(double));
    table->r = malloc(nr * sizeof(double));
    table->values = malloc(count * sizeof(double));
    if (!table->alpha || !table->r || !table->values) {
        status = AXION_ERR_MEMORY; goto cleanup;
    }
    for (i = 0; i < na; ++i) table->alpha[i] = rows[i * nr].alpha;
    for (j = 0; j < nr; ++j) table->r[j] = rows[j].r;
    for (i = 0; i < count; ++i) table->values[i] = rows[i].value;
    *out = table;
    table = NULL;
cleanup:
    free(rows);
    axion_table_free(table);
    return status;
}

int axion_table_dimensions(const AxionTable *table, size_t *nalpha, size_t *nr)
{
    if (!table || !nalpha || !nr) return AXION_ERR_ARGUMENT;
    *nalpha = table->nalpha; *nr = table->nr;
    return AXION_OK;
}

/* Inclusive endpoints, binary search, no out-of-range array access.
 * The one-alpha case is permitted only at its exact stored alpha. */
static int bracket(const double *grid, size_t n, double x, size_t *lo)
{
    size_t left = 0, right = n - 1;
    if (!isfinite(x) || x < grid[0] || x > grid[n - 1])
        return AXION_ERR_DOMAIN;
    if (n == 1) { *lo = 0; return AXION_OK; }
    if (x == grid[n - 1]) { *lo = n - 2; return AXION_OK; }
    while (right - left > 1) {
        size_t mid = left + (right - left) / 2;
        if (x < grid[mid]) right = mid;
        else left = mid;
    }
    *lo = left;
    return AXION_OK;
}

/* Rescale only if an axis subtraction overflows (finite opposite-sign ends). */
static double fraction(double x, double lo, double hi)
{
    double d = hi - lo;
    if (x == lo) return 0;
    if (x == hi) return 1;
    if (isfinite(d)) return (x - lo) / d;
    return (x * 0.5 - lo * 0.5) / (hi * 0.5 - lo * 0.5);
}

static double lerp(double lo, double hi, double weight)
{
    /* R is nonnegative: hi-lo cannot overflow, even at DBL_MAX. */
    if (weight == 0) return lo;
    if (weight == 1) return hi;
    return lo + weight * (hi - lo);
}

int axion_table_interpolate(const AxionTable *table, double alpha, double r,
                            double *value)
{
    size_t i, j;
    double wr, wa, low, high, result;
    int status;
    if (!table || !value) return AXION_ERR_ARGUMENT;
    status = bracket(table->alpha, table->nalpha, alpha, &i);
    if (status) return status;
    status = bracket(table->r, table->nr, r, &j);
    if (status) return status;
    wr = fraction(r, table->r[j], table->r[j + 1]);
    low = lerp(table->values[i * table->nr + j],
               table->values[i * table->nr + j + 1], wr);
    result = low;
    if (table->nalpha > 1) {
        high = lerp(table->values[(i + 1) * table->nr + j],
                    table->values[(i + 1) * table->nr + j + 1], wr);
        wa = fraction(alpha, table->alpha[i], table->alpha[i + 1]);
        result = lerp(low, high, wa);
    }
    if (!isfinite(result) || result < 0) return AXION_ERR_NUMERIC;
    *value = result;
    return AXION_OK;
}

static int model_parameters(const AxionModel *m)
{
    if (!isfinite(m->norm) || m->norm < 0 ||
        !isfinite(m->alpha) || m->alpha <= 0 ||
        !isfinite(m->omega) || m->omega <= 0 ||
        !isfinite(m->phase) || m->phase < 0 || m->phase > 15 ||
        floor(m->phase) != m->phase ||
        !isfinite(m->spin) || fabs(m->spin) >= 1 || !isfinite(m->hslope))
        return AXION_ERR_ARGUMENT;
    if (m->coords != AXION_BL && m->coords != AXION_KS &&
        m->coords != AXION_MKSBHAC) return AXION_ERR_COORDINATES;
    return AXION_OK;
}

int axion_model_init(AxionModel *model, const char *table_path,
                     double alpha, double omega, double norm, double phase,
                     double spin, double hslope, AxionCoordinates coords)
{
    AxionModel next = {0};
    AxionTable *table = NULL;
    size_t index;
    int status;
    if (!model) return AXION_ERR_ARGUMENT;
    memset(model, 0, sizeof(*model));
    if (!isfinite(norm) || norm < 0) return AXION_ERR_ARGUMENT;
    if (norm == 0) return AXION_OK;
    next.alpha = alpha; next.omega = omega; next.norm = norm;
    next.phase = phase; next.spin = spin; next.hslope = hslope; next.coords = coords;
    status = model_parameters(&next);
    if (status) return status;
    status = axion_table_load(table_path, &table);
    if (status) return status;
    status = bracket(table->alpha, table->nalpha, alpha, &index);
    if (status) { axion_table_free(table); return status; }
    next.table = table;
    *model = next;
    return AXION_OK;
}

void axion_model_free(AxionModel *model)
{
    if (!model) return;
    axion_table_free(model->table);
    memset(model, 0, sizeof(*model));
}

int axion_field_value(const AxionModel *model, const double X[4], double *value)
{
    double r, theta, t, phi, s, rp, rm, window, radial, angle, result;
    int k, status;
    if (!model || !value) return AXION_ERR_ARGUMENT;
    /* Baseline must not even dereference X or inspect table/spin/coords. */
    if (model->norm == 0) { *value = 0; return AXION_OK; }
    status = model_parameters(model);
    if (status) return status;
    if (!X || !model->table) return AXION_ERR_ARGUMENT;
    for (k = 0; k < 4; ++k)
        if (!isfinite(X[k])) return AXION_ERR_ARGUMENT;
    t = X[0]; theta = X[2]; phi = X[3];
    r = model->coords == AXION_MKSBHAC ? exp(X[1]) : X[1];
    if (!isfinite(r)) return AXION_ERR_NUMERIC;
    s = sqrt((1 - model->spin) * (1 + model->spin));
    rp = 1 + s; rm = 1 - s;
    if (!(r > rp)) return AXION_ERR_DOMAIN;
    /* Historical user-selected W25,15. Only exact underflow bypasses lookup. */
    window = exp(-pow(r / 25.0, 15.0));
    if (window == 0) { *value = 0; return AXION_OK; }
    if (model->coords == AXION_MKSBHAC)
        theta += 0.5 * model->hslope * sin(2 * X[2]);
    if (model->coords != AXION_BL) {
        /* M=1 ingoing KS -> BL, with the same integration constants as
         * GRmath.c. Avoid forming delta=r*r-2*r+a*a or a near-unity ratio. */
        double log_ratio = r - rp < 2 * s
            ? log(r - rp) - log(r - rm)
            : log1p(-2 * s / (r - rm));
        double time_shift = log(r - rp) + log(r - rm) + log_ratio / s;
        double phi_shift = model->spin * log_ratio / (2 * s);
        t -= time_shift;
        phi -= phi_shift;
    }
    if (!isfinite(theta) || !isfinite(t) || !isfinite(phi))
        return AXION_ERR_NUMERIC;
    status = axion_table_interpolate(model->table, model->alpha, r, &radial);
    if (status) return status;
    angle = -phi + (1022.0 + t) * model->omega + model->phase * AXION_PI / 8.0;
    if (!isfinite(angle)) return AXION_ERR_NUMERIC;
    result = model->norm / (2.0 * AXION_PI) * sin(theta) * radial *
             cos(angle) * window;
    if (!isfinite(result)) return AXION_ERR_NUMERIC;
    *value = result;
    return AXION_OK;
}
