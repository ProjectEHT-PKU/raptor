# Explicit axion input and startup contract

## Parameter file

`RAPTOR model.in GRMHD_FILE TIME_INIT` requires exactly three arguments after
its executable name. `TIME_INIT` must be finite and the GRMHD path at most 255
bytes. No GRMHD simulation is needed for the isolated input tests below.

The reader accepts whitespace-delimited `LABEL UNIT VALUE` triples. Newlines
are optional, ordering is arbitrary, and every field below must occur exactly
once, with the exact case-sensitive label and unit. Unknown/duplicate fields,
missing values, extra tokens, invalid units, nonfinite numbers, trailing numeric
text and tokens over 255 bytes are errors. Comments are not parameter triples
and are not accepted. Integer fields require integer spelling, not `1.0`.

| Label | Unit | Allowed value |
| --- | --- | --- |
| `MBH` | `(Msun)` | positive, finite; cgs conversion must remain finite |
| `DISTANCE` | `(kpc)` | positive, finite; cgs conversion must remain finite |
| `M_UNIT` | `(g)` | positive, finite |
| `R_LOW`, `R_HIGH` | `(-)` | positive, finite; no new ordering assumption |
| `INCLINATION` | `(deg)` | finite, 0 through 180 inclusive |
| `IMG_WIDTH`, `IMG_HEIGHT` | `(pixels)` | equal positive integers, multiples of compiled `num_pixels_1d` (normally 10); block count fits `int` |
| `CAM_SIZE_X`, `CAM_SIZE_Y` | `(Rg)` | positive, finite |
| `FREQS_PER_DEC` | `(-)` | positive, finite |
| `FREQ_MIN` | `(Hz)` | positive, finite |
| `STEPSIZE` | `(-)` | positive, finite |
| `AXION_ALPHA` | `(-)` | independent radial-table mass parameter; positive when enabled |
| `AXION_OMEGA` | `(-)` | independent true real frequency `Re(M omega)`; positive when enabled |
| `AXION_NORM` | `(-)` | finite nonnegative reference field strength `N` |
| `AXION_PHASE` | `(-)` | integer phase index 0 through 15; phase is index times pi/8 |
| `AXION_TABLE` | `(-)` | a single filename/path token, at most 255 bytes, no spaces |
| `MAX_LEVEL` | `(-)` | positive integer |

All five `AXION_*` fields are required even in an unpolarized build and even
when `AXION_NORM=0`. For that disabled baseline, alpha and omega may be zero
(but not negative), and the table need not exist. The core example
`run/model.in` is such a baseline. All existing non-axion example values are
unchanged. `POL=0` with nonzero norm is rejected, rather than silently ignoring
rotation.

Relative table paths are relative to the **process working directory**, not
the input file's directory. The parser records but does not open the table.

## Migration: alpha is not omega

Old inputs supplied one `AXION_OMEGA` value for both the radial mass parameter
and temporal frequency. They are deliberately rejected. A complete old file,
including an old script's concatenated/no-newline form, reports missing
`AXION_ALPHA` with a migration explanation. Do not simply rename the old field
or copy its number into both new fields.

Supply `AXION_ALPHA` for the chosen radial table and compute an independent
`AXION_OMEGA` from the **user-selected physical mode/formula**. The core neither
chooses an n=2/n=4 approximation nor computes omega from alpha. The old
`Run.py`/input generators need a future redesign; this change does not modify
any actual running repository, script pipeline, upstream solver, or simulation.

Example enabled triples (numbers are independent illustrative inputs, **not**
a prescribed frequency formula; the table must cover alpha=0.2):

```
AXION_ALPHA (-) 0.2
AXION_OMEGA (-) 0.197
AXION_NORM  (-) 1
AXION_PHASE (-) 3
AXION_TABLE (-) axion_radial.csv
```

## Lifetime and coordinates

`AxionModel axion_model` is global. Startup calls `init_model()` and
`set_constants()` first, so model globals `a` and `hslope` are available. It
then calls exactly once, on the startup thread:

```c
axion_model_init(&axion_model, axion_table, axion_alpha, axion_omega,
                 axion_norm, phase_index, a, hslope, axion_coords);
```

`AXION_OK` (zero) means success. A failure is reported with `axion_error()` and
aborts before camera/refinement/ray work. Enabled configurations support only
`metric==BL`, `KS`, or `MKSBHAC`, mapped respectively to `AXION_BL`, `AXION_KS`,
or `AXION_MKSBHAC`. Other metrics fail before tracing; disabled baselines work
without loading a table even on unsupported metrics. The table is headerless
`alpha,r,R` CSV and uses the validation/domain contract in `src/axion.h`.
An enabled model also validates spin, slope, alpha coverage and table shape.

The radial profile and oscillatory phase are defined in **Boyer–Lindquist**
coordinates: linear radius r, theta, `t_BL`, `phi_BL`. KS/MKSBHAC positions
must be converted by the axion model; do not use KS time or phi directly in a
BL phase. The MKSBHAC mapping assumes the model's supported R0=0 convention.
The table's physical provenance/spin consistency is the caller's responsibility.

The historical time reference **1022**, window `exp(-(r/25)^15)`, phase-index
convention, and post-plasma Lie splitting remain unchanged. Endpoints and their
rotation sign are corrected together: `delta_chi=A(p-1)-A(p)` and
`P -> exp(2i*delta_chi) P`, following the draft's explicit observer-minus-emitter
endpoint convention. This input/startup integration makes no upstream solver,
EOS, SIGMA, or plasma-coefficient changes; a separate deterministic screen
rotation is added for axion propagation without a fluid tetrad. Polarized transfer now receives
`const AxionModel *axion, int block, int pixel` instead of three legacy scalars.
Every pixel gets `&axion_model`. Queries are immutable, table loading never
occurs inside ray/OpenMP work, and `axion_model_free()` runs after the image is
finished. No lazy per-ray initialization is permitted.

## Build and isolated tests

The core-owned `run/makefile` adds `axion.c`, emits `-MMD -MP` dependencies,
and uses order-only build-directory creation for parallel make. Header edits
therefore rebuild dependent objects. It still expects the normal configured
model headers/source and HDF5/GSL toolchain; it does not redesign a run project.

`tests/test_input.c` is a standalone synthetic parser/startup regression harness.
Its compile/run instructions are at the top of that file. It uses the real
core/main and axion implementation but stubs camera, GRMHD and transfer work;
it is not a scientific run. It checks full and newline-free files, old/missing/
malformed inputs, exact labels/units/ranges, argv handling, disabled missing-table
behavior, and single startup initialization before any parallel query. Compile
variants exercise POL=0 and unsupported coordinates as well.
