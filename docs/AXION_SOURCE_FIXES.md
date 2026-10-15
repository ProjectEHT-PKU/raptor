# Axion source fixes: scope, conventions and validation

## Baseline and boundaries

This independent clone retains the ProjectEHT-PKU/raptor history. Baseline is
`b6df33c69e650f3617845f4324ba642bdb878582`, branch `fix/axion-source`.
The original mounted engineering, its run repository, scientific inputs and
outputs are unchanged. No saved RAPTOR binary was executed. No scientific
image was recomputed, and nothing was pushed.

The authorized scope is the fork-added axion defects numbered 1–7 in the
review: ray endpoints, the fork's OR fluid gate, thread-unsafe radial loading,
lookup bounds, zero-normalization isolation, alpha/frequency separation and
BL/KS field consistency. This is not an upstream RAPTOR upgrade.

The post-plasma Lie split, the `exp(-(r/25)^15)` window and the `1022` phase
reference are retained as the user's choices. They are not silently removed
or promoted to high-order/fully converged physics. Normalization and all
existing radial CSV bytes remain unchanged. Electron-temperature closure,
SIGMA_CUT, thermal emissivity/Faraday fits and the original stiffness/clip
algorithms are not changed. The original model/ and tetrad/metric/geodesic
source files remain byte-identical to the fork baseline.

## Changes

### 1. Forward endpoints and sign

The stored ray has records `[0,steps)`. Radiation propagates from `p=steps-1`
toward the observer at record zero. Segment endpoints are **p and p-1**, with
step length from p-1. No access to the unstored record `lightpath[steps]` remains.
Empty/single-record rays return zero without dereferencing missing records.

Use the explicit observer-minus-emitter endpoint convention:

```
A = g*a = N/(2*pi) * sin(theta) * R(alpha,r)
          * cos(-phi_BL + (1022+t_BL)*omega + phase*pi/8) * W25,15(r)
delta_chi = A(forward) - A(current)
P_out = exp(+2i*delta_chi) P_plasma
```

This agrees with the draft's explicit TotalDerivative equation. The draft's
transfer matrix combined with its printed `rho'_V=rho_V-2g da/ds` yields the
opposite sign; those draft expressions are internally inconsistent. We do not
claim the code agrees with every printed draft equation. Endpoint and rotation
sign were corrected together, rather than mechanically replacing +1 with -1.

### 2. Plasma and axion propagation are separate

The AND→OR change was introduced by this fork's axion patch, not inherited
from jordydavelaar RAPTOR. Restore valid-fluid AND plasma-domain gating, with
finite-positive checks for the current synchrotron coefficient/tetrad domain.
Zero-density and zero-field states are treated as no plasma interaction;
this is not a general statement about arbitrary other radiation mechanisms.

The original plasma solver no longer applies a second axis rotation. After
plasma, a deterministic local screen rotates the polarization vector once per
forward segment, even in vacuum, sigma-cut or pitch-aligned regions. It uses
a timelike normal observer and the same screen handedness as the existing
observer tetrad; it needs no stale plasma velocity/magnetic field. For a valid
null/transverse state it preserves I, V and |P| and is gauge-independent under
`f -> f+c*k`. No RNG or lazy global state is used.

Fluid lookup receives a coordinate copy, so its existing wrapping/offset does
not change the endpoints used by the field and GR transport. The geometric
wavevector/step parameterization is restored after the plasma's reciprocal
unit rescaling; this is equivalent affine parameterization, not a claim that
upstream vacuum transport units were wrong. The original plasma coefficient
formulas and GR transport functions are unchanged.

**The disabled baseline may differ from the old fork**, because the bad OR
plasma gate was also active at norm=0. Tests comparing a disabled model with
no model establish equality in the fixed chain, not historical image equality.

### 3. Immutable dynamic radial tables

`src/axion.c` loads and validates an entire headerless `alpha,r,R` CSV on the
startup thread, after the GRMHD model has supplied spin/coordinates and before
any refinement/OpenMP ray work. Queries are read-only. The loader accepts both
current 23x2400 and 15x2400 tables without hardcoded dimension assumptions;
rectangular shape, finite/nonnegative R, ordered axes and shared radial grid
are checked. Missing/short/malformed tables and out-of-range queries fail
explicitly. No renormalization or extrapolation is performed.

`norm=0` does not open a radial file, inspect its coordinates or query a table.
The module's error strings are immutable and indexed by status, not shared
thread-writable diagnostic state. Model free/reinitialization must happen
outside query regions.

The historical window is still evaluated. Only when it is exactly zero in
floating point can a field query bypass a radial table outside its radial
support; no new approximate small-field threshold is introduced. Other
out-of-domain requests fail. Supporting a wider cloud with a different window
requires an explicitly supported table/domain, not endpoint extrapolation.

### 4. Alpha and omega are independent explicit inputs

`AXION_ALPHA` selects the radial slice. `AXION_OMEGA` is the real oscillation
frequency Re(M omega). No n=2/n=4/table frequency formula is silently chosen.
All axion inputs are explicitly required and validated; omega-only legacy
files are rejected before model loading. The future run redesign must supply
the agreed frequency and record the selected table's spin/provenance.

The core-owned `run/model.in` is repaired to be a disabled example; no actual
production run repository is modified. The core-owned makefile adds axion.c
and header dependency files, so users can actually build the new source.

### 5. BL field represented in KS coordinates

The separated radial envelope and phase are defined in BL. KS/MKSBHAC inputs
are transformed with background metric spin a and M=1:

```
s = sqrt(1-a*a), r_plus=1+s, r_minus=1-s
F = log((r-r_plus)*(r-r_minus))
    + log((r-r_plus)/(r-r_minus))/s
G = a*log((r-r_plus)/(r-r_minus))/(2*s)
t_BL = t_KS-F
phi_BL = phi_KS-G
```

The implementation evaluates stable log differences/log1p rather than forming
a cancellation-prone delta/near-unity ratio. MKSBHAC uses r=exp(x1) and
`theta=x2+hslope*sin(2*x2)/2`, retaining the true angular mapping. BL, linear KS
and MKSBHAC(R0=0) are supported; other coordinate systems and extremal spins
are rejected for enabled axions rather than silently misinterpreted. Disabled
baselines do not load an axion model even on other coordinate systems.

This represents the **chosen real radial-envelope approximation in BL**. It
is not a claim to reconstruct the missing complex radial eigenfunction or its
phase, nor a change to radial eigenstate/normalization science.

## Validation performed

Reproduce the standalone suite from the repository root:

```bash
OMP_NUM_THREADS=8 bash tests/run_tests.sh
```

Optionally append existing radial CSV paths for read-only load/shape tests.
No GRMHD dump is accessed by these tests.

| Test | Result |
|---|---|
| Strict standalone C99 module, -Wall/-Wextra/-Wpedantic/-Werror | passed |
| Dynamic CSV, first/last nodes, single-alpha, malformed/short/nonrectangular/nonfinite input, domain errors | passed |
| 100000 immutable parallel queries, 8 OpenMP threads | passed |
| norm0 missing-table isolation, independent alpha/omega, 16 phases, fixed window and t_ref | passed |
| BL/KS/MKSBHAC same-event field covariance, background spin and near-horizon domain | passed |
| Existing reference CSV read-only load | 23x2400 and 15x2400 passed |
| Named input/startup, BHAC/BL/KS | 128 checks each passed |
| POL0 and unsupported metric input/startup | 127 checks each passed |
| Actual curved MKSBHAC metric/tetrad screen | 396 spin/slope/radius/theta/ray-direction cases passed |
| Curved screen spin-2 sign, circular component, norm, reverse and gauge invariance | passed |
| Actual transfer/tetrad/axion chain with synthetic flat geometry and fluid/coeff stubs | passed |
| steps0/1/exact max_steps, invalid fluid with untouched/garbage state, mutated lookup coordinates, sigma/ne0/pitch gates, endpoint telescoping | passed |
| ASan/UBSan/leak checks for module/input/curved screen/transfer | passed |
| Full BHAC source compilation and HDF5/GSL/OpenMP linking in a private /tmp build | passed; binary not scientifically run |

Logs: `tests/validation/`. Full-build compiler warnings inherited from the
baseline (geodesic theta history, photon declaration bounds, unchecked model
IO, unused model rc) remain; a successful build is not a declaration that the
unmodified upstream portions are warning-free or scientifically validated.

The container's Nix HDF5 h5cc wrapper contains a nonexistent compiler path.
Full build used the normal GCC with the Nix HDF5/GSL environment and explicit
`-lhdf5`; no system file or wrapper was edited. The exact invocation is in
`tests/validation/BUILD.md`.

An independent GPT-6.1 Sol read-only review found no remaining confirmed
blocking axion defect after checking the sign, screen handedness/gauge,
endpoints, coordinates and fork ownership of the OR change. This and the tests
are implementation checks, not real-data accuracy estimates.

## Not performed / future run redesign

- No real snapshot or image rerun; no estimate of historical image/limit changes.
- No replacement/relabeling of old CSVs, cubes or caches.
- No update to upstream solvers, fits, EOS, floors or diagnostics.
- No decision to remove the window, change 1022, change splitting order, choose
  a frequency formula, change physical amplitude convention or grid/spin model.
- No redesign of ProjectEHT-PKU/EHT-axion-pixel-by-pixel or its submodule pin.

Future runners must update the core pin and build sources, provide all new
inputs and isolate outputs. Old wrappers cannot be used unchanged; deliberate
input rejection prevents them from silently running the old mixed-alpha model.
