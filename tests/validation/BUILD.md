# Full-source build check (not a scientific execution)

Base fork: b6df33c69e650f3617845f4324ba642bdb878582.
The source changes are on `fix/axion-source`. Run this from an independent
core checkout; create a **new** temporary directory, never use a production
run or replace saved binaries.

The tested build copied only the core's BHAC model files and makefile into a
private temporary directory. It did not use setup defaults or modify any
model header. Sources came from the independent core checkout:

```bash
root=$(pwd -P)
build=$(mktemp -d /tmp/raptor-fixed-build-XXXXXX)
cp "$root"/model/bhac/* "$build/"
cp "$root/run/makefile" "$build/Makefile"
export RAPTOR="$root"
cd "$build"
nix-shell -p hdf5 gsl gcc gnumake --run \
  'make -j2 CC="gcc -I$RAPTOR/src -I$PWD" \
      LDFLAGS="-fopenmp -lm -lgsl -lgslcblas -lhdf5"'
```

Dependencies: GCC 15.2.0, HDF5 1.14.6, GSL 2.8, GNU Make 4.4.1, Linux x86_64.
The Nix h5cc wrapper in this container references an absent compiler store
path; the explicit GCC link is an environment workaround, not a source or
system-wrapper modification. In a normal HDF5 toolchain the supplied makefile
continues to use h5cc.

`full-build.log` records the complete successful compilation/link, including
inherited compiler warnings. `incremental-build.log` records recompilation and
linking after the final polarized-transfer error-message refinement. No new
binary was run on scientific data. Build output lives only under `/tmp`.

`suite.log` records standalone tests using synthetic fixtures, the actual
production transport/geometry functions, ASan/UBSan and read-only radial
reference CSV loading. The reference paths are local audit provenance, not
required default test dependencies.
