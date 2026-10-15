#!/usr/bin/env bash
# Real production transfer/tetrad/axion, no HDF5 or simulation-data dependency.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build=$(mktemp -d /tmp/raptor-transfer.XXXXXX)
trap 'rm -rf -- "$build"' EXIT
cc=${CC:-gcc}

# Keep the BHAC config/model headers unmodified; all overrides/dependencies
# live in this private build directory or tests/stubs. In particular, do not
# compile against a stale run/definitions.h or any old RAPTOR binary.
cp "$root/model/bhac/definitions.h" "$root"/model/bhac/model_*.h "$build/"
printf '0.2,3,1\n0.2,9,2\n' > "$build/radial.csv"
flags=(-std=c99 -fopenmp -O1 -g -ffunction-sections -fdata-sections
       -fno-omit-frame-pointer -I"$build" -I"$root/tests/stubs" -I"$root/src")
link=(-Wl,--gc-sections -Wl,--wrap=create_observer_tetrad
      -Wl,--wrap=axion_field_value -lm)
sources=("$root/tests/test_transfer.c" "$root/src/pol_rte_integrator.c"
         "$root/src/tetrad.c" "$root/src/axion.c")

printf '%s\n' 'Building and running production-chain transfer tests (C99/OpenMP)'
"$cc" "${flags[@]}" "${sources[@]}" "${link[@]}" -o "$build/transfer"
"$build/transfer" "$build/radial.csv"
printf '%s\n' 'Building and running production-chain transfer tests (ASan/UBSan)'
# Non-PIE avoids intermittent ASan shadow-map startup conflicts in containers.
"$cc" "${flags[@]}" -fsanitize=address,undefined -fno-pie -no-pie \
    "${sources[@]}" "${link[@]}" -o "$build/transfer-asan"
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
    "$build/transfer-asan" "$build/radial.csv"
