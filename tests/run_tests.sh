#!/usr/bin/env bash
# Synthetic regression suite. Never executes a saved RAPTOR or opens GRMHD data.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build=$(mktemp -d /tmp/raptor-axion-suite.XXXXXX)
trap 'rm -rf -- "$build"' EXIT
cd "$root"
cc=${CC:-gcc}
common=(-std=c99 -O1 -g -fopenmp -fno-omit-frame-pointer)
strict=(-Wall -Wextra -Wpedantic -Werror)
includes=(-Isrc -Imodel/bhac -Itests/stubs)

"$cc" "${common[@]}" "${strict[@]}" -Isrc src/axion.c tests/test_axion.c -lm -o "$build/axion"
mkdir "$build/fixtures"
"$build/axion" "$build/fixtures" "$@"
for variant in bhac bl ks pol0 unsupported; do
    case "$variant" in
        bhac) extra=();;
        bl) extra=(-DTEST_METRIC=1);;
        ks) extra=(-DTEST_METRIC=3);;
        pol0) extra=(-DTEST_POL=0);;
        unsupported) extra=(-DTEST_METRIC=8);;
    esac
    "$cc" "${common[@]}" -Wall -Wextra -Werror "${includes[@]}" \
        "${extra[@]}" tests/test_input.c src/axion.c -lm -o "$build/input-$variant"
    "$build/input-$variant"
done
"$cc" "${common[@]}" "${strict[@]}" -fsanitize=address,undefined -fno-pie -no-pie \
    -Isrc src/axion.c tests/test_axion.c -lm -o "$build/axion-asan"
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$build/axion-asan" "$build/fixtures" "$@"
"$cc" "${common[@]}" -Wall -Wextra -Werror "${includes[@]}" \
    -fsanitize=address,undefined -fno-pie -no-pie tests/test_input.c src/axion.c -lm -o "$build/input-asan"
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$build/input-asan"

# Actual curved MKSBHAC geometry; section GC discards unrelated fluid/HDF code.
geometry=(tests/test_screen_curved.c src/pol_rte_integrator.c src/tetrad.c src/GRmath.c src/metric.c)
"$cc" "${common[@]}" "${includes[@]}" -ffunction-sections -fdata-sections \
    "${geometry[@]}" -Wl,--gc-sections -lm -o "$build/curved"
"$build/curved"
"$cc" "${common[@]}" "${includes[@]}" -ffunction-sections -fdata-sections \
    -fsanitize=address,undefined -fno-pie -no-pie \
    "${geometry[@]}" -Wl,--gc-sections -lm -o "$build/curved-asan"
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$build/curved-asan"
bash tests/run_transfer_tests.sh
printf '%s\n' 'PASS all axion source regressions'
