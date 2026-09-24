#!/usr/bin/env bash
# Generate ./cuda11_gcc11_fix/bits/std_function.h to work around a bug where
# CUDA <= 11.5 nvcc (cudafe++) fails to parse gcc-11's libstdc++ <functional>:
#     bits/std_function.h: error: parameter packs not expanded with '...'
#
# It copies your own system header and strips the SFINAE constraints from the
# std::function(_Functor&&) constructor / operator= signatures that cudafe
# mis-parses (std::function is only pulled in transitively and is not used in
# xcorr_cc's compute path, so relaxing the overload constraints is harmless).
# The Makefile then auto-prepends this directory (see NVCC_EXTRA_INC).
#
# ONLY needed for CUDA <= 11.5 combined with gcc >= 11. On CUDA >= 11.6 this is
# unnecessary -- do not run it, and delete cuda11_gcc11_fix/ if it exists.
# The generated header is a machine-specific copy of your libstdc++ header and
# is intentionally git-ignored (not committed).
set -euo pipefail
cd "$(dirname "$0")/.."

SYS=$(echo '#include <functional>' | "${CXX:-g++}" -E -x c++ - 2>/dev/null \
        | grep -m1 -oE '/[^"]*bits/std_function\.h' || true)
[ -z "${SYS:-}" ] && SYS=$(ls /usr/include/c++/*/bits/std_function.h 2>/dev/null | sort -V | tail -1 || true)
[ -z "${SYS:-}" ] && { echo "ERROR: could not locate the system bits/std_function.h" >&2; exit 1; }
echo "source header: $SYS"

mkdir -p cuda11_gcc11_fix/bits
python3 - "$SYS" cuda11_gcc11_fix/bits/std_function.h <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
lines = open(src).read().split('\n')
out = []
rc = rn = rr = 0
for ln in lines:
    st = ln.strip()
    if st == 'typename _Constraints = _Requires<_Callable<_Functor>>>':
        for j in range(len(out) - 1, -1, -1):
            if out[j].strip() == 'template<typename _Functor,':
                out[j] = out[j].replace('template<typename _Functor,', 'template<typename _Functor>')
                break
        rc += 1
        continue
    if st.startswith('noexcept(') and '_S_nothrow_init' in st:
        rn += 1
        continue
    if st == '_Requires<_Callable<_Functor>, function&>':
        out.append(ln.replace('_Requires<_Callable<_Functor>, function&>', 'function&'))
        rr += 1
        continue
    out.append(ln)
open(dst, 'w').write('\n'.join(out))
print(f"patched: constraint={rc} noexcept={rn} return={rr}")
if rc == 0 and rr == 0:
    sys.stderr.write("WARNING: expected patterns not found; your gcc version may differ. "
                     "If the build still fails, patch bits/std_function.h by hand.\n")
PY

echo "generated cuda11_gcc11_fix/bits/std_function.h -- now just run 'make' (or ./build_local.sh)"
