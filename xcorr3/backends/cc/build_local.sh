#!/usr/bin/env bash
# Build helper for "system" CUDA installs where nvcc lives in /usr/bin and there
# is no /usr/local/cuda (e.g. Ubuntu's `nvidia-cuda-toolkit` package). It points
# the Makefile at the system CUDA headers/libs; the GPU architecture is
# auto-detected by the Makefile unless CUDA_ARCH is exported.
#
# For a standard /usr/local/cuda install you do NOT need this -- just run `make`.
# For CUDA <= 11.5 + gcc >= 11 also run scripts/make_gcc11_compat.sh first
# (see README, "老工具链 / 系统版 CUDA").
#
# Override any of CUDA_NVCC / CUDA_INC / CUDA_LIB / CUDA_ARCH via the environment,
# e.g.  CUDA_ARCH=sm_86 ./build_local.sh
set -euo pipefail
cd "$(dirname "$0")"

: "${CUDA_NVCC:=$(command -v nvcc || echo /usr/bin/nvcc)}"
: "${CUDA_INC:=/usr/include}"
: "${CUDA_LIB:=/usr/lib/x86_64-linux-gnu}"

ARGS=(CUDA_NVCC="$CUDA_NVCC" CUDA_INC="$CUDA_INC" CUDA_LIB="$CUDA_LIB" CUDA_RPATH="$CUDA_LIB")
[ -n "${CUDA_ARCH:-}" ] && ARGS+=(CUDA_ARCH="$CUDA_ARCH")   # else Makefile auto-detects

exec make "${ARGS[@]}" "$@"
