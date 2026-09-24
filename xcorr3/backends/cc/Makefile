CC ?= gcc
CXX ?= g++
PKG_CONFIG ?= pkg-config

CUDA_HOME ?= /usr/local/cuda
CUDA_NVCC ?= $(CUDA_HOME)/bin/nvcc
CUDA_INC ?= $(CUDA_HOME)/include
CUDA_LIB ?= $(CUDA_HOME)/lib64
CUDA_RPATH ?= $(CUDA_LIB)
CUFFT_LIB ?= -lcufft
CUDART_LIB ?= -lcudart

CUDA_BIN_DIR := $(patsubst %/,%, $(dir $(CUDA_NVCC)))
CUDA_NVVM_BIN ?= $(abspath $(CUDA_BIN_DIR)/../nvvm/bin)
NVCC_RUN = PATH=$(CUDA_BIN_DIR):$(CUDA_NVVM_BIN):$$PATH $(CUDA_NVCC)

GLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags glib-2.0)
GLIB_LIBS := $(shell $(PKG_CONFIG) --libs glib-2.0)

# GPU architecture: auto-detect from first visible GPU, fall back to sm_70
# Override with: make CUDA_ARCH=sm_89
CUDA_ARCH ?= $(shell $(CUDA_NVCC) detect_arch.cu -o /dev/null --run 2>/dev/null || \
               nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | \
               head -1 | sed 's/\.//' | sed 's/^/sm_/' || echo sm_70)

# Local workaround for CUDA <= 11.5 nvcc (cudafe++) failing to parse gcc-11's
# libstdc++ <functional> ("std_function.h: parameter packs not expanded").
# If ./cuda11_gcc11_fix/bits/std_function.h exists, prepend that dir so the
# patched header is found before the system one. Harmless on newer toolchains.
GCC11_FIX_DIR := $(wildcard $(CURDIR)/cuda11_gcc11_fix)
NVCC_EXTRA_INC ?= $(if $(GCC11_FIX_DIR),-I$(GCC11_FIX_DIR),)

CFLAGS += -O3 -std=gnu99 -I. $(GLIB_CFLAGS)
NVCCFLAGS += -O3 -std=c++14 -arch=$(CUDA_ARCH) $(NVCC_EXTRA_INC) -I. -I$(CUDA_INC) $(GLIB_CFLAGS)
LDLIBS_CUDA += $(GLIB_LIBS) -L$(CUDA_LIB) $(CUFFT_LIB) $(CUDART_LIB) -Xlinker -rpath -Xlinker $(CUDA_RPATH)

all: xcorr_cc

prm_helper.o: prm_helper.c xcorr2.h
	$(CC) -c -o $@ $< $(CFLAGS)

xcorr2_args.o: xcorr2_args.c xcorr2.h xcorr2_args.h
	$(CC) -c -o $@ $< $(CFLAGS)

xcorr_postprocess.o: xcorr_postprocess.cpp xcorr_postprocess.h xcorr_io.h xcorr2.h xcorr2_args.h
	$(CXX) -c -o $@ $< -O3 -std=c++14 -I. $(GLIB_CFLAGS) $(CXXFLAGS)

xcorr_cc.o: xcorr_cc.cu xcorr2.h xcorr2_args.h xcorr_io.h xcorr_postprocess.h
	$(NVCC_RUN) -c -o $@ $< $(NVCCFLAGS)

xcorr_cc: xcorr_cc.o prm_helper.o xcorr2_args.o xcorr_postprocess.o
	$(NVCC_RUN) -o $@ $^ $(NVCCFLAGS) $(LDLIBS_CUDA)

.PHONY: all clean print-config

print-config:
	@echo "CUDA_NVCC  = $(CUDA_NVCC)"
	@echo "CUDA_INC   = $(CUDA_INC)"
	@echo "CUDA_LIB   = $(CUDA_LIB)"
	@echo "CUDA_ARCH  = $(CUDA_ARCH)"
	@echo "NVCC_EXTRA_INC = $(NVCC_EXTRA_INC)"

clean:
	rm -f *.o xcorr_cc
