#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <thread>
#include <csignal>
#include "xcorr_io.h"
#include "xcorr_postprocess.h"

#include <sys/statvfs.h>
#include <unistd.h>

extern "C" {
#include "xcorr2_args.h"
}

#include <cuda_runtime.h>
#include <cufft.h>
#include <cub/cub.cuh>

#define CUDA_CHECK(call) do { \
    cudaError_t err__ = (call); \
    if (err__ != cudaSuccess) { \
        throw std::runtime_error(std::string(#call) + ": " + cudaGetErrorString(err__)); \
    } \
} while (0)

static const char *cufft_error_string(cufftResult code) {
    switch (code) {
        case CUFFT_SUCCESS: return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
        default: return "CUFFT_UNKNOWN_ERROR";
    }
}

#define CUFFT_CHECK(call) do { \
    cufftResult err__ = (call); \
    if (err__ != CUFFT_SUCCESS) { \
        throw std::runtime_error(std::string(#call) + ": " + cufft_error_string(err__)); \
    } \
} while (0)

/* Software fallback for double-precision atomicAdd on pre-Pascal GPUs
 * (compute capability < 6.0), where the hardware intrinsic is unavailable.
 * On sm_60+ the built-in overload is used and this block is skipped, so it is
 * zero-cost on modern GPUs. Verbatim from the CUDA C Programming Guide. */
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__ double atomicAdd(double *address, double val) {
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

int div_up(int value, int block) {
    return value / block + (value % block != 0);
}

uint64_t checked_mul_u64(uint64_t a, uint64_t b, const char *label) {
    if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) {
        throw std::runtime_error(std::string("Integer overflow while computing ") + label);
    }
    return a * b;
}

uint64_t checked_add_u64(uint64_t a, uint64_t b, const char *label) {
    if (a > std::numeric_limits<uint64_t>::max() - b) {
        throw std::runtime_error(std::string("Integer overflow while computing ") + label);
    }
    return a + b;
}

bool is_power_of_two_int(int value) {
    return value > 0 && (value & (value - 1)) == 0;
}

std::string format_bytes(uint64_t bytes) {
    static const char *units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1024.0 && unit < 4) {
        value /= 1024.0;
        ++unit;
    }

    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f %s", value, units[unit]);
    return std::string(buffer);
}

template <typename T>
class DeviceBuffer {
public:
    DeviceBuffer() : data_(nullptr), size_(0) {}

    explicit DeviceBuffer(size_t count) : data_(nullptr), size_(0) {
        allocate(count);
    }

    ~DeviceBuffer() {
        release();
    }

    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;

    DeviceBuffer(DeviceBuffer &&other) noexcept : data_(other.data_), size_(other.size_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }

    DeviceBuffer &operator=(DeviceBuffer &&other) noexcept {
        if (this != &other) {
            release();
            data_ = other.data_;
            size_ = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    void allocate(size_t count) {
        if (count == size_ && data_ != nullptr) {
            return;
        }
        release();
        size_ = count;
        if (size_ > 0) {
            CUDA_CHECK(cudaMalloc(&data_, size_ * sizeof(T)));
        }
    }

    void zero() {
        if (size_ > 0) {
            CUDA_CHECK(cudaMemset(data_, 0, size_ * sizeof(T)));
        }
    }

    void copy_from_host(const std::vector<T> &host) {
        allocate(host.size());
        if (!host.empty()) {
            CUDA_CHECK(cudaMemcpy(data_, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice));
        }
    }

    void copy_from_host(const T *host, size_t count, size_t offset = 0) {
        if (offset + count > size_) {
            throw std::runtime_error("DeviceBuffer::copy_from_host out of bounds");
        }
        if (count > 0) {
            CUDA_CHECK(cudaMemcpy(data_ + offset, host, count * sizeof(T), cudaMemcpyHostToDevice));
        }
    }

    void copy_to_host(std::vector<T> &host) const {
        host.resize(size_);
        if (size_ > 0) {
            CUDA_CHECK(cudaMemcpy(host.data(), data_, size_ * sizeof(T), cudaMemcpyDeviceToHost));
        }
    }

    void copy_to_host(T *host, size_t count, size_t offset = 0) const {
        if (offset + count > size_) {
            throw std::runtime_error("DeviceBuffer::copy_to_host out of bounds");
        }
        if (count > 0) {
            CUDA_CHECK(cudaMemcpy(host, data_ + offset, count * sizeof(T), cudaMemcpyDeviceToHost));
        }
    }

    void copy_from_device(const DeviceBuffer<T> &other) {
        allocate(other.size_);
        if (size_ > 0) {
            CUDA_CHECK(cudaMemcpy(data_, other.data_, size_ * sizeof(T), cudaMemcpyDeviceToDevice));
        }
    }

    T *get() { return data_; }
    const T *get() const { return data_; }
    size_t size() const { return size_; }

private:
    void release() {
        if (data_ != nullptr) {
            cudaFree(data_);
            data_ = nullptr;
            size_ = 0;
        }
    }

    T *data_;
    size_t size_;
};

template <typename T>
class PinnedHostBuffer {
public:
    PinnedHostBuffer() : data_(nullptr), size_(0) {}

    explicit PinnedHostBuffer(size_t count) : data_(nullptr), size_(0) {
        allocate(count);
    }

    ~PinnedHostBuffer() {
        release();
    }

    PinnedHostBuffer(const PinnedHostBuffer &) = delete;
    PinnedHostBuffer &operator=(const PinnedHostBuffer &) = delete;

    void allocate(size_t count) {
        if (count == size_ && data_ != nullptr) {
            return;
        }
        release();
        size_ = count;
        if (size_ > 0) {
            CUDA_CHECK(cudaMallocHost(&data_, size_ * sizeof(T)));
        }
    }

    T *data() { return data_; }
    const T *data() const { return data_; }
    size_t size() const { return size_; }

private:
    void release() {
        if (data_ != nullptr) {
            cudaFreeHost(data_);
            data_ = nullptr;
            size_ = 0;
        }
    }

    T *data_;
    size_t size_;
};

struct CufftPlan2D {
    cufftHandle handle;

    CufftPlan2D(int height, int width) : handle(0) {
        CUFFT_CHECK(cufftPlan2d(&handle, height, width, CUFFT_C2C));
    }

    ~CufftPlan2D() {
        if (handle != 0) {
            cufftDestroy(handle);
        }
    }

    CufftPlan2D(const CufftPlan2D &) = delete;
    CufftPlan2D &operator=(const CufftPlan2D &) = delete;
};

std::string format_duration(double seconds) {
    long total = static_cast<long>(seconds + 0.5);
    long hours = total / 3600;
    long minutes = (total % 3600) / 60;
    long secs = total % 60;
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%02ld:%02ld:%02ld", hours, minutes, secs);
    return std::string(buffer);
}

void render_progress(
        size_t processed,
        size_t total,
        const std::chrono::steady_clock::time_point &start_time,
        bool done) {
    using namespace std::chrono;

    auto now = steady_clock::now();
    double elapsed = duration_cast<duration<double>>(now - start_time).count();
    double ratio = total > 0 ? static_cast<double>(processed) / static_cast<double>(total) : 1.0;
    ratio = std::max(0.0, std::min(1.0, ratio));
    double eta = (processed > 0 && processed < total)
        ? elapsed * static_cast<double>(total - processed) / static_cast<double>(processed)
        : 0.0;

    const int bar_width = 32;
    int filled = static_cast<int>(ratio * bar_width + 0.5);
    if (filled > bar_width) {
        filled = bar_width;
    }

    char bar[bar_width + 1];
    for (int i = 0; i < bar_width; ++i) {
        bar[i] = (i < filled) ? '#' : '-';
    }
    bar[bar_width] = '\0';

    std::fprintf(stderr,
                 "\r[%s] %6.2f%% %zu/%zu elapsed %s eta %s",
                 bar,
                 ratio * 100.0,
                 processed,
                 total,
                 format_duration(elapsed).c_str(),
                 format_duration(eta).c_str());
    if (done) {
        std::fprintf(stderr, "\n");
    }
    std::fflush(stderr);
}

__device__ cufftComplex complex_mul(cufftComplex a, cufftComplex b) {
    return make_cuFloatComplex(
        a.x * b.x - a.y * b.y,
        a.x * b.y + a.y * b.x);
}

__device__ cufftComplex complex_conj(cufftComplex a) {
    return make_cuFloatComplex(a.x, -a.y);
}

__global__ void crop_complex_kernel(
        const cufftComplex *src,
        int src_width,
        int x0,
        int width,
        int height,
        cufftComplex *dst) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    dst[y * width + x] = src[y * src_width + x0 + x];
}

__global__ void complex_abs_kernel(const cufftComplex *src, float *dst, int count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) {
        return;
    }

    dst[idx] = hypotf(src[idx].x, src[idx].y);
}

/* Bug #4 fix: extract real part of complex array (matches original GMTSAR
 * highres_corr.c which uses .r for peak finding, not Cabs). */
__global__ void complex_real_kernel(const cufftComplex *src, float *dst, int count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) {
        return;
    }

    dst[idx] = src[idx].x;
}

__global__ void center_and_mask_kernel(
        float *c1r,
        float *c2r,
        int width,
        int height,
        int xsearch,
        int ysearch,
        float mean1,
        float mean2) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    int idx = y * width + x;
    c1r[idx] -= mean1;

    float value = c2r[idx] - mean2;
    if (y < ysearch || y >= height - ysearch || x < xsearch || x >= width - xsearch) {
        value = 0.0f;
    }
    c2r[idx] = value;
}

__global__ void real_to_complex_kernel(const float *src, cufftComplex *dst, int count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) {
        return;
    }

    dst[idx] = make_cuFloatComplex(src[idx], 0.0f);
}

__global__ void freq_product_kernel(
        const cufftComplex *a,
        const cufftComplex *b,
        int height,
        int width,
        cufftComplex *out) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    int idx = y * width + x;
    float sign = ((x + y) & 1) ? -1.0f : 1.0f;
    cufftComplex value = complex_mul(a[idx], complex_conj(b[idx]));
    out[idx] = make_cuFloatComplex(sign * value.x, sign * value.y);
}

/* Nyquist handling during FFT zero-padding for sub-pixel interpolation, chosen
 * at runtime via `split` (set from -nyquist_split; default 0):
 *   split==0: bit-compatible with GMTSAR fft_arrange_interpolate() -- the
 *     Nyquist bin is mapped only to the negative-frequency slot at full
 *     amplitude (no split). Reproduces the reference implementation.
 *   split==1: split the Nyquist bin symmetrically (0.5 to +N/2, 0.5 to -N/2),
 *     keeping the interpolated correlation surface real/symmetric. */
__global__ void pad_fft_kernel(
        const cufftComplex *in_fft,
        int height,
        int width,
        cufftComplex *out_fft,
        int out_height,
        int out_width,
        int split) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    cufftComplex value = in_fft[y * width + x];

    if (!split) {
        /* GMTSAR-compatible: match fft_arrange_interpolate() -- no Nyquist split.
         * Positive freqs keep their index; negative freqs (and the Nyquist bin,
         * treated as the first negative freq) shift to the top of the spectrum. */
        int oy = (y < height / 2) ? y : y + (out_height - height);
        int ox = (x < width / 2)  ? x : x + (out_width  - width);
        out_fft[oy * out_width + ox] = value;
        return;
    }

    /* Y Nyquist needs splitting only when we are actually expanding in Y */
    bool y_nyq = (y == height / 2) && (out_height > height);
    bool x_nyq = (x == width / 2)  && (out_width  > width);

    /* Determine base output y coordinate(s) */
    int oy0, oy1;  /* oy1 is only used when y_nyq */
    if (y_nyq) {
        oy0 = height / 2;                  /* positive-frequency Nyquist row */
        oy1 = out_height - height / 2;     /* negative-frequency Nyquist row */
    } else {
        oy0 = (y < height / 2) ? y : y + (out_height - height);
        oy1 = oy0;  /* unused but keep compiler happy */
    }

    /* Apply halving factors for Nyquist bins */
    float y_scale = y_nyq ? 0.5f : 1.0f;
    float x_scale = x_nyq ? 0.5f : 1.0f;
    float total_scale = y_scale * x_scale;
    cufftComplex sv = make_cuFloatComplex(value.x * total_scale,
                                          value.y * total_scale);

    /* Determine output x coordinate(s) */
    if (x_nyq) {
        int ox_lo = width / 2;
        int ox_hi = out_width - width / 2;
        out_fft[oy0 * out_width + ox_lo] = sv;
        out_fft[oy0 * out_width + ox_hi] = sv;
        if (y_nyq) {
            out_fft[oy1 * out_width + ox_lo] = sv;
            out_fft[oy1 * out_width + ox_hi] = sv;
        }
    } else if (x < width / 2) {
        out_fft[oy0 * out_width + x] = sv;
        if (y_nyq) {
            out_fft[oy1 * out_width + x] = sv;
        }
    } else {
        int mapped_x = x + out_width - width;
        out_fft[oy0 * out_width + mapped_x] = sv;
        if (y_nyq) {
            out_fft[oy1 * out_width + mapped_x] = sv;
        }
    }
}

__global__ void scale_complex_kernel(cufftComplex *data, float scale, int count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) {
        return;
    }

    data[idx].x *= scale;
    data[idx].y *= scale;
}

__global__ void scaled_abs_crop_kernel(
        const cufftComplex *src,
        int src_width,
        int x0,
        int y0,
        int width,
        int height,
        float scale,
        float *dst) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    cufftComplex value = src[(y0 + y) * src_width + (x0 + x)];
    dst[y * width + x] = hypotf(value.x, value.y) * scale;
}

__global__ void crop_scale_pow_quarter_kernel(
        const float *src,
        int src_width,
        int x0,
        int y0,
        int width,
        int height,
        float scale,
        float *dst) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) {
        return;
    }

    float value = src[(y0 + y) * src_width + (x0 + x)] * scale;
    dst[y * width + x] = powf(value, 0.25f);
}

__global__ void sum_two_arrays_kernel(
        const float *a,
        const float *b,
        int count,
        double *sums) {
    extern __shared__ double shared[];
    double *sum_a = shared;
    double *sum_b = shared + blockDim.x;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    int stride = blockDim.x * gridDim.x;

    double local_a = 0.0;
    double local_b = 0.0;
    for (int i = idx; i < count; i += stride) {
        local_a += static_cast<double>(a[i]);
        local_b += static_cast<double>(b[i]);
    }

    sum_a[tid] = local_a;
    sum_b[tid] = local_b;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            sum_a[tid] += sum_a[tid + offset];
            sum_b[tid] += sum_b[tid + offset];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(&sums[0], sum_a[0]);
        atomicAdd(&sums[1], sum_b[0]);
    }
}

/* Reduce one array to (sum, sum-of-squares) -> used for the correlation-surface
 * peak-significance SNR = (peak - mean)/std that flags spurious (gross) peaks. */
__global__ void sum_sumsq_kernel(
        const float *x,
        int count,
        double *out) {
    extern __shared__ double shared[];
    double *s_sum = shared;
    double *s_sq = shared + blockDim.x;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    int stride = blockDim.x * gridDim.x;

    double local_sum = 0.0;
    double local_sq = 0.0;
    for (int i = idx; i < count; i += stride) {
        double v = static_cast<double>(x[i]);
        local_sum += v;
        local_sq += v * v;
    }

    s_sum[tid] = local_sum;
    s_sq[tid] = local_sq;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            s_sum[tid] += s_sum[tid + offset];
            s_sq[tid] += s_sq[tid + offset];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(&out[0], s_sum[0]);
        atomicAdd(&out[1], s_sq[0]);
    }
}

__global__ void corr_metrics_kernel(
        const float *c1r,
        const float *c2r,
        int nx_win,
        int xsearch,
        int ysearch,
        int nx_corr,
        int ny_corr,
        int peak_x,
        int peak_y,
        double *accum) {
    extern __shared__ double shared[];
    double *num_shared = shared;
    double *denom1_shared = shared + blockDim.x;
    double *denom2_shared = shared + 2 * blockDim.x;

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    int total = nx_corr * ny_corr;
    int stride = blockDim.x * gridDim.x;

    double numerator = 0.0;
    double denom1 = 0.0;
    double denom2 = 0.0;

    for (int linear = idx; linear < total; linear += stride) {
        int y = linear / nx_corr;
        int x = linear % nx_corr;
        int row1 = (ysearch + peak_y + y) * nx_win;
        int row2 = (ysearch + y) * nx_win;
        float a = c1r[row1 + xsearch + peak_x + x];
        float b = c2r[row2 + xsearch + x];
        numerator += static_cast<double>(a) * static_cast<double>(b);
        denom1 += static_cast<double>(a) * static_cast<double>(a);
        denom2 += static_cast<double>(b) * static_cast<double>(b);
    }

    num_shared[tid] = numerator;
    denom1_shared[tid] = denom1;
    denom2_shared[tid] = denom2;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            num_shared[tid] += num_shared[tid + offset];
            denom1_shared[tid] += denom1_shared[tid + offset];
            denom2_shared[tid] += denom2_shared[tid + offset];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(&accum[0], num_shared[0]);
        atomicAdd(&accum[1], denom1_shared[0]);
        atomicAdd(&accum[2], denom2_shared[0]);
    }
}

void read_slc_rows_chunk(
        std::ifstream &fin,
        int start,
        int n_rows,
        int nx,
        int16_t *io_buffer,
        cufftComplex *out_rows) {
    std::streamoff offset = static_cast<std::streamoff>(nx) * start * sizeof(int16_t) * 2;
    size_t complex_count = static_cast<size_t>(n_rows) * nx;
    size_t int16_count = complex_count * 2;

    fin.clear();
    fin.seekg(offset, std::ios::beg);
    if (!fin) {
        throw std::runtime_error("Failed to seek SLC file");
    }

    fin.read(reinterpret_cast<char *>(io_buffer), static_cast<std::streamsize>(int16_count * sizeof(int16_t)));
    if (fin.gcount() != static_cast<std::streamsize>(int16_count * sizeof(int16_t))) {
        throw std::runtime_error("Failed to read data from SLC file");
    }

    for (size_t i = 0; i < complex_count; ++i) {
        out_rows[i] = make_cuFloatComplex(
            static_cast<float>(io_buffer[2 * i]),
            static_cast<float>(io_buffer[2 * i + 1]));
    }
}

class SlcRowCache {
public:
    SlcRowCache(std::ifstream &fin, int nx, int n_rows)
        : fin_(fin), nx_(nx), n_rows_(n_rows), current_start_(std::numeric_limits<int>::min()), active_(0) {
        size_t row_count = static_cast<size_t>(nx_) * n_rows_;
        buffers_[0].allocate(row_count);
        buffers_[1].allocate(row_count);
        io_buffer_.allocate(row_count * 2);
        host_rows_.allocate(row_count);
    }

    const cufftComplex *data() const {
        return buffers_[active_].get();
    }

    void ensure_window(int start_row) {
        if (current_start_ == std::numeric_limits<int>::min()) {
            full_reload(start_row);
            return;
        }

        int delta = start_row - current_start_;
        if (delta == 0) {
            return;
        }

        if (delta > 0 && delta < n_rows_) {
            slide_forward(delta);
            return;
        }

        full_reload(start_row);
    }

private:
    void full_reload(int start_row) {
        read_slc_rows_chunk(fin_, start_row, n_rows_, nx_, io_buffer_.data(), host_rows_.data());
        buffers_[active_].copy_from_host(host_rows_.data(), static_cast<size_t>(nx_) * n_rows_);
        current_start_ = start_row;
    }

    void slide_forward(int delta_rows) {
        int keep_rows = n_rows_ - delta_rows;
        size_t keep_count = static_cast<size_t>(keep_rows) * nx_;
        DeviceBuffer<cufftComplex> &src = buffers_[active_];
        DeviceBuffer<cufftComplex> &dst = buffers_[1 - active_];

        if (keep_count > 0) {
            CUDA_CHECK(cudaMemcpy(dst.get(),
                                  src.get() + static_cast<size_t>(delta_rows) * nx_,
                                  keep_count * sizeof(cufftComplex),
                                  cudaMemcpyDeviceToDevice));
        }

        read_slc_rows_chunk(fin_, current_start_ + n_rows_, delta_rows, nx_, io_buffer_.data(), host_rows_.data());
        dst.copy_from_host(host_rows_.data(), static_cast<size_t>(delta_rows) * nx_, keep_count);

        active_ = 1 - active_;
        current_start_ += delta_rows;
    }

    std::ifstream &fin_;
    int nx_;
    int n_rows_;
    int current_start_;
    int active_;
    DeviceBuffer<cufftComplex> buffers_[2];
    PinnedHostBuffer<int16_t> io_buffer_;
    PinnedHostBuffer<cufftComplex> host_rows_;
};

void initialize_cuda(const st_xcorr_args &args) {
    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA devices available");
    }

    CUDA_CHECK(cudaSetDevice(0));

    if (args.device == XCORR2_DEVICE_OPENCL) {
        std::fprintf(stderr, "xcorr3: -af opencl is accepted for compatibility and mapped to CUDA\n");
    } else if (args.device == XCORR2_DEVICE_CPU) {
        std::fprintf(stderr, "xcorr3: -af cpu is accepted for compatibility and mapped to CUDA\n");
    }
}

void dft_interpolate(
        const DeviceBuffer<cufftComplex> &input,
        int height,
        int width,
        int scale_h,
        int scale_w,
        int split,
        CufftPlan2D &input_plan,
        CufftPlan2D &output_plan,
        DeviceBuffer<cufftComplex> &fft_input,
        DeviceBuffer<cufftComplex> &fft_padded,
        DeviceBuffer<cufftComplex> &output) {
    int out_height = height * scale_h;
    int out_width = width * scale_w;
    int input_count = height * width;
    int output_count = out_height * out_width;

    CUDA_CHECK(cudaMemcpy(fft_input.get(), input.get(), static_cast<size_t>(input_count) * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemset(fft_padded.get(), 0, static_cast<size_t>(output_count) * sizeof(cufftComplex)));

    CUFFT_CHECK(cufftExecC2C(input_plan.handle, fft_input.get(), fft_input.get(), CUFFT_FORWARD));

    dim3 block(16, 16);
    dim3 grid(div_up(width, block.x), div_up(height, block.y));
    pad_fft_kernel<<<grid, block>>>(fft_input.get(), height, width, fft_padded.get(), out_height, out_width, split);
    CUDA_CHECK(cudaGetLastError());

    CUFFT_CHECK(cufftExecC2C(output_plan.handle, fft_padded.get(), output.get(), CUFFT_INVERSE));

    int threads = 256;
    int blocks = div_up(output_count, threads);
    scale_complex_kernel<<<blocks, threads>>>(output.get(), 1.0f / (height * width), output_count);
    CUDA_CHECK(cudaGetLastError());
}

void crop_center_columns(
        const DeviceBuffer<cufftComplex> &src,
        int src_width,
        int dst_width,
        int height,
        DeviceBuffer<cufftComplex> &dst) {
    dim3 block(16, 16);
    dim3 grid(div_up(dst_width, block.x), div_up(height, block.y));
    int x0 = src_width / 2 - dst_width / 2;
    crop_complex_kernel<<<grid, block>>>(src.get(), src_width, x0, dst_width, height, dst.get());
    CUDA_CHECK(cudaGetLastError());
}

cub::KeyValuePair<int, float> argmax_to_host(
        const DeviceBuffer<float> &input,
        int count,
        DeviceBuffer<cub::KeyValuePair<int, float>> &output,
        DeviceBuffer<unsigned char> &temp_storage,
        size_t temp_storage_bytes) {
    CUFFT_CHECK(CUFFT_SUCCESS);
    cub::DeviceReduce::ArgMax(temp_storage.get(), temp_storage_bytes, input.get(), output.get(), count);
    CUDA_CHECK(cudaGetLastError());
    cub::KeyValuePair<int, float> result;
    output.copy_to_host(&result, 1);
    return result;
}

float compute_max_corr_gpu(
        const DeviceBuffer<float> &c1r,
        const DeviceBuffer<float> &c2r,
        int nx_win,
        int xsearch,
        int ysearch,
        int nx_corr,
        int ny_corr,
        int peak_x,
        int peak_y,
        DeviceBuffer<double> &metrics) {
    metrics.zero();
    int total = nx_corr * ny_corr;
    int threads = 256;
    int blocks = std::min(div_up(total, threads), 256);
    corr_metrics_kernel<<<blocks, threads, static_cast<size_t>(threads) * 3 * sizeof(double)>>>(
        c1r.get(), c2r.get(), nx_win, xsearch, ysearch, nx_corr, ny_corr, peak_x, peak_y, metrics.get());
    CUDA_CHECK(cudaGetLastError());

    double host_metrics[3] = {0.0, 0.0, 0.0};
    metrics.copy_to_host(host_metrics, 3);
    if (host_metrics[1] == 0.0 || host_metrics[2] == 0.0) {
        return 0.0f;
    }

    return static_cast<float>(100.0 * std::fabs(host_metrics[0] / std::sqrt(host_metrics[1] * host_metrics[2])));
}

struct XcorrWorkspace {
    explicit XcorrWorkspace(const st_xcorr &xcorr)
        : nx_win(xcorr.xsearch * 4),
          ny_win(xcorr.ysearch * 4),
          nx_corr(xcorr.xsearch * 2),
          ny_corr(xcorr.ysearch * 2),
          sample_count(nx_win * ny_win),
          corr_count(nx_corr * ny_corr),
          range_interp_width(nx_win * std::max(1, xcorr.ri)),
          range_output_count(ny_win * range_interp_width),
          hi_small_count(std::max(2, xcorr.n2x) * std::max(2, xcorr.n2y)),
          hi_big_count(std::max(2, xcorr.n2x * std::max(1, xcorr.interp_factor)) *
                       std::max(2, xcorr.n2y * std::max(1, xcorr.interp_factor))),
          c1(sample_count),
          c2(sample_count),
          c1r(sample_count),
          c2r(sample_count),
          c1_fft(sample_count),
          c2_fft(sample_count),
          c3_fft(sample_count),
          corr(corr_count),
          sums(2),
          metrics(3),
          corr_argmax(1) {
        if (xcorr.ri > 1) {
            range_fft_input.allocate(sample_count);
            range_fft_padded.allocate(range_output_count);
            range_output.allocate(range_output_count);
        }

        if (xcorr.interp_factor > 1) {
            corr2_real.allocate(hi_small_count);
            corr2_complex.allocate(hi_small_count);
            hi_fft_input.allocate(hi_small_count);
            hi_fft_padded.allocate(hi_big_count);
            hi_corr_complex.allocate(hi_big_count);
            hi_corr_real.allocate(hi_big_count);
            hi_argmax.allocate(1);
        }

        size_t temp_bytes = 0;
        cub::DeviceReduce::ArgMax(nullptr, temp_bytes, corr.get(), corr_argmax.get(), corr_count);
        corr_argmax_temp.allocate(temp_bytes);
        corr_argmax_temp_bytes = temp_bytes;

        if (xcorr.interp_factor > 1) {
            temp_bytes = 0;
            cub::DeviceReduce::ArgMax(nullptr, temp_bytes, hi_corr_real.get(), hi_argmax.get(), hi_big_count);
            hi_argmax_temp.allocate(temp_bytes);
            hi_argmax_temp_bytes = temp_bytes;
        } else {
            hi_argmax_temp_bytes = 0;
        }
    }

    int nx_win;
    int ny_win;
    int nx_corr;
    int ny_corr;
    int sample_count;
    int corr_count;
    int range_interp_width;
    int range_output_count;
    int hi_small_count;
    int hi_big_count;

    DeviceBuffer<cufftComplex> c1;
    DeviceBuffer<cufftComplex> c2;
    DeviceBuffer<float> c1r;
    DeviceBuffer<float> c2r;
    DeviceBuffer<cufftComplex> c1_fft;
    DeviceBuffer<cufftComplex> c2_fft;
    DeviceBuffer<cufftComplex> c3_fft;
    DeviceBuffer<float> corr;

    DeviceBuffer<cufftComplex> range_fft_input;
    DeviceBuffer<cufftComplex> range_fft_padded;
    DeviceBuffer<cufftComplex> range_output;

    DeviceBuffer<float> corr2_real;
    DeviceBuffer<cufftComplex> corr2_complex;
    DeviceBuffer<cufftComplex> hi_fft_input;
    DeviceBuffer<cufftComplex> hi_fft_padded;
    DeviceBuffer<cufftComplex> hi_corr_complex;
    DeviceBuffer<float> hi_corr_real;

    DeviceBuffer<double> sums;
    DeviceBuffer<double> metrics;
    DeviceBuffer<cub::KeyValuePair<int, float>> corr_argmax;
    DeviceBuffer<cub::KeyValuePair<int, float>> hi_argmax;
    DeviceBuffer<unsigned char> corr_argmax_temp;
    DeviceBuffer<unsigned char> hi_argmax_temp;
    size_t corr_argmax_temp_bytes;
    size_t hi_argmax_temp_bytes;
};

struct SamplingAxis {
    std::vector<int> positions;
};

struct SamplingGrid {
    SamplingAxis x;
    SamplingAxis y;
};

/* Keep the original truncation rule, but never narrow a nonfinite/out-of-range
 * stretch or add signed ints before proving that the final center fits. */
int slave_azimuth_center(int center, int offset, double stretch) {
    const double shift = std::trunc(center * stretch);
    if (!std::isfinite(shift) || shift < std::numeric_limits<int>::min() || shift > std::numeric_limits<int>::max())
        throw std::runtime_error("PRF stretch exceeds the azimuth index range");
    const int64_t result = static_cast<int64_t>(center) + offset + static_cast<int>(shift);
    if (result < std::numeric_limits<int>::min() || result > std::numeric_limits<int>::max())
        throw std::runtime_error("ashift/PRF stretch produces an out-of-range azimuth index");
    return static_cast<int>(result);
}

SamplingAxis build_sampling_axis_original(
        int requested_count,
        int win_size,
        int search_half,
        int master_size,
        bool is_x_axis) {
    if (requested_count <= 0) {
        throw std::runtime_error(std::string(is_x_axis ? "range" : "azimuth") + " sample count must be positive");
    }

    const int corr_size = search_half * 2;
    const int inc = (master_size - 2 * (search_half + corr_size)) /
                    (requested_count + (is_x_axis ? 3 : 1));
    if (inc <= 0) {
        throw std::runtime_error(std::string("Invalid sampling increment for ") +
                                 (is_x_axis ? "range" : "azimuth") +
                                 "; requested nx/ny is too dense for the scene size and search window");
    }

    SamplingAxis axis = {};
    axis.positions.resize(static_cast<size_t>(requested_count));
    for (int k = 0; k < requested_count; ++k) {
        const int multiplier = is_x_axis ? (k + 2) : (k + 1);
        axis.positions[static_cast<size_t>(k)] = win_size + multiplier * inc;
    }
    return axis;
}

SamplingAxis build_sampling_axis_overlap(
        int requested_count,
        int window_size,
        int master_size,
        int slave_size,
        long double scale,
        int offset,
        const char *axis_name) {
    if (requested_count <= 0) {
        throw std::runtime_error(std::string(axis_name) + " sample count must be positive");
    }
    if (window_size > master_size || window_size > slave_size) {
        throw std::runtime_error(std::string("Search window is larger than the valid ") + axis_name + " image extent");
    }

    const long double half_window = static_cast<long double>(window_size) / 2.0L;
    const long double master_min = half_window;
    const long double master_max = static_cast<long double>(master_size) - half_window;
    const long double slave_min = (half_window - static_cast<long double>(offset)) / scale;
    const long double slave_max = (static_cast<long double>(slave_size) - half_window - static_cast<long double>(offset)) / scale;
    const long double effective_min = std::max(master_min, slave_min);
    const long double effective_max = std::min(master_max, slave_max);

    const int min_center = static_cast<int>(std::ceil(effective_min - 1e-12L));
    const int max_center = static_cast<int>(std::floor(effective_max + 1e-12L));
    if (min_center > max_center) {
        throw std::runtime_error(
            std::string("No valid overlapping ") + axis_name +
            " sampling region remains after applying the search window and PRM shift/stretch");
    }

    SamplingAxis axis = {};
    axis.positions.resize(static_cast<size_t>(requested_count));
    if (requested_count == 1) {
        axis.positions[0] = min_center + (max_center - min_center) / 2;
        return axis;
    }

    const unsigned long long unique_count = static_cast<unsigned long long>(max_center - min_center) + 1ULL;
    const bool require_unique = static_cast<unsigned long long>(requested_count) <= unique_count;
    const long double span = static_cast<long double>(max_center - min_center);
    for (int i = 0; i < requested_count; ++i) {
        const long double t = static_cast<long double>(i) / static_cast<long double>(requested_count - 1);
        int center = static_cast<int>(std::llround(static_cast<long double>(min_center) + t * span));

        const int remaining_slots = requested_count - 1 - i;
        const int min_allowed = (i == 0)
            ? min_center
            : (require_unique ? axis.positions[static_cast<size_t>(i) - 1] + 1
                              : axis.positions[static_cast<size_t>(i) - 1]);
        const int max_allowed = require_unique ? (max_center - remaining_slots) : max_center;
        if (center < min_allowed) center = min_allowed;
        if (center > max_allowed) center = max_allowed;
        axis.positions[static_cast<size_t>(i)] = center;
    }
    return axis;
}

bool positions_fit_overlap(
        const std::vector<int> &positions,
        int window_size,
        int master_size,
        int slave_size,
        long double scale,
        int offset) {
    const long double half_window = static_cast<long double>(window_size) / 2.0L;
    for (int center : positions) {
        const long double slave_center = scale * static_cast<long double>(center) + static_cast<long double>(offset);
        if (static_cast<long double>(center) - half_window < 0.0L) return false;
        if (static_cast<long double>(center) + half_window > static_cast<long double>(master_size)) return false;
        if (slave_center - half_window < 0.0L) return false;
        if (slave_center + half_window > static_cast<long double>(slave_size)) return false;
    }
    return true;
}

SamplingGrid build_sampling_grid(const st_xcorr &xcorr) {
    const long double scale = 1.0L + static_cast<long double>(xcorr.astretcha);
    if (scale <= 0.0L) {
        throw std::runtime_error("Invalid azimuth stretch factor; resulting slave coordinate scale is not positive");
    }

    SamplingGrid grid = {};
    grid.x = build_sampling_axis_original(xcorr.nxl, xcorr.xsearch * 4, xcorr.xsearch, xcorr.m_nx, true);
    grid.y = build_sampling_axis_original(xcorr.nyl, xcorr.ysearch * 4, xcorr.ysearch, xcorr.m_ny, false);

    /* Bug #1 fix: range (x) axis uses scale=1.0 (no astretcha), only offset */
    if (!positions_fit_overlap(grid.x.positions, xcorr.xsearch * 4, xcorr.m_nx, xcorr.s_nx, 1.0L, xcorr.x_offset)) {
        throw std::runtime_error(
            "Original GMTSAR range sampling grid would exceed the slave image bounds with the current PRM shifts/stretch. "
            "xcorr3 no longer changes the sampling grid silently. Reduce nx/ny, reduce xsearch/ysearch, or rerun with -noshift "
            "for pixel tracking on already aligned SLCs.");
    }
    if (!positions_fit_overlap(grid.y.positions, xcorr.ysearch * 4, xcorr.m_ny, xcorr.s_ny, scale, xcorr.y_offset)) {
        throw std::runtime_error(
            "Original GMTSAR azimuth sampling grid would exceed the slave image bounds with the current PRM shifts/stretch. "
            "xcorr3 no longer changes the sampling grid silently. Reduce nx/ny, reduce xsearch/ysearch, or rerun with -noshift "
            "for pixel tracking on already aligned SLCs.");
    }
    return grid;
}

struct ResourceEstimate {
    uint64_t device_bytes;
    uint64_t host_pinned_bytes;
    uint64_t output_bytes;
    uint64_t total_patches;
    size_t corr_plan_work_bytes;
    size_t range_plan_work_bytes;
    size_t hi_in_plan_work_bytes;
    size_t hi_out_plan_work_bytes;

    uint64_t total_device_bytes() const {
        uint64_t total = device_bytes;
        total = checked_add_u64(total, static_cast<uint64_t>(corr_plan_work_bytes), "cuFFT plan workspace");
        total = checked_add_u64(total, static_cast<uint64_t>(range_plan_work_bytes), "cuFFT plan workspace");
        total = checked_add_u64(total, static_cast<uint64_t>(hi_in_plan_work_bytes), "cuFFT plan workspace");
        total = checked_add_u64(total, static_cast<uint64_t>(hi_out_plan_work_bytes), "cuFFT plan workspace");
        return total;
    }
};

size_t estimate_cufft_plan_bytes(int height, int width) {
    size_t work_bytes = 0;
    CUFFT_CHECK(cufftEstimate2d(height, width, CUFFT_C2C, &work_bytes));
    return work_bytes;
}

ResourceEstimate estimate_resources(const st_xcorr &xcorr) {
    const uint64_t nx_win = checked_mul_u64(static_cast<uint64_t>(xcorr.xsearch), 4, "nx_win");
    const uint64_t ny_win = checked_mul_u64(static_cast<uint64_t>(xcorr.ysearch), 4, "ny_win");
    const uint64_t nx_corr = checked_mul_u64(static_cast<uint64_t>(xcorr.xsearch), 2, "nx_corr");
    const uint64_t ny_corr = checked_mul_u64(static_cast<uint64_t>(xcorr.ysearch), 2, "ny_corr");
    const uint64_t sample_count = checked_mul_u64(nx_win, ny_win, "sample_count");
    const uint64_t corr_count = checked_mul_u64(nx_corr, ny_corr, "corr_count");
    const uint64_t master_row_count = checked_mul_u64(static_cast<uint64_t>(xcorr.m_nx), ny_win, "master row cache");
    const uint64_t slave_row_count = checked_mul_u64(static_cast<uint64_t>(xcorr.s_nx), ny_win, "slave row cache");
    const uint64_t total_patches = checked_mul_u64(static_cast<uint64_t>(xcorr.nxl), static_cast<uint64_t>(xcorr.nyl), "total patches");

    uint64_t device_bytes = 0;
    uint64_t host_pinned_bytes = 0;

    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "workspace c1"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "workspace c2"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(float), "workspace c1r"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(float), "workspace c2r"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "workspace c1_fft"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "workspace c2_fft"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "workspace c3_fft"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(corr_count, sizeof(float), "workspace corr"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(master_row_count, sizeof(cufftComplex) * 2ULL, "master cache device"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, checked_mul_u64(slave_row_count, sizeof(cufftComplex) * 2ULL, "slave cache device"), "workspace device");
    device_bytes = checked_add_u64(device_bytes, sizeof(double) * 5ULL, "workspace device");
    device_bytes = checked_add_u64(device_bytes, sizeof(cub::KeyValuePair<int, float>) * 2ULL, "workspace device");

    host_pinned_bytes = checked_add_u64(host_pinned_bytes, checked_mul_u64(master_row_count, sizeof(cufftComplex), "master host rows"), "workspace host");
    host_pinned_bytes = checked_add_u64(host_pinned_bytes, checked_mul_u64(master_row_count, sizeof(int16_t) * 2ULL, "master host io"), "workspace host");
    host_pinned_bytes = checked_add_u64(host_pinned_bytes, checked_mul_u64(slave_row_count, sizeof(cufftComplex), "slave host rows"), "workspace host");
    host_pinned_bytes = checked_add_u64(host_pinned_bytes, checked_mul_u64(slave_row_count, sizeof(int16_t) * 2ULL, "slave host io"), "workspace host");

    if (xcorr.ri > 1) {
        const uint64_t range_interp_width = checked_mul_u64(nx_win, static_cast<uint64_t>(xcorr.ri), "range interp width");
        const uint64_t range_output_count = checked_mul_u64(ny_win, range_interp_width, "range output count");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(sample_count, sizeof(cufftComplex), "range fft input"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(range_output_count, sizeof(cufftComplex), "range fft padded"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(range_output_count, sizeof(cufftComplex), "range output"), "workspace device");
    }

    if (xcorr.interp_factor > 1) {
        const uint64_t hi_small_count = checked_mul_u64(static_cast<uint64_t>(std::max(2, xcorr.n2x)),
                                                        static_cast<uint64_t>(std::max(2, xcorr.n2y)),
                                                        "hi corr small count");
        const uint64_t hi_big_count = checked_mul_u64(static_cast<uint64_t>(std::max(2, xcorr.n2x * std::max(1, xcorr.interp_factor))),
                                                      static_cast<uint64_t>(std::max(2, xcorr.n2y * std::max(1, xcorr.interp_factor))),
                                                      "hi corr big count");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_small_count, sizeof(float), "corr2 real"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_small_count, sizeof(cufftComplex), "corr2 complex"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_small_count, sizeof(cufftComplex), "hi fft input"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_big_count, sizeof(cufftComplex), "hi fft padded"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_big_count, sizeof(cufftComplex), "hi corr complex"), "workspace device");
        device_bytes = checked_add_u64(device_bytes, checked_mul_u64(hi_big_count, sizeof(float), "hi corr real"), "workspace device");
    }

    ResourceEstimate estimate = {};
    estimate.device_bytes = device_bytes;
    estimate.host_pinned_bytes = host_pinned_bytes;
    estimate.output_bytes = checked_mul_u64(total_patches, 64ULL, "estimated output bytes");
    estimate.total_patches = total_patches;
    estimate.corr_plan_work_bytes = estimate_cufft_plan_bytes(static_cast<int>(ny_win), static_cast<int>(nx_win));
    estimate.range_plan_work_bytes = (xcorr.ri > 1)
        ? estimate_cufft_plan_bytes(static_cast<int>(ny_win), static_cast<int>(nx_win * static_cast<uint64_t>(xcorr.ri)))
        : 0;
    estimate.hi_in_plan_work_bytes = (xcorr.interp_factor > 1)
        ? estimate_cufft_plan_bytes(std::max(2, xcorr.n2y), std::max(2, xcorr.n2x))
        : 0;
    estimate.hi_out_plan_work_bytes = (xcorr.interp_factor > 1)
        ? estimate_cufft_plan_bytes(std::max(2, xcorr.n2y * std::max(1, xcorr.interp_factor)),
                                    std::max(2, xcorr.n2x * std::max(1, xcorr.interp_factor)))
        : 0;
    return estimate;
}

void validate_runtime_configuration(const st_xcorr &xcorr, const SamplingGrid &sampling) {
    if (xcorr.nxl <= 0 || xcorr.nyl <= 0) {
        throw std::runtime_error("nx and ny must both be positive");
    }
    if (!is_power_of_two_int(xcorr.xsearch) || !is_power_of_two_int(xcorr.ysearch)) {
        throw std::runtime_error("xsearch and ysearch must both be powers of two");
    }
    if (xcorr.ri <= 0 || !is_power_of_two_int(xcorr.ri)) {
        throw std::runtime_error("range_interp must be a positive power of two");
    }
    if (xcorr.interp_factor <= 0) {
        throw std::runtime_error("interp factor must be positive");
    }

    if (sampling.x.positions.empty() || sampling.y.positions.empty()) {
        throw std::runtime_error("Sampling grid is empty after overlap validation");
    }
}

/* Host memory available for new allocations.  Prefer MemAvailable from
 * /proc/meminfo: it counts reclaimable page cache, which the kernel frees on
 * demand.  Bare _SC_AVPHYS_PAGES (~MemFree) rejects runs on machines whose
 * cache is warm (e.g. right after copying the input SLCs) even though the
 * allocation would succeed.  Fall back to sysconf when /proc is unreadable. */
static uint64_t host_available_bytes() {
    std::ifstream f("/proc/meminfo");
    std::string key, unit;
    uint64_t kb = 0;
    while (f >> key >> kb >> unit) {
        if (key == "MemAvailable:")
            return kb * 1024ULL;
    }
    const long host_pages = sysconf(_SC_AVPHYS_PAGES);
    const long page_size = sysconf(_SC_PAGE_SIZE);
    if (host_pages > 0 && page_size > 0)
        return static_cast<uint64_t>(host_pages) * static_cast<uint64_t>(page_size);
    return 0;
}

void check_resource_limits(const ResourceEstimate &estimate) {
    size_t gpu_free = 0;
    size_t gpu_total = 0;
    CUDA_CHECK(cudaMemGetInfo(&gpu_free, &gpu_total));

    const uint64_t required_gpu_bytes = estimate.total_device_bytes();
    const uint64_t safe_gpu_budget = static_cast<uint64_t>(gpu_free) * 85ULL / 100ULL;
    if (required_gpu_bytes > safe_gpu_budget) {
        throw std::runtime_error(
            "Estimated GPU memory requirement " + format_bytes(required_gpu_bytes) +
            " exceeds the safe budget of currently free GPU memory " + format_bytes(static_cast<uint64_t>(gpu_free)) +
            ". Reduce xsearch/ysearch/interp settings or free GPU memory.");
    }

    const uint64_t host_available = host_available_bytes();
    if (host_available > 0) {
        const uint64_t safe_host_budget = host_available * 60ULL / 100ULL;
        if (estimate.host_pinned_bytes > safe_host_budget) {
            throw std::runtime_error(
                "Estimated pinned host memory requirement " + format_bytes(estimate.host_pinned_bytes) +
                " exceeds 60% of available host memory (MemAvailable " + format_bytes(host_available) + ").");
        }
    }

    struct statvfs fs = {};
    if (statvfs(".", &fs) == 0) {
        const uint64_t disk_available = static_cast<uint64_t>(fs.f_bavail) * static_cast<uint64_t>(fs.f_frsize);
        const uint64_t safe_disk_budget = disk_available * 95ULL / 100ULL;
        if (estimate.output_bytes > safe_disk_budget) {
            throw std::runtime_error(
                "Estimated freq_xcorr.dat size " + format_bytes(estimate.output_bytes) +
                " exceeds the free disk space in the current directory (" + format_bytes(disk_available) +
                "). Reduce nx/ny or run in a filesystem with more space.");
        }
    }
}

int main(int argc, char **argv) {
    try {
        std::signal(SIGXFSZ, SIG_IGN); // Report file-size/write failures through checked I/O.
        st_xcorr_args args;
        st_xcorr xcorr;

        parse_opts(&args, argc, argv);
        apply_args(&args, &xcorr);
        std::unique_ptr<char, decltype(&std::free)> master_path(xcorr.m_path, &std::free);
        std::unique_ptr<char, decltype(&std::free)> slave_path(xcorr.s_path, &std::free);
        std::vector<std::string> products = {"freq_xcorr.dat"};
        if (xcorr.do_geocode) products.insert(products.end(), {"azi_offset.grd", "rng_offset.grd", "azi_offset_ll.grd", "rng_offset_ll.grd"});
        OutputTransaction outputs(products);
        SamplingGrid sampling = build_sampling_grid(xcorr);
        validate_runtime_configuration(xcorr, sampling);
        for (int y : sampling.y.positions) {
            const int64_t center = slave_azimuth_center(y, xcorr.y_offset, xcorr.astretcha);
            if (center - 2LL * xcorr.ysearch < 0 || center + 2LL * xcorr.ysearch > xcorr.s_ny)
                throw std::runtime_error("ashift/PRF stretch places the azimuth window outside the secondary SLC");
        }
        initialize_cuda(args);
        cudaDeviceProp device_properties;
        CUDA_CHECK(cudaGetDeviceProperties(&device_properties, 0));
        if (div_up(4 * xcorr.xsearch, 16) > device_properties.maxGridSize[0] ||
            div_up(4 * xcorr.ysearch, 16) > device_properties.maxGridSize[1])
            throw std::runtime_error("xsearch/ysearch exceeds this GPU's two-dimensional launch grid limits");
        ResourceEstimate estimate = estimate_resources(xcorr);
        check_resource_limits(estimate);

        std::ifstream fmaster(xcorr.m_path, std::ios::binary);
        std::ifstream fslave(xcorr.s_path, std::ios::binary);
        if (!fmaster) {
            throw std::runtime_error(std::string("Failed to open master SLC image: ") + xcorr.m_path);
        }
        if (!fslave) {
            throw std::runtime_error(std::string("Failed to open slave SLC image: ") + xcorr.s_path);
        }

        const int xsearch = xcorr.xsearch;
        const int ysearch = xcorr.ysearch;
        const int nx_corr = xcorr.xsearch * 2;
        const int ny_corr = xcorr.ysearch * 2;
        const int nx_win = nx_corr * 2;
        const int ny_win = ny_corr * 2;
        const int sample_count = nx_win * ny_win;
        const std::vector<int> &x_positions = sampling.x.positions;
        const std::vector<int> &y_positions = sampling.y.positions;
        const double slave_scale = 1.0 + xcorr.astretcha;
        const size_t total_patches = static_cast<size_t>(xcorr.nxl) * static_cast<size_t>(xcorr.nyl);
        size_t processed_patches = 0;
        auto start_time = std::chrono::steady_clock::now();
        auto last_progress = start_time;

        CufftPlan2D corr_plan(ny_win, nx_win);
        std::unique_ptr<CufftPlan2D> range_out_plan;
        if (xcorr.ri > 1) {
            range_out_plan.reset(new CufftPlan2D(ny_win, nx_win * std::max(1, xcorr.ri)));
        }
        std::unique_ptr<CufftPlan2D> hi_in_plan;
        std::unique_ptr<CufftPlan2D> hi_out_plan;
        if (xcorr.interp_factor > 1) {
            hi_in_plan.reset(new CufftPlan2D(std::max(2, xcorr.n2y), std::max(2, xcorr.n2x)));
            hi_out_plan.reset(new CufftPlan2D(std::max(2, xcorr.n2y * std::max(1, xcorr.interp_factor)),
                                              std::max(2, xcorr.n2x * std::max(1, xcorr.interp_factor))));
        }

        XcorrWorkspace workspace(xcorr);
        SlcRowCache master_rows(fmaster, xcorr.m_nx, ny_win);
        SlcRowCache slave_rows(fslave, xcorr.s_nx, ny_win);

        CheckedFile output(outputs.path("freq_xcorr.dat"));
        FILE *fout = output.get();

        render_progress(0, total_patches, start_time, false);

        dim3 block2d(16, 16);
        int block1d = 256;
        int sum_blocks = std::min(div_up(sample_count, block1d), 256);

        for (int j = 0; j < xcorr.nyl; ++j) {
            int loc_y = y_positions[static_cast<size_t>(j)];
            /* Bug #1 fix: astretcha only applies to azimuth (y), not range (x).
             * Match original GMTSAR: ishft = (int)(loc_y * astretcha)
             *   slave_y = loc_y + y_offset + ishft
             *   slave_x = loc_x + x_offset              (no astretcha)       */
            int slave_loc_y = slave_azimuth_center(loc_y, xcorr.y_offset, xcorr.astretcha);
            master_rows.ensure_window(loc_y - ny_win / 2);
            slave_rows.ensure_window(slave_loc_y - ny_win / 2);

            for (int i = 0; i < xcorr.nxl; ++i) {
                int loc_x = x_positions[static_cast<size_t>(i)];
                int slave_loc_x = loc_x + xcorr.x_offset;

                dim3 patch_grid(div_up(nx_win, block2d.x), div_up(ny_win, block2d.y));
                crop_complex_kernel<<<patch_grid, block2d>>>(
                    master_rows.data(), xcorr.m_nx, loc_x - nx_win / 2, nx_win, ny_win, workspace.c1.get());
                crop_complex_kernel<<<patch_grid, block2d>>>(
                    slave_rows.data(), xcorr.s_nx, slave_loc_x - nx_win / 2, nx_win, ny_win, workspace.c2.get());
                CUDA_CHECK(cudaGetLastError());

                if (xcorr.ri > 1) {
                    dft_interpolate(workspace.c1, ny_win, nx_win, 1, xcorr.ri, xcorr.nyquist_split,
                                    corr_plan, *range_out_plan,
                                    workspace.range_fft_input, workspace.range_fft_padded, workspace.range_output);
                    crop_center_columns(workspace.range_output, workspace.range_interp_width, nx_win, ny_win, workspace.c1);

                    dft_interpolate(workspace.c2, ny_win, nx_win, 1, xcorr.ri, xcorr.nyquist_split,
                                    corr_plan, *range_out_plan,
                                    workspace.range_fft_input, workspace.range_fft_padded, workspace.range_output);
                    crop_center_columns(workspace.range_output, workspace.range_interp_width, nx_win, ny_win, workspace.c2);
                    }

                complex_abs_kernel<<<div_up(sample_count, block1d), block1d>>>(workspace.c1.get(), workspace.c1r.get(), sample_count);
                complex_abs_kernel<<<div_up(sample_count, block1d), block1d>>>(workspace.c2.get(), workspace.c2r.get(), sample_count);
                CUDA_CHECK(cudaGetLastError());

                workspace.sums.zero();
                sum_two_arrays_kernel<<<sum_blocks, block1d, static_cast<size_t>(block1d) * 2 * sizeof(double)>>>(
                    workspace.c1r.get(), workspace.c2r.get(), sample_count, workspace.sums.get());
                CUDA_CHECK(cudaGetLastError());

                double host_sums[2] = {0.0, 0.0};
                workspace.sums.copy_to_host(host_sums, 2);
                float mean1 = static_cast<float>(host_sums[0] / sample_count);
                float mean2 = static_cast<float>(host_sums[1] / sample_count);

                center_and_mask_kernel<<<patch_grid, block2d>>>(
                    workspace.c1r.get(), workspace.c2r.get(), nx_win, ny_win, xsearch, ysearch, mean1, mean2);
                CUDA_CHECK(cudaGetLastError());

                real_to_complex_kernel<<<div_up(sample_count, block1d), block1d>>>(workspace.c1r.get(), workspace.c1_fft.get(), sample_count);
                real_to_complex_kernel<<<div_up(sample_count, block1d), block1d>>>(workspace.c2r.get(), workspace.c2_fft.get(), sample_count);
                CUDA_CHECK(cudaGetLastError());

                CUFFT_CHECK(cufftExecC2C(corr_plan.handle, workspace.c1_fft.get(), workspace.c1_fft.get(), CUFFT_FORWARD));
                CUFFT_CHECK(cufftExecC2C(corr_plan.handle, workspace.c2_fft.get(), workspace.c2_fft.get(), CUFFT_FORWARD));
                freq_product_kernel<<<patch_grid, block2d>>>(workspace.c1_fft.get(), workspace.c2_fft.get(), ny_win, nx_win, workspace.c3_fft.get());
                CUDA_CHECK(cudaGetLastError());
                CUFFT_CHECK(cufftExecC2C(corr_plan.handle, workspace.c3_fft.get(), workspace.c3_fft.get(), CUFFT_INVERSE));

                dim3 corr_grid(div_up(nx_corr, block2d.x), div_up(ny_corr, block2d.y));
                scaled_abs_crop_kernel<<<corr_grid, block2d>>>(
                    workspace.c3_fft.get(), nx_win, xsearch, ysearch, nx_corr, ny_corr,
                    1.0f / (nx_win * ny_win), workspace.corr.get());
                CUDA_CHECK(cudaGetLastError());

                cub::KeyValuePair<int, float> corr_peak = argmax_to_host(
                    workspace.corr, workspace.corr_count,
                    workspace.corr_argmax, workspace.corr_argmax_temp, workspace.corr_argmax_temp_bytes);

                int corr_max_idx = corr_peak.key;
                float cmax = corr_peak.value;
                int peak_x = corr_max_idx % nx_corr - xsearch;
                int peak_y = corr_max_idx / nx_corr - ysearch;
                float max_corr = compute_max_corr_gpu(
                    workspace.c1r, workspace.c2r,
                    nx_win, xsearch, ysearch, nx_corr, ny_corr,
                    peak_x, peak_y, workspace.metrics);

                /* Peak-significance SNR (scheme #1): z-score of the correlation
                 * peak above the surface background = (cmax - mean)/std.  A
                 * sharp, unique peak scores high; a spurious peak lost among
                 * comparable side-lobes scores low, letting downstream filtering
                 * discard gross mismatches that the plain correlation value keeps.
                 * Scale-invariant (cmax/mean/std share the FFT scale factor). */
                workspace.sums.zero();
                sum_sumsq_kernel<<<std::min(div_up(workspace.corr_count, block1d), 256), block1d,
                                   static_cast<size_t>(block1d) * 2 * sizeof(double)>>>(
                    workspace.corr.get(), workspace.corr_count, workspace.sums.get());
                CUDA_CHECK(cudaGetLastError());
                double corr_stats[2] = {0.0, 0.0};
                workspace.sums.copy_to_host(corr_stats, 2);
                double cave = corr_stats[0] / workspace.corr_count;
                double cvar = corr_stats[1] / workspace.corr_count - cave * cave;
                double cstd = cvar > 0.0 ? std::sqrt(cvar) : 0.0;
                float peak_snr = (cstd > 0.0)
                    ? static_cast<float>((static_cast<double>(cmax) - cave) / cstd)
                    : 0.0f;

                float xfrac = 0.0f;
                float yfrac = 0.0f;

                if (xcorr.interp_factor > 1) {
                    int factor = xcorr.interp_factor;
                    int nx_corr2 = xcorr.n2x;
                    int ny_corr2 = xcorr.n2y;

                    if (peak_y + ysearch < ny_corr2 / 2) {
                        peak_y = ny_corr2 / 2 - ysearch;
                    } else if (peak_y + ysearch >= ny_corr - ny_corr2 / 2) {
                        peak_y = ny_corr - ny_corr2 / 2 - ysearch - 1;
                    }

                    if (peak_x + xsearch < nx_corr2 / 2) {
                        peak_x = nx_corr2 / 2 - xsearch;
                    } else if (peak_x + xsearch >= nx_corr - nx_corr2 / 2) {
                        peak_x = nx_corr - nx_corr2 / 2 - xsearch - 1;
                    }

                    float scale = (cmax != 0.0f) ? (max_corr / cmax) : 0.0f;
                    dim3 corr2_grid(div_up(nx_corr2, block2d.x), div_up(ny_corr2, block2d.y));
                    crop_scale_pow_quarter_kernel<<<corr2_grid, block2d>>>(
                        workspace.corr.get(), nx_corr,
                        peak_x + xsearch - nx_corr2 / 2,
                        peak_y + ysearch - ny_corr2 / 2,
                        nx_corr2, ny_corr2, scale,
                        workspace.corr2_real.get());
                    CUDA_CHECK(cudaGetLastError());

                    real_to_complex_kernel<<<div_up(workspace.hi_small_count, block1d), block1d>>>(
                        workspace.corr2_real.get(), workspace.corr2_complex.get(), workspace.hi_small_count);
                    CUDA_CHECK(cudaGetLastError());

                    dft_interpolate(workspace.corr2_complex, ny_corr2, nx_corr2, factor, factor, xcorr.nyquist_split,
                                    *hi_in_plan, *hi_out_plan,
                                    workspace.hi_fft_input, workspace.hi_fft_padded, workspace.hi_corr_complex);
                    /* Bug #4 fix: use real part for peak finding, matching
                     * original GMTSAR highres_corr.c line 53: cd_exp[].r */
                    complex_real_kernel<<<div_up(workspace.hi_big_count, block1d), block1d>>>(
                        workspace.hi_corr_complex.get(), workspace.hi_corr_real.get(), workspace.hi_big_count);
                    CUDA_CHECK(cudaGetLastError());

                    cub::KeyValuePair<int, float> hi_peak = argmax_to_host(
                        workspace.hi_corr_real, workspace.hi_big_count,
                        workspace.hi_argmax, workspace.hi_argmax_temp, workspace.hi_argmax_temp_bytes);

                    int nx_hi = nx_corr2 * factor;
                    int ny_hi = ny_corr2 * factor;
                    int peak_x2 = hi_peak.key % nx_hi - nx_hi / 2;
                    int peak_y2 = hi_peak.key / nx_hi - ny_hi / 2;
                    xfrac = peak_x2 / static_cast<float>(factor);
                    yfrac = peak_y2 / static_cast<float>(factor);
                }

                /* Bug #1 fix (output): match original truncation of astretcha to int */
                float xoff = xcorr.x_offset - ((peak_x + xfrac) / xcorr.ri);
                float yoff = xcorr.y_offset - (peak_y + yfrac) + static_cast<float>(static_cast<int>(loc_y * xcorr.astretcha));
                if (std::fprintf(fout, " %d %6.3lf %d %6.3lf %6.2lf %6.2lf \n",
                             loc_x, static_cast<double>(xoff),
                             loc_y, static_cast<double>(yoff),
                             static_cast<double>(max_corr),
                             static_cast<double>(peak_snr)) < 0)
                    throw io_error("write correlation output");

                processed_patches++;
                if (processed_patches == total_patches || processed_patches % 64 == 0) {
                    auto now = std::chrono::steady_clock::now();
                    if (processed_patches == total_patches || now - last_progress >= std::chrono::milliseconds(500)) {
                        render_progress(processed_patches, total_patches, start_time, processed_patches == total_patches);
                        last_progress = now;
                    }
                }

                /* GPU load headroom: on a single GPU that also drives the display,
                 * back-to-back kernels starve the desktop compositor (screen
                 * freezes) and hold the card at max power/heat. Drain the queue
                 * and rest briefly every throttle_every patches so the display
                 * gets GPU time and thermals/power ease off. Skipped entirely
                 * when throttle_ms == 0 (e.g. display on the iGPU). */
                if (xcorr.throttle_ms > 0 &&
                    processed_patches != total_patches &&
                    processed_patches % static_cast<size_t>(xcorr.throttle_every) == 0) {
                    CUDA_CHECK(cudaDeviceSynchronize());
                    std::this_thread::sleep_for(std::chrono::milliseconds(xcorr.throttle_ms));
                }
            }
        }

        CUDA_CHECK(cudaDeviceSynchronize());
        output.close();
        fmaster.close();
        fslave.close();
        run_postprocess(xcorr, args.m_prm, sampling.x.positions, sampling.y.positions, outputs);
        outputs.publish();
        return 0;
    } catch (const std::exception &ex) {
        std::fprintf(stderr, "xcorr3: %s\n", ex.what());
        return 1;
    }
}
