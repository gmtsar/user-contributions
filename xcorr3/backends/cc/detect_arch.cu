// Print the compute capability of the first visible CUDA device as "sm_XY".
// Used by the Makefile to auto-select -arch for the local GPU.
#include <cstdio>
#include <cuda_runtime.h>
int main(void) {
    int dev = 0;
    cudaDeviceProp prop;
    if (cudaGetDevice(&dev) != cudaSuccess) return 1;
    if (cudaGetDeviceProperties(&prop, dev) != cudaSuccess) return 1;
    std::printf("sm_%d%d\n", prop.major, prop.minor);
    return 0;
}
