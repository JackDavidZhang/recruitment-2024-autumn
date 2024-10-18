#include "cuda.h"

void __global__ cuda_cal1(int *dst, unsigned *src, int n, int num);
void __global__ cuda_cal2(int *dst, unsigned *src, int n, int num1, int num2);

void __host__ cuda_caluc1(int *dst, unsigned *src, int n, int num) {
  int *d_dst;
  unsigned *d_src;
  cudaMalloc((void **)&d_dst, n * sizeof(int));
  cudaMalloc((void **)&d_src, n * sizeof(unsigned));
  cudaMemcpy(d_src, src, n * sizeof(unsigned), cudaMemcpyHostToDevice);
  cuda_cal1<<<(n / 1024) + 1, 1024>>>(d_dst, d_src, n, num);
  cudaDeviceSynchronize();
  cudaMemcpy(dst, d_dst, n * sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_dst);
  cudaFree(d_src);
}

void __host__ cuda_caluc2(int *dst, unsigned *src, int n, int num1, int num2) {
  int *d_dst;
  unsigned *d_src;
  cudaMalloc((void **)&d_dst, n * sizeof(int));
  cudaMalloc((void **)&d_src, n * sizeof(unsigned));
  cudaMemcpy(d_src, src, n * sizeof(unsigned), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dst, dst, n * sizeof(unsigned), cudaMemcpyHostToDevice);
  cuda_cal2<<<(n / 1024) + 1, 1024>>>(d_dst, d_src, n, num1, num2);
  cudaDeviceSynchronize();
  cudaMemcpy(dst, d_dst, n * sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_dst);
  cudaFree(d_src);
}

void __global__ cuda_cal1(int *dst, unsigned *src, int n, int num) {
  unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < n) {
    dst[index] = src[index] + num;
  }
}

void __global__ cuda_cal2(int *dst, unsigned *src, int n, int num1, int num2) {
  unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < n) {
    dst[index] = dst[index] * num1 + num2 + src[index];
  }
}