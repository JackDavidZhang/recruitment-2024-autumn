#include "cuda.h"

void __host__ cuda_caluc1(int * dst, unsigned * src, int n,int num) {
  int *d_dst,*d_src;
  cudaMalloc((void**)&d_dst,n*sizeof(int));
  cudaMalloc((void**)&d_src,n*sizeof(unsigned));
  cudaMemcpy(d_src,src,n*sizeof(unsigned),cudaMemcpyHostToDevice);
  cuda_cal1<<<n/1024+1,1024>>>(d_dst,d_src,n,num);
  cudaMemcpy(dst,d_dst,n*sizeof(int),cudaMemcpyDeviceToHost);
}

void __global__ cuda_cal1(int * dst, unsigned * src, int n,int num) {
   int index = blockIdx.x*blockDim.x+threadIdx.x;
   if(index<n) {
     dst[index] = src[index]+num;
    }
}