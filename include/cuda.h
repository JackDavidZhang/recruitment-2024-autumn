//
// Created by ztsubaki on 24-10-18.
//

#ifndef CUDA_H
#define CUDA_H

void __host__ cuda_caluc1(int * dst, unsigned * src, int n,int num);
void __global__ cuda_cal1(int * dst, unsigned * src, int n,int num);

#endif //CUDA_H
