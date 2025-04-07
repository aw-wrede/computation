#include "marrayCuda.cuh"

#include <cstdio>
#include <cuda_runtime.h>
#include <iostream>
#include <ostream>


/**
 * Cuda kernel for dot product
 * @param result pointer to the result
 * @param a pointer to data of matrix a
 * @param b pointer to data of matrix b
 * @param aRows amount of rows matrix a has
 * @param aCols amount of columns matrix a has
 * @param bCols amount of columns matrix b has
*/
__global__ void matMulKernel(double *result, const double *a, const double *b, const int aRows, const int aCols, const int bCols) {
    const unsigned int row = blockIdx.y * blockDim.y + threadIdx.y; // row index for result
    const unsigned int col = blockIdx.x * blockDim.x + threadIdx.x; // column index for result

    if (row >= aRows || col >= bCols) {
        return;
    }

    double sum = 0.0;

    for (int k = 0; k < aCols; k++) {
        sum += a[row * aCols + k] * b[k * bCols + col];
    }

    result[row * bCols + col] = sum;
}


namespace marrayCuda {

    void free(marray *m) {
        if (m == nullptr) return;

        cudaFree(m);
    }


    void print(const marray *m) {
        if (m == nullptr) {
            return;
        }

        for (int i = 0; i < m->rows; i++) {
            for (int j = 0; j < m->cols; j++) {
                printf("%lf\t", m->data[i * m->cols + j]);
            }

            printf("\n");
        }
    }


    marray *zeroes(const int rows, const int cols) {
        // alloc mem for structure and data
        // calloc sets the entire memory area to 0
        marray *m;
        const cudaError_t err = cudaMallocManaged(&m, sizeof(marray) + rows * cols * sizeof(double));

        // check if mallocManaged is supported or other errors occur
        // TODO add fallback cudaMalloc
        if (err != cudaSuccess) {
            std::cerr << "Failed to allocate marray memory!\n" << cudaGetLastError() << "\n";
            return nullptr;
        }

        cudaMemset(m, 0, sizeof(marray) + rows * cols * sizeof(double));

        // return null if not enough memory is free
        if (m == nullptr) {
            return nullptr;
        }

        // apply dimensions
        m->rows = rows;
        m->cols = cols;

        // set data pointer to memory after structure
        m->data = reinterpret_cast<double *>(m + 1);

        return m;
    }


    void dotGeneral(const marray *result, const marray *a, const marray *b) {
        if (a->cols != b->rows) {
            return;
        }

        // TODO add support for fallback cudaMalloc
        //const size_t size_a = a->rows * a->cols * sizeof(double);
        //const size_t size_b = b->rows * b->cols * sizeof(double);
        //const size_t size_result = a->rows * b->cols * sizeof(double);

        // memory for matrices in cuda ram
        //double *d_a, *d_b, *d_result;
        //cudaMalloc(&d_a, size_a);
        //cudaMalloc(&d_b, size_b);
        //cudaMalloc(&d_result, size_result);

        // copy matrices data to cuda memory
        //cudaMemcpy(d_a, a->data, size_a, cudaMemcpyHostToDevice);
        //cudaMemcpy(d_b, b->data, size_b, cudaMemcpyHostToDevice);

        constexpr dim3 blockDim(16, 16);
        dim3 gridDim( (b->cols + blockDim.x - 1) / blockDim.x, (a->rows + blockDim.y - 1) / blockDim.y );

        // call cuda kernel
        matMulKernel<<<gridDim, blockDim>>>(result->data, a->data, b->data, a->rows, a->cols, b->cols);

        // wait for completion
        cudaDeviceSynchronize();

        // copy results
        //cudaMemcpy(result->data, d_result, size_result, cudaMemcpyDeviceToHost);

        // free cuda memory
        //cudaFree(d_a);
        //cudaFree(d_b);
        //cudaFree(d_result);
    }

}



