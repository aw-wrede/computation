#ifndef MARRAY_CUDA_CUH
#define MARRAY_CUDA_CUH

namespace marrayCuda {

    typedef struct {
        int rows, cols;
        double *data;
    } marray;


    /**
     * Frees the memory of the given matrix
     * @param m pointer to a matrix
     */
    void free(marray *m);


    /**
     * Prints the given matrix to console
     * @param m pointer to a matrix
     */
    void print(const marray *m);


    /**
     * Creates a new matrix with given dimension (rows x cols) with value zero in each element
     * @param rows amount of rows the matrix will have
     * @param cols amount of columns the matrix will have
     * @return pointer to the created matrix with given dimensions and value zero in each element
     */
    marray *zeroes(int rows, int cols);


    /**
     * Calculates the dot product of two matrices
     * @param result pointer to the result of the matrix dot product a*b
     * @param a pointer to a matrix
     * @param b pointer to a matrix
     */
    void dotGeneral(const marray *result, const marray *a, const marray *b);

}

#endif //MARRAY_CUDA_CUH
