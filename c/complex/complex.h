#ifndef COMPLEX_H
#define COMPLEX_H

#include <complex.h>

typedef struct {
    int rows, cols;
    double complex *data;
} carray;

/*
Frees the memory of the given matrix

Input:
    carray *m: pointer to a matrix

Output:
    None
*/
void cmatrix_free(carray *m);

/*
Prints the given matrix to console

Input:
    carray *m: pointer to a matrix

Output:
    None, but the matrix is displayed in the console
*/
void cmatrix_print(const carray *m);

/*
Prints the given matrix to console with given precision

Input:
    carray *m: pointer to a matrix
    int precision: amount of decimals should be printed

Output:
    None, but the matrix is displayed in the console with given decimal precision
*/
void cmatrix_print_precision(const carray *m, int precision);

/*
Creates a new complex matrix with given dimension (rows x cols) with value zero in each element

Input:
    int rows: amount of rows the matrix will have
    int cols: amount of columns the matrix will have

Output:
    marray*: pointer to the created complex matrix with given dimensions and value zero in each element
*/
carray *cmatrix_zeroes(int rows, int cols);

/*
Creates a new complex matrix with one on diagonal elements and zero elsewhere

Input:
    int n: dimensions for quadratic matrix

Output:
    carray*: pointer to the created identity matrix
*/
carray *cmatrix_identity(int n);
#endif //COMPLEX_H
