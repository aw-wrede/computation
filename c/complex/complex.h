#ifndef COMPLEX_H
#define COMPLEX_H

#include <complex.h>
#include <stdbool.h>

#define COMPLEX_NORM int
#define COMPLEX_NORM_UNSCALED 0
#define COMPLEX_NORM_ORTHO 1

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
    carray*: pointer to the created complex matrix with given dimensions and value zero in each element
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

/*
Creates a copy of the complex matrix

Input:
    carray *m: pointer to the matrix to copy

Output:
    carray*: pointer to the copy of the input matrix
*/
carray *cmatrix_copy(const carray *m);

/*
Creates the transposed of the given complex matrix

Input:
    carray *m: matrix from which the transposed is to be created.

Output:
    carray*: pointer to the transposed
*/
carray *cmatrix_transposed(const carray *m);

/*
Creates the complex-conjugate of the given matrix

Input:
    carray *m: matrix from which the complex-conjugate is to be created.

Output:
    carray*: pointer to the complex-conjugate
*/
carray * cmatrix_conj(const carray *m);

/*
Creates the adjoint (conjugate transpose) of the given matrix

Input:
    carray *m: matrix from which the adjoint is to be created.

Output:
    carray*: pointer to the adjoint
*/
carray *cmatrix_adjoint(const carray *m);

/*
Adds the value b to the complex matrix a

Input:
    carray *a: matrix to which the value b is added
    double complex b: complex value that is added to the matrix a

Output:
    None, but the matrix a is the result of the addition
*/
void cmatrix_addi_val(const carray *a, double complex b);

/*
Calculates the result of the complex matrix a, to which value b is added

Input:
    carray *a: matrix to which the value b is added
    double complex b: complex value that is added to the matrix a

Output:
    carray*: pointer to the result of the addition of matrix a and value b
*/
carray *cmatrix_add_val(const carray *a, double complex b);

/*
Adds the second matrix b to the first matrix a

Input:
    carray *a: matrix to which the matrix b is added
    carray *b: matrix that is added to the matrix a

Output:
    None, but the matrix a is the result of the matrix addition of a and b
*/
void cmatrix_addi(const carray *a, const carray *b);

/*
Calculates the complex matrix addition of a and b

Input:
    carray *a: matrix to which the matrix b is added
    carray *b: matrix that is added to the matrix a

Output:
    carray*: pointer to the result of the matrix addition of a and b
*/
carray *cmatrix_add(const carray *a, const carray *b);

/*
Subtracts the second matrix b from the first matrix a

Input:
    carray *a: matrix to which the matrix b is subtracted
    carray *b: matrix that is subtracted from the matrix a

Output:
    None, but the matrix a is the result of the matrix subtraction of a and b
*/
void cmatrix_subi(const carray *a, const carray *b);

/*
Calculates the complex matrix subtraction of a and b

Input:
    carray *a: matrix to which the matrix b is subtracted
    carray *b: matrix that is subtracted from the matrix a

Output:
    carray*: pointer to the result of the matrix addition of a and b
*/
carray *cmatrix_sub(const carray *a, const carray *b);

/*
Multiples each element of the complex matrix a with b

Input:
    carray *a: matrix which will be multiplied with b
    double complex b: complex value that the matrix a will be multiplied with

Output:
    None, but the matrix a is the result of the multiplication
*/
void cmatrix_muli_val(const carray *a, double complex b);

/*
Multiples each element of matrix a with b

Input:
    carray *a: matrix which will be multiplied with b
    double complex b: value that the matrix a will be multiplied with

Output:
    carray*: pointer to the matrix containing the result of the multiplication
*/
carray *cmatrix_mul_val(const carray *a, double complex b);

/*
Creates a partition or sub-area of the given matrix

Input:
    carray **dest: pointer in which the new matrix partition pointer is to be saved
    carray *m: matrix from which a partition is to be created
    row_start: row index to start (included)
    row_end: row index to end (excluded)
    col_start: column index to start (included)
    col_end: column index to end (excluded)
Output:
    None, but the given pointer dest contains the pointer in which the new matrix partition pointer is to be saved
*/
void cmatrix_get_partition(carray **dest, const carray *m, int row_start, int row_end, int col_start, int col_end);

/*
Calculates the dot product of two matrices

Input:
    carray *a: a pointer to a matrix
    carray *b: a pointer to a matrix

Output:
    carray*: pointer to the result of the matrix dot product a*b
*/
carray *cmatrix_dot(const carray *a, const carray *b);

/*
Checks whether the individual elements of both matrices are close to each other and returns the result per element

Input:
    carray *a: first matrix
    carray *b: second matrix, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    carray*: matrix containing the result of each element
*/
carray *cmatrix_close(const carray *a, const carray *b, double rtol, double atol);

/*
Checks whether the individual elements of both matrices are close to each other

Input:
    carray *a: first matrix
    carray *b: second matrix, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    bool: result of whether the matrices are close to each other
*/
bool cmatrix_close_all(const carray *a, const carray *b, double rtol, double atol);

/*
Creates the dft matrix with given dimension

Input:
    COMPLEX_NORM n: dimension of the dft matrix

Output:
    carray*: pointer to the created dft matrix of size n*n
*/
carray *dft_matrix(int n, COMPLEX_NORM norm);

/*
Checks if a matrix is unitary

Input:
    carray *m: matrix to be checked
    double rtol: relative tolerance to the unit matrix
    double atol: absolute tolerance to the unit matrix

Output:
    bool: whether the matrix is unitary
*/
bool cmatrix_is_unitary(const carray *m, double rtol, double atol);

/*
Performs the Discrete Fourier Transform on on real values

Input:
    carray *m: real data to be transformed
    COMPLEX_NORM norm: Normalization mode, default is unscaled

Output:
    carray*: pointer to the fourier transformed data
*/
carray *cmatrix_dft(const carray *m, COMPLEX_NORM norm);

/*
Performs the Fourier Transform on real values
Is fastest when the input length is a power of two

Input:
    carray *data: real data to be transformed
    COMPLEX_NORM: Normalization mode

Output:
    carray*: pointer to the fourier transformed data
*/
carray *cmatrix_fft(const carray *data, COMPLEX_NORM norm);

#endif //COMPLEX_H
