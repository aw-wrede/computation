#ifndef DTYPE_MATRIX_H
#define DTYPE_MATRIX_H

#include <stdbool.h>

#define MATRIX_DTYPE char /* type of matrix data */
#define MATRIX_DTYPE_DEFAULT 2 /* default data type of matrix (double) */
#define MATRIX_DTYPE_FLOAT 0 /* float */
#define MATRIX_DTYPE_DOUBLE 1 /* double */
#define MATRIX_DTYPE_LONG_DOUBLE 2 /* long double */
#define MATRIX_DTYPE_COMPLEX_FLOAT 3 /* float complex */
#define MATRIX_DTYPE_COMPLEX_DOUBLE 4 /* double complex */
#define MATRIX_DTYPE_COMPLEX_LONG_DOUBLE 5 /* long double complex */

#define MATRIX_PRINT_PRECISION_DEFAULT 2

/* data structure of matrix */
typedef struct {
    int rows, cols; /* dimensions of matrix */
    MATRIX_DTYPE dtype; /* data type of matrix data */
    void *data; /* pointer to data */
} tarray;

/*
Frees the memory of the given matrix

Input:
    marray *m: pointer to a matrix

Output:
    None
*/
void tarray_free(tarray *m);

/*
Prints the given matrix to console

Input:
    marray *m: pointer to a matrix

Output:
    None, but the matrix is displayed in the console
*/
void tarray_print(const tarray *m);

/*
Prints the given matrix to console with given precision

Input:
    marray *m: pointer to a matrix
    int precision: amount of decimals should be printed

Output:
    None, but the matrix is displayed in the console with given decimal precision
*/
void tarray_print_precision(const tarray *m, int precision);

/*
Creates a new matrix with given dimension (rows x cols) with value zero in each element

Input:
    int rows: amount of rows the matrix will have
    int cols: amount of columns the matrix will have

Output:
    marray*: pointer to the created matrix with given dimensions and value zero in each element
*/
tarray *tarray_zeroes(int rows, int cols, MATRIX_DTYPE dtype);

/*
Creates a new matrix with one on diagonal elements and zero elsewhere

Input:
    int n: dimensions for quadratic matrix

Output:
    marray*: pointer to the created identity matrix
*/
tarray *tarray_identity(int n, MATRIX_DTYPE dtype);

/*
Creates a copy of the matrix

Input:
    marray *m: pointer to the matrix to copy

Output:
    marray*: pointer to the copy of the input matrix
*/
tarray *tarray_copy(const tarray *m);

/*
Creates the transposed of the given matrix

Input:
    marray *m: matrix from which the transposed is to be created.

Output:
    marray*: pointer to the transposed
*/
tarray *tarray_transposed(const tarray *m);

/*
Creates a matrix from data saved in a file

Input:
    char *filename: file path

Output:
    tarray*: pointer to the matrix containing the data in the file
*/
tarray *tarray_from_file(const char *filename);


/*
Saves the matrix into a file

Input:
    tarray *m: matrix to save
    char *filename: file path

Output:
    bool: result of save operation
*/
bool tarray_to_file(const tarray *m, const char *filename);

/*
Adds the value b to the matrix a

Input:
    marray *a: matrix to which the value b is added
    double b: value that is added to the matrix a

Output:
    None, but the matrix a is the result of the addition
*/
void tarray_addi_val(const tarray *a, double b);

/*
Calculates the result of the matrix a, to which value b is added

Input:
    marray *a: matrix to which the value b is added
    double b: value that is added to the matrix a

Output:
    marray*: pointer to the result of the addition of matrix a and value b
*/
tarray *tarray_add_val(const tarray *a, double b);

/*
Adds the second matrix b to the first matrix a

Input:
    marray *a: matrix to which the matrix b is added
    marray *b: matrix that is added to the matrix a

Output:
    None, but the matrix a is the result of the matrix addition of a and b
*/
void tarray_addi(const tarray *a, const tarray *b);

/*
Calculates the matrix addition of a and b

Input:
    marray *a: matrix to which the matrix b is added
    marray *b: matrix that is added to the matrix a

Output:
    marray*: pointer to the result of the matrix addition of a and b
*/
tarray *tarray_add(const tarray *a, const tarray *b);

/*
Subtracts the second matrix b from the first matrix a

Input:
    marray *a: matrix to which the matrix b is subtracted
    marray *b: matrix that is subtracted from the matrix a

Output:
    None, but the matrix a is the result of the matrix subtraction of a and b
*/
void tarray_subi(const tarray *a, const tarray *b);

/*
Calculates the matrix subtraction of a and b

Input:
    marray *a: matrix to which the matrix b is subtracted
    marray *b: matrix that is subtracted from the matrix a

Output:
    marray*: pointer to the result of the matrix addition of a and b
*/
tarray *tarray_sub(const tarray *a, const tarray *b);

/*
Multiples each element of matrix a with b

Input:
    marray *a: matrix which will be multiplied with b
    double b: value that the matrix a will be multiplied with

Output:
    None, but the matrix a is the result of the multiplication
*/
void tarray_muli_val(const tarray *a, double b);

/*
Multiples each element of matrix a with b

Input:
    marray *a: matrix which will be multiplied with b
    double b: value that the matrix a will be multiplied with

Output:
    marray*: pointer to the matrix containing the result of the multiplication
*/
tarray *tarray_mul_val(const tarray *a, double b);

/*
Creates a partition or sub-area of the given matrix

Input:
    marray **dest: pointer in which the new matrix partition pointer is to be saved
    marray *m: matrix from which a partition is to be created
    row_start: row index to start (included)
    row_end: row index to end (excluded)
    col_start: column index to start (included)
    col_end: column index to end (excluded)
Output:
    None, but the given pointer dest contains the pointer in which the new matrix partition pointer is to be saved
*/
void tarray_get_partition(tarray **dest, const tarray *m, int row_start, int row_end, int col_start, int col_end);

/*
Calculates the dot product of two matrices

Input:
    marray *a: a pointer to a matrix
    marray *b: a pointer to a matrix

Output:
    marray*: pointer to the result of the matrix dot product a*b
*/
tarray *tarray_dot(const tarray *a, const tarray *b);

/*
Checks whether the individual elements of both matrices are close to each other and returns the result per element

Input:
    marray *a: first matrix
    marray *b: second matrix, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    marray*: matrix containing the result of each element
*/
tarray *tarray_close(const tarray *a, const tarray *b, const double rtol, const double atol);

/*
Checks whether the individual elements of both matrices are close to each other

Input:
    marray *a: first matrix
    marray *b: second matrix, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    bool: result of whether the matrices are close to each other
*/
bool tarray_close_all(const tarray *a, const tarray *b, const double rtol, const double atol);


#endif //DTYPE_MATRIX_H
