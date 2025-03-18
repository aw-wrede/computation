#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
  int rows, cols;
  double *data;
} marray;

/*
Frees the memory of the given matrix

Input:
  marray* m: pointer to a matrix

Output:
  None
*/
void matrix_free(marray* m);

/*
Prints the given matrix to console

Input:
  marray* m: pointer to a matrix

Output:
  None, but the matrix is displayed in the console
*/
void matrix_print(const marray* m);

/*
Creates a new matrix with given dimension (rows x cols) with value zero in each element

Input:
  int rows: amount of rows the matrix will have
  int cols: amount of columns the matrix will have

Output:
  marray*: pointer to the created matrix with given dimensions and value zero in each element
*/
marray* matrix_zeroes(int rows, int cols);

/*
Creates a new matrix with one on diagonal elements and zero elsewhere

Input:
  int n: dimensions for quadratic matrix

Output:
  marray*: pointer to the created identity matrix
*/
marray* matrix_identity(const int n);

/*
Creates a copy of the matrix

Input:
  marray* m: pointer to the matrix to copy

Output:
  marray*: pointer to the copy of the input matrix
*/
marray* matrix_copy(const marray* m);

/*
Creates the transposed of the given matrix

Input:
  marray* m: matrix from which the transposed is to be created.

Output:
  marray*: pointer to the transposed
*/
marray* matrix_transposed(const marray* m);

/*
Adds the value b to the matrix a

Input:
  marray* a: matrix to which the matrix b is added
  double b: value that is added to the matrix a

Output:
  None, but the matrix a is the result of the addition
*/
void matrix_addi_val(const marray* a, double b);

/*
Calculates the result of the matrix a, to which value b is added

Input:
  marray* a: matrix to which the value b is added
  double b: value that is added to the matrix a

Output:
  marray*: pointer to the result of the addition of matrix a and value b
*/
marray* matrix_add_val(const marray* a, double b);

/*
Adds the second b matrix to the first matrix a

Input:
  marray* a: matrix to which the matrix b is added
  marray* b: matrix that is added to the matrix a

Output:
  None, but the matrix a is the result of the matrix addition of a and b
*/
void matrix_addi(const marray* a, const marray* b);

/*
Calculates the matrix addition of a and b

Input:
  marray* a: matrix to which the matrix b is added
  marray* b: matrix that is added to the matrix a

Output:
  marray*: pointer to the result of the matrix addition of a and b
*/
marray* matrix_add(const marray* a, const marray* b);

/*
Subtracts the second b matrix to the first matrix a

Input:
  marray* a: matrix to which the matrix b is subtracted
  marray* b: matrix that is subtracted from the matrix a

Output:
  None, but the matrix a is the result of the matrix subtraction of a and b
*/
void matrix_subi(const marray* a, const marray* b);

/*
Calculates the matrix subtraction of a and b

Input:
  marray* a: matrix to which the matrix b is subtracted
  marray* b: matrix that is subtracted from the matrix a

Output:
  marray*: pointer to the result of the matrix addition of a and b
*/
marray* matrix_sub(const marray* a, const marray* b);

/*
Calculates the dot product of two matrices

Input:
  marray* a: a pointer to a matrix
  marray* b: a pointer to a matrix

Output:
  marray*: pointer to the result of the matrix dot product a*b
*/
marray* matrix_dot(const marray* a, const marray* b);

#endif //MATRIX_H
