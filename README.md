# Calculation
This repository contains algorithms in Python and C to calculate basic mathematical problems.

## C Documentaton

### Marray

__marray__ is a matrix data structure that stores the individual values as a double. \
It works with the functions in linalg to solve linear systems of equations, for example.

__Current features:__
* Creating zero matrices *(marray_zero)*, unit matrices *(marray_identity)*, copies *(marray_copy)*, transposed matrices *(marray_transposed)*
* Addition *(marray_add_val, marray_addi_val)*, subtraction, multiplication *(marray_mul_val, marray_muli_val)* of matrices with a single value, both as new matrices and on the first one 
* Addition, subtraction, multiplication of matrices, both as new matrices and on the first one 
* Creating a sub-matrix *(marray_get_partition)*
* Matrix multiplication *(marray_dot)*
* Evaluate whether matrices are close to each other *(marray_close, marray_close_all)*

```
marray *a = marray_zeroes(4, 4);    // creates a new matrix with zeros and dimension 4x4 (rows x columns)
marray_addi_val(a, 2.1);            // adds 2.1 to every element of the matrix
```

### Carray

__carray__ is a matrix data structure that stores the individual values as a double complex and thus enables complex calculations. \
It offers functions such as the Fourier transform and the inverse.

__Current features:__
* All functions from marray *(carray_...)*
* Conjugate of a matrix *(carray_zero)*
* Adjoint (conjugate transpose) of a matrix *(carray_adjoint)*
* DFT matrix *(carray_dft_matrix)*
* Discrete Fourier Transform *(carray_dft)* and Inverse *(carray_idft)*
* Fast Fourier Transform *(carray_fft)* and Inverse *(carray_ifft)*

```
carray *signal = marray_zeroes(8, 1);                       // creates a new complex matrix with zeros and dimension 8x1 (rows x columns)
carray *freq = cmatrix_fft(signal, COMPLEX_NORM_BACKWARD);  // calculates the frequencies of the signal with the Fast Fourier Transform
```

### Tarray in dtype_matrix

__tarray__ is a matrix data structure that supports different data types, real and complex. \
_However, it is uncertain to what extent this structure will be further developed in this repository_


### linalg

linalg works with marray matrices to solve linear algebra problems.

__Current features:__
* Solving systems of linear equations (with full rank) *(linalg_solve)*
* Gaussian elimination *(linalg_gaussian_elimination)*
* Back substitution *(linalg_back_substitution)*

#### Example
$$
\begin{align}
I:  & 2x &- y &= 1 \\
II: & x &+ y &= 5
\end{align}
$$

```math
\begin{bmatrix} 2 & -1 \\ 1 & 1 \end{bmatrix}
\cdot
\begin{bmatrix} x \\ y \end{bmatrix}
=
\begin{bmatrix} 1 \\ 5 \end{bmatrix}
```

```
marray *a = matrix_zeroes(2, 1);    // matrix for equations
marray *b = matrix_zeroes(2, 1);    // matrix for solutions

marray *x = linalg_solve(a, b);     // values for x=x[0] and y=x[1]

```
