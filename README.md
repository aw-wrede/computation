# Calculation
This repository contains algorithms in Python and C to calculate basic mathematical problems.

## C Documentaton

### mArray

__marray__ is a matrix data structure that stores the individual values as a double. \
It works with the functions in linalg to solve linear systems of equations, for example.

__Current features:__
* Creating zero matrices *(marray_zero)*, unit matrices *(marray_identity)*, copies *(marray_copy)*, transposed matrices *(marray_transposed)*
* Addition *(marray_add_val, marray_addi_val)*, subtraction, multiplication *(marray_mul_val, marray_muli_val)* of matrices with a single value, both as new matrices and on the first one 
* Addition, subtraction, multiplication of matrices, both as new matrices and on the first one (immediate)
* Creating a sub-matrix *(marray_get_partition)*
* Matrix multiplication *(marray_dot)*
* Evaluate whether matrices are close to each other *(marray_close, marray_close_all)*

```C
marray *a = marray_zeroes(4, 4);    // creates a new matrix with zeros and dimension 4x4 (rows x columns)
marray_addi_val(a, 2.1);            // adds 2.1 to every element of the matrix
```

### cArray

__carray__ is a matrix data structure that stores the individual values as a double complex and thus enables complex calculations. \
It offers functions such as the Fourier transform and the inverse.

__Current features:__
* All functions from marray *(carray_...)*
* Conjugate of a matrix *(carray_zero)*
* Adjoint (conjugate transpose) of a matrix *(carray_adjoint)*
* Creating DFT matrix *(carray_dft_matrix)*
* Discrete Fourier Transform *(carray_dft)* and Inverse *(carray_idft)*
* Fast Fourier Transform *(carray_fft)* and Inverse *(carray_ifft)*

```C
carray *signal = marray_zeroes(8, 1);                       // creates a new complex matrix with zeros and dimension 8x1 (rows x columns)

// add your signal to the matrix

carray *freq = cmatrix_fft(signal, COMPLEX_NORM_BACKWARD);  // calculates the frequencies of the signal with the Fast Fourier Transform

// change singal (e.g. frequency filter)

carray *transformed_signal = cmatrix_ifft(freq, COMPLEX_NORM_BACKWARD);

```

### tArray in dtype_matrix

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
\times
\begin{bmatrix} x \\ y \end{bmatrix}
=
\begin{bmatrix} 1 \\ 5 \end{bmatrix}
``` 

```C
marray *a = matrix_zeroes(2, 1);    // matrix for equations
marray *b = matrix_zeroes(2, 1);    // matrix for solutions

marray *x = linalg_solve(a, b);     // values for x=x[0] and y=x[1]

```

## Mathematical background

### Matrix dot product

The dot product of two matrices is defined as follows for
$A \in \mathbb{C}^{m,n}, B \in \mathbb{C}^{n,p}$: \
$A \times B = C \in \mathbb{C}^{n,p}$ with 
$c_{i,j} = \sum_{k=1}^{m} a_{i,k} \cdot b_{k,j}$

This results in the runtime $O(n^3)$.

If both matrices are quadratic with dimension $n \times n$ and $n$ is even 
($n = 2k, k \in \mathbb{N}$),
the matrix multiplication can be simplified as follows 
*(Strassen algorithm)*:
```math
\begin{bmatrix} A_{1,1} & A_{1,2} \\ A_{2,1} & A_{2,2} \end{bmatrix}
\times
\begin{bmatrix} B_{1,1} & B_{1,2} \\ B_{2,1} & B_{2,2} \end{bmatrix}
=
\begin{bmatrix} C_{1,1} & C_{1,2} \\ A_{2,1} & C_{2,2} \end{bmatrix}
```

$$
\begin{align}
C_{1,1} &= A_{1,1} \times B_{1,1} + A_{1,2} \times B_{2,1} \\
C_{1,2} &= A_{1,1} \times B_{1,2} + A_{1,2} \times B_{2,2} \\
C_{2,1} &= A_{2,1} \times B_{1,1} + A_{2,2} \times B_{2,1} \\
C_{2,2} &= A_{2,1} \times B_{1,2} + A_{2,2} \times B_{2,2}
\end{align}
$$

However, this still leads to the following runtime:
$T(n) = 8 \cdot T(n/2) + O(n^2) \Rightarrow T(n) = O(n^3)$

One multiplication can be saved by clever algebraic transformation:

$$
\begin{align}
M_1 &= ( A_{1,1} + A_{2,2} ) \times ( B_{1,1} + B_{2,2} ) \\
M_2 &= ( A_{2,1} + A_{2,2} ) \times B_{1,1} \\
M_3 &= A_{1,1} \times ( B_{1,2} - B_{2,2} ) \\
M_4 &= A_{2,2} \times ( B_{2,1} - B_{1,1} ) \\
M_5 &= ( A_{1,1} + A_{1,2} ) \times B_{2,2} \\
M_6 &= ( A_{2,1} - A_{1,1} ) \times ( B_{1,1} + B_{1,2} ) \\
M_7 &= ( A_{1,2} - A_{2,2} ) \times ( B_{2,1} + B_{2,2} )
\end{align}
$$

This leads to the runtime
$T(n) = 7 \cdot T(n/2) + O(n^2) \Rightarrow T(n) = O(n^{\log_2 7}) \approx O(n^{2.81})$.

This is not the fastest algorithm. The Coppersmith-Winograd algorithm has the runtime
$O(n^{2.376})$. \
*(It's planned to implement better algorithms, also for non-square matrices.)*

### Fourier Transform

#### Unit roots
We are looking for all solutions $z \in \mathbb{C}$ to the equation $z^n = 1$. \
In Euler form this results in $z^n = (r \cdot e^{-i\alpha})^n$
$\Leftrightarrow r=1$ and $\alpha = \frac{2 \pi}{n} \cdot k$ for $k = 0, \ldots, n-1$. \
Solutions: $z_k = e^{i \cdot \frac{2 \pi}{n}  \cdot k}$

The $n$-th unit root is defined as $\omega_n = e^{-i \cdot \frac{2 \pi}{n} }$
and therefore $z_k = {\omega_k}^n$

In addition: $e^{i \cdot \frac{2 \pi}{n}} \cdot k = e^{-i \cdot \frac{2 \pi}{n} \cdot (k \mod n)}$

*To ensure that the functions work like numpy, the unit root is defined with a negative sign in the exponent, so the rotation is clockwise.*

#### DFT matrix
The matrix $\Omega_n$ is a unitary $n \times n$ matrix with $\Omega_n = \frac{1}{\sqrt{n}} \cdot ({\omega_n}^{k \cdot j})_{0 \leq j,k \leq n-1}$ \
The matrix is normalized with $\frac{1}{\sqrt{n}}$.

Because the matrix is symmetrical and unitary, the following applies to the inverse: \
${\Omega_n}^{-1} = \overline{\Omega_n}^T = \overline{\Omega_n}$

```C
// C
#include <complex.h>

const int n = 4; // dimension			
carray *dft_matrix = carray_dft_matrix(n, COMPLEX_NORM_ORTHO); // creates unitary dft matrix
```

```Python
# Python
import complex


n = 4 # dimension
f4 = complex.dft_matrix(n, norm='ortho') # creates unitary dft matrix
```

#### Discrete Fourier Transform
The DFT can be used to analyze the frequencies of a given audio signal. \
For this purpose, the continuous audio signal is discretized (it is sampled at n equidistant sampling points) and the signal $z \in \mathbb{R}^n$ is generated.

Local space to Frequency space: $\Omega_n \times z = \hat{z}$ \
Frequency space to Local space: ${\Omega_n}^{-1} \times \hat{z} = z$

The real parts of the transformed signal represent the cosine oscillations and the imaginary parts the sine oscillations.

The following symmetry exists for real signals: \
$\Re(\hat{z}_k) = \Re(\hat{z}_{n-k})$
for
$1 \leq k \leq \lfloor\frac{n}{2}\rfloor$ \
$\Im(\hat{z}_k) = - \Im(\hat{z}_{n-k})$
for
$1 \leq k \leq \lfloor\frac{n}{2}\rfloor$

This leads to the conclusion that no distinction can be made between frequencies $k$ and $n-k$ and therefore the highest analyzable frequency is in the middle of the frequency spectrum.

However, due to the matrix multiplication, this leads to the runtime $O(n^3)$.

#### Fast Fourier Transform
The general idea of the FFT is to represent $\Omega_{2n}$ with $\Omega_n$
```math
\Omega_{2n} = \frac{1}{\sqrt{2}} \cdot
\begin{bmatrix}
	\Omega_n & F_n \times \Omega_n \\ 
	\Omega_n & - F_n \times \Omega_n
\end{bmatrix}
\times P_{2n}
```
$F_n = diag({\omega_{2n}}^0, {\omega_{2n}}^1, \ldots, {\omega_{2n}}^{n-1})$ 
and $P_{2n}$ is the Permutation matrix *(Bit-Reversal-Matrix)* and permutes the order of the entries so that the decomposition into even and odd indices is mapped correctly.

This leads to an algorithm to calculate $\Omega_n \times z$ with the runtime $O(n \cdot \log n)$

Step 1) Re-sorting the signal *(Bit reversing)*\
Step 2) Step-by-Step combination *(A more detailed explanation may be provided later)*