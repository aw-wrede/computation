#include "complex.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** PRIVATE FUNCTIONS **/

/*
Performs the dot multiplication with the general method

Input:
    carray *result: pointer to a matrix to store the result
    carray *a: a pointer to a matrix
    carray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
void carray_dot_general(const carray *result, const carray *a, const carray *b) {
    for (int i = 0; i < a->rows; i++) {
        const int row = i * a->cols;

        for (int j = 0; j < b->cols; j++) {
            double complex val = 0;

            for (int k = 0; k < a->cols; k++) {
                val += a->data[row + k] * b->data[k * b->cols + j];
            }

            result->data[i * result->cols + j] = val;
        }
    }
}


/*
Splits the matrix into 4 partitions

Input:
    carray *a: pointer to the matrix from which the partitions are created
    carray **a11: upper left partition
    carray **a12: upper right partition
    carray **a21: bottom left partition
    carray **a22: bottom right partition

Output:
    None, but the 4 partitions are stored in input pointers
*/
int carray_get_partitions(const carray *a, carray **a11, carray **a12, carray **a21, carray **a22) {
    carray_get_partition(a11, a, 0, a->rows / 2, 0, a->cols / 2);
    carray_get_partition(a12, a, 0, a->rows / 2, a->cols / 2, a->cols);
    carray_get_partition(a21, a, a->rows / 2, a->rows, 0, a->cols / 2);
    carray_get_partition(a22, a, a->rows / 2, a->rows, a->cols / 2, a->cols);

    if (*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
        // free possibly created matrices
        carray_free(*a11);
        carray_free(*a12);
        carray_free(*a21);
        carray_free(*a22);
        return 0;
    }

    return 1;
}


/*
Inserts the given partition into the matrix

Input:
    carray *dest: destination matrix
    carray *m: partition matrix that is applied
    int row: start row
    int col: start col

Output:
    None, but the partition is applied to the target matrix
*/
void carray_apply_partition(const carray *dest, const carray *m, const int row, const int col) {
    const size_t data_size = sizeof(double complex);

    for (int i = 0; i < m->rows; i++) {
        const int dest_offset = (i + row) * dest->cols + col;
        const int src_row = i * m->cols;

        memcpy(dest->data + dest_offset,
            m->data + src_row,
            m->cols * data_size);
    }
}


/*
Inserts all 4 partitions into the destination matrix

Input:
    carray *dest: destination matrix
    carray **a11: upper left partition
    carray **a12: upper right partition
    carray **a21: bottom left partition
    carray **a22: bottom right partition

Output:
    None, but the destination matrix is the combination of all 4 partitions
*/
void carray_from_partitions(const carray *dest, const carray *a11, const carray *a12, const carray *a21, const carray *a22) {
    const int n = dest->rows;

    carray_apply_partition(dest, a11, 0, 0);
    carray_apply_partition(dest, a12, 0, n / 2);
    carray_apply_partition(dest, a21, n / 2, 0);
    carray_apply_partition(dest, a22, n / 2, n / 2);
}


/*
Calculates the dot product of two quadratic matrices with even dimension

Input:
    carray *result: pointer to a matrix to store the result
    carray *a: a pointer to a matrix
    carray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
int carray_dot_quadratic(const carray *result, const carray *a, const carray *b) {
    //TODO Think about better solution to free matrices if errors occur

    carray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
    carray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

    // check if partitions created successfully
    if (!carray_get_partitions(a, &a11, &a12, &a21, &a22) || !carray_get_partitions(b, &b11, &b12, &b21, &b22)) {
        return 0;
    }

    // intermediate steps
    carray *m1_a = carray_add(a11, a22), *m1_b = carray_add(b11, b22);
    carray *m2_a = carray_add(a21, a22);
    carray *m3_b = carray_sub(b12, b22);
    carray *m4_b = carray_sub(b21, b11);
    carray *m5_a = carray_add(a11, a12);
    carray *m6_a = carray_sub(a21, a11), *m6_b = carray_add(b11, b12);
    carray *m7_a = carray_sub(a12, a22), *m7_b = carray_add(b21, b22);

    // check for successful matrix operations
    if (m1_a == NULL || m1_b == NULL || m2_a == NULL || m3_b == NULL || m4_b == NULL || m5_a == NULL || m6_a == NULL ||
        m6_b == NULL || m7_a == NULL || m7_b == NULL) {
        carray_free(a11);
        carray_free(a12);
        carray_free(a21);
        carray_free(a22);

        carray_free(b11);
        carray_free(b12);
        carray_free(b21);
        carray_free(b22);

        carray_free(m1_a);
        carray_free(m1_b);
        carray_free(m2_a);
        carray_free(m3_b);
        carray_free(m4_b);
        carray_free(m5_a);
        carray_free(m6_a);
        carray_free(m6_b);
        carray_free(m7_a);
        carray_free(m7_b);

        return 0;
    }

    carray *m1 = carray_dot(m1_a, m1_b);
    carray *m2 = carray_dot(m2_a, b11);
    carray *m3 = carray_dot(a11, m3_b);
    carray *m4 = carray_dot(a22, m4_b);
    carray *m5 = carray_dot(m5_a, b22);
    carray *m6 = carray_dot(m6_a, m6_b);
    carray *m7 = carray_dot(m7_a, m7_b);

    // free up memory that is no longer required
    carray_free(a11);
    carray_free(a12);
    carray_free(a21);
    carray_free(a22);

    carray_free(b11);
    carray_free(b12);
    carray_free(b21);
    carray_free(b22);

    carray_free(m1_a);
    carray_free(m1_b);
    carray_free(m2_a);
    carray_free(m3_b);
    carray_free(m4_b);
    carray_free(m5_a);
    carray_free(m6_a);
    carray_free(m6_b);
    carray_free(m7_a);
    carray_free(m7_b);

    // check for successful matrix operations
    if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
        carray_free(m1);
        carray_free(m2);
        carray_free(m3);
        carray_free(m4);
        carray_free(m5);
        carray_free(m6);
        carray_free(m7);
        return 0;
    }

    // calculate final partitions
    carray *c11 = carray_add(m1, m4);
    carray_subi(c11, m5);
    carray_addi(c11, m7);

    carray *c12 = carray_add(m3, m5);

    carray *c21 = carray_add(m2, m4);

    carray *c22 = carray_sub(m1, m2);
    carray_addi(c22, m3);
    carray_addi(c22, m6);

    // free up memory that is no longer required
    carray_free(m1);
    carray_free(m2);
    carray_free(m3);
    carray_free(m4);
    carray_free(m5);
    carray_free(m6);
    carray_free(m7);

    // check for successful matrix operations
    if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
        carray_free(c11);
        carray_free(c12);
        carray_free(c21);
        carray_free(c22);
        return 0;
    }

    carray_from_partitions(result, c11, c12, c21, c22);

    carray_free(c11);
    carray_free(c12);
    carray_free(c21);
    carray_free(c22);

    return 1;
}


/** matrix fft functions **/

/*
Mirrors a binary number with fixed length

Input:
    int number: number from which the binary representation is to be mirrored
    int length: fixed number of bits that are set during mirroring

Output:
    int: mirrored number
*/
int carray_bin_mirror(const int number, const int length) {
    if (number == 0) {
        return 0;
    }

    int binary[length];

    // get binary representation of number (inverted)
    for (int i = 0; i < length; i++) {
        binary[i] = (number >> i) & 1;
    }

    int number_mirrored = 0;
    for (int i = 0; i < length; i++) {
        number_mirrored = (number_mirrored << 1) | binary[i];
    }

    return number_mirrored;
}


/*
Shuffles the elements of data using bit-reversal of list index.

Input:
    carray *data: data to be shuffled

Output:
    carray: Shuffled data array
*/
carray *carray_shuffle_bit_reversed_order(const carray *data) {
    // calculate bits required to store the largest index
    // a.e. 0001 -> 1000, filling zeroes are needed for mirroring
    const int bit_length = (int)log2(data->rows);

    carray *shuffled_data = carray_zeroes(data->rows, 1);

    // mirror every index and put entry to this position
    for (int i = 0; i < data->rows; i++) {
        const int i_mirrored = carray_bin_mirror(i, bit_length);
        shuffled_data->data[i_mirrored] = data->data[i];
    }

    return shuffled_data;
}


/* ========================================================================================= */


/** PUBLIC HEADER FUNCTIONS **/

void carray_free(carray *m) {
    if (m == NULL) {
        return;
    }

    free(m);
}


/** matrix print functions **/

void carray_print(const carray *m) {
    if (m == NULL) {
        return;
    }

    for (int i = 0; i < m->rows; i++) {
        const int row = i * m->cols;

        for (int j = 0; j < m->cols; j++) {
            const double complex z = m->data[row + j];
            printf("%lf %+lfi\t", creal(z), cimag(z));
        }

        printf("\n");
    }
}


void carray_print_precision(const carray *m, const int precision) {
    if (m == NULL) {
        return;
    }

    for (int i = 0; i < m->rows; i++) {
        const int row = i * m->cols;

        for (int j = 0; j < m->cols; j++) {
            const double complex z = m->data[row + j];
            printf("%.*lf %+.*lfi\t", precision, creal(z), precision, cimag(z));
        }

        printf("\n");
    }
}


/** matrix initialize functions **/

carray *carray_zeroes(const int rows, const int cols) {
    // alloc mem for structure and data
    // calloc sets the entire memory area to 0
    carray *m = calloc(sizeof(carray) + rows * cols * sizeof(double complex), 1);

    // return null if not enough memory is free
    if (m == NULL) {
        return NULL;
    }

    // apply dimensions
    m->rows = rows;
    m->cols = cols;

    // set data pointer to memory after structure
    m->data = (double complex *) (m + 1);

    return m;
}


carray *carray_identity(const int n) {
    carray *m = carray_zeroes(n, n);

    if (m == NULL) {
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }

    return m;
}


carray *carray_copy(const carray *m) {
    carray *copy = carray_zeroes(m->rows, m->cols);

    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double complex));

    return copy;
}


carray *carray_transposed(const carray *m) {
    carray *t = carray_zeroes(m->cols, m->rows);

    if (t == NULL) {
        return NULL;
    }

    for (int i = 0; i < t->rows; i++) {
        const int row = i * t->cols;

        for (int j = 0; j < t->cols; j++) {
            t->data[row + j] = m->data[j * m->cols + i];
        }
    }

    return t;
}


carray * carray_conj(const carray *m) {
    carray *c = carray_zeroes(m->rows, m->cols);

    if (c == NULL) {
        return NULL;
    }

    for (int i = 0; i < c->rows; i++) {
        const int row = i * c->cols;

        for (int j = 0; j < c->cols; j++) {
            c->data[row + j] = conj(m->data[row + j]);
        }
    }

    return c;
}


carray *carray_adjoint(const carray *m) {
    carray *t = carray_zeroes(m->cols, m->rows);

    if (t == NULL) {
        return NULL;
    }

    for (int i = 0; i < t->rows; i++) {
        const int row = i * t->cols;

        for (int j = 0; j < t->cols; j++) {
            t->data[row + j] = conj(m->data[j * m->cols + i]);
        }
    }

    return t;
}


/** matrix calculation functions **/

void carray_addi_val(const carray *a, const double complex b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b;
    }
}

carray * carray_add_val(const carray *a, const double complex b) {
    carray *result = carray_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    // copy data to result matrix
    memcpy(result->data, a->data, sizeof(double complex) * a->rows * a->cols);

    // add b to result matrix
    carray_addi_val(result, b);

    return result;
}

void carray_addi(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b->data[i];
    }
}

carray * carray_add(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    carray *m = carray_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] + b->data[i];
    }

    return m;
}

void carray_subi(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] - b->data[i];
    }
}

carray * carray_sub(const carray *a, const carray *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    carray *m = carray_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] - b->data[i];
    }

    return m;
}

void carray_muli_val(const carray *a, const double complex b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }
}

carray * carray_mul_val(const carray *a, const double complex b) {
    carray *result = carray_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }

    return result;
}


/** matrix dot functions **/

void carray_get_partition(carray **dest, const carray *m, const int row_start, const int row_end, const int col_start,
    const int col_end) {
    *dest = carray_zeroes(row_end - row_start, col_end - col_start);

    for (int i = 0; i < row_end - row_start; i++) {
        // dest: i-th row of p
        // src: i + row_start of matrix a, skip rows behind start row
        memcpy(&(*dest)->data[i * (*dest)->cols], &m->data[(i + row_start) * m->cols + col_start],
               sizeof(double complex) * (*dest)->cols);
    }
}


carray * carray_dot(const carray *a, const carray *b) {
    // TODO Think about (required for linalg back substitution)
    if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
        carray *m = carray_zeroes(1, 1);
        m->data[0] = 0;

        return m;
    }

    // alloc matrix to store result
    carray *m = carray_zeroes(a->rows, b->cols);

    if (m == NULL) {
        return NULL;
    }

    // check if a and b are quadratic
    if (a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
        if (a->cols % 2 == 0) {
            if (carray_dot_quadratic(m, a, b)) {
                return m;
            }

            // dot_quadratic was not successful -> free matrix and return false
            carray_free(m);

            return NULL;
        }
    }

    carray_dot_general(m, a, b);

    return m;
}


/** matrix close functions **/

carray * carray_close(const carray *a, const carray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    // alloc matrix to store result
    carray *results = carray_zeroes(a->rows, a->cols);

    if (results == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    // get result for each entry and store separately
    for (int i = 0; i < elems; i++) {
        results->data[i] = cabs(a->data[i] - b->data[i]) <= atol + rtol * cabs(b->data[i]);
    }

    return results;
}


bool carray_close_all(const carray *a, const carray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }

    const int elems = a->rows * a->cols;

    // check each entry, return false if one mismatches
    for (int i = 0; i < elems; i++) {
        if (! (cabs(a->data[i] - b->data[i]) <= atol + rtol * cabs(b->data[i]))) {
            return false;
        }
    }

    return true;
}


/** matrix fft functions **/

carray *carray_dft_matrix(const int n, const COMPLEX_NORM norm) {
    carray *f = carray_zeroes(n, n);

    if (f == NULL) {
        return NULL;
    }

    // calculate omega
    const double complex omega = cexp(-1 * I * 2 * M_PI / n);

    // calculate all different unit roots
    carray *omegas = carray_zeroes(1, n);

    if (omegas == NULL) {
        // free f if omegas matrix could not be created
        carray_free(f);
        return NULL;
    }

    // normalize based on selected mode
    switch (norm) {
        // normalize unit roots with 1/sqrt(n) if norm is ortho
        case COMPLEX_NORM_ORTHO:
            const double n_sqrt = sqrt(n);
            for (int i = 0; i < n; i++) {
                // normalize unit root with 1/sqrt(n)
                omegas->data[i] = cpow(omega, i) / n_sqrt;
            }
            break;

        // normalize unit roots with 1/n if norm is forward
        case COMPLEX_NORM_FORWARD:
            for (int i = 0; i < n; i++) {
                // normalize unit root with 1/n
                omegas->data[i] = cpow(omega, i) / n;
            }
            break;

        // otherwise just calculate unit roots
        default:
            for (int i = 0; i < n; i++) {
                omegas->data[i] = cpow(omega, i);
            }
            break;
    }

    // fill matrix with unit roots
    for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
            // there are only n different unit roots -> calculate required one
            const double complex result = omegas->data[k * l % n];

            // dft matrix is symmetric -> set symmetry pairs
            f->data[k * f->cols + l] = result;
            f->data[l * f->cols + k] = result;
        }
    }

    // free temp matrix
    carray_free(omegas);

    return f;
}


bool carray_is_unitary(const carray *m, const double rtol, const double atol) {
    // check if matrix is quadratic
    if (m->rows != m->cols) {
        return false;
    }

    // calculate the adjoint (conjugate transpose) of the matrix to calculate the dot product
    carray *m_adj = carray_adjoint(m);
    carray *dot = carray_dot(m, m_adj);

    carray_free(m_adj);

    carray *unit = carray_identity(m->rows);

    // the product of the matrix with its adjoint gives the unit matrix if the matrix is unitary
    const bool result = carray_close_all(dot, unit, rtol, atol);

    // free temp matrices
    carray_free(dot);
    carray_free(unit);

    return result;
}


carray *carray_dft(const carray *m, const COMPLEX_NORM norm) {
    if (m == NULL) {
        return NULL;
    }

    // compute Omega
    carray *o = carray_dft_matrix(m->rows, norm);
    if (o == NULL) {
        return NULL;
    }

    // perform discrete fourier transform (Omega * m)
    carray *freq = carray_dot(o, m);

    // free temp matrix
    carray_free(o);

    return freq;
}


carray *carray_idft(const carray *m, const COMPLEX_NORM norm) {
    if (m == NULL) {
        return NULL;
    }

    carray *m_conj = carray_conj(m);
    carray *result = NULL;

    if (m_conj == NULL) {
        return NULL;
    }

    switch (norm) {
        case COMPLEX_NORM_ORTHO:
            result = carray_dft(m_conj, COMPLEX_NORM_ORTHO);
            break;
        case COMPLEX_NORM_FORWARD:
            result = carray_dft(m_conj, COMPLEX_NORM_BACKWARD);
            break;
        default:
            result = carray_dft(m_conj, COMPLEX_NORM_FORWARD);
            break;
    }

    carray_free(m_conj);

    if (result == NULL) {
        return NULL;
    }

    // allows complex input on dft/fft
    carray *result_conj = carray_conj(result);

    carray_free(result);

    return result_conj;
}


carray *carray_fft(const carray *data, const COMPLEX_NORM norm) {
    if (data == NULL) {
        return NULL;
    }

    // check if data is empty
    if (data->rows == 0) {
        return NULL;
    }

    // check if data is one-dimensional
    if (data->cols != 1) {
        return NULL;
    }

    // check if input length is power of two
    // otherwise perform slow discrete fourier transform
    if (data->rows & (data->cols - 1) != 0) {
        return carray_dft(data, norm);
    }

    // first step of fft: shuffle data
    carray *data_shuffled = carray_shuffle_bit_reversed_order(data);

    // second step: iteratively merge transforms
    const int steps = (int) log2(data->rows);

    for (int step = 0; step < steps; step++) {
        // k = 2^(step+1)
        const int k = 1 << (step + 1);

        // unit square root
        const double complex omega = cexp(-1 * I * 2 * M_PI / k);

        // for each block
        // each step has 2^(steps-step-1) blocks
        const int blocks = 1 << (steps - step - 1);
        for (int i = 0; i < blocks ; i++) {
            // offset or index of first element in block: i * 2^(step+1)
            const int block_start = i * (1 << (step + 1));

            // for each pair of elements in a block
            // every block has 2^(step) elements
            const int elems = 1 << step;
            for (int j = 0; j < elems; j++) {
                // index of first element to calculate
                const int j_0 = block_start + j;
                // index of second element to calculate
                const int j_1 = block_start + j + (1 << step);

                const double complex s_0 = data_shuffled->data[block_start + j] +
                        cpow(omega, j) * data_shuffled->data[(1 << step) + block_start + j];

                const double complex s_1 = data_shuffled->data[block_start + j] +
                    cpow(omega, j) * (-1) * data_shuffled->data[(1 << step) + block_start + j];

                data_shuffled->data[j_0] = s_0;
                data_shuffled->data[j_1] = s_1;
            }

        }
    }

    // normalize fft signal with 1/sqrt(n) if norm is ortho
    // normalize fft signal with 1/n if norm is forward
    // otherwise return unmodified
    if (norm == COMPLEX_NORM_ORTHO) {
        carray_muli_val(data_shuffled, 1/sqrt(data->rows));
    } else if (norm == COMPLEX_NORM_FORWARD) {
        carray_muli_val(data_shuffled, (double) 1/data->rows);
    }

    return data_shuffled;
}


carray *carray_ifft(const carray *data, const COMPLEX_NORM norm) {
    if (data == NULL) {
        return NULL;
    }

    carray *data_conj = carray_conj(data);
    carray *result = NULL;

    if (data_conj == NULL) {
        return NULL;
    }

    switch (norm) {
        case COMPLEX_NORM_ORTHO:
            result = carray_fft(data_conj, COMPLEX_NORM_ORTHO);
            break;
        case COMPLEX_NORM_FORWARD:
            result = carray_fft(data_conj, COMPLEX_NORM_BACKWARD);
            break;
        default:
            result = carray_fft(data_conj, COMPLEX_NORM_FORWARD);
            break;
    }

    carray_free(data_conj);

    if (result == NULL) {
        return NULL;
    }

    // allows complex input on dft/fft
    carray *result_conj = carray_conj(result);

    carray_free(result);

    return result_conj;
}