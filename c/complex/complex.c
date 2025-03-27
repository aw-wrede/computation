#include "complex.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//TODO Think about combining matrix and complex with dtype property in marray

void cmatrix_free(carray *m) {
    if (m == NULL) {
        return;
    }

    free(m);
}

void cmatrix_print(const carray *m) {
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

void cmatrix_print_precision(const carray *m, const int precision) {
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

carray *cmatrix_zeroes(const int rows, const int cols) {
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

carray *cmatrix_identity(const int n) {
    carray *m = cmatrix_zeroes(n, n);

    if (m == NULL) {
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }

    return m;
}

carray *cmatrix_copy(const carray *m) {
    carray *copy = cmatrix_zeroes(m->rows, m->cols);

    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double complex));

    return copy;
}

carray *cmatrix_transposed(const carray *m) {
    carray *t = cmatrix_zeroes(m->cols, m->rows);

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

carray * cmatrix_conj(const carray *m) {
    carray *c = cmatrix_zeroes(m->rows, m->cols);

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

carray *cmatrix_adjoint(const carray *m) {
    carray *t = cmatrix_zeroes(m->cols, m->rows);

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

void cmatrix_addi_val(const carray *a, const double complex b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b;
    }
}

carray * cmatrix_add_val(const carray *a, const double complex b) {
    carray *result = cmatrix_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    // copy data to result matrix
    memcpy(result->data, a->data, sizeof(double complex) * a->rows * a->cols);

    // add b to result matrix
    cmatrix_addi_val(result, b);

    return result;
}

void cmatrix_addi(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b->data[i];
    }
}

carray * cmatrix_add(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    carray *m = cmatrix_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] + b->data[i];
    }

    return m;
}

void cmatrix_subi(const carray *a, const carray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] - b->data[i];
    }
}

carray * cmatrix_sub(const carray *a, const carray *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    carray *m = cmatrix_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] - b->data[i];
    }

    return m;
}

void cmatrix_muli_val(const carray *a, const double complex b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }
}

carray * cmatrix_mul_val(const carray *a, const double complex b) {
    carray *result = cmatrix_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }

    return result;
}

void cdot_general(const carray *result, const carray *a, const carray *b) {
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

void cmatrix_get_partition(carray **dest, const carray *m, const int row_start, const int row_end, const int col_start,
    const int col_end) {
    *dest = cmatrix_zeroes(row_end - row_start, col_end - col_start);

    for (int i = 0; i < row_end - row_start; i++) {
        // dest: i-th row of p
        // src: i + row_start of matrix a, skip rows behind start row
        memcpy(&(*dest)->data[i * (*dest)->cols], &m->data[(i + row_start) * m->cols + col_start],
               sizeof(double complex) * (*dest)->cols);
    }
}

int cget_partitions(const carray *a, carray **a11, carray **a12, carray **a21, carray **a22) {
    cmatrix_get_partition(a11, a, 0, a->rows / 2, 0, a->cols / 2);
    cmatrix_get_partition(a12, a, 0, a->rows / 2, a->cols / 2, a->cols);
    cmatrix_get_partition(a21, a, a->rows / 2, a->rows, 0, a->cols / 2);
    cmatrix_get_partition(a22, a, a->rows / 2, a->rows, a->cols / 2, a->cols);

    if (*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
        // free possibly created matrices
        cmatrix_free(*a11);
        cmatrix_free(*a12);
        cmatrix_free(*a21);
        cmatrix_free(*a22);
        return 0;
    }

    return 1;
}

void capply_partition(const carray *dest, const carray *m, const int row, const int col) {
    for (int i = 0; i < m->rows; i++) {
        const int dest_offset = (i + row) * dest->cols + col;
        const int src_row = i * m->cols;

        for (int j = 0; j < m->cols; j++) {
            dest->data[dest_offset + j] = m->data[src_row + j];
        }
    }
}

void cfrom_partitions(const carray *dest, const carray *a11, const carray *a12, const carray *a21, const carray *a22) {
    const int n = dest->rows;

    capply_partition(dest, a11, 0, 0);
    capply_partition(dest, a12, 0, n / 2);
    capply_partition(dest, a21, n / 2, 0);
    capply_partition(dest, a22, n / 2, n / 2);
}

int cdot_quadratic(const carray *result, const carray *a, const carray *b) {
    //TODO Think about better solution to free matrices if errors occur

    carray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
    carray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

    // check if partitions created successfully
    if (!cget_partitions(a, &a11, &a12, &a21, &a22) || !cget_partitions(b, &b11, &b12, &b21, &b22)) {
        return 0;
    }

    // intermediate steps
    carray *m1_a = cmatrix_add(a11, a22), *m1_b = cmatrix_add(b11, b22);
    carray *m2_a = cmatrix_add(a21, a22);
    carray *m3_b = cmatrix_sub(b12, b22);
    carray *m4_b = cmatrix_sub(b21, b11);
    carray *m5_a = cmatrix_add(a11, a12);
    carray *m6_a = cmatrix_sub(a21, a11), *m6_b = cmatrix_add(b11, b12);
    carray *m7_a = cmatrix_sub(a12, a22), *m7_b = cmatrix_add(b21, b22);

    // check for successful matrix operations
    if (m1_a == NULL || m1_b == NULL || m2_a == NULL || m3_b == NULL || m4_b == NULL || m5_a == NULL || m6_a == NULL ||
        m6_b == NULL || m7_a == NULL || m7_b == NULL) {
        cmatrix_free(a11);
        cmatrix_free(a12);
        cmatrix_free(a21);
        cmatrix_free(a22);

        cmatrix_free(b11);
        cmatrix_free(b12);
        cmatrix_free(b21);
        cmatrix_free(b22);

        cmatrix_free(m1_a);
        cmatrix_free(m1_b);
        cmatrix_free(m2_a);
        cmatrix_free(m3_b);
        cmatrix_free(m4_b);
        cmatrix_free(m5_a);
        cmatrix_free(m6_a);
        cmatrix_free(m6_b);
        cmatrix_free(m7_a);
        cmatrix_free(m7_b);

        return 0;
    }

    carray *m1 = cmatrix_dot(m1_a, m1_b);
    carray *m2 = cmatrix_dot(m2_a, b11);
    carray *m3 = cmatrix_dot(a11, m3_b);
    carray *m4 = cmatrix_dot(a22, m4_b);
    carray *m5 = cmatrix_dot(m5_a, b22);
    carray *m6 = cmatrix_dot(m6_a, m6_b);
    carray *m7 = cmatrix_dot(m7_a, m7_b);

    // free up memory that is no longer required
    cmatrix_free(a11);
    cmatrix_free(a12);
    cmatrix_free(a21);
    cmatrix_free(a22);

    cmatrix_free(b11);
    cmatrix_free(b12);
    cmatrix_free(b21);
    cmatrix_free(b22);

    cmatrix_free(m1_a);
    cmatrix_free(m1_b);
    cmatrix_free(m2_a);
    cmatrix_free(m3_b);
    cmatrix_free(m4_b);
    cmatrix_free(m5_a);
    cmatrix_free(m6_a);
    cmatrix_free(m6_b);
    cmatrix_free(m7_a);
    cmatrix_free(m7_b);

    // check for successful matrix operations
    if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
        cmatrix_free(m1);
        cmatrix_free(m2);
        cmatrix_free(m3);
        cmatrix_free(m4);
        cmatrix_free(m5);
        cmatrix_free(m6);
        cmatrix_free(m7);
        return 0;
    }

    // calculate final partitions
    carray *c11 = cmatrix_add(m1, m4);
    cmatrix_subi(c11, m5);
    cmatrix_addi(c11, m7);

    carray *c12 = cmatrix_add(m3, m5);

    carray *c21 = cmatrix_add(m2, m4);

    carray *c22 = cmatrix_sub(m1, m2);
    cmatrix_addi(c22, m3);
    cmatrix_addi(c22, m6);

    // free up memory that is no longer required
    cmatrix_free(m1);
    cmatrix_free(m2);
    cmatrix_free(m3);
    cmatrix_free(m4);
    cmatrix_free(m5);
    cmatrix_free(m6);
    cmatrix_free(m7);

    // check for successful matrix operations
    if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
        cmatrix_free(c11);
        cmatrix_free(c12);
        cmatrix_free(c21);
        cmatrix_free(c22);
        return 0;
    }

    cfrom_partitions(result, c11, c12, c21, c22);

    // check for successful matrix operations
    cmatrix_free(c11);
    cmatrix_free(c12);
    cmatrix_free(c21);
    cmatrix_free(c22);

    return 1;
}

carray * cmatrix_dot(const carray *a, const carray *b) {
    // TODO Think about (required for linalg back substitution)
    if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
        carray *m = cmatrix_zeroes(1, 1);
        m->data[0] = 0;

        return m;
    }

    // alloc matrix to store result
    carray *m = cmatrix_zeroes(a->rows, b->cols);

    if (m == NULL) {
        return NULL;
    }

    // check if a and b are quadratic
    if (a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
        if (a->cols % 2 == 0) {
            if (cdot_quadratic(m, a, b)) {
                return m;
            }

            // dot_quadratic was not successful -> free matrix and return false
            cmatrix_free(m);

            return NULL;
        }
    }

    cdot_general(m, a, b);

    return m;
}

carray * cmatrix_close(const carray *a, const carray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    // alloc matrix to store result
    carray *results = cmatrix_zeroes(a->rows, a->cols);

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

bool cmatrix_close_all(const carray *a, const carray *b, const double rtol, const double atol) {
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

carray *dft_matrix(const int n, const COMPLEX_NORM norm) {
    carray *f = cmatrix_zeroes(n, n);

    if (f == NULL) {
        return NULL;
    }

    // calculate omega
    const double complex omega = cexp(-1 * I * 2 * M_PI / n);

    // calculate all different unit roots
    carray *omegas = cmatrix_zeroes(1, n);

    if (omegas == NULL) {
        // free f if omegas matrix could not be created
        cmatrix_free(f);
        return NULL;
    }

    // Normalize unit roots with 1/sqrt(n) if norm is ortho
    // Otherwise just calculate unit roots
    if (norm == COMPLEX_NORM_ORTHO) {
        const double n_sqrt = sqrt(n);
        for (int i = 0; i < n; i++) {
            // normalize unit root with 1/sqrt(n)
            omegas->data[i] = cpow(omega, i) / n_sqrt;
        }
    } else {
        for (int i = 0; i < n; i++) {
            omegas->data[i] = cpow(omega, i);
        }
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
    cmatrix_free(omegas);

    return f;
}

bool cmatrix_is_unitary(const carray *m, const double rtol, const double atol) {
    // check if matrix is quadratic
    if (m->rows != m->cols) {
        return false;
    }

    // calculate the adjoint (conjugate transpose) of the matrix to calculate the dot product
    carray *m_adj = cmatrix_adjoint(m);
    carray *dot = cmatrix_dot(m, m_adj);

    cmatrix_free(m_adj);

    carray *unit = cmatrix_identity(m->rows);

    // the product of the matrix with its adjoint gives the unit matrix if the matrix is unitary
    const bool result = cmatrix_close_all(dot, unit, rtol, atol);

    // free temp matrices
    cmatrix_free(dot);
    cmatrix_free(unit);

    return result;
}


/* FOURIER TRANSFORMATION */

carray *cmatrix_dft(const carray *m, const COMPLEX_NORM norm) {
    if (m == NULL) {
        return NULL;
    }

    // compute Omega
    carray *o = dft_matrix(m->rows, norm);
    if (o == NULL) {
        return NULL;
    }

    // perform discrete fourier transform (Omega * m)
    carray *freq = cmatrix_dot(o, m);

    // free temp matrix
    cmatrix_free(o);

    return freq;
}

/*
Mirrors a binary number with fixed length

Input:
    int number: number from which the binary representation is to be mirrored
    int length: fixed number of bits that are set during mirroring

Output:
    int: mirrored number
*/
int bin_mirror(const int number, const int length) {
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
carray *shuffle_bit_reversed_order(const carray *data) {
    // calculate bits required to store the largest index
    // a.e. 0001 -> 1000, filling zeroes are needed for mirroring
    const int bit_length = (int)log2(data->rows);

    carray *shuffled_data = cmatrix_zeroes(data->rows, 1);

    // mirror every index and put entry to this position
    for (int i = 0; i < data->rows; i++) {
        const int i_mirrored = bin_mirror(i, bit_length);
        shuffled_data->data[i_mirrored] = data->data[i];
    }

    return shuffled_data;
}

carray *cmatrix_fft(const carray *data, const COMPLEX_NORM norm) {
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
        return cmatrix_dft(data, norm);
    }

    // first step of fft: shuffle data
    carray *data_shuffled = shuffle_bit_reversed_order(data);

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
    if (norm == COMPLEX_NORM_ORTHO) {
        cmatrix_muli_val(data_shuffled, 1/sqrt(data->rows));
    }

    return data_shuffled;
}