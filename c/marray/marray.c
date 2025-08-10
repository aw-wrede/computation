#include "marray.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** PRIVATE FUNCTIONS **/

/*
Frees many matrices at once

Input:
    marray *m[]: list of all matrices to free
    int n: amount of elements in the array

Output:
    None
*/
void marray_free_many(marray *m[], const int n) {
    for (int i = 0; i < n; i++) {
        marray_free(m[i]);
    }
}


/*
Performs the dot multiplication with the general method

Input:
    marray *result: pointer to a matrix to store the result
    marray *a: a pointer to a matrix
    marray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
void marray_dot_general(const marray *result, const marray *a, const marray *b) {
    for (int i = 0; i < a->rows; i++) {
        const int row = i * a->cols;

        for (int j = 0; j < b->cols; j++) {
            double val = 0;

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
    marray *a: pointer to the matrix from which the partitions are created
    marray **a11: upper left partition
    marray **a12: upper right partition
    marray **a21: bottom left partition
    marray **a22: bottom right partition

Output:
    None, but the 4 partitions are stored in input pointers
*/
int marray_get_partitions(const marray *a, marray **a11, marray **a12, marray **a21, marray **a22) {
    marray_get_partition(a11, a, 0, a->rows / 2, 0, a->cols / 2);
    marray_get_partition(a12, a, 0, a->rows / 2, a->cols / 2, a->cols);
    marray_get_partition(a21, a, a->rows / 2, a->rows, 0, a->cols / 2);
    marray_get_partition(a22, a, a->rows / 2, a->rows, a->cols / 2, a->cols);

    if (*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
        // free possibly created matrices
        marray_free(*a11);
        marray_free(*a12);
        marray_free(*a21);
        marray_free(*a22);
        return 0;
    }

    return 1;
}


/*
Inserts the given partition into the matrix

Input:
    marray *dest: destination matrix
    marray *m: partition matrix that is applied
    int row: start row
    int col: start col

Output:
    None, but the partition is applied to the target matrix
*/
void marray_apply_partition(const marray *dest, const marray *m, const int row, const int col) {
    const size_t data_size = sizeof(double);

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
    marray *dest: destination matrix
    marray **a11: upper left partition
    marray **a12: upper right partition
    marray **a21: bottom left partition
    marray **a22: bottom right partition

Output:
    None, but the destination matrix is the combination of all 4 partitions
*/
void marray_from_partitions(const marray *dest, const marray *a11, const marray *a12, const marray *a21, const marray *a22) {
    const int n = dest->rows;

    marray_apply_partition(dest, a11, 0, 0);
    marray_apply_partition(dest, a12, 0, n / 2);
    marray_apply_partition(dest, a21, n / 2, 0);
    marray_apply_partition(dest, a22, n / 2, n / 2);
}


/*
Calculates the dot product of two quadratic matrices with even dimension

Input:
    marray *result: pointer to a matrix to store the result
    marray *a: a pointer to a matrix
    marray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
int marray_dot_quadratic(const marray *result, const marray *a, const marray *b) {
    // TODO Improve (dot general is way faster)

    marray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
    marray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

    // check if partitions created successfully
    if (!marray_get_partitions(a, &a11, &a12, &a21, &a22) || !marray_get_partitions(b, &b11, &b12, &b21, &b22)) {
        return 0;
    }

    // intermediate steps
    marray *m1_a = marray_add(a11, a22), *m1_b = marray_add(b11, b22);
    marray *m2_a = marray_add(a21, a22);
    marray *m3_b = marray_sub(b12, b22);
    marray *m4_b = marray_sub(b21, b11);
    marray *m5_a = marray_add(a11, a12);
    marray *m6_a = marray_sub(a21, a11), *m6_b = marray_add(b11, b12);
    marray *m7_a = marray_sub(a12, a22), *m7_b = marray_add(b21, b22);

    // check for successful matrix operations
    if (m1_a == NULL || m1_b == NULL || m2_a == NULL || m3_b == NULL || m4_b == NULL || m5_a == NULL || m6_a == NULL ||
        m6_b == NULL || m7_a == NULL || m7_b == NULL) {
        marray *arrays_1[] = {
            a11, a12, a21, a22, // 4
            b11, b12, b21, b22, // 4
            m1_a, m1_b, m2_a, m3_b, m4_b, m5_a, m6_a, m6_b, m7_a, m7_b}; // 10

        marray_free_many(arrays_1, 18);

        return 0;
    }

    marray *m1 = marray_zeroes(a11->rows, a11->cols);
    marray *m2 = marray_zeroes(a11->rows, a11->cols);
    marray *m3 = marray_zeroes(a11->rows, a11->cols);
    marray *m4 = marray_zeroes(a11->rows, a11->cols);
    marray *m5 = marray_zeroes(a11->rows, a11->cols);
    marray *m6 = marray_zeroes(a11->rows, a11->cols);
    marray *m7 = marray_zeroes(a11->rows, a11->cols);

    marray_dot_general(m1, m1_a, m1_b);
    marray_dot_general(m2, m2_a, b11);
    marray_dot_general(m3, a11, m3_b);
    marray_dot_general(m4, a22, m4_b);
    marray_dot_general(m5, m5_a, b22);
    marray_dot_general(m6, m6_a, m6_b);
    marray_dot_general(m7, m7_a, m7_b);

    // free up memory that is no longer required
    marray *arrays_1[] = {
        a11, a12, a21, a22, // 4
        b11, b12, b21, b22, // 4
        m1_a, m1_b, m2_a, m3_b, m4_b, m5_a, m6_a, m6_b, m7_a, m7_b}; // 10

    marray_free_many(arrays_1, 18);

    // check for successful matrix operations
    if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
        marray *arrays_2[] = {m1, m2, m3, m4, m5, m6, m7};
        marray_free_many(arrays_2, 7);

        return 0;
    }

    // calculate final partitions
    marray *c11 = marray_add(m1, m4);
    marray_subi(c11, m5);
    marray_addi(c11, m7);

    marray *c12 = marray_add(m3, m5);

    marray *c21 = marray_add(m2, m4);

    marray *c22 = marray_sub(m1, m2);
    marray_addi(c22, m3);
    marray_addi(c22, m6);

    // free up memory that is no longer required
    marray *arrays_2[] = {m1, m2, m3, m4, m5, m6, m7};
    marray_free_many(arrays_2, 7);

    // check for successful matrix operations
    if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
        marray *arrays_3[] = {c11, c12, c21, c22};
        marray_free_many(arrays_3, 4);

        return 0;
    }

    marray_from_partitions(result, c11, c12, c21, c22);

    marray *arrays_3[] = {c11, c12, c21, c22};
    marray_free_many(arrays_3, 4);

    return 1;
}


/* ========================================================================================= */


/** PUBLIC HEADER FUNCTIONS **/

void marray_free(marray *m) {
    if (m == NULL) {
        return;
    }

    free(m);
}


/** matrix print functions **/

void marray_print(const marray *m) {
    if (m == NULL) {
        return;
    }

    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            printf("%lf\t", m->data[i * m->cols + j]);
        }

        printf("\n");
    }
}


void marray_print_precision(const marray *m, const int precision) {
    if (m == NULL) {
        return;
    }

    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            printf("%.*lf\t", precision, m->data[i * m->cols + j]);
        }

        printf("\n");
    }
}


/** matrix initialize functions **/

marray *marray_zeroes(const int rows, const int cols) {
    // alloc mem for structure and data
    // calloc sets the entire memory area to 0
    marray *m = calloc(sizeof(marray) + rows * cols * sizeof(double), 1);

    // return null if not enough memory is free
    if (m == NULL) {
        return NULL;
    }

    // apply dimensions
    m->rows = rows;
    m->cols = cols;

    // set data pointer to memory after structure
    m->data = (double *) (m + 1);

    return m;
}


marray *marray_identity(const int n) {
    marray *m = marray_zeroes(n, n);

    if (m == NULL) {
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }

    return m;
}


marray *marray_copy(const marray *m) {
    marray *copy = marray_zeroes(m->rows, m->cols);

    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double));

    return copy;
}


marray *marray_transposed(const marray *m) {
    marray *t = marray_zeroes(m->cols, m->rows);

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


/** matrix calculation functions **/

void marray_addi_val(const marray *a, const double b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b;
    }
}


marray *marray_add_val(const marray *a, const double b) {
    marray *result = marray_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    // copy data to result matrix
    memcpy(result->data, a->data, sizeof(double) * a->rows * a->cols);

    // add b to result matrix
    marray_addi_val(result, b);

    return result;
}


void marray_addi(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b->data[i];
    }
}


marray *marray_add(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    marray *m = marray_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] + b->data[i];
    }

    return m;
}


void marray_subi(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] - b->data[i];
    }
}


marray *marray_sub(const marray *a, const marray *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    marray *m = marray_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] - b->data[i];
    }

    return m;
}


void marray_muli_val(const marray *a, const double b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }
}


marray *marray_mul_val(const marray *a, const double b) {
    marray *result = marray_zeroes(a->rows, a->cols);

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

void marray_get_partition(marray **dest, const marray *m, const int row_start, const int row_end, const int col_start,
                          const int col_end) {
    *dest = marray_zeroes(row_end - row_start, col_end - col_start);

    for (int i = 0; i < row_end - row_start; i++) {
        // dest: i-th row of p
        // src: i + row_start of matrix a, skip rows behind start row
        memcpy(&(*dest)->data[i * (*dest)->cols], &m->data[(i + row_start) * m->cols + col_start],
               sizeof(double) * (*dest)->cols);
    }
}


marray *marray_dot(const marray *a, const marray *b) {
    // TODO Think about (required for linalg back substitution)
    if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
        marray *m = marray_zeroes(1, 1);
        m->data[0] = 0;

        return m;
    }

    // alloc matrix to store result
    marray *m = marray_zeroes(a->rows, b->cols);

    if (m == NULL) {
        return NULL;
    }

    // check if a and b are quadratic
    if (a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
        if (a->cols % 2 == 0) {
            if (marray_dot_quadratic(m, a, b)) {
                return m;
            }

            // dot_quadratic was not successful -> free matrix and return false
            marray_free(m);

            return NULL;
        }
    }

    marray_dot_general(m, a, b);

    return m;
}


/** matrix sum, exp functions **/

double marray_sum(const marray *m) {
    double sum = 0.0;

    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        sum += m->data[i];
    }

    return sum;
}


marray *marray_exp(const marray *m) {
    marray *result = marray_zeroes(m->rows, m->cols);

    if (result == NULL) {
        return NULL;
    }

    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        result->data[i] = exp(m->data[i]);
    }

    return result;
}


/** matrix close functions **/

marray *marray_close(const marray *a, const marray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    // alloc matrix to store result
    marray *results = marray_zeroes(a->rows, a->cols);

    if (results == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    // get result for each entry and store separately
    for (int i = 0; i < elems; i++) {
        results->data[i] = fabs(a->data[i] - b->data[i]) <= atol + rtol * fabs(b->data[i]);
    }

    return results;
}


bool marray_close_all(const marray *a, const marray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }

    const int elems = a->rows * a->cols;

    // check each entry, return false if one mismatches
    for (int i = 0; i < elems; i++) {
        if (! (fabs(a->data[i] - b->data[i]) <= atol + rtol * fabs(b->data[i]))) {
            return false;
        }
    }

    return true;
}