#include "matrix.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void matrix_free(marray *m) {
    if (m == NULL) {
        return;
    }

    free(m);
}

void matrix_print(const marray *m) {
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

marray *matrix_zeroes(const int rows, const int cols) {
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

marray *matrix_identity(const int n) {
    marray *m = matrix_zeroes(n, n);

    if (m == NULL) {
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        m->data[i * n + i] = 1.0;
    }

    return m;
}

marray *matrix_copy(const marray *m) {
    marray *copy = matrix_zeroes(m->rows, m->cols);

    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double));

    return copy;
}

marray *matrix_transposed(const marray *m) {
    marray *t = matrix_zeroes(m->cols, m->rows);

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

void matrix_addi_val(const marray *a, const double b) {
    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b;
    }
}

marray *matrix_add_val(const marray *a, const double b) {
    marray *result = matrix_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    // copy data to result matrix
    memcpy(result->data, a->data, sizeof(double) * a->rows * a->cols);

    // add b to result matrix
    matrix_addi_val(result, b);

    return result;
}

void matrix_addi(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] + b->data[i];
    }
}

marray *matrix_add(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    marray *m = matrix_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] + b->data[i];
    }

    return m;
}

void matrix_subi(const marray *a, const marray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] - b->data[i];
    }
}

marray *matrix_sub(const marray *a, const marray *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    marray *m = matrix_zeroes(a->rows, a->cols);

    if (m == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        m->data[i] = a->data[i] - b->data[i];
    }

    return m;
}

void matrix_muli_val(const marray *a, const double b) {
    for (int i = 0; i < a->rows * a->cols; i++) {
        a->data[i] = a->data[i] * b;
    }
}

marray *matrix_mul_val(const marray *a, const double b) {
    marray *result = matrix_zeroes(a->rows, a->cols);

    if (result == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    for (int i = 0; i < elems; i++) {
        a->data[i] = a->data[i] * b;
    }

    return result;
}

void dot_general(const marray *result, const marray *a, const marray *b) {
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

void matrix_get_partition(marray **dest, const marray *m, const int row_start, const int row_end, const int col_start,
                          const int col_end) {
    *dest = matrix_zeroes(row_end - row_start, col_end - col_start);

    for (int i = 0; i < row_end - row_start; i++) {
        // dest: i-th row of p
        // src: i + row_start of matrix a, skip rows behind start row
        memcpy(&(*dest)->data[i * (*dest)->cols], &m->data[(i + row_start) * m->cols + col_start],
               sizeof(double) * (*dest)->cols);
    }
}

int get_partitions(const marray *a, marray **a11, marray **a12, marray **a21, marray **a22) {
    matrix_get_partition(a11, a, 0, a->rows / 2, 0, a->cols / 2);
    matrix_get_partition(a12, a, 0, a->rows / 2, a->cols / 2, a->cols);
    matrix_get_partition(a21, a, a->rows / 2, a->rows, 0, a->cols / 2);
    matrix_get_partition(a22, a, a->rows / 2, a->rows, a->cols / 2, a->cols);

    if (*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
        return 0;
    }

    return 1;
}

void apply_partition(const marray *dest, const marray *m, const int row, const int col) {
    for (int i = 0; i < m->rows; i++) {
        const int dest_offset = (i + row) * dest->cols + col;
        const int src_row = i * m->cols;

        for (int j = 0; j < m->cols; j++) {
            dest->data[dest_offset + j] = m->data[src_row + j];
        }
    }
}

void from_partitions(const marray *dest, const marray *a11, const marray *a12, const marray *a21, const marray *a22) {
    const int n = dest->rows;

    apply_partition(dest, a11, 0, 0);
    apply_partition(dest, a12, 0, n / 2);
    apply_partition(dest, a21, n / 2, 0);
    apply_partition(dest, a22, n / 2, n / 2);
}

int dot_quadratic(const marray *result, const marray *a, const marray *b) {
    //TODO Think about better solution to free matrices if errors occur

    marray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
    marray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

    // check if partitions created successfully
    if (!get_partitions(a, &a11, &a12, &a21, &a22) || !get_partitions(b, &b11, &b12, &b21, &b22)) {
        return 0;
    }

    // intermediate steps
    marray *m1_a = matrix_add(a11, a22), *m1_b = matrix_add(b11, b22);
    marray *m2_a = matrix_add(a21, a22);
    marray *m3_b = matrix_sub(b12, b22);
    marray *m4_b = matrix_sub(b21, b11);
    marray *m5_a = matrix_add(a11, a12);
    marray *m6_a = matrix_sub(a21, a11), *m6_b = matrix_add(b11, b12);
    marray *m7_a = matrix_sub(a12, a22), *m7_b = matrix_add(b21, b22);

    // check for successful matrix operations
    if (m1_a == NULL || m1_b == NULL || m2_a == NULL || m3_b == NULL || m4_b == NULL || m5_a == NULL || m6_a == NULL ||
        m6_b == NULL || m7_a == NULL || m7_b == NULL) {
        matrix_free(a11);
        matrix_free(a12);
        matrix_free(a21);
        matrix_free(a22);

        matrix_free(b11);
        matrix_free(b12);
        matrix_free(b21);
        matrix_free(b22);

        matrix_free(m1_a);
        matrix_free(m1_b);
        matrix_free(m2_a);
        matrix_free(m3_b);
        matrix_free(m4_b);
        matrix_free(m5_a);
        matrix_free(m6_a);
        matrix_free(m6_b);
        matrix_free(m7_a);
        matrix_free(m7_b);

        return 0;
    }

    marray *m1 = matrix_dot(m1_a, m1_b);
    marray *m2 = matrix_dot(m2_a, b11);
    marray *m3 = matrix_dot(a11, m3_b);
    marray *m4 = matrix_dot(a22, m4_b);
    marray *m5 = matrix_dot(m5_a, b22);
    marray *m6 = matrix_dot(m6_a, m6_b);
    marray *m7 = matrix_dot(m7_a, m7_b);

    // free up memory that is no longer required
    matrix_free(a11);
    matrix_free(a12);
    matrix_free(a21);
    matrix_free(a22);

    matrix_free(b11);
    matrix_free(b12);
    matrix_free(b21);
    matrix_free(b22);

    matrix_free(m1_a);
    matrix_free(m1_b);
    matrix_free(m2_a);
    matrix_free(m3_b);
    matrix_free(m4_b);
    matrix_free(m5_a);
    matrix_free(m6_a);
    matrix_free(m6_b);
    matrix_free(m7_a);
    matrix_free(m7_b);

    // check for successful matrix operations
    if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
        matrix_free(m1);
        matrix_free(m2);
        matrix_free(m3);
        matrix_free(m4);
        matrix_free(m5);
        matrix_free(m6);
        matrix_free(m7);
        return 0;
    }

    // calculate final partitions
    marray *c11 = matrix_add(m1, m4);
    matrix_subi(c11, m5);
    matrix_addi(c11, m7);

    marray *c12 = matrix_add(m3, m5);

    marray *c21 = matrix_add(m2, m4);

    marray *c22 = matrix_sub(m1, m2);
    matrix_addi(c22, m3);
    matrix_addi(c22, m6);

    // free up memory that is no longer required
    matrix_free(m1);
    matrix_free(m2);
    matrix_free(m3);
    matrix_free(m4);
    matrix_free(m5);
    matrix_free(m6);
    matrix_free(m7);

    // check for successful matrix operations
    if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
        matrix_free(c11);
        matrix_free(c12);
        matrix_free(c21);
        matrix_free(c22);
        return 0;
    }

    from_partitions(result, c11, c12, c21, c22);

    // check for successful matrix operations
    matrix_free(c11);
    matrix_free(c12);
    matrix_free(c21);
    matrix_free(c22);

    return 1;
}

marray *matrix_dot(const marray *a, const marray *b) {
    // TODO Think about (required for linalg back substitution)
    if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
        marray *m = matrix_zeroes(1, 1);
        m->data[0] = 0;

        return m;
    }

    // alloc matrix to store result
    marray *m = matrix_zeroes(a->rows, b->cols);

    if (m == NULL) {
        return NULL;
    }

    // check if a and b are quadratic
    if (a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
        if (a->cols % 2 == 0) {
            if (dot_quadratic(m, a, b)) {
                return m;
            }

            // dot_quadratic was not successful -> free matrix and return false
            matrix_free(m);

            return NULL;
        }
    }

    dot_general(m, a, b);

    return m;
}

marray *matrix_close(const marray *a, const marray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    // alloc matrix to store result
    marray *results = matrix_zeroes(a->rows, a->cols);

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

bool matrix_close_all(const marray *a, const marray *b, const double rtol, const double atol) {
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