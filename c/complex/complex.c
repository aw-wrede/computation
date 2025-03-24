#include "complex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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