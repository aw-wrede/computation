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