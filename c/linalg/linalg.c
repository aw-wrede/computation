#include "linalg.h"
#include <stdio.h>
#include <string.h>

int find_largest_pivot_element(const marray* m, const int col) {

  int pivot_index = 0;
  const double* pivot = NULL;

  // go from start row till end to find the largest pivot element
  for (int i = col; i < m->rows; i++) {

    // if no pivot element exists or
    // current is larger than the old one and current is not zero
    // (prevent unsolvable if largest element would be zero because other candidates are negative)
    if (pivot == NULL || (m->data[i * m->cols + col] > *pivot && m->data[i * m->cols + col] != 0)) {
      pivot = &m->data[i * m->cols + col];
      pivot_index = i;
    }

  }

  return pivot_index;
}

bool linalg_gaussian_elimination(const marray *a, const marray *b) {
  for (int i = 0; i < a->rows; i++) {

    const int max_row = find_largest_pivot_element(a, i);
    // return false if no pivot element is non-zero
    if (a->data[i * a->cols + max_row] == 0.0) { //TODO Discuss use of close function instead
      return 0;
    }

    // swap rows of matrix a
    double origin[a->cols];

    // save original row i in temp
    memcpy(origin, &a->data[i * a->cols], a->cols * sizeof(double));

    // copy new pivot row (max_row) to i-th row
    memcpy(&a->data[i * a->cols], &a->data[max_row * a->cols], a->cols * sizeof(double));

    // copy temp (old i-th) row to max_row
    memcpy(&a->data[max_row * a->cols], origin, a->cols * sizeof(double));

    // swap rows of matrix b
    const double temp = b->data[i];
    b->data[i] = b->data[max_row];
    b->data[max_row] = temp;

    const double pivot = a->data[i * a->cols + i];

    // eliminate all fields below pivot element
    for (int row = i + 1; row < a->rows; row++) {
      const double factor = a->data[row * a->cols + i] / pivot;

      for (int col = i; col < a->cols; col++) {
        a->data[row * a->cols + col] -= a->data[i * a->cols + col] * factor;
      }

      b->data[row] -= factor * b->data[i];

    }

  }

  return 1;
}

marray* linalg_solve(const marray* a, const marray* b) {

  if (a == NULL || b == NULL) {
    return NULL;
  }

  // check dimensions
  // only allow well-determined matrix a (full rank)
  if (a->rows != a->cols || a->rows != b->rows || b->cols != 1) {
    return NULL;
  }

  //TODO

  return NULL;
}