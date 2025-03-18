#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void matrix_free(marray* m) {
  if (m == NULL) {
    return;
  }

  if(m->data != NULL) {
    free(m->data);
  }

  free(m);
}

void matrix_print(const marray* m) {
  if (m == NULL) {
    return;
  }

  for(int i = 0; i < m->rows; i++) {
    for(int j = 0; j < m->cols; j++) {
      printf("%lf ", m->data[i * m->cols + j]);
    }
    printf("\n");
  }
}

marray* matrix_zeroes(const int rows, const int cols) {

  // alloc mem for structure
  marray* m = malloc(sizeof(marray));

  // return null if not enough memory is free
  if(m == NULL) {
    return NULL;
  }

  // apply dimensions
  m->rows = rows;
  m->cols = cols;

  // alloc memory for array
  m->data = (double*) calloc(m->rows * m->cols, sizeof(double));

  // return null if not enough memory is free
  if(m->data == NULL) {
    free(m);
    return NULL;
  }

  return m;
}

marray* matrix_identity(const int n) {
  marray* m = matrix_zeroes(n, n);

  if(m == NULL) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    m->data[i*n + i] = 1.0;
  }

  return m;
}

marray* matrix_copy(const marray* m) {

  marray* copy = matrix_zeroes(m->rows, m->cols);

  if(copy == NULL) {
    return NULL;
  }

  memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double));

  return copy;
}

marray* matrix_transposed(const marray* m) {
  marray* t = matrix_zeroes(m->cols, m->rows);

  if(t == NULL) {
    return NULL;
  }

  for (int i = 0; i < t->rows; i++) {
    for (int j = 0; j < t->cols; j++) {
      t->data[i * t->cols + j] = m->data[j * m->cols + i];
    }
  }

  return t;
}

void matrix_addi_val(const marray* a, const double b) {
  for(int i = 0; i < a->rows * a->cols; i++) {
    a->data[i] = a->data[i] + b;
  }
}

marray* matrix_add_val(const marray* a, const double b) {
  marray* result = matrix_zeroes(a->rows, a->cols);

  if(result == NULL) {
    return NULL;
  }

  // copy data to result matrix
  memcpy(result->data, a->data, sizeof(double) * a->rows * a->cols);

  // add b to result matrix
  matrix_addi_val(result, b);

  return result;
}

void matrix_addi(const marray* a, const marray* b) {
  for(int i = 0; i < a->rows * a->cols; i++) {
    a->data[i] = a->data[i] + b->data[i];
  }
}

marray* matrix_add(const marray* a, const marray* b) {

  // Check dimensions
  if(a->rows != b->rows || a->cols != b->cols) {
    return NULL;
  }

  marray* m = matrix_zeroes(a->rows, a->cols);

  if(m == NULL) {
    return NULL;
  }

  for(int i = 0; i < a->rows * a->cols; i++) {
    m->data[i] = a->data[i] + b->data[i];
  }

  return m;
}

void matrix_subi(const marray* a, const marray* b) {

  // Check dimensions
  if(a->rows != b->rows || a->cols != b->cols) {
    return;
  }

  for(int i = 0; i < a->rows * a->cols; i++) {
    a->data[i] = a->data[i] - b->data[i];
  }
}

marray* matrix_sub(const marray* a, const marray* b) {

  if(a->rows != b->rows || a->cols != b->cols) {
    return NULL;
  }

  marray* m = matrix_zeroes(a->rows, a->cols);
  if(m == NULL) {
    return NULL;
  }

  for(int i = 0; i < a->rows * a->cols; i++) {
    m->data[i] = a->data[i] - b->data[i];
  }

  return m;
}

void dot_general(const marray* result, const marray* a, const marray* b) {
  for(int i = 0; i < a->rows; i++) {
    for(int j = 0; j < b->cols; j++) {
      double val = 0;
      for(int k = 0; k < a->rows; k++) {
        val += a->data[i * a->cols + k] * b->data[k * b->cols + j];
      }

      result->data[i * result->rows + j] = val;
    }
  }
}

void get_partition(marray** p, const marray* a,  const int rowStart, const int rowEnd, const int colStart, const int colEnd) {
  *p = matrix_zeroes(rowEnd - rowStart, colEnd - colStart);
  for(int i = 0; i < rowEnd - rowStart; i++) {
    // dest: i-th row of p
    // src: i + row_start of matrix a, skip rows behind start row
    memcpy( &(*p)->data[i * (*p)->cols], &a->data[(i + rowStart) * a->cols + colStart], sizeof(double) * (*p)->cols );
  }
}

int get_partitions(const marray* a, marray** a11, marray** a12, marray** a21, marray** a22) {
  get_partition(a11, a, 0,         a->rows/2, 0,         a->cols/2);
  get_partition(a12, a, 0,         a->rows/2, a->cols/2, a->cols);
  get_partition(a21, a, a->rows/2, a->rows,   0,         a->cols/2);
  get_partition(a22, a, a->rows/2, a->rows,   a->cols/2, a->cols);

  if(*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
    return 0;
  }

  return 1;
}

void apply_partition(const marray* dest, const marray* m, const int row, const int col) {
  for(int i = 0; i < m->rows; i++) {
    for(int j = 0; j < m->cols; j++) {
      dest->data[(i + row) * dest->cols + j + col] = m->data[i * m->cols + j];
    }
  }
}

void from_partitions(const marray* dest, const marray* a11, const marray* a12, const marray* a21, const marray* a22) {
  const int n = dest->rows;

  apply_partition(dest, a11, 0, 0);
  apply_partition(dest, a12, 0, n/2);
  apply_partition(dest, a21, n/2, 0);
  apply_partition(dest, a22, n/2, n/2);
}

int dot_quadratic(marray* result, const marray* a, const marray* b) {
  marray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
  marray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

  // check if partitions created successfully
  if(!get_partitions(a, &a11, &a12, &a21, &a22) || !get_partitions(b, &b11, &b12, &b21, &b22)) {
    return 0;
  }

  marray *m1_a = matrix_add(a11, a22), *m1_b = matrix_add(b11, b22);
  marray *m2_a = matrix_add(a21, a22);
  marray *m3_b = matrix_sub(b12, b22);
  marray *m4_b = matrix_sub(b21, b11);
  marray *m5_a = matrix_add(a11, a12);
  marray *m6_a = matrix_sub(a21, a11), *m6_b = matrix_add(b11, b12);
  marray *m7_a = matrix_sub(a12, a22), *m7_b = matrix_add(b21, b22);

  marray *m1 = matrix_dot(m1_a, m1_b);
  marray *m2 = matrix_dot(m2_a, b11);
  marray *m3 = matrix_dot(a11, m3_b);
  marray *m4 = matrix_dot(a22, m4_b);
  marray *m5 = matrix_dot(m5_a, b22);
  marray *m6 = matrix_dot(m6_a, m6_b);
  marray *m7 = matrix_dot(m7_a, m7_b);

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

  if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
    return 0;
  }


  marray *c11 = matrix_add(m1, m4);
  matrix_subi(c11, m5);
  matrix_addi(c11, m7);

  marray *c12 = matrix_add(m3, m5);

  marray *c21 = matrix_add(m2, m4);

  marray *c22 = matrix_sub(m1, m2);
  matrix_addi(c22, m3);
  matrix_addi(c22, m6);

  matrix_free(m1);
  matrix_free(m2);
  matrix_free(m3);
  matrix_free(m4);
  matrix_free(m5);
  matrix_free(m6);
  matrix_free(m7);

  if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
    matrix_free(c11);
    matrix_free(c12);
    matrix_free(c21);
    matrix_free(c22);
    return 0;
  }

  from_partitions(result, c11, c12, c21, c22);

  matrix_free(c11);
  matrix_free(c12);
  matrix_free(c21);
  matrix_free(c22);

  return 1;
}

marray* matrix_dot(const marray* a, const marray* b) {

  // alloc matrix to store result
  marray* m = matrix_zeroes(a->rows, a->cols);

  if (m == NULL) {
    return NULL;
  }

  // check if a and b are quadratic
  if(a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
    if(a->cols % 2 == 0) {
      dot_quadratic(m, a, b);
      return m;
    }
  }

  dot_general(m, a, b);

  return m;
}
