#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
  int rows, cols;
  double *data;
} marray;

void matrix_free(marray* m);

void matrix_print(const marray* m);

marray* matrix_zeroes(int rows, int cols);

void matrix_addi_val(const marray* a, double b);

marray* matrix_add_val(const marray* a, double b);

void matrix_addi(const marray* a, const marray* b);

marray* matrix_add(const marray* a, const marray* b);

void matrix_subi(const marray* a, const marray* b);

marray* matrix_sub(const marray* a, const marray* b);

marray* matrix_dot(const marray *a, const marray *b);

#endif //MATRIX_H
