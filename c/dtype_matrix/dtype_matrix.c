#include "dtype_matrix.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** PRIVATE FUNCTIONS **/

/*
Frees many matrices at once

Input:
    tarray *m[]: list of all matrices to free
    int n: amount of elements in the array

Output:
    None
*/
void tarray_free_many(tarray *m[], const int n) {
    for (int i = 0; i < n; i++) {
        tarray_free(m[i]);
    }
}


/** matrix_print sub functions based on tarray dtype **/

typedef void (*print_func_t)(const void *data, int index, int precision);

void print_matrix_float(const void *data, const int index, const int precision) {
    printf("%.*f\t", precision, ((const float*)data)[index]);
}

void print_matrix_double(const void *data, const int index, const int precision) {
    printf("%.*lf\t", precision, ((const double*)data)[index]);
}

void print_matrix_long_double(const void *data, const int index, const int precision) {
    printf("%.*Lf\t", precision, ((const long double*)data)[index]);
}

void print_matrix_complex_float(const void *data, const int index, const int precision) {
    const float complex z = ((const float complex*)data)[index];
    printf("%.*f %+.*fi\t\t", precision, crealf(z), precision, cimagf(z));
}

void print_matrix_complex_double(const void *data, const int index, const int precision) {
    const double complex z = ((const double complex*)data)[index];
    printf("%.*lf %+.*lfi\t\t", precision, creal(z), precision, cimag(z));
}

void print_matrix_complex_long_double(const void *data, const int index, const int precision) {
    const long double complex z = ((const long double complex*)data)[index];
    printf("%.*Lf %+.*Lfi\t\t", precision, creall(z), precision, cimagl(z));
}

/*
Returns the print function based on given data type

Input:
    MATRIX_DTYPE dtype: data type

Output:
    print_func_t: pointer to the matrix print function
*/
print_func_t get_print_func(const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            return print_matrix_float;
        case MATRIX_DTYPE_DOUBLE:
            return print_matrix_double;
        case MATRIX_DTYPE_LONG_DOUBLE:
            return print_matrix_long_double;
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return print_matrix_complex_float;
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return print_matrix_complex_double;
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return print_matrix_complex_long_double;
        default:
            return NULL;
    }
}


/** matrix assign functions based on dtype **/

typedef void (*assign_func_t)(const void *data, int index, const void *value);

void assign_matrix_float(const void *data, const int index, const void *value) {
    ((float*)data)[index] = *(float*)value;
}

void assign_matrix_double(const void *data, const int index, const void *value) {
    ((double*)data)[index] = *(double*)value;
}

void assign_matrix_long_double(const void *data, const int index, const void *value) {
    ((long double*)data)[index] = *(long double*)value;
}

void assign_matrix_complex_float(const void *data, const int index, const void *value) {
    ((float complex*)data)[index] = *(float complex*)value;
}

void assign_matrix_complex_double(const void *data, const int index, const void *value) {
    ((double complex*)data)[index] = *(double complex*)value;
}

void assign_matrix_complex_long_double(const void *data, const int index, const void *value) {
    ((long double complex*)data)[index] = *(long double complex*)value;
}

/*
Returns the assign function based on given data type

Input:
    MATRIX_DTYPE dtype: data type

Output:
    assign_func_t: pointer to the matrix assign function
*/
assign_func_t get_assign_func(const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            return assign_matrix_float;
        case MATRIX_DTYPE_DOUBLE:
            return assign_matrix_double;
        case MATRIX_DTYPE_LONG_DOUBLE:
            return assign_matrix_long_double;
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return assign_matrix_complex_float;
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return assign_matrix_complex_double;
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return assign_matrix_complex_long_double;
        default:
            return NULL;
    }
}

/*
Returns the size of a single data element

Input:
    MATRIX_DTYPE: data type of the matrix

Output:
    size_t: memory space in bytes of a single data element
*/
size_t dtype_data_size(const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            return sizeof(float);
        case MATRIX_DTYPE_DOUBLE:
            return sizeof(double);
        case MATRIX_DTYPE_LONG_DOUBLE:
            return sizeof(long double);
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return sizeof(float complex);
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return sizeof(double complex);
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return sizeof(long double complex);
        default:
            return 0;
    }
}


/*
Calculates the required memory space to store matrix data

Input:
    int elems: amount of elements inside matrix (rows * cols)
    MATRIX_DTYPE: data type of matrix elements

Output:
    size_t: required memory space in bytes to store matrix data
*/
size_t matrix_data_size(const int elems, const MATRIX_DTYPE dtype) {
    return elems * dtype_data_size(dtype);
}


/** add functions **/

#define MATRIX_OP unsigned char
#define MATRIX_OP_ADD 0
#define MATRIX_OP_SUB 1
#define MATRIX_OP_MUL 2
#define MATRIX_OP_DIV 3

#define MATRIX_OPS_AMOUNT 4

typedef void (*op_val_double_func_t)(MATRIX_OP op, const tarray *target, const tarray *data, double value);

typedef void (*op_func_t)(const void *target, const void *a, const void *b);

void float_op_add(const void *target, const void *a, const void *b) {
    *(float*)target = *(float*)a + *(float*)b;
}

void float_op_sub(const void *target, const void *a, const void *b) {
    *(float*)target = *(float*)a - *(float*)b;
}

void float_op_mul(const void *target, const void *a, const void *b) {
    *(float*)target = *(float*)a * *(float*)b;
}

void float_op_div(const void *target, const void *a, const void *b) {
    *(float*)target = *(float*)a / *(float*)b;
}

void double_op_add(const void *target, const void *a, const void *b) {
    *(double*)target = *(double*)a + *(double*)b;
}

void double_op_sub(const void *target, const void *a, const void *b) {
    *(double*)target = *(double*)a - *(double*)b;
}

void double_op_mul(const void *target, const void *a, const void *b) {
    *(double*)target = *(double*)a * *(double*)b;
}

void double_op_div(const void *target, const void *a, const void *b) {
    *(double*)target = *(double*)a / *(double*)b;
}

void long_double_op_add(const void *target, const void *a, const void *b) {
    *(long double*)target = *(long double*)a + *(long double*)b;
}

void long_double_op_sub(const void *target, const void *a, const void *b) {
    *(long double*)target = *(long double*)a - *(long double*)b;
}

void long_double_op_mul(const void *target, const void *a, const void *b) {
    *(long double*)target = *(long double*)a * *(long double*)b;
}

void long_double_op_div(const void *target, const void *a, const void *b) {
    *(long double*)target = *(long double*)a / *(long double*)b;
}

void complex_float_op_add(const void *target, const void *a, const void *b) {
    *(float complex*)target = *(float complex*)a + *(float complex*)b;
}

void complex_float_op_sub(const void *target, const void *a, const void *b) {
    *(float complex*)target = *(float complex*)a - *(float complex*)b;
}

void complex_float_op_mul(const void *target, const void *a, const void *b) {
    *(float complex*)target = *(float complex*)a * *(float complex*)b;
}

void complex_float_op_div(const void *target, const void *a, const void *b) {
    *(float complex*)target = *(float complex*)a / *(float complex*)b;
}

void complex_double_op_add(const void *target, const void *a, const void *b) {
    *(double complex*)target = *(double complex*)a + *(double complex*)b;
}

void complex_double_op_sub(const void *target, const void *a, const void *b) {
    *(double complex*)target = *(double complex*)a - *(double complex*)b;
}

void complex_double_op_mul(const void *target, const void *a, const void *b) {
    *(double complex*)target = *(double complex*)a * *(double complex*)b;
}

void complex_double_op_div(const void *target, const void *a, const void *b) {
    *(double complex*)target = *(double complex*)a / *(double complex*)b;
}

void complex_long_double_op_add(const void *target, const void *a, const void *b) {
    *(long double complex*)target = *(long double complex*)a + *(long double complex*)b;
}

void complex_long_double_op_sub(const void *target, const void *a, const void *b) {
    *(long double complex*)target = *(long double complex*)a - *(long double complex*)b;
}

void complex_long_double_op_mul(const void *target, const void *a, const void *b) {
    *(long double complex*)target = *(long double complex*)a * *(long double complex*)b;
}

void complex_long_double_op_div(const void *target, const void *a, const void *b) {
    *(long double complex*)target = *(long double complex*)a / *(long double complex*)b;
}

/*
Returns the operation function for the given operation based on the data type provided

Input:
    MATRIX_OP op: operation
    MATRIX_DTYPE dtype: data type

Output:
    op_func_t: function pointer of requested operation
*/
op_func_t get_op(const MATRIX_OP op, const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            const op_func_t op_funcs_float[] = {float_op_add, float_op_sub, float_op_mul, float_op_div};
            return  op_funcs_float[op];

        case MATRIX_DTYPE_DOUBLE:
            const op_func_t op_funcs_double[] = {double_op_add, double_op_sub, double_op_mul, double_op_div};
            return op_funcs_double[op];

        case MATRIX_DTYPE_LONG_DOUBLE:
            const op_func_t op_funcs_long_double[] = {long_double_op_add, long_double_op_sub, long_double_op_mul, long_double_op_div};
            return  op_funcs_long_double[op];

        case MATRIX_DTYPE_COMPLEX_FLOAT:
            const op_func_t op_funcs_complex_float[] = {complex_float_op_add, complex_float_op_sub, complex_float_op_mul, complex_float_op_div};
            return op_funcs_complex_float[op];

        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            const op_func_t op_funcs_complex_double[] = {complex_double_op_add, complex_double_op_sub, complex_double_op_mul, complex_double_op_div};
            return  op_funcs_complex_double[op];

        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            const op_func_t op_funcs_complex_long_double[] = {complex_long_double_op_add, complex_long_double_op_sub, complex_long_double_op_mul, complex_long_double_op_div};
            return op_funcs_complex_long_double[op];
        default:
            return NULL;
    }
}

/*
Performs the given operation on the matrix with the value

Input:
    MATRIX_OP op: operation code
    tarray *target: target matrix
    tarray *m: source matrix
    double val: value to add/sub/mul/div with

Output:
    None, but the target matrix will be modified
*/
void op(const MATRIX_OP op, const tarray *target, const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    if (op > MATRIX_OPS_AMOUNT) {
        return;
    }

    op_func_t op_func;

    switch (target->dtype) {
        case MATRIX_DTYPE_FLOAT:
            const op_func_t op_funcs_float[] = {float_op_add, float_op_sub, float_op_mul, float_op_div};
            op_func = op_funcs_float[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((float*)target->data)[i], &((float*)m->data)[i], (float*)&val);
            }

            return;

        case MATRIX_DTYPE_DOUBLE:
            const op_func_t op_funcs_double[] = {double_op_add, double_op_sub, double_op_mul, double_op_div};
            op_func = op_funcs_double[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((double*)target->data)[i], &((double*)m->data)[i], (double*)&val);
            }

            return;

        case MATRIX_DTYPE_LONG_DOUBLE:
            const op_func_t op_funcs_long_double[] = {long_double_op_add, long_double_op_sub, long_double_op_mul, long_double_op_div};
            op_func = op_funcs_long_double[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((long double*)target->data)[i], &((long double*)m->data)[i], (long double*)&val);
            }

            return;

        case MATRIX_DTYPE_COMPLEX_FLOAT:
            const op_func_t op_funcs_complex_float[] = {complex_float_op_add, complex_float_op_sub, complex_float_op_mul, complex_float_op_div};
            op_func = op_funcs_complex_float[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((float complex*)target->data)[i], &((float complex*)m->data)[i], (float complex*)&val);
            }

            return;

        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            const op_func_t op_funcs_complex_double[] = {complex_double_op_add, complex_double_op_sub, complex_double_op_mul, complex_double_op_div};
            op_func = op_funcs_complex_double[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((double complex*)target->data)[i], &((double complex*)m->data)[i], (double complex*)&val);
            }

            return;

        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            const op_func_t op_funcs_complex_long_double[] = {complex_long_double_op_add, complex_long_double_op_sub, complex_long_double_op_mul, complex_long_double_op_div};
            op_func = op_funcs_complex_long_double[op];

            for (int i = 0; i < elems; i++) {
                op_func(&((long double complex*)target->data)[i], &((long double complex*)m->data)[i], (long double complex*)&val);
            }

            return;
        default:
            return;
    }
}


/** matrix add matrix function **/

typedef void (*op_matrix_func_t)(const tarray *target, const tarray *a, const tarray *b);

void add_matrix_float(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((float*)target->data)[i] = ((float*)a->data)[i] + ((float*)b->data)[i];
    }
}

void add_matrix_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((double*)target->data)[i] = ((double*)a->data)[i] + ((double*)b->data)[i];
    }
}

void add_matrix_long_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((long double*)target->data)[i] = ((long double*)a->data)[i] + ((long double*)b->data)[i];
    }
}

void add_matrix_complex_float(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((float complex*)target->data)[i] = ((float complex*)a->data)[i] + ((float complex*)b->data)[i];
    }
}

void add_matrix_complex_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((double complex*)target->data)[i] = ((double complex*)a->data)[i] + ((double complex*)b->data)[i];
    }
}

void add_matrix_complex_long_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((long double complex*)target->data)[i] = ((long double complex*)a->data)[i] + ((long double complex*)b->data)[i];
    }
}

void sub_matrix_float(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((float*)target->data)[i] = ((float*)a->data)[i] - ((float*)b->data)[i];
    }
}

void sub_matrix_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((double*)target->data)[i] = ((double*)a->data)[i] - ((double*)b->data)[i];
    }
}

void sub_matrix_long_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((long double*)target->data)[i] = ((long double*)a->data)[i] - ((long double*)b->data)[i];
    }
}

void sub_matrix_complex_float(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((float complex*)target->data)[i] = ((float complex*)a->data)[i] - ((float complex*)b->data)[i];
    }
}

void sub_matrix_complex_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((double complex*)target->data)[i] = ((double complex*)a->data)[i] - ((double complex*)b->data)[i];
    }
}

void sub_matrix_complex_long_double(const tarray *target, const tarray *a, const tarray *b) {
    const int elems = target->rows * target->cols;

    for (int i = 0; i < elems; i++) {
        ((long double complex*)target->data)[i] = ((long double complex*)a->data)[i] - ((long double complex*)b->data)[i];
    }
}

/*
Returns the matrix add function based on given data type

Input:
    MATRIX_DTYPE dtype: data type

Output:
    assign_func_t: pointer to the matrix add function
*/
op_matrix_func_t get_matrix_add_func(const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            return add_matrix_float;
        case MATRIX_DTYPE_DOUBLE:
            return add_matrix_double;
        case MATRIX_DTYPE_LONG_DOUBLE:
            return add_matrix_long_double;
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return add_matrix_complex_float;
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return add_matrix_complex_double;
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return add_matrix_complex_long_double;
        default:
            return NULL;
    }
}

/*
Returns the matrix add function based on given data type

Input:
    MATRIX_DTYPE dtype: data type

Output:
    assign_func_t: pointer to the matrix add function
*/
op_matrix_func_t get_matrix_sub_func(const MATRIX_DTYPE dtype) {
    switch (dtype) {
        case MATRIX_DTYPE_FLOAT:
            return sub_matrix_float;
        case MATRIX_DTYPE_DOUBLE:
            return sub_matrix_double;
        case MATRIX_DTYPE_LONG_DOUBLE:
            return sub_matrix_long_double;
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return sub_matrix_complex_float;
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return sub_matrix_complex_double;
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return sub_matrix_complex_long_double;
        default:
            return NULL;
    }
}


typedef void (*mul_matrix_func_t)(const tarray *a, double b);

void float_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((float*)m->data)[i] *= (float)val;
    }
}

void double_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((double*)m->data)[i] *= val;
    }
}

void long_double_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((long double*)m->data)[i] *= val;
    }
}

void complex_float_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((float complex*)m->data)[i] *= (float)val;
    }
}

void complex_double_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((double complex*)m->data)[i] *= val;
    }
}

void complex_long_double_mul(const tarray *m, const double val) {
    const int elems = m->rows * m->cols;

    for (int i = 0; i < elems; i++) {
        ((long double complex*)m->data)[i] *= val;
    }
}


/** matrix dot **/

/*
Performs the dot multiplication with the general method

Input:
    tarray *result: pointer to a matrix to store the result
    tarray *a: a pointer to a matrix
    tarray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
void tarray_dot_general(const tarray *result, const tarray *a, const tarray *b) {
    const op_func_t op_add = get_op(MATRIX_OP_ADD, result->dtype);
    const op_func_t op_mul = get_op(MATRIX_OP_MUL, result->dtype);
    const assign_func_t assign_func = get_assign_func(result->dtype);

    const size_t data_size = dtype_data_size(result->dtype);
    void *val = calloc(data_size, 3);
    void *temp_val = val + data_size;
    void *zero_val = val + 2*data_size;

    // check for any errors
    if (op_add == NULL || op_mul == NULL || assign_func == NULL || val == NULL) {
        return;
    }

    for (int i = 0; i < a->rows; i++) {
        const int row = i * a->cols;

        for (int j = 0; j < b->cols; j++) {
            assign_func(val, 0, &zero_val);

            for (int k = 0; k < a->cols; k++) {
                //assign_func(temp_val, 0, &zero_val);
                op_mul(temp_val, a->data + (row + k) * data_size, b->data + (k * b->cols + j) * data_size);
                op_add(val, val, temp_val);
                //val += a->data[row + k] * b->data[k * b->cols + j];
            }

            assign_func(result->data, i * result->cols + j, val);
            //result->data[i * result->cols + j] = val;
        }
    }

    free(val);
}


/*
Splits the matrix into 4 partitions

Input:
    tarray *a: pointer to the matrix from which the partitions are created
    tarray **a11: upper left partition
    tarray **a12: upper right partition
    tarray **a21: bottom left partition
    tarray **a22: bottom right partition

Output:
    None, but the 4 partitions are stored in input pointers
*/
int tarray_get_partitions(const tarray *a, tarray **a11, tarray **a12, tarray **a21, tarray **a22) {
    tarray_get_partition(a11, a, 0, a->rows / 2, 0, a->cols / 2);
    tarray_get_partition(a12, a, 0, a->rows / 2, a->cols / 2, a->cols);
    tarray_get_partition(a21, a, a->rows / 2, a->rows, 0, a->cols / 2);
    tarray_get_partition(a22, a, a->rows / 2, a->rows, a->cols / 2, a->cols);

    if (*a11 == NULL || *a12 == NULL || *a21 == NULL || *a22 == NULL) {
        // free possibly created matrices
        tarray_free(*a11);
        tarray_free(*a12);
        tarray_free(*a21);
        tarray_free(*a22);
        return 0;
    }

    return 1;
}


/*
Inserts the given partition into the matrix

Input:
    tarray *dest: destination matrix
    tarray *m: partition matrix that is applied
    int row: start row
    int col: start col

Output:
    None, but the partition is applied to the target matrix
*/
void tarray_apply_partition(const tarray *dest, const tarray *m, const int row, const int col) {
    const size_t data_size = dtype_data_size(dest->dtype);

    if (data_size == 0) {
        return;
    }

    for (int i = 0; i < m->rows; i++) {
        const int dest_offset = (i + row) * dest->cols + col;
        const int src_row = i * m->cols;

        memcpy(dest->data + dest_offset * data_size,
            m->data + src_row * data_size,
            m->cols * data_size);
    }
}


/*
Inserts all 4 partitions into the destination matrix

Input:
    tarray *dest: destination matrix
    tarray **a11: upper left partition
    tarray **a12: upper right partition
    tarray **a21: bottom left partition
    tarray **a22: bottom right partition

Output:
    None, but the destination matrix is the combination of all 4 partitions
*/
void tarray_from_partitions(const tarray *dest, const tarray *a11, const tarray *a12, const tarray *a21, const tarray *a22) {
    const int n = dest->rows;

    tarray_apply_partition(dest, a11, 0, 0);
    tarray_apply_partition(dest, a12, 0, n / 2);
    tarray_apply_partition(dest, a21, n / 2, 0);
    tarray_apply_partition(dest, a22, n / 2, n / 2);
}


/*
Calculates the dot product of two quadratic matrices with even dimension

Input:
    tarray *result: pointer to a matrix to store the result
    tarray *a: a pointer to a matrix
    tarray *b: a pointer to a matrix

Output:
    None, but the result is stored in *result
*/
bool tarray_dot_quadratic(const tarray *result, const tarray *a, const tarray *b) {
    //TODO Think about better solution to free matrices if errors occur

    tarray *a11 = NULL, *a12 = NULL, *a21 = NULL, *a22 = NULL;
    tarray *b11 = NULL, *b12 = NULL, *b21 = NULL, *b22 = NULL;

    // check if partitions created successfully
    if (!tarray_get_partitions(a, &a11, &a12, &a21, &a22) || !tarray_get_partitions(b, &b11, &b12, &b21, &b22)) {
        return 0;
    }

    // intermediate steps
    tarray *m1_a = tarray_add(a11, a22), *m1_b = tarray_add(b11, b22);
    tarray *m2_a = tarray_add(a21, a22);
    tarray *m3_b = tarray_sub(b12, b22);
    tarray *m4_b = tarray_sub(b21, b11);
    tarray *m5_a = tarray_add(a11, a12);
    tarray *m6_a = tarray_sub(a21, a11), *m6_b = tarray_add(b11, b12);
    tarray *m7_a = tarray_sub(a12, a22), *m7_b = tarray_add(b21, b22);

    // check for successful matrix operations
    if (m1_a == NULL || m1_b == NULL || m2_a == NULL || m3_b == NULL || m4_b == NULL || m5_a == NULL || m6_a == NULL ||
        m6_b == NULL || m7_a == NULL || m7_b == NULL) {
        tarray *arrays_1[] = {
            a11, a12, a21, a22, // 4
            b11, b12, b21, b22, // 4
            m1_a, m1_b, m2_a, m3_b, m4_b, m5_a, m6_a, m6_b, m7_a, m7_b}; // 10

        tarray_free_many(arrays_1, 18);

        return 0;
    }

    tarray *m1 = tarray_dot(m1_a, m1_b);
    tarray *m2 = tarray_dot(m2_a, b11);
    tarray *m3 = tarray_dot(a11, m3_b);
    tarray *m4 = tarray_dot(a22, m4_b);
    tarray *m5 = tarray_dot(m5_a, b22);
    tarray *m6 = tarray_dot(m6_a, m6_b);
    tarray *m7 = tarray_dot(m7_a, m7_b);

    // free up memory that is no longer required
    tarray *arrays_1[] = {
        a11, a12, a21, a22, // 4
        b11, b12, b21, b22, // 4
        m1_a, m1_b, m2_a, m3_b, m4_b, m5_a, m6_a, m6_b, m7_a, m7_b}; // 10

    tarray_free_many(arrays_1, 18);

    // check for successful matrix operations
    if (m1 == NULL || m2 == NULL || m3 == NULL || m4 == NULL || m5 == NULL || m6 == NULL || m7 == NULL) {
        tarray *arrays_2[] = {m1, m2, m3, m4, m5, m6, m7};
        tarray_free_many(arrays_2, 7);

        return 0;
    }

    // calculate final partitions
    tarray *c11 = tarray_add(m1, m4);
    tarray_subi(c11, m5);
    tarray_addi(c11, m7);

    tarray *c12 = tarray_add(m3, m5);

    tarray *c21 = tarray_add(m2, m4);

    tarray *c22 = tarray_sub(m1, m2);
    tarray_addi(c22, m3);
    tarray_addi(c22, m6);

    // free up memory that is no longer required
    tarray *arrays_2[] = {m1, m2, m3, m4, m5, m6, m7};
    tarray_free_many(arrays_2, 7);

    // check for successful matrix operations
    if (c11 == NULL || c12 == NULL || c21 == NULL || c22 == NULL) {
        tarray *arrays_3[] = {c11, c12, c21, c22};
        tarray_free_many(arrays_3, 4);

        return 0;
    }

    tarray_from_partitions(result, c11, c12, c21, c22);

    tarray *arrays_3[] = {c11, c12, c21, c22};
    tarray_free_many(arrays_3, 4);

    return 1;
}


/*
Checks whether the individual elements of both matrices are close to each other

Input:
    tarray *a: first matrix
    tarray *b: second matrix, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    bool: result of whether the matrices are close to each other
*/
bool entry_close(const tarray *a, const tarray *b, const int index, const double rtol, const double atol) {
    switch (a->dtype) {
        case MATRIX_DTYPE_FLOAT:
            return fabsf(((float*)a->data)[index] - ((float*)b->data)[index]) <= atol + rtol * fabsf(((float*)b->data)[index]);
        case MATRIX_DTYPE_DOUBLE:
            return fabs(((double*)a->data)[index] - ((double*)b->data)[index]) <= atol + rtol * fabs(((double*)b->data)[index]);
        case MATRIX_DTYPE_LONG_DOUBLE:
            return fabsl(((long double*)a->data)[index] - ((long double*)b->data)[index]) <= atol + rtol * fabsl(((long double*)b->data)[index]);
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            return cabsf(((float complex*)a->data)[index] - ((float complex*)b->data)[index]) <= atol + rtol * cabsf(((float complex*)b->data)[index]);
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            return cabs(((double complex*)a->data)[index] - ((double complex*)b->data)[index]) <= atol + rtol * cabs(((double complex*)b->data)[index]);
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            return cabsl(((long double complex*)a->data)[index] - ((long double complex*)b->data)[index]) <= atol + rtol * cabsl(((long double complex*)b->data)[index]);
        default:
            return 0;
    }
}


/* ========================================================================================= */


/** PUBLIC HEADER FUNCTIONS **/

void tarray_free(tarray *m) {
    if (m == NULL) {
        return;
    }

    free(m);
}


/** matrix print functions **/

void tarray_print_precision(const tarray *m, const int precision) {
    if (m == NULL) {
        return;
    }

    const print_func_t print_func = get_print_func(m->dtype);

    for (int i = 0; i < m->rows; i++) {
        const int row = i * m->cols;

        for (int j = 0; j < m->cols; j++) {
            print_func(m->data, row + j, precision);
        }
        printf("\n");
    }
}


void tarray_print(const tarray *m) {
    tarray_print_precision(m, MATRIX_PRINT_PRECISION_DEFAULT);
}


/** matrix initialize functions **/

tarray *tarray_zeroes(const int rows, const int cols, const MATRIX_DTYPE dtype) {
    // calc space based on selected dtype
    const size_t data_size = matrix_data_size(rows * cols, dtype);

    // return NULL if dtype is unknown
    if (data_size == 0) {
        return NULL;
    }

    // alloc mem for structure and data
    // calloc sets the entire memory area to 0
    // ReSharper disable CppDFAMemoryLeak
    tarray *m = calloc(sizeof(tarray) + data_size, 1);

    // return null if not enough memory is free
    if (m == NULL) {
        return NULL;
    }

    // apply dimensions and dtype
    m->rows = rows;
    m->cols = cols;
    m->dtype = dtype;

    // set data pointer to memory after structure
    m->data = m + 1;

    return m;
}


tarray *tarray_identity(const int n, const MATRIX_DTYPE dtype) {
    tarray *m = tarray_zeroes(n, n, dtype);

    if (m == NULL) {
        return NULL;
    }

    const size_t data_size = dtype_data_size(dtype);
    if (data_size == 0) {
        free(m);
        return NULL;
    }

    void *value = malloc(data_size);

    // pointer to assign function
    assign_func_t assign_func = NULL;

    // select assign function based on matrix dtype
    switch (m->dtype) {
        case MATRIX_DTYPE_FLOAT:
            assign_func = assign_matrix_float;
            *(float*)value = 1.0f;
        break;
        case MATRIX_DTYPE_DOUBLE:
            assign_func = assign_matrix_double;
            *(double*)value = 1.0;
        break;
        case MATRIX_DTYPE_LONG_DOUBLE:
            assign_func = assign_matrix_long_double;
            *(long double*)value = 1.0;
        break;
        case MATRIX_DTYPE_COMPLEX_FLOAT:
            assign_func = assign_matrix_complex_float;
            *(float complex*)value = 1.0f;
        break;
        case MATRIX_DTYPE_COMPLEX_DOUBLE:
            assign_func = assign_matrix_complex_double;
            *(double complex*)value = 1.0;
        break;
        case MATRIX_DTYPE_COMPLEX_LONG_DOUBLE:
            assign_func = assign_matrix_complex_long_double;
            *(long double complex*)value = 1.0;
        break;
        default:
            return NULL;
    }

    for (int i = 0; i < n; i++) {
        assign_func(m->data, i * n + i, value);
    }

    free(value);

    return m;
}


tarray *tarray_copy(const tarray *m) {
    tarray *copy = tarray_zeroes(m->rows, m->cols, m->dtype);

    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, m->data, matrix_data_size(m->rows * m->cols, m->dtype));

    return copy;
}


tarray *tarray_transposed(const tarray *m) {
    tarray *t = tarray_zeroes(m->cols, m->rows, m->dtype);

    if (t == NULL) {
        return NULL;
    }

    const assign_func_t assign_func = get_assign_func(m->dtype);
    const size_t data_size = dtype_data_size(m->dtype);
    if (assign_func == NULL) {
        tarray_free(t);
        return NULL;
    }

    for (int i = 0; i < t->rows; i++) {
        const int row = i * t->cols;

        for (int j = 0; j < t->cols; j++) {
            assign_func(t->data, j * m->cols + i, m->data + (row + j) * data_size);
        }
    }

    return t;
}


/** matrix calculation functions **/

void tarray_addi_val(const tarray *a, const double b) {
    if (a == NULL) {
        return;
    }

    op(MATRIX_OP_ADD, a, a, b);
}


tarray *tarray_add_val(const tarray *a, const double b) {
    tarray *result = tarray_zeroes(a->rows, a->cols, a->dtype);

    if (result == NULL) {
        return NULL;
    }

    op(MATRIX_OP_ADD, result, a, b);

    return result;
}


void tarray_addi(const tarray *a, const tarray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const op_matrix_func_t add_func = get_matrix_add_func(a->dtype);
    if (add_func == NULL) {
        return;
    }

    add_func(a, a, b);
}


tarray *tarray_add(const tarray *a, const tarray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    tarray *m = tarray_zeroes(a->rows, a->cols, a->dtype);

    if (m == NULL) {
        return NULL;
    }

    const op_matrix_func_t add_func = get_matrix_add_func(a->dtype);
    if (add_func == NULL) {
        return NULL;
    }

    add_func(m, a, b);

    return m;
}


void tarray_subi(const tarray *a, const tarray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return;
    }

    const op_matrix_func_t sub_func = get_matrix_sub_func(a->dtype);
    if (sub_func == NULL) {
        return;
    }

    sub_func(a, a, b);
}


tarray *tarray_sub(const tarray *a, const tarray *b) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    tarray *m = tarray_zeroes(a->rows, a->cols, a->dtype);

    if (m == NULL) {
        return NULL;
    }

    const op_matrix_func_t sub_func = get_matrix_sub_func(a->dtype);
    if (sub_func == NULL) {
        return NULL;
    }

    sub_func(m, a, b);

    return m;
}


void tarray_muli_val(const tarray *a, const double b) {
    op(MATRIX_OP_MUL, a, a, b);
}


tarray *tarray_mul_val(const tarray *a, const double b) {
    tarray *result = tarray_zeroes(a->rows, a->cols, a->dtype);

    if (result == NULL) {
        return NULL;
    }

    op(MATRIX_OP_MUL, result, a, b);

    return result;
}


void tarray_get_partition(tarray **dest, const tarray *m, const int row_start, const int row_end, const int col_start,
                          const int col_end) {
    *dest = tarray_zeroes(row_end - row_start, col_end - col_start, m->dtype);

    if (*dest == NULL) {
        return;
    }

    const size_t data_size = dtype_data_size(m->dtype);

    for (int i = 0; i < row_end - row_start; i++) {
        // dest: i-th row of p
        // src: i + row_start of matrix a, skip rows behind start row
        memcpy((*dest)->data + i * (*dest)->cols * data_size,
            m->data + ((i + row_start) * m->cols + col_start) * data_size,
               data_size * (*dest)->cols);
    }
}


tarray *tarray_dot(const tarray *a, const tarray *b) {
    // TODO Think about (required for linalg back substitution)
    if (a->rows == 0 || a->cols == 0 || b->rows == 0 || b->cols == 0) {
        tarray *m = tarray_zeroes(1, 1, a->dtype);

        return m;
    }

    // alloc matrix to store result
    tarray *m = tarray_zeroes(a->rows, b->cols, a->dtype);

    if (m == NULL) {
        return NULL;
    }

    // check if a and b are quadratic
    if (a->rows == b->rows && a->cols == b->cols && a->rows == a->cols) {
        if (a->cols % 2 == 0) {
            if (tarray_dot_quadratic(m, a, b)) {
                return m;
            }

            // dot_quadratic was not successful -> free matrix and return NULL
            tarray_free(m);

            return NULL;
        }
    }

    tarray_dot_general(m, a, b);

    return m;
}


tarray *tarray_close(const tarray *a, const tarray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return NULL;
    }

    // alloc matrix to store result
    tarray *results = tarray_zeroes(a->rows, a->cols, MATRIX_DTYPE_FLOAT);

    if (results == NULL) {
        return NULL;
    }

    const int elems = a->rows * a->cols;

    // get result for each entry and store separately
    for (int i = 0; i < elems; i++) {
        ((float*)results->data)[i] = (float)entry_close(a, b, i, rtol, atol);
    }

    return results;
}


bool tarray_close_all(const tarray *a, const tarray *b, const double rtol, const double atol) {
    // Check dimensions
    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }

    const int elems = a->rows * a->cols;

    // check each entry, return false if one mismatches
    for (int i = 0; i < elems; i++) {
        if (! entry_close(a, b, i, rtol, atol)) {
            return false;
        }
    }

    return true;
}
