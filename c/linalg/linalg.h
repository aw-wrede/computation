#ifndef LINALG_H
#define LINALG_H

#include <stdbool.h>
#include "../matrix/matrix.h"

/*
Prints the given linear equation system to console

Input:
  marray* a: matrix containing the linear equations
  marray* b: matrix containing the results of linear equations

Output:
  None, but the linear equation system is displayed in the console
*/
void linalg_print_system(const marray* a, const marray* b);

/*
Performs Gaussian elimination on the input matrices a and b

Input:
  marray* a: matrix containing the linear equations
  marray* b: matrix containing the results of linear equations

Output:
  bool: success of Gaussian elimination
  modified matrices a and b
*/
bool linalg_gaussian_elimination(const marray* a, const marray* b);

/*
Performs back substitution on given upper triangular matrix a with b

Input:
  marray* a: upper triangular matrix containing the linear equations
  marray* b: matrix containing the results of linear equations

Output:
  marray*: solution vector
*/
marray* linalg_back_substitution(const marray* a, const marray* b);

/*
Solves the given system of linear equations

Input:
  marray* a: matrix containing the linear equations
  marray* b: matrix containing the results of linear equations

Output:
  marray*: solution vector of given linear equation system
*/
marray* linalg_solve(const marray* a, const marray* b);

#endif //LINALG_H
