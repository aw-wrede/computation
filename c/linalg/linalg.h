#ifndef LINALG_H
#define LINALG_H

#include <stdbool.h>
#include "../matrix/matrix.h"

void linalg_print_system(const marray* a, const marray* b);

bool linalg_gaussian_elimination(const marray* a, const marray* b);

marray* linalg_back_substitution(const marray* a, const marray* b);

marray* linalg_solve(const marray* a, const marray* b);

#endif //LINALG_H
