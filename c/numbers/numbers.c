#include "numbers.h"

#include <math.h>

double numbers_compute_epsilon() {
    double eps = 1.0;
    const double one = 1.0;

    while (one + eps != one) {
        eps = eps / 2.0;
    }

    eps *= 2;

    return eps;
}

bool numbers_is_close(const double a, const double b, const double rtol, const double atol) {
    return fabs(a - b) <= atol + rtol * fabs(b);
}