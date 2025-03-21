#ifndef NUMBERS_H
#define NUMBERS_H
#include <stdbool.h>

#endif //NUMBERS_H

/*
Checks whether two given numbers are close to each other within the given tolerances

Input:
    double a: first number
    double b: second number, used as reference
    double rtol: relative tolerance, which depends on the value b (default 1e-05)
    double atol: absolute tolerance (default 1e-08)

Output:
    bool: Whether the numbers are close to each other
*/
bool numbers_close(double a, double b, double rtol, double atol);