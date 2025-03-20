from encodings import undefined

import numpy as np


def gaussian_elimination(a: np.ndarray, b:np.ndarray):
    """
    Performs Gaussian elimination on the input matrices a and b
    :param a: Matrix containing the linear equations
    :param b: Matrix containing the results of the linear equations
    :return: matrices a and b after Gaussian elimination
    """

    def find_largest_pivot_element(a: np.ndarray, col:int) -> int or None:
        pivot_index = None
        pivot = None

        for i in range(col, a.shape[0]):
            if pivot is None or (a[i, col] > pivot and a[i, col] != 0):
                pivot = a[i, col]
                pivot_index = i

        return pivot_index

    n, m = a.shape
    if n != m:
        raise ValueError("Matrix a must be square!")

    if b.shape[0] != n:
        raise ValueError("Vector b must have same number of rows as matrix A!")

    # Perform gaussian elimination
    for i in range(n):
        # Find row with the largest pivot element
        # (prevent unsolvable if largest element would be zero because other candidates are negative)
        max_row = find_largest_pivot_element(a, i)

        if a[max_row, i] is None:
            raise ValueError("Matrix a is singular and not solvable!")

        # Swap rows in A and b
        if max_row != i:
            a[[i, max_row]] = a[[max_row, i]]
            b[[i, max_row]] = b[[max_row, i]]

        pivot = a[i, i]

        # Eliminate all fields below pivot element
        for row in range(i + 1, n):
            factor = a[row, i] / pivot

            if a[row, i] == 0:
                continue

            for column in range(i, n):
                a[row, column] -= factor * a[i, column]

            b[row] -= factor * b[i]

    return a, b