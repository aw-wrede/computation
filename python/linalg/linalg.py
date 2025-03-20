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

def back_substitution(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Performs back substitution on given upper triangular matrix a with b
    :param a: Upper triangular matrix containing the linear equations
    :param b: Matrix containing the results of the linear equations
    :return: Solution vector
    """
    n, m = a.shape
    if n != m:
        raise ValueError("Matrix a must be square!")

    if b.shape[0] != n:
        raise ValueError("Vector b must have same number of rows as matrix A!")

    # Initialize solution vector with proper size
    x = np.zeros(m)

    for i in range(m - 1, -1, -1):
        pivot = a[i, i]

        if pivot == 0:
            raise ValueError("Matrix a has zero at pivot element [{}, {}] and is not solvable!".format(i, i))

        x[i] = ( b[i] - np.dot(a[i, i+1:], x[i+1:]) ) / pivot
        print("i={}".format(i))
        print(a[i, i+1:])
        print(x[i+1:])

    return x

def solve(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solves the given system of linear equations
    :param a: Matrix containing the linear equations
    :param b: Matrix containing the results of the linear equations
    :return: Solution vector
    """
    # Create copies of input matrix and vector to leave them unmodified
    a = a.copy()
    b = b.copy()

    gaussian_elimination(a, b)
    return back_substitution(a, b)
