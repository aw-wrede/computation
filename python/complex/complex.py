import numpy as np


def dft_matrix(n: int, dtype:np.dtype = np.complex128) -> np.ndarray:
    """
    Constructs the DFT matrix of size n*n
    :param n: size of DFT matrix
    :param dtype: data type of the matrix
    :return: DFT matrix of size n*n
    """
    f = np.ones((n, n), dtype=dtype)

    # Calculate omega (counter clockwise -> positive sign)
    omega = np.exp(-1j * 2 * np.pi / n)

    # Calculate all different unit roots
    omegas = np.zeros(n, dtype=dtype)
    n_sqrt = np.sqrt(n)
    for i in range(n):
        # Normalize unit root with 1/sqrt(n)
        omegas[i] = (omega ** i) / n_sqrt

    # Fill matrix with unit roots
    for k in range(n):
        for l in range(k, n):
            result = omegas[(k*l) % n]
            f[k, l] = result
            f[l, k] = result

    return f