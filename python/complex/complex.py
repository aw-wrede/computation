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


def is_unitary(m: np.ndarray) -> bool:
    """
    Checks if a matrix is unitary
    :param m: Matrix to be checked
    :return: Whether the matrix is unitary
    """

    # Check if matrix is quadratic
    if m.shape[0] != m.shape[1]:
        return False

    # Calculate adjoint (conjugate transpose) of matrix
    m_a = np.conjugate(m).T

    # The product of the matrix with its adjoint gives the unit matrix if the matrix is unitary
    return np.allclose(np.dot(m, m_a), np.eye(m.shape[0]))