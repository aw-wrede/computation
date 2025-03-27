import numpy as np


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
