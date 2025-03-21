import numpy as np

def dot(a:np.ndarray ,b:np.ndarray) -> np.ndarray:

    def get_partitions(a:np.ndarray) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
        """
        Split input matrix a into 4 partitions
        :param a: Input matrix
        :return: matrix a split into 4 partitions
        """
        # a = [ a11 a12
        #       a21 a22 ]
        h = int(a.shape[0] / 2)
        a11 = a[:h, :h]
        a12 = a[:h, h:]
        a21 = a[h:, :h]
        a22 = a[h:, h:]
        return a11, a12, a21, a22


    def from_partitions(a11: np.ndarray, a12:np.ndarray, a21:np.ndarray, a22:np.ndarray) -> np.ndarray:
        """
        Combine partitions into one matrix
        :param a11: partition
        :param a12: partition
        :param a21: partition
        :param a22: partition
        :return: Matrix
        """
        h = a11.shape[0]
        n = h + a21.shape[0]
        b = np.zeros((n, n), a11.dtype)
        b[:h, :h] = a11
        b[:h, h:] = a12
        b[h:, :h] = a21
        b[h:, h:] = a22

        return b


    def dot_general(a:np.ndarray, b:np.ndarray) -> np.ndarray:
        """
        Calculate dot product of two matrices with general formula
        :param a: Matrix
        :param b: Matrix
        :return: Dot product of two matrices
        """

        c = np.zeros((a.shape[0], b.shape[1]), a.dtype)

        for i in range(a.shape[0]):
            for j in range(b.shape[1]):
                val = 0
                for k in range(a.shape[1]):
                    val += a[i][k] * b[k][j]
                c[i][j] = val

        return c

    def dot_quadratic(a:np.ndarray , b:np.ndarray) -> np.ndarray:
        """
        Calculate dot product of two matrices n x n with n is even
        :param a: Matrix
        :param b: Matrix
        :return: dot product of two matrices
        """

        a11, a12, a21, a22 = get_partitions(a)
        b11, b12, b21, b22 = get_partitions(b)

        m1 = dot(a11 + a22 , b11 + b22)
        m2 = dot(a21 + a22,  b11)
        m3 = dot(a11,        b12 - b22)
        m4 = dot(a22,        b21 - b11)
        m5 = dot(a11 + a12,  b22)
        m6 = dot(a21 - a11,  b11 + b12)
        m7 = dot(a12 - a22,  b21 + b22)

        c11 = m1 + m4 - m5 + m7
        c12 = m3 + m5
        c21 = m2 + m4
        c22 = m1 - m2 + m3 + m6

        return from_partitions(c11, c12, c21, c22)


    m, n_a = a.shape
    n_b, p = b.shape

    # Check if matrices are compatible for multiplication
    if n_a != n_b:
        raise ValueError("Matrices must have compatible size")

    # Check if matrix is quadratic and size is even
    # Use partition method then to calculate dot product
    if a.shape == b.shape and m == n_a and m % 2 == 0:
        return dot_quadratic(a, b)

    return dot_general(a, b)

def close(a:np.ndarray, b:np.ndarray, rtol=1e-05, atol=1e-08) -> np.ndarray:
    """
    Checks whether the individual elements of both matrices are close to each other and returns the result per element
    :param a: first matrix
    :param b: second matrix, used as reference
    :param rtol: relative tolerance, which depends on the value b
    :param atol: absolute tolerance
    :return: matrix containing the result of each element
    """
    if a.shape != b.shape:
        raise ValueError("Matrices must have compatible size")

    results = (np.abs(a - b) <= atol + rtol * abs(b))

    return results

def close_all(a:np.ndarray, b:np.ndarray, rtol=1e-05, atol=1e-08) -> bool:
    """
        Checks whether the individual elements of both matrices are close to each other
        :param a: first matrix
        :param b: second matrix, used as reference
        :param rtol: relative tolerance, which depends on the value b
        :param atol: absolute tolerance
        :return: result of whether the matrices are close to each other
        """
    if a.shape != b.shape:
        raise ValueError("Matrices must have compatible size")

    return np.all(np.abs(a - b) <= (atol + rtol * np.abs(b)))