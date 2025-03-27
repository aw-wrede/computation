from typing import Literal

import numpy as np


def dft_matrix(n: int, norm: Literal['unscaled', 'ortho'] = 'unscaled', dtype:np.dtype = np.complex128) -> np.ndarray:
    """
    Constructs the DFT matrix of size n*n
    :param n: size of DFT matrix
    :param norm: Normalization mode, default is unscaled
    :param dtype: data type of the matrix
    :return: DFT matrix of size n*n
    """
    f = np.zeros((n, n), dtype=dtype)

    # Calculate omega (counter clockwise -> positive sign)
    omega = np.exp(-1j * 2 * np.pi / n)

    # Calculate all different unit roots
    omegas = np.zeros(n, dtype=dtype)

    # Normalize unit roots with 1/sqrt(n) if norm is ortho
    # Otherwise just calculate unit roots
    if norm == 'ortho':
        n_sqrt = np.sqrt(n)
        for i in range(n):
            omegas[i] = (omega ** i) / n_sqrt
    else:
        for i in range(n):
            omegas[i] = (omega ** i)

    # Fill matrix with unit roots
    for k in range(n):
        for l in range(k, n):
            result = omegas[(k*l) % n]
            f[k, l] = result
            f[l, k] = result

    return f

def dft(data: np.ndarray, norm: Literal['unscaled', 'ortho'] = 'unscaled', dtype: np.dtype = np.complex128) -> np.ndarray:
    """
    Performs the Discrete Fourier Transform on a real numpy array
    :param data: Real data to be transformed
    :param norm: Normalization mode, default is unscaled
    :param dtype: Data type of the fourier transformed data
    :return: Fourier transformed data
    """

    # Compute Omega
    o = dft_matrix(len(data), norm=norm, dtype=dtype)

    # Perform Discrete Fourier Transform (Omega * data)
    return np.dot(o, data)


def bin_mirror(number: int, length: int) -> int:
    """
    Mirrors a binary number with fixed length
    :param number: Number to be mirrored
    :param length: Fixed length
    :return: Mirrored number
    """

    if number == 0:
        return 0

    bin_number = format(number, f'0{length}b')
    bin_mirrored = bin_number[::-1]

    return int(bin_mirrored, 2)


def shuffle_bit_reversed_order(data: np.ndarray) -> np.ndarray:
    """
    Shuffles the  elements of data using bit-reversal of list index.
    :param data: Data to be transformed
    :return: Shuffled data array
    """

    n = data.shape[0]
    bit_length = (n - 1).bit_length()

    shuffled_data = np.zeros(n, data.dtype)

    # Shuffling by reversing index bits
    for i in range(n):
        mirrored = bin_mirror(i, bit_length)
        shuffled_data[mirrored] = data[i]

    return shuffled_data


def fft(data: np.ndarray, norm: Literal['unscaled', 'ortho'] = 'unscaled', dtype:np.dtype = np.complex128) -> np.ndarray:
    """
    Performs the Fourier Transform on a real numpy array.
    :param data: Real data to be transformed
    :param norm: Normalization mode, default is unscaled
    :param dtype: Data type of the fourier transformed data
    :return: Fourier transformed data
    """

    fdata = np.asarray(data, dtype=dtype)
    n = fdata.size

    # Check if array length is not zero, prevent possible power of 2 error
    if n == 0:
        raise ValueError('Data is empty')

    # Check if input length is power of two
    if (n & (n - 1)) != 0:
        return dft(fdata, norm=norm, dtype=dtype)

    # First step of FFT: shuffle data
    shuffled_data = shuffle_bit_reversed_order(fdata)

    # Second step, recursively merge transforms
    steps = int(np.log2(n))

    for step in range(steps):
        # unit square root
        k = 2 ** (step + 1)
        omega = np.exp(-1j * 2 * np.pi / k)

        # For each block
        for i in range(2 ** (steps - step - 1)):
            # Offset or index of first element in block
            block_start = i * 2 ** (step + 1)

            # For each pair of elements in block
            # (2^step is the amount of elements in each block)
            for j in range(2 ** step):
                j_0 = block_start + j               # index of first element to calculate
                j_1 = block_start + j + 2 ** step   # index of second element to calculate

                s_0 = shuffled_data[block_start + j] + (omega ** j) * shuffled_data[
                    2 ** step + block_start + j]

                s_1 = shuffled_data[block_start + j] + (omega ** j) * (-1) * shuffled_data[
                    2 ** step + block_start + j]

                shuffled_data[j_0] = s_0
                shuffled_data[j_1] = s_1


    # Normalize fft signal with 1/sqrt(n) if norm is ortho
    if norm == 'ortho':
        return shuffled_data / np.sqrt(n)

    return shuffled_data