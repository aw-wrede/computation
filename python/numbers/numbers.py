import numpy as np

def close(a:np.number, b:np.number, rtol=1e-05, atol=1e-08):
    """
    Checks whether two given numbers are close to each other within the given tolerances
    :param a: First number
    :param b: Second number, used as reference
    :param rtol: Relative tolerance, which depends on the value b
    :param atol: Absolute tolerance
    :return: Whether the numbers are close to each other
    """
    return abs(a - b) <= atol + rtol * abs(b)