import unittest

import numpy as np

import complex as c


class MyTestCase(unittest.TestCase):
    f_4 = c.fft.dft_matrix(4, norm='ortho')
    f_5 = c.fft.dft_matrix(5, norm='ortho')
    f_10 = c.fft.dft_matrix(10, norm='ortho')

    eye_4 = np.eye(4)

    rnd_matrix = np.random.randn(4, 3)

    def test_dft_matrix(self):
        self.assertTrue(np.allclose(self.f_4, np.fft.fft(np.eye(4), norm='ortho')))
        self.assertTrue(np.allclose(self.f_5, np.fft.fft(np.eye(5), norm='ortho')))
        self.assertTrue(np.allclose(self.f_10, np.fft.fft(np.eye(10), norm='ortho')))

    def test_unitary(self):
        self.assertTrue(c.is_unitary(self.f_4))
        self.assertTrue(c.is_unitary(self.f_5))
        self.assertTrue(c.is_unitary(self.f_10))

        self.assertTrue(c.is_unitary(self.eye_4))
        self.assertFalse(c.is_unitary(self.rnd_matrix))

    def test_fft(self):
        random = np.random.randn(64)
        self.assertTrue(np.allclose(c.fft.fft(random), np.fft.fft(random)))
        self.assertTrue(np.allclose(c.fft.fft(random, norm='ortho'), np.fft.fft(random, norm='ortho')))

    def test_dft(self):
        random = np.random.randn(63)
        self.assertTrue(np.allclose(c.fft.dft(random), np.fft.fft(random)))
        self.assertTrue(np.allclose(c.fft.dft(random, norm='ortho'), np.fft.fft(random, norm='ortho')))

if __name__ == '__main__':
    unittest.main()
