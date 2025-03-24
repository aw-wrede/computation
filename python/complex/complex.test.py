import unittest

import complex as c


class MyTestCase(unittest.TestCase):
    def test_something(self):
        f_4 = c.dft_matrix(4)

        # TODO add test cases
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
