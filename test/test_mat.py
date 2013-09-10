# -*- coding: utf-8 -*-
import os
import time
import unittest
import emc

# testing matrix operations
class EMC_Test(unittest.TestCase):
    def assertListAlmostEqual(self, list1, list2, places, msg):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(float(a), float(b), places, msg)

    def test_transpose(self):
        m = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]
        m_t = emc.T(m)
        for i in range(len(m)):
            row = [ m[j][i] for j in range(len(m)) ]
            self.assertListAlmostEqual(m_t[i], row, places=5, msg='Matrix transpose failed for row %d' % i)

    def test_mat_n_vec(self):
        m = [[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]
        v = [-1., -2., -3.]
        p = emc.MAT_m_VEC(m, v)
        self.assertListAlmostEqual(p, [-14., -32., -50.], places=5, msg='MAT_n_VEC failed')

    def test_normal(self):
        m = [-1., 2., -3.]
        m_n = emc.N(m)
        self.assertListAlmostEqual(m_n, [ item/(-3.) for item in m ], places=5, msg='Normal function failed')

    def test_det(self):
        m = [[-1., 2., 3.], [4., 5., -6.], [7., -8., 9.]]
        det = emc.DET_3X3(m)
        self.assertAlmostEqual(det, -354., places=2, msg='Determinant failed')

    def test_eigenval_diag(self):
        m_diag = [[9., 0., 0.], [0., -4., 0.], [0., 0., 1.]]
        eigvec, eigvals = emc.jacobi(m_diag)
        self.assertListAlmostEqual(eigvals, [-4., 1., 9.], places=1, msg='Eigenvals for diagonal symmetric matrix failed')

    def test_eigenval(self):
        m = [[1.,2.,30.],[2.,5.,10.],[30.,10.,9.]]
        eigvec, eigvals = emc.jacobi(m)
        self.assertListAlmostEqual(eigvals, [-26.13065, 3.443361, 37.68729], places=5, msg='Eigenvals for symmetric matrix failed')

if __name__ == '__main__':
    unittest.main()

