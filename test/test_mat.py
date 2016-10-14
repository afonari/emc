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

    def test_get_eff_masses(self):
        m = [[9.,5.,1.], [5.,7.,3.], [1.,3.,6.]]
        basis = [[6.,2.,1.], [3.,9.,4.], [1.,3.,8.]]
        em, vecs_cart, vecs_frac, vecs_n = emc.get_eff_masses(m, basis)
        self.assertListAlmostEqual(em, [0.48965,0.16844,0.07132], places=5, msg='Effective mass value failed')
        self.assertListAlmostEqual(vecs_cart[0], [-0.47715, 0.75420,-0.45113], places=5, msg='Cartesian eigenvectors failed')
        self.assertListAlmostEqual(vecs_cart[1], [-0.53266, 0.16010, 0.83104], places=5, msg='Cartesian eigenvectors failed')
        self.assertListAlmostEqual(vecs_cart[2], [ 0.69900, 0.63683, 0.32534], places=5, msg='Cartesian eigenvectors failed')
        self.assertListAlmostEqual(vecs_frac[0], [-0.13660, 0.15271,-0.11567], places=5, msg='Fractional eigenvectors failed')
        self.assertListAlmostEqual(vecs_frac[1], [-0.10988, 0.00360, 0.11581], places=5, msg='Fractional eigenvectors failed')
        self.assertListAlmostEqual(vecs_frac[2], [ 0.09126, 0.04887, 0.00483], places=5, msg='Fractional eigenvectors failed')
        self.assertListAlmostEqual(vecs_n[0], [-0.89450, 1.00000,-0.75745], places=5, msg='Normalized eigenvectors failed')
        self.assertListAlmostEqual(vecs_n[1], [-0.94877, 0.03110, 1.00000], places=5, msg='Normalized eigenvectors failed')
        self.assertListAlmostEqual(vecs_n[2], [ 1.00000, 0.53551, 0.05287], places=5, msg='Normalized eigenvectors failed')

if __name__ == '__main__':
    unittest.main()

