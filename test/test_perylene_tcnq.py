# -*- coding: utf-8 -*-
import os
import time
import unittest
import emc

# testing cubic system
class EMC_Test(unittest.TestCase):
    def assertListAlmostEqual(self, list1, list2, places, msg):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(float(a), float(b), places, msg)

    def setUp(self):
        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
        inpcar_path = os.path.join(script_dir, "Perylene-TCNQ-B3LYP-6-21G-CRY09/INPCAR")
        self.inpcar_fh = open(inpcar_path, 'r')

        eigenval_path = os.path.join(script_dir, "Perylene-TCNQ-B3LYP-6-21G-CRY09/EIGENVAL")
        self.eigenval_fh = open(eigenval_path, 'r')

        kpoints_path = os.path.join(script_dir, "Perylene-TCNQ-B3LYP-6-21G-CRY09/KPOINTS")
        self.kpoints_fh = open(kpoints_path, 'r')

    def tearDown(self):
        self.inpcar_fh.close()
        self.eigenval_fh.close()
        self.kpoints_fh.close()

    def test_parse_inpcar(self):
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need stepsize later
        self.assertListEqual(kpt, [0.0, 0.0, 0.0], msg='Failed to parse K-point')
        self.assertEquals(stepsize, 0.01, msg='Failed to parse stepsize')
        self.assertEquals(band, 1, msg='Failed to parse band')
        self.assertEquals(prg, 'C', msg='Failed to parse program identifier')
        self.assertListEqual(basis, [[13.8324582, -0.0965703, 0.0],[0.0, 27.4955152, 0.0],[0.0, 0.0, 20.5602203]], msg='Failed to parse basis')

    def test_calculate_effmass(self):
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need band and stepsize later

        energies = emc.parse_EIGENVAL_VASP(self.eigenval_fh, band, len(emc.st3), debug=False)
        m = emc.fd_effmass_st3(energies, stepsize)
        self.assertListAlmostEqual(m[0], [-1.23535, 0.04464, 0.0], places=5, msg='Failed to calculate effective mass tensor')
        self.assertListAlmostEqual(m[1], [0.04464, -0.45875, 0.0], places=5, msg='Failed to calculate effective mass tensor')
        self.assertListAlmostEqual(m[2], [0.0, 0.0, -0.67875], places=5, msg='Failed to calculate effective mass tensor')

        print ''
        print 'Effective mass values and directions for Perylene-TCNQ crystal'
        print 'are tested against data computed in the JLB group'
        print ''
        masses, vecs_cart, vecs_frac, vecs_n = emc.get_eff_masses(m, basis)
        self.assertListAlmostEqual(masses, [-0.808, -1.473, -2.192], places=3, msg='Effective mass calculations failed')
        self.assertListAlmostEqual(vecs_n[0], [1.0, -0.02531, 0.0], places=5, msg='Effective mass direction 1 failed')
        self.assertListAlmostEqual(vecs_n[1], [0.0, 0.0, 1.0], places=5, msg='Effective mass direction 2 failed')
        self.assertListAlmostEqual(vecs_n[2], [0.11385, 1.0, 0.0], places=5, msg='Effective mass direction 3 failed')

    def test_kpoints(self):
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need stepsize later
        kpts = emc.generate_kpoints(kpt, emc.st3, stepsize, prg, basis)

        self.kpoints_fh.readline() # title
        nkpt = int(self.kpoints_fh.readline()) # Reciprocal
        self.assertEquals(nkpt, len(kpts), msg='Length of the list is not equal to the number from KPOINTS')
        self.kpoints_fh.readline() # Reciprocal

        for i in range(len(kpts)):
            kp = [float(x) for x in self.kpoints_fh.readline().split()[0:3]] # kx ky kz w
            self.assertListAlmostEqual(kpts[i], kp, places=5, msg='K-point %d is not equal to that from the KPOINTS file' % i)

if __name__ == '__main__':
    unittest.main()

