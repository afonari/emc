# -*- coding: utf-8 -*-
import os
import time
import unittest
import emc

# testing cubic system
class EMC_Test(unittest.TestCase):
    def setUp(self):
        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
        inpcar_path = os.path.join(script_dir, "GaAs-5.648-PBE-VASP/INPCAR")
        self.inpcar_fh = open(inpcar_path, 'r')

        eigenval_path = os.path.join(script_dir, "GaAs-5.648-PBE-VASP/EIGENVAL")
        self.eigenval_fh = open(eigenval_path, 'r')

    def test_parse_inpcar(self):
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need stepsize later
        self.assertListEqual(kpt, [0.0, 0.0, 0.0], msg='Failed to parse K-point')
        self.assertEquals(stepsize, 0.01, msg='Failed to parse stepsize')
        self.assertEquals(band, 16, msg='Failed to parse band')
        self.assertEquals(prg, 'V', msg='Failed to parse program identifier')
        self.assertListEqual(basis, [[5.648, 0.0, 0.0],[0.0, 5.648, 0.0],[0.0, 0.0, 5.648]], msg='Failed to parse basis')

    def test_calculate_effmass(self):
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need stepsize later

        energies = emc.parse_EIGENVAL_VASP(self.eigenval_fh, band, len(emc.diff_d2), debug=False)
        #m = emc.fd_effmass(energies, stepsize)
        #self.assertListEqual(m, [[-2.90687, 0.0, 0.0],[0.0, -2.90687, 0.0],[0.0, 0.0, -2.90687]], msg='Failed to parse basis')

if __name__ == '__main__':
    unittest.main()

