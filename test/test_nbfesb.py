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

    def test_calculate_effmass(self):
        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in

        # See NbFeSb-CASTEP/emcpy.out_electrons_band19
        inpcar_fn = os.path.join(script_dir, 'NbFeSb-CASTEP', 'input_el')
        inpcar_fh = open(inpcar_fn, 'r')
        eigenval_fn = os.path.join(script_dir, 'NbFeSb-CASTEP', 'NbFeSb.electrons.bands')
        eigenval_fh = open(eigenval_fn, 'r')

        kpt, stepsize, band, prg, basis = emc.parse_inpcar(inpcar_fh, debug=False) # will need band and stepsize later

        energies = emc.parse_bands_CASTEP(eigenval_fh, band, len(emc.st3), debug=False)
        m = emc.fd_effmass_st3(energies, stepsize)
        self.assertListAlmostEqual(m[0], [3.459576, 0.0, 0.0], places=5, msg='Failed to calculate effective mass tensor')
        self.assertListAlmostEqual(m[1], [0.0, 3.459576, 0.0], places=5, msg='Failed to calculate effective mass tensor')
        self.assertListAlmostEqual(m[2], [0.0, 0.0, 1.759904], places=5, msg='Failed to calculate effective mass tensor')

if __name__ == '__main__':
    unittest.main()

