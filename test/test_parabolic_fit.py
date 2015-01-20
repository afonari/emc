# -*- coding: utf-8 -*-
import os
import time
import unittest
import emc

# testing matrix operations
class EMC_Test(unittest.TestCase):
    #
    def setUp(self):
        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
        inpcar_path = os.path.join(script_dir, "Parabolic_Fit/INPCAR")
        self.inpcar_fh = open(inpcar_path, 'r')
    #
    def tearDown(self):
        self.inpcar_fh.close()
    #
    def test_parabolic_fit(self):
        print ""
        kpt, stepsize, band, prg, basis = emc.parse_inpcar(self.inpcar_fh, debug=False) # will need stepsize later
        st = []
        st.append([-5.0, -5.0, 0.0])
        st.append([-4.0, -4.0, 0.0])
        st.append([-3.0, -3.0, 0.0])
        st.append([-2.0, -2.0, 0.0])
        st.append([-1.0, -1.0, 0.0])
        st.append([0.0, 0.0, 0.0])
        st.append([1.0, 1.0, 0.0])
        st.append([2.0, 2.0, 0.0])
        st.append([3.0, 3.0, 0.0])
        st.append([4.0, 4.0, 0.0])
        st.append([5.0, 5.0, 0.0])
        #
        kpts = emc.generate_kpoints(kpt, st, stepsize, prg, basis)
        #print kpts
        kpoints_fh = open('KPOINTS', 'w')
        kpoints_fh.write("EMC \n")
        kpoints_fh.write("%d\n" % len(kpts))
        kpoints_fh.write("Reciprocal\n")
        #
        for i, kpt in enumerate(kpts):
            kpoints_fh.write( '%15.10f %15.10f %15.10f 0.01\n' % (kpt[0], kpt[1], kpt[2]) )
        #
        kpoints_fh.close()
#
if __name__ == '__main__':
    unittest.main()

