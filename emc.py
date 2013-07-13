#!/usr/bin/env python

EMC_VERSION='2.0b1'

diff_d2 = []
diff_d2.append([0.0, 0.0, 0.0]); # 0
diff_d2.append([-1.0, 0.0, 0.0]); diff_d2.append([1.0, 0.0, 0.0]);  # dx  1-2
diff_d2.append([0.0, -1.0, 0.0]); diff_d2.append([0.0, 1.0, 0.0])   # dy  3-4
diff_d2.append([0.0, 0.0, -1.0]); diff_d2.append([0.0, 0.0, 1.0])   # dz  5-6
diff_d2.append([-1.0, -1.0, 0.0]); diff_d2.append([1.0, 1.0, 0.0]); diff_d2.append([1.0, -1.0, 0.0]); diff_d2.append([-1.0, 1.0, 0.0]); # dxdy 7-10
diff_d2.append([-1.0, 0.0, -1.0]); diff_d2.append([1.0, 0.0, 1.0]); diff_d2.append([1.0, 0.0, -1.0]); diff_d2.append([-1.0, 0.0, 1.0]); # dxdz 11-14
diff_d2.append([0.0, -1.0, -1.0]); diff_d2.append([0.0, 1.0, 1.0]); diff_d2.append([0.0, 1.0, -1.0]); diff_d2.append([0.0, -1.0, 1.0]); # dydz 15-18

Bohr = 0.5291772

def fd_effmass(e, stepsize, order=2):
    m = [[0.0 for i in range(3)] for j in range(3)]
    m[0][0] = (e[1] - 2.0*e[0] + e[2])/stepsize**2
    m[1][1] = (e[3] - 2.0*e[0] + e[4])/stepsize**2
    m[2][2] = (e[5] - 2.0*e[0] + e[6])/stepsize**2

    m[0][1] = (e[7] + e[8] - e[9] - e[10])/(4*stepsize**2)
    m[0][2] = (e[11] + e[12] - e[13] - e[14])/(4*stepsize**2)
    m[1][2] = (e[15] + e[16] - e[17] - e[18])/(4*stepsize**2)

    # symmetrize
    m[1][0] = m[0][1]
    m[2][0] = m[0][2]
    m[2][1] = m[1][2]

    return m

def get_kpoints(kpt, stepsize, prg, basis, debug=False):
    import numpy as np
    import sys

    basis_r_np = np.linalg.inv(np.array(basis).T)* 2*np.pi
    #else:
    #    basis_r_np = np.linalg.inv(np.array(basis))

    kpt_r_cart = np.dot(np.array(kpt), basis_r_np)
    if debug: print kpt_r_cart

    if prg == 'V':
        stepsize = stepsize*(1/Bohr) # [1/A]

    kpoints = []
    for i in range(len(diff_d2)):
        k_cart = kpt_r_cart + np.array(diff_d2[i])*stepsize
        k_frac = np.dot(k_cart, np.linalg.inv(basis_r_np))
        kpoints.append( [k_frac[0], k_frac[1], k_frac[2]] )

    return kpoints

def parse_EIGENVAL_VASP(eigenval_fh, band, diff2_size):
    import sys
    import re

    ev2h = 1.0/27.21138505
    eigenval_fh.seek(0) # just in case
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()

    nelec, nkpt, nband = [int(s) for s in eigenval_fh.readline().split()]
    print 'Reading __filename__:\n'
    print '  Number of the valence band is %d (NELECT/2)' % (nelec/2)
    if band > nband:
        print 'Requested band (%d) is larger than total number of the calculated bands (%d), exiting...' % (band, nband)
        sys.exit(0)

    energies = []
    for i in range(diff2_size):
        eigenval_fh.readline() # empty line
        eigenval_fh.readline() # k point coordinates
        for j in range(1, nband+1):
            line = eigenval_fh.readline()
            if band == j:
                energies.append(float(line.split()[1])*ev2h)

    return energies

def parse_inpcar(filename="INPCAR", debug=False):
    import sys
    import re

    kpt = []       # k-point at which eff. mass in reciprocal reduced coords (3 floats)
    stepsize = 0.0 # stepsize for finite difference (1 float) in Bohr
    band = 0       # band for which eff. mass is computed (1 int)
    prg = ''       # program identifier (1 char)
    basis = []     # basis vectors in cartesian coords (3x3 floats), units depend on the program identifier

    inpcar = open(filename, 'r')

    p = re.search(r'^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', inpcar.readline())
    if p:
        kpt = [float(p.group(1)), float(p.group(2)), float(p.group(3))]
        print "Found k point in the reduced reciprocal space: %5.3f %5.3f %5.3f" % (kpt[0], kpt[1], kpt[2])
    else:
        print "Was expecting k point on the line 0 (3 floats), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\d+\.\d+)', inpcar.readline())
    if p:
        stepsize = float(p.group(1))
        print "Found stepsize of: %5.3f (1/Bohr)" % stepsize
    else:
        print "Was expecting a stepsize on line 1 (1 float), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\d+)', inpcar.readline())
    if p:
        band = int(p.group(1))
        print "Requested band is : %5d" % band
    else:
        print "Was expecting band number on line 2 (1 int), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\w)', inpcar.readline())
    if p:
        prg = p.group(1)
        print "Program identifier is: %5c" % prg
    else:
        print "Was expecting program identifier on line 3 (1 char), didn't get it, exiting..."
        sys.exit(1)

    for i in range(3):
        p = re.search(r'^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', inpcar.readline())
        if p:
            basis.append([float(p.group(1)), float(p.group(2)), float(p.group(3))])

    print "Basis:"
    print basis

    return kpt, stepsize, band, prg, basis

if __name__ == '__main__':
    import sys
    import re
    import numpy as np

    kpt, stepsize, band, prg, basis = parse_inpcar()

    output_filename = ''
    if len(sys.argv) > 1:
        output_filename = sys.argv[1]
    else:
        if prg.upper() == 'V':
            output_filename = 'EIGENVAL'
        else:
            output_filename = "OUTCAR"

    output_fh = 0
    output_exists = False
    try:
        output_fh = open(output_filename, 'r')
        output_exists = True
    except IOError:
        pass

    energies = []
    if output_exists:
        print 'Successfully opened '+output_filename+', preparing to parse it...\n'

        if prg.upper() == 'V':
            energies = parse_EIGENVAL_VASP(output_fh, band, len(diff_d2))
            m = fd_effmass(energies, stepsize)

        print 'Original matrix:'
        for i in range(len(m)):
            print '%7.5f %7.5f %7.5f' % (m[i][0], m[i][1], m[i][2])
        print '\n'

        eigval, eigvec = np.linalg.eigh(np.array(m))
        print 'Principle effective masses and directions:\n'
        for i in range(len(m)):
            eigenvec_norm = eigvec[:,i]/np.max(eigvec[:,i])
            print 'Effective mass (%d): %12.3f' % (i, 1.0/eigval[i])
            print 'Cartesian coordinates: %7.5f %7.5f %7.5f\n' % (eigenvec_norm[0], eigenvec_norm[1], eigenvec_norm[2])

    else:
        print 'No '+output_filename+' file found, entering the Generation regime...\n'
        kpoints = get_kpoints(kpt, stepsize, prg, basis, debug=False)

        kpoints_fh = open('KPOINTS', 'w')
        kpoints_fh.write("EMC "+EMC_VERSION+"\n")
        kpoints_fh.write("%d\n" % len(diff_d2))
        kpoints_fh.write("Reciprocal\n")

        for i, kpt in enumerate(kpoints):
            kpoints_fh.write( '%7.5f %7.5f %7.5f 0.01\n' % (kpt[0], kpt[1], kpt[2]) )

        kpoints_fh.close()
        print 'KPOINTS file has been generated in the current directory...\n'






