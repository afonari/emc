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

def MAT_m_VEC(m, v):
    p = [ 0.0 for i in range(len(v)) ]
    for i in range(len(m)):
        assert len(v) == len(m[i]), 'Length of the matrix row is not equal to the length of the vector'
        p[i] = sum( [ m[i][j]*v[j] for j in range(len(v)) ] )
    return p

def T(m):
    p = [[ m[i][j] for i in range(len( m[j] )) ] for j in range(len( m )) ]
    return p

def N(v):
    max_ = 0.
    for item in v:
        if abs(item) > abs(max_): max_ = item

    return [ item/max_ for item in v ]
def DET_3X3(m):
    assert len(m) == 3, 'Matrix should be of the size 3 by 3'
    return m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2] - \
           m[0][2]*m[1][1]*m[2][0] - m[2][1]*m[1][2]*m[0][0] - m[2][2]*m[0][1]*m[1][0]

def IS_SYMMETRIC(m):
    for i in range(len(m)):
        for j in range(len(m[i])):
            if m[i][j] != m[j][i]: return False # automatically checks square-shape

    return True

def EIGVALUES_SYMMETRIC_3X3(m):
    # follows http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    import math

    assert IS_SYMMETRIC(m), 'Supplied matrix is not symmetric'

    eigs = [ 0.0 for i in range(3) ]
    n = [[0.0 for i in range(3)] for j in range(3)]
    identity = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]

    p = m[0][1]**2 + m[0][2]**2 + m[1][2]**2
    if p == 0.0:
        # m is diagonal.
        eigs = sorted( [ m[i][i] for i in range(3) ] ) # ascending order
    else:
        q = sum([ m[i][i] for i in range(3) ])/3.0
        p = (m[0][0] - q)**2 + (m[1][1] - q)**2 + (m[2][2] - q)**2 + 2.0*p
        p = math.sqrt(p/6.0)

        for i in range(3):
            for j in range(3):
                n[i][j] = (1.0/p) * (m[i][j] - q*identity[i][j])

        r = DET_3X3(n)/2.0

        if r <= -1.:
            phi = math.pi/3.0
        elif r >= 1.:
            phi = 0.0
        else:
            phi = math.acos(r)/3.0

        eigs[0] = q + 2.0*p*math.cos(phi)
        eigs[1] = q + 2.0*p*math.cos(phi + math.pi * (2.0/3.0))
        eigs[2] = 3.0*q - eigs[0] - eigs[1]

        eigs = sorted(eigs) # I'm aware of the list.sort() method

    return eigs

def cart2frac(basis, v):
    import numpy as np
    return np.dot(v, np.linalg.inv(basis)).tolist()

def fd_effmass(e, stepsize, order=2, debug=False):
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

    if debug:
        print 'Assembling effective mass tensor...'
        for i in range(len(m)):
            print '%7.5f %7.5f %7.5f' % (m[i][0], m[i][1], m[i][2])

    if debug: print ''
    return m

def generate_kpoints(kpt_frac, stepsize, prg, basis, debug=False):
    import numpy as np
    import sys

    # working in the reciprocal space
    basis_r = (np.linalg.inv( T(basis) )* 2*np.pi).tolist()

    k_c = MAT_m_VEC(T(basis_r), kpt_frac)
    if debug: print kpt_r_cart

    if prg == 'V':
        stepsize = stepsize*(1/Bohr) # [1/A]

    kpoints = []
    for i in range(len(diff_d2)):
        k_c_ = [ k_c[j] + diff_d2[i][j]*stepsize for j in range(3) ] # getting displaced k points in Cartesian coordinates
        k_f = cart2frac(basis_r, k_c_)
        kpoints.append( [k_f[0], k_f[1], k_f[2]] )

    return kpoints

def parse_EIGENVAL_VASP(eigenval_fh, band, diff2_size, debug=False):
    import sys

    ev2h = 1.0/27.21138505
    eigenval_fh.seek(0) # just in case
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()
    eigenval_fh.readline()

    nelec, nkpt, nband = [int(s) for s in eigenval_fh.readline().split()]
    if debug: print 'From EIGENVAL: Number of the valence band is %d (NELECT/2)' % (nelec/2)
    if band > nband:
        print 'Requested band (%d) is larger than total number of the calculated bands (%d), exiting...' % (band, nband)
        sys.exit(1)

    energies = []
    for i in range(diff2_size):
        eigenval_fh.readline() # empty line
        eigenval_fh.readline() # k point coordinates
        for j in range(1, nband+1):
            line = eigenval_fh.readline()
            if band == j:
                energies.append(float(line.split()[1])*ev2h)

    if debug: print ''
    return energies

def parse_inpcar(inpcar_fh, debug=False):
    import sys
    import re

    kpt = []       # k-point at which eff. mass in reciprocal reduced coords (3 floats)
    stepsize = 0.0 # stepsize for finite difference (1 float) in Bohr
    band = 0       # band for which eff. mass is computed (1 int)
    prg = ''       # program identifier (1 char)
    basis = []     # basis vectors in cartesian coords (3x3 floats), units depend on the program identifier

    inpcar_fh.seek(0) # just in case
    p = re.search(r'^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', inpcar_fh.readline())
    if p:
        kpt = [float(p.group(1)), float(p.group(2)), float(p.group(3))]
        if debug: print "Found k point in the reduced reciprocal space: %5.3f %5.3f %5.3f" % (kpt[0], kpt[1], kpt[2])
    else:
        print "Was expecting k point on the line 0 (3 floats), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\d+\.\d+)', inpcar_fh.readline())
    if p:
        stepsize = float(p.group(1))
        if debug: print "Found stepsize of: %5.3f (1/Bohr)" % stepsize
    else:
        print "Was expecting a stepsize on line 1 (1 float), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\d+)', inpcar_fh.readline())
    if p:
        band = int(p.group(1))
        if debug: print "Requested band is : %5d" % band
    else:
        print "Was expecting band number on line 2 (1 int), didn't get it, exiting..."
        sys.exit(1)

    p = re.search(r'^\s*(\w)', inpcar_fh.readline())
    if p:
        prg = p.group(1)
        if debug: print "Program identifier is: %5c" % prg
    else:
        print "Was expecting program identifier on line 3 (1 char), didn't get it, exiting..."
        sys.exit(1)

    for i in range(3):
        p = re.search(r'^\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', inpcar_fh.readline())
        if p:
            basis.append([float(p.group(1)), float(p.group(2)), float(p.group(3))])

    if debug: 
        print "Real space basis:"
        for i in range(len(basis)):
            print '%9.7f %9.7f %9.7f' % (basis[i][0], basis[i][1], basis[i][2])

    if debug: print ''

    return kpt, stepsize, band, prg, basis

if __name__ == '__main__':
    import sys
    import re
    import datetime
    import numpy as np

    print '\nEffective mass calculator '+EMC_VERSION
    print 'License: MIT'
    print 'Developed by: Alexandr Fonari'
    print 'Started at: '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+'\n'

    inpcar_fh = 0
    try:
        inpcar_fh = open('INPCAR', 'r')
    except IOError:
        print "Couldn't open INPCAR file, exiting...\n"
        sys.exit(1)

    kpt, stepsize, band, prg, basis = parse_inpcar(inpcar_fh)

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

    if output_exists:
        print 'Successfully opened '+output_filename+', preparing to parse it...\n'

        energies = []
        if prg.upper() == 'V' or prg.upper() == 'C':
            energies = parse_EIGENVAL_VASP(output_fh, band, len(diff_d2))
            m = fd_effmass(energies, stepsize)

        eigval, eigvec = np.linalg.eigh(np.array(m))
        print 'Principle effective masses and directions:\n'
        for i in range(len(m)):
            vec_cart = eigvec[:,i].tolist()
            vec_frac = cart2frac(basis, vec_cart)
            #vec_real = MAT_m_VEC(T(basis), vec)
            vec_n = N(vec_frac)
            print 'Effective mass (%d): %12.3f' % (i, 1.0/eigval[i])
            print 'Original eigenvectors: %7.5f %7.5f %7.5f\n' % (vec_cart[0], vec_cart[1], vec_cart[2])
            print 'Normal fractional coordinates: %7.5f %7.5f %7.5f\n' % (vec_n[0], vec_n[1], vec_n[2])

    else:
        print 'No '+output_filename+' file found, entering the Generation regime...\n'

        if prg.upper() == 'C' and band != 1:
            print 'Band should be set to 1 for CRYSTAL calculations,'
            print 'desired band number is set as a parameter (-b) for cry-getE.pl script.\n'
            sys.exit(1)

        kpoints = generate_kpoints(kpt, stepsize, prg, basis, debug=False)

        kpoints_fh = open('KPOINTS', 'w')
        kpoints_fh.write("EMC "+EMC_VERSION+"\n")
        kpoints_fh.write("%d\n" % len(diff_d2))
        kpoints_fh.write("Reciprocal\n")

        for i, kpt in enumerate(kpoints):
            kpoints_fh.write( '%7.5f %7.5f %7.5f 0.01\n' % (kpt[0], kpt[1], kpt[2]) )

        kpoints_fh.close()
        print 'KPOINTS file has been generated in the current directory...\n'






