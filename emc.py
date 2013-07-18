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

# def get_eigensystem_3x3(m): # SHOULD BE SYMMETRIC!
#     # en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
#     eigvals = [0.0 for i in range(3)]
#     b = [[0.0 for i in range(3)] for j in range(3)]
#     identity = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
# 
#     p = m[0][1]**2 + m[0][2]**2 + m[1][2]**2
#     if (p == 0.0):
#         # m is diagonal
#         eigvals[0] = m[0][0]
#         eigvals[1] = m[1][1]
#         eigvals[2] = m[2][2]
#     else:
#         q = (m[0][0] + m[1][1] + m[2][2])/3
#         p = (m[0][0] - q)**2 + (m[1][1] - q)**2 + (m[2][2] - q)**2 + 2*p
#         p = sqrt(p/6.0)
# 
#         # B = (1 / p) * (A - q * I)    I is the identity matrix
#         for i in range(3):
#             for j in range(3):
#                 b[i][j] = (1.0/p) * (m[i][j] - q*identity[i][j])
#                 r = det(B)/2.0
#  
#    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
#    % but computation error can leave it slightly outside this range.
#    if (r <= -1) 
#       phi = pi / 3
#    elseif (r >= 1)
#       phi = 0
#    else
#       phi = acos(r) / 3
#    end
#  
#    % the eigenvalues satisfy eig3 <= eig2 <= eig1
#    eig1 = q + 2 * p * cos(phi)
#    eig3 = q + 2 * p * cos(phi + pi * (2/3))
#    eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
# end
# 
# def det_3by3(m):
#     return m[0][0]*m[1][1]*m[2][2] + m[0][1]*$g[6]*$g[7] + $g[3]*$g[8]*$g[4] - $g[7]*$g[5]*$g[3] - $g[4]*$g[2]*$g[9] - $g[1]*$g[6]*$g[8]; 

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

def MAT_m_VEC(m, v):
    p = [ 0.0 for i in range(len(v)) ]
    for i in range(len(m)):
        assert len(v) == len(m[i]), 'Length of the matrix row is not equal to the length of the vector'
        p[i] = sum( [ m[i][j]*v[j] for j in range(len(v)) ] )
    return p

def T(m):
    p = [[ m[i][j] for i in range(len( m[j] )) ] for j in range(len( m )) ]
    return p
    
def generate_kpoints(kpt, stepsize, prg, basis, debug=False):
    import numpy as np
    import sys

    basis_r_np = (np.linalg.inv(np.array(T(basis)))* 2*np.pi).tolist()

    kpt_r_cart = MAT_m_VEC(T(basis_r_np), kpt)
    if debug: print kpt_r_cart

    if prg == 'V':
        stepsize = stepsize*(1/Bohr) # [1/A]

    kpoints = []
    for i in range(len(diff_d2)):
        k_cart = kpt_r_cart + np.array(diff_d2[i])*stepsize
        k_frac = np.dot(k_cart, np.linalg.inv(basis_r_np))
        kpoints.append( [k_frac[0], k_frac[1], k_frac[2]] )

    return kpoints

def parse_EIGENVAL_VASP(eigenval_fh, band, diff2_size, debug=False):
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
            print m

        eigval, eigvec = np.linalg.eigh(np.array(m))
        print 'Principle effective masses and directions:\n'
        for i in range(len(m)):
            max_i = np.argmax(np.abs(eigvec[:,i]))
            eigenvec_norm = eigvec[:,i]/eigvec[max_i,i]
            print 'Effective mass (%d): %12.3f' % (i, 1.0/eigval[i])
            print 'Cartesian coordinates: %7.5f %7.5f %7.5f\n' % (eigenvec_norm[0], eigenvec_norm[1], eigenvec_norm[2])

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






