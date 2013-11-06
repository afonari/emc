! Copyright (c) "2012, Alexandr Fonari
!                URL: https://github.com/alexandr-fonari/emc
!                License: MIT License"

program emc_calc ! version 1.0
use emc_functions
implicit none
real(kind=8), parameter :: b2a = 0.52917721092d0
real(kind=8), parameter :: pi = 4.d0*DATAN(1.d0)
real(kind=8), parameter :: ev2h = 1.d0/27.21138505d0
integer(kind=4), parameter :: nkpoints = 61
integer(kind=4), parameter :: iunt = 10, ilog = 11
character(len=8), parameter :: version_number = '1.5f'

real(kind=8) :: get_next_eigeval, get_2nd_deriv, get_mixed_deriv1
real(kind=8) :: f(3,3), g(3,3), ev(3,3), v(3), len
real(kind=8) :: E(-2:2,-2:2,-2:2), m(3,3), b(3), WORK(100)
real(kind=8) :: kp(3), dk
integer i, i1, j, j1, ok
integer(kind=4) :: itrash, count, nbands, band, A(4)
character(len=1) :: prg
external :: DSYEV

open(unit=ilog,file='emc_calc.log',form='formatted')
write(ilog,*) "Effective Mass Calculator calculator ", version_number
call print_time(ilog)
write(ilog,*)

! read input ########################################################

open(unit=iunt,file='inp',form='formatted')
    read(iunt,fmt=*) (kp(j),j=1,3)
    read(iunt,fmt=*) dk
    read(iunt,fmt=*) band
    read(iunt,fmt=*) prg
    read(iunt,*) ((f(i,j),j=1,3),i=1,3)
close(iunt)

if(prg .eq. 'C') then
    write(ilog,*) "It is Crystal - setting band to 1"
    band=1
end if

! g = inverse(transpose(f)), will do it in two steps
g=f
g=transpose(g)
call inverse(g, 3)
g=g*2.d0*pi

write(ilog,*) "direct lattice vectors              reciprocal lattice vectors"
do i=1,3
    write(ilog,"(3F10.6,A,3F10.6)") , (f(i,j), j=1,3), "     ", (g(i,j), j=1,3)
end do
write(ilog,*)

! read EIGENVAL ###########################################

open(unit=iunt,file='EIGENVAL',form='formatted')

read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*) itrash, itrash, nbands

write(ilog,"(A,2I5,F10.6)") "nbands, band, dk: ", nbands, band, dk
write(ilog,"(A,3F10.6)") "k-point in reciprocal space:", (kp(j), j=1,3)
write(ilog,*)

E(0,0,0) = get_next_eigeval(iunt, band, nbands)

! http://stackoverflow.com/questions/9791001/loop-in-fortran-from-a-list
A = (/ -2, -1, 1, 2 /)
! x
do i = 1, size(A)
    i1 = A(i)
    E(i1,0,0) = get_next_eigeval(iunt, band, nbands)
end do

! y
do i = 1, size(A)
    i1 = A(i)
    E(0,i1,0) = get_next_eigeval(iunt, band, nbands)
end do

! z
do i = 1, size(A)
    i1 = A(i)
    E(0,0,i1) = get_next_eigeval(iunt, band, nbands)
end do

! xy
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        E(j1,i1,0) = get_next_eigeval(iunt, band, nbands)
    end do
end do

! xz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        E(j1,0,i1) = get_next_eigeval(iunt, band, nbands)
    end do
end do

! yz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        E(0,j1,i1) = get_next_eigeval(iunt, band, nbands)
    end do
end do

close(iunt)

E = E*ev2h ! to Hartree, we want to deal with a.u.

!d2N
m(1,1) = get_2nd_deriv(E, dk, 1, 0, 0)
m(2,2) = get_2nd_deriv(E, dk, 0, 1, 0)
m(3,3) = get_2nd_deriv(E, dk, 0, 0, 1)

!dxdN
m(2,1) = get_mixed_deriv1(E, dk, 1, 0)
m(3,1) = get_mixed_deriv1(E, dk, 0, 1)

!dydz
m(3,2) = (-63.0D0*(E(0,1,-2) + E(0,2,-1) + E(0,-2,1) +  E(0,-1,2))&
         &+63.0D0*(E(0,-1,-2)+ E(0,-2,-1)+ E(0,1,2) +   E(0,2,1))&
         &+44.0D0*( E(0,2,-2)+ E(0, -2, 2)-  E(0,-2,-2)- E(0,2,2))&
         &+74.0D0*( E(0,-1,-1)+E(0, 1, 1) -  E(0,1,-1)-E(0,-1,1)))&
         &/(600.0D0*dk**2)

m(1,2) = m(2,1)
m(1,3) = m(3,1)
m(2,3) = m(3,2)

write(ilog,*) "Original matrix:"
write(ilog,"(3F15.8)") ((m(i,j), j=1,3),i=1,3)
write(ilog,*)

! At this point m can be either (a) inverted and diagonalized or (b) diagonalized and eigenvalues inverted.
! The problem with (a) is that the inverse procedure does NOT guarantee symmetric matrix as the result.
! (a):
!   call inverse(m, 3)
!   call DGEEV('N', 'V', 3, m, 3, b, bi, DUMMY, 1, ev, 3, WORK, 12, ok)
!
! (b):
call DSYEV( 'V', 'U', 3, m, 3, b, WORK, size(WORK), ok )
if (ok .eq. 0) then
    write(ilog,*) "Principle effective masses and directions:"
    do i=1,3
        write(ilog,"(A25,F10.3)") "Effective mass:", 1.0D0/b(i)
        write(ilog,"(A25,3F10.6)") "Cartesian coordinate:", (m(j,i), j=1,3)
        v = pureDGESV(f, m(:,i), 'T')
        v = normal(v, size(v))
        write(ilog,"(A25,3F10.3)") "Direct lattice vectors:", (v(j), j=1,3)
        write(ilog,*)
    end do
else
    write (ilog,*) "INFO from DGEEV: ", ok
endif

end program

real(kind=8) function get_2nd_deriv(E, dk, x, y, z)
    implicit none
    integer(kind=4) :: x, y, z
    real(kind=8) :: E(-2:2,-2:2,-2:2), dk

    get_2nd_deriv = (-E(-2*x,-2*y,-2*z) + 16.0D0*E(-1*x,-1*y,-1*z) - 30.0D0*E(0,0,0)&
         &+16.0D0*E(1*x,1*y,1*z) - E(2*x,2*y,2*z))/(12.0D0*dk**2.0D0)
    return
end function

real(kind=8) function get_mixed_deriv1(E, dk, y, z)
    implicit none
    integer(kind=4) :: x, y, z
    real(kind=8) :: E(-2:2,-2:2,-2:2), dk

    get_mixed_deriv1=&
    &(-63.0D0*(E(1,-2*y,-2*z) + E(2,-1*y,-1*z) + E(-2,1*y,1*z) +  E(-1,2*y,2*z))&
     &+63.0D0*(E(-1,-2*y,-2*z)+ E(-2,-1*y,-1*z)+  E(1,2*y,2*z) +   E(2,1*y,1*z))&
     &+44.0D0*( E(2,-2*y,-2*z)+ E(-2, 2*y, 2*z)-  E(-2,-2*y,-2*z)- E(2,2*y,2*z))&
    &+74.0D0*( E(-1,-1*y,-1*z)+E(1, 1*y, 1*z) -  E( 1,-1*y,-1*z)-E(-1,1*y,1*z)))&
    &/(600.0D0*dk**2.0D0)
    return
end function

real(kind=8) function get_next_eigeval(iunt, band, nbands)
    implicit none
    integer(kind=4) :: iunt, band, nbands, j, itrash

    read(iunt,fmt=*)
    read(iunt,fmt=*)
    do j=1,nbands
        if(j == band) then
            read(iunt,fmt=*) itrash, get_next_eigeval
            ! write(*,*) get_next_eigeval
        else
            read(iunt,fmt=*)
        end if
    end do
    return
end function
