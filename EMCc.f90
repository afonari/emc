!##########################################++##########################
!####                                                               
!#### Effective Mass Calculator
!#### c.alculator program
!#### licensed under MIT LICENSE
!####                                                               
!######################################################################

program EMCg ! version 1.0
implicit none
real(kind=8) :: eigenval, get_next_eigeval, get_2nd_deriv, get_mixed_deriv1, ft(3,3), ev(3,3)
real(kind=8) :: kp(3), dk, E(-2:2,-2:2,-2:2), A(4), m(3,3), b(3), bi(3), WORK(100), DUMMY(1,1)
real(kind=8) :: trash1(3), trash2(3), tmp(3)

real(kind=8), parameter :: a2b = 1.0D0/0.52917721092D0
real(kind=8), parameter :: b2a = 0.52917721092D0
real(kind=8), parameter :: pi = 3.14159265358979324D0
real(kind=8), parameter :: ev2h = 1.0D0/27.21138505D0

integer h,i,i1,j,j1,k,l, LWORK, ok
integer(kind=4) :: iunt, nkpoints, itrash, count, nbands, band
character*3 version_number
character*20 wc1
character(len=80) :: cha
character(len=1) :: prg
external DSYEV

version_number='1.0'
nkpoints = 61

write(*,*) "Effective Mass Calculator calculator ", version_number
write(*,*)

! read input ########################################################

iunt=10
open(unit=iunt,file='inp',form='formatted')
    read(iunt,fmt=*) cha
    read(iunt,fmt=*) dk
    read(iunt,fmt=*) band    
    read(iunt,fmt=*) prg
close(iunt)

!if(prg .eq. 'V') then
!    write(*,*) "dk will be converted to VASP units (dk*2*Pi/a2b)"
!    dk = dk/(2.0D0*pi*b2a)
!end if
write(*,*) "band, dk: ", band, dk

! read OUTCAR ###########################################

!open(unit=iunt,file='OUTCAR',form='formatted')
!    cha="direct lattice vectors"
!    call search_key(cha, iunt)

    !reading TRANSPOSE f
!    do i=1,3
!        read(iunt, fmt=*)ft(i,1),ft(i,2),ft(i,3)
!    end do
!close(iunt)

! call inverse(ft, 3)
!write(*,*) "Inverse transposed f matrix:"
!do i=1,3
!    write(*,*) ft(1,i), ft(2,i), ft(3,i)
!enddo

! read EIGENVAL ###########################################

open(unit=iunt,file='EIGENVAL',form='formatted')

read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*)
read(iunt,fmt=*) itrash, itrash, nbands

write(*,*) "NBands: ", nbands

! x
do i=-2,2
    E(i,0,0) = get_next_eigeval(iunt, band, nbands)
end do

! http://stackoverflow.com/questions/9791001/loop-in-fortran-from-a-list
A = (/ -2, -1, 1, 2 /)

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

write(*,*) "Original matrix:"
write(*,"(3F15.8)") ((m(i,j), j=1,3),i=1,3)
write(*,*)

! At this point m can be either (a) inverted and diagonalized or (b) diagonalized and eigenvalues inverted.
! The problem with (a) is that the inverse procedure does NOT guarantee symmetric matrix as the result.
! (a):
!   call inverse(m, 3)
!   call DGEEV('N', 'V', 3, m, 3, b, bi, DUMMY, 1, ev, 3, WORK, 12, ok)
!
! (b):
call DSYEV( 'V', 'U', 3, m, 3, b, WORK, 100, ok )
if (ok .eq. 0) then
    write(*,*) "Principle effective masses and directions:"
    do i=1, 3
        write(*,"(A25,F10.2)") "Effective mass:", 1.0D0/b(i)
        write(*,"(A25,3F10.6)") "Cartesian coordinate:", (m(j,i), j=1,3)
    !       call cart2fract(vabc,emabc(1,i),em(1,i))
    !       call normal(emabc(1,i),3)
    !       write(*,"(A25,3F10.3)") "Direct lattice vectors:", (emabc(j,i), j=1,3)
    end do
else
    write (*,*) "An error occured while diagonalizing matrix in DGEEV()"
endif

!write(*,*) "Directions:"
!do i=1,3
!    ! trash1 = (/ ev(1,i), ev(2,i), ev(3,i) /)
!    call DGEMV('N', 3, 3, 1.0d0, ft, 3, ev(:,i), 1, 0.0d0, trash2, 1)
!    call normalize(trash2, 3)
!    write(*,"(F10.5,F10.5,F10.5)") trash2(1), trash2(2), trash2(3)
!end do

end program

subroutine set_next_eigeval(iunt, kp, i, j, k, w, dk)
    implicit none
    integer(kind=4) :: iunt, i, j, k
    real(kind=8) :: kp(3), w, dk

    write(unit=iunt,fmt='(F15.10,F15.10,F15.10,F5.1)') kp(1)+i*dk, kp(2)+j*dk, kp(3)+k*dk, w
    return
end subroutine

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

! http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
subroutine inverse(M, n)
    implicit none
    integer(kind=4) :: n, LWORK, ok
    real(kind=8) :: M(n,n)
    integer(kind=4),allocatable :: ipiv(:)
    real(kind=8),allocatable :: WORK(:)
    external DGETRF, DGETRI

    LWORK = n*n
    allocate(WORK(LWORK))
    allocate(ipiv(n+1))

    call DGETRF(n, n, M, n, ipiv, ok)
    if (ok > 0) then
        write(*,*) 'Matrix is singular in: ', ok
        return
    end if
    
    call DGETRI(n, M, n, ipiv, WORK, LWORK, ok)

    if (ok > 0) then
        write(*,*) 'Matrix is singular in: ', ok
        return
    end if

    deallocate(ipiv)
    deallocate(WORK)
    return

end subroutine

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

subroutine normalize(v, n)
    implicit none
    real(kind=8) :: v(n), norm
    integer(kind=4) :: i, n

    norm = 0.0D0
    do i=1,n
        norm = norm + v(i)**2
    end do

    v = v/SQRT(norm)
    return
end subroutine

subroutine search_key(key_word,iunt)
    character(len=80) :: key_word
    integer(kind=4) :: iunt,n
    logical(kind=4) :: yesno
    character(len=30) :: cha
    n=len_trim(key_word)
    yesno=.true.
    do while(yesno)
       read(iunt,100) cha
       cha=trim(adjustl(cha))
       if(cha(1:n)==key_word) yesno=.false.
    end do
100 format(A30)
    return
end