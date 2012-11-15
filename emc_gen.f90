! Copyright (c) "2012, by Georgia Institute of Technology
!                Contributors: Alexandr Fonari
!                Affiliation: Dr. Bredas group
!                URL: https://github.com/alexandr-fonari/emc
!                License: MIT License"

program emc_gen ! version 1.0
use emc_functions
implicit none
real(kind=8), parameter :: b2a = 0.52917721092d0
real(kind=8), parameter :: pi = 4.d0*DATAN(1.d0)
integer(kind=4), parameter :: nkpoints = 61
integer(kind=4), parameter :: iunt = 10, ilog = 11, w = 1 ! k point weight for VASP
character(len=3), parameter :: version_number = '1.0'

real(kind=8) :: kp(3), kpr(3), dk, E(-2:2,-2:2,-2:2), A(4), f(3,3), g(3,3)
integer(kind=4) :: i, i1, j, j1, band
character prg

open(unit=ilog,file='emc_gen.log',form='formatted')
write(ilog,*) "Effective Mass Calculator generator ", version_number
call print_time(ilog)
write(ilog,*)

open(unit=iunt,file='inp',form='formatted')
    read(iunt,fmt=*) (kp(i),i=1,size(kp))
    read(iunt,fmt=*) dk
    read(iunt,fmt=*) band
    read(iunt,fmt=*) prg
    read(iunt,fmt=*) ((f(i,j),j=1,3),i=1,3)
close(iunt)

if(prg .eq. 'V') then
    write(ilog,*) "dk will be converted to VASP units (2Pi/A)"
    dk = dk/(2.0D0*pi*b2a)
end if

! g = inverse(transpose(f)), will do it in two steps
g=f
g=transpose(g)
call inverse(g, 3)

write(ilog,*) "direct lattice vectors              reciprocal lattice vectors"
do i=1,3
    write(ilog,"(3F10.6,A,3F10.6)") , (f(i,j), j=1,3), "     ", (g(i,j), j=1,3)
end do

kpr = fract2cart(g, kp)
write(ilog,*) "k-point in: reciprocal space       reciprocal Cartesian space"
write(ilog,"(3F10.6,A,3F10.6)") , (kp(j), j=1,3), "     ", (kpr(j), j=1,3)


! write KPOINTS file ###########################################

open(unit=iunt,file='KPOINTS',form='formatted')

write(unit=iunt,fmt='(A,F7.4,F7.4,F7.4,A,F7.5)') 'kp: ', (kpr(j), j=1,3), ', dk: ', dk
write(unit=iunt,fmt='(I2)') nkpoints
write(unit=iunt,fmt='(A)') 'Cartesian'

! x
do i=-2,2
    call set_next_eigeval(iunt, kpr, i, 0, 0, w, dk)
end do

! http://stackoverflow.com/questions/9791001/loop-in-fortran-from-a-list
A = (/ -2, -1, 1, 2 /)

! y
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, kpr, 0, i1, 0, w, dk)
end do

! z
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, kpr, 0, 0, i1, w, dk)
end do

! xy
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kpr, j1, i1, 0, w, dk)
    end do
end do

! xz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kpr, j1, 0, i1, w, dk)
    end do
end do

! yz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kpr, 0, j1, i1, w, dk)
    end do
end do

close(30)

end program

subroutine set_next_eigeval(iunt, kp, i, j, k, w, dk)
    implicit none
    integer(kind=4), intent(in) :: iunt, i, j, k, w
    real(kind=8), intent(in) :: kp(3), dk

    write(unit=iunt,fmt='(F15.10,F15.10,F15.10,F5.1)') kp(1)+i*dk, kp(2)+j*dk, kp(3)+k*dk, dble(w)
    ! write(unit=iunt,fmt='(I5,I5,I5)') i, j, k
    return
end subroutine