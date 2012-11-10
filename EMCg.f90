!##########################################++##########################
!####                                                               
!#### Effective Mass Calculator
!#### g.enerator program
!#### licensed under MIT LICENSE
!####                                                               
!######################################################################

program EMCg ! version 1.0
implicit none

real(kind=8) :: kp(3), dk, E(-2:2,-2:2,-2:2), count, A(4), w
integer(kind=4) :: h,i,i1,j,j1,k,l, iunt, nkpoints, band
character prg
character*3 version_number
character*20 wc1

real(kind=8), parameter :: a2b = 1.0D0/0.52917721092D0
real(kind=8), parameter :: b2a = 0.52917721092D0
real(kind=8), parameter :: pi = 3.14159265358979324D0

version_number='1.0'

write(*,*) "Effective Mass Calculator generator ", version_number
write(*,*)

count = 1
iunt = 10
nkpoints = 61
w = 1.0 ! k point weight for VASP

! read input ########################################################

open(unit=iunt,file='inp',form='formatted')
    read(iunt,fmt=*) (kp(i),i=1,size(kp))
    read(iunt,fmt=*) dk
    read(iunt,fmt=*) band
    read(iunt,fmt=*) prg
close(iunt)

if(prg .eq. 'V') then
    write(*,*) "dk will be converted to VASP units (2Pi/A)"
    dk = dk/(2.0D0*pi*b2a)
end if
  write(*,*) "k-point: ", (kp(i),i=1,size(kp)), "dk: ", dk


! write KPOINTS file ###########################################

open(unit=iunt,file='KPOINTS',form='formatted')

write(unit=iunt,fmt='(A,F7.3,F7.3,F7.3,A,F8.5)') 'central point: ', kp(1), kp(2), kp(3), ', dk: ', dk
write(unit=iunt,fmt='(I2)') nkpoints
write(unit=iunt,fmt='(A)') 'Cartesian'

! x
do i=-2,2
    call set_next_eigeval(iunt, kp, i, 0, 0, w, dk)
end do

! http://stackoverflow.com/questions/9791001/loop-in-fortran-from-a-list
A = (/ -2, -1, 1, 2 /)

! y
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, kp, 0, i1, 0, w, dk)
end do

! z
do i = 1, size(A)
    i1 = A(i)
    call set_next_eigeval(iunt, kp, 0, 0, i1, w, dk)
end do

! xy
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kp, j1, i1, 0, w, dk)
    end do
end do

! xz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kp, j1, 0, i1, w, dk)
    end do
end do

! yz
do i = 1, size(A)
    i1 = A(i)
    do j=1, size(A)
        j1 = A(j)
        call set_next_eigeval(iunt, kp, 0, j1, i1, w, dk)
    end do
end do

close(30)

end program

subroutine set_next_eigeval(iunt, kp, i, j, k, w, dk)
    implicit none
    integer(kind=4) :: iunt, i, j, k
    real(kind=8) :: kp(3), w, dk

    write(unit=iunt,fmt='(F15.10,F15.10,F15.10,F5.1)') kp(1)+i*dk, kp(2)+j*dk, kp(3)+k*dk, w
    ! write(unit=iunt,fmt='(I5,I5,I5)') i, j, k
    return
end subroutine