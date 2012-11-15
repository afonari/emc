! http://www.tek-tips.com/viewthread.cfm?qid=1516435
module emc_functions
    implicit none
contains
function cart2fract(f, cart_coords)
    implicit none
    integer(kind=4), parameter :: size = 3
    real(kind=8), intent(in) :: f(size,size), cart_coords(size)
    real(kind=8) :: cart2fract(size)
    integer(kind=4) :: IPIV(size), INFO, i, j
    external :: DGETRS

    cart2fract = cart_coords
    call DGETRS( 'N', size, 1, f, size, IPIV, cart2fract, size, INFO )
    if (INFO /= 0) then
        write(*,*) "INFO from cart2fract function: ", INFO
        cart2fract = 0.d0
    else
        cart2fract = cart_coords
    endif
    return
end function cart2fract

function fract2cart(f, frac_coords)
    implicit none
    integer(kind=4), parameter :: size = 3
    real(kind=8), intent(in) :: f(size,size), frac_coords(size)
    real(kind=8) :: fract2cart(size), ft(size,size)
    external :: DGEMV

    call DGEMV('t', size, size, 1.d0, f, size, frac_coords, 1, 1.d0, fract2cart, 1)
    return
end function fract2cart

subroutine normalize(v, n)
    implicit none
    integer(kind=4) :: i, n
    real(kind=8) :: v(n), norm

    norm = 0.d0
    do i=1,n
        norm = norm + v(i)**2
    end do

    v = v/dsqrt(norm)
    return
end subroutine

subroutine normal(a,n)
    implicit none
    integer(kind=4) :: n, i
    real(kind=8) :: a(n), amax

    amax=0.0d0
    do i=1, n
       if (dabs(amax)<dabs(a(i))) amax=a(i)
    end do
    do i=1, n
       a(i)=a(i)/amax
    end do
    return
end

subroutine print_time(iunt)
    implicit none
    integer(kind=4), intent(in) :: iunt
    integer(kind=4) :: today(3), now(3)

    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
    write(iunt,"('Date ',i2.2,'/',i2.2,'/',i4.4,';time ',i2.2,':',i2.2,':',i2.2 )")today(2),today(1),today(3),now
    return
end subroutine print_time

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

end module emc_functions
