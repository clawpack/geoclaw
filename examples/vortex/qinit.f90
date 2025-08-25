
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: g => grav
    use amr_module, only: t0
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i, j
    real(kind=8) :: x, y, u, v

    ! Constants
    real(kind=8), parameter :: pi = 4.d0 * atan (1.d0)
    real(kind=8), parameter :: c(2) = [0.04d0, 0.02d0]
    real(kind=8), parameter :: x0(2) = [-20.d0, -10.d0]
    real(kind=8), parameter :: M = 0.5d0
    real(kind=8), parameter :: alpha = pi / 6.d0

    do i=1-mbc, mx+mbc
        x = xlower + (i - 0.5d0) * dx
        do j=1-mbc, my+mbc
            y = ylower + (j - 0.5d0) * dy

            q(1, i, j) = 1.d0 - c(1)**2 / (4.d0 * c(2) * g) * exp(2.d0 * f(x, y, t0))
            u = M * cos(alpha) + c(1) * (y - x0(2) - M * t0 * sin(alpha)) * exp(f(x, y, t0))
            v = M * sin(alpha) - c(1) * (x - x0(1) - M * t0 * cos(alpha)) * exp(f(x, y, t0))
            q(2, i, j) = u * q(1, i, j)
            q(3, i, j) = v * q(1, i, j)
        end do
    end do

contains

    real(kind=8) pure function f(x, y, t)
        implicit none
        real(kind=8), intent(in) :: x, y, t

        f = -c(2) * (  (x - x0(1) - M * t * cos(alpha))**2           &
                     + (y - x0(2) - M * t * sin(alpha))**2)
    end function f

end subroutine qinit
