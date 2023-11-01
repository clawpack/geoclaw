subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use grid_module, only: xcell

    use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i

    real(kind=8) :: eta, width, x0, ampl

    width = 2.d3    ! controls width of Gaussian
    x0 = -85.d3   ! initial location of Gaussian
    ampl = 2.d0  ! amplitude

    do i=1,mx
      eta = ampl * exp(-((xcell(i)-x0)/width)**2)
      q(1,i) = max(sea_level, eta - aux(1,i))
      q(2,i) = eta*sqrt(grav*q(1,i))  ! right-going

   enddo


end subroutine qinit
