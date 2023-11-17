subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use grid_module, only: xcell

    ! uncomment if any of these needed...
    !use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i

    real(kind=8) :: eta, width, x0, ampl, wavelength, pi, freq

    width = 100.    ! controls width of Gaussian
    x0 = 0.d0   ! initial location of Gaussian
    ampl = 5.0d0  ! amplitude
    !wavelength = 10d3
    !freq = 0.d0  !1.d0/wavelength
    pi = acos(-1d0)

    do i=1,mx
      !eta = ampl * exp(-((xcell(i)-x0)/width)**2) * cos(2*pi*xcell(i)*freq)
      eta = ampl * exp(-((xcell(i)-x0)/width)**2)
      q(1,i) = max(sea_level, eta - aux(1,i))
      q(2,i) = 0.d0  !eta*sqrt(grav*q(1,i))  ! right-going

   enddo


end subroutine qinit
