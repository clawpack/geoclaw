subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use grid_module, only: xcell

    ! uncomment if any of these needed...
    use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: x0,width,eta

    x0 = -10.d3   ! initial location of step

    do i=1,mx
      if (xcell(i) > x0) then 
          eta = 0.d0
        else
          eta = 1.d0
        endif
      q(1,i) = max(0.0, eta - aux(1,i))
      q(2,i) = eta*sqrt(grav*q(1,i))  ! right-going

   enddo


end subroutine qinit
