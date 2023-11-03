
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level, grav
    use geoclaw_module, only: coordinate_system, earth_radius
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m
    real(kind=8) :: x0,y0,x,y,dx0,dy0,dsigma,r,hU
    real(kind=8) :: cos_theta,sin_theta,a1,a2,u,theta,xx,yy
    integer :: k

    real(kind=8) :: eta, width, ampl

    width = 10.d3    ! controls width of Gaussian
    ampl = 20.0d0  ! amplitude

    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall

    x0 = 0.d0
    y0 = 0.d0

    do i=1,mx
      x = (xlower + (i-0.5d0)*dx)
      dx0 = x - x0
      do j=1,my
          if (q(1,i,j) > 0.d0) then
              y = (ylower + (j-0.5d0)*dy)
              dy0 = y - y0
              r = sqrt(dx0**2 + dy0**2)
              eta = ampl * exp(-((r/width)**2))
              
              q(1,i,j) = max(0.d0, eta - aux(1,i,j))
              q(2,i,j) = 0.d0
              q(3,i,j) = 0.d0

              endif
          enddo
       enddo

    
end subroutine qinit
