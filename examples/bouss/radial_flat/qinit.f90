
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j
    real(kind=8) :: r0,x,y,r

    real(kind=8) :: eta, width, ampl

    width = 100.d0    ! controls width of Gaussian
    r0 = 0.d0   ! initial radius of Gaussian
    ampl = 5.0d0  ! amplitude

    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall


    do i=1,mx
      x = (xlower + (i-0.5d0)*dx)
      do j=1,my
          if (q(1,i,j) > 0.d0) then
              y = (ylower + (j-0.5d0)*dy)
              r = sqrt(x**2 + y**2)
              eta = ampl * exp(-(((r-r0)/width)**2))
              
              q(1,i,j) = max(0.d0, eta - aux(1,i,j))
              q(2,i,j) = 0.d0
              q(3,i,j) = 0.d0

              endif
          enddo
       enddo

    
end subroutine qinit
