
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level
    use amr_module, only: t0
    use qinit_module, only: qinit_type,add_perturbation
    use qinit_module, only: variable_eta_init
    use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    use qinit_module, only: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry
    use qinit_module, only: dx_fdry, dy_fdry
    use qinit_module, only: tend_force_dry
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m, ii,jj
    real(kind=8) :: x,y
    real(kind=8) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: ddxy
    
    q = 0.d0   ! initialize all elements to 0
    
    if (variable_eta_init) then
        ! Set initial surface eta based on eta_init
        call set_eta_init(mbc,mx,my,xlower,ylower,dx,dy,t0,veta)
      else
        veta = sea_level  ! same value everywhere
      endif

    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, veta(i,j) - aux(1,i,j))
    end forall

    if (use_force_dry .and. (t0 <= tend_force_dry)) then
     ! only use the force_dry if it specified on a grid that matches the 
     ! resolution of this patch, since we only check the cell center:
     ddxy = max(abs(dx-dx_fdry), abs(dy-dy_fdry))
     if (ddxy < 0.01d0*min(dx_fdry,dy_fdry)) then
       do i=1,mx
          x = xlower + (i-0.5d0)*dx
          ii = int((x - xlow_fdry + 1d-7) / dx_fdry)
          do j=1,my
              y = ylower + (j-0.5d0)*dy
              jj = int((y - ylow_fdry + 1d-7) / dy_fdry)
              jj = my_fdry - jj  ! since index 1 corresponds to north edge
              if ((ii>=1) .and. (ii<=mx_fdry) .and. &
                  (jj>=1) .and. (jj<=my_fdry)) then
                  ! grid cell lies in region covered by force_dry,
                  ! check if this cell is forced to be dry 
                  ! Otherwise don't change value set above:                  
                  if (force_dry(ii,jj) == 1) then
                      q(1,i,j) = 0.d0
                      endif
                  endif
          enddo ! loop on j
       enddo ! loop on i
       endif ! dx and dy agree with dx_fdry, dy_fdry
    endif ! use_force_dry

    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif
    
end subroutine qinit
