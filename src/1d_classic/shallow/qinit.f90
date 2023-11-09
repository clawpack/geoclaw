subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version simply sets eta = max(h + B, sea_level)

    ! For more specific initial conditions
    !  copy this to an application directory and
    !  loop over all grid cells to set values of q(1:meqn, 1:mx).

    use geoclaw_module, only: sea_level

    ! uncomment if any of these needed...
    !use geoclaw_module, only: dry_tolerance, grav
    !use grid_module, only: xcell,xp_edge,mx_edge
    !use topo_module, only: zcell0 ! initial topo in each cell

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: eta

    do i=1,mx
      eta = sea_level
      q(1,i) = max(0.d0, eta - aux(1,i))
      q(2,i) = 0.d0
   enddo


end subroutine qinit
