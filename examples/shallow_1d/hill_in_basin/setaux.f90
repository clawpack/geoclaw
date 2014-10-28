subroutine setaux(mbc,mx,xlower,dx,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This default version sets aux(1,:) to b(x) for shallow flow.
    ! this version simply sets b(x) = 0.0 ...for specific problems
    !   copy to an application directory

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    !use geoclaw_module, only: grav  !uncomment if needed

    implicit none
    integer, intent(in) :: mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xcell


    do i=1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         !call mapc2p(xcell,xp)
         !xp = xcell

         aux(1,i)=1.8d0*exp(-10.d-2*xcell**2) + .5d0*1.0d0*exp(-.5d-2*(xcell-0.5d2)**2) &
            &  + 1.8d0*exp(-10.d-2*(xcell-1.d2)**2)

    enddo

end subroutine setaux
