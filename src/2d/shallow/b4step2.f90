! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,actualstep)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use geoclaw_module, only: speed_limit
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET
    use amr_module, only: outunit

    use storm_module, only: set_storm_fields

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    logical, intent (in) :: actualstep

    ! Local storage
    integer :: index,i,j,k,dummy
    real(kind=8) :: h,u,v, sratio, xs,ys, s

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    end forall

    ! Check for fluid speed sqrt(u**2 + v**2) > speed_limit
    ! and reset by scaling (u,v) down to this value (preserving direction)
    ! Note: similar check is done in getmaxspeed
    ! This helps avoid too many dt reductions when flow off very steep topo
    ! with delta B larger than fluid depth gives big speeds in Riemann solution 
    ! (shallow water equations aren't valid for flow off a cliff)

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (q(1,i,j) > 0.d0) then
                s = sqrt((q(2,i,j)**2 + q(3,i,j)**2)) / q(1,i,j)
                if (s > speed_limit) then
                    sratio = speed_limit / s
                    q(2,i,j) = q(2,i,j) * sratio
                    q(3,i,j) = q(3,i,j) * sratio

                    if (.true.) then
                        ! write out info useful for investigating topo:
                        xs = xlower + (i-0.5d0)*dx
                        ys = ylower + (j-0.5d0)*dy
                        write(6,604) t, i,j, mx,my, dx
                        write(outunit,604) t, i,j, mx,my, dx
                        write(6,603) s,q(1,i,j),aux(1,i,j),xs,ys
                        write(outunit,603) s,q(1,i,j),aux(1,i,j),xs,ys

 604                    format('b4step2 at t =',f10.2, '  i,j,mx,my:',4i4, &
                               '  dx = ',f10.7)
 603                    format('     reset s =',f9.2,'  h=',e11.3, ' B=',f8.2,&
                               '  x,y = ', f11.6,',',f10.6)
                    endif
               endif
            endif
        enddo
    enddo



    if (aux_finalized < 2 .and. actualstep) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    ! Set wind and pressure aux variables for this grid
    call set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

end subroutine b4step2
