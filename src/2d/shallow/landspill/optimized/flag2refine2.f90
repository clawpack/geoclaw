! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for GeoClaw for tsunami applications and related problems
!
!
! The logical function allowflag(x,y,t) is called to
! check whether further refinement at this level is allowed in this cell
! at this time.
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)
!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags)

    use amr_module, only: mxnest, t0, DOFLAG, UNSET
    use topo_module, only: tlowtopo,thitopo,xlowtopo,xhitopo,ylowtopo,yhitopo
    use topo_module, only: minleveltopo,mtopofiles
    use landspill_module, only: point_sources

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp

    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Flagging
    real(kind=8), intent(in out) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)

    logical :: allowflag
    external allowflag

    ! Generic locals
    integer :: i,j,m
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi

    ! Don't initialize flags, since they were already
    ! flagged by flagregions2
    ! amrflags = DONTFLAG

    ! if t == t0, flag the cells containing point sources
    if (t .eq. t0) then ! should we use abs(...) .lt. some tolerance?
        call point_sources%flag_cells( &
            mbuff, mx, my, xlower, ylower, dx, dy, amrflags)
    endif

    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy

        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

            ! Check to see if refinement is forced in any topography file region:
            do m=1,mtopofiles
                if (level < minleveltopo(m) .and. t >= tlowtopo(m) .and. t <= thitopo(m)) then
                    if (  x_hi > xlowtopo(m) .and. x_low < xhitopo(m) .and. &
                          y_hi > ylowtopo(m) .and. y_low < yhitopo(m) &
                           .and. amrflags(i,j) == UNSET) then

                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                endif
            enddo

            ! -----------------------------------------------------------------
            ! Refinement not forced, so check if it is allowed 
            ! and if the flag is still UNSET. If so,
            ! check if there is a reason to flag this point.
            ! If flag == DONTFLAG then refinement is forbidden by a region,
            ! if flag == DOFLAG checking is not needed
            if (allowflag(x_c,y_c,t,level) .and. amrflags(i,j) == UNSET) then
                if (q(1,i,j) > 0D0) then
                    amrflags(i,j) = DOFLAG
                endif
            endif

        enddo x_loop
    enddo y_loop
end subroutine flag2refine2
