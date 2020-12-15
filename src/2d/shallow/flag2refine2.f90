! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for GeoClaw for tsunami applications and related problems
!
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
! This routine only checks cells in which amrflags(i,j) = UNSET,
! If the value is already DONTFLAG or DOFLAG then this we determined by
! flagregions or other criteria and is not modified by this routine.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags)

    use amr_module, only: mxnest, t0, DOFLAG, UNSET
    use geoclaw_module, only:dry_tolerance, sea_level
    use geoclaw_module, only: spherical_distance, coordinate_system


    use storm_module, only: storm_specification_type, wind_refine, R_refine
    use storm_module, only: storm_location, wind_forcing, wind_index, wind_refine

    use regions_module, only: num_regions, regions
    use refinement_module, only: wave_tolerance, speed_tolerance

    use adjoint_module, only: totnum_adjoints,innerprod_index, &
                              adjoint_flagging,select_snapshots
    use adjointsup_module, only: calculate_innerproduct

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp

    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Flagging
    real(kind=8), intent(in out) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)

    ! Generic locals
    integer :: i,j,m,r
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: speed, eta, ds

    ! Storm specific variables
    real(kind=8) :: R_eye(2), wind_speed

    ! Adjoint method specific variables
    logical mask_selecta(totnum_adjoints)


    if(adjoint_flagging) then
        aux(innerprod_index,:,:) = 0.0
        call select_snapshots(t,mask_selecta)

        ! Loop over adjoint snapshots
        aloop: do r=1,totnum_adjoints

            ! Consider only snapshots that are within the desired time range
            if (mask_selecta(r)) then
                ! Calculate inner product with current snapshot
                call calculate_innerproduct(q,r,mx,my,xlower,   &
                        ylower,dx,dy,meqn,mbc,maux,aux)
            endif

        enddo aloop
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
            
            if (amrflags(i,j) .ne. UNSET) then
                ! do not consider this cell further
                cycle x_loop
            endif

            ! ************* Storm Based Refinement ****************
            ! Check to see if we are some specified distance from the eye of
            ! the storm and refine if we are
            if (storm_specification_type /= 0) then
                R_eye = storm_location(t)
                do m=1,size(R_refine,1)
                    if (coordinate_system == 2) then
                        ds = spherical_distance(x_c, y_c, R_eye(1), R_eye(2))
                    else
                        ds = sqrt((x_c - R_eye(1))**2 + (y_c - R_eye(2))**2)
                    end if
                    
                    if (ds < R_refine(m) .and. level <= m) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                enddo
                
                ! Refine based on wind speed
                if (wind_forcing) then
                    wind_speed = sqrt(aux(wind_index,i,j)**2 + aux(wind_index+1,i,j)**2)
                    do m=1,size(wind_refine,1)
                        if ((wind_speed > wind_refine(m)) .and. (level <= m)) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif
            endif

            ! ********* criteria applied only to wet cells:

            if (q(1,i,j) > dry_tolerance) then

                if(adjoint_flagging) then
                    if(aux(innerprod_index,i,j) > tolsp) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                else
                
                    ! Check wave criteria
                    eta = q(1,i,j) + aux(1,i,j)
                    if (abs(eta - sea_level) > wave_tolerance) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif

                    ! Check speed criteria, note that it might be useful to
                    ! also have a per layer criteria since this is not
                    ! gradient based
                    speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                    do m=1,min(size(speed_tolerance),mxnest)
                        if (speed > speed_tolerance(m) .and. level <= m) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif
            endif


        enddo x_loop
    enddo y_loop
end subroutine flag2refine2
