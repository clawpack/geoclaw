! ::::::::::::::::::::: flagregions ::::::::::::::::::::::::::::::::::
!
!> Modify array of flagged points to respect minlevels and maxlevels
!! specified by regions. This subroutine processes ONE grid associated with this amrflags
!!
!! Second version with outer loop on regions, should be faster.
!!
!! This routine may change flags only in cells that are (partially)
!! covered by one or more regions. amrflags will be later modified 
!! by Richardson extrapolation and/or flag2refine routine, as requested,
!! which will only add DOFLAG points to cells that are still UNSET
!!
!! If any part of a grid cell is covered by one or more regions, then
!! refinement is *required* to at least the max of all region min_levels
!! and is *allowed* to at most to the max of all region max_levels.
!!
!! Note that buffering is done *after* this, so additional cells may
!! be refined in areas covered by regions that do not allow it!

!! amrflags  = array to be flagged with either the value
!!             DONTFLAG (no refinement needed)  or
!!             DOFLAG   (refinement desired)    
!!
!! \param mx number of cells in *i* direction
!! \param my number of cells in *j* direction
!! \param mbuff width of buffer region
!! \param maux number of auxiliary variables
!! \param xlower x-coordinate of left physical boundary
!! \param ylower y-coordinate of lower physical boundary
!! \param dx spacing in *i* direction
!! \param dy spacing in *j* direction
!! \param level AMR level of this grid
!! \param t simulation time on this grid
!! \param amrflags array to be flagged with either the value **DONTFLAG** or **DOFLAG** for each cell. 
!!        It is enlarged from grid size to include buffer regions around the grid. 
!! \param DONTFLAG value to be assigned to amrflags for cells that need no refinement
!! \param DOFLAG value to be assigned to amrflags for cells that do need refinement
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine flagregions2(mx,my,mbuff,xlower,ylower,dx,dy,level,t, &
                            amrflags,mptr)

    use regions_module
    use amr_module, only : DOFLAG, UNSET, DONTFLAG

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,level,mbuff,mptr
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    
    type(ruled_region_type), pointer :: rr

    ! Locals
    integer :: i,j,m,i1,i2,j1,j2, k,k1,k2
    real(kind=8) :: x_low,y_low,x_hi,y_hi, xupper,yupper
    integer :: minlevel(mx,my), maxlevel(mx,my)
    !integer, allocatable :: minlevel(:,:), maxlevel(:,:)
    integer :: min_current_minlevel, min_current_maxlevel

    real(kind=8) :: x_lower_y_low, x_lower_y_hi, x_upper_y_low, x_upper_y_hi, &
                    x_lower_s, x_upper_s, alpha_y_low, alpha_y_hi, &
                    y_lower_x_low, y_lower_x_hi, y_upper_x_low, y_upper_x_hi, &
                    y_lower_s, y_upper_s, alpha_x_low, alpha_x_hi

#ifdef WHERE_AM_I
      write(*,*)"starting flagregions2"
#endif

    !write(58,*) 'level, lower: ',level,xlower,ylower  ! +++
    !allocate(minlevel(mx,my), maxlevel(mx,my))
    
    minlevel = 0
    maxlevel = 0

    xupper = xlower + mx*dx
    yupper = ylower + my*dy

    ! mark ghost cells with DONTFLAG, in case they were UNSET:
    do i=1-mbuff,mx+mbuff
        do j=1-mbuff,0
            amrflags(i,j) = DONTFLAG
            enddo
        do j=my+1,my+mbuff
            amrflags(i,j) = DONTFLAG
            enddo
        enddo

    do j=1-mbuff,my+mbuff
        do i=1-mbuff,0
            amrflags(i,j) = DONTFLAG
            enddo
        do i=mx+1,mx+mbuff
            amrflags(i,j) = DONTFLAG
            enddo
        enddo

    ! loop on old-style regions, still supported for now:
    
    rloop: do m=1,num_regions
        if (t < regions(m)%t_low .or. t > regions(m)%t_hi) then
            cycle rloop  ! no intersection
        endif
        
        if (xlower >= regions(m)%x_hi .or. xupper <= regions(m)%x_low) then
            cycle rloop  ! no intersection
        else
            i1 = max(floor((regions(m)%x_low - xlower) / dx) + 1, 1)
            i2 = min(floor((regions(m)%x_hi -xlower) / dx) + 1, mx)
        endif

        if (ylower >= regions(m)%y_hi .or. yupper <= regions(m)%y_low) then
            cycle rloop  ! no intersection
        else
            j1 = max(floor((regions(m)%y_low - ylower) / dy) + 1, 1)
            j2 = min(floor((regions(m)%y_hi - ylower) / dy) + 1, my)
        endif

        do j=j1,j2
            do i=i1,i2
                minlevel(i,j) = max(minlevel(i,j), regions(m)%min_level)
                maxlevel(i,j) = max(maxlevel(i,j), regions(m)%max_level)
            enddo
         enddo
    enddo rloop


    ! Loop on new "ruled rectangle" flagregions:
    
    rrloop: do m=1,num_rregions
    
        rr => rregions(m)

        if (t < rr%t_low .or. t > rr%t_hi) then
            cycle rrloop  ! no intersection
        endif
        
        ! compute intersection of patch with rr bounding box:
        if (xlower >= rr%x2bb .or. xupper <= rr%x1bb) then
            cycle rrloop  ! no intersection
        else
            i1 = max(floor((rr%x1bb - xlower) / dx) + 1, 1)
            i2 = min(floor((rr%x2bb -xlower) / dx) + 1, mx)
        endif

        if (ylower >= rr%y2bb .or. yupper <= rr%y1bb) then
            cycle rrloop  ! no intersection
        else
            j1 = max(floor((rr%y1bb - ylower) / dy) + 1, 1)
            j2 = min(floor((rr%y2bb - ylower) / dy) + 1, my)
        endif

        !write(58,*) 'ixy, ds: ',rr%ixy, rr%ds  ! +++
        ! patch overlaps bounding box, so need to check:
        
        ! first check if this rregion could affect anything already set:    
!        min_current_minlevel = minval(minlevel(i1:i2,j1:j2))
!        min_current_maxlevel = minval(maxlevel(i1:i2,j1:j2))
!        if ((min_current_minlevel >= rr%minlevel) .and. &
!            (min_current_maxlevel >= rr%maxlevel)) then
!                cycle rrloop ! this rregion won't change anything
!                endif

        !write(82,*) '+++ t = ',t, ' level =', level
        
        if (rr%ixy == 1) then
            ! rr%s corresponds to x, while lower,upper are in y
            ! outer loop on x, compute k1,k2 once for each col of cells
            iloop1: do i=i1,i2
                x_low = xlower + (i - 1) * dx
                x_hi = xlower + i * dx
                if (rr%ds > 0) then
                    k1 = floor((x_low-rr%s(1) + 1d-6) / rr%ds) + 1
                    k1 = max(k1, 1)
                    k2 = floor((x_hi-rr%s(1) + 1d-6) / rr%ds) + 1
                    k2 = min(k2, rr%nrules-1)
                else
                    ! rr%ds <= 0 means s point are not equally spaced
                    do k=1,rr%nrules
                        if (x_low < rr%s(k)) exit
                        enddo
                    k1 = max(k-1, 1)
                    do k=1,rr%nrules
                        if (x_hi < rr%s(k)) exit
                        enddo
                    k2 = min(k-1, rr%nrules-1)
                    endif ! rr%ds <= 0
                !write(58,*) '+++ ds = ',rr%ds,rr%nrules,k1,k2
                    
                    
                jloop1: do j=j1,j2
                    y_low = ylower + (j - 1) * dy
                    y_hi = ylower + j * dy
                    
                    if (rr%method == 0) then
                        ! pw constant interpolation (rectangles):
                        do k=k1,k2
                            if ((y_hi > rr%lower(k)) .and. &
                                (y_low < rr%upper(k))) then
                                minlevel(i,j) = max(minlevel(i,j), &
                                                    rr%min_level)
                                maxlevel(i,j) = max(maxlevel(i,j), &
                                                    rr%max_level)
                                cycle jloop1
                                endif
                            enddo
                    else 
                        ! rr%method == 1, linear interpolation (trapezoids):
                        do k=k1,k2
                            alpha_x_low = (x_low - rr%s(k)) &
                                        / (rr%s(k+1) - rr%s(k))
                            alpha_x_hi = (x_hi - rr%s(k)) &
                                        / (rr%s(k+1) - rr%s(k))
                            y_lower_x_low = (1-alpha_x_low) * rr%lower(k) &
                                         + alpha_x_low * rr%lower(k+1)
                            y_upper_x_low = (1-alpha_x_low) * rr%upper(k) &
                                          + alpha_x_low * rr%upper(k+1)
                            y_lower_x_hi = (1-alpha_x_hi) * rr%lower(k) &
                                         + alpha_x_hi * rr%lower(k+1)
                            y_upper_x_hi = (1-alpha_x_hi) * rr%upper(k) &
                                          + alpha_x_hi * rr%upper(k+1)
                            y_lower_s = min(y_lower_x_low, y_lower_x_hi)
                            y_upper_s = max(y_upper_x_low, y_upper_x_hi)
                            if ((y_hi > y_lower_s) .and. &
                                (y_low < y_upper_s)) then
                                  minlevel(i,j) = max(minlevel(i,j), &
                                                      rr%min_level)
                                  maxlevel(i,j) = max(maxlevel(i,j), &
                                                      rr%max_level)
                                  cycle jloop1
                                  endif
                             enddo  ! k                                  
                        endif ! rr%method
                    enddo jloop1
                enddo iloop1
                
            else ! rr%ixy = 2
                ! rr%s corresponds to y, while lower,upper are in x
                ! outer loop on y, compute k1,k2 once for each row of cells
                jloop2: do j=j1,j2
                    y_low = ylower + (j - 1) * dy
                    y_hi = ylower + j * dy
                    if (rr%ds > 0) then
                        k1 = floor((y_low-rr%s(1) + 1d-6) / rr%ds) + 1
                        k1 = max(k1, 1)
                        k2 = floor((y_hi-rr%s(1) + 1d-6) / rr%ds) + 1
                        k2 = min(k2, rr%nrules-1)
                    else
                        ! rr%ds <= 0 means s point are not equally spaced
                        do k=1,rr%nrules
                            if (y_low < rr%s(k)) exit
                            enddo
                        k1 = max(k-1, 1)
                        do k=1,rr%nrules
                            if (y_hi < rr%s(k)) exit
                            enddo
                        k2 = min(k-1, rr%nrules-1)
                        endif ! rr%ds <= 0
                    !write(58,*) '+++ ds = ',rr%ds,rr%nrules,k1,k2

                        
                    iloop2: do i=i1,i2
                        x_low = xlower + (i - 1) * dx
                        x_hi = xlower + i * dx
                        
                        if (rr%method == 0) then
                            ! pw constant interpolation (rectangles):
                            do k=k1,k2
                                if ((x_hi > rr%lower(k)) .and. &
                                    (x_low < rr%upper(k))) then
                                    minlevel(i,j) = max(minlevel(i,j), &
                                                        rr%min_level)
                                    maxlevel(i,j) = max(maxlevel(i,j), &
                                                        rr%max_level)
                                    cycle iloop2
                                    endif
                                enddo
                        else 
                            ! rr%method == 1, linear interpolation (trapezoids):
                            do k=k1,k2
                                alpha_y_low = (y_low - rr%s(k)) &
                                            / (rr%s(k+1) - rr%s(k))
                                alpha_y_hi = (y_hi - rr%s(k)) &
                                            / (rr%s(k+1) - rr%s(k))
                                x_lower_y_low = (1-alpha_y_low) * rr%lower(k) &
                                             + alpha_y_low * rr%lower(k+1)
                                x_upper_y_low = (1-alpha_y_low) * rr%upper(k) &
                                              + alpha_y_low * rr%upper(k+1)
                                x_lower_y_hi = (1-alpha_y_hi) * rr%lower(k) &
                                             + alpha_y_hi * rr%lower(k+1)
                                x_upper_y_hi = (1-alpha_y_hi) * rr%upper(k) &
                                              + alpha_y_hi * rr%upper(k+1)
                                x_lower_s = min(x_lower_y_low, x_lower_y_hi)
                                x_upper_s = max(x_upper_y_low, x_upper_y_hi)
                                if ((x_hi > x_lower_s) .and. &
                                    (x_low < x_upper_s)) then
                                      minlevel(i,j) = max(minlevel(i,j), &
                                                          rr%min_level)
                                      maxlevel(i,j) = max(maxlevel(i,j), &
                                                          rr%max_level)
                                      cycle iloop2
                                      endif
                              enddo  ! k                                  
                            endif ! rr%method
                        enddo iloop2
                    enddo jloop2
                endif ! rr%ixy == 2
        enddo rrloop
    
    
    do j=1,my
        do i=1,mx
         if (minlevel(i,j) > maxlevel(i,j)) then
              write(6,*) '*** Error: this should never happen!'
              write(6,*) '*** minlevel > maxlevel in flagregions'
              stop
         endif

         if (maxlevel(i,j) /= 0) then
             ! this point lies in at least one region, so may need
             ! to modify the exisiting flag at this point...
             if (level < minlevel(i,j)) then
                 ! Require refinement of this cell:
                 amrflags(i,j) = DOFLAG
             else if (level >= maxlevel(i,j)) then
                 ! Do not refine of this cell:
                 amrflags(i,j) = DONTFLAG
             ! else leave amrflags(i,j) alone.
             endif
         endif

        enddo 
    enddo 

#ifdef WHERE_AM_I
      write(*,*)"ending   flagregions2"
#endif

end subroutine flagregions2
