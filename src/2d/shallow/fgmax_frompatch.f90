

subroutine fgmax_frompatch(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
           xlower,ylower,level,time,time_afterstep)

    ! Do the new fixed grid stuff on all fgrids, updating 
    ! based on the patch passed in.

    use fgmax_module
    use geoclaw_module, only: sea_level, dry_tolerance

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux,level
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower,time,time_afterstep

    real(kind=8), allocatable, dimension(:,:) :: fg_values
    !integer, allocatable, dimension(:) :: fg_klist ! replaced by fg%klist
    type(fgrid), pointer :: fg
    integer :: ifg,k,mv, fg_klist_length, indexk
    real(kind=8) :: h,B,eta
    real(kind=8) :: x1,x2,y1,y2,x,y,xupper,yupper,eps
    integer :: i1,i2,j1,j2
    integer :: clock_start, clock_finish, clock_rate
    integer :: mythread, omp_get_thread_num
    integer :: i_fg,i1_fg,i2_fg,  j_fg,j1_fg,j2_fg
    
      mythread = 0 ! if not using OpenMP
!$    mythread = omp_get_thread_num()

    !write(6,*) '++++ mythread = ', mythread

    !print *, '+++ FG_DEBUG = ', FG_DEBUG

    if (FG_num_fgrids == 0) then
        return
        endif 

    eps = min(dx,dy)*0.01d0  ! buffer width to avoid rounding error issues
    
    !write(61,*) '++++ In frompatch, level,mx,my,xlower,ylower'
    !write(61,*) level,mx,my,xlower,ylower
    
    ifg_loop: do ifg=1,FG_num_fgrids
        fg => FG_fgrids(ifg)
        !write(6,*) '++++ in frompatch, tstart,tend: ',fg%tstart_max, fg%tend_max
        !write(6,*) '++++ ifg = ',ifg
        !write(6,*) '++++ fg%xNbb = ',fg%x1bb,fg%x2bb
        !write(61,*) '++++ frompatch xNbb', fg%x1bb, fg%x2bb

        if (FG_DEBUG) then
            write(61,61) ifg,level,time
 61         format('---------- In fgmax_frompatch ----------',/, &
               'ifg = ',i2,' level = ',i1,' time = ',d16.6)
            if (time_afterstep <= minval(fg%t_last_updated)+fg%dt_check) then
                write(61,68) time, minval(fg%t_last_updated)
     68         format('++++ Skipping update at t = ', e20.11,' min t_last = ',&
                       e20.11)
            else
                write(61,67) time, minval(fg%t_last_updated)
     67         format('++++ Doing update at t = ', e20.11,' min t_last = ', &
                       e20.11)
                endif
            endif

!        if ((time >= fg%tstart_max) .and. (time <= fg%tend_max) .and. &
!                (level >= fg%min_level_check) .and. &
!                (level >= minval(fg%levelmax)) .and. &
!                (time_afterstep > minval(fg%t_last_updated)+fg%dt_check)) &
!                then
!            ! Otherwise this level won't update any fg%valuemax elements.
!            !write(61,*) '+++ Setting fg_values'

!        if ((time < fg%tstart_max) .or. (time > fg%tend_max) .or. &
!                (time_afterstep < fg%t_last_updated(level)+fg%dt_check) .or. &
!                (level < fg%min_level_check) .or. &
!                (level < fg%min_levelmax_checked)) then 
                
        if (.not. fg%update_now(level)) then
        
!            if (level == 6) then
!                write(6,601) fg%min_level_check, &
!                 fg%min_levelmax_checked, time_afterstep, &
!                 fg%t_last_updated(level)+fg%dt_check
! 601            format('++++ not now 6: ',2i3,2d26.17)
!                endif
                
            !write(6,*) '++++ nothing:time,level,mythread:',time,level,mythread
            !write(6,*) '++++ tstart,tend: ',fg%tstart_max, fg%tend_max
            !write(6,*) '++++ ',fg%min_level_check, minval(fg%levelmax)
            
            cycle ifg_loop  ! nothing to do, go to next fgmax grid
            endif
            
            
        !===  unindented...
        call system_clock(clock_start,clock_rate)
        
        ! Determine intersection of this patch with bounding box of fgrid:

        xupper = xlower+mx*dx
        yupper = ylower+my*dy
        x1 = max(xlower, fg%x1bb) - eps
        x2 = min(xupper, fg%x2bb) + eps
        y1 = max(ylower, fg%y1bb) - eps
        y2 = min(yupper, fg%y2bb) + eps
        if (FG_DEBUG) then
            write(61,*) 'xlower,xupper: ',xlower,xupper
            write(61,*) 'ylower,yupper: ',ylower,yupper
            write(61,*) 'x1bb, x2bb: ',fg%x1bb, fg%x2bb
            write(61,*) 'y1bb, y2bb: ',fg%y1bb, fg%y2bb
            write(61,*) 'Intersection x1, x2:',x1,x2
            write(61,*) 'Intersection y1, y2:',y1,y2
            endif

        ! If this grid does not intersect the bounding box of fg, do nothing:
        if ((x1 > x2) .or. (y1 > y2)) then
            if (FG_DEBUG) then
                write(61,*) '+++ No intersection found!'
                endif
            !return
            !write(6,*) '++++ No intersection'
            cycle ifg_loop ! go to next fgmax grid
            endif
                    
        
        ! No longer create mask, just use i1,i2, j1,j2 directly
        ! only values(:,i1:i2,j1:j2) on this patch are needed for fgmax
        ! adding 0.5d0 cellwidth necessary for bilinear interp option
        i1 = max(int((x1 - xlower + 0.5d0*dx) / dx), 0)
        i2 = min(int((x2 - xlower + 0.5d0*dx) / dx) + 1, mx+1)
        j1 = max(int((y1 - ylower + 0.5d0*dy) / dy), 0)
        j2 = min(int((y2 - ylower + 0.5d0*dy) / dy) + 1, my+1)

        if (FG_DEBUG) then
            write(61,*) 'patch intersecting fgrid: i1,i2: ',i1,i2
            write(61,*) 'patch intersecting fgrid: j1,j2: ',j1,j2
            endif

        ! Create list of k values for which (fg%x(k),fg%y(k)) is in patch:
        
        fg_klist_length = 0

        if ((fg%point_style == 2) .or. (fg%point_style == 4)) then
            ! fgmax points are on masked grid giving index of each point,
            ! so we don't need to search through all fg%npts points.
            ! instead just loop over part of masked grid intersecting patch:

            ! June 2020: RJL Fixed index error arising when x1 close to fg%xll
            !    due to fact that int(A) = sgn(A)*floor(abs(A)) returns 0
            !    when A is between -1 and +1 regardless of sign.
            !    Hence i1_fg might have been 2 when it should have been 1.
            !    And similarly in y direction.
            
            ! reset eps since now we are computing ranges of indices 
            ! on the fgmax grid that lie within this patch:
            eps = min(fg%dx, fg%dy) * 0.01d0
            
            if (x1 <= fg%xll) then
                i1_fg = 1
            else
                i1_fg = max(int((x1 - fg%xll) / fg%dx + 2), 1)
            endif
            i2_fg = min(int((x2 - fg%xll) / fg%dx + 1), fg%nx)

            if (y1 <= fg%yll) then
                j1_fg = 1
            else
                j1_fg = max(int((y1 - fg%yll) / fg%dy + 2), 1)
            endif
            j2_fg = min(int((y2 - fg%yll) / fg%dy + 1), fg%ny)
            
            
            !write(6,602) x1,x2,y1,y2
            !write(6,602) xlower,xupper,ylower,yupper
            !write(6,602) fg%x1bb, fg%x2bb, fg%y1bb, fg%y2bb
            !write(6,602) fg%xll, fg%xll+(fg%nx-1)*fg%dx, &
            !             fg%yll, fg%yll+(fg%ny-1)*fg%dy
            !write(6,603) i1,i2,j1,j2, i1_fg, i2_fg, j1_fg, j2_fg
 602        format('+++ xy: ',4f13.6)
 603        format('+++ ij: ',8i7)

 
            do j_fg=j1_fg,j2_fg
                do i_fg=i1_fg,i2_fg
                    if (fg%index(i_fg,j_fg) > 0) then
                        k = fg%index(i_fg,j_fg)
                        fg_klist_length = fg_klist_length+1
                        fg%klist(fg_klist_length,mythread) = k
                        !write(6,*) '++++ k, fg%x(k), fg%y(k): ',k, fg%x(k), fg%y(k)

                        endif
                    enddo
                enddo
        else
            ! we must loop over all fgmax points to look for those on patch:
            do k=1,fg%npts
                if ((fg%x(k) >= x1) .and. (fg%x(k) <= x2) .and. &
                    (fg%y(k) >= y1) .and. (fg%y(k) <= y2)) then
                        fg_klist_length = fg_klist_length+1
                        fg%klist(fg_klist_length,mythread) = k
                        endif
                enddo
            endif
        
        !write(6,*) '+++ klength = ',   fg_klist_length
         
        if (fg_klist_length == 0) then
            !write(6,*) '++++ t with klength=0: ', time
            cycle ifg_loop  ! go to next fgmax grid
            endif

        ! now fg_values is indexed by indexk where fg%klist(indexk,mythread) = k
        ! and k is the index into the fgmax grid
        ! so fg_values can be much shorter than fg%npts:
        allocate(fg_values(FG_NUM_VAL, fg_klist_length))
        
        !write(6,*) '++++ t, klength, kfirst, klast: ', &
        !        time, fg_klist_length, fg_klist(1), fg_klist(fg_klist_length)
        
    
        ! Interpolate from q on patch to desired values fg_values on fgrid.
        call fgmax_interpolate(mx,my,meqn,mbc,maux,q,aux, &
             dx,dy,xlower,ylower,ifg,level,fg_values, &
             i1,i2,j1,j2,fg%klist(1:fg%npts,mythread), &
             fg_klist_length,fg%npts)

        if (FG_DEBUG) then
            write(61,*) '+++ updating level: ',level
            write(61,*) '+++ in fgmax_frompatch, count = ',fg_klist_length 
            endif

    
        do indexk=1,fg_klist_length
            k = fg%klist(indexk,mythread)
            ! fg_values is set only at points k where the fgrid intersects the
            ! patch, now stored in this list
            ! Also: only update valuemax if the current patch is at least
            ! as high a level as the patch last used to update:
            if ((level >= fg%levelmax(k))) then
                !fg%t_last_updated(k) = time
                do mv=1,FG_NUM_VAL
                    !print *,'+++ updating fg%valuemax(mv,k),fg_values(mv,indexk):',&
                    !      fg%valuemax(mv,k),fg_values(mv,indexk)
                    if ((level > fg%levelmax(k)) .or. &
                          (fg_values(mv,indexk) > fg%valuemax(mv,k))) then
                        fg%valuemax(mv,k) = fg_values(mv,indexk)
                        ! also keep track of time maximum happened:
                        fg%tmax(mv,k) = time  
                        endif
                    enddo
                fg%levelmax(k) = level
                endif

            ! keep track of arrival time...
            ! This might not be correct when there is subsidence or uplift?
            h = fg_values(1,indexk)
            B = fg%aux(level,1,k)
            eta = h+B
            if ((level >= fg%min_level_check) .and. &
               (h>dry_tolerance) .and. &
               (eta > sea_level + fg%arrival_tol) .and. &
               (fg%arrival_time(k) == FG_NOTSET)) then
                    fg%arrival_time(k) = time
                endif

            enddo
            
        deallocate(fg_values)
            
        !call system_clock(clock_finish,clock_rate)
        !write(6,*) '+++ t, cpu time , klength, mythread : ', &
        !    time,(clock_finish - clock_start), fg_klist_length, mythread

        enddo ifg_loop

end subroutine fgmax_frompatch
