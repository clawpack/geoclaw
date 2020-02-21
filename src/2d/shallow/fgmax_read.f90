
subroutine fgmax_read(fgmax_unit,ifg)

    ! Read in data file describing any fixed grids.
    ! Modified for v5.7.0 to read all data from fgmax_grids.data
    ! This code now looks for the data for fgmax grid number ifg.
    
    ! The information is assumed to have the form:
    ! 
    ! tstart_max  # start time for monitoring this fgrid.
    ! tend_max    # end time for monitoring this fgrid.
    ! dt_check    # desired maximum time increment between updating max values.
    ! min_level_check # minimum level to check for monitoring values/arrivals
    ! arrival_tol           # tolerance for identifying arrival. 
    ! point_style   # indicating how fgmax points are specified on this grid
    
    ! if point_style==0 this is followed by:
    !   npts        # number of grid points
    !   x(1), y(1)  # first grid point
    !   ...
    !   x(npts), y(npts)  # last grid point
    
    ! OR if npts==0, instead of the list the next line after npts should be
    !   xy_fname  # path to file containing npts followed by list of points.

    ! if point_style==1:
    !   npts       # desired number of points on transect
    !   x1, y1     # first point
    !   x2, y2     # last point
    ! if point_style==2:
    !   nx, ny     # desired number of points in x and y
    !   x1, y1     # lower left corner of cartesian grid
    !   x2, y2     # upper right corner of cartesian grid
    ! if point_style==3:
    !   n12, n23   # number of points on sides 1-2 (also 3-4) and 2-3 (also 4-1)
    !   x1, y1     # first corner of quadrilateral
    !   x2, y2     # second corner
    !   x3, y3     # third corner
    !   x4, y4     # fourth corner
    ! if point_style==4:
    !   xy_fname   # a file with topo_type==3 format specifying Z = 0 or 1,
    !              # with 1 at fgmax points

    use fgmax_module
    use amr_module, only: mxnest
    use topo_module, only: read_topo_header

    implicit none
    integer, intent(in) :: ifg, fgmax_unit
    
    ! local
    character(150) :: fname
    integer :: k,i,j,point_style,n12,n23
    real(kind=8) :: x1,x2,y1,y2,yj
    real(kind=8) :: x3,x4,y3,y4,x14,y14,x23,y23,xi,eta
    type(fgrid), pointer :: fg
    logical :: foundFile
    character(len=150) :: fname2
    integer omp_get_max_threads, maxthreads
    integer :: clock_start, clock_finish, clock_rate
    real(kind=8), allocatable :: fg_row(:)
    real(kind=8) :: fg_y
    integer :: fg_npts_max, jj

    
    fg => FG_fgrids(ifg)   ! point to next element of array of fgrids
    read(fgmax_unit,*) fg%fgno
    read(fgmax_unit,*) fg%tstart_max
    read(fgmax_unit,*) fg%tend_max
    read(fgmax_unit,*) fg%dt_check  
    read(fgmax_unit,*) fg%min_level_check
    read(fgmax_unit,*) fg%arrival_tol
    read(fgmax_unit,*) fg%interp_method
    read(fgmax_unit,*) point_style
    fg%point_style = point_style
    
    if (point_style == 0) then
        read(fgmax_unit,*) fg%npts
        
        if (fg%npts .ne. 0) then
            ! points are also in fgmax_grids.data file:
            allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
            do k=1,fg%npts
                read(fgmax_unit,*) fg%x(k), fg%y(k)
                enddo
        else
            ! in this case fname2 is read, separate file containing x,y values:
            read(fgmax_unit,*) fname2
            write(6,*) 'Reading fgmax points from '
            write(6,*) '    ',trim(fname2)
            open(unit=FG_UNIT,file=trim(fname2),status='old')
            read(FG_UNIT,*) fg%npts
            write(6,*) 'npts = ',fg%npts

            allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
            do k=1,fg%npts
                read(FG_UNIT,*) fg%x(k), fg%y(k)
                enddo
            close(FG_UNIT)
        endif

    else if (point_style == 1) then
        read(fgmax_unit,*) fg%npts
        read(fgmax_unit,*) x1,y1
        read(fgmax_unit,*) x2,y2
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        do k=1,fg%npts
            fg%x(k) = x1 + (k-1)*(x2-x1)/(fg%npts - 1)
            fg%y(k) = y1 + (k-1)*(y2-y1)/(fg%npts - 1)
            enddo
    else if (point_style == 2) then
        read(fgmax_unit,*) fg%nx, fg%ny
        read(fgmax_unit,*) x1,y1
        read(fgmax_unit,*) x2,y2
        fg%xll = x1
        fg%yll = y1
        fg%dx = (x2-x1)/(fg%nx - 1)
        fg%dy = (y2-y1)/(fg%ny - 1)
        fg%npts = fg%nx*fg%ny
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        allocate(fg%index(fg%nx, fg%ny)) ! to store index into list of pts
        k = 0
        do j=1,fg%ny
            yj = y1 + (j-1)*(y2-y1)/(fg%ny - 1)
            do i=1,fg%nx
                k = k+1
                fg%x(k) = x1 + (i-1)*(x2-x1)/(fg%nx - 1)
                fg%y(k) = yj
                fg%index(i,j) = k
                enddo
            enddo
    else if (point_style == 3) then
        read(fgmax_unit,*) n12, n23
        read(fgmax_unit,*) x1,y1
        read(fgmax_unit,*) x2,y2
        read(fgmax_unit,*) x3,y3
        read(fgmax_unit,*) x4,y4
        fg%npts = n12*n23
        allocate(fg%x(1:fg%npts), fg%y(1:fg%npts))
        
        k = 0
        do j=1,n23
            eta = (j-1)/(n23-1.d0)
            x14 = x1 + eta*(x4-x1)
            y14 = y1 + eta*(y4-y1)
            x23 = x2 + eta*(x3-x2)
            y23 = y2 + eta*(y3-y2)
            do i=1,n12
                k = k+1
                xi = (i-1)/(n12-1.d0)
                fg%x(k) = x14 + xi*(x23 - x14)
                fg%y(k) = y14 + xi*(y23 - y14)
                enddo
            enddo
    else if (point_style == 4) then
        ! in this case read in a file name for a file in 
        ! topotype=3 format that contains 0 at points not to be used,
        ! and 1 (or maybe an index > 0) at the desired fgmax points
        read(fgmax_unit,*) fname2
        !close(unit=FG_UNIT)
        write(6,*) 'Reading fgmax points from '
        write(6,*) '    ',trim(fname2)
        
        call read_topo_header(fname2,3,fg%nx,fg%ny,fg%xll,fg%yll, &
                         fg%xhi,fg%yhi,fg%dx,fg%dy)

        open(unit=FG_UNIT,file=trim(fname2),status='old')
        ! Skip over header lines read above:
        do i=1,6
            read(FG_UNIT,*)
            enddo
        
        allocate(fg%index(fg%nx, fg%ny)) ! to store index into list of pts
        allocate(fg_row(fg%nx))          ! temporary for one row
        fg_npts_max = fg%nx * fg%ny
        allocate(fg%x(1:fg_npts_max), fg%y(1:fg_npts_max))  ! list of pts
        
        k = 0
        fg%index = 0
        do j=1,fg%ny
            jj = fg%ny - j + 1  ! since topo-style file goes from N to S
            read(FG_UNIT,*) (fg_row(i), i=1,fg%nx)
            fg_y = fg%yll + (jj-1)*fg%dy
            !write(6,*) '+++ fg_y = ',fg_y
            do i=1,fg%nx
                if (fg_row(i) > 0) then
                    k = k+1
                    fg%index(i,jj) = k
                    fg%x(k) = fg%xll + (i-1)*fg%dx
                    fg%y(k) = fg_y
                    if ((fg_row(i)>1) .and. (fg_row(i) /= k)) then
                        write(6,*) '*** fg index mismatch, ',fg_row(i),k
                        endif
                    endif
                enddo
            enddo
        fg%npts = k
        write(6,*) 'npts = ',fg%npts
            
    else
        write(6,*) '*** Unexpected value of point_style in fgmax_read: ',point_style
        write(6,*) '*** Note that format of fgmax file has changed, see documentation'
        stop
        endif


    ! needed to allocate fg%klist array for use on each patch by all threads:
    maxthreads = 1  ! if not using OpenMP
!$  maxthreads = omp_get_max_threads()

    ! needed for timing if printing at end of this routine:
    !call system_clock(clock_start,clock_rate)

    ! allocate and initialize arrays
    allocate(fg%valuemax(1:FG_NUM_VAL, 1:fg%npts))
    allocate(fg%levelmax(1:fg%npts))
    allocate(fg%aux(1:mxnest, 1:FG_NUM_AUX, 1:fg%npts))
    allocate(fg%auxdone(1:mxnest))
    allocate(fg%update_now(1:mxnest))
    allocate(fg%tmax(1:FG_NUM_VAL, 1:fg%npts))
    allocate(fg%t_last_updated(1:mxnest))
    allocate(fg%arrival_time(1:fg%npts))
    allocate(fg%klist(1:fg%npts, 0:maxthreads-1))

    fg%valuemax = FG_NOTSET
    fg%levelmax = 0
    fg%aux = FG_NOTSET
    fg%auxdone = .false.
    fg%update_now = .false.
    fg%tmax = FG_NOTSET
    fg%t_last_updated = fg%tstart_max - fg%dt_check - 1.d0
    fg%arrival_time = FG_NOTSET

    ! Set corners of bounding box.
    ! note that when point_style==4 these arrays are longer than fg%npts
    fg%x1bb = minval(fg%x(1:fg%npts))
    fg%x2bb = maxval(fg%x(1:fg%npts))
    fg%y1bb = minval(fg%y(1:fg%npts))
    fg%y2bb = maxval(fg%y(1:fg%npts))

    write(6,*) '+++ fgbb: ', fg%x1bb,fg%x2bb,fg%y1bb,fg%y2bb
    write(6,*) '+++ fgll: ', fg%xll, fg%yll, fg%nx, fg%ny
    write(6,*) '+++ fg%x(npts), fg%y(1): ',fg%x(fg%npts), fg%y(1)

    close(FG_UNIT)

    !call system_clock(clock_finish,clock_rate)
    !write(6,*) '+++ Done with fgmax_read, cpu time, clock_rate: ', &
    !        (clock_finish - clock_start), clock_rate
    !write(6,*) '+++ maxthreads, fg%npts, fg%klist_last: ', maxthreads, &
    !        fg%npts, fg%klist(fg%npts, maxthreads-1)

end subroutine fgmax_read
