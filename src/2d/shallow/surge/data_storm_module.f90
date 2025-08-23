! ==============================================================================
! data_storm_module
!
! Module contains routines for constructing a wind and pressure field based on a
! provided set of data files.
!
! ==============================================================================
!                   Copyright (C) Clawpack Developers 2017
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================
module data_storm_module
        
    implicit none
    save

    logical, private :: module_setup = .false.
    logical, private :: DEBUG = .false.

    ! Model storm type definition
    type data_storm_type

        ! Total number of wind/pressure fields
        ! integer :: num_casts
        integer :: num_regions

        ! Time offset in string - parsed for times by not used directly
        character(len=19) :: time_offset

        ! Storm data, wind velocity in x and y, pressure and wind speed
        real(kind=8), allocatable :: pressure(:,:,:)
        real(kind=8), allocatable :: wind_u(:,:,:)
        real(kind=8), allocatable :: wind_v(:,:,:)

        ! Wind field dimension arrays
        real(kind=8), allocatable :: latitude(:)
        real(kind=8), allocatable :: longitude(:)
        integer, allocatable :: time(:)

        ! Paths to data
        character(len=256), allocatable :: paths(:)

    end type data_storm_type

    integer, private :: last_storm_index

    ! Ramping settings
    integer :: window_type
    real(kind=8) :: window(4), ramp_width

    ! Time tracking tolerance allowance - allows for the beginning of the storm
    ! track to be close to but not equal the start time of the simulation
    real(kind=8), parameter :: TRACKING_TOLERANCE = 1d-10

contains

    ! ==========================================================================
    !  Setup and Data I/O
    ! ==========================================================================

    ! === set_storm ============================================================
    ! Initializes the storm type for an HWRF type data derived storm that is
    ! saved as a netcdf format
    subroutine set_storm(storm_data_path, storm, storm_spec_type, log_unit)
        
#ifdef NETCDF
        use netcdf
#endif

        use amr_module, only: t0, tfinal
        use utility_module, only: to_upper, check_netcdf_error
        use utility_module, only: seconds_from_epoch, ISO_time_format

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! General data
        integer, parameter :: data_unit = 10
        integer :: i, io_status, file_format, num_data_files
        integer :: mt, mx, my, time(6, 2)

        ! ASCII
        real(kind=8) :: sw_lon, sw_lat, dx, dy
        
        ! NetCDF
        integer :: nc_fid, num_dims, num_vars, var_id
        character(len=16) :: name, dim_names(3), var_names(3)

        if (.not.module_setup) then
            ! Open data file
            print *,'Reading storm data file ', storm_data_path
            open(unit=data_unit, file=storm_data_path, status='old',        &
                 action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening storm data file. status = ", io_status
                stop
            endif

            read(data_unit, *) ! Comment line
            read(data_unit, "(a)") storm%time_offset
            read(data_unit, "(i2)") file_format
            read(data_unit, "(i2)") num_data_files
            read(data_unit, "(i2)") window_type
            read(data_unit, *) ramp_width

            write(log_unit, "('time_offset = ',a)") storm%time_offset
            write(log_unit, "('Format = ',i2)") file_format
            write(log_unit, "('Num files = ',i2)") num_data_files
            write(log_unit, "('window_type = ',i2)") window_type
            write(log_unit, "('ramp_width = ',d16.8)") ramp_width

            ! If not provided, then we use the file extents after being read
            if (window_type > 0) then
                read(data_unit, *) window
            else
                read(data_unit, *)
            end if
            read(data_unit, *)
            
            if (DEBUG) then
                print "('time_offset = ',a)", storm%time_offset
                print "('Format = ',i2)", file_format
                print "('Num files = ',i2)", num_data_files
            end if

            ! Format specific
            select case(file_format)
                case(1) ! ASCII/OWI
                    read(data_unit, *)
                    read(data_unit, *)
                case(2) ! NetCDF
                    read(data_unit, *)
                    read(data_unit, *) dim_names
                    read(data_unit, *) var_names
                    read(data_unit, *)
                    write(log_unit, "('dims = ',a)") dim_names
                    write(log_unit, "('vars = ',a)") var_names

                    if (DEBUG) then
                        print "('dims = ',a)", dim_names
                        print "('vars = ',a)", var_names
                    end if
                case default
                    print *, "Storm data format", file_format," not available."
                    stop
            end select

            ! Read rest of data file and close
            allocate(storm%paths(num_data_files))
            read(data_unit, *)
            write(log_unit, *) "File Paths = "
            do i=1, num_data_files
                read(data_unit, "(a)") storm%paths(i)
                write(log_unit, "(' ', a)") storm%paths(i)
                if (DEBUG) then
                    print "(' ', a)", storm%paths(i)
                end if
            end do
            close(data_unit)

            ! Calculate number of seconds from epoch - needed below
            read(storm%time_offset, ISO_time_format) time(:, 1)

            ! Read actual data
            select case(file_format)
                case(1)
                    ! Load headers for setting up the rest of the data
                    call read_OWI_ASCII_header(storm%paths(1), mx, my, mt,   &
                                                               sw_lon, sw_lat, &
                                                               dy, dx)
                    write(log_unit, *) "Storm header info:"
                    write(log_unit, "('  mx,mx,mt = ', 3i4)") mx, my, mt
                    write(log_unit, "('  sw = ', 2D25.16)") sw_lon, sw_lat
                    write(log_unit, "('  delta = ', 2D25.16)") dx, dy

                    if (DEBUG) then
                        print *, "Storm header info:"
                        print "('  mx,mx,mt = ', 3i4)", mx, my, mt
                        print "('  sw = ', 2D25.16)", sw_lon, sw_lat
                        print "('  delta = ', 2D25.16)", dx, dy
                    end if

                    ! allocate arrays in storm object
                    allocate(storm%pressure(mx, my, mt))
                    allocate(storm%wind_u(mx, my, mt))
                    allocate(storm%wind_v(mx, my, mt))
                    allocate(storm%longitude(mx))
                    allocate(storm%latitude(my))
                    allocate(storm%time(mt))

                    ! Fill out variable data/info
                    print *, "Reading pressure file ", storm%paths(1)
                    print *, "Reading wind file ", storm%paths(2)
                    call read_OWI_ASCII(storm%paths(2), storm%paths(1),        &
                                        mx, my, mt, sw_lat, sw_lon, dy, dx,    &
                                        storm, seconds_from_epoch(time(:, 1)))

                case(2)
#ifdef NETCDF                    
                    ! Open file and get file ID
                    ! :TODO: Only read in times that are between t0 and tfinal
                    print *, "Reading storm NetCDF file ", storm%paths(1)
                    call check_netcdf_error(nf90_open(storm%paths(1), nf90_nowrite, nc_fid))
                    ! Check dim/var number
                    call check_netcdf_error(nf90_inquire(nc_fid, num_dims, num_vars)) 
                    if (num_dims /= 3 .and. num_vars /= 3) then
                        print *, "Invalid number of dimensions/variables."
                        print *, "  num_dims = ", num_dims
                        print *, "  num_vars = ", num_vars
                        stop
                    end if

                    ! Get dimensions
                    call check_netcdf_error(nf90_inq_dimid(nc_fid, dim_names(1), var_ID))
                    call check_netcdf_error(nf90_inquire_dimension(nc_fid, var_ID, name, mx))
                    call check_netcdf_error(nf90_inq_dimid(nc_fid, dim_names(2), var_ID))
                    call check_netcdf_error(nf90_inquire_dimension(nc_fid, var_ID, name, my))
                    call check_netcdf_error(nf90_inq_dimid(nc_fid, dim_names(3), var_ID))
                    call check_netcdf_error(nf90_inquire_dimension(nc_fid, var_ID, name, mt))

                    ! allocate arrays in storm object
                    allocate(storm%pressure(mx, my, mt))
                    allocate(storm%wind_u(mx, my, mt))
                    allocate(storm%wind_v(mx, my, mt))
                    allocate(storm%longitude(mx))
                    allocate(storm%latitude(my))
                    allocate(storm%time(mt))
                    
                    write(log_unit, *) "Storm header info:"
                    write(log_unit, "('  mx,mx,mt = ', 3i4)") mx, my, mt

                    if (DEBUG) then
                        print *, "Storm header info:"
                        print "('  mx,mx,mt = ', 3i4)", mx, my, mt
                    end if

                    ! Dimensions
                    call check_netcdf_error(nf90_inq_varid(nc_fid, dim_names(1), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%longitude))
                    call check_netcdf_error(nf90_inq_varid(nc_fid, dim_names(2), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%latitude))
                    call check_netcdf_error(nf90_inq_varid(nc_fid, dim_names(3), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%time))

                    ! Convert time to seconds from offset
                    storm%time(:) = storm%time(:) - seconds_from_epoch(time(:, 1))

                    ! Wind
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_names(1), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%wind_u))
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_names(2), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%wind_v))
                    ! ! Pressure
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_names(3), var_ID))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_ID, storm%pressure))

                    ! Close file to stop corrupting the netcdf files
                    call check_netcdf_error(nf90_close(nc_fid))
#else
                    print *, "GeoClaw was not compiled with NetCDF support."
                    stop
#endif
            end select

            ! Make sure first time step in setrun is inside the times in the data
            if (t0 - storm%time(1) < -TRACKING_TOLERANCE) then
                print *, 'Start time', t0, " is outside of the tracking"
                print *, 'tolerance range with the track start'
                print *, storm%time(1), '.'
                stop
            end if

            last_storm_index = 2
            last_storm_index = storm_index(t0, storm)
            if (last_storm_index == -1) then
                print *, 'Forecast not found for time ', t0, '.'
                stop
            end if

            ! Set bounds on window to file bounds if not provided
            if (window_type == 0) then
                window = [storm%longitude(1), storm%latitude(1),            &
                          storm%longitude(mx), storm%latitude(my)]
            end if
            write(log_unit, "('window = ',4d16.8)") window

            module_setup = .true.
        endif

    end subroutine set_storm


    ! === read_OWI_ASCII_header ================================================
    ! Opens ASCII formatted OWI file and reads the header for use in parsing the
    ! data structure
    subroutine read_OWI_ASCII_header(pressure_file, mx, my, mt, swlon, swlat, &
                                                    dy, dx)
        
        use utility_module, only: seconds_from_epoch

        implicit none

        ! Subroutine I/O
        character(len=*), intent(in) :: pressure_file

        ! Output arguments
        integer, intent(out) ::  mx, my, mt
        real(kind=8), intent(out) :: swlat, swlon, dx, dy

        ! Format strings
        character(len=*), parameter :: header_format =          &
                                            "(t56,i4,i2,i2,i2,t71,i4,i2,i2,i2)"
        character(len=*), parameter :: full_info_format =       &
                                            "(t6,i4,t16,i4,t23,f6.0,t32,"   // &
                                            "f6.0,t44,f8.0,t58,f8.0,t69,"   // &
                                            "i4,i2,i2,i2,i2)"
        character(len=*), parameter :: time_info_format =       &
                                            "(t69,i4,i2,i2,i2,i2)"
        
        ! Local
        integer, parameter :: OWI_unit = 156
        integer :: i, time(5, 2), total_time, dt


        ! Read in start and end dates of file from first header line
        open(OWI_unit, file=pressure_file, status='old', action='read')
       
        ! Total number of time in seconds of the dataset
        read(OWI_unit, header_format) time(1:4, 1), time(1:4, 2) 
        total_time =   seconds_from_epoch(time(1:4, 2))           &
                     - seconds_from_epoch(time(1:4, 1))
        
        ! Read second header from file to obtain array sizes and first timestep
        read(OWI_unit, full_info_format) my, mx, dx, dy, swlat, swlon, time(:, 1)

        ! ! Would need to multiply by 2 if using the wind file
        do i=1, (mx * my / 8 + merge(0, 1, mod(mx * my, 8) == 0))
            read(OWI_unit, *)
        end do
        read(OWI_unit, time_info_format) time(:, 2)
        
        close(OWI_unit)

        dt = seconds_from_epoch(time(:, 2)) - seconds_from_epoch(time(:, 1))
        mt = (total_time / dt) + 1

    end subroutine read_OWI_ASCII_header  
    

    ! === read_OWI_ASCII =======================================================
    ! reads the data files and fills out the storm object and it's dataarrays
    subroutine read_OWI_ASCII(wind_file, pressure_file, mx, my, mt,         &
                                                        swlat, swlon,       &
                                                        dx, dy, storm,      &
                                                        seconds_from_offset)
        
        use utility_module, only: seconds_from_epoch
        
        implicit none
        
        ! Input arguments
        type(data_storm_type) :: storm
        character(len=*), intent(in) :: wind_file, pressure_file
        integer, intent(in) :: my, mx, mt, seconds_from_offset 
        real(kind=8), intent(in) :: swlat, swlon, dx, dy
        
        ! Local storage
        integer, parameter :: wind_unit = 700, pressure_unit = 800
        integer :: i, j, n, time(5)

        ! Open both storm forcing files
        open(wind_unit, file=wind_file, status='old', action='read')
        open(pressure_unit, file=pressure_file, status='old', action='read')

        ! Skip file headers
        read(wind_unit, *)
        read(pressure_unit, *)
        ! Loop over all timesteps 
        do n = 1, mt
            ! Read each time from the next array
            read(wind_unit, '(t69, i4,i2,i2,i2, i2)') time
            storm%time(n) = seconds_from_epoch(time) - seconds_from_offset

            read(wind_unit, '(8f10.0)') ((storm%wind_u(i,j, n),i=1,mx),j=1,my)
            read(wind_unit, '(8f10.0)') ((storm%wind_v(i,j, n),i=1,mx),j=1,my)
            
            read(pressure_unit, *) ! Skip header line since we have it from above
            read(pressure_unit, '(8f10.0)') ((storm%pressure(i,j, n),i=1,mx),j=1,my) 
            ! Convert from Pa to hPa
            storm%pressure(:,:,n) = storm%pressure(:,:,n) * 100.0
        end do
        close(wind_unit)
        close(pressure_unit)

        ! Calculate lat/lon from SW corner and resolution
        ! :TODO: do we need to store these?
        do i = 1, my
            storm%latitude(i) = swlat + i * dy
        end do
        do j = 1, mx
            storm%longitude(j) = swlon + j * dx
        end do
        
    end subroutine read_OWI_ASCII


    ! ==========================================================================
    !  Management functions
    ! ==========================================================================

    ! === storm_index(t,storm) =================================================
    ! Finds the index of the next storm data point
    integer function storm_index(t, storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        ! Figure out where we are relative to the last time we checked for the
        ! index (stored in last_storm_index)
        ! Check if we are already beyond the end of the last forecast time
        if (last_storm_index == size(storm%time) + 1) then
            index = size(storm%time) + 1
        else
            t0 = storm%time(last_storm_index - 1)
            t1 = storm%time(last_storm_index)

            ! Check to see if we are close enough to the current index to just
            ! use that, tolerance is based on TRACKING_TOLERANCE
            if ((abs(t0 - t) < TRACKING_TOLERANCE) .or.   &
                (abs(t1 - t) < TRACKING_TOLERANCE) .or.   &
                (t0 < t .and. t < t1)) then

                index = last_storm_index
            else if ( t1 < t ) then
                found = .false.
                do index=last_storm_index+1,size(storm%time)
                    if (t < storm%time(index)) then
                        found = .true.
                        exit
                    endif
                enddo
                ! Assume we have gone past last forecast time
                if (.not. found) then
                    index = size(storm%time) + 1
                endif
            else
                ! Fail gracefully
                if (last_storm_index == 2) then
                    index = -1
                else
                    do index=last_storm_index-1,2,-1
                        if (storm%time(index-1) < t) exit
                    enddo
                endif
            endif
        endif

    end function storm_index


    ! === storm_location(t,storm) ==============================================
    ! Interpolate location of hurricane in the current time interval
    function storm_location(t,storm) result(location)
        ! used in other potential data storm types, requires eye location
        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(inout) :: storm

        ! Output
        real(kind=8) :: location(2)

        stop "Storm location for data storms is not supported."

    end function storm_location


    ! === storm_direction ======================================================
    ! Angle off of due north that the storm is traveling
    real(kind=8) function storm_direction(t, storm) result(theta)
        ! Used for other types of data storm types, requires eye location
        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in) :: storm

        stop "Storm direction for data storms is not supported."

    end function storm_direction


    ! === set_data_fields ======================================================
    ! Fills out data for current time step and current patch
    subroutine set_data_fields(maux, mbc, mx, my, xlower, ylower,    &
                               dx, dy, t, aux, wind_index,           &
                               pressure_index, storm)

        use geoclaw_module, only: Pa => ambient_pressure, pi
        use utility_module, only: point_in_rectangle, rects_intersect

        implicit none
        ! subroutine i/o
        integer, intent(in) :: maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t

        ! Storm descrption
        type(data_storm_type), intent(inout) :: storm

        ! Array storing wind and pressure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: u_value, v_value, p_value
        integer :: i, j
        real(kind=8) :: x, y

        ! Temporally interpolated fields
        real(kind=8) :: wind_u(size(storm%longitude),size(storm%latitude))
        real(kind=8) :: wind_v(size(storm%longitude), size(storm%latitude))
        real(kind=8) :: pressure(size(storm%longitude),size(storm%latitude))

        ! Ramping window support
        real(kind=8) :: d, ramp, xhi, yhi

        ! Check to see if this anything on this grid is inside of the window
        ! Patch boundaries [xlower, xhi] x [ylower, yhi] \in window + ramp_width
        xhi = xlower + (mx + mbc - 0.5d0) * dx
        yhi = ylower + (my + mbc - 0.5d0) * dy
        if (rects_intersect([xlower, ylower, xhi, yhi],                        &
                        [window(1) - ramp_width, window(2) - ramp_width,       &
                         window(3) + ramp_width, window(4) + ramp_width])) then

            ! Get the wind and pressure arrays for the current timestep
            call interp_time(storm, t, wind_u, wind_v, pressure)

            ! Loop over every point in the patch and fill in the data
            do j=1-mbc, my+mbc
                y = ylower + (j-0.5d0) * dy
                do i=1-mbc, mx+mbc
                    x = xlower + (i-0.5d0) * dx

                    ! Inside of window + ramp_width
                    if (point_in_rectangle([x, y], [window(1) - ramp_width,   &
                                                    window(2) - ramp_width,   &
                                                    window(3) + ramp_width,   &
                                                    window(4) + ramp_width])) then

                        aux(wind_index, i, j) = spatial_interp(storm, x, y, wind_u)
                        aux(wind_index + 1, i, j) = spatial_interp(storm, x, y, wind_v)
                        aux(pressure_index, i, j) = spatial_interp(storm, x, y, pressure)

                        ! Ramp function, hardcoded right now
                        d = min(x - window(1), window(3) - x,         &
                                y - window(2), window(4) - y)
                        if (-ramp_width < d .and. d < 0.d0) then
                            ramp = 0.5d0 * (1.d0 + sin(-pi / ramp_width * d + pi / 2.0))
                        else if (d > 0.d0) then
                            ramp = 1.d0
                        else
                            ramp = 0.d0
                        end if
                        aux(pressure_index,i,j) =                              &
                                    Pa + (aux(pressure_index,i,j) - Pa) * ramp
                        aux(wind_index:wind_index+1,i,j) =                     &
                                    aux(wind_index:wind_index+1,i,j) * ramp

                    else
                        aux(wind_index:wind_index + 1, i, j) = 0.d0
                        aux(pressure_index, i, j) = Pa
                    end if
                end do
            end do
        else
            aux(wind_index:wind_index+1, :, :) = 0.d0
            aux(pressure_index, :, :) = Pa
        end if
    end subroutine set_data_fields

    ! === set_ascii_fields =====================================================
    ! Stub that points to the generic set_data_fields
    subroutine set_ascii_fields(maux, mbc, mx, my, xlower, ylower,    &
                                dx, dy, t, aux, wind_index,           &
                                pressure_index, storm)
        implicit none
        integer, intent(in) :: maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        call set_data_fields(maux, mbc, mx, my, xlower, ylower, dx, dy, t,     &
                             aux, wind_index, pressure_index, storm)

    end subroutine set_ascii_fields


    ! === set_netcdf_fields =====================================================
    ! Stub that points to the generic set_data_fields
    subroutine set_netcdf_fields(maux, mbc, mx, my, xlower, ylower,    &
                                 dx, dy, t, aux, wind_index,           &
                                 pressure_index, storm)
        implicit none
        integer, intent(in) :: maux, mbc, mx, my
        real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        call set_data_fields(maux, mbc, mx, my, xlower, ylower, dx, dy, t,     &
                             aux, wind_index, pressure_index, storm)

    end subroutine set_netcdf_fields


    ! === interp_time ==========================================================
    ! Interpolate storm arrays to get current time
    subroutine interp_time(storm, t, wind_u, wind_v, pressure)
        
        implicit none
        ! Subroutine IO
        type(data_storm_type), intent(in) :: storm
        real(kind=8), intent(in) :: t
        real(kind=8), intent(inout) :: wind_u(size(storm%longitude),        &
                                              size(storm%latitude))
        real(kind=8), intent(inout) :: wind_v(size(storm%longitude),        &
                                              size(storm%latitude))
        real(kind=8), intent(inout) :: pressure(size(storm%longitude),      &
                                                size(storm%latitude))

        ! Local storage
        integer :: i
        real(kind=8) :: tn, tnm, weight

        ! Get indices of bracketing storm fields
        last_storm_index = 2
        i = storm_index(t, storm)
        last_storm_index = i

        ! List of possible error conditions
        if (i <= 1) then
            if (i == 0) then
                print *,"Invalid storm forecast requested for t = ",t
                print *,"Time requested is before any forecast data."
                print *,"    first time = ",storm%time(1)
                print *,"   ERROR = ",i
                stop
            else if (i > size(storm%time) + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif
        if (i == size(storm%time) + 1) then
            i = i -1
            wind_u = storm%wind_u(:,:,i)
            wind_v = storm%wind_v(:,:,i)
            pressure = storm%pressure(:,:,i)
        else
            tn = storm%time(i)
            tnm = storm%time(i - 1)
            weight = (t - tnm) / (tn - tnm)
            wind_u = (storm%wind_u(:,:,i) - storm%wind_u(:,:,i - 1)) * weight + &
                      storm%wind_u(:,:, i-1)
            wind_v = (storm%wind_v(:,:,i) - storm%wind_v(:,:,i - 1)) * weight + &
                      storm%wind_v(:,:, i-1)
            pressure = (storm%pressure(:,:,i) - storm%pressure(:,:,i - 1)) * weight + &
                        storm%pressure(:,:,i - 1)
        end if

    end subroutine interp_time


    ! === spatial_intrp ========================================================
    ! Spatially interpolate interp_array onto (x,y) point
    pure real(kind=8) function spatial_interp(storm, x, y, interp_array) result(value)
        
        implicit none
        
        ! Subroutine I/O
        type(data_storm_type), intent(in) :: storm
        real(kind=8), intent(in) :: x, y, interp_array(size(storm%longitude),  &
                                                       size(storm%latitude))

        ! Local storage
        real(kind=8) :: q11, q12, q21, q22
        real(kind=8) :: llon, llat, ulon, ulat
        real(kind=8) :: storm_dx, storm_dy ! initial wind field dx and dy
        integer :: xidx_low, xidx_high, yidx_low, yidx_high

        ! Get the data resolution
        storm_dx = storm%longitude(2) - storm%longitude(1)
        storm_dy = storm%latitude(2) - storm%latitude(1)

        ! Add a edge case handling that reflects the boundary
        if (x < minval(storm%longitude) .or. x > maxval(storm%longitude) .or. &
            y < minval(storm%latitude) .or. y > maxval(storm%latitude)) then
                ! Set safe defaults in case of out-of-bounds
                xidx_low = 1
                xidx_high = 1
                yidx_low = 1
                yidx_high = 1

                if (x < minval(storm%longitude)) then
                    xidx_low = 1
                    xidx_high = 1
                elseif (x > maxval(storm%longitude)) then
                    xidx_low = size(storm%longitude)
                    xidx_high = size(storm%longitude)
                endif
                if (y < minval(storm%latitude)) then
                    yidx_low = 1
                    yidx_high = 1
                elseif (y > maxval(storm%latitude)) then
                    yidx_low = size(storm%latitude)
                    yidx_high = size(storm%latitude)
                endif

                value = interp_array(xidx_low, yidx_low)
        
        else
            call find_nearest(storm, x - storm_dx, y - storm_dy, llon, llat,  xidx_low, yidx_low)
            call find_nearest(storm, x + storm_dx, y + storm_dy, ulon, ulat,  xidx_high, yidx_high)
            ! Find the values at the corners
            q11 = interp_array(xidx_low, yidx_low)
            q12 = interp_array(xidx_low, yidx_high)
            q21 = interp_array(xidx_high, yidx_low)
            q22 = interp_array(xidx_high, yidx_high)
                        
            ! Calculate the value at the center of the box using bilinear interpolation
            value = (q11 * (ulon - x) * (ulat - y) + &
                     q21 * (x - llon) * (ulat - y) + &
                     q12 * (ulon - x) * (y - llat) + &
                     q22 * (x - llon) * (y - llat)) / ((ulon - llon) * (ulat - llat) + 0.00)
            
        end if
        
    end function spatial_interp

    ! === find_nearest =========================================================
    ! Finds nearest value to x and y for interpolation points. Finds the nearest
    ! actual point to the patch values using minimum distance
    pure subroutine find_nearest(storm, x, y, lon, lat,  xidx, yidx)
        implicit none
        ! Subroutine I/O
        type(data_storm_type), intent(in) :: storm
        real(kind=8), intent(in) :: x, y
        real(kind=8), intent(out) :: lon, lat
        integer, intent(out) :: xidx, yidx

        xidx = minloc(abs(storm%longitude - x), dim=1)
        yidx = minloc(abs(storm%latitude - y), dim=1)
        lon = storm%longitude(xidx)
        lat = storm%latitude(yidx)
    end subroutine find_nearest

end module data_storm_module