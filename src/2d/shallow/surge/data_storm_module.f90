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
        integer :: num_casts

        ! size of latitude and longitude arrays
        integer :: mx, my

        ! Storm data, wind velocity in x and y, pressure and wind speed
        real(kind=8), allocatable :: pressure(:,:,:)
        real(kind=8), allocatable :: wind_u(:,:,:)
        real(kind=8), allocatable :: wind_v(:,:,:)

        ! Wind field latitude/longitude arrays
        real(kind=8), allocatable :: latitude(:)
        real(kind=8), allocatable :: longitude(:)
        ! Time steps from wind/pressure files in seconds
        integer, allocatable :: time(:)

    end type data_storm_type

    integer, private :: last_storm_index

    ! Time tracking tolerance allowance - allows for the beginning of the storm
    ! track to be close to but not equal the start time of the simulation
    real(kind=8), parameter :: TRACKING_TOLERANCE = 1d-10

contains
    ! ==========================================================================
    !  set_storm(storm_data_path, storm, storm_spec_type, log_unit)
    !    Initializes the storm type for an Oceanweather, Inc type data derived 
    !    storm and calls subroutines based on input filetype ie (nc or ascii)
    ! ==========================================================================
    subroutine set_storm(storm_data_path, storm, storm_spec_type, log_unit)
        implicit none
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        if (-3 <= storm_spec_type .and. storm_spec_type < 0) then
            select case(storm_spec_type)
            case(-2) ! netcdf
                call set_netcdf_storm(storm_data_path, storm, storm_spec_type, &
                                  log_unit)
            case(-3) ! ascii fixed width
                call set_data_storm(storm_data_path, storm, storm_spec_type, &
                                    log_unit)
            end select
        end if
    end subroutine set_storm

    ! ==========================================================================
    !  set_netcdf_storm(storm_data_path, storm, storm_spec_type, log_unit)
    !    Initializes the storm type for an Oceanweather, Inc type data derived 
    !    storm that is saved as a netcdf format
    ! ==========================================================================
    subroutine set_netcdf_storm(storm_data_path, storm, storm_spec_type, log_unit)
#ifdef NETCDF
    use netcdf
#endif
        use amr_module, only: t0, rinfinity
        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! Local Storage
        ! Set up for dimension info, dims are lat/lon/time
        integer :: num_dims, x_dim_id, y_dim_id, t_dim_id, dim_ids(3)
        character (len=10) :: x_dim_name, y_dim_name, t_dim_name

        ! Length of dims, mx, my, mt
        integer :: mx, my, mt

        ! Set up for variable info, variables are pressure, lat, lon, time
        integer :: nvars, var_type, var_id
        character (len=13) :: var_name

        ! Netcdf file id and counter for looping
        integer :: nc_fid, n
#ifdef NETCDF
        if (.not. module_setup) then

            ! Open data file
            print *, 'Reading storm data file ', storm_data_path

            ! Open file and get file ID
            call check_netcdf_error(nf90_open(storm_data_path, nf90_nowrite, nc_fid))
            ! Read number of dimensions and variables in nc file
            call check_netcdf_error(nf90_inquire(nc_fid, num_dims, nvars))

            ! Get the dimension names and sizes from the file
            call get_dim_info(nc_fid, num_dims, x_dim_id, x_dim_name, mx, &
            y_dim_id, y_dim_name, my, t_dim_id, t_dim_name, mt)

            ! allocate arrays in storm object
            allocate(storm%pressure(mx, my, mt))
            allocate(storm%wind_u(mx, my, mt))
            allocate(storm%wind_v(mx, my, mt))
            allocate(storm%longitude(mx))
            allocate(storm%latitude(my))
            allocate(storm%time(mt))
            ! Set up lat/lon lengths in storm object for use in linear interpolation of time
            storm%mx = mx
            storm%my = my

            ! Number of time steps in data
            storm%num_casts = mt

            ! Fill out variable data/info
            do n = 1, nvars
                ! read file for one variable and parse the data in storm object
                call check_netcdf_error(nf90_inquire_variable(nc_fid, n, var_name, var_type, num_dims, dim_ids))
                if ('PRESSURE' == Upper(var_name)) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%pressure))
                elseif(ANY((/'LON      ','LONGITUDE'/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%longitude))
                elseif(ANY((/'LAT     ', 'LATITUDE'/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%latitude))
                elseif(ANY((/'TIME', 'T   '/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%time))
                elseif(ANY((/'U     ', 'WIND_U', 'UU    '/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%wind_u))
                elseif(ANY((/'V     ', 'WIND_V', 'VV    '/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%wind_v))
                end if
            end do

        ! Close file to stop corrupting the netcdf files
        call check_netcdf_error(nf90_close(nc_fid))
        end if
#endif
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

        ! Write out a surge file after discussing the data with kyle
        if (.not. module_setup) then

            module_setup = .true.
        end if
    end subroutine set_netcdf_storm

    ! ==========================================================================
    !  set_data_storm(storm_data_path, storm, storm_spec_type, log_unit)
    !    Initializes the storm type for an Oceanweather, Inc type data derived 
    !    storm that is saved as a fixed width fortran format
    ! ==========================================================================
    subroutine set_data_storm(storm_data_path, storm, storm_spec_type, log_unit)
        use amr_module, only: t0, rinfinity
        implicit none

        ! Subroutine I/O
        character(len=*), intent(in) :: storm_data_path
        character(len=(len(storm_data_path))) :: WIND_FILE, PRESSURE_FILE
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! Local Storage
        ! Set up for dimension info, dims are lat/lon/time
        integer :: initial_time, mt
        integer :: mx, my, yr, mo, da, hr, minute
        real(kind=8) :: swlat, swlon, dx, dy
        integer :: wunit=7, punit=8, ext_pos
        character(len=3) :: wind_ext='WIN', pres_ext='PRE', file_ext
        logical :: file_exists
        ! Set up for variable info, variables are pressure, lat, lon, time
        integer :: nvars, var_type, var_id
        character (len=13) :: var_name
        
        if (.not. module_setup) then
            ext_pos = scan(trim(storm_data_path), ".", BACK=.true.)
            file_ext = Upper(storm_data_path(ext_pos+1:len_trim(storm_data_path)))
           
            if (ANY((/'WIN', 'WND'/) == file_ext)) then
                WIND_FILE = storm_data_path
                PRESSURE_FILE = storm_data_path(1:ext_pos)//'PRE'
            else if (file_ext == 'PRE') then
                PRESSURE_FILE = storm_data_path
                WIND_FILE = storm_data_path(1:ext_pos)//'WIN'
                inquire(file=trim(WIND_FILE), exist=file_exists)
                if (.NOT. file_exists) then
                    WIND_FILE = storm_data_path(1:ext_pos)//'WND'
                end if
            end if

            ! Open data file
            print *, 'Reading storm data file ', storm_data_path, WIND_FILE, PRESSURE_FILE

            ! Load headers for setting up the rest of the data
            call initialize_storm_data(WIND_FILE, my, mx, mt, swlat, &
                              swlon, dy, dx, initial_time)

            ! ! Determine array sizes from the file
            ! call get_array_sizes(WIND_FILE, mt, yr, mo, da, hr, minute, &
            !                      my, mx, total_time)
            
            ! allocate arrays in storm object
            allocate(storm%pressure(mx, my, mt))
            allocate(storm%wind_u(mx, my, mt))
            allocate(storm%wind_v(mx, my, mt))
            allocate(storm%longitude(mx))
            allocate(storm%latitude(my))
            allocate(storm%time(mt))
            ! Set up lat/lon lengths in storm object for use in linear interpolation of time
            storm%mx = mx
            storm%my = my

            ! Number of time steps in data
            storm%num_casts = mt

            ! Fill out variable data/info
            call fill_data_arrays(WIND_FILE, PRESSURE_FILE, my, mx, mt, swlat, swlon,&
                                 dy, dx, storm, initial_time)
        end if

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

        ! Write out a surge file after discussing the data with kyle
        if (.not. module_setup) then

            module_setup = .true.
        end if
    end subroutine set_data_storm

    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function storm_location(t,storm) result(location)
        ! used in other potential data storm types, requires eye location
        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(inout) :: storm

        ! Output
        real(kind=8) :: location(2)

        stop "Data-derived storm are not yet implemented!"

    end function storm_location

    ! ==========================================================================
    !  storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function storm_direction(t, storm) result(theta)
        ! Used for other types of data storm types, requires eye location
        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in) :: storm

        stop "Data-derived storm are not yet implemented!"

    end function storm_direction


    ! ==========================================================================
    !  Use the 1980 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_HWRF_fields(maux, mbc, mx, my, xlower, ylower,    &
                          dx, dy, t, aux, wind_index,           &
                          pressure_index, storm)

  
        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(data_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        stop "HWRF data input is not yet implemented!"

    end subroutine set_HWRF_fields

    ! ==========================================================================
    ! read_headers() Opens fixed width format file and reads the header
    ! for use in parsing the data structure
    ! ==========================================================================
    subroutine initialize_storm_data(WIND_FILE, my, mx, mt, swlat, swlon, dy, &
                                     dx,  initial_time)
        implicit none
        character(len=*), intent(in) :: WIND_FILE

        ! Output arguments, time in seconds
        integer, intent(out) ::  my, mx, mt
        real(kind=8), intent(out) :: swlat, swlon, dx, dy

        ! Local Storage
        integer, parameter :: iunit=9
        integer :: start_yr, start_mo, start_da, start_hr, start_min, dt
        integer :: end_yr, end_mo, end_da, end_hr, end_min
        integer :: yr1, mo1, da1, hr1, minute1, total_time, initial_time, final_time
        integer :: yr2, mo2, da2, hr2, minute2, timestep2, timestep1
        integer :: i, j

        ! Read in start and end dates of file from first header line
        open(iunit, file=WIND_FILE, status='old', action='read')
       
        read(iunit, "(t56, i4, i2, i2, i2, t71, i4, i2, i2, i2)") &
                    start_yr, start_mo, start_da, start_hr, end_yr, end_mo, &
                    end_da, end_hr
        initial_time = seconds_from_epoch(start_yr, start_mo, start_da, start_hr, 00)
        final_time = seconds_from_epoch(end_yr, end_mo, end_da, end_hr, 00)
        total_time = final_time - initial_time
        
        ! Read second header from file to obtain array sizes and first timestep
        read(iunit, '(t6,i4,t16,i4,t23,f6.0,t32,f6.0,t44,f8.0, t58, f8.0,t69, i4,i2,i2,i2, i2)') &
             my, mx, dx, dy, swlat, swlon, yr1, mo1, da1, hr1, minute1

        ! Skip the wind velocity matrices first u then v
        print *, "Skipping first matrix (u)..."
        do i = 1, (mx * my / 8)
            read(iunit, *) ! Skip first matrix
        end do
        print *, "First matrix (u) skipped."
        if (mod(mx * my,8) /= 0) then
            print *, "Reading line after 1st matrix..."
            read(iunit, *) 
        else
            print *, 'proceeding to reading the matrix'
        end if
        print *, "Skipping second matrix (v)..."
        do j = 1, (mx * my / 8)
            read(iunit, *) ! Skip second matrix
        end do
        print *, "Second matrix (v) skipped."
        if (mod(mx * my,8) /= 0) then
            print *, "Reading line after 2nd matrix..."
            
        else
            print *, 'proceeding to reading the header'
        end if
        read(iunit, '(t69, i4,i2,i2,i2, i2)') yr2, mo2, da2, hr2, minute2
        close(iunit)

        timestep1 = seconds_from_epoch(yr1, mo1, da1, hr1, minute1)
        timestep2 = seconds_from_epoch(yr2, mo2, da2, hr2, minute2)

        dt = timestep2 - timestep1
        mt = (total_time/dt) + 1


    end subroutine initialize_storm_data       
    
    ! ==========================================================================
    ! get_array_sizes() Reads the 2nd header data from the datafile and 
    ! saves the array sizes for allocation inside the data_storm_type
    ! ==========================================================================
    ! subroutine get_array_sizes(WIND_FILE, mt, yr, mo, da, hr, minute, &
    !                            my, mx, total_time)
    !     implicit none
    !     ! Input arguments
    !     integer, intent(in) :: yr, mo, da, hr, minute
    !     integer, intent(in) :: my, mx, total_time
    !     character(len=*), intent(in) :: WIND_FILE

    !     ! Output argument
    !     integer, intent(out) :: mt

    !     ! Local storage
    !     integer :: iunit = 8
    !     integer :: yr2, mo2, da2, hr2, min2
    !     integer :: time1, time2, i, j, dt
    !     integer, dimension(8) :: values

    !     ! Error Handling
    !     integer :: ios
    !     character(len=100) :: line_after_2nd_matrix
    !     open(iunit, file=WIND_FILE, status='old', action='read', iostat=ios)
    !     if (ios /= 0) then
    !         print *, "Error opening file:", ios
    !         stop
    !     end if
    !     print *, "File opened successfully."

    !     ! Read time data from first two headers/timestep
    !     print *, "Reading first header..."
    !     read(iunit, *) 
    !     print *, "First header read."

    !     print *, "Reading second header..."
    !     read(iunit, *)
    !     print *, "Second header read."

    !     ! Skip the wind velocity matrices first u then v
    !     print *, "Skipping first matrix (u)..."
    !     do i = 1, (mx * my / 8)
    !         read(iunit, *) ! Skip first matrix
    !     end do
    !     print *, "First matrix (u) skipped."
    !     if (mod(mx * my,8) /= 0) then
    !         print *, "Reading line after 1st matrix..."
    !         read(iunit, *) 
    !     else
    !         print *, 'proceeding to reading the matrix'
    !     end if
    !     ! print *, "Reading line separating matrices..."
    !     ! read(iunit, *) ! Skip line separating matrices
    !     ! print *, "Line separating matrices read."

    !     print *, "Skipping second matrix (v)..."
    !     do j = 1, (mx * my / 8)
    !         read(iunit, *) ! Skip second matrix
    !     end do
    !     print *, "Second matrix (v) skipped."
    !     if (mod(mx * my,8) /= 0) then
    !         print *, "Reading line after 2nd matrix..."
    !         read(iunit, *) line_after_2nd_matrix
    !         print *, "Line after 2nd matrix:", line_after_2nd_matrix
    !     else
    !         print *, 'proceeding to reading the header'
    !     end if

    !     ! Read time data from second header/timestep
    !     print *, "ABOUT TO READ 1st HEADER"
    !     ! read(iunit, *) line_after_2nd_matrix
    !     read(iunit, '(t69, i4,i2,i2,i2, i2)') yr2, mo2, da2, hr2, min2
    !     ! print *, "Second header/timestep read: ", line_after_2nd_matrix
    !     !++++++++++++ WHY NOT JUST CALCULATE THE NUMBER OF DAYS BETWEEN TIME STEPS AND THEN DIVIDE BY HOURS?????+++++++++++++++
    !     close(iunit)
    !     print *, "File closed."

    !     ! open(iunit, file=WIND_FILE, status='old', action='read', iostat=ios)
        
    !     ! ! Read time data from first two headers/timestep
    !     ! read(iunit, *) 
    !     ! read(iunit, *)
    !     ! ! Skip the wind velocity matrices first u then v
    !     ! do i = 1, (mx * my / 8)
    !     !     read(iunit, *) ! Skip first matrix
    !     ! end do
    !     ! read(iunit, *) ! Skip line separating matrices
    !     ! do j = 1, (mx * my / 8)
    !     !     read(iunit, *) ! Skip second matrix
    !     ! end do
    !     ! read (iunit, *) ! Skip line after 2nd matrix
      
    !     ! ! Read time data from second header/timestep
    !     !  print *, 'ABOUT TO READ 1st HEADER'
    !     ! read(iunit, '(t69, i4,i2,i2,i2, i2)') yr2, mo2, da2, hr2, min2
    !     ! close(iunit)
       

    !     dt =  parse_dates(yr, mo, da, hr, minute, yr2, mo2, da2, hr2, min2)
    !     print *, 'DT and TOTAL TIME', dt, total_time
    !     mt = (total_time/dt) + 1 ! total number of timesteps
    ! end subroutine get_array_sizes
    
    ! ==========================================================================
    ! fill_data_arrays() reads the data files and fills out the storm object
    ! and it's dataarrays
    ! ==========================================================================
    subroutine fill_data_arrays(WIND_FILE, PRESSURE_FILE, my, mx, mt, swlat, swlon,  &
                                dy, dx, storm,  initial_time)
        ! Read in the data file and parse the arrays to be passed to other calculations
        implicit none
        ! Input arguments
        type(data_storm_type) :: storm
        character(len=*), intent(in) :: WIND_FILE, PRESSURE_FILE
        integer, intent(in) :: mt, my, mx, initial_time 
        real(kind=8), intent(in) :: swlat, swlon, dx, dy
        
        ! Local storage
        integer :: num_lat, num_lon, i, j, n, current_timestep
        integer :: yr, mo, da, hr, minute
        integer :: wunit = 700, punit = 800
        open(wunit, file=WIND_FILE, status='old', action='read')
        open(punit, file=PRESSURE_FILE, status='old', action='read')
        ! Skip file headers
        read(wunit, *)
        read(punit, *)
        ! allocate(storm%time(mt))
        do n = 1, mt
            
            read(wunit, '(t69, i4,i2,i2,i2, i2)') yr, mo, da, hr, minute
            
            current_timestep = seconds_from_epoch(yr, mo, da, hr, minute) - initial_time
            print *, '##### Current time ####',n, current_timestep
            storm%time(n) = current_timestep
            
            read(wunit, '(8f10.0)') ((storm%wind_u(i,j, n),i=1,mx),j=1,my)
            read(wunit, '(8f10.0)') ((storm%wind_v(i,j, n),i=1,mx),j=1,my)
            
        end do
        print *, 'Finished WIND, starting PRESSURE'
        
        do n = 1, mt
            read(punit, *) ! Skip header line since we have it from above
            read(punit, '(8f10.0)') ((storm%pressure(i,j, n),i=1,mx),j=1,my)
            ! Multiply each pressure value by 100 to convert from Pa to hPa
            storm%pressure(:,:,n) = storm%pressure(:,:,n) * 100.0
        end do
        
        do i = 1, my
            storm%latitude(i) = swlat + i * dy
        end do
        do j = 1, mx
            storm%longitude(j) = swlon + j * dx
        end do
        
        
    end subroutine fill_data_arrays

    ! ==========================================================================
    ! get_dim_info() pulls dim names and lengths from the nc file
    !
    ! ==========================================================================

    subroutine get_dim_info(nc_file, ndims, x_dim_id, x_dim_name, mx, &
        y_dim_id, y_dim_name, my, t_dim_id, t_dim_name, mt)
#ifdef NETCDF
    use netcdf
#endif
        implicit none
        integer, intent(in) :: nc_file, ndims
        integer, intent(out) :: x_dim_id, y_dim_id, mx, my, t_dim_id, mt
        character (len = *), intent(out) :: x_dim_name, y_dim_name, t_dim_name
        integer :: m_tmp, n
        character(20) :: dim_name_tmp
#ifdef NETCDF
        ! get indices to start at for reading netcdf within domain
        do n=1, ndims
            call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                n, dim_name_tmp, m_tmp))
            if (ANY((/ 'LON      ','LONGITUDE','X        ' /) == Upper(dim_name_tmp))) then
                x_dim_name = dim_name_tmp
                mx = m_tmp
                x_dim_id = n
            else if (ANY((/ 'LAT     ','LATITUDE','Y       ' /) == Upper(dim_name_tmp))) then
                y_dim_name = dim_name_tmp
                my = m_tmp
                y_dim_id = n
            else if (ANY((/'TIME', 'T   '/) == Upper(dim_name_tmp))) then
                t_dim_name = dim_name_tmp
                mt = m_tmp
                t_dim_id = n
            end if
        end do
#endif
    end subroutine get_dim_info

    ! ==========================================================================
    !  storm_index(t,storm)
    !    Finds the index of the next storm data point
    ! ==========================================================================
    integer pure function storm_index(t, storm) result(index)

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
        if (last_storm_index == storm%num_casts + 1) then
            index = storm%num_casts + 1
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
                do index=last_storm_index+1,storm%num_casts
                    if (t < storm%time(index)) then
                        found = .true.
                        exit
                    endif
                enddo
                ! Assume we have gone past last forecast time
                if (.not. found) then
                    index = storm%num_casts + 1
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


    ! ==========================================================================
    ! set_owi_fields()
    ! Fills out data for current time step and current patch
    ! ==========================================================================
    subroutine set_owi_fields(maux, mbc, mx, my, xlower, ylower,    &
                              dx, dy, t, aux, wind_index,           &
                              pressure_index, storm)
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
        logical :: convert_height=.true.

        ! wind and pressure arrays for current time step
        real(kind=8) :: wind_tu(storm%mx,storm%my), wind_tv(storm%mx, storm%my), pressure_t(storm%mx,storm%my)
        ! get the wind and pressure arrays for the current timestep
        call get_storm_time(storm, t, wind_tu, wind_tv, pressure_t, mx, my)
        ! Loop over every point in the patch and fill in the data
        do j=1-mbc, my+mbc
            y = ylower + (j-0.5d0) * dy
            do i=1-mbc, mx+mbc
                x = xlower + (i-0.5d0) * dx
                call interp_array_data(storm, x, y, wind_tu, u_value)
                aux(wind_index, i, j) = u_value
                call interp_array_data(storm, x, y, wind_tv, v_value)
                aux(wind_index + 1, i, j) = v_value
                call interp_array_data(storm, x, y, pressure_t, p_value)
                aux(pressure_index, i, j) = p_value
            end do
        end do
    end subroutine set_owi_fields


    ! ==========================================================================
    ! get_storm_time()
    ! Calculates the wind and pressure fields for the current time step
    ! ==========================================================================
    subroutine get_storm_time(storm, t, wind_tu, wind_tv,  pressure_t,  mx, my)
        implicit none
        ! Subroutine IO
        type(data_storm_type), intent(in) :: storm
        integer, intent(in) :: mx, my
        real(kind=8), intent(in) :: t
        real(kind=8), intent(inout) :: wind_tu(storm%mx,storm%my), wind_tv(storm%mx, storm%my),pressure_t(storm%mx,storm%my)

        ! Local storage
        integer :: i, j, k, mx_orig, my_orig
        real(kind=8) :: tn, tnm, weight
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
            else if (i > storm%num_casts + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif
        if (i == storm%num_casts + 1) then
            i = i -1
            wind_tu = storm%wind_u(:,:,i)
            wind_tv = storm%wind_v(:,:,i)
            pressure_t = storm%pressure(:,:,i)
        else
            tn = storm%time(i)
            tnm = storm%time(i - 1)
            weight = (t - tnm) / (tn - tnm)
            wind_tu = (storm%wind_u(:,:,i) - storm%wind_u(:,:,i - 1)) * weight + &
                      storm%wind_u(:,:, i-1)
            wind_tv = (storm%wind_v(:,:,i) - storm%wind_v(:,:,i - 1)) * weight + &
                      storm%wind_v(:,:, i-1)
            pressure_t = (storm%pressure(:,:,i) - storm%pressure(:,:,i - 1)) * weight + &
                          storm%pressure(:,:,i - 1)
        end if

    end subroutine get_storm_time


    ! ==========================================================================
    ! interp_array_data() Obtains wind and pressure data for each patch
    ! Uses bilinear interpolation to get the data
    ! ==========================================================================
    subroutine interp_array_data(storm, x, y, interp_array, value)
        implicit none
        ! Subroutine I/O
        type(data_storm_type), intent(inout) :: storm
        real(kind=8), intent(inout) :: x, y
        real(kind=8), intent(in) :: interp_array(size(storm%longitude), size(storm%latitude))
        real(kind=8), intent(out) :: value

        ! Local storage
        real(kind=8) :: q11, q12, q21, q22
        real(kind=8) :: llon, llat, ulon, ulat
        real(kind=8) :: storm_dx, storm_dy ! initial wind field dx and dy
        integer :: xidx_low, xidx_high, yidx_low, yidx_high

        ! Get the data resolution
        storm_dx = storm%longitude(2) - storm%longitude(1)
        storm_dy = storm%latitude(2) - storm%latitude(1)

        ! Duplicate the boundary at the bottom corner
        if (x < minval(storm%longitude) .and. y < maxval(storm%latitude)) then
            value = interp_array(1, 1)
        ! Duplicate the boundary at the top corner
        elseif (x > maxval(storm%longitude) .and. y > maxval(storm%latitude)) then
            value = interp_array(storm%mx, storm%my)
        ! if x is to the west of the boundary but y is inside
        elseif (x < minval(storm%longitude)) then
            call find_nearest(x, y, llon, llat, storm, xidx_low, yidx_low)
            value = interp_array (xidx_low, yidx_low)
        ! if x is to the east of the boundary but y is inside
        elseif (x > maxval(storm%longitude)) then
            call find_nearest(x, y, llon, llat, storm, xidx_low, yidx_low)
            value = interp_array(xidx_low, yidx_low)
        ! if y is south of the boundary but x is inside
        elseif (y < minval(storm%latitude)) then
            call find_nearest(x, y, llon, llat, storm, xidx_low, yidx_low)
            value = interp_array(xidx_low, yidx_low)
        ! if y is north of the boundary but x is inside
        elseif (y > maxval(storm%latitude)) then
            call find_nearest(x, y, llon, llat, storm, xidx_low, yidx_low)
            value = interp_array(xidx_low, yidx_low)
        ! x and y both inside boundary
        else
            call find_nearest(x - storm_dx, y - storm_dy, llon, llat, storm, xidx_low, yidx_low)
            call find_nearest(x + storm_dx, y + storm_dy, ulon, ulat, storm, xidx_high, yidx_high)
            ! Find the values at the corners
            q11 = interp_array(xidx_low, yidx_low)
            q12 = interp_array(xidx_low, yidx_high)
            q21 = interp_array(xidx_high, yidx_low)
            q22 = interp_array(xidx_high, yidx_high)

            ! Calculate the value at the center of the box using bilinear interpolation
            value = (q11 * (ulon - x) * (ulat -y) + &
                     q21 * (x - llon) * (ulat - y) + &
                     q12 * (ulon - x) * (y - llat) + &
                     q22 * (x - llon) * (y - llat)) / ((ulon - llon) * (ulat - llat) + 0.00)
        end if

    end subroutine interp_array_data

    ! ==========================================================================
    ! find_nearest() finds nearest value to x and y for interpolation points
    ! Finds the nearest actual point to the patch values using minimum distance
    ! ==========================================================================
    pure subroutine find_nearest(x, y, lon, lat, storm, xidx, yidx)
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

    subroutine check_netcdf_error(ios)
#ifdef NETCDF
    use netcdf
#endif
        implicit none

        integer, intent(in) :: ios
#ifdef NETCDF
        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if
#endif
    end subroutine check_netcdf_error
    
    function seconds_from_epoch(year, month, day, hour, minute) result(seconds)
        integer, intent(in) :: year, month, day, hour, minute
        integer, dimension(12) :: days_in_month=[31, 28, 31, 30, 31, &
                                                 30, 31, 31, 30, 31, 30, 31]
        integer :: seconds, leap_days, total_days
        
        total_days = (year - 1970) * 365

        leap_days = (year - 1968)/4 - (year - 1900)/4 + (year - 1600)/400

        total_days = total_days + leap_days
        if (mod(year, 4) == 0 .and. (mod(year,100) /= 0 &
            .or. mod(year, 400) == 0).and.month >2) THEN
                total_days = total_days + 1
        endif
        total_days = total_days + sum(days_in_month(1:month-1)) + day

        seconds = (total_days*86400) + (hour*3600) + (minute*60)
       
       
    end function seconds_from_epoch
    

    function Upper(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i

        DO i = 1,LEN(s1)
           ch = s1(i:i)
           IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
           s2(i:i) = ch
        END DO
    END function Upper
end module data_storm_module


   ! ==========================================================================
   !
   !
   ! ==========================================================================



