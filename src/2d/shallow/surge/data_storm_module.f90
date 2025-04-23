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
        integer :: num_regions

        ! Landfall in string - parsed for times by not used directly
        character(len=19) :: landfall

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
    !    storm and calls subroutines based on input filetype  (nc or ascii)
    ! ==========================================================================
    subroutine set_storm(storm_data_path, storm, storm_spec_type, log_unit)

        implicit none
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit
        
        if (-2 <= storm_spec_type .and. storm_spec_type < 0) then
            select case(storm_spec_type)
            case(-1) ! HWRF
                call set_HWRF_storm(storm_data_path, storm, storm_spec_type, log_unit)
            case(-2) ! OWI
                call set_OWI_storm(storm_data_path, storm, storm_spec_type, log_unit)
            end select
        end if
    end subroutine set_storm

    ! ==========================================================================
    !  set_HWRF_storm(storm_data_path, storm, storm_spec_type, log_unit)
    !    Initializes the storm type for an HWRF type data derived storm that
    !    is saved as a netcdf format
    ! ==========================================================================
    subroutine set_HWRF_storm(storm_data_path, storm, storm_spec_type, log_unit)
        
        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        stop "HWRF storm support not implemented."

    end subroutine set_HWRF_storm

    ! ==========================================================================
    !  set_owi_storm(storm_data_path, storm, storm_spec_type, log_unit)
    !    Initializes the storm type for an Oceanweather, Inc type data derived 
    !    storm that is saved as a netcdf format
    ! ==========================================================================
    subroutine set_OWI_storm(storm_data_path, storm, storm_spec_type, log_unit)

#ifdef NETCDF
        use netcdf
#endif
        use amr_module, only: t0, rinfinity

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! Locals
        integer :: file_format, io_status
        integer, parameter :: data_unit = 10
        ! integer, parameter :: OWI_unit = 156s

        ! Locals
        integer :: i, time(6, 2), dt, total_time, seconds_from_landfall
        integer :: mt, mx, my

        ! ASCII / NWS12
        character(len=*), parameter :: ISO_time_Format = "(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)"
        ! character(len=*), parameter :: header_format = "(t56,i4,i2,i2,i2,t71,i4,i2,i2,i2)"
        ! character(len=*), parameter :: full_info_format = "(t6,i4,t16,i4,t23,f6.0,t32,f6.0,t44,f8.0,t58,f8.0,t69,i4,i2,i2,i2,i2)"
        ! character(len=*), parameter :: time_info_format = "(t69,i4,i2,i2,i2,i2)"
        ! real(kind=8) :: sw_lon, sw_lat, dx, dy
        character(len=256), allocatable :: wind_files(:), pressure_files(:)

        ! Temp conversion
        ! integer :: yr, mo, da, hr, minute, seconds
        
        real(kind=8) :: swlat, swlon, dx, dy

        ! Local Storage
        ! Set up for dimension info, dims are lat/lon/time
        integer :: num_dims, x_dim_id, y_dim_id, t_dim_id, dim_ids(3)
        character (len=10) :: x_dim_name, y_dim_name, t_dim_name

        ! Set up for variable info, variables are pressure, lat, lon, time
        integer :: nvars, var_type, var_id
        character (len=13) :: var_name

        ! Netcdf file id and counter for looping
        integer :: nc_fid, n
        character(len=256) :: nc_path

        if (.not. module_setup) then
            ! Open data file
            print *,'Reading storm data file ', storm_data_path
            open(unit=data_unit, file=storm_data_path, status='old',        &
                 action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening storm data file. status = ", io_status
                stop
            endif

            read(data_unit, *) ! Comment line
            read(data_unit, "(a)") storm%landfall
            read(data_unit, "(i2)") file_format
            read(data_unit, *)
            read(data_unit, *)
            write(log_unit, "('Landfall = ',a)") storm%landfall
            write(log_unit, "('Format = ',i1)") file_format

            ! ASCII / NWS12 input files
            if (file_format == 1) then
                read(data_unit, '(i2)') storm%num_regions
                write(log_unit, "('Num regions = ',i2)") storm%num_regions
                if (storm%num_regions > 1) then
                    print *, "More than 1 regions is not implemented."
                    stop
                end if
                
                ! Read in file paths
                allocate(wind_files(storm%num_regions))
                allocate(pressure_files(storm%num_regions))
                do i=1, storm%num_regions
                    read(data_unit, "(a)") pressure_files(i)
                    read(data_unit, "(a)") wind_files(i)
                end do
                close(data_unit)

                ! Calculate number of seconds from epoch
                read(storm%landfall, ISO_time_format) time(:, 1)

                ! Load headers for setting up the rest of the data
                call read_OWI_ASCII_header(pressure_files(1), mx, my, mt,   &
                                                              swlon, swlat, &
                                                              dy, dx)
                write(log_unit, *) "Storm header info:"
                write(log_unit, *) mx, my, mt
                write(log_unit, *) swlon, swlat, dx, dy

                ! Set up lat/lon lengths in storm object for use in linear interpolation of time
                storm%mx = mx
                storm%my = my

                ! Number of time steps in data
                storm%num_casts = mt

                ! allocate arrays in storm object
                allocate(storm%pressure(mx, my, mt))
                allocate(storm%wind_u(mx, my, mt))
                allocate(storm%wind_v(mx, my, mt))
                allocate(storm%longitude(mx))
                allocate(storm%latitude(my))
                allocate(storm%time(mt))

                ! Fill out variable data/info
                print *, "Reading wind file ", wind_files(1)
                print *, "Reading pressure file ", pressure_files(1)
                call read_OWI_ASCII(wind_files(1), pressure_files(1),       &
                                        mx, my, mt, swlat, swlon, dy, dx,   &
                                        storm, seconds_from_epoch(time(1:5, 1)))

            ! NetCDf / NWS13 input file
            else if (file_format == 2) then
#ifdef NETCDF
                ! Read rest of data file and close
                read(data_unit, "(a)") nc_path
                close(data_unit)

                ! Open file and get file ID
                print *, "Reading storm NetCDF file ", nc_path
                call check_netcdf_error(nf90_open(nc_path, nf90_nowrite, nc_fid))
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
                    if ('PRESSURE' == upper(var_name)) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%pressure))
                    elseif(ANY((/'LON      ','LONGITUDE'/) == upper(var_name))) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%longitude))
                    elseif(ANY((/'LAT     ', 'LATITUDE'/) == upper(var_name))) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%latitude))
                    elseif(ANY((/'TIME', 'T   '/) == upper(var_name))) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%time))
                    elseif(ANY((/'U     ', 'WIND_U', 'UU    '/) == upper(var_name))) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%wind_u))
                    elseif(ANY((/'V     ', 'WIND_V', 'VV    '/) == upper(var_name))) then
                        call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                        call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%wind_v))
                    end if
                end do

                ! Close file to stop corrupting the netcdf files
                call check_netcdf_error(nf90_close(nc_fid))
#else
                print *, "GeoClaw was not compiled with NetCDF support needed for"
                print *, "OWI NetCDF support."
                stop
#endif                
            else
                print *, "Invalid file format ", file_format,"."
                stop
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

            module_setup = .true.
        end if

    end subroutine set_OWI_storm

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
    ! read_OWI_ASCII_header() Opens ASCII formatted OWI file and reads the 
    ! header for use in parsing the data structure
    ! ==========================================================================
    subroutine read_OWI_ASCII_header(pressure_file, mx, my, mt, swlon, swlat, &
                                                    dy, dx)
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
    
    ! ==========================================================================
    ! read_OWI_ASCII() reads the data files and fills out the storm object
    ! and it's dataarrays
    ! ==========================================================================
    subroutine read_OWI_ASCII(wind_file, pressure_file, mx, my, mt,         &
                                                        swlat, swlon,       &
                                                        dx, dy, storm,      &
                                                        seconds_from_landfall)
        
        implicit none
        
        ! Input arguments
        type(data_storm_type) :: storm
        character(len=*), intent(in) :: wind_file, pressure_file
        integer, intent(in) :: my, mx, mt, seconds_from_landfall 
        real(kind=8), intent(in) :: swlat, swlon, dx, dy
        
        ! Local storage
        integer, parameter :: wind_unit = 700, pressure_unit = 800
        integer :: num_lat, num_lon, i, j, n, time(5)

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
            storm%time(n) = seconds_from_epoch(time) - seconds_from_landfall

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
    ! get_dim_info() pulls dim names and lengths from the nc file
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
            if (ANY((/ 'LON      ','LONGITUDE','X        ' /) == upper(dim_name_tmp))) then
                x_dim_name = dim_name_tmp
                mx = m_tmp
                x_dim_id = n
            else if (ANY((/ 'LAT     ','LATITUDE','Y       ' /) == upper(dim_name_tmp))) then
                y_dim_name = dim_name_tmp
                my = m_tmp
                y_dim_id = n
            else if (ANY((/'TIME', 'T   '/) == upper(dim_name_tmp))) then
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
    ! integer pure function storm_index(t, storm) result(index)
    integer function storm_index(t, storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        integer :: i

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

        stop "Storm location for data stroms is not supported."

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

        stop "Storm direction for data stroms is not supported."

    end function storm_direction

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
        real(kind=8), intent(inout) :: wind_tu(storm%mx,storm%my)
        real(kind=8), intent(inout) :: wind_tv(storm%mx, storm%my)
        real(kind=8), intent(inout) :: pressure_t(storm%mx,storm%my)

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
        
    end subroutine interp_array_data

    ! ==========================================================================
    ! find_nearest() finds nearest value to x and y for interpolation points
    ! Finds the nearest actual point to the patch values using minimum distance
    ! ==========================================================================
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

   ! ==========================================================================
   ! Check for netcdf file errors when loading data, only active if the
   ! NETCDF FFLAGS are in the Makefile
   ! ==========================================================================
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

    ! ==========================================================================
    ! seconds_from_epoch() Calculates seconds from 1970 from a datetime
    ! Returns the total seconds from the epoch includes leap years and days
    ! ==========================================================================
    pure function seconds_from_epoch(time) result(seconds)
        implicit none

        integer, intent(in) :: time(:) ! year, month, day, hour, minutes (optional)
        integer :: seconds
        integer, parameter, dimension(12) :: days_in_month=[31, 28, 31, 30, &
                                                            31, 30, 31, 31, &
                                                            30, 31, 30, 31]
        integer :: leap_days, days, minutes

        days = (time(1) - 1970) * 365
        leap_days = (time(1) - 1968)/4 - (time(1) - 1900)/4 + (time(1) - 1600)/400
        days = days + leap_days
        if (mod(time(1), 4) == 0 .and. (mod(time(1),100) /= 0 &
            .or. mod(time(1), 400) == 0).and.time(2) >2) THEN
                days = days + 1
        endif
        days = days + sum(days_in_month(1:time(2)-1)) + time(3)
        minutes = merge(0, time(size(time)), size(time) == 4)
        seconds = (days*86400) + (time(4)*3600) + (minutes*60)

    end function seconds_from_epoch
    

    ! ==========================================================================
    ! function upper() Returns the Upper case characater in response to input
    ! character
    ! ==========================================================================
    pure function upper(s1) result (s2)
        implicit none
        character(len=*), intent(in):: s1
        character(len=len(s1)) :: s2
        character :: ch
        integer, parameter :: duc = ichar('A') - ichar('a')
        integer :: i

        do i = 1, len(s1)
           ch = s1(i:i)
           if (ch >= 'a' .and. ch <= 'z') ch = char(ichar(ch)+duc)
           s2(i:i) = ch
        end do
    end function upper

end module data_storm_module
