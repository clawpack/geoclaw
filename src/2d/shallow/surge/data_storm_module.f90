! ==============================================================================
! model_storm_module 
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

        ! Storm data, wind velocity in x and y, pressure and wind speed
        real(kind=8), allocatable :: pressure(:,:,:)
        real(kind=8), allocatable :: wind_speed(:,:,:)

        ! Wind field latitude/longitude arrays
        real(kind=8), allocatable :: latitude(:)
        real(kind=8), allocatable :: longitude(:)

        ! Time steps from wind/pressure files in seconds
        integer, allocatable :: time(:)

    end type data_storm_type

contains

    ! Setup routine for data type storms
    subroutine set_storm(storm_data_path, storm, model_type, log_unit)
        use netcdf
        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(data_storm_type), intent(inout) :: storm
        integer, intent(in) :: model_type, log_unit

        ! Set up for dimension info, dims are lat/lon/time
        integer :: num_dims, x_dim_id, y_dim_id, t_dim_id, dim_ids(3)
        character (len=10) :: x_dim_name, y_dim_name, t_dim_name
        ! Length of dims, mx, my, mt
        integer :: mx, my, mt

        ! Set up for variable info, variables are speed, pressure, lat, lon, time
        integer :: nvars, var_type, var_id
        character (len=10) :: var_name

        ! Local storage
        integer :: nc_fid, n

        if (.not. module_setup) then
            ! Open file and get file ID
            call check_netcdf_error(nf90_open(storm_data_path, nf90_nowrite, nc_fid))
            call check_netcdf_error(nf90_inquire(nc_fid, num_dims, nvars))

            ! Get the dimension names, sizes from the file
            call get_dim_info(nc_fid, num_dims, x_dim_id, x_dim_name, mx, &
            y_dim_id, y_dim_name, my, t_dim_id, t_dim_name, mt)

            ! allocate arrays in storm object
            allocate(storm%wind_speed(mx, my, mt))
            allocate(storm%pressure(mx, my, mt))
            allocate(storm%longitude(mx))
            allocate(storm%latitude(my))
            allocate(storm%time(mt))

            ! Fill out variable data/info
            do n = 1, nvars
                call check_netcdf_error(nf90_inquire_variable(nc_fid, n, var_name, var_type, num_dims, dim_ids))
                if (ANY((/'SPEED     ', 'WIND SPEED'/) == Upper(var_name))) then
                    call check_netcdf_error(nf90_inq_varid(nc_fid, var_name, var_id))
                    call check_netcdf_error(nf90_get_var(nc_fid, var_id, storm%wind_speed))
                elseif ('PRESSURE' == Upper(var_name)) then
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
                end if
            end do


        end if

        stop "Data-derived storm are not yet implemented!"

        if (.not. module_setup) then

            module_setup = .true.
        end if

    end subroutine set_storm


    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(data_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        stop "Data-derived storm are not yet implemented!"

    end function storm_location

    ! ==========================================================================
    !  storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function storm_direction(t, storm) result(theta)

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
        type(data_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        stop "HWRF data input is not yet implemented!"

    end subroutine set_HWRF_fields

    subroutine check_netcdf_error(ios)

        use netcdf

        implicit none

        integer, intent(in) :: ios

        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if

    end subroutine check_netcdf_error

    subroutine get_dim_info(nc_file, ndims, x_dim_id, x_dim_name, mx, &
        y_dim_id, y_dim_name, my, t_dim_id, t_dim_name, mt)
        use netcdf
        implicit none
        integer, intent(in) :: nc_file, ndims
        integer, intent(out) :: x_dim_id, y_dim_id, mx, my, t_dim_id, mt
        character (len = *), intent(out) :: x_dim_name, y_dim_name, t_dim_name
        integer :: m_tmp, n
        character(20) :: dim_name_tmp

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
    end subroutine get_dim_info


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



