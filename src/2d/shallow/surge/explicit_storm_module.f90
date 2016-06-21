! ==============================================================================
! explicit_storm_module 
!
! Module contains routines for returning wind and pressure field based on
! specified values in data files.
!
! The data files are expected to be of the form:
! *lat.dat
! *lon.dat
! *u10.dat
! *v10.dat
! *pmsl.dat
! where * is specified in setrun.py or empty by default.
! 
! ==============================================================================
module explicit_storm_module

    implicit none
    save

    ! Explicit storm type definition
    ! Specified wind & pressure field 
    type explicit_storm_type
        ! Size of spatial grids,
        !  corresponds to lengths of lat/lon arrays
        integer :: num_rows
        integer :: num_cols

        ! Location of storm field values
        ! longitude and latitude arrays
        ! start from SW corner
        real(kind=8), allocatable :: lat(:)
        real(kind=8), allocatable :: lon(:)

        ! We keep two time snapshots in memory 
        !  for interpolation, with t_next > t > t_prev
        real(kind=8) :: t_next, t_prev
        real(kind=8), allocatable :: u_next(:,:)
        real(kind=8), allocatable :: v_next(:,:)
        real(kind=8), allocatable :: p_next(:,:)
        real(kind=8), allocatable :: u_prev(:,:)
        real(kind=8), allocatable :: v_prev(:,:)
        real(kind=8), allocatable :: p_prev(:,:)
        ! The values will be updated incrementally:
        !  when t>t_next, replace older states with newer states
        !  and read in next snapshot.
        
        ! These will be used during the interpolation step
        !  for times between t_prev and t_next
        real(kind=8) :: t
        real(kind=8), allocatable :: u(:,:)
        real(kind=8), allocatable :: v(:,:)
        real(kind=8), allocatable :: p(:,:)

        ! Keep track of how many snapshots have been read
        integer :: last_storm_index

        ! for compatibility with parent storm module (?)
        ! Storm physics
        real(kind=8) :: ambient_pressure = 101.3d3 ! 101300 Pascals
        real(kind=8) :: rho_air = 1.3d0

        ! Store the storm data file for repeated reading
        character(len=90) :: data_path_root

    end type explicit_storm_type

    ! Internal tracking variables for storm
    logical, private :: DEBUG = .true. 
    !logical, private :: DEBUG = .false. 

contains

    ! Setup routine for the explicit model
    ! Open data files, get parameters, allocate memory,
    !  and read in first two time snapshots of storm data
    subroutine set_explicit_storm(storm_data_path_root, storm, log_unit)

        use geoclaw_module, only: coordinate_system
        use amr_module, only: t0

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path_root
        type(explicit_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: l_file = 701
        integer :: i, j, k, io_status, num_rows, num_cols
        real(kind=8) :: forecast_time,last_time,x(2),y(2),ds,dt,dx,dy
        character(len=100) :: storm_data_path

        ! Reading buffer variables
        character(len=20000) :: dummy_line

        ! Storm type only works on lat-long coordinate systems
        if (coordinate_system /= 2) then
            stop "explicit storm type does only works on lat-long coordinates."
        endif

        ! We need to count two things:
        !   number of latitude coords (ny)
        !   number of longitude coords (nx)
        ! We don't need number of time snapshots (nt)
        !   since we will load snapshots only as needed
        ! The datafiles for lat & lon contain nx*ny values
        ! The datafiles for u, v, p contain nx*ny*nt values

        ! Open latitudes data file
        if (present(storm_data_path_root)) then
            storm%data_path_root = storm_data_path_root
        else
            storm%data_path_root = './'
        endif
        storm_data_path = trim(storm%data_path_root) // "lat.dat"
        print *,'Reading latitudes data file ',storm_data_path
        open(unit=l_file,file=storm_data_path,status='old', &
                action='read',iostat=io_status)
        if (io_status /= 0) then
            print "(a,i2)", "Error opening latitudes data file. status = ", io_status
            stop 
        endif            

        ! Count number of data columns
        ! by reading the first line and 
        ! counting the decimal points
        read( l_file, '( a )',iostat=io_status) dummy_line
        num_cols = count( [ ( dummy_line( i:i ), i = 1, len( dummy_line ) ) ] == '.' )
        storm%num_cols = num_cols

        ! Count number of data rows
        num_rows = 1
        do
            read (l_file, *, iostat=io_status)
            ! Exit loop if we ran into an error or we reached the end of the file
            if (io_status /= 0) then
                exit
            else
                num_rows = num_rows + 1
            endif
        end do
        rewind(l_file) ! closed later
        storm%num_rows = num_rows

        ! save lat & lon coords, times, and u/v/p fields
        allocate(storm%lat(num_rows))
        allocate(storm%lon(num_cols))
        allocate(storm%u_prev(num_cols,num_rows))
        allocate(storm%v_prev(num_cols,num_rows))
        allocate(storm%p_prev(num_cols,num_rows))
        allocate(storm%u_next(num_cols,num_rows))
        allocate(storm%v_next(num_cols,num_rows))
        allocate(storm%p_next(num_cols,num_rows))
        allocate(storm%u(num_cols,num_rows))
        allocate(storm%v(num_cols,num_rows))
        allocate(storm%p(num_cols,num_rows))

        ! Read in latitude coords
        do i = 1, num_rows
            read (l_file, *, iostat=io_status) storm%lat(i)
            if (io_status /= 0) exit
        end do
        close(l_file)

        ! Open longitudes data file
        storm_data_path = trim(storm%data_path_root) // "lon.dat"
        print *,'Reading longitudes data file ',storm_data_path
        open(unit=l_file,file=storm_data_path,status='old', &
                action='read',iostat=io_status)
        if (io_status /= 0) then
            print "(a,i2)", "Error opening longitudes data file. status = ", io_status
            stop 
        endif            
        ! Read in longitude coords. Just need first line.
        read (l_file, *, iostat=io_status) storm%lon
        close(l_file)

        ! This is used to speed up searching for correct storm data
        storm%last_storm_index = 0

        ! Read in the first two storm data snapshots
        call read_explicit_storm(storm,t0)
        call read_explicit_storm(storm,t0)
        ! storm%last_storm_index will be 2

        if (t0 < storm%t_prev) then
            print *,t0,storm%t_prev
            stop "Simulation start time preceeds storm data."
        endif

    end subroutine set_explicit_storm


    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function date_to_seconds(year,months,days,hours,minutes, &
                                               seconds) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: year, months, days, hours, minutes
        real(kind=8), intent(in) :: seconds

        ! Local storage
        integer :: total_days

        ! Count number of days
        total_days = days

        ! Add days for months that have already passed
        if (months > 1) total_days = total_days + 31
        if (months > 2) then
            if (int(year / 4) * 4 == year) then
                total_days = total_days + 29
            else
                total_days = total_days + 28
            endif
        endif
        if (months > 3)  total_days = total_days + 30
        if (months > 4)  total_days = total_days + 31
        if (months > 5)  total_days = total_days + 30
        if (months > 6)  total_days = total_days + 31
        if (months > 7)  total_days = total_days + 30
        if (months > 8)  total_days = total_days + 31
        if (months > 9)  total_days = total_days + 30
        if (months > 10) total_days = total_days + 31
        if (months > 11) total_days = total_days + 30

        ! Convert everything to seconds since the beginning of the year
        time = real((total_days - 1) * 86400 + hours * 3600 + minutes * 60,kind=8)
        time = time + seconds

    end function date_to_seconds

    ! ==========================================================================
    !  real(kind=8) pure ymdh_to_seconds(ymdh)
    !    Convert time from year, month, day, hour (YEARMODAHR)
    !    to seconds since the beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function ymdh_to_seconds(ymdh) result(time)
      
        implicit none

        ! Input
        integer, intent(in) :: ymdh

        ! Local storage
        integer :: year, month, day, hour

        ! collect datetime (input format is YEARMODAHR)
        hour = mod(ymdh,100)
        day = mod(ymdh/100,100)
        month = mod(ymdh/10000,100)
        year = ymdh / 1000000

        ! process as seconds from start of year
        time = date_to_seconds(year,month,day,hour,0,0.d0)

    end function ymdh_to_seconds

    ! ==========================================================================
    !  read_explicit_storm_data_file()
    !    Opens storm data file and reads next storm entry
    !    Currently only for ASCII file
    ! ==========================================================================
    subroutine read_explicit_storm_file(data_path,storm_array,num_rows,last_storm_index,timestamp)

        implicit none

        ! Subroutine I/O
        real(kind=8), intent(in out) :: storm_array(:,:)
        !real(kind=8), intent(in) :: t
        character(len=*), intent(in) :: data_path
        integer, intent(in) :: num_rows, last_storm_index
        integer, intent(inout) :: timestamp

        ! Local storage
        integer :: j, k, datafile, iostatus

        ! Open the input file
        !
        open(unit=datafile,file=data_path,status='old', &
                action='read',iostat=iostatus)
        if (iostatus /= 0) then
            print "(a,i2)", "Error opening data file: ",data_path
            print "(a,i2)", "Status = ", iostatus
            stop 
        endif            
        ! Advance to the next time step to be read in
        ! Skip entries based on total number previously read
        do k = 1, last_storm_index
            do j = 1, num_rows
                read(datafile, *)
                ! Exit loop if we ran into an error or we reached the end of the file
                if (iostatus /= 0) then
                    print "(a,i2)", "Unexpected end-of-file reading ",data_path
                    print "(a,i2)", "Status = ", iostatus
                    exit
                endif
            enddo
        enddo
        ! Read in next time snapshot 
        do j = 1, num_rows
            read(datafile, *) timestamp, storm_array(:,j) 
            ! Exit loop if we ran into an error or we reached the end of the file
            if (iostatus /= 0) then
                print "(a,i2)", "Unexpected end-of-file reading ",data_path
                print "(a,i2)", "Status = ", iostatus
                exit
            endif
        enddo
        close(datafile) 

        timestamp = ymdh_to_seconds(timestamp)

    end subroutine read_explicit_storm_file

    ! ==========================================================================
    !  read_explicit_storm_data()
    !    Reads storm fields for next time snapshot
    !    Currently only for ASCII files
    ! ==========================================================================

    subroutine read_explicit_storm(storm,t)

        implicit none

        ! Subroutine I/O
        type(explicit_storm_type), intent(in out) :: storm
        real(kind=8), intent(in) :: t

        ! Local storage
        integer :: j, k
        character(len=100) :: data_path

        ! Reading buffer variables
        integer :: timestamp

        ! Overwrite older storm states with newer storm states
        storm%t_prev = storm%t_next
        storm%u_prev = storm%u_next 
        storm%v_prev = storm%v_next 
        storm%p_prev = storm%p_next 
        
        ! Current time t currently unused in favor of storm%last_storm_index.
        ! This should probably be changed in the future.

        ! Read the u-velocity file
        !
        data_path = trim(storm%data_path_root) // "u10.dat"
        call read_explicit_storm_file(data_path,storm%u_next,storm%num_rows,storm%last_storm_index,timestamp)

        ! Read v-velocity file
        !
        data_path = trim(storm%data_path_root) // "v10.dat"
        call read_explicit_storm_file(data_path,storm%v_next,storm%num_rows,storm%last_storm_index,timestamp)

        ! Read pressure file
        !
        data_path = trim(storm%data_path_root) // "pmsl.dat"
        call read_explicit_storm_file(data_path,storm%p_next,storm%num_rows,storm%last_storm_index,timestamp)

        ! Convert pressure units: mbar to Pa
        storm%p_next = storm%p_next * 1.0e2

        ! Save timestamp (sec) of next snapshot
        storm%t_next = timestamp

        ! Update number of storm snapshots read in
        storm%last_storm_index = storm%last_storm_index + 1
        if (DEBUG) print *, "last_storm_index=", storm%last_storm_index


    end subroutine read_explicit_storm

    ! ==========================================================================
    !  integer pure get_lat_index(lat)
    !    Returns index of latitude array of the storm data
    !    corresponding to input lat.
    ! ==========================================================================
    pure integer function get_lat_index(lat,storm) result(i)
      
        implicit none

        ! Input
        real(kind=8), intent(in) :: lat
        type(explicit_storm_type), intent(in) :: storm

        ! Local storage
        real(kind=8) :: dy

        ! Out-of-bound conditions:
        !  no error messages generated to optimize performace
        !  (pure function specification)
        if (lat < storm%lat(1)) then
            i = 1
            !print *, "Warning: latitude ", lat, " out of bounds. Using ", &
                !storm%lat(1), " instead."
        else if (lat > storm%lat(storm%num_rows)) then
            i = storm%num_rows
            !print *, "Warning: latitude ", lat, " out of bounds. Using ", &
                !storm%lat(storm%num_rows), " instead."
        else
            ! Find spacing between latitude values
            dy = (storm%lat(storm%num_rows) - storm%lat(1)) / storm%num_rows
            ! Determine index based on spacing
            i = 1 + (lat - storm%lat(1)) / dy
        endif

    end function get_lat_index

    ! ==========================================================================
    !  integer pure get_lon_index(lon)
    !    Returns index of lontitude array of the storm data
    !    corresponding to input lon.
    ! ==========================================================================
    pure integer function get_lon_index(lon,storm) result(i)
      
        implicit none

        ! Input
        real(kind=8), intent(in) :: lon
        type(explicit_storm_type), intent(in) :: storm

        ! Local storage
        real(kind=8) :: dx

        ! Out-of-bound conditions:
        !  no error messages generated to optimize performace
        !  (pure function specification)
        if (lon < storm%lon(1)) then
            i = 1
            !print *, "Warning: Longitude ", lon, " out of bounds. Using ", &
                !storm%lon(1), " instead."
        else if (lon > storm%lon(storm%num_cols)) then
            i = storm%num_cols
            !print *, "Warning: Longitude ", lon, " out of bounds. Using ", &
                !storm%lon(storm%num_cols), " instead."
        else
            ! Find spacing between longitude values
            dx = (storm%lon(storm%num_cols) - storm%lon(1)) / storm%num_cols
            ! Determine index based on spacing
            i = 1 + (lon - storm%lon(1)) / dx
        endif

    end function get_lon_index

    ! ==========================================================================
    !  set_explicit_storm_fields()
    ! ==========================================================================
    subroutine set_explicit_storm_fields(maux,mbc,mx,my,xlower, &
                                    ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, storm)

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(explicit_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y
        integer :: i,j,k,l

        if (t < storm%t_prev) then
            print *, "Simulation time precedes storm data in memory."
            print *, "t=",t,"< t_prev=",storm%t_prev,"t_next=",storm%t_next
        endif

        if (t > storm%t_next) then
            ! Load two snapshots into memory, at times t_next and t_prev
            !$OMP CRITICAL (READ_STORM)
            do while (t > storm%t_next)
            ! update all storm data, including value of t_next
                if (DEBUG) print *,"loading new storm snapshot ","t=",t,"t_next=",storm%t_next
                call read_explicit_storm(storm,t)
                ! TODO: simulation should keep running even if storm data ends
            enddo
            !$OMP END CRITICAL (READ_STORM)
        endif
        
        ! Interpolate storm data in time
        ! t_prev <= t <= t_next
        if (t > storm%t) then
            ! Bring interpolated storm data up to speed
            !$OMP CRITICAL (INTERP_STORM)
            if (t > storm%t) then
                storm%t = t
                storm%u = (storm%u_prev*(storm%t_next-t) + storm%u_next*(t-storm%t_prev)) / &
                        (storm%t_next-storm%t_prev)
                storm%v = (storm%v_prev*(storm%t_next-t) + storm%v_next*(t-storm%t_prev)) / &
                        (storm%t_next-storm%t_prev)
                storm%p = (storm%p_prev*(storm%t_next-t) + storm%p_next*(t-storm%t_prev)) / &
                        (storm%t_next-storm%t_prev)
            endif
            !$OMP END CRITICAL (INTERP_STORM)
        endif

        ! Set fields
        ! Determine lat/long of each cell in layer,
        !  determine corresponding storm cell indices
        !  (unfortunately the lat/lons aren't evenly spaced),
        !  then get value of corresponding storm data cell.
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            k = get_lat_index(y,storm) ! storm index of latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude
                l = get_lon_index(x,storm) ! storm index of longitude
                ! Set pressure field
                aux(pressure_index,i,j) = storm%p(l,k)
                ! Set velocity components of storm 
                aux(wind_index,i,j)   = storm%u(l,k)
                aux(wind_index+1,i,j) = storm%v(l,k)
            enddo
        enddo

    end subroutine set_explicit_storm_fields

end module explicit_storm_module

