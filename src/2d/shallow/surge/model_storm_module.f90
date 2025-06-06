! ==============================================================================
! model_storm_module
!
! Module contains routines for constructing a wind and pressure field based on
! the a parameterized model of the wind and pressure fields.
!
! ==============================================================================
!                   Copyright (C) Clawpack Developers 2017
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ==============================================================================
module model_storm_module

    implicit none
    save

    logical, private :: module_setup = .false.
    logical, private :: DEBUG = .false.

    ! Model storm type definition
    type model_storm_type
        ! Fore/hindcast size and current position
        integer :: num_casts

        ! Landfall - This is not used explicitly (t0 = landfall ideally)
        character(len=19) :: landfall

        ! These parameters are located at time points but are interpolated in
        ! time and space when the relevant fields are requested.

        ! Location of storm
        ! Track is a triplet with (time,longitude,latitude)
        real(kind=8), allocatable :: track(:,:)

        ! Storm parameterization
        real(kind=8), allocatable :: max_wind_radius(:)
        real(kind=8), allocatable :: max_wind_speed(:)
        real(kind=8), allocatable :: central_pressure(:)

        ! This is not always provided and is often determined by the last
        ! closed iso-bar of the storm
        real(kind=8), allocatable :: radius(:)

        ! Approximate velocity of storm, approximated via the track points
        ! using a first order difference on the sphere
        real(kind=8), allocatable :: velocity(:, :)

        ! Approximate pressure change of storm, approximated using first order differences
        real(kind=8), allocatable :: central_pressure_change(:)

    end type model_storm_type

    ! How to deterimine which way a storm should be made to spin
    ! The default defers simply to assuming y is a latitude
    ! Here either an integer or bool can be returned but as implemented 1
    ! refers to the Northern Hemisphere and therefore causes the storm to spin
    ! in a counter-clockwise direction.
    abstract interface
        logical pure function rotation_def(x, y)
            implicit none
            real(kind=8), intent(in) :: x, y
        end function rotation_def
    end interface
    procedure(rotation_def), pointer :: rotation

    ! Interal tracking variables for storm
    integer, private :: last_storm_index

    ! Atmospheric boundary layer, input variable in ADCIRC but always is
    ! set to the following value
    real(kind=8), parameter :: atmos_boundary_layer = 0.9d0

    ! Sampling adjustment from 1 min to 10 min winds
    real(kind=8), parameter :: sampling_time = 0.88d0

    ! Storm field ramping width - Represents crudely the ramping radial area
    real(kind=8), parameter :: RAMP_WIDTH = 100.0d3

    ! Time tracking tolerance allowance - allows for the beginning of the storm
    ! track to be close to but not equal the start time of the simulation
    real(kind=8), parameter :: TRACKING_TOLERANCE = 1d-10

contains


    ! Setup routine for model storms
    subroutine set_storm(storm_data_path, storm, storm_spec_type, log_unit)

        use geoclaw_module, only: deg2rad, spherical_distance, coordinate_system
        use amr_module, only: t0, rinfinity

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(model_storm_type), intent(inout) :: storm
        integer, intent(in) :: storm_spec_type, log_unit

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: i, k, io_status
        real(kind=8) :: x(2), y(2), ds, dt

        if (.not. module_setup) then

            ! Open data file
            print *,'Reading storm date file ', storm_data_path
            open(unit=data_file, file=storm_data_path, status='old',        &
                 action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening storm data file. status = ", io_status
                stop
            endif

            read(data_file, "(i4)") storm%num_casts
            read(data_file, "(a)") storm%landfall
            read(data_file, *)

            write(log_unit, "('Data length = ',i3)") storm%num_casts

            ! Allocate storm parameter file variables
            allocate(storm%track(3, storm%num_casts))
            allocate(storm%max_wind_speed(storm%num_casts))
            allocate(storm%max_wind_radius(storm%num_casts))
            allocate(storm%central_pressure(storm%num_casts))
            allocate(storm%radius(storm%num_casts))

            ! Now read in the storm data - note that the units are expected to
            ! be consistent with:
            ! max_wind_speed = m/s
            ! max_wind_radius = m
            ! central_pressure = Pa
            ! radius = m
            do i=1, storm%num_casts
                read(data_file, *) storm%track(:, i), &
                                   storm%max_wind_speed(i), &
                                   storm%max_wind_radius(i), &
                                   storm%central_pressure(i), &
                                   storm%radius(i)
            enddo
            close(data_file)

            ! Calculate storm speed
            allocate(storm%velocity(2, storm%num_casts))
            do i=1,storm%num_casts - 1
                ! Calculate velocity based on great circle distance between

                ! locations of storm
                x = storm%track(2:3,i)
                y = storm%track(2:3,i+1)

                dt = storm%track(1,i + 1) - storm%track(1,i)

                if (coordinate_system == 2) then
                    ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
                                            y(1), 0.5d0 * (x(2) + y(2)))
                    storm%velocity(1,i) = sign(ds / dt, y(1) - x(1))

                    ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
                                            0.5d0 * (x(1) + y(1)), y(2))
                    storm%velocity(2, i) = sign(ds / dt, y(2) - x(2))
                else
                    storm%velocity(1, i) = (y(1) - x(1)) / dt
                    storm%velocity(2, i) = (y(2) - x(2)) / dt
                end if
            end do

            ! Use last approximation for velocity point going forward
            storm%velocity(:, storm%num_casts) = storm%velocity(:,  &
                                                            storm%num_casts - 1)

            ! Calculate central pressure change
            allocate(storm%central_pressure_change(storm%num_casts))
            do i=1,storm%num_casts - 1
                ds = storm%central_pressure(i + 1) - storm%central_pressure(i)
                dt = storm%track(1,i + 1) - storm%track(1,i)
                storm%central_pressure_change(i) = ds / dt
            end do
            storm%central_pressure_change(storm%num_casts) = &
                storm%central_pressure_change(storm%num_casts - 1)

            if (t0 - storm%tracK(1, 1) < -TRACKING_TOLERANCE) then
                print *, "Start time", t0, " is outside of the tracking"
                print *, "tolerance range with the track start"
                print *, storm%track(1, 1), "."
                stop
            endif

            ! This is used to speed up searching for correct storm data
            last_storm_index = 2
            last_storm_index = storm_index(t0, storm)
            if (last_storm_index == -1) then
                print *,"Forecast not found for time ",t0,'.'
                stop
            endif

            ! Log everything to the surge log file
            write(log_unit,*) ""
            write(log_unit,*) "Storm Track and Strength"
            write(log_unit,*) ""
            do i=1, storm%num_casts
                write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
                                            (storm%velocity(k,i),k=1,2), &
                                             storm%max_wind_radius(i), &
                                             storm%max_wind_speed(i),  &
                                             storm%central_pressure(i)
            enddo

            module_setup = .true.
        end if

    end subroutine set_storm

    ! ==========================================================================
    !  storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function storm_location(t, storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Junk storage
        real(kind=8) :: junk(7)

        call get_storm_data(t, storm, location,                 &
                            junk(1:2), junk(3), junk(4), junk(5), junk(6), junk(7))

    end function storm_location

    ! ==========================================================================
    !  storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function storm_direction(t, storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: junk(7), velocity(2)

        ! Fetch velocity of storm which has direction encoded in it
        call get_storm_data(t, storm, junk(1:2), velocity, junk(3), junk(4),   &
                            junk(5), junk(6), junk(7))

        ! Unit directional vector
        theta = atan2(velocity(2),velocity(1))

    end function storm_direction

    ! ==========================================================================
    ! set_pressure
    !   Set pressure at a radius r
    ! ==========================================================================
    real(kind=8) pure function set_pressure(Pc, r, dp, mwr, B) result(pres)
        implicit none

        real(kind=8), intent(in) :: Pc, r, dp, mwr, B

        ! for any situation that could raise an underflow error, we set the
        ! second term to 0
        if ((mwr / r)**B > 100) then
            pres = Pc
        else
            pres = Pc + dp * exp(-(mwr / r)**B)
        endif

    end function set_pressure

    ! ==========================================================================
    !  storm_index(t,storm)
    !    Finds the index of the next storm data point
    !    This duplicates data_storm_module:storm_index
    ! ==========================================================================
    integer pure function storm_index(t, storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(model_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        ! Figure out where we are relative to the last time we checked for the
        ! index (stored in last_storm_index)

        ! Check if we are already beyond the end of the last forecast time
        if (last_storm_index == storm%num_casts + 1) then
            index = storm%num_casts + 1
        else
            t0 = storm%track(1,last_storm_index - 1)
            t1 = storm%track(1,last_storm_index)

            ! Check to see if we are close enough to the current index to just
            ! use that, tolerance is based on TRACKING_TOLERANCE
            if ((abs(t0 - t) < TRACKING_TOLERANCE) .or.   &
                (abs(t1 - t) < TRACKING_TOLERANCE) .or.   &
                (t0 < t .and. t < t1)) then

                index = last_storm_index
            else if ( t1 < t ) then
                found = .false.
                do index=last_storm_index+1,storm%num_casts
                    if (t < storm%track(1,index)) then
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
                        if (storm%track(1,index-1) < t) exit
                    enddo
                endif
            endif
        endif

    end function storm_index


    ! ==========================================================================
    !  get_storm_data()
    !    Interpolates in time and returns storm data.
    ! ==========================================================================
    subroutine get_storm_data(t, storm, location, velocity, max_wind_radius,  &
                              max_wind_speed, central_pressure,               &
                              radius, central_pressure_change)

        use geoclaw_module, only: deg2rad, latlon2xy, xy2latlon, coordinate_system

        implicit none

        ! Input
        real(kind=8), intent(in) :: t                 ! Current time
        type(model_storm_type), intent(in) :: storm   ! Storm

        ! Output
        real(kind=8), intent(out) :: location(2), velocity(2)
        real(kind=8), intent(out) :: max_wind_radius, max_wind_speed
        real(kind=8), intent(out) :: central_pressure, radius, central_pressure_change

        ! Local
        real(kind=8) :: fn(9), fnm(9), weight, tn, tnm, x(2)
        integer :: i

        ! Increment storm data index if needed and not at end of forecast
        i = storm_index(t,storm)
        last_storm_index = i

        ! List of possible error conditions
        if (i <= 1) then
            if (i == 0) then
                print *,"Invalid storm forecast requested for t = ",t
                print *,"Time requested is before any forecast data."
                print *,"    first time = ",storm%track(1,1)
                print *,"   ERROR = ",i
                stop
            else if (i > storm%num_casts + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif

        ! Interpolate in time for all parameters
        if (i == storm%num_casts + 1) then
            i = i - 1
            ! At last forecast, use last data for storm strength parameters and
            ! velocity, location uses last velocity and constant motion forward

            ! Convert coordinates temporarily to meters so that we can use
            ! the pre-calculated m/s velocities from before
            if (coordinate_system == 2) then
                x = latlon2xy(storm%track(2:3,i),storm%track(2:3,i))
                x = x + (t - storm%track(1,i)) * storm%velocity(:,i)
                fn = [xy2latlon(x,storm%track(2:3,i)), &
                      storm%velocity(:,i), storm%max_wind_radius(i), &
                      storm%max_wind_speed(i), storm%central_pressure(i), &
                      storm%radius(i), storm%central_pressure_change(i)]
            else
                x = x + (t - storm%track(1,i)) * storm%velocity(:,i)
                fn = [x, &
                      storm%velocity(:,i), storm%max_wind_radius(i), &
                      storm%max_wind_speed(i), storm%central_pressure(i), &
                      storm%radius(i), storm%central_pressure_change(i)]
            end if
        else
            ! Inbetween two forecast time points (the function storm_index
            ! ensures that we are not before the first data point, i.e. i > 1)
            tn = storm%track(1,i)
            tnm = storm%track(1,i-1)
            weight = (t - tnm) / (tn - tnm)
            fn = [storm%track(2:3,i),storm%velocity(:,i), &
                  storm%max_wind_radius(i),storm%max_wind_speed(i), &
                  storm%central_pressure(i), storm%radius(i), &
                  storm%central_pressure_change(i)]
            fnm = [storm%track(2:3,i - 1),storm%velocity(:,i - 1), &
                   storm%max_wind_radius(i - 1),storm%max_wind_speed(i - 1), &
                   storm%central_pressure(i - 1), storm%radius(i - 1), &
                   storm%central_pressure_change(i - 1)]
            fn = weight * (fn - fnm) + fnm
        endif

        ! Set output variables
        location = fn(1:2)
        velocity = fn(3:4)
        max_wind_radius = fn(5)
        max_wind_speed = fn(6)
        central_pressure = fn(7)
        radius = fn(8)
        central_pressure_change = fn(9)

    end subroutine get_storm_data


    ! ==========================================================================
    ! A set of common functions applied before calculating the circular component of
    ! wind speed using a particular wind model.

    ! 1. Adjust max wind speed by removing translational speed (bounding at 0).

    ! 2. If specified, convert surface winds to PBL winds
    ! ==========================================================================

    pure subroutine adjust_max_wind(tv, mws, mod_mws, convert_height)

        real (kind=8), intent(inout) :: tv(2)
        real (kind=8), intent(in) :: mws
        logical, intent(in) :: convert_height
        real (kind=8), intent(out) :: mod_mws

        real (kind=8) :: trans_speed, trans_mod

        trans_speed = sqrt(tv(1)**2 + tv(2)**2)
        mod_mws = mws - trans_speed

        ! Bound this at 0
        if (mod_mws < 0) then
            trans_mod = mws / trans_speed
            tv = tv * trans_mod
            mod_mws = 0
        end if

        ! Convert wind speed (10 m) to top of atmospheric boundary layer
        if (convert_height) mod_mws = mod_mws / atmos_boundary_layer

    end subroutine adjust_max_wind


    ! ==========================================================================
    ! A set of common functions applied after the circular component of wind speed
    ! has been calculated given a particular wind model:

    ! 1. Determine translation speed that should be added to final storm wind speed.
    !    This is tapered to zero as the storm wind tapers to zero toward the eye of the
    !    storm and at long distances from the storm. Avoid divide-by-0 error

    ! 2. If specified, convert wind velocity from top of atmospheric boundary layer to
    !    wind velocity at 10 m above the earth's surface.

    ! 3. Convert from 1 minute averaged winds to 10 minute averaged winds

    ! 4. Add back in translational velocity

    ! 5. Apply distance ramp to limit scope
    ! ==========================================================================
    pure subroutine post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
        wind_index, pressure_index, r, radius, tv, mod_mws, theta, convert_height, &
        n_hemisphere)

        use geoclaw_module, only: Pa => ambient_pressure

        integer, intent(in) :: maux, mbc, mx, my, i, j
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc), wind

        real (kind=8), intent(in) :: tv(2), mod_mws, theta, r, radius
        integer, intent(in) :: wind_index, pressure_index

        logical, intent(in) :: convert_height, n_hemisphere

        real (kind=8) :: trans_speed_x, trans_speed_y, ramp

        if (mod_mws > 0) then
            trans_speed_x = (abs(wind) / mod_mws) * tv(1)
            trans_speed_y = (abs(wind) / mod_mws) * tv(2)
        else
            trans_speed_x = 0
            trans_speed_y = 0
        end if

        if (convert_height) wind = wind * atmos_boundary_layer

        wind = wind * sampling_time

        ! Velocity components of storm (assumes perfect vortex shape)
        ! including addition of translation speed, and adjustment of orientation
        ! according to hemisphere (counter-clockwise on northern hemisphere).
        aux(wind_index,i,j)   = wind * merge(-1, 1, n_hemisphere) * sin(theta) + trans_speed_x
        aux(wind_index+1,i,j) = wind * merge(1, -1, n_hemisphere) * cos(theta) + trans_speed_y

        ! Apply distance ramp down(up) to fields to limit scope
        ramp = 0.5d0 * (1.d0 - tanh((r - radius) / RAMP_WIDTH))
        aux(pressure_index,i,j) = Pa + (aux(pressure_index,i,j) - Pa) &
                                * ramp
        aux(wind_index:wind_index+1,i,j) =                        &
                                aux(wind_index:wind_index+1,i,j)  &
                                * ramp

    end subroutine post_process_wind_estimate

    ! ==========================================================================
    ! Calculate storm centric polar coordinate location of grid cell center
    ! ==========================================================================
    pure subroutine calculate_polar_coordinate(x, y, sloc, r, theta)

        use geoclaw_module, only: deg2rad, coordinate_system
        use geoclaw_module, only: spherical_distance

        real(kind=8), intent(in) :: x, y, sloc(2)
        real(kind=8), intent(out) :: r, theta

        if (coordinate_system == 2) then
            ! lat-long coordinates, uses Haversine formula
            r = spherical_distance(x, y, sloc(1), sloc(2))
            theta = atan2((y - sloc(2)),(x - sloc(1)))
        else
            ! Strictly cartesian
            r = sqrt( (x - sloc(1))**2 + (y - sloc(2))**2)
            theta = atan2((y - sloc(2)), (x - sloc(1)))
        end if

    end subroutine calculate_polar_coordinate

    ! ==========================================================================
    ! Calculate central pressure difference
    ! ==========================================================================
    pure function get_pressure_diff(Pc) result(dp)
        use geoclaw_module, only: rho_air, ambient_pressure

        real (kind=8), intent(in) :: Pc
        real (kind=8) :: dp


        dp = ambient_pressure - Pc

        ! Limit central pressure deficit due to bad ambient pressure,
        ! really should have better ambient pressure...
        if (dp < 100.d0) dp = 100.d0

    end function get_pressure_diff

    ! ==========================================================================
    ! Calculate Holland parameters and limit the result
    ! ==========================================================================
    pure function get_holland_b(mod_mws, dp) result(B)
        use geoclaw_module, only: rho_air

        real (kind=8), intent(in) :: mod_mws, dp
        real (kind=8) :: B

        B = rho_air * exp(1.d0) * (mod_mws**2) / dp
        if (B <  1.d0) B = 1.d0
        if (B > 2.5d0) B = 2.5d0

    end function get_holland_b

    ! ==========================================================================
    ! Calculate B parameter according to Holland 2008 and limit the result
    ! ==========================================================================
    pure function get_holland_b_2008(dp, dpdt, slat, tv) result(B)
        real (kind=8), intent(in) :: dp, dpdt, slat, tv(2)
        real (kind=8) :: B, x, trans_speed

        trans_speed = sqrt(tv(1)**2 + tv(2)**2)
        x = 0.6d0 * (1.d0 - dp / 215.d0)
        B = -4.4d-5 * dp**2 + 0.01d0 * dp + 0.03d0 * dpdt &
            - 0.014d0 * abs(slat) + 0.15d0 * trans_speed**x + 1.d0
        if (B <  1.d0) B = 1.d0
        if (B > 2.5d0) B = 2.5d0

    end function get_holland_b_2008

    ! ==========================================================================
    !  Use the 1980 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_holland_1980_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: unsigned_coriolis

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: f, mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! Holland 1980 requires converting surface to gradient winds
        logical :: convert_height=.true.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! remove translational component
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            f = unsigned_coriolis(y)
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Set speed of wind at this point, handling case of grid cell centroid
                ! at eye of storm
                if ((mwr / r)**B > 100) then
                    wind = 0
                else
                    wind = sqrt((mwr / r)**B &
                        * exp(1.d0 - (mwr / r)**B) * mod_mws**2.d0 &
                        + (r * f)**2.d0 / 4.d0) - r * f / 2.d0
                endif

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                    wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                    convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_holland_1980_fields

    ! ==========================================================================
    !  Use the 2008 Holland model to set the storm fields
    ! ==========================================================================
    subroutine set_holland_2008_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: coriolis

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: f, mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! Holland 2008 uses surface winds exclusively
        logical :: convert_height=.false.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        dp = get_pressure_diff(Pc)
        B = get_holland_b_2008(dp, dPc_dt, sloc(2), tv)

        ! Use estimated mws instead of recorded mws
        mod_mws = sqrt(B / (rho_air * exp(1.d0)) * dp)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            f = coriolis(y)
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Set speed of wind at this point, handling case of grid cell centroid
                ! at eye of storm
                if ((mwr / r)**B > 100) then
                    wind = 0
                else
                    wind = sqrt((mwr / r)**B &
                        * exp(1.d0 - (mwr / r)**B) * mod_mws**2.d0 &
                        + (r * f)**2.d0 / 4.d0) - r * f / 2.d0
                endif

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                    wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                    convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_holland_2008_fields

    ! ==========================================================================
    !  Use the 2010 Holland model to set the storm fields
    !
    !  As in the original publication, this implementation computes
    !  cyclostrophic winds, i.e., it does not include the Coriolis force.
    ! ==========================================================================
    subroutine set_holland_2010_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)


        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        real(kind=8) :: rn, vn, xn, xx
        real(kind=8) :: dg, rg
        integer :: i,j

        ! Holland 2010 does not require converting surface to gradient winds
        logical :: convert_height=.false.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        dp = get_pressure_diff(Pc)
        B = get_holland_b_2008(dp, dPc_dt, sloc(2), tv)

        ! Use estimated mws instead of recorded mws
        mod_mws = sqrt(B / (rho_air * exp(1.d0)) * dp)

        ! Additional Holland parameters needed
        ! Since we don't have a second wind speed measurement, we assume that the
        ! wind speed is 10 m/s at a radius of 500 km. For comparison: The original
        ! publication assumed 17 m/s at a radius of 300 km.
        vn = 10
        rn = 500000
        dg = (mwr/rn)**B
        rg = dg * exp(1.d0-dg)
        xn = log(vn/mws)/log(rg)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Set speed of wind at this point, handling case of grid cell centroid
                ! at eye of storm
                if ((mwr / r)**B > 100) then
                    wind = 0
                else
                    ! # HOLLAND 2010 WIND SPEED CALCULATION
                    if (r <= mwr) then
                        xx = 0.5
                    else
                        xx = 0.5 + (r - mwr) * (xn - 0.5) / (rn - mwr)
                    endif
                    wind = mod_mws * ((mwr / r)**B * exp(1.d0 - (mwr / r)**B))**xx
                endif

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                    wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                    convert_height, rotation(x, y))
            enddo
        enddo

    end subroutine set_holland_2010_fields


    ! ==========================================================================
    !  Use the Chavas, Lin, and Emmanuel 2016 model
    ! ==========================================================================
    subroutine set_CLE_fields(maux, mbc, mx, my, xlower, ylower,    &
                              dx, dy, t, aux, wind_index,           &
                              pressure_index, storm)

        use geoclaw_module, only: rho_air, Pa => ambient_pressure
        use geoclaw_module, only: unsigned_coriolis, deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(in out) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2)
        real(kind=8) :: f, mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius
        real(kind=8) :: mod_mws, trans_speed, ramp, trans_mod
        integer :: i,j

        ! Variables for CLE model calculations
        integer, PARAMETER :: res=2000
        real(kind=8), save :: v_vec(1:res), p_vec(1:res)
        real(kind=8), save :: m_out_vec(1:res)
        real(kind=8), save :: last_time=-1
        real(kind=8), save :: dr, r_a, r_0
        integer :: io_status
        integer, save :: in_res, out_res


        ! Additional local storage for CLE Model
        real(kind=8) :: v, next_loc(2)
        integer :: r_index
        real(kind=8) :: r_factor, t_factor

        ! CLE does not require surface wind to be converted to gradient wind (?)
        logical :: convert_height=.false.

        ! Output the pressure and velocity fields.
        ! Determines the radial wind ONLY IF called at a new time.

        !First, interpolate all of the relevant tracked storm parameters.
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        ! (It appears the azimuthal winds are not used in the algorithm)
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        ! Fill array with CLE model parameters
        if (last_time /= t) then

            f = unsigned_coriolis(sloc(2))

            ! Determine the wind profile, but only if the profile isn't already
            ! saved. In v_vec and p_vec.
            last_time = t
            call solve_hurricane_wind_parameters(f,mwr,mws,res,r_0,r_a)
            dr = real(r_0) / (res - 1)
            !out_res = number of points in outer model/in_res = points in inner model
            out_res = ceiling(real((res - 1) * (r_0 - r_a)) / r_0)
            in_res = res - out_res
            !Determine angular momentum profile for outer model
            call integrate_m_out(f,r_0,r_0 - out_res*dr,out_res,m_out_vec)
            do j= 1, out_res
               !convert angular momentum to azimuthal velocity.
               r = r_0 - (j-1)*r_0/res
               m_out_vec(out_res + 1 - j) = (m_out_vec(out_res + 1 - j)- .5*f*r**2)/r
            end do

            p_vec(1) = 0
            v_vec(1) = 0

            do j=2, res
                ! Combine inner and outer model into single radial wind profile.
                if(j <= in_res) then
                   r = (j - 1)*dr
                   v_vec(j) = evaluate_v_a(f, mwr, mws, r)
                else
                   v_vec(j) = m_out_vec(j - in_res)
                end if
                ! Pressure is determined using the gradient wind equation.
                p_vec(j) = p_vec(j - 1) + (v_vec(j)**2)/r + f*v_vec(j)
            end do

            do j=1, res
                ! Normalize pressure to match measurements.
                p_vec(j) = (Pa - Pc)*p_vec(j)/p_vec(res) + Pc
            end do

            open(unit=4242, file='test.txt',status='replace', &
                 action='write',iostat=io_status)
            write(4242,*) v_vec
        end if

        ! Set fields
        do j=1-mbc,my+mbc
           ! Assign a clockwise, radially symmetric wind field by piece-wise linear
           ! interpolation of v_vec.  Likewise with p_vec.
            y = ylower + (j-0.5d0) * dy ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)
                r_index = floor(r/dr)
                r_factor = (r - r_index * dr) / dr
                r_index = r_index + 1

                if (r_index < res) then
                    aux(pressure_index, i, j) = (1 - r_factor) * p_vec(r_index) &
                                        + r_factor * p_vec(r_index + 1)
                    aux(wind_index, i, j) = -0.9 * ((1 - r_factor) * v_vec(r_index) &
                                        + r_factor * v_vec(r_index + 1)) * sin(theta)
                    aux(wind_index+1, i, j) = 0.9 * ((1 - r_factor) * v_vec(r_index) &
                                        + r_factor * v_vec(r_index + 1)) * cos(theta)
                else
                    aux(pressure_index,i,j) = Pa
                    aux(wind_index,i,j) = 0
                    aux(wind_index+1,i,j) = 0
                end if
            end do
        end do

    end subroutine set_CLE_fields

    ! ==========================================================================
    !  inner_derivative(f,r_m,v_m,r_a)
    !   Evaluate the derivative of the inner model with respect to the
    !   merge point radius. Fortran can not handle the computation of the
    !   derivative in one line, thus it is calculated by breaking it up.
    !   The analytical solution of the derivative is very messy, but possible
    !   to solve.
    ! ==========================================================================
    real(kind=8) function evaluate_inner_derivative(f,r_m,v_m,r_a,v_a) result(dMa)

        ! Input
        real(kind=8), intent(in) :: f, r_m, v_m, r_a, v_a

        ! Variables
        real(kind=8) :: denominator, M_a
        real(kind=8), parameter :: alpha=1.0

        ! Calculate necessary components of the derivative to make
        ! calculations easier
        denominator = 2 - alpha + alpha * (r_a / r_m)**2
        M_a = 0.5 * f * r_a**2 + r_a * v_a
        dMa = (2 * M_a / r_a) * (1/denominator)

    end function evaluate_inner_derivative

    ! ==========================================================================
    !  evaluate_v_a(f,r_m,v_m,r_a)
    !    Evaluates the inner velocity at the merge point that allows the inner
    !    and outer model to be continuous.
    ! ==========================================================================
    real(kind=8) function evaluate_v_a(f,r_m,v_m,r_a) result(v_a)

        ! Input
        real(kind=8), intent(in) :: f, r_m, v_m, r_a
        real(kind=8), parameter :: alpha=1.0

        v_a = ((2.0*(r_a/r_m)**2)/(2.0-alpha+alpha*(r_a/r_m)**2))**(1.0/(2.0-alpha))
        v_a = v_a*(0.5*f*r_m**2 + r_m*v_m)
        v_a = v_a - 0.5*f*r_a**2
        v_a = v_a/r_a

    end function evaluate_v_a


    ! ==========================================================================
    !  solve_r_0(f,r_a,v_a,res,r_0, r_guess, m_out)
    !    Determines the wind profile for a hurricane using
    !    the inner and outer model for the Chavas model.
    ! ==========================================================================
    subroutine solve_r_0(f, r_a, v_a, res, r_0, r_guess)

        ! Input
        real(kind=8), intent(in) :: f, r_a, v_a
        integer, intent(in) :: res

        ! Variables with less strict i/o
        real(kind=8) :: r_0, r_guess, m_diff, r_step, dr,r_m
        real(kind=8) :: m_out(1:res)
        integer :: i

        ! Parameters
        real(kind=8), PARAMETER :: R_TOL=0.1
        integer :: STEP_FLAG=1

        ! Intialize the guess
        r_0 = r_guess
        r_step = r_guess

        ! Golden section method to determine r_0, the radius
        ! at which the wind's speed drops to zero.
        do
            call integrate_m_out(f, r_0, r_a, res, m_out)
            if (m_out(1) - 0.5*f*r_a**2 - r_a*v_a < 0) then
                r_step = r_step/STEP_FLAG
                r_0 = r_0 + r_step
            else
                r_step = r_step/2.0
                r_0 = r_0 - r_step
                STEP_FLAG = 2
            end if

            ! Check the tolerance value
            if (r_step < R_TOL) then
                r_guess = r_0
                exit
            end if
        end do

    end subroutine solve_r_0

    ! ==========================================================================
    !  integrate_m_out(f, r_0, r_a, res, m_out)
    !   Integrates from r_0 to r_a to determine the outer wind profile.
    ! ==========================================================================
    subroutine integrate_m_out(f,r_0,r_a,res,m_out)

        ! Input variables
        real(kind=8), intent(in) :: f, r_0, r_a
        integer, intent(in) :: res
        real(kind=8), dimension(1:res) :: m_out

        ! Parameters and Other variables
        real(kind=8) :: V_m, dr, r_guess, r_p
        integer :: i
        real(kind=8), parameter :: chi=1.0

        ! Intialize
        m_out(res) = 0.5*f*r_0**2
        dr = (r_0 - r_a)/res
        r_p = r_0 - dr
        m_out(res-1) = 0.5*f*r_p**2 + r_p*f*(r_0-r_p)

        ! Integrates the outer model equation from r_0 to r_a
        do i = res-2, 1, -1
            r_guess = r_p
            r_p = r_guess - dr
            V_m = (m_out(i+1) - 0.5*f*r_p**2)/r_p
            m_out(i) = m_out(i+1) - chi*dr*((r_p*V_m)**2)/(r_0**2 - r_p**2)
        end do

    end subroutine integrate_m_out


    ! ==========================================================================
    !  solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)
    !   Determine r_0, r_a, and v_a.
    !   Determine the merge point when given the f value. At this current
    !   moment we are using values of chi = 1.0 and alpha = 1.0.
    ! ==========================================================================
    subroutine solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)

        ! Input
        real(kind=8), intent(in) :: f, r_m, v_m
        integer, intent(in) :: res

        ! Variables
        real(kind=8) :: r_0, r_0_guess, r_a, r_step, v_a
        real(kind=8) :: slope_difference
        real(kind=8) :: m_out(1:res)

        ! Parameters
        real(kind=8), parameter :: r_tol = 0.1
        integer :: step_flag = 1, i

        ! Parameters
        real(kind=8) :: dr, r
        real(kind=8) :: inner_res, outer_res
        real(kind=8), dimension(1:res) :: v_out, v_in
        real(kind=8), parameter :: chi=1.0
        real(kind=8), parameter :: alpha=1.0

        ! Initialize guesses for merge and r_0
        r_a = 2.0*r_m
        r_0_guess = 5.0*r_m
        r_step = r_a

        do
            v_a = evaluate_v_a(f,r_m,v_m,r_a)
            !print *, "v_a: ", v_a
            if (v_a < 0) then
                r_step = r_step/2
                r_a = r_a - r_step
                cycle
            end if

            call solve_r_0(f,r_a,v_a,res,r_0,r_0_guess)

            slope_difference = evaluate_inner_derivative(f,r_m,v_m,r_a,v_a) - &
                                chi*((r_a*v_a)**2)/(r_0**2 - r_a**2)
            if (slope_difference < 0) then
                r_step = r_step/2.0
                r_a = r_a - r_step
                step_flag = 1
            else
                r_step = r_step/step_flag
                r_a = r_a + r_step
            end if

            if (r_step < r_tol) then
                exit
            end if
        end do

    end subroutine solve_hurricane_wind_parameters

    ! ==========================================================================
    !  solve_wind_profile(f, r_m, v_m, res)
    ! ==========================================================================
    subroutine solve_wind_profile(f, r_m, v_m, res)

        ! Input
        real(kind=8), intent(in) :: f, r_m, v_m
        integer, intent(in) :: res
        real(kind=8), save :: r_0, r_a, v_a

        ! Parameters
        real(kind=8) :: dr, r
        real(kind=8) :: inner_res, outer_res
        real(kind=8), dimension(1:res) :: v_out, v_in, m_out
        integer :: i

        call solve_hurricane_wind_parameters(f, r_m, v_m, res, r_0, r_a)

        ! Determine the inner model velocity vector
        dr = (r_a/res)
        r = dr
        v_in(1) = 0


        do i = 2, res, 1
            v_in(i) = evaluate_v_a(f, r_m, v_m, r)
            r = r + dr
        end do

        ! Determine the outer model velocity vector
        dr = (r_0 - r_a)/res
        r = r_0 - dr
        v_out(res) = 0

        call integrate_m_out(f, r_0, r_a, res, m_out)

        do i = res-1, 1, -1
            v_out = (m_out(i) - 0.5*f*r**2)/r
            r = r - dr
        end do

    end subroutine solve_wind_profile

    ! ==========================================================================
    !  Use the SLOSH model to set the storm fields
    !  based on info from: https://www.mdpi.com/2073-4433/10/4/193/pdf
    !  and https://www.hsdl.org/?abstract&did=779560
    ! ==========================================================================
    subroutine set_SLOSH_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius
        real(kind=8) :: mod_mws, trans_speed_x, trans_speed_y
        integer :: i,j

        ! SLOSH model requires converting surface to gradient winds (for calculating
        ! Holland B)
        logical :: convert_height=.true.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! slosh winds use a 10-min wind speed so we apply correction factor upstream
        mod_mws = mod_mws * sampling_time

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Speed of wind at this point
                ! Note that in reality SLOSH does not directly input mws but instead
                ! estimates it from central pressure and mwr using a lookup table and
                ! iterative procedure (see p. 14 of Jelesnianski et al. 1992)
                ! It's assumed to be equal to 10-min winds, so we apply the correction
                ! upstream
                wind = (2.d0 * mws * mwr * r) / (mwr**2.d0 + r**2.d0)

                !! SLOSH has custom addition of translational velocity
                ! so we do not use `post_process_wind_estimate` function
                trans_speed_x = tv(1) * mwr * r / (mwr**2.d0 + r**2.d0)
                trans_speed_y = tv(2) * mwr * r / (mwr**2.d0 + r**2.d0)
                
                aux(wind_index,i,j)   = wind * merge(-1, 1, rotation(x, y)) * sin(theta) + trans_speed_x
                aux(wind_index+1,i,j) = wind * merge(1, -1, rotation(x, y)) * cos(theta) + trans_speed_y
            enddo
        enddo

    end subroutine set_SLOSH_fields


    ! ==========================================================================
    !  Use the Rankine model to set the wind field and H80 for pressure
    ! ==========================================================================
    subroutine set_rankine_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! Rankine requires converting surface to gradient winds
        logical :: convert_height=.true.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Set speed of wind at this point, handling case of grid cell centroid
                ! at eye of storm
                if ((mwr / r)**B > 100) then
                    wind = 0
                else
                    ! Speed of wind at this point
                    if (r < mwr) then
                        wind = mws * (r / mwr)
                    else
                        wind = mws * (mwr / r)
                    endif
                endif

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_rankine_fields


    ! ==========================================================================
    !  Use the Modified Rankine model to set the wind field and H80 for pressure
    ! ==========================================================================
    subroutine set_modified_rankine_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! Rankine requires converting surface to gradient winds
        logical :: convert_height=.true.

        ! additional decay parameter
        real(kind=8) :: alpha

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! set modified rankine decay param (from Mallen, Montgomery and Wang 2005)
        if ((mws) < 30.d0) then
            alpha = 0.31d0
        else if (mws > 50.d0) then
            alpha = 0.48d0
        else
            alpha = 0.35d0
        end if

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Set speed of wind at this point, handling case of grid cell centroid
                ! at eye of storm
                if ((mwr / r)**B > 100) then
                    wind = 0
                else
                    ! Speed of wind at this point
                    if (r < mwr) then
                        wind = mws * (r / mwr)
                    else
                        wind = mws * (mwr / r)**alpha
                    endif
                endif

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_modified_rankine_fields


    ! ==========================================================================
    !  Use the DeMaria model to set the wind fields, using Holland 1980 for
    !  pressures.
    !  original paper: deMaria 1987, Monthly Weather Review
    ! ==========================================================================
    subroutine set_deMaria_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! deMaria requires converting surface to gradient winds (just for Holland B)
        logical :: convert_height=.true.

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Speed of wind at this point
                ! Commented out version has 2 more free parameters that we would need
                ! to set (chosen to be equal to 1 for the time being)
                !wind = (r * mws / mwr) &
                !        * exp((1.d0 - (r / mwr)**b) * (1.d0/d))
                wind = mod_mws * (r / mwr) &
                        * exp(1.d0 - (r / mwr))

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                    wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                    convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_deMaria_fields


    ! ==========================================================================
    !  Use the Willoughby 2006 model to set the wind field, using H80 for pressure
    !  Original paper https://doi.org/10.1175/MWR3106.1
    ! ==========================================================================
    subroutine set_willoughby_fields(maux, mbc, mx, my, xlower, ylower,    &
                                       dx, dy, t, aux, wind_index,           &
                                       pressure_index, storm)

        use geoclaw_module, only: rho_air
        use geoclaw_module, only: deg2rad
        use geoclaw_module, only: spherical_distance

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(model_storm_type), intent(inout) :: storm

        ! Array storing wind and presure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2), B
        real(kind=8) :: mwr, mws, Pc, dPc_dt, dp, wind, tv(2), radius, mod_mws
        integer :: i,j

        ! Willoughby requires converting surface to gradient winds
        logical :: convert_height=.true.

        ! other Willoughby params
        real(kind=8) :: n, X1, X2, A, xi, W

        ! Get interpolated storm data
        call get_storm_data(t, storm, sloc, tv, mwr, mws, Pc, radius, dPc_dt)

        ! pre-process max wind speed and translational velocity
        call adjust_max_wind(tv, mws, mod_mws, convert_height)

        dp = get_pressure_diff(Pc)
        B = get_holland_b(mod_mws, dp)

        ! Get willoughby coeffs
        ! (See Section 4)
        n = 2.134d0 + 0.0077d0*mod_mws - 0.4522d0*log(mwr/1.d3) - 0.0038*abs(sloc(2))
        X1 = (317.1d0 - 2.026d0*mod_mws + 1.915d0*abs(sloc(2))) * 1.d3
        X2 = 25.d3
        A = max(0.5913d0 + 0.0029d0*mod_mws - 0.1361*log(mwr/1.d3) - 0.0042d0*abs(sloc(2)), &
                0.d0)

        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                call calculate_polar_coordinate(x, y, sloc, r, theta)

                ! Set pressure field
                aux(pressure_index,i,j) = set_pressure(Pc, r, dp, mwr, B)

                ! Speed of wind at this point
                ! (.9 and 1.1 are chosen as R1 and R2 rather than fit )
                ! Note this may mess up the location of Vmax slightly since it does not
                ! have to lie at xi=0.5 (see Section 2 of paper)
                if (r < 0.9d0 * mwr) then
                    wind = mod_mws * (r / mwr)**n
                else if (r > 1.1d0 * mwr) then
                    wind = mod_mws * ( (1-A) * exp(-(r-mwr)/X1) + A * exp(-(r-mwr)/X2) )
                else
                    xi = (r - 0.9d0*mwr)/(0.2d0*mwr)
                    W = 126.d0*xi**5 - 420.d0*xi**6 + 540.d0*xi**7 - 315.d0*xi**8 + &
                        70.d0*xi**9
                    wind = mod_mws * (r / mwr)**n * (1-W) + &
                            mod_mws * &
                                ( (1-A) * exp(-(r-mwr)/X1) + A * exp(-(r-mwr)/X2) ) * W
                end if

                call post_process_wind_estimate(maux, mbc, mx, my, i, j, wind, aux, &
                    wind_index, pressure_index, r, radius, tv, mod_mws, theta, &
                    convert_height, rotation(x, y))

            enddo
        enddo

    end subroutine set_willoughby_fields
end module model_storm_module
