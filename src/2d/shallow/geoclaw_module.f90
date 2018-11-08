! ============================================================================
!  File:        geoclaw_mod
! ============================================================================
!    Copyright (C) 2010-04-21 Clawpack Developers http://www.clawpack.org
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module geoclaw_module

    implicit none
    save

    logical, private :: module_setup = .false.
    
    ! ========================================================================
    !  Constants
    ! ========================================================================
    integer, parameter :: GEO_PARM_UNIT = 78
    integer, parameter :: KAPPA_UNIT = 42
    real(CLAW_REAL), parameter :: pi = 4.d0*datan(1.d0)
    real(CLAW_REAL), parameter :: DEG2RAD = pi / 180.d0
    real(CLAW_REAL), parameter :: RAD2DEG = 180.d0 / pi
    
    ! ========================================================================
    !  Physics
    ! ========================================================================
    real(CLAW_REAL) :: grav, rho_air, ambient_pressure, earth_radius, sea_level
    ! Water density can be an array to handle multiple layers
    real(CLAW_REAL), allocatable :: rho(:)
    integer :: coordinate_system

    ! Rotational velocity of Earth
    real(CLAW_REAL), parameter :: omega = 2.0d0 * pi / 86164.2d0
    
    ! Forcing
    logical :: coriolis_forcing ! True then coriolis terms included in src
    real(CLAW_REAL) :: theta_0 ! Used if using the beta-plane approximation
    logical :: friction_forcing ! Friction forcing will be applied
    real(CLAW_REAL), dimension(:),allocatable :: manning_coefficient, manning_break
    integer :: num_manning
    real(CLAW_REAL) :: friction_depth
    
    ! Method parameters    
    real(CLAW_REAL) :: dry_tolerance

#ifdef CUDA
    real(kind=CLAW_REAL) :: grav_d, dry_tolerance_d, rho_d, earth_radius_d, deg2rad_d
    attributes(constant) :: grav_d, dry_tolerance_d, rho_d, earth_radius_d, deg2rad_d

    interface
    subroutine set_dry_tolerance(val) bind(C,name='set_dry_tolerance')
        use iso_c_binding
        implicit none
#if (CLAW_REAL == 8)
        real(kind=c_double), value, intent(in) :: val
#else
        real(kind=c_float), value, intent(in) :: val
#endif
    end subroutine set_dry_tolerance
    end interface
#endif


contains

    ! ========================================================================
    !  set_geo(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_geo(file_name)

        use utility_module, only: get_value_count
        
        use amr_module, only: mcapa, rinfinity

        implicit none

        ! Input
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: unit = 127

        integer :: i
        character(len=128) :: line

        if (.not.module_setup) then
            open(unit=GEO_PARM_UNIT,file='fort.geo',status="unknown",action="write")

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Physics Parameters:'
            write(GEO_PARM_UNIT,*) '-------------------'

            ! Read user parameters from setgeo.data
            if (present(file_name)) then
                call opendatafile(unit, file_name)
            else
                call opendatafile(unit, 'geoclaw.data')
            endif

            read(unit,*) grav
            read(unit,'(a)') line
            allocate(rho(get_value_count(line)))
            read(line,*) (rho(i), i=1, size(rho, 1))
            read(unit,*) rho_air
            read(unit,*) ambient_pressure
            read(unit,*) earth_radius
            read(unit,*) coordinate_system
            read(unit,*) sea_level
            read(unit,*)
            read(unit,*) coriolis_forcing
            if (coordinate_system == 1 .and. coriolis_forcing) then
                read(unit,*) theta_0
            else
                theta_0 = 0.d0
            endif
            read(unit,*) friction_forcing
            if (friction_forcing) then
                read(unit,*) num_manning
                allocate(manning_coefficient(num_manning))
                allocate(manning_break(num_manning))
                manning_break(num_manning) = rinfinity
                read(unit,*) manning_coefficient(:)
                read(unit,*) manning_break(1:num_manning-1)
                read(unit,*) friction_depth
            else
                friction_depth = rinfinity
            endif
            read(unit,*)
            read(unit,*) dry_tolerance
            
            close(unit)

            ! coordinate_system = 1 means Cartesian grid in meters
            ! coordinate_system = 2 means lat-long grid on sphere
            ! Check that coordinate_system is consistent with mcapa:
            if ((coordinate_system > 1) .and. (mcapa == 0)) then
                print *, 'ERROR in setgeo:  if coordinate_system > 1 then'
                print *, '      mcapa should be nonzero'
                stop
            endif
            if ((coordinate_system == 1) .and. (mcapa > 0)) then
                print *, 'ERROR in setgeo:  if coordinate_system = 1 then'
                print *, '      mcapa should be zero'
                stop
            endif

            write(GEO_PARM_UNIT,*) '   gravity:',grav
            write(GEO_PARM_UNIT,*) '   density water:',rho
            write(GEO_PARM_UNIT,*) '   density air:',rho_air
            write(GEO_PARM_UNIT,*) '   ambient pressure:',ambient_pressure
            write(GEO_PARM_UNIT,*) '   earth_radius:',earth_radius
            write(GEO_PARM_UNIT,*) '   coordinate_system:',coordinate_system
            write(GEO_PARM_UNIT,*) '   sea_level:',sea_level
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '   coriolis_forcing:',coriolis_forcing
            write(GEO_PARM_UNIT,*) '   theta_0:',theta_0
            write(GEO_PARM_UNIT,*) '   friction_forcing:',friction_forcing
            if (friction_forcing) then
                write(GEO_PARM_UNIT,*) '   manning_coefficient:', manning_coefficient
                write(GEO_PARM_UNIT,*) '   friction_depth:',friction_depth
            else
                write(GEO_PARM_UNIT,*) '   manning_coefficient: not used'
                write(GEO_PARM_UNIT,*) '   friction_depth: not used'
            end if

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '   dry_tolerance:', dry_tolerance

            module_setup = .true.
        end if
#ifdef CUDA
        grav_d = grav
        dry_tolerance_d = dry_tolerance
        rho_d = rho(1) ! assume only one layer of water
        earth_radius_d = earth_radius
        deg2rad_d = deg2rad
        call set_dry_tolerance(dry_tolerance)
#endif

    end subroutine set_geo


    ! ==========================================================================
    !  Calculate the coriolis constant f
    !   If coordinate_system == 1 then
    !       A beta-plane approximation is used and y should be in meters
    !   if coordinate_system == 2 then
    !       Grid is in lat-long and y should be in degrees which is converted
    !       to radians
    ! ==========================================================================
    real(CLAW_REAL) pure function coriolis(y)

        implicit none
        
        ! Input
        real(CLAW_REAL), intent(in) :: y
        
        ! Locals
        real(CLAW_REAL) :: theta
        
        ! Assume beta plane approximation and y is in meters    
        if (coordinate_system == 1) then
            theta = y / 111d3 * DEG2RAD + theta_0
            coriolis = 2.d0 * omega * (sin(theta_0) + (theta - theta_0)     &
                                                    * cos(theta_0))
        else if (coordinate_system == 2) then        
            coriolis = 2.d0 * omega * sin(y * DEG2RAD)
        else
            ! Unknown coordinate system, return nothing
            coriolis = 0.d0
        endif
    end function coriolis

    ! ==========================================================================
    !  Calculate the distance along a sphere
    !    real(CLAW_REAL) spherical_distance(x1,y1,x2,y2)
    !       x1 = (long,lat)
    !       x2 = (long,lat)
    ! ==========================================================================
    real(CLAW_REAL) pure function spherical_distance(x1,y1,x2,y2) result(distance)

        implicit none

        ! Input
        real(CLAW_REAL), intent(in) :: x1,y1,x2,y2

        ! Locals
        real(CLAW_REAL) :: dx ,dy

        dx = (x2 - x1) * DEG2RAD
        dy = (y2 - y1) * DEG2RAD

        distance = earth_radius * 2.d0 * asin(sqrt(sin(0.5d0*dy)**2 &
                                   + cos(y1 * DEG2RAD)*cos(y2 * DEG2RAD) &
                                   * sin(0.5d0*dx)**2))

    end function spherical_distance

    !=================================================================
    ! Transform long,lat --> (x,y) coordinates.
    !
    ! On input:
    !    coords(2) = (longitude (E/W),latitude (N/S))
    !    projection_center(2) = (longitude (E/W),latitude (N/S)) - coordinates 
    !                        where projection is true
    !
    ! On output:
    !    x(2)          (meters)
    !=================================================================
    pure function latlon2xy(coords,projection_center) result(x)

        real(CLAW_REAL), intent(in) :: coords(2), projection_center(2)
        real(CLAW_REAL) :: x(2)

        x(1) = deg2rad * earth_radius * (coords(1) - &
                    projection_center(1)) * cos(deg2rad * projection_center(2))
        x(2) = deg2rad * earth_radius * coords(2)

    end function latlon2xy

    !=================================================================
    ! Transform (x,y) --> (lat,lon) coordinates.
    !
    ! On input:
    !    x(2) = (meters)          
    !    projection_center(2) = (longitude (E/W),latitude (N/S)) - coordinates 
    !                        where projection is true
    !
    ! On output:
    !    coords(2) = (longitude,latitude)
    !=================================================================
    pure function xy2latlon(x,projection_center) result(coords)

        real(CLAW_REAL), intent(in) :: x(2), projection_center(2)
        real(CLAW_REAL) :: coords(2)

        coords(1) = projection_center(1) + x(1) &
                / (deg2rad * earth_radius * cos(deg2rad * projection_center(2)))
        coords(2) = x(2) / (deg2rad * earth_radius)
    end function xy2latlon

end module geoclaw_module
