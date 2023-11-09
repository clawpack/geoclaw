! ============================================================================
!  File:        geoclaw_mod
! ============================================================================
!    Copyright (C) 2010-04-21 Clawpack Developers http://www.clawpack.org
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

! Adapted from 2d GeoClaw by deleting some thing that are unused in 1d.
! Other things still remain that are note used in 1d, but they are written to
! data files by default and so still read in here.

module geoclaw_module

    implicit none
    save

    logical, private :: module_setup = .false.
    
    ! ========================================================================
    !  Constants
    ! ========================================================================
    integer, parameter :: GEO_PARM_UNIT = 78
    integer, parameter :: KAPPA_UNIT = 42
    real(kind=8), parameter :: pi = 4.d0*datan(1.d0)
    real(kind=8), parameter :: DEG2RAD = pi / 180.d0
    real(kind=8), parameter :: RAD2DEG = 180.d0 / pi
    
    ! ========================================================================
    !  Physics
    ! ========================================================================
    real(kind=8) :: grav, rho_air, ambient_pressure, earth_radius, sea_level
    ! Water density can be an array to handle multiple layers
    real(kind=8), allocatable :: rho(:)
    integer :: coordinate_system, sphere_source

    ! Rotational velocity of Earth
    real(kind=8), parameter :: omega = 2.0d0 * pi / 86164.2d0
    
    ! Forcing
    logical :: coriolis_forcing ! True then coriolis terms included in src
    real(kind=8) :: theta_0 ! Used if using the beta-plane approximation
    logical :: friction_forcing ! Friction forcing will be applied
    real(kind=8), dimension(:),allocatable :: manning_coefficient, manning_break
    integer :: num_manning
    real(kind=8) :: friction_depth
    real(kind=8) :: friction_coefficient ! used in src2 rather than manning
    
    ! Method parameters    
    real(kind=8) :: dry_tolerance

contains

    ! ========================================================================
    !  set_geo(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_geo(file_name)

        use utility_module, only: get_value_count
        
        !use amr_module, only: mcapa, rinfinity

        implicit none

        ! Input
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: unit = 127

        integer :: i
        character(len=128) :: line

        real(kind=8), parameter :: rinfinity = 10.e32

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
            read(unit,*) sphere_source
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
                if (num_manning > 1) then
                    write(6,*) '*** num_manning > 1 not yet implemented'
                    stop
                    endif
                allocate(manning_coefficient(num_manning))
                allocate(manning_break(num_manning))
                manning_break(num_manning) = rinfinity
                read(unit,*) manning_coefficient(:)
                read(unit,*) manning_break(1:num_manning-1)
                read(unit,*) friction_depth
                friction_coefficient = manning_coefficient(1) ! used in src2
            else
                friction_depth = rinfinity
            endif
            read(unit,*)
            read(unit,*) dry_tolerance
            
            close(unit)


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

    end subroutine set_geo

end module geoclaw_module
