! ============================================================================
!  File:        geoclaw_module (for 1D)
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
    real(kind=8) :: grav, earth_radius, sea_level
    integer :: coordinate_system

    ! Forcing
    ! friction
    real(kind=8) :: friction_depth, friction_coefficient, friction_forcing

    ! Method parameters
    real(kind=8) :: dry_tolerance

contains

    ! ========================================================================
    !  set_geo(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_geo(file_name)

        implicit none

        ! Input
        character(len=*), intent(in), optional :: file_name

        ! Locals
        integer, parameter :: unit = 127

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
        read(unit,*) earth_radius
        read(unit,*) coordinate_system
        read(unit,*) sea_level
        read(unit,*) friction_forcing
        read(unit,*) friction_coefficient
        read(unit,*) friction_depth
        read(unit,*) dry_tolerance

        close(unit)

        write(GEO_PARM_UNIT,*) '   gravity:',grav
        write(GEO_PARM_UNIT,*) '   earth_radius:',earth_radius
        write(GEO_PARM_UNIT,*) '   coordinate_system:',coordinate_system
        write(GEO_PARM_UNIT,*) '   sea_level:',sea_level
        write(GEO_PARM_UNIT,*) '   friction_forcing:',friction_forcing
        write(GEO_PARM_UNIT,*) '   friction_coefficient:', friction_coefficient
        write(GEO_PARM_UNIT,*) '   friction_depth:',friction_depth
        write(GEO_PARM_UNIT,*) '   dry_tolerance:',dry_tolerance

    end subroutine set_geo


end module geoclaw_module



