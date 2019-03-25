!> @file landspill_module.f90
!! @brief Top-level module for land-spill simulations.
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-17

module landspill_module
    use point_source_collection_module
    use darcy_weisbach_module
    use hydro_feature_collection_module
    use evap_module
    implicit none
    save
    public
    private:: get_kinematic_viscosity, landspill_log

    !> @brief The state of this module.
    logical, private:: module_setup = .false.

    !> @brief Kinematic viscosity (m^2/s).
    real(kind=8):: nu = 0D0

    !> @brief Reference dynamic viscosity (mPa-s = 1e-3 kg/s/m = cP).
    real(kind=8):: ref_mu = 0D0

    !> @brief Reference temperature (Celsius).
    real(kind=8):: ref_temperature = 0D0

    !> @brief Ambient temperature (Celsius).
    real(kind=8):: ambient_temperature = 0D0

    !> @brief Density at reference temperature.
    real(kind=8):: density = 0D0


    !> @brief File name of point source settings.
    character(len=:), allocatable:: point_source_file

    !> @brief Object for a collection of point sources.
    type(PointSourceCollection):: point_sources


    !> @brief File name of Darcy-Weisbach settings.
    character(len=:), allocatable:: darcy_weisbach_file

    !> @brief Object for Darcy-Weisbach.
    type(DarcyWeisbach):: darcy_weisbach


    !> @brief File name of hydrological feature settings.
    character(len=:), allocatable:: hydro_feature_file

    !> @brief Object for a collection of hydro features.
    type(HydroFeatureCollection):: hydro_features


    !> @brief File name of evaporation settings.
    character(len=:), allocatable:: evaporation_file

    !> @brief Object for a collection of hydro features.
    type(EvapModel):: evaporation

contains

    !> @brief Initialize landspill module
    !! @param[in] landspill_file an optional arg; file for landspill module
    subroutine set_landspill(landspill_file)
        use geoclaw_module, only: geo_module_setup => coordinate_system
        use geoclaw_module, only: geo_friction => friction_forcing
        use geoclaw_module, only: geo_rho => rho

        ! arguments
        character(len=*), intent(in), optional:: landspill_file
        integer(kind=4), parameter:: funit = 200

        ! We ASSUME that coordinate_system in GeoClaw is initially zero.
        ! So if the geoclaw module is set, then coordinate_system should be > 0.
        if (geo_module_setup == 0) then
            print *, "Error: geoclaw module should be set up " // &
                "before setting landspill module"
            stop
        end if

        ! open module configuration file
        if (present(landspill_file)) then
            call opendatafile(funit, landspill_file)
        else
            call opendatafile(funit, "landspill.data")
        endif

        ! read data
        read(funit, *) ref_mu
        read(funit, *) ref_temperature
        read(funit, *) ambient_temperature
        read(funit, *) density
        read(funit, *) point_source_file
        read(funit, *) darcy_weisbach_file
        read(funit, *) hydro_feature_file
        read(funit, *) evaporation_file

        ! close file
        close(funit)

        ! overwrite the density in geoclaw module
        geo_rho = density

        ! get kinematic viscosity at ambient temperature
        nu = get_kinematic_viscosity( &
            ref_mu, ref_temperature, ambient_temperature, density)


        ! point source collection
        call point_sources%init(point_source_file)

        ! Darcy-Weisbach
        call darcy_weisbach%init(nu, darcy_weisbach_file)

        ! Disable Manning's friction in GeoClaw
        ! TODO: this is just a temporary workaround. Need to integrate to GeoClaw.
        if (darcy_weisbach%get_type() .gt. 0) geo_friction = .false.

        ! hydro feature collection
        call hydro_features%init(hydro_feature_file)

        ! evaporation
        call evaporation%init(ambient_temperature, evaporation_file)
        call evaporation%reset_release_profile(point_sources)

        ! set module_setup
        module_setup = .true.

        call landspill_log()

    end subroutine set_landspill

    !> @brief Calculate kinematic viscosity.
    function get_kinematic_viscosity(mu_ref, T_ref, T, rho) result(nu)
        implicit none
        real(kind=8), intent(in):: mu_ref ! cP, i.e., 1e-3 kg / m /s
        real(kind=8), intent(in):: T_ref ! Celsius
        real(kind=8), intent(in):: T ! Celsius
        real(kind=8), intent(in):: rho ! kg / m^3
        real(kind=8):: nu

        ! get dynamic viscosity at ambient temperature (unit: cP)
        nu = mu_ref**(-0.2661D0) + (T - T_ref) / 233D0
        nu = - dlog(nu) / 0.2661D0
        nu = dexp(nu)

        ! convert to kg / s / m
        nu = nu * 1D-3

        ! get kinematic viscosity (m^2 / s)
        nu = nu / rho
    end function get_kinematic_viscosity

    !> @brief output a log
    subroutine landspill_log()
        implicit none
        integer(kind=4), parameter:: funit = 199
        integer(kind=4):: i, npts, nt
        real(kind=8), allocatable, dimension(:):: times, rates

        open(unit=funit, file="landspill.log", action="write")

        write(funit, *) "Reference Dynamic Viscosity (cP): ", ref_mu
        write(funit, *) "Reference Temperature (Celsius): ", ref_temperature
        write(funit, *) "Ambient Temperature (Celsius): ", ambient_temperature
        write(funit, *) "Density (kg / m^3): ", density
        write(funit, *) "Calculated kinematic viscosity (m^2 / s): ", nu

        write(funit, *) "Evaporation type: ", evaporation%get_type()
        write(funit, *) "Remained percentage from 60min to 61min: ", &
            evaporation%remained_percentage(36D2, 6D1)

        npts = point_sources%get_n_points()
        write(funit, *) "Number of point sources: ", npts

        do i = 1, npts
            nt = point_sources%get_n_stages(i)
            write(funit, *) "Number of stages of the point ", i, ": ", nt

            allocate(times(nt), rates(nt))
            call point_sources%get_times(i, times)
            call point_sources%get_v_rates(i, rates)
            write(funit, *) "Stage times of the point ", i, ": ", times
            write(funit, *) "Stage rates of the point ", i, ": ", rates
            deallocate(times, rates)
        end do
        close(funit)
    end subroutine landspill_log

end module landspill_module
