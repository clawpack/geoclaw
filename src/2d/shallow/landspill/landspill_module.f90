!> @file landspill_module.f90
!! @brief Top-level module for land-spill simulations. 
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-17

module landspill_module
    use point_source_collection_module
    use darcy_weisbach_module
    use hydro_feature_collection_module
    implicit none
    save
    public

    !> @brief The state of this module.
    logical, private:: module_setup = .false.

    !> @brief Object for a collection of point sources.
    type(PointSourceCollection):: point_sources

    !> @brief Object for Darcy-Weisbach.
    type(DarcyWeisbach):: darcy_weisbach

    !> @brief Object for a collection of hydro features.
    type(HydroFeatureCollection):: hydro_features

contains

    !> @brief Initialize landspill module
    !! @param[in] point_source_file an optional arg; file for point sources
    !! @param[in] darcy_weisbach_file an optional arg; file for Darcy-Weisbach
    !! @param[in] hydro_feature_file an optional arg; file for hydro features
    subroutine set_landspill(&
        point_source_file, darcy_weisbach_file, hydro_feature_file)
        use geoclaw_module, only: geo_friction => friction_forcing 

        ! arguments
        character(len=*), intent(in), optional:: point_source_file 
        character(len=*), intent(in), optional:: darcy_weisbach_file 
        character(len=*), intent(in), optional:: hydro_feature_file 

        ! point source collection
        if (present(point_source_file)) then
            call point_sources%init(point_source_file)
        else
            call point_sources%init("point_source.data")
        endif

        ! Darcy-Weisbach
        if (present(darcy_weisbach_file)) then
            call darcy_weisbach%init(darcy_weisbach_file)
        else
            call darcy_weisbach%init("darcy_weisbach.data")
        endif

        ! Disable Manning's friction in GeoClaw
        ! Warning: this assumes geoclaw module has been set up!
        ! TODO: this is just a temporary workaround. Need to integrate to GeoClaw.
        if (darcy_weisbach%get_type() > 0) then
            geo_friction = .false.
        endif

        ! hydro feature collection
        if (present(hydro_feature_file)) then
            call hydro_features%init(hydro_feature_file)
        else
            call hydro_features%init("hydro_feature.data")
        endif

        ! set module_setup
        module_setup = .true.

    end subroutine set_landspill

end module landspill_module
