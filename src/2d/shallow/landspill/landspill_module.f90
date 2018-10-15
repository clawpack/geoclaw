!> @file landspill_module.f90
!! @brief Top-level module for land-spill simulations. 
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-17

module landspill_module
    use point_source_collection_module
    implicit none
    save
    public

    !> @brief The state of this module.
    logical, private:: module_setup = .false.

    !> @brief Object for a collection of point sources.
    type(PointSourceCollection):: point_sources

contains

    !> @brief Initialize landspill module
    !! @param[in] point_source_file an optional arg; file for point sources
    subroutine set_landspill(point_source_file)
        character(len=*), intent(in), optional:: point_source_file 

        ! open data file for point source collection
        if (present(point_source_file)) then
            call point_sources%init(point_source_file)
        else
            call point_sources%init("point_source.data")
        endif

        ! set module_setup
        module_setup = .true.

    end subroutine set_landspill

end module landspill_module
