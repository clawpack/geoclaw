!> @file landspill_module.f90
!! @brief Top-level module for land-spill simulations. 
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-17

module landspill_module
    use point_source_module, only: point_source
    implicit none
    save
    public

    !> @brief The state of this module.
    logical, private:: module_setup = .false.

    !> @brief Number of point sources
    integer(kind=4):: nptss

    !> @brief Point source onjects
    type(point_source), allocatable, dimension(:):: ptss

contains

    subroutine set_landspill(filename)
        use geoclaw_module, only: coordinate_system
        character(len=*), intent(in), optional:: filename 
        integer(kind=4), parameter:: funit = 256 ! local file unit
        integer(kind=4):: i, id

        ! open data file
        if (present(filename)) then
            call opendatafile(funit, filename)
        else
            call opendatafile(funit, "point_source.data")
        endif

        ! read number of point sources
        read(funit, *) nptss

        ! if there's no point source, exit
        if (nptss .eq. 0) then
            module_setup = .true.
            close(funit)
            return
        endif

        ! so far, we only support xy coordinates (Cartesian)
        if (coordinate_system .ne. 1) then
            write(*, *) "Point source functionality now only works with &
                Cartesian coordinates."
            stop
        endif

        ! allocate the array of point source instances
        allocate(ptss(nptss))

        ! read data for each point source
        do i = 1, nptss
            read(funit, *) id
            read(funit, "(DT)", advance='no') ptss(i)
        enddo

        ! close the file
        close(funit)

    end subroutine set_landspill

end module landspill_module
