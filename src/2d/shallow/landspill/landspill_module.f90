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

    !> @brief Initialize landspill module
    !! @param[in] filename an optional arg; the input file
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

    !> @brief Add point source to the RHS of continuity equation.
    !! @param[in] meqn number of equations (the 1st dimension of variable q)
    !! @param[in] mbc number of ghost cell layers
    !! @param[in] mx number of cells in the x direction
    !! @param[in] my number of cells in the y direction
    !! @param[in] xlower the x-coordinate of the bottom-left corner of the mesh
    !! @param[in] ylower the y-coordinate of the bottom-left corner of the mesh
    !! @param[in] dx the cell size in x direction
    !! @param[in] dy the cell size in y direction
    !! @param[in] q the array holding values
    !! @param[in] t the current time
    !! @param[in] dt time-step size
    subroutine apply_point_sources( &
        meqn, mbc, mx, my, xlower, ylower, dx, dy, q, t, dt)

        ! declarations
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        integer(kind=4):: pti, i, j
        real(kind=8):: d

        ! code
        do pti = 1, nptss
            ! get indices of this point source on provided mesh
            call ptss(pti)%cell_id(mx, my, xlower, ylower, dx, dy, i, j)

            ! if this point source located outside this mesh, we skip this point
            if (i .eq. -999) cycle

            d = ptss(pti)%d_rate(t, dx, dy) ! get depth increment
            q(1, i, j) = q(1, i, j) + dt * d ! add to the continuity equation
        enddo

    end subroutine apply_point_sources

end module landspill_module
