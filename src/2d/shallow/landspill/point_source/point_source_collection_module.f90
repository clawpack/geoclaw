!> @file point_source_collection_module.f90
!! @brief PointSourceCollection class.
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-10-03


!> @brief Definition and implementations of PointSourceCollection class.
module point_source_collection_module
    use point_source_module
    implicit none
    private
    public:: PointSourceCollection

    !> @brief A class for a collection of point sources
    type:: PointSourceCollection
        private ! member variables
        !> @brief Number of point sources in the collection.
        integer(kind=4):: npts = -1
        !> @brief Array of point source objects.
        type(PointSource), allocatable, dimension(:):: pts

        contains ! member functions
        !> @brief Initialization with an input file.
        procedure:: init
        !> @brief Underlying output driver.
        procedure, private:: write_data
        !> @brief Underlying input driver.
        procedure, private:: read_data
        !> @brief Apply point sources to the RHS of equations.
        procedure:: apply_point_sources
        !> @brief Apply tiny amount depth to IC.
        procedure:: apply_trivial_ic
        !> @brief Flag cells to be refined
        procedure:: flag_cells
        !> @brief Get number of point sources.
        procedure:: get_n_points
        !> @brief Get number of segments of a point.
        procedure:: get_n_stages
        !> @brief Get time information of releasing fluids of a point.
        procedure:: get_times
        !> @brief Set new times of a point.
        procedure:: set_times
        !> @brief Get volumetric rates of a point.
        procedure:: get_v_rates
        !> @brief Set volumetric rates of a point.
        procedure:: set_v_rates
        !> @brief Overload intrinsic write.
        generic:: write(formatted) => write_data
        !> @brief Overload intrinsic write.
        generic:: read(formatted) => read_data
        !> @brief Destructor.
        final:: destructor
    end type PointSourceCollection

    !> @brief Constructor of PointSourceCollection.
    interface PointSourceCollection
        procedure constructor
    end interface PointSourceCollection

contains

    ! implementation of constructor
    function constructor(filename)
        type(PointSourceCollection):: constructor
        character(len=*), intent(in), optional:: filename

        if (present(filename)) then
            call constructor%init(filename)
        else
            call constructor%init()
        endif
    end function constructor

    ! implementation of init
    subroutine init(this, filename)
        use geoclaw_module, only: coordinate_system
        class(PointSourceCollection), intent(inout):: this
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
        read(funit, *) this%npts

        ! if there's no point source, exit
        if (this%npts == 0) then
            close(funit)
            return
        endif

        ! allocate the array of point source instances
        allocate(this%pts(this%npts))

        ! read data for each point source
        do i = 1, this%npts
            call this%pts(i)%init(funit)
        enddo

        ! close the file
        close(funit)

        ! so far, we only support xy coordinates (Cartesian)
        if (coordinate_system .ne. 1) then
            print *, "Point source functionality now only works with &
                Cartesian coordinates."
            stop
        endif
    end subroutine init

    ! implementation of destructor
    subroutine destructor(this)
        type(PointSourceCollection), intent(inout):: this

        this%npts = -1
        if (allocated(this%pts)) deallocate(this%pts)
    end subroutine destructor

    subroutine write_data(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) n
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%npts, t, t, t, t, "=: n_point_sources # Number of point sources"

        do i = 1, this%npts
            write(iounit, *, iostat=stat, iomsg=msg) n
            write(iounit, "(DT)", iostat=stat, iomsg=msg) this%pts(i)
        enddo

    end subroutine write_data

    subroutine read_data(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(PointSourceCollection), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        print *, "Direct read of this object is prohibited!"
        stop

    end subroutine read_data

    !> @brief Add point source to the RHS of continuity equation.
    !! @param[in[ this a PointSourceCollection object.
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
    subroutine apply_point_sources(this, &
        meqn, mbc, mx, my, xlower, ylower, dx, dy, q, t, dt)

        ! declarations
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        integer(kind=4):: pti, i, j
        real(kind=8):: d

        ! code
        do pti = 1, this%npts
            ! get indices of this point source on provided mesh
            call this%pts(pti)%cell_id(mbc, mx, my, xlower, ylower, dx, dy, i, j)

            ! if this point source located outside this mesh, we skip this point
            if (i == -999) cycle

            d = this%pts(pti)%d_rate(t, dx, dy) ! get depth increment
            q(1, i, j) = q(1, i, j) + dt * d ! add to the continuity equation
        enddo

    end subroutine apply_point_sources

    !> @brief Add a trivial amount depth to IC.
    !! @param[in[ this a PointSourceCollection object.
    !! @param[in] meqn number of equations (the 1st dimension of variable q)
    !! @param[in] mbc number of ghost cell layers
    !! @param[in] mx number of cells in the x direction
    !! @param[in] my number of cells in the y direction
    !! @param[in] xlower the x-coordinate of the bottom-left corner of the mesh
    !! @param[in] ylower the y-coordinate of the bottom-left corner of the mesh
    !! @param[in] dx the cell size in x direction
    !! @param[in] dy the cell size in y direction
    !! @param[in] q the array holding values
    subroutine apply_trivial_ic(this, &
        meqn, mbc, mx, my, xlower, ylower, dx, dy, q)
        use geoclaw_module, only: dry_tolerance

        ! declarations
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        integer(kind=4):: pti, i, j

        ! code
        do pti = 1, this%npts
            ! get indices of this point source on provided mesh
            call this%pts(pti)%cell_id(mbc, mx, my, xlower, ylower, dx, dy, i, j)

            ! if this point source located outside this mesh, we skip this point
            if (i == -999) cycle

            if ((q(1, i, j) / dry_tolerance) < 0.9) then
                q(1, i, j) = dry_tolerance * 0.9
            endif
        enddo

    end subroutine apply_trivial_ic

    ! flag_cells
    subroutine flag_cells(this, &
        mbuff, mx, my, xlower, ylower, dx, dy, amrflags)
        use amr_module, only: DOFLAG

        ! declarations
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: mbuff, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy
        real(kind=8), intent(inout):: amrflags(1-mbuff:mx+mbuff, 1-mbuff:my+mbuff)
        integer(kind=4):: pti, i, j

        ! code
        do pti = 1, this%npts
            ! get indices of this point source on provided mesh
            call this%pts(pti)%cell_id(mbuff, mx, my, xlower, ylower, dx, dy, i, j)

            ! if this point source located outside this mesh, we skip this point
            if (i == -999) cycle

            ! flag the cell
            amrflags(i, j) = DOFLAG
        enddo
    end subroutine flag_cells

    ! get_n_points.
    function get_n_points(this) result(ans)
        class(PointSourceCollection), intent(in):: this
        integer(kind=4):: ans
        ans = this%npts
    end function get_n_points

    ! get_n_stages
    function get_n_stages(this, i) result(ans)
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: i
        integer(kind=4):: ans
        ans = this%pts(i)%get_n_stages()
    end function get_n_stages

    ! get_times
    subroutine get_times(this, i, times)
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: i
        real(kind=8), dimension(:), intent(inout):: times
        call this%pts(i)%get_times(times)
    end subroutine get_times

    ! set_times
    subroutine set_times(this, i, times)
        class(PointSourceCollection), intent(inout):: this
        integer(kind=4), intent(in):: i
        real(kind=8), dimension(:), intent(in):: times
        call this%pts(i)%set_times(times)
    end subroutine set_times

    ! get_v_rates
    subroutine get_v_rates(this, i, rates)
        class(PointSourceCollection), intent(in):: this
        integer(kind=4), intent(in):: i
        real(kind=8), dimension(:), intent(inout):: rates
        call this%pts(i)%get_v_rates(rates)
    end subroutine get_v_rates

    ! set_v_rates
    subroutine set_v_rates(this, i, rates)
        class(PointSourceCollection), intent(inout):: this
        integer(kind=4), intent(in):: i
        real(kind=8), dimension(:), intent(in):: rates
        call this%pts(i)%set_v_rates(rates)
    end subroutine set_v_rates

end module point_source_collection_module
