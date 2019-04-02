!> @file point_source_module.f90
!! @brief PointSource class.
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-12

!> @brief Module for type declaration and implementation of point_source class.
!!
!! Only the point_source class is exposed to public. All others are private.
!!
module point_source_module
    implicit none
    private
    public:: PointSource

    !> @brief Type declaration of PointSource class.
    type:: PointSource
        private
        !> @brief A tag for this point source for future use.
        integer(kind=4):: id = -1
        !> @brief XY-coordinate of this point source.
        real(kind=8):: coord(2) = (/ 0D0, 0D0 /)
        !> @brief The number of time segments.
        integer(kind=4):: nt = -1
        !> @brief The end time of each time segment.
        real(kind=8), dimension(:), allocatable:: t
        !> @brief The volumetric value in each time segment.
        real(kind=8), dimension(:), allocatable:: v_rate

        contains
        !> @brief Reinitialize this instance directly.
        procedure, private:: init_direct
        !> @brief Reinitialize this instance with file unit.
        procedure, private:: init_from_funit
        !> @brief Calculate and return depth rate.
        procedure:: d_rate
        !> @brief Obtain the cell indices
        procedure:: cell_id
        !> @brief Get number of profile stages.
        procedure:: get_n_stages
        !> @brief Get time information of releasing fluids.
        procedure:: get_times
        !> @brief Set new times.
        procedure:: set_times
        !> @brief Get volumetric rates.
        procedure:: get_v_rates
        !> @brief Set volumetric rates.
        procedure:: set_v_rates
        !> @brief Write the information and data of this instance.
        procedure, private:: writef
        !> @brief Read the information and data of this instance.
        procedure, private:: readf
        !> @brief A generic interface for initialization
        generic:: init => init_direct, init_from_funit
        !> @brief Overriding intrinsic write function.
        generic:: write(formatted) => writef
        !> @brief Overriding intrinsic read function.
        generic:: read(formatted) => readf
        !> @brief Destructor.
        final:: destructor
    end type PointSource

    !> @brief An overloading of C++ style constructor.
    interface PointSource
        procedure constructor_direct
        procedure constructor_funit
    end interface PointSource

contains

    !> @brief A mimic to C++ stype constructor.
    !! @param[in] coord a 1D real(kind=8) array with size 2, representing x and y coordinates
    !! @param[in] nt number of time segments
    !! @param[in] t a 1D array of reals with size nt; times at the end of segments
    !! @param[in] v_rate a 1D array of reals with size nt; values
    !! @return a new instance with provided setup.
    function constructor_direct(id, coord, nt, t, v_rate)
        ! variable declaration
        type(PointSource):: constructor_direct
        integer(kind=4), intent(in):: id
        real(kind=8), dimension(2), intent(in):: coord
        integer(kind=4), intent(in):: nt
        real(kind=8), dimension(nt), intent(in):: t
        real(kind=8), dimension(nt), intent(in):: v_rate

        call init_direct(constructor_direct, id, coord, nt, t, v_rate)
    end function constructor_direct

    !> @brief A mimic to C++ stype constructor.
    !! @param[in] funit the unit of the file.
    !! @return a new instance with provided setup.
    function constructor_funit(funit)
        ! variable declaration
        type(PointSource):: constructor_funit
        integer(kind=4), intent(in):: funit

        call init_from_funit(constructor_funit, funit)
    end function constructor_funit

    !> @brief Finalize a PointSource instance.
    !!
    !! Deallocate dynamic arrays and finalize the instance.
    !!
    !! @param[in, out] this an instance of PointSource
    subroutine destructor(this)
        ! declaration
        type(PointSource), intent(inout):: this

        ! code
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%v_rate)) deallocate(this%v_rate)

        this%coord = (/ 0D0, 0D0 /)
        this%nt = -1
        this%id = -1

    end subroutine destructor

    !> @brief Initialize a PointSource instance.
    !!
    !! Initializ a PointSource instance.
    !!
    !! @param[in, out] this an instance of PointSource
    !! @param[in] id ID/tag of this object.
    !! @param[in] coord a 1D real(kind=8) array with size 2, representing x and y coordinates
    !! @param[in] nt number of time segments
    !! @param[in] t a 1D array of reals with size nt; times at the end of segments
    !! @param[in] v_rate a 1D array of reals with size nt; values
    subroutine init_direct(this, id, coord, nt, t, v_rate)
        ! variable declaration
        class(PointSource), intent(inout):: this
        integer(kind=4), intent(in):: id
        real(kind=8), dimension(2), intent(in):: coord
        integer(kind=4), intent(in):: nt
        real(kind=8), dimension(nt), intent(in):: t
        real(kind=8), dimension(nt), intent(in):: v_rate

        ! =====================================================================
        ! Note: t array should be already sorted. Using GeoClaw's preprocessing
        !       mechanism should have already checked this.
        ! =====================================================================

        ! code
        call destructor(this)

        this%id = id
        this%coord = coord
        this%nt = nt
        this%t = t ! implicit allocation
        this%v_rate = v_rate ! implicit allocation

    end subroutine init_direct

    !> @brief Initialize a PointSource instance with a file unit provided.
    !!
    !! Initializ a PointSource instance.
    !!
    !! @param[in, out] this an instance of PointSource
    !! @param[in] funit unit of the file.
    subroutine init_from_funit(this, funit)
        ! variable declaration
        class(PointSource), intent(inout):: this
        integer(kind=4), intent(in):: funit

        ! code
        call destructor(this)

        read(funit, *) this%id
        read(funit, *) this%coord
        read(funit, *) this%nt

        allocate(this%t(this%nt), this%v_rate(this%nt))

        read(funit, *) this%t
        read(funit, *) this%v_rate

        ! TODO: make sure the t array is sorted and v_rate is sorted based on t

    end subroutine init_from_funit

    !> @brief Print information of an instance of PointSource.
    !!
    !! Print information of an instance of PointSource to specified file unit.
    !!
    !! @param[in] this an instance of PointSource
    !! @param[in] iounit see intrinsic write function
    !! @param[in] iotype see intrinsic write function
    !! @param[in] v_list see intrinsic write function
    !! @param[out] stat see intrinsic write function
    !! @param[in, out] msg see intrinsic write function
    subroutine writef(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(PointSource), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        ! code
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%id, t, t, t, t, "=: id # ID of this point source"
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%coord, "=: coord # coordinates"
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%nt, t, t, t, t, "=: n_times # number of time segments"
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%t, "=: end_times # end times of segments"
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%v_rate, "=: vol_rates # volumetric rates of segments"

    end subroutine writef

    !> @brief Read formatted data block and re-initialize the instance.
    !!
    !! Read formatted data block and re-initialize the instance.
    !!
    !! @param[in] this an instance of PointSource
    !! @param[in] iounit see intrinsic write function
    !! @param[in] iotype see intrinsic write function
    !! @param[in] v_list see intrinsic write function
    !! @param[out] stat see intrinsic write function
    !! @param[in, out] msg see intrinsic write function
    subroutine readf(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(PointSource), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        ! code
        print *, "Direct read of this object is prohibited!"
        stop

    end subroutine readf

    !> @brief Return a depth rate according to time, dx and dy.
    !!
    !! Calculate depth rate based on given time, dx and dy.
    !!
    !! @param[in] this an instance of PointSource
    !! @param[in] t time
    !! @param[in] dx cell size in x direction
    !! @param[in] dy cell size in y direction
    !! @return depth rate
    function d_rate(this, t, dx, dy) result(rate)
        ! declarations
        class(PointSource), intent(in):: this
        real(kind=8), intent(in):: t, dx, dy
        real(kind=8):: rate
        integer(kind=4):: i

        ! code
        i = count(this%t<t) + 1

        if (i > this%nt) then
            rate = 0D0
        else
            rate = this%v_rate(i) / (dx * dy)
        endif

    end function d_rate

    !> @brief Find the (i, j) of this point source in a provided mesh.
    !!
    !! Find the index (i, j) in a provided mesh where this point source located
    !! in. If the returned i and j are -999, this means this point source is not
    !! located in this mesh.
    !!
    !! @param[in] this an instance of PointSource
    !! @param[in] mbuff the number of buffer layers.
    !! @param[in] mx the number of cells in x direction.
    !! @param[in] my the number of cells in y direction.
    !! @param[in] xlower the lower limit in x direction of this mesh.
    !! @param[in] ylower the lower limit in y direction of this mesh.
    !! @param[in] dx the cell size in x direction.
    !! @param[in] dy the cell size in y direction.
    !! @param[out] i the cell index in x direction.
    !! @param[out] j the cell index in y direction.
    subroutine cell_id(this, mbuff, mx, my, xlower, ylower, dx, dy, i, j)
        ! declaration
        class(PointSource), intent(in):: this
        integer(kind=4), intent(in):: mbuff, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy
        integer(kind=4), intent(out):: i, j
        real(kind=8):: xleft, xright, ybot, ytop

        ! code
        i = -999
        j = -999

        xleft = xlower - mbuff * dx - dx * 1d-6
        xright = xlower + (mx + mbuff) * dx + dx * 1d-6
        ybot = ylower - mbuff * dy - dx * 1d-6
        ytop = ylower + (my + mbuff) * dy + dx * 1d-6

        if ((this%coord(1) >= xleft) .and. (this%coord(1) < xright)) then
            if ((this%coord(2) >= ybot) .and. (this%coord(2) < ytop)) then
                i = int((this%coord(1)-xleft)/dx) + 1 - mbuff
                j = int((this%coord(2)-ybot)/dy) + 1 - mbuff
            endif
        endif

    end subroutine cell_id

    ! get_n_stages
    function get_n_stages(this) result(ans)
        class(PointSource), intent(in):: this
        integer(kind=4):: ans
        ans = this%nt
    end function get_n_stages

    ! get_times
    subroutine get_times(this, times)
        class(PointSource), intent(in):: this
        real(kind=8), dimension(this%nt), intent(inout):: times
        times = this%t
    end subroutine get_times

    ! set_times
    subroutine set_times(this, times)
        class(PointSource), intent(inout):: this
        real(kind=8), dimension(this%nt), intent(in):: times
        this%t = times
    end subroutine set_times

    ! get_v_rates
    subroutine get_v_rates(this, rates)
        class(PointSource), intent(in):: this
        real(kind=8), dimension(this%nt), intent(inout):: rates
        rates = this%v_rate
    end subroutine get_v_rates

    ! set_v_rates
    subroutine set_v_rates(this, rates)
        class(PointSource), intent(inout):: this
        real(kind=8), dimension(this%nt), intent(in):: rates
        this%v_rate = rates
    end subroutine set_v_rates

end module point_source_module
