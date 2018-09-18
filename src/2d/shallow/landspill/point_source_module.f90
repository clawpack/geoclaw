!> @file point_source_module.f90
!! @brief Point source class. 
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
    public:: point_source

    !> @brief Type declaration of point_source class.
    type point_source
        !> @brief XY-coordinate of this point source.
        real(kind=8):: coord(2) = (/ 0D0, 0D0 /)
        !> @brief The number of time segments.
        integer(kind=4):: nt = -1
        !> @brief The end time of each time segment.
        real(kind=8), dimension(:), allocatable:: t
        !> @brief The volumetric value in each time segment.
        real(kind=8), dimension(:), allocatable:: v_rate

        contains
        !> @brief Reinitialize this instance.
        procedure:: reinit => init
        !> @brief Manually clear this instance.
        procedure:: clear
        !> @brief Calculate and return depth rate.
        procedure:: d_rate
        !> @brief Write the information and data of this instance.
        procedure:: writef
        !> @brief Read the information and data of this instance.
        procedure:: readf
        !> @brief Overriding intrinsic write function.
        generic:: write(formatted) => writef
        !> @brief Overriding intrinsic read function.
        generic:: read(formatted) => readf
        !> @brief Destructor.
        final:: destructor
    end type point_source

    !> @brief An overloading of C++ style constructor.
    interface point_source
        procedure constructor
    end interface point_source

contains

    !> @brief A mimic to C++ stype constructor.
    !! @param[in] coord a 1D real(kind=8) array with size 2, representing x and y coordinates
    !! @param[in] nt number of time segments
    !! @param[in] t a 1D array of reals with size nt; times at the end of segments
    !! @param[in] v_rate a 1D array of reals with size nt; values
    !! @return a new instance with provided setup.
    function constructor(coord, nt, t, v_rate)
        ! variable declaration
        type(point_source):: constructor
        real(kind=8), dimension(2), intent(in):: coord
        integer(kind=4), intent(in):: nt
        real(kind=8), dimension(nt), intent(in):: t
        real(kind=8), dimension(nt), intent(in):: v_rate

        call init(constructor, coord, nt, t, v_rate)
    end function constructor

    !> @brief Finalize a point_source instance.
    !!
    !! Deallocate dynamic arrays and finalize the instance.
    !!
    !! @param[in] this an instance of point_source
    subroutine destructor(this)
        ! declaration
        type(point_source), intent(inout):: this

        ! code
        call this%clear()

    end subroutine destructor

    !> @brief Initialize a point_source instance.
    !!
    !! Initializ a point_source instance.
    !!
    !! @param[in, out] this an instance of point_source
    !! @param[in] coord a 1D real(kind=8) array with size 2, representing x and y coordinates
    !! @param[in] nt number of time segments
    !! @param[in] t a 1D array of reals with size nt; times at the end of segments
    !! @param[in] v_rate a 1D array of reals with size nt; values
    subroutine init(this, coord, nt, t, v_rate)
        ! variable declaration
        class(point_source), intent(inout):: this
        real(kind=8), dimension(2), intent(in):: coord
        integer(kind=4), intent(in):: nt
        real(kind=8), dimension(nt), intent(in):: t
        real(kind=8), dimension(nt), intent(in):: v_rate

        ! code
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%v_rate)) deallocate(this%t)

        this%coord = coord
        this%nt = nt
        this%t = t ! implicit allocation
        this%v_rate = v_rate ! implicit allocation

        ! TODO: make sure the t array is sorted and v_rate is sorted based on t

    end subroutine init

    !> @brief Manually clear a point_source instance.
    !!
    !! Deallocate dynamic arrays and finalize the instance.
    !!
    !! @param[in] this an instance of point_source
    subroutine clear(this)
        ! declaration
        class(point_source), intent(inout):: this

        ! code
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%v_rate)) deallocate(this%v_rate)

        this%coord = (/ 0D0, 0D0 /)
        this%nt = -1

    end subroutine clear

    !> @brief Print information of an instance of point_source.
    !!
    !! Print information of an instance of point_source to specified file unit.
    !!
    !! @param[in] this an instance of point_source
    !! @param[in] iounit see intrinsic write function
    !! @param[in] iotype see intrinsic write function
    !! @param[in] v_list see intrinsic write function
    !! @param[out] stat see intrinsic write function
    !! @param[in, out] msg see intrinsic write function
    subroutine writef(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(point_source), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        character:: c

        ! code
        write(iounit, *, iostat=stat, iomsg=msg) new_line(c), this%coord, &
            new_line(c), this%nt, new_line(c), this%t, new_line(c), this%v_rate

    end subroutine writef

    !> @brief Read formatted data block and re-initialize the instance.
    !!
    !! Read formatted data block and re-initialize the instance.
    !!
    !! @param[in] this an instance of point_source
    !! @param[in] iounit see intrinsic write function
    !! @param[in] iotype see intrinsic write function
    !! @param[in] v_list see intrinsic write function
    !! @param[out] stat see intrinsic write function
    !! @param[in, out] msg see intrinsic write function
    subroutine readf(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(point_source), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        ! code
        read(iounit, *, iostat=stat, iomsg=msg) this%coord
        read(iounit, *, iostat=stat, iomsg=msg) this%nt

        allocate(this%t(this%nt), this%v_rate(this%nt))

        read(iounit, *, iostat=stat, iomsg=msg) this%t
        read(iounit, *, iostat=stat, iomsg=msg) this%v_rate

    end subroutine readf

    !> @brief Return a depth rate according to time, dx and dy.
    !!
    !! Calculate depth rate based on given time, dx and dy.
    !!
    !! @param[in] this an instance of point_source
    !! @param[in] t time
    !! @param[in] dx cell size in x direction
    !! @param[in] dy cell size in y direction
    !! @return depth rate
    function d_rate(this, t, dx, dy) result(rate)
        ! declarations
        class(point_source), intent(in):: this
        real(kind=8), intent(in):: t, dx, dy
        real(kind=8):: rate
        integer(kind=4):: i

        ! code
        i = count(this%t<t) + 1

        if (i .ge. this%nt) then
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
    !! @param[in] this an instance of point_source
    !! @param[in] mx the number of cells in x direction.
    !! @param[in] my the number of cells in y direction.
    !! @param[in] xlower the lower limit in x direction of this mesh.
    !! @param[in] ylower the lower limit in y direction of this mesh.
    !! @param[in] dx the cell size in x direction.
    !! @param[in] dy the cell size in y direction.
    !! @param[out] i the cell index in x direction.
    !! @param[out] j the cell index in y direction.
    subroutine cell_id(this, mx, my, xlower, ylower, dx, dy, i, j)
        ! declaration
        class(point_source), intent(in):: this
        integer(kind=4), intent(in):: mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy
        integer(kind=4), intent(out):: i, j
        real(kind=8):: xupper, yupper

        ! code
        i = -999
        j = -999

        xupper = xlower + mx * dx
        yupper = ylower + my * dy

        if ((this%coord(1) .ge. xlower) .and. (this%coord(1) .lt. xupper)) then
            if ((this%coord(2) .ge. ylower) .and. (this%coord(2) .lt. yupper)) then
                i = int(this%coord(1)/dx) + 1
                j = int(this%coord(2)/dy) + 1
            endif
        endif

    end subroutine cell_id

end module point_source_module
