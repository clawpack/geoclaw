!> @file point_source.f90
!! @brief Point source class. 
!! @author Pi-Yueh Chuang
!! @version alpha
!! @date 2018-09-12

!> @brief Module for type declaration and implementation of point_source class.
!!
!! Only the point_source class is exposed to public. All others are private.
!!
module point_source_mod
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
        !> @brief Calculate and return depth rate.
        procedure:: d_rate
        !> @brief Write the information and data of this instance.
        procedure:: writef
        !> @brief Overriding intrinsic write function.
        generic:: write(formatted) => writef
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
        this%coord = coord
        this%nt = nt
        this%t = t ! implicit allocation
        this%v_rate = v_rate ! implicit allocation

        ! TODO: make sure the t array is sorted and v_rate is sorted based on t

    end subroutine init

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

        ! code
        write(iounit, *, iostat=stat, iomsg=msg) this%coord
        write(iounit, "(/)", iostat=stat, iomsg=msg)
        write(iounit, *, iostat=stat, iomsg=msg) this%nt
        write(iounit, "(/)", iostat=stat, iomsg=msg)
        write(iounit, *, iostat=stat, iomsg=msg) this%t
        write(iounit, "(/)", iostat=stat, iomsg=msg)
        write(iounit, *, iostat=stat, iomsg=msg) this%v_rate

    end subroutine writef

    !> @brief Finalize a point_source instance.
    !!
    !! Deallocate dynamic arrays and finalize the instance.
    !!
    !! @param[in] this an instance of point_source
    subroutine destructor(this)
        ! declaration
        type(point_source), intent(inout):: this

        ! code
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%v_rate)) deallocate(this%v_rate)

        this%coord = (/ 0D0, 0D0 /)
        this%nt = -1

    end subroutine destructor

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

end module point_source_mod
