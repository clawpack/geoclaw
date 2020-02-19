!
! darcy_weisbach_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!

!> @brief The top-level class wrapper for different types of Darcy-Weisbach.
module darcy_weisbach_module
    use darcy_weisbach_abstract_module
    use darcy_weisbach_constant_module
    use darcy_weisbach_block_constants_module
    use darcy_weisbach_cells_module
    use darcy_weisbach_three_regimes_module
    use darcy_weisbach_churchill_module
    use darcy_weisbach_two_regimes_module
    implicit none
    private
    public:: DarcyWeisbach

    type:: DarcyWeisbach
        ! member variables
        private
        integer(kind=4):: type = -1
        class(DarcyWeisbachBase), pointer:: ptr => null()

        ! member functions
        contains
        procedure:: init
        !> @brief Outputing.
        procedure:: write_data
        !> @brief Inputting.
        procedure:: read_data
        !> @brief Get the coefficient of a single cell.
        procedure:: get_coefficient
        !> @brief Add forcing terms to a grid.
        procedure:: apply_to_grid
        !> @brief An interface to get type.
        procedure:: get_type
        !> @brief Overriding generic write
        generic:: write(formatted) => write_data
        !> @brief Overriding generic read
        generic:: read(formatted) => read_data
        !> @brief Destructor.
        final:: destructor
    end type DarcyWeisbach

    interface DarcyWeisbach
        procedure constructor
    end interface DarcyWeisbach

contains

    ! implementation of constrtuctor
    function constructor(kin_vis, filename) 
        type(DarcyWeisbach):: constructor
        character(len=*), intent(in), optional:: filename
        real(kind=8), intent(in):: kin_vis

        if (present(filename)) then
            call constructor%init(kin_vis, filename)
        else
            call constructor%init(kin_vis)
        endif
    end function constructor

    ! implementation of init
    subroutine init(this, kin_vis, filename)
        class(DarcyWeisbach), intent(inout):: this
        character(len=*), intent(in), optional:: filename
        real(kind=8), intent(in):: kin_vis
        integer(kind=4), parameter:: funit = 254

        if (present(filename)) then
            call opendatafile(funit, filename)
        else
            call opendatafile(funit, "darcy_weisbach.data")
        endif

        ! in case this instance has been initialized
        call destructor(this)

        read(funit, *) this%type

        select case (this%type)
        case (0) ! a trivial null object
            allocate(DarcyWeisbachNull::this%ptr)
        case (1)
            allocate(DarcyWeisbachConstant::this%ptr)
        case (2)
            allocate(DarcyWeisbachBlockConstants::this%ptr)
        case (3)
            allocate(DarcyWeisbachCells::this%ptr)
        case (4)
            allocate(DarcyWeisbachThreeRegimes::this%ptr)
        case (5)
            allocate(DarcyWeisbachChurchill::this%ptr)
        case (6)
            allocate(DarcyWeisbachTwoRegimes::this%ptr)
        case default
            write(*, *) "Invalid Darcy-Weisbach friction type."
            stop
        end select

        ! TODO: maybe design a dummy DarcyWeisbach object for type=0 case?

        call this%ptr%init_from_funit(funit, kin_vis)

        close(funit)
    end subroutine init

    ! implementation of destructor
    subroutine destructor(this)
        type(DarcyWeisbach), intent(inout):: this

        if (associated(this%ptr)) deallocate(this%ptr)

        this%type = -1
    end subroutine destructor

    ! implementation of write_data
    subroutine write_data(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbach), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        ! write type
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%type, t, t, "=: type of Darcy-Weisbach factor"

        ! use underlying object's method
        write(iounit, *, iostat=stat, iomsg=msg) this%ptr
    end subroutine write_data

    ! implementation of read_data
    subroutine read_data(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbach), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        ! use underlying object's method
        read(iounit, *, iostat=stat, iomsg=msg) this%ptr
    end subroutine read_data

    ! implementation fo get_coefficient
    function get_coefficient(this, x, y, q) result(coef)
    class(DarcyWeisbach), intent(in):: this
    real(kind=8), intent(in):: x, y, q(3)
    real(kind=8):: coef
        ! use underlying object's method
        coef = this%ptr%get_coefficient(x, y, q)
    end function get_coefficient

    ! implementation of apply_to_grid
    subroutine apply_to_grid(this, meqn, mbc, mx, my, &
        xlower, ylower, dx, dy, q, dt)
        ! passed variables
        class(DarcyWeisbach), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! local variables
        integer(kind=4)::i, j
        real(kind=4):: dgamma

        ! use underlying object's method
        call this%ptr%apply_to_grid(meqn, mbc, mx, my, &
            xlower, ylower, dx, dy, q, dt)
    end subroutine apply_to_grid

    ! implementation of get_type
    function get_type(this)
        class(DarcyWeisbach), intent(in):: this
        integer(kind=4):: get_type
        get_type = this%type
    end function get_type
end module darcy_weisbach_module
