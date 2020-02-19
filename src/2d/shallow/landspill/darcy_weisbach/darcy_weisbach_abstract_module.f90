!
! darcy_weisbach_abstract_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Abstract and null class for Darcy-Weisbach classes.
module darcy_weisbach_abstract_module
    implicit none
    private
    public:: DarcyWeisbachBase, DarcyWeisbachNull

    !> @brief The abstract class for friction objects.
    type, abstract:: DarcyWeisbachBase
        private ! variables
        !> @brief The name/type of the actual class type of an instance.
        character(len=15), public:: name = "None"
        !> @brief The friction model will only apply to depth smaller to this value.
        real(kind=8), public:: friction_tol = 1D6
        !> @brief The tolerance for momentum to be considered zero
        real(kind=8), public:: dry_tol = 1D-30
        !> @brief Kinematic viscosity (m^2 / s).
        real(kind=8), public:: nu = 0D0

        contains ! member functions
        !> @brief Initialize.
        procedure(init_from_funit_base), deferred:: init_from_funit
        !> @brief Underlying output function.
        procedure(write_data_base), deferred:: write_data
        !> @brief The function to return the coefficient of a cell.
        procedure(get_coefficient_base), deferred:: get_coefficient
        !> @brief The function to add forcing term to a grid.
        procedure:: apply_to_grid
        !> @brief Underlying read function.
        procedure:: read_base
        !> @brief Overriding generic write.
        generic, public:: write(formatted) => write_data
        !> @brief Overriding generic read.
        generic, public:: read(formatted) => read_base
    end type DarcyWeisbachBase

    !> @brief Interface of abstract (i.e., virtual) functions.
    abstract interface
        ! init_from_funit_base
        subroutine init_from_funit_base(this, funit, kin_vis)
            import 
            implicit none
            class(DarcyWeisbachBase), intent(inout):: this
            integer(kind=4), intent(in):: funit
            real(kind=8), intent(in):: kin_vis
        end subroutine init_from_funit_base

        ! write_data_base
        subroutine write_data_base(this, iounit, iotype, v_list, stat, msg)
            import
            implicit none
            class(DarcyWeisbachBase), intent(in):: this
            integer(kind=4), intent(in):: iounit
            character(*), intent(in)::iotype
            integer(kind=4), intent(in):: v_list(:)
            integer(kind=4), intent(out):: stat
            character(*), intent(inout):: msg
        end subroutine write_data_base

        ! get_coefficient_base
        function get_coefficient_base(this, x, y, q) result(coef)
            import
            implicit none
            class(DarcyWeisbachBase), intent(in):: this
            real(kind=8), intent(in):: x, y
            real(kind=8), intent(in):: q(3)
            real(kind=8):: coef
        end function get_coefficient_base
    end interface

    type, extends(DarcyWeisbachBase):: DarcyWeisbachNull
        private

        contains
        !> @brief Initialize.
        procedure:: init_from_funit => init_from_funit_null
        !> @brief Underlying outputing.
        procedure:: write_data => write_data_null
        !> @brief Getting the coefficient of a single cell.
        procedure:: get_coefficient => get_coefficient_null
        !> @brief Trivial function.
        procedure:: apply_to_grid => apply_to_grid_null
        !> @bried Destructor
        final:: destructor_null
    end type DarcyWeisbachNull

contains

    !> @brief Kernel function for Darcy-Weisbach friction force.
    function kernel(q, coefficient) result(gamma)
        implicit none

        ! passed-in variables
        real(kind=8), intent(in):: q(3), coefficient
        real(kind=8):: gamma
        
        ! code
        gamma = 0.125D0 * coefficient * dsqrt(q(2)**2+q(3)**2) / (q(1)**2)

    end function kernel      

    !> @brief Calculate and add friction force to the whole array of q.
    subroutine apply_to_grid(this, meqn, mbc, mx, my, &
        xlower, ylower, dx, dy, q, dt)
        implicit none

        ! passed-in variables
        class(DarcyWeisbachBase), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! local variables
        integer(kind=4)::i, j
        real(kind=8):: x, y, dgamma
        
        ! code
        do j = 1, my
            y = ylower + j * dy - dy / 2D0 ! y coordinate
            do i = 1, mx
                x = xlower + i * dx - dx / 2D0 ! x coordinate
                ! skip this cell if depth is larger than the friction tolerance
                if (q(1, i, j) > this%friction_tol) cycle

                ! set momentum to be zero if depth is smaller than dry_tol
                if (q(1, i, j) < this%dry_tol) then
                    q(2:3, i, j) = 0D0
                    cycle
                endif

                dgamma = 1D0 + dt * kernel(q(:, i, j), &
                    this%get_coefficient(x, y, q(1:3, i, j)))

                q(2, i, j) = q(2, i, j) / dgamma
                q(3, i, j) = q(3, i, j) / dgamma
            enddo
        enddo
    end subroutine apply_to_grid      

    ! implementation of read_base
    subroutine read_base(this, iounit, iotype, v_list, stat, msg)
        implicit none

        ! variable declaration
        class(DarcyWeisbachBase), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        write(*, *) "Error: direct read of this object is prohibited."
        stop
    end subroutine read_base

    ! implementation of init_from_funit_null
    subroutine init_from_funit_null(this, funit, kin_vis)
        implicit none
        class(DarcyWeisbachNull), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: kin_vis

        this%name = "Trivial Darcy-Weisbach (Do nothing object)"

        ! do not read from file. Force to safety numbers
        this%friction_tol = 0D0
        this%dry_tol = 1D7
        this%nu = 0D0
    end subroutine init_from_funit_null

    ! implementation of write_data_null
    subroutine write_data_null(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbachNull), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%friction_tol, t, t, "=: ", "friction_tol", &
            n, this%dry_tol, t, t, "=: ", "dry_tol"
    end subroutine write_data_null

    ! implementation of get_coefficient_null
    function get_coefficient_null(this, x, y, q) result(coef)
        class(DarcyWeisbachNull), intent(in):: this
        real(kind=8), intent(in):: x, y, q(3)
        real(kind=8):: coef

        coef = 0D0
    end function get_coefficient_null

    ! implementation of apply_to_grid_null
    subroutine apply_to_grid_null(this, meqn, mbc, mx, my, &
        xlower, ylower, dx, dy, q, dt)

        ! passed-in variables
        class(DarcyWeisbachNull), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my
        real(kind=8), intent(in):: xlower, ylower, dx, dy, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        return
    end subroutine apply_to_grid_null

    ! implementation of destructor_null
    subroutine destructor_null(this)
        type(DarcyWeisbachNull), intent(inout):: this
        this%name = ''
    end subroutine destructor_null

end module darcy_weisbach_abstract_module
