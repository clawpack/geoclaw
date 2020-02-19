!
! darcy_weisbach_constant_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Derived class from abstract Darcy-Weisbach class for single constant.
module darcy_weisbach_constant_module
    use darcy_weisbach_abstract_module
    implicit none
    private
    public:: DarcyWeisbachConstant

    !> @brief Darcy-Weisbach with a constant coefficient.
    type, extends(DarcyWeisbachBase):: DarcyWeisbachConstant
        private ! variables
        !> @brief The constant coefficient.
        real(kind=8):: coefficient

        contains ! member functions
        !> @brief Initialize.
        procedure:: init_from_funit => init_from_funit_constant
        !> @brief Underlying outputing.
        procedure:: write_data => write_data_constant
        !> @brief Getting the coefficient of a single cell.
        procedure:: get_coefficient => get_coefficient_constant
        !> @bried Destructor
        final:: destructor_constant
    end type DarcyWeisbachConstant

contains

    ! implementation of init_from_funit_constant
    subroutine init_from_funit_constant(this, funit, kin_vis)
        class(DarcyWeisbachConstant), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: kin_vis

        this%name = "Constant Darcy-Weisbach"
        this%nu = kin_vis
        read(funit, *) this%friction_tol
        read(funit, *) this%dry_tol
        read(funit, *) this%coefficient
    end subroutine init_from_funit_constant

    ! implementation of write_data_constant
    subroutine write_data_constant(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbachConstant), intent(in):: this
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
            n, this%dry_tol, t, t, "=: ", "dry_tol", &
            n, this%coefficient, t, t, "=: ", "coefficient"
    end subroutine write_data_constant

    ! implementation of get_coefficient_constant
    function get_coefficient_constant(this, x, y, q) result(coef)
        class(DarcyWeisbachConstant), intent(in):: this
        real(kind=8), intent(in):: x, y, q(3)
        real(kind=8):: coef

        coef = this%coefficient
    end function get_coefficient_constant

    ! implementation of destructor_constant
    subroutine destructor_constant(this)
        type(DarcyWeisbachConstant), intent(inout):: this
        this%coefficient = -999
        this%name = ''
    end subroutine destructor_constant

end module darcy_weisbach_constant_module
