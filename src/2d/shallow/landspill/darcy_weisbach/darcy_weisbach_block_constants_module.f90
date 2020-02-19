!
! darcy_weisbach_constant_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Derived class from abstract Darcy-Weisbach class for block constants.
module darcy_weisbach_block_constants_module
    use darcy_weisbach_abstract_module
    implicit none
    private
    public:: DarcyWeisbachBlockConstants

    !> @brief Darcy-Weisbach with block constant coefficients.
    type, extends(DarcyWeisbachBase):: DarcyWeisbachBlockConstants
        private ! variables
        !> @brief Number of blocks
        integer(kind=4):: n_blocks = 0
        !> @brief Default coefficient for regions uncovered by the blocks
        real(kind=8):: default_coefficient = 0D0
        !> @brief Coefficients of blocks
        real(kind=8), allocatable, dimension(:):: coefficients
        !> @brief xlower of each block
        real(kind=8), allocatable, dimension(:):: xlowers
        !> @brief xupper of each block
        real(kind=8), allocatable, dimension(:):: xuppers
        !> @brief ylower of each block
        real(kind=8), allocatable, dimension(:):: ylowers
        !> @brief yupper of each block
        real(kind=8), allocatable, dimension(:):: yuppers

        contains ! member functions
        !> @brief Initialize.
        procedure:: init_from_funit => init_from_funit_block_constants
        !> @brief Underlying outputing.
        procedure:: write_data => write_data_block_constants
        !> @brief Getting the coefficient of a single cell.
        procedure:: get_coefficient => get_coefficient_block_constants
        !> @bried Destructor
        final:: destructor_block_constants
    end type DarcyWeisbachBlockConstants

contains

    ! implementation of init_from_funit_block_constant
    subroutine init_from_funit_block_constants(this, funit, kin_vis)
        class(DarcyWeisbachBlockConstants), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: kin_vis

        this%name = "Block Constant Darcy-Weisbach"
        this%nu = kin_vis
        read(funit, *) this%friction_tol
        read(funit, *) this%dry_tol
        read(funit, *) this%default_coefficient
        read(funit, *) this%n_blocks

        allocate(this%coefficients(this%n_blocks))
        allocate(this%xlowers(this%n_blocks), this%xuppers(this%n_blocks))
        allocate(this%ylowers(this%n_blocks), this%yuppers(this%n_blocks))

        read(funit, *) this%xlowers
        read(funit, *) this%xuppers
        read(funit, *) this%ylowers
        read(funit, *) this%yuppers
        read(funit, *) this%coefficients
    end subroutine init_from_funit_block_constants

    ! implementation of write_data_block_constants
    subroutine write_data_block_constants(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbachBlockConstants), intent(in):: this
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
            n, this%default_coefficient, t, t, "=: ", "default_coefficientl", &
            n, this%n_blocks, t, t, "=: n_blocks", t, t, "# (number of blocks)", &
            n, this%xlowers, t, t, "=: xlowers", &
            n, this%xuppers, t, t, "=: xuppers", &
            n, this%ylowers, t, t, "=: ylowers", &
            n, this%yuppers, t, t, "=: yuppers", &
            n, this%coefficients, t, t, "=: coefficients"

    end subroutine write_data_block_constants

    ! implementation of get_coefficient_block_constants
    function get_coefficient_block_constants(this, x, y, q) result(coef)
        class(DarcyWeisbachBlockConstants), intent(in):: this
        real(kind=8), intent(in):: x, y, q(3)
        real(kind=8):: coef
        logical:: found = .false.
        integer(kind=4):: i

        ! TODO: think about a better algorithm
        do i = 0, this%n_blocks
            if ((x >= this%xlowers(i)) .and. (x < this%xuppers(i))) then
                if ((y >= this%ylowers(i)) .and. (y < this%yuppers(i))) then
                    coef = this%coefficients(i)
                    found = .true.
                    exit ! exit the loop
                endif
            endif
        enddo

        if (.not. found) coef = this%default_coefficient

    end function get_coefficient_block_constants

    ! implementation of destructor_block_constants
    subroutine destructor_block_constants(this)
        type(DarcyWeisbachBlockConstants), intent(inout):: this

        this%n_blocks = 0
        if (allocated(this%coefficients)) deallocate(this%coefficients)
        if (allocated(this%xlowers)) deallocate(this%xlowers)
        if (allocated(this%xuppers)) deallocate(this%xuppers)
        if (allocated(this%ylowers)) deallocate(this%ylowers)
        if (allocated(this%yuppers)) deallocate(this%yuppers)
        this%name = ''
    end subroutine destructor_block_constants

end module darcy_weisbach_block_constants_module
