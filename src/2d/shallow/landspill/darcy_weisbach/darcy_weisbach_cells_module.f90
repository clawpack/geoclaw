!
! darcy_weisbach_constant_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Derived class from abstract Darcy-Weisbach class for single constant.
module darcy_weisbach_cells_module
    use darcy_weisbach_abstract_module
    implicit none
    private
    public:: DarcyWeisbachCells

    !> @brief Darcy-Weisbach with cell-wide coefficients
    type, extends(DarcyWeisbachBase):: DarcyWeisbachCells
        private ! variables
        !> @brief Keep the underlying coefficient file name.
        character(len=255):: filename
        !> @brief Number of cells in x direction
        integer(kind=4):: mx
        !> @brief Number of cells in y direction
        integer(kind=4):: my
        !> @brief xlower of the coefficient grid provided.
        real(kind=8):: xlower
        !> @brief ylower of the coefficient grid provided.
        real(kind=8):: ylower
        !> @brief xupper of the coefficient grid provided.
        real(kind=8):: xupper
        !> @brief yupper of the coefficient grid provided.
        real(kind=8):: yupper
        !> @brief Cell size (assume uniform in both directions).
        real(kind=8):: cellsize
        !> @brief No data value.
        real(kind=8):: nodatavalue
        !> @brief Default coefficient for regions uncovered by the file
        real(kind=8):: default_coefficient = 0D0
        !> @brief Coefficient array.
        real(kind=8), allocatable, dimension(:, :):: coefficients

        contains ! member functions
        !> @brief Initialize.
        procedure:: init_from_funit => init_from_funit_cells
        !> @brief Underlying outputing.
        procedure:: write_data => write_data_cells
        !> @brief Getting the coefficient of a single cell.
        procedure:: get_coefficient => get_coefficient_cells
        !> @brief Read real coefficient file (in Esri ASCII format)
        procedure:: read_coefficient_file
        !> @bried Destructor
        final:: destructor_cells
    end type DarcyWeisbachCells

contains

    ! implementation of init_from_funit_cells
    subroutine init_from_funit_cells(this, funit, kin_vis)
        class(DarcyWeisbachCells), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: kin_vis

        this%name = "Cell-wide Darcy-Weisbach"
        this%nu = kin_vis
        read(funit, *) this%friction_tol
        read(funit, *) this%dry_tol
        read(funit, *) this%filename
        read(funit, *) this%default_coefficient
        call this%read_coefficient_file()
    end subroutine init_from_funit_cells

    ! implementation of write_data_cells
    subroutine write_data_cells(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbachCells), intent(in):: this
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
            n, this%filename, t, t, "=: filename", t, t, "# (coefficient file)", &
            n, this%default_coefficient, t, t, "=: default_coefficien", t, t, &
            "# (coefficient for cells uncovered by the file)"

    end subroutine write_data_cells

    ! implementation of get_coefficient_cells
    function get_coefficient_cells(this, x, y, q) result(coef)
        class(DarcyWeisbachCells), intent(in):: this
        real(kind=8), intent(in):: x, y, q(3)
        real(kind=8):: coef
        integer(kind=4):: i, j

        ! when the coordinate is covered by the provided coefficient file
        if ((x >= this%xlower) .and. (x < this%xupper)) then
            if ((y >= this%ylower) .and. (y < this%yupper)) then
                i = int((x-this%xlower)/this%cellsize) + 1
                j = int((y-this%ylower)/this%cellsize) + 1
                coef = this%coefficients(i, j)
                ! TODO: should we use at least linear interpolation?

                return ! function return earlier
            endif
        endif

        ! for the case that the cell is not covered
        coef = this%default_coefficient

    end function get_coefficient_cells

    ! implementation of read_coefficient_file
    subroutine read_coefficient_file(this, filename)
        use utility_module, only: parse_values
        class(DarcyWeisbachCells), intent(inout):: this
        character(len=*), intent(in), optional:: filename
        integer(kind=4), parameter:: funit = 253
        character(len=255):: line
        integer(kind=4):: n_values, i, j
        real(kind=8):: values(10)

        if (present(filename)) then
            open(unit=funit, file=filename, action="read", &
                status="old", form="formatted")
        else
            open(unit=funit, file=this%filename, action="read", &
                status="old", form="formatted")
        endif

        ! mx
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%mx = idnint(values(1))

        ! my
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%my = idnint(values(1))

        ! xlower
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%xlower = values(1)

        ! ylower
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%ylower = values(1)

        ! cellsize
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%cellsize = values(1)

        ! nodatavalue
        read(funit, "(A)") line
        call parse_values(line, n_values, values)
        this%nodatavalue = values(1)

        ! calculate xupper and yupper
        this%xupper = this%xlower + this%cellsize * this%mx
        this%yupper = this%ylower + this%cellsize * this%my

        ! allocate coefficients
        allocate(this%coefficients(this%mx, this%my))
        
        ! read coefficients
        do j = 1, this%my
            read(funit, *) this%coefficients(:, j)
        enddo

        ! handle missing data
        if (any(this%coefficients == this%nodatavalue)) then
            where(this%coefficients == this%nodatavalue) this%coefficients = 0D0

            write(*, *) "WARNING: missing data found in the Darcy-Weisbach &
                coefficient file. Set these data to zero automatically."
        endif

        close(funit)

        ! TODO: NetCDF file
    end subroutine read_coefficient_file

    ! implementation of destructor_cells
    subroutine destructor_cells(this)
        type(DarcyWeisbachCells), intent(inout):: this
        this%name = ""
        this%filename = ""
        this%mx = 0
        this%my = 0
        this%xlower = 0D0
        this%ylower = 0D0
        this%xupper = 0D0
        this%yupper = 0D0
        this%cellsize = 0D0
        this%nodatavalue = 0D0
        this%default_coefficient = 0D0

        if (allocated(this%coefficients)) deallocate(this%coefficients)
    end subroutine destructor_cells

end module darcy_weisbach_cells_module
