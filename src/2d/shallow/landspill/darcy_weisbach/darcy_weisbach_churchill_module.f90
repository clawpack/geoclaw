!
! darcy_weisbach_churchill_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Calculate coefficients with Churchill's model.
module darcy_weisbach_churchill_module
    use darcy_weisbach_abstract_module
    implicit none
    private
    public:: DarcyWeisbachChurchill

    !> @brief Darcy-Weisbach with cell-wide coefficients
    type, extends(DarcyWeisbachBase):: DarcyWeisbachChurchill
        private ! variables
        !> @brief Keep the underlying roughness file name.
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
        !> @brief Default roughness for regions uncovered by the file
        real(kind=8):: default_roughness = 0D0
        !> @brief Coefficient array.
        real(kind=8), allocatable, dimension(:, :):: roughness

        contains ! member functions
        !> @brief Initialize.
        procedure:: init_from_funit => init_from_funit_churchill
        !> @brief Underlying outputing.
        procedure:: write_data => write_data_churchill
        !> @brief Getting the coefficient of a single cell.
        procedure:: get_coefficient => get_coefficient_churchill
        !> @brief Read real coefficient file (in Esri ASCII format)
        procedure:: read_roughness_file
        !> @bried Destructor
        final:: destructor_churchill
    end type DarcyWeisbachChurchill

contains

    ! implementation of init_from_funit_churchill
    subroutine init_from_funit_churchill(this, funit, kin_vis)
        class(DarcyWeisbachChurchill), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: kin_vis

        this%name = "Churchill Darcy-Weisbach"
        this%nu = kin_vis
        read(funit, *) this%friction_tol
        read(funit, *) this%dry_tol
        read(funit, *) this%filename
        read(funit, *) this%default_roughness
        call this%read_roughness_file()
    end subroutine init_from_funit_churchill

    ! implementation of write_data_churchill
    subroutine write_data_churchill(this, iounit, iotype, v_list, stat, msg)
        class(DarcyWeisbachChurchill), intent(in):: this
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
            n, this%default_roughness, t, t, "=: default_coefficien", t, t, &
            "# (coefficient for cells uncovered by the file)"

    end subroutine write_data_churchill

    ! implementation of get_coefficient_churchill
    function get_coefficient_churchill(this, x, y, q) result(coef)
        class(DarcyWeisbachChurchill), intent(in):: this
        real(kind=8), intent(in):: x, y, q(3)
        real(kind=8):: coef, roughness, Re, theta1, theta2
        integer(kind=4):: i, j

        ! *********************************************************************
        ! Note: the following code is based on Re being defined by hydraulic
        !       radius, not hydraulic diameter
        ! *********************************************************************

        ! calculate local Reynolds number (defined by hydraulic radius, i.e., h)
        Re = dsqrt(q(2)**2+q(3)**2) / this%nu ! depth h is included in q(2) & q(3)

        ! initialize roughness with default value
        roughness = this%default_roughness

        ! when the coordinate is covered by the provided roughness file
        if ((x >= this%xlower) .and. (x < this%xupper)) then
            if ((y >= this%ylower) .and. (y < this%yupper)) then
                ! TODO: should we use at least linear interpolation?

                i = int((x-this%xlower)/this%cellsize) + 1
                j = int((y-this%ylower)/this%cellsize) + 1
                roughness = this%roughness(i, j)
            endif
        endif

        ! Churchill's model (note we use hydraulic radius, not hydraulic diameter)
        theta1 = (1.75D0 / Re)**0.9 + 0.0675 * roughness / q(1)
        theta1 = 2.457D0 * dlog(1D0/theta1)
        theta1 = theta1**16

        theta2 = (9382.5D0 / Re)**16

        coef = (theta1 + theta2)**3
        coef = dsqrt(coef)
        coef = 1D0 /coef
        coef = (coef + (2D0 / Re)**12)**(1D0/12D0)
        coef = 8D0 * coef

    end function get_coefficient_churchill

    ! implementation of read_roughness_file
    subroutine read_roughness_file(this, filename)
        use utility_module, only: parse_values
        class(DarcyWeisbachChurchill), intent(inout):: this
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
        allocate(this%roughness(this%mx, this%my))
        
        ! read coefficients
        do j = this%my, 1, -1
            read(funit, *) this%roughness(:, j)
        enddo

        ! handle missing data
        if (any(this%roughness == this%nodatavalue)) then
            where(this%roughness == this%nodatavalue) this%roughness = 0D0

            write(*, *) "WARNING: missing data found in the roughness file. &
                Set these data to zero automatically."
        endif

        close(funit)

        ! TODO: NetCDF file
    end subroutine read_roughness_file

    ! implementation of destructor_churchill
    subroutine destructor_churchill(this)
        type(DarcyWeisbachChurchill), intent(inout):: this
        this%name = ''
        this%filename = ''
        this%mx = 0
        this%my = 0
        this%xlower = 0D0
        this%ylower = 0D0
        this%xupper = 0D0
        this%yupper = 0D0
        this%cellsize = 0D0
        this%nodatavalue = 0D0
        this%default_roughness = 0D0

        if (allocated(this%roughness)) deallocate(this%roughness)
    end subroutine destructor_churchill

end module darcy_weisbach_churchill_module
