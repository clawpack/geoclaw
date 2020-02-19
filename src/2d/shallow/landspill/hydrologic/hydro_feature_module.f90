!
! hydro_feature_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!


!> @brief Definition and implementation of HydroFeature
module hydro_feature_module
    implicit none
    private
    public:: HydroFeature

    !> @brief Hydro feature from a single file
    type:: HydroFeature
        private
        !> @brief The name of the underlying file name.
        character(len=:), allocatable:: filename
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
        !> @brief Array indicating whether a cell belongs to hydrologic heature.
        integer(kind=1), allocatable, dimension(:, :):: hydro_cell

        contains ! member function
        !> @brief Init from file
        procedure:: init
        !> @brief Change the value of hydro feature boundary
        procedure:: process_boundary
        !> @brief Read hydro data file (in Esri ASCII format)
        procedure, private:: read_file
        !> @brief Underlying outputing.
        procedure, private:: write_data
        !> @brief Return if a cell is a hydro cell.
        procedure:: cell_type
        !> @brief Overriding intrinsic write function.
        generic:: write(formatted) => write_data
        !> @bried Destructor
        final:: destructor
    end type HydroFeature

contains

    ! implementation of init_from_funit
    subroutine init(this, filename)
        class(HydroFeature), intent(inout):: this
        character(len=*), intent(in):: filename

        this%filename = filename
        call this%read_file()
        call this%process_boundary()
    end subroutine init

    ! implementation of destructor
    subroutine destructor(this)
        type(HydroFeature), intent(inout):: this
        this%filename = ''
        this%mx = 0
        this%my = 0
        this%xlower = 0D0
        this%ylower = 0D0
        this%xupper = 0D0
        this%yupper = 0D0
        this%cellsize = 0D0
        this%nodatavalue = 0D0

        if (allocated(this%hydro_cell)) deallocate(this%hydro_cell)
    end subroutine destructor

    ! implementation of write_data
    subroutine write_data(this, iounit, iotype, v_list, stat, msg)
        class(HydroFeature), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%filename, t, t, "=: filename"

    end subroutine write_data
    
    ! implementation of read_file
    subroutine read_file(this, filename)
        use utility_module, only: parse_values
        class(HydroFeature), intent(inout):: this
        character(len=*), intent(in), optional:: filename
        integer(kind=4), parameter:: funit = 251
        integer(kind=4):: n_values, i, j
        real(kind=8):: values(10)
        real(kind=8), allocatable, dimension(:):: temp
        character(len=255):: line

        if (present(filename)) then
            print *, "Reading hydro feature:", filename
            open(unit=funit, file=filename, action="read", &
                status="old", form="formatted")
        else
            print *, "Reading hydro feature:", this%filename
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
        allocate(this%hydro_cell(this%mx, this%my))
        this%hydro_cell = 0

        ! allocate temporary array for real/integer numbers of every line
        allocate(temp(this%mx))
        
        ! read coefficients
        do j = this%my, 1, -1
            read(funit, *) temp
            where(temp .ne. this%nodatavalue) this%hydro_cell(:, j) = 1
        enddo

        close(funit)
        deallocate(temp)

        ! TODO: NetCDF file
    end subroutine read_file

    ! implementation of process_boundary
    subroutine process_boundary(this)
        class(HydroFeature), intent(inout):: this

        integer(kind=4):: i, j

        ! we don't care the cells on computational domain boundary, really!
        do j = 2, this%my-1
            do i = 2, this%mx-1
                ! skip interation if this is not a hydro cell
                if (this%hydro_cell(i, j) == 0) cycle

                if (any(this%hydro_cell(i-1:i+1, j-1:j+1) == 0)) then
                    this%hydro_cell(i, j) = 2
                end if
            end do
        end do
    end subroutine process_boundary

    ! implementation of is_hydro_cell
    function cell_type(this, x_cell_lower, x_cell_higher, &
        y_cell_lower, y_cell_higher)
        
        ! function argument
        class(HydroFeature), intent(in):: this
        real(kind=8), intent(in):: x_cell_lower, x_cell_higher 
        real(kind=8), intent(in):: y_cell_lower, y_cell_higher 
        integer(kind=1):: cell_type

        ! local variables
        integer(kind=4):: i_lower, i_higher
        integer(kind=4):: j_lower, j_higher

        ! initialize cell_type
        cell_type = 0

        ! completely out of domain, return false
        if ((x_cell_higher <= this%xlower) .or. &
            (x_cell_lower >= this%xupper) .or. &
            (y_cell_higher <= this%ylower) .or. &
            (y_cell_lower >= this%yupper)) then

            return
        endif

        if (x_cell_lower < this%xlower) then
            i_lower = 1
        else
            i_lower = int((x_cell_lower - this%xlower) / this%cellsize) + 1
        endif

        if (x_cell_higher >= this%xupper) then
            i_higher = this%mx
        else
            i_higher = int((x_cell_higher - this%xlower) / this%cellsize) + 1
        endif

        if (y_cell_lower < this%ylower) then
            j_lower = 1
        else
            j_lower = int((y_cell_lower - this%ylower) / this%cellsize) + 1
        endif

        if (y_cell_higher >= this%yupper) then
            j_higher = this%my
        else
            j_higher = int((y_cell_higher - this%ylower) / this%cellsize) + 1
        endif

        ! see if any cell covered in the block is a hydro boundary cell
        if (any(this%hydro_cell(i_lower:i_higher, j_lower:j_higher) == 2)) then
            cell_type = 2
            return
        endif

        ! see if this is a interior hydro cell
        if (any(this%hydro_cell(i_lower:i_higher, j_lower:j_higher) == 1)) then
            cell_type = 1
        endif

    end function cell_type

end module hydro_feature_module
