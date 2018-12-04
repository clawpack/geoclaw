!
! hydro_feature_module.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
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
        logical, allocatable, dimension(:, :):: hydro_cell
        !> @brief Number of the cells laying on the outlines of the features
        integer(kind=4):: n_ol_cells
        !> @brief x-coordinates of the outline cell
        real(kind=8), allocatable, dimension(:):: x_ol_cells
        !> @brief y-coordinates of the outline cell
        real(kind=8), allocatable, dimension(:):: y_ol_cells
        !> @brief Accumulated removed fluid on the outline at each AMR level
        real(kind=8), allocatable, dimension(:, :):: rmvd_ol_cells

        contains ! member function
        !> @brief Init from file
        procedure:: init
        !> @brief Read hydro data file (in Esri ASCII format)
        procedure:: read_file
        !> @brief Underlying outputing.
        procedure:: write_data
        !> @brief Return if a cell is a hydro cell.
        procedure:: is_hydro_cell
        !> @brief Identify cells located on the outlines of the features
        procedure, private:: identify_outline
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
        call this%identify_outline()
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
        allocate(this%hydro_cell(this%mx, this%my))
        this%hydro_cell = .false.

        ! allocate temporary array for real/integer numbers of every line
        allocate(temp(this%mx))
        
        ! read coefficients
        do j = this%my, 1, -1
            read(funit, *) temp
            where(temp .ne. this%nodatavalue) this%hydro_cell(:, j) = .true.
        enddo

        close(funit)
        deallocate(temp)

        ! TODO: NetCDF file
    end subroutine read_file

    ! implementation of identify_outline
    subroutine identify_outline(this)
        use:: amr_module, only: mxnest
        
        ! function argument
        class(HydroFeature), intent(inout):: this

        ! definitions of point. Used only in this subroutine
        type:: point
            real(kind=8):: x
            real(kind=8):: y
            type(point), pointer:: prev => null()
            type(point), pointer:: next => null()
        end type point

        ! local variables
        type(point), pointer:: bg => null()
        type(point), pointer:: ed => null()
        type(point), pointer:: head => null()
        integer(kind=4):: i, j
        real(kind=8):: xlower_shift
        real(kind=8):: ylower_shift

        ! code
        xlower_shift = this%xlower - this%cellsize / 2D0
        ylower_shift = this%ylower - this%cellsize / 2D0
        this%n_ol_cells = 0

        ! create a linked list for temporarily recording outlines
        do i = 1, this%mx
            do j = 1, this%my
                ! if this is not a hydro cell, skip
                if (.not. this%hydro_cell(i, j)) cycle
                
                ! if all 3x3 cells are hydro cells, this is not on outline
                if (all(this%hydro_cell(max(1, i-1):min(this%mx, i+1), &
                    max(1, j-1):min(this%my, j+1)))) cycle

                ! if this is the first point found
                if (.not. associated(bg)) then
                    allocate(bg)
                    head => bg
                else ! else, append to the list
                    allocate(head%next)
                    head%next%prev => head
                    head => head%next
                end if

                head%x = i * this%cellsize + xlower_shift
                head%y = j * this%cellsize + ylower_shift
                this%n_ol_cells = this%n_ol_cells + 1
            end do
        end do

        ! the last element of the list
        ed => head
        nullify(head)

        allocate(this%x_ol_cells(this%n_ol_cells))
        allocate(this%y_ol_cells(this%n_ol_cells))
        allocate(this%rmvd_ol_cells(mxnest, this%n_ol_cells))

        ! init/copy values
        this%rmvd_ol_cells = 0D0
        
        head => bg
        do i = 1, this%n_ol_cells
            this%x_ol_cells(i) = head%x
            this%y_ol_cells(i) = head%y
            head => head%next
        end do
        nullify(head)

        ! deallocate linked list
        head => ed%prev
        do while(associated(head))
            deallocate(head%next)
            head => head%prev
        end do
        nullify(head)
        deallocate(bg)

    end subroutine identify_outline

    ! implementation of is_hydro_cell
    function is_hydro_cell(this, x_cell_lower, x_cell_higher, &
        y_cell_lower, y_cell_higher)
        
        ! function argument
        class(HydroFeature), intent(in):: this
        real(kind=8), intent(in):: x_cell_lower, x_cell_higher 
        real(kind=8), intent(in):: y_cell_lower, y_cell_higher 
        logical:: is_hydro_cell

        ! local variables
        integer(kind=4):: i_lower, i_higher
        integer(kind=4):: j_lower, j_higher

        ! initialize is_hydro_cell
        is_hydro_cell = .false.

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
            i_higher = int((x_cell_lower - this%xlower) / this%cellsize) + 1
        endif

        if (y_cell_lower < this%ylower) then
            j_lower = 1
        else
            j_lower = int((y_cell_lower - this%ylower) / this%cellsize) + 1
        endif

        if (y_cell_higher >= this%yupper) then
            j_higher = this%my
        else
            j_higher = int((y_cell_lower - this%ylower) / this%cellsize) + 1
        endif

        ! see if any cell covered in the block is a hydro cell
        if (any(this%hydro_cell(i_lower:i_higher, j_lower:j_higher))) then
            is_hydro_cell = .true.
        endif

    end function is_hydro_cell

end module hydro_feature_module
