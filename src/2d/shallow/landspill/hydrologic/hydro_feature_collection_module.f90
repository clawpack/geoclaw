!
! hydro_feature_collection_module.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
!

module hydro_feature_collection_module
    use hydro_feature_module
    implicit none
    private
    public:: HydroFeatureCollection

    !> brief A class for collection of hydro features
    type:: HydroFeatureCollection
        private
        !> @brief Number of features in the collection.
        integer(kind=4):: nfeats = -1
        !> @brief Array of HydroFeature objects.
        type(HydroFeature), allocatable, dimension(:):: feats

        contains ! member functions
        !> @brief Initialization
        procedure:: init
        !> @brief Underlying output driver.
        procedure, private:: write_data
        !> @brief Underlying input driver.
        procedure, private:: read_data
        !> @brief Return if a cell is a hydro cell.
        procedure:: is_hydro_cell
        !> @brief Overload intrinsic write.
        generic:: write(formatted) => write_data
        !> @brief Overload intrinsic write.
        generic:: read(formatted) => read_data
        !> @brief Destructor
        final:: destructor
    end type HydroFeatureCollection


    !> @brief C++ style constructor.
    interface HydroFeatureCollection
        procedure:: constructor
    end interface HydroFeatureCollection
contains
    
    ! implementation of init
    subroutine init(this, filename)
        use geoclaw_module, only: coordinate_system

        ! arguments
        class(HydroFeatureCollection), intent(inout):: this
        character(len=*), intent(in), optional:: filename

        ! local variables
        integer(kind=4), parameter:: funit=252
        integer(kind=4):: i
        character(len=255):: feat_file

        ! so far, we only support xy coordinates (Cartesian)
        if (coordinate_system .ne. 1) then
            print *, "Hydrological functionality now only works with &
                Cartesian coordinates."
            stop
        endif

        ! open data file
        if (present(filename)) then
            call opendatafile(funit, filename)
        else
            call opendatafile(funit, "hydro_feature.data")
        endif

        ! read the number of hydro feature files
        read(funit, *) this%nfeats

        ! if no features, exit
        if (this%nfeats == 0) then
            close(funit)
            return
        endif

        ! allocating
        allocate(this%feats(this%nfeats))

        ! initialize each HydroFeature object
        do i = 1, this%nfeats
            read(funit, *) feat_file
            call this%feats(i)%init(trim(feat_file))
        enddo

        ! close the file
        close(funit)
    end subroutine init

    ! implementation of C++ stype constructor
    function constructor(filename)
        character(len=*), intent(in), optional:: filename
        type(HydroFeatureCollection):: constructor

        if (present(filename)) then
            call constructor%init(filename)
        else
            call constructor%init()
        endif
    end function constructor

    ! implementation destructor
    subroutine destructor(this)
        type(HydroFeatureCollection), intent(inout):: this
        this%nfeats = -1
        if (allocated(this%feats)) deallocate(this%feats)
    end subroutine destructor

    ! implementation of write_data
    subroutine write_data(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i
        character:: n, t

        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) n
        write(iounit, *, iostat=stat, iomsg=msg) &
            n, this%nfeats, t, t, t, t, "=: n_files # Number of hydro files"

        do i = 1, this%nfeats
            write(iounit, *, iostat=stat, iomsg=msg) n
            write(iounit, "(DT)", iostat=stat, iomsg=msg) this%feats(i)
        enddo

    end subroutine write_data

    ! implementation of read_data
    subroutine read_data(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(HydroFeatureCollection), intent(inout):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        print *, "Direct read of this object is prohibited!"
        stop

    end subroutine read_data

    ! implementation of is_hydro_cell
    function is_hydro_cell(this, x_cell_lower, x_cell_higher, &
        y_cell_lower, y_cell_higher)
        
        ! function argument
        class(HydroFeatureCollection), intent(in):: this
        real(kind=8), intent(in):: x_cell_lower, x_cell_higher 
        real(kind=8), intent(in):: y_cell_lower, y_cell_higher 
        logical:: is_hydro_cell

        ! local variables
        integer(kind=4):: i

        do i = 1, this%nfeats
            is_hydro_cell = this%feats(i)%is_hydro_cell(&
                x_cell_lower, x_cell_higher, y_cell_lower, y_cell_higher)

            if (is_hydro_cell) exit
        enddo

    end function is_hydro_cell
end module hydro_feature_collection_module
