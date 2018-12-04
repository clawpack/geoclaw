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
        !> @brief The index in aux array for hydro cell indicator
        integer(kind=4):: hydro_index = -1
        !> @brief Number of features in the collection.
        integer(kind=4):: nfeats = -1
        !> @brief Array of HydroFeature objects.
        type(HydroFeature), allocatable, dimension(:):: feats

        contains ! member functions
        !> @brief Initialization
        procedure:: init
        !> @brief Set hydro cell indicator index in aux
        procedure, private:: set_index
        !> @brief Underlying output driver.
        procedure, private:: write_data
        !> @brief Underlying input driver.
        procedure, private:: read_data
        !> @brief Return if a cell is a hydro cell.
        procedure:: is_hydro_cell
        !> @brief Update aux.
        procedure:: update_aux
        !> @brief Remove working fluid from hydro cells.
        procedure:: remove_fluid
        !> @brief Get a copy of the number of features.
        procedure:: get_n_features
        !> @brief Get a copy of the aux index.
        procedure:: get_hydro_index
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

        ! set index in aux
        call this%set_index()
        
    end subroutine init

    ! implementation of set_index
    subroutine set_index(this)
        use:: geoclaw_module, only: coordinate_system
        use:: friction_module, only: variable_friction
        use:: storm_module, only: wind_forcing, pressure_forcing
        use:: multilayer_module, only: num_layers
        use:: amr_module, only: auxtype
        class(HydroFeatureCollection), intent(inout):: this

        if (this%nfeats == -1) then
            print *, "Error: Hydro feature collection class is not yet initialized."
            stop
        else if (this%nfeats == 0) then
            this%hydro_index = -1
            return
        endif

        ! topo data
        this%hydro_index = 1

        ! coordinate system
        if (coordinate_system == 2) this%hydro_index = this%hydro_index + 2

        ! friction coefficients
        if (variable_friction) this%hydro_index = this%hydro_index + 1

        ! from wind module
        if (wind_forcing) this%hydro_index = this%hydro_index + 2

        ! from wind pressure
        if (pressure_forcing) this%hydro_index = this%hydro_index + 1

        ! from multi-layer
        if (num_layers > 1) this%hydro_index = this%hydro_index + num_layers

        ! finally, the index of hydro cell indicators
        this%hydro_index = this%hydro_index + 1

        ! check is the size of aux array is correct
        if (size(auxtype) /= this%hydro_index) then
            print *, "The index for hydro cell indicator,", this%hydro_index, &
                " is not correct. Maybe the number of aux is  not set &
                correctly in setrun.py"
            stop
        endif
    end subroutine set_index

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

    ! implementation of update_aux
    subroutine update_aux(this, mbc, mx, my, xlow, ylow, dx, dy, maux, aux)
        ! arguments
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4), intent(in):: mbc, mx, my, maux
        real(kind=8), intent(in):: xlow, ylow, dx, dy
        real(kind=8), intent(inout):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! local variables
        integer(kind=4):: i, j
        logical:: hydro_cell
        real(kind=8):: xl, xr, yb, yt

        ! if there's no hydrolic feature, exit the subroutine
        if (this%nfeats == 0) return

        xl = xlow - mbc * dy
        yb = ylow - mbc * dy

        do j=1-mbc, my+mbc
            yt = ylow + j * dy
            do i=1-mbc, mx+mbc
                xr = xlow + i * dx
                hydro_cell = this%is_hydro_cell(xl, xr, yb, yt)
                if (hydro_cell) aux(this%hydro_index, i, j) = 2D0
                xl = xr
            enddo
            yb = yt
        enddo
    end subroutine update_aux

    ! implementation of remove_fluid
    subroutine remove_fluid(this, meqn, mbc, mx, my, q, maux, aux)
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4), intent(in):: meqn, mbc, mx, my, maux
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! if there's no hydrolic feature, exit the subroutine
        if (this%nfeats == 0) return

        where(abs(aux(this%hydro_index, :, :)-2D0) < 1e-6) 
            q(1, :, :) = 0D0
            q(2, :, :) = 0D0
            q(3, :, :) = 0D0
        end where
    end subroutine remove_fluid

    ! implementation of get_n_features
    function get_n_features(this)
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4):: get_n_features
        get_n_features = this%nfeats 
    end function get_n_features

    ! implementation of get_hydro_index
    function get_hydro_index(this)
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4):: get_hydro_index
        get_hydro_index = this%hydro_index 
    end function get_hydro_index

end module hydro_feature_collection_module
