!
! hydro_feature_collection_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!

module hydro_feature_collection_module
    use hydro_feature_module
    use SPM_module
    implicit none
    private
    public:: HydroFeatureCollection

    !> @brief A class for collection of hydro features
    type:: HydroFeatureCollection
        private
        !> @brief The index in aux array for hydro cell indicator
        integer(kind=4):: hydro_index = -1
        !> @brief Number of features in the collection.
        integer(kind=4):: nfeats = -1
        !> @brief Array of HydroFeature objects.
        type(HydroFeature), allocatable, dimension(:):: feats
        !> @brief A removed-fluid tracer represented by CSR.
        type(MultiLayerCSR):: tracer
        !> @brief xlower of the mesh tracing removed fluid.
        real(kind=8):: tracer_xlower
        !> @brief ylower of the mesh tracing removed fluid.
        real(kind=8):: tracer_ylower
        !> @brief xupper of the mesh tracing removed fluid.
        real(kind=8):: tracer_xupper
        !> @brief yupper of the mesh tracing removed fluid.
        real(kind=8):: tracer_yupper
        !> @brief Number of cell in x-dir of the mesh tracing removed fluid.
        integer(kind=4):: tracer_mx
        !> @brief Number of cell in y-dir of the mesh tracing removed fluid.
        integer(kind=4):: tracer_my
        !> @brief Cell size in x-dir of the mesh tracing removed fluid.
        real(kind=8):: tracer_dx
        !> @brief Cell size in x-dir of the mesh tracing removed fluid.
        real(kind=8):: tracer_dy

        contains ! member functions
        !> @brief Initialization
        procedure:: init
        !> @brief Set hydro cell indicator index in aux
        procedure, private:: set_aux_index
        !> @brief Initialize tracer.
        procedure, private:: init_rmvd_fluid_tracer
        !> @brief Underlying output driver.
        procedure, private:: write_data
        !> @brief Underlying input driver.
        procedure, private:: read_data
        !> @brief Return the id of the hydro feature owning this cell.
        procedure:: cell_type
        !> @brief Update aux.
        procedure:: update_aux
        !> @brief Remove working fluid from hydro cells.
        procedure:: remove_fluid
        !> @brief Evaporate fluid removed.
        procedure:: evap_fluid
        !> @brief Get a copy of the number of features.
        procedure:: get_n_features
        !> @brief Get a copy of the aux index.
        procedure:: get_hydro_index
        !> @brief Output data of removed fluid.
        procedure:: output_removed_fluid
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

        ! so far, we only support xy coordinates (Cartesian)
        if (coordinate_system .ne. 1) then
            print *, "Hydrological functionality now only works with &
                Cartesian coordinates."
            stop
        endif

        ! set index in aux
        call this%set_aux_index()

        ! set tracer tracing removed fluid on feature boundary
        call this%init_rmvd_fluid_tracer()
        
    end subroutine init

    ! implementation of set_aux_index
    subroutine set_aux_index(this)
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
    end subroutine set_aux_index

    ! implementation init_rmvd_fluid_tracer
    subroutine init_rmvd_fluid_tracer(this)
        use:: amr_module, only: xlower, ylower, xupper, yupper
        use:: amr_module, only: mxnest, hxposs, hyposs, intratx, intraty
        class(HydroFeatureCollection), intent(inout):: this

        integer(kind=4):: i, j
        integer(kind=1):: tracer_cell_type
        real(kind=8):: xl, xh, yl, yh
        type(COO):: temp_mtx


        this%tracer_xlower = xlower
        this%tracer_ylower = ylower
        this%tracer_xupper = xupper
        this%tracer_yupper = yupper

        this%tracer_dx = hxposs(1)
        this%tracer_dy = hyposs(1)

        do i = 1, mxnest-1
            this%tracer_dx = this%tracer_dx / intratx(i)
            this%tracer_dy = this%tracer_dy / intraty(i)
        end do

        this%tracer_mx = idnint((xupper-xlower)/this%tracer_dx)
        this%tracer_my = idnint((yupper-ylower)/this%tracer_dy)

        ! initialize temporary COO matrix
        call temp_mtx%init(this%tracer_mx, this%tracer_my)

        ! loop
        do j = 1, this%tracer_my

            yl = (j - 1) * this%tracer_dy + this%tracer_ylower
            yh = yl + this%tracer_dy

            do i = 1, this%tracer_mx

                xl = (i - 1) * this%tracer_dx + this%tracer_xlower
                xh = xl + this%tracer_dx

                tracer_cell_type = this%cell_type(xl, xh, yl, yh)

                if (tracer_cell_type == 2) then
                    call temp_mtx%append(IndexSet(i, j, 0D0))
                end if
            end do
        end do

        ! transform to CSR
        call this%tracer%init(3, temp_mtx)

        ! let the time recorder be -1 at the beginning
        call this%tracer%set_all(1, -1D0)

        ! destroy temporary COO matrix
        call temp_mtx%destroy()

    end subroutine init_rmvd_fluid_tracer

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

    ! implementation of cell_type
    function cell_type(this, x_cell_lower, x_cell_higher, &
                       y_cell_lower, y_cell_higher)
        
        ! function argument
        class(HydroFeatureCollection), intent(in):: this
        real(kind=8), intent(in):: x_cell_lower, x_cell_higher 
        real(kind=8), intent(in):: y_cell_lower, y_cell_higher 
        integer(kind=1):: cell_type

        ! local variables
        integer(kind=4):: i

        do i = 1, this%nfeats
            cell_type = this%feats(i)%cell_type(x_cell_lower, &
                x_cell_higher, y_cell_lower, y_cell_higher)

            ! no need to go through other hydro feature
            if (cell_type /= 0) exit
        enddo

    end function cell_type

    ! implementation of update_aux
    subroutine update_aux(this, mbc, mx, my, xlow, ylow, dx, dy, maux, aux)
        use:: amr_module, only: xlower, ylower

        ! arguments
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4), intent(in):: mbc, mx, my, maux
        real(kind=8), intent(in):: xlow, ylow, dx, dy
        real(kind=8), intent(inout):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! local variables
        integer(kind=4):: i, j, ilo, jlo
        logical:: hydro_cell
        real(kind=8):: xl, xr, yb, yt

        ! if there's no hydrolic feature, exit the subroutine
        if (this%nfeats == 0) return

        ilo = floor((xlow-xlower+.05d0*dx)/dx)
        jlo = floor((ylow-ylower+.05d0*dy)/dy)

        do j=1-mbc, my+mbc
            yb = ylower + real(jlo+j-1, 8) * dy
            yt = yb + dy
            do i=1-mbc, mx+mbc
                xl = xlower + real(ilo+i-1, 8) * dx
                xr = xl + dx
                aux(this%hydro_index, i, j) = this%cell_type(xl, xr, yb, yt)
            enddo
        enddo
    end subroutine update_aux

    ! implementation of remove_fluid
    subroutine remove_fluid(this, level, meqn, mbc, mx, my, xlow, &
                            ylow, dx, dy, q, maux, aux, time)
        use:: amr_module, only: mxnest

        ! input arguments
        class(HydroFeatureCollection), intent(inout):: this
        integer(kind=4), intent(in):: level, meqn, mbc, mx, my, maux
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: xlow, ylow, dx, dy, time

        ! local variables
        integer(kind=4):: i, j, nnz, k
        integer(kind=4):: rowl, rowh, coll, colh
        integer(kind=4):: idxs((int(dx/this%tracer_dx)+2)*(int(dy/this%tracer_dy)+2))
        real(kind=8):: val

        ! if there's no hydrolic feature, exit the subroutine
        if (this%nfeats == 0) return

        ! not the finest grid, just remove fluids
        if (level /= mxnest) then
            do j = 1-mbc, my+mbc
                do i = 1-mbc, mx+mbc
                    if (idnint(aux(this%hydro_index, i, j)) /= 0) q(:, i, j) = 0D0
                end do
            end do

        ! the finest grid, remove fluids and also record the volume
        else
            do j = 1-mbc, my+mbc
                do i = 1-mbc, mx+mbc
                    ! skip, if this is not a boundary of hydor features
                    if (idnint(aux(this%hydro_index, i, j)) == 0) cycle

                    ! only boundary cells with non-zero depth undergoes this
                    if (q(1, i, j) /= 0D0) then

                        rowl = int(((i - 1) * dx + xlow - &
                            this%tracer_xlower) / this%tracer_dx) + 1
                        rowh = int((i * dx + xlow - &
                            this%tracer_xlower) / this%tracer_dx) + 1
                        coll = int(((j - 1) * dy + ylow - &
                            this%tracer_ylower) / this%tracer_dy) + 1
                        colh = int((j * dy + ylow - &
                            this%tracer_ylower) / this%tracer_dy) + 1

                        call this%tracer%get_indices(rowl, rowh, coll, colh, nnz, idxs)
                        val = q(1, i, j) * dx * dy / nnz

                        call this%tracer%add_multiples(nnz, idxs, 2, val)
                        call this%tracer%add_multiples(nnz, idxs, 3, val)

                        ! record the time of first contact
                        ! TODO: this may be somehow inefficient!!
                        do k = 1, nnz
                            if (this%tracer%get_value(idxs(k), 1) < 0D0) &
                                call this%tracer%set(idxs(k), 1, time)
                        end do
                    end if
                    
                    ! zero the cell quantities
                    q(:, i, j) = 0D0
                end do
            end do
        end if

    end subroutine remove_fluid

    ! implementation of evap_fluid
    subroutine evap_fluid(this, remained_rate, tracker)
        class(HydroFeatureCollection), intent(inout):: this
        real(kind=8), intent(in):: remained_rate
        real(kind=8), intent(inout), optional:: tracker
        real(kind=8):: temp
        
        if (this%nfeats == 0) return

        if (present(tracker)) then ! we need to track the volume evaporated
            tracker = tracker + this%tracer%sum(3)
            call this%tracer%mult_all(3, remained_rate)
            tracker = tracker - this%tracer%sum(3)
        else
            call this%tracer%mult_all(3, remained_rate)
        end if
    end subroutine evap_fluid

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

    ! implementation of output_removed_fluid
    subroutine output_removed_fluid(this)
        class(HydroFeatureCollection), intent(in):: this
        integer(kind=4):: i, j, x_prec, y_prec, idx
        real(kind=8):: val(3), x, y
        character(len=255):: fmt_str

        ! calculate the required precision for coordinates
        if (this%tracer_dx >= 1D0) then
            x_prec = 1
        else
            x_prec = iabs(floor(dlog10(this%tracer_dx))) + 1
        end if

        if (this%tracer_dy >= 1D0) then
            y_prec = 1
        else
            y_prec = iabs(floor(dlog10(this%tracer_dy))) + 1
        end if

        ! prepare the format for each line of output
        write(fmt_str, "('(F0.', I0, ', '', '', F0.', &
            I0, ', '', '', F0.1, 2('', '', ES23.15E3))')") x_prec, y_prec

        ! write to file
        open(unit=995, file="removed_fluid.csv", action="write", recl=9999)
        do i = 1, this%tracer_mx
            do j = 1, this%tracer_my

                idx = this%tracer%get_index(i, j)
                if (idx == 0) cycle ! not a non-zero element

                val(2) = this%tracer%get_value(idx, 2)

                if (val(2) /= 0D0) then
                    x = (i - 0.5) * this%tracer_dx + this%tracer_xlower
                    y = (j - 0.5) * this%tracer_dy + this%tracer_ylower
                    val(1) = this%tracer%get_value(idx, 1)
                    val(3) = this%tracer%get_value(idx, 3)

                    write(995, fmt_str) x, y, val
                end if
            end do
        end do
        close(995)
    end subroutine output_removed_fluid

end module hydro_feature_collection_module
