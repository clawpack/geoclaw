!
! evap_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!

!> @brief The top-level evaporation module.
module evap_module
    use evap_base_module
    use fingas1996_module
    implicit none
    private
    public:: EvapModel

    !> @brief A top-level class for evaporation models.
    type:: EvapModel
        private
        !> @brief The type of underlying evaporation model.
        integer(kind=4):: type = -1
        !> @brief The pointer pointing to the real model.
        class(EvapBase), pointer:: ptr => null()

        contains
        !> @brief Initialization.
        procedure:: init
        !> @brief Remove evaporated fluid.
        procedure:: apply_to_grid
        !> @brief Return a copy of the evaporation type.
        procedure:: get_type
        !> @brief Reset the releasing profile of a point source.
        procedure:: reset_release_profile
        !> @brief Evaporate the fluid in hydrological features.
        procedure:: evap_hydro_fluid
        !> @brief Return percentage remained.
        procedure:: remained_percentage
        !> @brief Output evaporated volume.
        procedure:: output
        !> @brief Destructor.
        final:: destructor
    end type EvapModel

contains

    ! init
    subroutine init(this, T, filename)
        class(EvapModel), intent(inout):: this
        real(kind=8), intent(in):: T
        character(len=*), intent(in), optional:: filename

        ! local variables
        integer(kind=4), parameter:: funit = 250

        ! open datafile
        if (present(filename)) then
            call opendatafile(funit, filename)
        else
            call opendatafile(funit, "evaporation.data")
        endif

        ! in case this instance has been initialized
        call destructor(this)

        ! read type of evaporation model
        read(funit, *) this%type

        ! create actual underlying evaporation model
        select case (this%type)
        case (0) ! a trivial null object
            allocate(EvapNull::this%ptr)
        case (1)
            allocate(EvapFingas1996Log::this%ptr)
        case (2)
            allocate(EvapFingas1996SQRT::this%ptr)
        case default
            print *, "Invalid evaporation type."
            stop
        end select

        ! initialize the underlying model
        call this%ptr%init_with_funit(funit, T)

        close(funit)
    end subroutine init

    ! apply_to_grid
    subroutine apply_to_grid(this, meqn, mbc, mx, my, xlower, &
                             ylower, dx, dy, q, maux, aux, t, dt)
        implicit none
        class(EvapModel), intent(inout):: this
        integer(kind=4), intent(in):: meqn, mbc, mx,my, maux
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        call this%ptr%apply_to_grid(meqn, mbc, mx, my, xlower, &
                                    ylower, dx, dy, q, maux, aux, t, dt)

    end subroutine apply_to_grid

    ! get_type
    function get_type(this) result(ans)
        class(EvapModel), intent(in):: this
        integer(kind=4):: ans
        ans = this%type
    end function get_type

    ! remained_percentage
    function remained_percentage(this, t, dt) result(ans)
        class(EvapModel), intent(in):: this
        real(kind=8), intent(in):: t, dt
        real(kind=8):: ans

        ans = this%ptr%remained_kernel(t, dt)
    end function remained_percentage

    ! destructor
    subroutine destructor(this)
        type(EvapModel), intent(inout):: this
        if (associated(this%ptr)) deallocate(this%ptr)
        this%type = -1
    end subroutine destructor

    ! modify the time of releasing profile
    subroutine reset_release_profile(this, pts, tol)
        use point_source_collection_module
        class(EvapModel), intent(in):: this
        type(PointSourceCollection), intent(inout):: pts
        real(kind=8), intent(in), optional:: tol

        ! local variables
        integer(kind=4):: i
        integer(kind=4):: npts ! number of point sources
        integer(kind=4):: nt ! number of stages in release profile
        real(kind=8):: Vp ! remained volume in pipeline
        real(kind=8):: dt, T, V0
        real(kind=8), allocatable, dimension(:):: rates, times

        if (this%type == 0) return

        if (present(tol)) then
            dt = tol
        else
            dt = 1D-3
        end if

        npts = pts%get_n_points()
        if (npts > 1) then
            print *, "ERROR: Evaporation model currently does not support &
                multiple point sources."
            stop
        end if

        nt = pts%get_n_stages(1) ! number of stages of 1st point
        allocate(times(nt), rates(nt))
        call pts%get_times(1, times)
        call pts%get_v_rates(1, rates)

        Vp = rates(1) * times(1)
        do i = 2, nt
            Vp = Vp + rates(i) * (times(i) - times(i-1))
        end do
        V0 = Vp

        ! find the time that no fuild remained in pipeline
        do while (Vp > 0D0)
            i = count(times < T) + 1 ! current stage
            Vp = Vp - rates(i) * dt ! minus the amount release to the field
            Vp = Vp * this%remained_percentage(T, dt)
            T = T + dt
        end do

        i = count(times < T) + 1
        times(i) = T
        rates(i+1:) = 0D0

        call pts%set_times(1, times)
        call pts%set_v_rates(1, rates)

        ! now Vp is the initial volume in the pipe AFTER modifying times
        Vp = rates(1) * times(1)
        do i = 2, nt
            Vp = Vp + rates(i) * (times(i) - times(i-1))
        end do
        V0 = V0 - Vp ! now V0 is the different between original and modified volume

        Vp = this%ptr%get_evaporated_volume() ! now Vp: current evaporated volume
        Vp = Vp + V0 ! now Vp: add the volume evaporated from the pipe
        call this%ptr%set_evaporated_volume(Vp) ! set the evaporated volume back

        deallocate(times, rates)
    end subroutine reset_release_profile

    ! evaporate the fluid already in hydrological features
    subroutine evap_hydro_fluid(this, hydro, t, dt)
        use hydro_feature_collection_module
        class(EvapModel), intent(in):: this
        type(HydroFeatureCollection), intent(inout):: hydro
        real(kind=8), intent(in):: t, dt

        call this%ptr%evap_hydro_fluid(hydro, t, dt)
    end subroutine evap_hydro_fluid

    ! output
    subroutine output(this)
        class(EvapModel), intent(in):: this

        open(unit=895, file="evaporated_fluid.dat", action="write")
        write(895, *) this%ptr%get_evaporated_volume()
        close(895)
    end subroutine output

end module evap_module
