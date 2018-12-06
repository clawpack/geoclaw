!
! evap_module.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
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
        procedure:: init => init_evapmodel
        !> @brief Remove evaporated fluid.
        procedure:: apply_to_grid => apply_to_grid_evapmodel
        !> @brief Return a copy of the evaporation type.
        procedure:: get_type => get_type_evapmodel
        !> @brief Return percentage remained.
        procedure:: remained_percentage
        !> @brief Destructor.
        final:: destructor_evapmodel
    end type EvapModel

contains

    ! init_evapmodel
    subroutine init_evapmodel(this, T, filename)
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
        call destructor_evapmodel(this)

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
    end subroutine init_evapmodel

    ! apply_to_grid_evapmodel
    subroutine apply_to_grid_evapmodel(this, meqn, mbc, mx, my, xlower, &
                             ylower, dx, dy, q, maux, aux, t, dt)
        use geoclaw_module, only: dry_tolerance
        implicit none
        class(EvapModel), intent(inout):: this
        integer(kind=4), intent(in):: meqn, mbc, mx,my, maux
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        call this%ptr%apply_to_grid(meqn, mbc, mx, my, xlower, &
                                    ylower, dx, dy, q, maux, aux, t, dt)

    end subroutine apply_to_grid_evapmodel

    ! get_type_evapmodel
    function get_type_evapmodel(this) result(ans)
        class(EvapModel), intent(in):: this
        integer(kind=4):: ans
        ans = this%type
    end function get_type_evapmodel

    ! remained_percentage
    function remained_percentage(this, t, dt) result(ans)
        class(EvapModel), intent(in):: this
        real(kind=8), intent(in):: t, dt
        real(kind=8):: ans

        ans = this%ptr%remained_kernel(t, dt)
    end function remained_percentage

    ! destructor_evapmodel
    subroutine destructor_evapmodel(this)
        type(EvapModel), intent(inout):: this
        if (associated(this%ptr)) deallocate(this%ptr)
        this%type = -1
    end subroutine destructor_evapmodel

end module evap_module
