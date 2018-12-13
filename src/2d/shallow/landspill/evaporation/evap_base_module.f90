!
! evap_base_module.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
!

!> @brief Evaporation abstract class module.
module evap_base_module
    implicit none
    private
    public:: EvapBase, EvapNull

    !> @brief Abstract class for all evaporation classes.
    type, abstract:: EvapBase
        private
        !> @brief The name of the model used.
        character(len=15), public:: model_name = "none"
        !> @brief The ambient_temperature.
        real(kind=8), public:: ambient_temperature = 0D0
        !> @brief A tracer to track evaporated volume.
        real(kind=8), public:: evap_volume_tracker = 0D0

        contains
        !> @brief Initialization with an opened file unit.
        procedure(init_funit_base), deferred:: init_with_funit
        !> @brief The kernel calculating remained percentage from t to t+dt.
        procedure(remained_kernel_base), deferred:: remained_kernel
        !> @brief Remove evaporated volume from a grid.
        procedure:: apply_to_grid
        !> @brief Evaporate the fluid in hydrological features.
        procedure:: evap_hydro_fluid
    end type EvapBase

    ! virtual functions/subroutines
    abstract interface

        ! init_funit_base
        subroutine init_funit_base(this, funit, T)
            import
            implicit none
            class(EvapBase), intent(inout):: this
            integer(kind=4), intent(in):: funit
            real(kind=8), intent(in):: T
        end subroutine init_funit_base

        ! remained_kernel_base
        function remained_kernel_base(this, t, dt) result(remained_percent)
            import
            implicit none
            class(EvapBase), intent(in):: this
            real(kind=8), intent(in):: t
            real(kind=8), intent(in):: dt
            real(kind=8):: remained_percent
        end function remained_kernel_base
    end interface

    !> @brief Do-nothing evaporation class.
    type, extends(EvapBase):: EvapNull
        private
        contains
        !> @brief Initialization with an opened file unit.
        procedure:: init_with_funit => init_with_funit_null
        !> @brief The kernel calculating remained percentage from t to t+dt.
        procedure:: remained_kernel => remained_kernel_null
        !> @brief A dummy function.
        procedure:: apply_to_grid => apply_to_grid_null
        !> @brief Evaporate the fluid in hydrological features.
        procedure:: evap_hydro_fluid => evap_hydro_fluid_null
    end type EvapNull

contains

    ! apply_to_grid
    subroutine apply_to_grid(this, meqn, mbc, mx, my, xlower, &
                             ylower, dx, dy, q, maux, aux, t, dt)
        use geoclaw_module, only: dry_tolerance
        implicit none
        class(EvapBase), intent(inout):: this
        integer(kind=4), intent(in):: meqn, mbc, mx,my, maux
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

        ! local variables
        real(kind=8):: remained, new_val
        integer(kind=4):: i, j

        ! get the percentage of evaporation between t to t+dt
        remained = this%remained_kernel(t, dt)

        do j = 1-mbc, my+mbc
            do i = 1-mbc, mx+mbc
                ! skip, because there's no fluid in this cell
                if (q(1, i, j) < 1e-9) cycle

                new_val = q(1, i, j) * remained

                ! update removed volume
                this%evap_volume_tracker = &
                    this%evap_volume_tracker + q(1, i, j) - new_val

                ! remove from grid
                q(1, i, j) = new_val

                ! lower than dry_tolerance, then fluid stops
                if (q(1, i, j) <= dry_tolerance) q(2:3, i, j) = 0D0

            end do
        end do
    end subroutine apply_to_grid

    ! evaporate the fluid already in hydrological features
    subroutine evap_hydro_fluid(this, hydro, t, dt)
        use hydro_feature_collection_module
        class(EvapBase), intent(in):: this
        type(HydroFeatureCollection), intent(inout):: hydro
        real(kind=8), intent(in):: t, dt
        real(kind=8):: remained_rate

        remained_rate = this%remained_kernel(t, dt)
        call hydro%evap_fluid(remained_rate)
    end subroutine evap_hydro_fluid

    ! init_with_funit_null
    subroutine init_with_funit_null(this, funit, T)
        class(EvapNull), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: T

        this%model_name = "Do-nothing evaporation."
        this%ambient_temperature = T
        this%evap_volume_tracker = 0D0

    end subroutine init_with_funit_null

    ! kernel_null
    function remained_kernel_null(this, t, dt) result(remained_percent)
        class(EvapNull), intent(in):: this
        real(kind=8), intent(in):: t
        real(kind=8), intent(in):: dt
        real(kind=8):: remained_percent

        remained_percent = 1D0
    end function remained_kernel_null

    ! apply_to_grid_null
    subroutine apply_to_grid_null(this, meqn, mbc, mx, my, xlower, &
                                  ylower, dx, dy, q, maux, aux, t, dt)
        use geoclaw_module, only: dry_tolerance
        implicit none
        class(EvapNull), intent(inout):: this
        integer(kind=4), intent(in):: meqn, mbc, mx,my, maux
        real(kind=8), intent(in):: xlower, ylower, dx, dy, t, dt
        real(kind=8), intent(inout):: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
        real(kind=8), intent(in):: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    end subroutine apply_to_grid_null

    ! evap_hydro_fluid_null
    subroutine evap_hydro_fluid_null(this, hydro, t, dt)
        use hydro_feature_collection_module
        class(EvapNull), intent(in):: this
        type(HydroFeatureCollection), intent(inout):: hydro
        real(kind=8), intent(in):: t, dt
    end subroutine evap_hydro_fluid_null

end module evap_base_module
