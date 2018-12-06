!
! fingas1996.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
!

!> @brief Implementation of Fingas 1996 models for non-volatile fluids.
module fingas1996_module
    use:: evap_base_module
    implicit none
    private
    public:: EvapFingas1996Log, EvapFingas1996SQRT

    !> @brief A model from Fingas 1996 based on natural log law.
    type, extends(EvapBase):: EvapFingas1996Log
        private
        !> @brief Final and combined coefficients used in Fingas 1996 model.
        real(kind=8):: final_coeffs

        contains
        !> @brief Initialization with an opened file unit.
        procedure:: init_with_funit => init_with_funit_fingas1996log
        !> @brief The kernel calculating remained percentage from t to t+dt.
        procedure:: remained_kernel => remained_kernel_fingas1996log
        !> @brief Destructor.
        final:: destructor_fingas1996log
    end type EvapFingas1996Log

    !> @brief A model from Fingas 1996 based on square-root law.
    type, extends(EvapBase):: EvapFingas1996SQRT
        private
        !> @brief Final and combined coefficients used in Fingas 1996 model.
        real(kind=8):: final_coeffs

        contains
        !> @brief Initialization with an opened file unit.
        procedure:: init_with_funit => init_with_funit_fingas1996sqrt
        !> @brief The kernel calculating remained percentage from t to t+dt.
        procedure:: remained_kernel => remained_kernel_fingas1996sqrt
        !> @brief Destructor.
        final:: destructor_fingas1996sqrt
    end type EvapFingas1996SQRT

contains

    ! init_with_funit_fingas1996log
    subroutine init_with_funit_fingas1996log(this, funit, T)
        class(EvapFingas1996Log), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: T

        real(kind=8):: C1, C2
        integer(kind=4):: n_coeffs

        this%model_name = "Fingas1996 Log"
        this%ambient_temperature = T
        this%evap_volume_tracker = 0D0

        read(funit, *) n_coeffs

        if (n_coeffs /= 2) then
            print *, "The number of coefficients in Fingas' model should be 2."
            stop
        end if

        read(funit, *) C1
        read(funit, *) C2

        this%final_coeffs = C1 + C2 * this%ambient_temperature
    end subroutine init_with_funit_fingas1996log

    ! remained_kernel_fingas1996log
    function remained_kernel_fingas1996log(this, t, dt) result(remained_percent)
        class(EvapFingas1996Log), intent(in):: this
        real(kind=8), intent(in):: t
        real(kind=8), intent(in):: dt
        real(kind=8):: remained_percent

        remained_percent = &
            (1D2 - this%final_coeffs * dlog((t+dt)/6D1)) / &
            (1D2 - this%final_coeffs * dlog(t/6D1))
    end function remained_kernel_fingas1996log

    ! destructor_fingas1996log
    subroutine destructor_fingas1996log(this)
        type(EvapFingas1996Log), intent(inout):: this
        this%model_name = "none"
        this%ambient_temperature = 0D0
        this%evap_volume_tracker = 0D0
        this%final_coeffs = 0D0
    end subroutine destructor_fingas1996log

    ! init_with_funit_fingas1996sqrt
    subroutine init_with_funit_fingas1996sqrt(this, funit, T)
        class(EvapFingas1996SQRT), intent(inout):: this
        integer(kind=4), intent(in):: funit
        real(kind=8), intent(in):: T

        real(kind=8):: C1, C2
        integer(kind=4):: n_coeffs

        this%model_name = "Fingas1996 Log"
        this%ambient_temperature = T
        this%evap_volume_tracker = 0D0

        if (n_coeffs /= 2) then
            print *, "The number of coefficients in Fingas' model should be 2."
            stop
        end if

        read(funit, *) C1
        read(funit, *) C2

        this%final_coeffs = C1 + C2 * this%ambient_temperature
    end subroutine init_with_funit_fingas1996sqrt

    ! remained_kernel_fingas1996sqrt
    function remained_kernel_fingas1996sqrt(this, t, dt) result(remained_percent)
        class(EvapFingas1996SQRT), intent(in):: this
        real(kind=8), intent(in):: t
        real(kind=8), intent(in):: dt
        real(kind=8):: remained_percent

        remained_percent = &
            (1D2 - this%final_coeffs * dsqrt((t+dt)/6D1)) / &
            (1D2 - this%final_coeffs * dsqrt(t/6D1))
    end function remained_kernel_fingas1996sqrt

    ! destructor_fingas1996sqrt
    subroutine destructor_fingas1996sqrt(this)
        type(EvapFingas1996SQRT), intent(inout):: this
        this%model_name = "none"
        this%ambient_temperature = 0D0
        this%evap_volume_tracker = 0D0
        this%final_coeffs = 0D0
    end subroutine destructor_fingas1996sqrt

end module fingas1996_module
