! =====================================================
subroutine fixdry(meqn, mbc, mx, my, q, maux, aux)
! =====================================================
!
! Boussinesq replacement for the GeoClaw library fixdry, see
! $CLAW/geoclaw/src/2d/shallow/fixdry.f90 for the override contract.
!
! Same rule as the single-layer routine, extended to the Boussinesq correction.
! For the SGN system meqn = 5, where q(4:5) holds psi, the correction applied as
! q(2:3) = q(2:3) + dt * q(4:5) (see implicit_update_bouss_2Calls.f90).  psi is a
! momentum rate rather than a conserved quantity, so it goes to zero along with
! the momenta it corrects; the sweep over 2:meqn covers both, and reduces to the
! single-layer rule when the build is run with meqn = 3.
!
! aux is unused here; it is in the signature to match the library routine.

    use geoclaw_module, only: dry_tolerance

    implicit none

    ! Arguments
    integer, intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Locals
    integer :: i, j

    do j = 1-mbc, my+mbc
        do i = 1-mbc, mx+mbc
            if (q(1,i,j) < dry_tolerance) then
                q(1,i,j) = max(q(1,i,j), 0.d0)
                q(2:meqn,i,j) = 0.d0
            end if
        end do
    end do

end subroutine fixdry
