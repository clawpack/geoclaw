! =====================================================
subroutine fixdry(meqn, mbc, mx, my, q, maux, aux)
! =====================================================
!
! Reset cells whose depth is below dry_tolerance to a consistent dry state:
! clamp a negative depth up to zero and zero the momenta.  Operates over the
! full patch, ghost cells included.
!
! Note what this routine does *not* do: it leaves every component past the
! momenta untouched.  Zeroing a conserved quantity at the dry state destroys
! conservation of that quantity, so a component beyond (h, hu, hv) is the
! application's business, not this routine's.  Where such a component does need
! handling, the right treatment is generally a positivity clamp mirroring the
! one applied to the depth, not a zero.
!
! This routine is an override point.  Applications and GeoClaw variants supply
! their own rule by placing a fixdry.f90 in the application SOURCES list (see
! the GeoClaw docs on replacing library routines); the multilayer and
! Boussinesq variants in this repository do exactly that.  Two constraints on
! any replacement:
!   * It must stay a plain external subroutine, not a module procedure
!   * It must be thread safe.  Every call site runs inside an !$OMP PARALLEL DO
!     over grids
!
! aux is unused here.  It is in the signature so that an override can reach
! topography and the capacity function without a signature change.

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
                q(2,i,j) = 0.d0
                q(3,i,j) = 0.d0
            end if
        end do
    end do

end subroutine fixdry
