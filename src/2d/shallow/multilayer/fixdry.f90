! =====================================================
subroutine fixdry(meqn, mbc, mx, my, q, maux, aux)
! =====================================================
!
! Multilayer replacement for the GeoClaw library fixdry, see
! $CLAW/geoclaw/src/2d/shallow/fixdry.f90 for the override contract.
!
! Same rule as the single-layer routine, applied per layer: clamp a negative
! depth up to zero and zero that layer's momenta.  Note that q stores rho*h
! rather than h, so the dry test divides by rho(k) before comparing against
! that layer's tolerance.
!
! aux is unused here; it is in the signature to match the library routine.

    use geoclaw_module, only: rho
    use multilayer_module, only: num_layers, dry_tolerance

    implicit none

    ! Arguments
    integer, intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Locals
    integer :: i, j, k

    do j = 1-mbc, my+mbc
        do i = 1-mbc, mx+mbc
            do k = 1, num_layers
                if (q(3*(k-1)+1,i,j) / rho(k) < dry_tolerance(k)) then
                    q(3*(k-1)+1,i,j) = max(q(3*(k-1)+1,i,j), 0.d0)
                    q(3*(k-1)+2,i,j) = 0.d0
                    q(3*(k-1)+3,i,j) = 0.d0
                end if
            end do
        end do
    end do

end subroutine fixdry
