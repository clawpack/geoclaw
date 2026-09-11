 ! =====================================================
 subroutine fixdry(meqn, mbc, mx, my, q, maux, aux)
 ! =====================================================
 !
 ! Reset cells whose depth is below dry_tolerance to a consistent dry state:
 ! zero negative depths to zero and zero the momenta.  Operates over the full
 ! patch, ghost cells included.
 !
 ! This routine is an override point.  Applications and GeoClaw variants may
 ! shadow it by placing their own fixdry.f90 in the application SOURCES list
 ! (see the GeoClaw docs on replacing library routines).  Two constraints on
 ! any replacement:
 !   * It must stay a plain external subroutine, not a module procedure
 !   * It must be thread safe.  Both call sites run inside an !$OMP PARALLEL DO
 !     over grids

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