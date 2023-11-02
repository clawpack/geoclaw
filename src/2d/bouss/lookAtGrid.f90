!
! -------------------------------------------------------------
!
subroutine lookAtGrid(val,valOther,dx,dy,xlowWithGhost,ylowWithGhost,nvar,mitot,mjtot,  &
                      ng,iu,doUpdate,lev)

   implicit none
   real(kind=8), intent(in) :: val(nvar,mitot,mjtot)
   real(kind=8), intent(in) :: valOther(nvar,mitot,mjtot)
   real(kind=8), intent(in) :: dx,dy,xlowWithGhost,ylowWithGhost
   integer, intent(in) :: nvar,mitot,mjtot,ng,iu,lev
   logical, intent(in) :: doUpdate
   real(kind=8) :: speed(mitot,mjtot),speedmax
   integer :: imax,jmax

   integer :: i,j,ivar
   real(kind=8) :: x,y

   speedmax = 0.d0
   imax = -1
   jmax = -1
   ! if ng=0 look at ghost cells. If ng = 2 are not looking at ghost cells.
   do i = ng+1, mitot-ng
      x = xlowWithGhost + (i-0.5d0)*dx
   do j = ng+1, mjtot-ng
      y = ylowWithGhost + (j-0.5d0)*dy
      if (val(1,i,j) .eq. 0.d0) cycle
      speed(i,j) = sqrt(val(2,i,j)**2+val(3,i,j)**2)/val(1,i,j)
      if (speed(i,j) .gt. speedmax) then
         speedmax = speed(i,j)
         imax = i
         jmax = j
      endif
      if (doUpdate) then
        write(iu,101)x,y,i,j,(val(ivar,i,j),ivar=1,3),  &
                         (valOther(ivar,i,j),ivar=4,5),speed(i,j)
      else
        write(iu,101)x,y,i,j,(val(ivar,i,j),ivar=1,nvar),speed(i,j)
      endif
 101  format(2e14.7,2i5,6e14.6)
   end do
   end do

   !write(*,102) lev,speedmax,imax,jmax
 102  format("from lookAtGrid level ",i3," speedmax = ",e15.7," at cell ",2i5)


end subroutine lookAtGrid

