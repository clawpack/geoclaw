!
!  :::::::::::::: BOUND :::::::::::::::::::::::::::::::::::::::::::
!     This routine sets the boundary values for a given grid 
!     at level level.
!     We are setting the values for a strip ng zones wide all
!     the way around the border, in 4 rectangular strips.
!
!     Outputs from this routine:
!     The values around the border of the grid are inserted
!     directly into the enlarged valbig array.
!
!     This routine calls the routine filpatch
!     which for any block of mesh points on a given level,
!     intersects that block with all grids on that level and with
!     the physical boundaries, copies the values into the
!     appropriate intersecting regions, and interpolates the remaining
!     cells from coarser grids as required.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine bound(time,nvar,ng,valbig,mitot,mjtot,mptr,aux,naux)

  use amr_module, only: rnode, node, hxposs, hyposs, cornxlo, cornxhi
  use amr_module, only: cornylo, cornyhi, ndilo, ndihi, ndjlo, ndjhi
  use amr_module, only: nestlevel, xlower, xupper, ylower, yupper
  use amr_module, only: xperdom, yperdom, spheredom

  implicit none

  ! Input
  integer, intent(in) :: nvar, ng, mitot, mjtot, mptr, naux
  real(kind=8), intent(in) :: time
  real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
  real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)

  ! Locals
  integer :: ilo, ihi, jlo, jhi, level  ! rect(4)
  real(kind=8) :: xleft, xright, ybot, ytop, hx, hy, xl, xr, yb, yt
  real(kind=8) :: xloWithGhost,  xhiWithGhost,  yloWithGhost, yhiWithGhost
  logical      :: patchOnly

  xleft  = rnode(cornxlo, mptr)
  xright = rnode(cornxhi, mptr)
  ybot   = rnode(cornylo, mptr)
  ytop   = rnode(cornyhi, mptr)
  ilo    = node(ndilo, mptr)
  ihi    = node(ndihi, mptr)
  jlo    = node(ndjlo, mptr)
  jhi    = node(ndjhi, mptr)
  level  = node(nestlevel, mptr)
  hx     = hxposs(level)
  hy     = hyposs(level)

  xloWithGhost = xleft  - ng*hx
  xhiWithGhost = xright + ng*hx
  yloWithGhost =  ybot  - ng*hy
  yhiWithGhost =  ytop  + ng*hy
  ! used in filaptch for bc2amr: for patches it is called. for full grids called from bound below
  patchOnly = .false.  




  ! left boundary
  xl = xleft - ng*hx
  xr = xleft
  yb = ybot 
  yt = ytop
  if ((xl < xlower) .and. xperdom) then
     call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,ng+1, &
          ilo-ng,ilo-1,jlo,jhi,ilo-ng,ihi+ng,jlo-ng,jhi+ng,patchOnly)
  else
     call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,ng+1,ilo-ng,ilo-1,jlo,jhi,patchOnly)
  endif

  ! right boundary
  xl = xright
  xr = xright + ng*hx
  yb = ybot
  yt = ytop

  if ((xr .gt. xupper) .and. xperdom) then
     call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
          mitot-ng+1,ng+1,ihi+1,ihi+ng,jlo,jhi,ilo-ng,ihi+ng,jlo-ng,jhi+ng,patchOnly)
  else
     call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
          mitot-ng+1,ng+1,ihi+1,ihi+ng,jlo,jhi,patchOnly)
  endif

  ! bottom boundary
  xl = xleft  - ng*hx
  xr = xright + ng*hx
  yb = ybot - ng*hy
  yt = ybot

  if ( ((yb < ylower) .and. (yperdom .or. spheredom)) .or. &
       ( ((xl < xlower) .or. (xr > xupper)) .and. xperdom) ) then
     call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,   &
                      1,1,ilo-ng,ihi+ng,jlo-ng,jlo-1,                &
                      ilo-ng,ihi+ng,jlo-ng,jhi+ng,patchOnly)
  else
     call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,1,ilo-ng,ihi+ng,jlo-ng,jlo-1,patchOnly)
  endif


  ! top boundary
  xl = xleft - ng*hx
  xr = xright + ng*hx
  yb = ytop
  yt = ytop + ng*hy

  if ( ((yt .gt. yupper) .and. (yperdom .or. spheredom)) .or. &
       (((xl .lt. xlower) .or. (xr .gt. xupper)) .and. xperdom) ) then
     call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
          1,mjtot-ng+1,ilo-ng,ihi+ng,jhi+1,jhi+ng,ilo-ng,ihi+ng,jlo-ng,jhi+ng,patchOnly)
  else
     call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
          1,mjtot-ng+1,ilo-ng,ihi+ng,jhi+1,jhi+ng,patchOnly)
  endif


  ! set all exterior (physical)  boundary conditions for this grid at once
  ! used to be done from filpatch, but now only for recursive calls with new patch
  ! where the info matches. more efficient to do whole grid at once, and avoid copying
       write(34,*) '+++ in bound, spheredom = ',spheredom

  call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,hx,hy,level,time,xloWithGhost,xhiWithGHost, &
       yloWithGhost,yhiWithGhost,xlower,ylower,xupper,yupper,xperdom,yperdom,spheredom)

end subroutine bound
