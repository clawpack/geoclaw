!
! ----------------------------------------------------------------
!
subroutine simpleBound(time,nvar,val,mitot,mjtot,mptr,aux,naux)

       use amr_module
       implicit none

       !input
       integer, intent(in) :: nvar, naux, mitot, mjtot, mptr
       real(kind=8), intent(in) :: time
       real(kind=8), intent(inout) :: val(nvar,mitot,mjtot)
       real(kind=8), intent(inout) :: aux(naux,mitot,mjtot)

       !locals
       integer :: ilo,ihi,jlo,jhi,iloadj,ihiadj,jloadj,jhiadj,i,j
       integer :: imlo,jmlo,imhi,jmhi,mgridmitot,locnew,ialloc,level
       integer :: ixlo,ixhi,jxlo,jxhi,bndLoc,bndNum,nextSpot,numg,icount,mgrid
       real(kind=8) :: hx,hy,xloWithGhost,xhiwithGhost,yloWithGhost,yhiWithGhost
       integer :: ist, iend, ivar, jst, jend

!      replace fine grid ghost cells with interior vals from other
!      grids at same level.  
!      ghost cells will already have been filled by coarser grids, 
!      so overwrite now if something better is available.

       level = node(nestlevel,mptr)
       hx = hxposs(level)
       hy = hyposs(level)
       xloWithGhost = rnode(cornxlo,mptr) - nghost*hx
       xhiWithGhost = rnode(cornxhi,mptr) + nghost*hx
       yloWithGhost = rnode(cornylo,mptr) - nghost*hy
       yhiWithGhost = rnode(cornyhi,mptr) + nghost*hy

       ! since filling ghost cells, adjust the integer indices to include them
       ilo = node(ndilo, mptr) 
       jlo = node(ndjlo, mptr) 
       ihi = node(ndihi, mptr) 
       jhi = node(ndjhi, mptr)

       iloadj = ilo - nghost
       jloadj = jlo - nghost
       ihiadj = ihi + nghost
       jhiadj = jhi + nghost

       !loop over bndry list to set ghost cells from grids at this level only
       ! then finish by calling bc2amr for domain bndry patches
       !
       ! NB this routine should overwrite and replace  the coarse interpolated 
       ! fine ghost cells from fillFine, called earlier,
       ! if there are neighboring fine grids 

       ! get list info of adjacent grids
       bndLoc = node(bndListSt, mptr)
       bndNum = node(bndListNum, mptr)
       nextSpot = node(bndListSt,mptr)
       numg = bndNum

       do icount = 1, numg
          mgrid = bndlist(nextspot, gridNbor)
          ! but use regular indices of neighboring grids as sources to fill
          imlo = node(ndilo, mgrid)
          jmlo = node(ndjlo, mgrid)
          imhi = node(ndihi, mgrid)
          jmhi = node(ndjhi, mgrid)
          mgridmitot = imhi - imlo + 1 + 2*nghost
          locnew = node(store1, mgrid)

          ixlo = max(imlo,iloadj)
          ixhi = min(imhi,ihiadj)
          jxlo = max(jmlo,jloadj)
          jxhi = min(jmhi,jhiadj)

          if (ixlo > ixhi  .or. jxlo > jxhi) go to 10

          ! only ghost cells should intersect, but do double loop since dont know which side
          do j = jxlo, jxhi
          do i = ixlo, ixhi
             ialloc = iadd(2, i-imlo+nghost+1, j-jmlo+nghost+1)
             val(2,i-ilo+nghost+1,j-jlo+nghost+1) = alloc(ialloc)
             val(3,i-ilo+nghost+1,j-jlo+nghost+1) = alloc(ialloc+1)
             if (nvar .eq. 5) then ! use fact that variables are consecutive in memory
                 val(4,i-ilo+nghost+1,j-jlo+nghost+1) = alloc(ialloc+2)
                 val(5,i-ilo+nghost+1,j-jlo+nghost+1) = alloc(ialloc+3)
             endif
          end do
          end do

 10       nextSpot = bndList(nextSpot,nextfree) ! get next grid

       end do

!  finish by calling bc2amr in case there are domain bndry touching patches
      call bc2amr(val,aux,mitot,mjtot,nvar,naux,hx,hy,level,time,                &
                  xloWithGhost,xhiWithGhost,yloWithGhost,yhiWithGhost)
contains

    integer pure function iadd(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar,i,j
        iadd = locnew + ivar-1 + nvar*((j-1)*mgridmitot+i-1)
     end function iadd

end subroutine simpleBound
