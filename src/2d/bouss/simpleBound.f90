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
       integer :: ilo,ihi,jlo,jhi,iloadj,ihiadj,jloadj,jhiadj,i,j,mtest
       integer :: imlo,jmlo,imhi,jmhi,mgridmitot,mgridmjtot
       integer :: locsrc,ialloc,level,levSt
       integer :: ixlo,ixhi,jxlo,jxhi,bndLoc,bndNum,nextSpot,numg,icount,mgrid
       real(kind=8) :: hx,hy,xloWithGhost,xhiwithGhost,yloWithGhost,yhiWithGhost
       integer :: ist, iend, ivar, jst, jend,nx_nbor,ny_nbor
       logical :: xlo_db,xhi_db,ylo_db,yhi_db
       logical :: xnbor_lo,xnbor_hi,ynbor_lo,ynbor_hi
       logical :: PATCH_TOUCHES_DOMAIN_BNDRY, periodic
       logical :: newlogic,newlogic2
       logical :: debug

       PATCH_TOUCHES_DOMAIN_BNDRY = (xlo_db .or. xhi_db .or. ylo_db .or. yhi_db)

!      replace fine grid ghost cells with interior vals from other
!      grids at same level.  
!      ghost cells will already have been filled by coarser grids, 
!      so overwrite now if something better is available.
!      Feb 2025 adding periodicity, but only checking this level, no recursive calls

       ! this routine uses fact that variables consecutive in alloc
       ! also assumes nvar 5 since wouldnt be in implicit update for SWE
       if (nvar .ne. 5) then
          write(*,*) "rewrite this routine for  nvar different than 5"
          stop
       endif
       periodic = xperdom .or. yperdom
       debug = .false.
       level = node(nestlevel,mptr)
       xlo_db = (node(ndilo,mptr) .eq. 0)
       xhi_db = (node(ndihi,mptr) .eq. (iregend(level)))
       ylo_db = (node(ndjlo,mptr) .eq. 0)
       yhi_db = (node(ndjhi,mptr) .eq. (jregend(level)))

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
          locsrc = node(store1, mgrid)

          ixlo = max(imlo,iloadj)
          ixhi = min(imhi,ihiadj)
          jxlo = max(jmlo,jloadj)
          jxhi = min(jmhi,jhiadj)

          if (ixlo > ixhi  .or. jxlo > jxhi) go to 10

          ! only ghost cells should intersect
          if (debug) write(*,*)" Interior ghost cells: mptr ",mptr," mgrid ",mgrid
          do j = jxlo, jxhi
          do i = ixlo, ixhi
             if (debug) then
               write(*,*)"indices of mptr ",i-ilo+nghost+1, j-jlo+nghost+1
               write(*,*)"indices of mgrid ",i-imlo+nghost+1, j-jmlo+nghost+1
             endif
             ialloc = iadd(2, i-imlo+nghost+1, j-jmlo+nghost+1)
             val(2:5,i-ilo+nghost+1,j-jlo+nghost+1) = alloc(ialloc:ialloc+3)
          end do
          end do

 10       nextSpot = bndList(nextSpot,nextfree) ! get next grid

       end do

       newlogic = PATCH_TOUCHES_DOMAIN_BNDRY
       !write(*,*) " newlogic 1  ",newlogic
       !write(*,*) " periodic ",periodic
       newlogic2 = (periodic .and. newlogic)
       !write(*,*) " newlogic 2 ",newlogic2
       !if (newlogic2) then
       if (.true.) then
       !if ((PATCH_TOUCHES_DOMAIN_BNDRY) .and. periodic) then
            ! if grid mptr is periodic and touches at least one bndry
            ! do slower loop over all grids to fill periodic bndries
            levSt = listStart(level)
            do icount = 1, numgrids(level)
               mtest = listOfGrids(levSt+icount-1)
!              check if grid mptr and mtest intersect
               imlo = node(ndilo,mtest)
               imhi = node(ndihi,mtest)
               jmlo = node(ndjlo,mtest)
               jmhi = node(ndjhi,mtest)
               mgridmitot = imhi - imlo + 1 + 2*nghost
               mgridmjtot = jmhi - jmlo + 1 + 2*nghost
               nx_nbor = imhi - imlo + 1
               ny_nbor = jmhi - jmlo + 1
               locsrc = node(store1, mtest)
               xnbor_lo = (node(ndilo,mtest) .eq. 0)
               xnbor_hi = (node(ndihi,mtest) .eq. (iregend(level)))
               ynbor_lo = (node(ndjlo,mtest) .eq. 0)
               ynbor_hi = (node(ndjhi,mtest) .eq. (jregend(level)))
               if (xperdom) then
                  ! adjust by num ghost cells to get corner cells if possible
                  ! this will only work if not touching bottom bndry since
                  ! nbor grid only uses interior cell numbering. 
                  ! corners handled below
                  jst  = max(jlo-nghost,jmlo)
                  jend = min(jhi+nghost,jmhi)
                  ! this doesnt include corner cells if stick out
                  ! would have too many cases. done after
                  ! see if any intersection and then copy
                  if (jst .le. jend) then 
                   if (xlo_db .and. xnbor_hi) then
                     if (debug) write(*,*)" Left bndry of mptr ",mptr," mtest ",mtest
                     do j = jst, jend
                     do i  = 1, nghost
                        ialloc = iadd(2, imhi-imlo+nghost-1+i, j-jmlo+nghost+1)
                        if (debug) then
                          !write(*,*)" loop indices ",i,j
                          write(*,*)"indices of mptr ",i, j-jlo+nghost+1
                          write(*,*)"indices of mtest ",imhi-imlo+nghost+i-1, j-jmlo+nghost+1
                        endif
                        val(2:5,i,j-jlo+nghost+1) = alloc(ialloc:ialloc+3)
                     end do
                     end do
                   endif ! mptr touches right and mtest touches left
                   if (xhi_db .and. xnbor_lo) then
                     if (debug)write(*,*)" Right bndry of mptr ", mptr," mtest ",mtest
                    do j = jst, jend
                    do i = ihi+1+nghost+1,ihi+1+2*nghost   ! mitot-nghost+1,mitot
                        if (debug) then
                          !write(*,*)" loop indices ",i,j
                          write(*,*)"indices of mptr ",i-ilo, j-jlo+nghost+1
                          write(*,*)"indices of mtest ",i-ihi-1, j-jmlo+nghost+1
                        endif
                       ialloc = iadd(2, i-ihi-1, j-jmlo+nghost+1)
                       val(2:5,i-ilo,j-jlo+nghost+1) = alloc(ialloc:ialloc+3)
                    end do
                    end do
                  endif ! mptr touches left and mtest touches right
                endif ! columns in range (jst <= jend)
               endif !periodic in x

               if (yperdom) then
                  ist  = max(ilo-nghost,imlo)
                  iend = min(ihi+nghost,imhi)
                  if (ist .le. iend) then
                   if (debug)  write(*,*)" Bot bndry of mptr ",mptr," mtest ",mtest
                   if (ylo_db .and. ynbor_hi) then
                    do j = 1,nghost   
                    do i = ist, iend
                       if (debug) then
                         !write(*,*)" loop indices ",i,j
                         write(*,*)"indices of mptr ",i-ilo+nghost+1,j
                         write(*,*)"indices of mtest ",i-imlo+nghost+1, jmhi-jmlo+nghost+j-1
                       endif
                       ialloc = iadd(2, i-imlo+nghost+1, jmhi-jmlo+nghost+j-1)
                       val(2:5,i-ilo+nghost+1,j) = alloc(ialloc:ialloc+3)
                    end do
                    end do
                   endif ! mptr touches top and mtest touches bot
                    if (debug) write(*,*)" Top bndry of mptr ",mptr," mtest ",mtest
                   if (yhi_db .and. ynbor_lo) then
                    do j = 1, nghost 
                    do i = ist, iend
                       if (debug) then
                         !write(*,*)" loop indices ",i,j
                         write(*,*)"indices of mptr ",i-ilo+nghost+1,mjtot-nghost+j
                         write(*,*)"indices of mtest ",i-imlo+nghost+1, nghost+j
                       endif
                       ialloc = iadd(2, i-imlo+nghost+1, nghost+j)
                       val(2:5,i-ilo+nghost+1,mjtot-nghost+j) = alloc(ialloc:ialloc+3)
                    end do
                    end do
                   endif ! mptr touches top and mtest touches bot
                  endif ! rows in range (ist <= iend)
               endif !periodic in y

               ! handle 4 corner cells if doubly periodic, otherwise leave as -1
               ! remember that numbering in mindex is 1 based, so 0 and nx/y+1
               ! are ghost cells
               if (xperdom .and. yperdom) then
                  ! need vals for corner cells if on periodic bndry
                  ! e.g. for single grid these are
                  ! (0,0),(0,reg_y+1),(reg_x+1,0),(reg_x+1,reg_y+1)
                  ! should come from (reg_x,reg_y),(reg_x,1),(1,reg_y),(1,1)
                  ! where reg_x,y are the iregend and jregend integer indices
                  ! in other words the source grid touches the appropriate bndry
                  ! grid indexing is 0 based, mindex is 1 based.
                  if (xlo_db .and. ylo_db .and. xnbor_hi .and. ynbor_hi) then
                      if (debug) write(*,*)"lower left corner mptr ",mptr," mtest ",mtest
                      do j = 1, nghost
                      do i = 1, nghost
                        ialloc = iadd(2,nx_nbor+i,ny_nbor+j)
                        if (debug) then
                          write(*,*)"indices of mptr ",i,j
                          write(*,*)"indices of mtest ",nx_nbor-nghost+i,ny_nbor-nghost+j
                        endif
                        val(2:5,i,j) = alloc(ialloc:ialloc+3)
                      end do
                      end do
                  endif
                  if (xlo_db .and. yhi_db .and. xnbor_hi .and. ynbor_lo) then
                      if (debug) write(*,*)"top left corner mptr ",mptr," mtest ",mtest
                      do j = 1, nghost
                      do i = 1, nghost
                        ialloc = iadd(2,nx_nbor+i,nghost+j)
                        if (debug) then
                          write(*,*)"indices of mptr ",i,mjtot-nghost+j
                          write(*,*)"indices of mtest ",nx_nbor+i,nghost+j
                        endif
                        val(2:5,i,mjtot-nghost+j) = alloc(ialloc:ialloc+3)
                      end do
                      end do
                  endif
                  if (xhi_db .and. ylo_db  .and. xnbor_lo .and. ynbor_hi) then
                      if (debug) write(*,*)"bottom right corner mptr ",mptr," mtest ",mtest
                      do j = 1, nghost
                      do i = 1, nghost
                        ialloc = iadd(2,nghost+i,ny_nbor+j)
                        if (debug) then
                          write(*,*)"indices of mptr ",mitot-nghost+i,j
                          write(*,*)"indices of mtest ",nghost+i,ny_nbor+j
                        endif
                        val(2:5,mitot-nghost+i,j) = alloc(ialloc:ialloc+3)
                      end do
                      end do
                  endif
                  if (xhi_db .and. yhi_db  .and. xnbor_lo .and. ynbor_lo) then
                      if (debug) write(*,*)"top right corner"
                      do j = 1, nghost
                      do i = 1, nghost
                        ialloc = iadd(2,nghost+i,nghost+j)
                        if (debug) then
                          write(*,*)"indices of mptr ",mitot-nghost+i,mjtot-nghost+j
                          write(*,*)"indices of mtest ",nghost+i,nghost+j
                        endif
                        val(2:5,mitot-nghost+i,mjtot-nghost+j) = alloc(ialloc:ialloc+3)
                      end do
                      end do
                  endif
               endif

            end do ! end loop over all grids at this level
       endif ! end periodic case

!  finish by calling bc2amr in case there are domain bndry touching patches
      call bc2amr(val,aux,mitot,mjtot,nvar,naux,hx,hy,level,time,                &
                  xloWithGhost,xhiWithGhost,yloWithGhost,yhiWithGhost)
contains

    integer pure function iadd(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar,i,j
        iadd = locsrc + ivar-1 + nvar*((j-1)*mgridmitot+i-1)
     end function iadd

end subroutine simpleBound
