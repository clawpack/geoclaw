!
! ----------------------------------------------------------
!
      subroutine setMatrixIndex(levelBouss)

      use amr_module
      use bouss_module
      implicit none

      integer, intent(in) :: levelBouss

      integer k,i,j,nb,nng,levSt, mptr,nx,ny,nx_nbor,ny_nbor
      integer ist,iend,jst,jend

      type(matrix_patchIndex), pointer :: mi
      type(matrix_patchIndex), pointer :: mi_nbor
      type(matrix_levInfo),  pointer :: minfo
      integer bndNum, imax,imin,jmax,jmin,icount,nextSpot,max_matrix_nelt
      integer ixlo,ixhi,jxlo,jxhi,imlo,imhi,jmlo,jmhi
      integer mtest,nbtest
      logical xlo_db,xhi_db,ylo_db,yhi_db,periodic
      logical xnbor_lo,xnbor_hi,ynbor_lo,ynbor_hi
      logical PATCH_TOUCHES_DOMAIN_BNDRY, debug

      PATCH_TOUCHES_DOMAIN_BNDRY = (xlo_db .or. xhi_db .or. ylo_db .or. yhi_db)

!
!     traverse grids in this exact order, after arrangeGrids has been called
!     to compute integer storage index for bouss matrix. 
!     also after makeBndryList called so can use it to find
!     neighboring grids at same level

!     Steps:
!      1. index array initialized to -1.
!      2. enumerate each grid, setting interior indices
!      3. go back and fill in overlapping ghost cell indices
!         cells at edges of grid  that still have a -1 in their stencil 
!         will be treated as SWE points.
!
      minfo => matrix_info_allLevs(levelBouss)
      periodic = xperdom .or. yperdom
      !debug = .true.
      debug = .false.

      k = 0  ! keep running total of number of unknowns in matrix system 
!     ## allocate array to store indices, only for levelBouss grids
      nb = 0  ! recount number of bouss grids to set pointers
      levSt = listStart(levelBouss)
      do nng = 1, numgrids(levelBouss)
         mptr = listOfGrids(levSt+nng-1)
         nb = nb + 1
         node(bpatchNum,mptr) = nb
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mi => minfo%matrix_indices(nb)
         allocate(mi%mindex(0:nx+1,0:ny+1))
         !mi%mindex = -1  ! indicates no matrix element set
         ! initialize on previous line doesnt work for allocate, need loop
         ! would be nicer to only go order cells but reset next loop
         do j = 0, ny+1
         do i = 0, nx+1
            mi%mindex(i,j) = -1  
         end do
         end do

         ! set interior with matrix index
         do j = 1, ny
         do i = 1, nx
            k = k + 1
            mi%mindex(i,j) = k
         end do
         end do
      end do ! end initializing the grids

!     save total number
      minfo%numBoussCells = k
      ! allocate overestimate since some cells will not have 5 entries at bndry
      ! remember there are twice num unknowns

      ! note do not generally need hv to left/right or hu above/below unless 
      ! introduced in place of a ghost cell, so max 14 elements per row, not 18
      if (ibouss .eq. 1) then
        max_matrix_nelt = 14 * minfo%numBoussCells 
      else ! ibouss  2 or 3
        max_matrix_nelt = 24 * minfo%numBoussCells 
      endif

      if (minfo%numBoussCells > 0) then
          if (.not. crs) then ! COO triplet format
            allocate(minfo%matrix_ia(max_matrix_nelt),     &
                     minfo%matrix_ja(max_matrix_nelt),     &
                     minfo%matrix_sa(max_matrix_nelt))
          else  ! CRS format
            allocate(minfo%cols(0:max_matrix_nelt),        &
                     minfo%vals(0:max_matrix_nelt),        &
                     minfo%rowPtr(0:2*minfo%numBoussCells))
          endif
      endif

! loop again to set boundary values from other grids
! this first part  is for non_periodic case, so it uses the bndList - which
! doesnt account for periodicity
      levSt = listStart(levelBouss)
      do nng = 1, numgrids(levelBouss)
         mptr = listOfGrids(levSt+nng-1)
!        set outer ghost cells by intersecting with other grids at this level
!        these grids only have one halo cell
         imin = node(ndilo,mptr)
         imax = node(ndihi,mptr)
         jmin = node(ndjlo,mptr)
         jmax = node(ndjhi,mptr)
         nx     = imax - imin + 1
         ny     = jmax - jmin + 1
         xlo_db = (node(ndilo,mptr) .eq. 0)
         xhi_db = (node(ndihi,mptr) .eq. (iregend(levelBouss)))
         ylo_db = (node(ndjlo,mptr) .eq. 0)
         yhi_db = (node(ndjhi,mptr) .eq. (jregend(levelBouss)))
         nb = node(bpatchNum,mptr)
         mi => minfo%matrix_indices(nb)

!        loop over this grids bndry list
         nextSpot = node(bndListSt,mptr)
         bndNum = node(bndListNum,mptr)
         do icount = 1, bndNum
           mtest = bndList(nextSpot,gridNbor)
!          check if grid mptr and mtest intersect
           imlo = node(ndilo,mtest) ! only test against INTERIOR of mtest
           imhi = node(ndihi,mtest)
           jmlo = node(ndjlo,mtest)
           jmhi = node(ndjhi,mtest)

           ixlo = max(imlo,imin-1)  ! our grid wants to fill its halo 
           ixhi = min(imhi,imax+1)  ! using other grids interior cells only
           jxlo = max(jmlo,jmin-1)
           jxhi = min(jmhi,jmax+1)

!          check if grids intersect, periodicity special case comes next
           if (ixlo > ixhi .or. jxlo > jxhi) cycle 

!          can only be one row or col, but dont know which so double loop
           nbtest = node(bpatchNum,mtest)
           mi_nbor => minfo%matrix_indices(nbtest)
           do j = jxlo,jxhi
           do i = ixlo,ixhi
              mi%mindex(i-imin+1,j-jmin+1) = mi_nbor%mindex(i-imlo+1,j-jmlo+1) 
           end do
           end do
           nextSpot = bndList(nextSpot,nextfree)
         end do  ! end loop over mptrs bndry grids 
      
         if ((PATCH_TOUCHES_DOMAIN_BNDRY) .and. periodic) then
            ! if grid mptr is periodic and touches at least one bndry
            ! do slower loop over all grids to fill periodic bndries
            levSt = listStart(levelBouss)
            do icount = 1, numgrids(levelBouss)
               mtest = listOfGrids(levSt+icount-1)
!              check if grid mptr and mtest intersect
               imlo = node(ndilo,mtest) ! only test against INTERIOR of mtest
               imhi = node(ndihi,mtest)
               jmlo = node(ndjlo,mtest)
               jmhi = node(ndjhi,mtest)
               nx_nbor  = imhi - imlo + 1
               ny_nbor  = jmhi - jmlo + 1
               nbtest = node(bpatchNum,mtest)
               mi_nbor => minfo%matrix_indices(nbtest)
               xnbor_lo = (node(ndilo,mtest) .eq. 0)
               xnbor_hi = (node(ndihi,mtest) .eq. (iregend(levelBouss)))
               ynbor_lo = (node(ndjlo,mtest) .eq. 0)
               ynbor_hi = (node(ndjhi,mtest) .eq. (jregend(levelBouss)))
  
               if (xperdom) then 
                  !jst  = max(jmin,jmlo)
                  !jend = min(jmax,jmhi)
                  ! adjust by 1 to get corner cells if possible
                  ! this will only work if not touching bottom bndry since
                  ! nbor grid only uses interior cell numbering. 
                  ! corners handled below
                  jst  = max(jmin-1,jmlo) 
                  jend = min(jmax+1,jmhi)
                  ! this doesnt include corner cells if stick out
                  ! would have too many cases. done after
                  if (xlo_db .and. xnbor_hi .and. jst .le. jend) then
                     ! see if any interesection and then copy
                     do j = jst, jend 
                        mi%mindex(0,j-jmin+1) = mi_nbor%mindex(nx_nbor,j-jmlo+1)
                     end do
                  endif
                  if (xhi_db .and. xnbor_lo) then
                    do j = jst, jend
                       mi%mindex(nx+1,j-jmin+1) = mi_nbor%mindex(1,j-jmlo+1)
                    end do
                  endif
               endif !periodic in x
               if (yperdom) then
                  ist  = max(imin-1,imlo)
                  iend = min(imax+1,imhi)
                  if (ylo_db .and. ynbor_hi .and. ist .le. iend) then
                    do i = ist, iend 
                       mi%mindex(i-imin+1,0) = mi_nbor%mindex(i-imlo+1,ny_nbor)
                    end do
                  endif
                  if (yhi_db .and. ynbor_lo) then
                    do i = ist, iend
                       mi%mindex(i-imin+1,ny+1) = mi_nbor%mindex(i-imlo+1,1)
                    end do
                  endif
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
                      mi%mindex(0,0) = mi_nbor%mindex(nx_nbor,ny_nbor) 
                  endif
                  if (xlo_db .and. yhi_db .and. xnbor_hi .and. ynbor_lo) then
                      mi%mindex(0,ny+1) = mi_nbor%mindex(nx_nbor,1) 
                  endif
                  if (xhi_db .and. ylo_db  .and. xnbor_lo .and. ynbor_hi) then
                      mi%mindex(nx+1,0) = mi_nbor%mindex(1,ny_nbor) 
                  endif
                  if (xhi_db .and. yhi_db  .and. xnbor_lo .and. ynbor_lo) then
                      mi%mindex(nx+1,ny+1) = mi_nbor%mindex(1,1)
                  endif
               endif
               
            end do  ! end loop over ALL grids at this level
         endif  ! end periodic case

         if (debug) then
            write(18,*)
            write(18,*)"grid ",mptr
            do j = ny+1,0,-1
               write(18,101)(mi%mindex(i,j),i=0,nx+1)
 101           format(300i7)
            end do
         endif
      end do  ! end loop over bouss grids

      return
      end subroutine setMatrixIndex
