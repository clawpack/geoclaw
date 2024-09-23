!
! ----------------------------------------------------------
!
      subroutine setMatrixIndex(levelBouss)

      use amr_module
      use bouss_module
      implicit none

      integer, intent(in) :: levelBouss

      integer k,i,j,nb,nng,levSt, mptr,nx,ny 

      type(matrix_patchIndex), pointer :: mi
      type(matrix_patchIndex), pointer :: mi_nbor
      type(matrix_levInfo),  pointer :: minfo
      integer bndNum, imax,imin,jmax,jmin,icount,nextSpot,max_matrix_nelt
      integer ixlo,ixhi,jxlo,jxhi,imlo,imhi,jmlo,jmhi,mtest,nbtest

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
      end do

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

!     loop again to set boundary values from other grids
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

           if (ixlo > ixhi .or. jxlo > jxhi) cycle !dont intersect
!          can only be one row or col, but dont know which so double loop
           nbtest = node(bpatchNum,mtest)
           mi_nbor => minfo%matrix_indices(nbtest)
           do j = jxlo,jxhi
           do i = ixlo,ixhi
              mi%mindex(i-imin+1,j-jmin+1) = mi_nbor%mindex(i-imlo+1,j-jmlo+1) 
           end do
           end do
           nextSpot = bndList(nextSpot,nextfree)
         end do  ! end loop over bndry grids for mptr
         if (0 .eq. 1) then
            write(*,*)
            write(*,*)"grid ",mptr
            do j = ny+1,0,-1
               write(*,101)(mi%mindex(i,j),i=0,nx+1)
 101           format(300i7)
            end do
         endif
      end do  ! end loop over bouss grids

      return
      end subroutine setMatrixIndex
