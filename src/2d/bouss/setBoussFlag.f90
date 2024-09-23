!
! ----------------------------------------------------------
!
      subroutine setBoussFlag(levelBouss,maux)

      use amr_module
      use bouss_module
      implicit none

      integer, intent(in) :: levelBouss

      integer k,m,nb,i,j,nng,levSt, mptr,nx,ny,mitot 

      type(matrix_patchIndex), pointer :: mi
      type(matrix_patchIndex), pointer :: mi_nbor
      type(matrix_levInfo),  pointer :: minfo
      integer :: locaux, iaddaux,maux
      integer :: trueCount
      real(kind=8) :: b_min,b_max,time
      logical :: debug
      integer :: ioff, joff
      real(kind=8) :: x,y,xlow,ylow,dx,dy

      iaddaux(m,i,j) = locaux +  maux*(i+nghost-1)   &
                     + maux*mitot*(j+nghost-1)

!
!     mark each cell in a bouss grid as to whether its deep enough
!     to use boussinesq and not SWE.

      minfo => matrix_info_allLevs(levelBouss)
      !debug = .true.
      debug = .false.

      k = 0  ! keep running total of number of unknowns in matrix system 
      levSt = listStart(levelBouss)
      do nng = 1, numgrids(levelBouss)
         mptr = listOfGrids(levSt+nng-1)
         time = rnode(timemult,mptr)
         nb = node(bpatchNum,mptr)
         mi => minfo%matrix_indices(nb)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         if (.not. allocated(mi%isBouss)) then
            allocate(mi%isBouss(0:nx+1,0:ny+1))
         endif
         mitot = nx + 2*nghost
         locaux = node(storeaux,mptr)
         trueCount = 0

         ! indexing to add ghost cells is in iaddaux
         do j = 0, ny+1
         do i = 0, nx+1
             b_max = max(alloc(iaddaux(1,i,j)),                                 &
                         alloc(iaddaux(1,i-1,j)),alloc(iaddaux(1,i+1,j)),       &
                         alloc(iaddaux(1,i,j-1)),alloc(iaddaux(1,i,j+1)),       &
                         alloc(iaddaux(1,i-1,j-1)),alloc(iaddaux(1,i-1,j+1)),   &
                         alloc(iaddaux(1,i+1,j-1)),alloc(iaddaux(1,i+1,j+1)))
             if (-boussMinDepth < b_max) then
               mi%isBouss(i,j) = .false.
             else
               mi%isBouss(i,j) = .true.
               trueCount = trueCount + 1
            endif
         end do
         end do

         if (debug) then
            write(55,*)
            write(55,108) mptr, levelBouss, trueCount, (nx+2)*(ny+2)
 108        format("grid ",i5," level ",i5," has ",i5," of ",i5," bouss cells")
            xlow = rnode(cornxlo,mptr)
            ylow = rnode(cornylo,mptr)
            dx = hxposs(levelBouss)
            dy = hyposs(levelBouss)
            do j = ny+1,0,-1
               y = ylow + (dfloat(j)-.5)*dy
               joff = j + node(ndjlo,mptr)-1
            do i = 0, nx+1 
               x = xlow + (dfloat(i)-.5)*dx
               ioff = i + node(ndilo,mptr)-1
               write(55,104)ioff,joff,x,y,mi%isBouss(i,j)
 101           format(300l1)
 104           format(2i6,2e12.4,l4)
            end do
            end do
            flush(55)
            if (1 .eq. 1) then
            write(56,*)
            write(56,*)"grid ",mptr
            !do j = ny+2,-1,-1
            do j = -1, ny+2
            do i = -1, nx+2
               !write(56,102)(alloc(iaddaux(1,i,j)),i=-1,nx+2)
               write(56,103)i+nghost,j+nghost,alloc(iaddaux(1,i,j))
 102           format(300e11.4)
 103           format(2i5,e12.4)
            end do
            end do
            endif
         endif
      end do  ! end loop over bouss grids

      return
      end subroutine setBoussFlag
