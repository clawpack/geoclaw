c
c -----------------------------------------------------------
c
      subroutine regrid  (nvar,lbase,cut,naux,start_time)
c
      use amr_module
      use bouss_module
      implicit double precision (a-h,o-z)
      integer newnumgrids(maxlv)
      integer(kind=8) :: clock_start2, clock_finish, clock_rate
      type(matrix_patchIndex), pointer :: mi
      type(matrix_levInfo),  pointer :: minfo

c
c :::::::::::::::::::::::::::: REGRID :::::::::::::::::::::::::::::::

!> Flag points on each grid with a level > = lbase.
!! cluster them, and fit new subgrids around the clusters.
!! the lbase grids stay fixed during regridding operation.
!! when a parent grid has its error estimated, add its kid grid
!! information to the error grid before clustering. (project)
!! order of grid examination - all grids at the same level, then
!! do the next coarser level.
c
c input parameters:
c     lbase  = highest level that stays fixed during regridding
c     cutoff = criteria for measuring goodness of rect. fit.

c local variables:
c     lcheck = the level being examined.
c     lfnew  = finest grid to be. will replace lfine.

c global
c    mstart  = start of very coarsest grids.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      lcheck    = min0(lfine,mxnest-1)
      lfnew     = lbase
      do  i   = 1, mxnest
        newnumgrids(i) = 0
        newstl(i) = 0
      end do
      time      = rnode(timemult, lstart(lbase))
c
 20   if (lcheck .lt. lbase) go to 50
          call grdfit(lbase,lcheck,nvar,naux,cut,time,start_time)
          if (newstl(lcheck+1) .eq. 0) go to 40
          lfnew = max0(lcheck + 1,lfnew)
 40       continue
          lcheck = lcheck - 1
c
      go to 20
 50   continue
c
c  end of level loop
c
c  remaining tasks left in regridding:
c  1.  count number of new grids at each level  
      maxnumnewgrids = 0   ! max over all levels. needed for dimensioning
      do lev = lbase+1,lfnew 
          ngridcount = 0
          mptr = newstl(lev)
 52       if (mptr .eq. 0) go to 55
             ngridcount = ngridcount + 1
             mptr = node(levelptr,mptr)
             go to 52

 55       newnumgrids(lev) = ngridcount
          maxnumnewgrids = max(maxnumnewgrids,ngridcount)
      end do
c
c  2. interpolate storage for the new grids.  the starting pointers
c  for each level are in newstl. also reclaim some space before new
c  allocations.
      call system_clock(clock_start2,clock_rate)
      call gfixup(lbase,lfnew,nvar,naux,newnumgrids,maxnumnewgrids)
      call system_clock(clock_finish,clock_rate)
      timeGrdfit2 = timeGrdfit2 + clock_finish - clock_start2
c
c  3. merge data structures (newstl and lstart )
c  finish storage allocation, reclaim space, etc. set up boundary
c  flux conservation arrays
c

c     get rid of old bouss storage
      if (max(lbase+1,minLevelBouss) .le. maxLevelBouss .and.
     &     ibouss .gt. 0) then  ! need to redo bouss stuff
        do lev = max(lbase+1,minLevelBouss), maxLevelBouss
          minfo => matrix_info_allLevs(lev)
          if (minfo%numBoussGrids .gt. 0) then
            do j = 1, minfo%numBoussGrids
               mi => minfo%matrix_indices(j)
               deallocate(mi%mindex)
               deallocate(mi%isBouss)
            end do
            deallocate(minfo%matrix_indices) ! get rid of old one 
            ! also deallocate matrix since new one will have diff size
            if (.not. crs) then ! COO triplet format
             deallocate(minfo%matrix_ia,minfo%matrix_ja,minfo%matrix_sa)
            else  ! CRS format
             deallocate(minfo%rowPtr, minfo%cols, minfo%vals)
            endif
            minfo%numBoussGrids = 0  ! reset for new counting
            minfo%numBoussCells = 0  ! reset for new counting
          endif
        end do
      endif

!     also get rid of pardiso storage, new levels need to new factorizations
      if (isolver .eq. 2) then
#ifdef HAVE_PARDISO
        phase = -1  ! release internal memory
        msglvl = 0
        do lev = lbase+1,mxnest
           newGrids(lev) = .true.
           if (lev.le.maxLevelBouss .and.  lev.ge.minLevelBouss) then
               !write(*,*)"regrid is releasing pardiso for lev ",lev
               CALL pardiso(pt(1,lev),maxfct,mnum,mtype,phase,nn,ddum,
     &                    idum, idum,idum,nrhs,iparm,msglvl,
     &                    ddum,ddum,error,dparm)
           endif
         end do
#endif
      else if (isolver .eq. 3) then
         newGrids(lbase+1:maxlv) = .true.
      endif

c     set deepest bathy for each grid
      if (lbase .lt. maxLevelBouss .and. ibouss.gt.0) then ! set up new bouss stuff
        do lev = lbase+1, lfine
           minfo => matrix_info_allLevs(lev)
           mptr = lstart(lev)
   74         continue
              nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              locaux = node(storeaux,mptr)
              !call setMaxMinDepth(alloc(locaux),mptr,naux,nx,ny,nghost)
              if (lev.le.maxLevelBouss .and. lev.ge.minLevelBouss) then
                    minfo%numBoussGrids = minfo%numBoussGrids+1
              endif
              mptr = node(levelptr, mptr)
              if (mptr .ne. 0) go to 74
           ! done with grids at this level
           if (minfo%numBoussGrids .gt. 0) then
              allocate(minfo%matrix_indices(minfo%numBoussGrids))
           endif
        end do
      endif

      ! initialize pardiso for new grids at this level
      if (isolver .eq. 2 .and. ibouss.gt.0) then
#ifdef HAVE_PARDISO
        do lev = lbase+1, mxnest
           if (lev.le.maxLevelBouss .and. lev.ge.minLevelBouss) then
                call pardisoinit(pt(1,lev),mtype,solver,iparm,dparm,
     &                           error)
                !write(*,*)"regrid initialized pardiso for lev ",lev
           endif
        end do
#endif
      endif
c
c  reset numgrids per level, needed for omp parallelization.
c  note that grids may have disappeared, so next loop resets to 0
c  if there are no grids from lfine+1 to mxnest
c
      do 72 levnew = lbase+1, mxnest
        mptr = lstart(levnew)
        ngridcount = 0
        ncells = 0
        do while (mptr .gt. 0)
           ngridcount = ngridcount + 1
           ncells = ncells + (node(ndihi,mptr)-node(ndilo,mptr)+1)
     .                     * (node(ndjhi,mptr)-node(ndjlo,mptr)+1)
           mptr = node(levelptr, mptr)
         end do
         if (ngridcount .ne. newnumgrids(levnew)) then
            write(*,*)"regrid grid count off at time",time
            stop
         endif
         numgrids(levnew) = ngridcount
         numcells(levnew) = ncells
         avenumgrids(levnew) = avenumgrids(levnew) + ngridcount
         iregridcount(levnew) = iregridcount(levnew) + 1
c        sort grids to first ones are the most work. this helps load
c        balancing, but doesn't help locality
         !commented out enxt line because now already sorted in prepnewgrids above
         ! changed 9/25/23 by mjb for better grid placement in filval par for loop 
         if (ngridcount .gt. 1) call arrangeGrids(levnew,ngridcount)

         if (verbosity_regrid .ge. levnew) then
           write(*,100) ngridcount,ncells,levnew
           write(outunit,100) ngridcount,ncells,levnew
 100       format("there are ",i6," grids with ",i10,
     &            " cells at level ", i3)
         endif
72     continue

      do 60 level = lbase, lfine-1
        call prepf(level+1,nvar,naux)
        call prepc(level,nvar)
 60   continue
c
c      set up array of grids instead of recomputing at each step
       call makeGridList(lbase)
       do levnew = lbase+1, lfine
          call makeBndryList(levnew)   ! does one level at a time
       end do
c
c     set up for Bouss matrices
      if (ibouss .gt. 0) then ! 0 is SWE, otherwise setup for bouss
         do lev = lbase+1,lfine
           if (lev.le. maxLevelBouss .and. lev .ge. minLevelBouss)then
              call setMatrixIndex(lev)
              call setBoussFlag(lev,naux)
           endif
         end do
      endif
c
      return
      end
c
c -------------------------------------------------------------------
c
!> Sort all grids at level **level**. 
!! Put the most expensive grid in **lstart(level)** 
!! Cost is measured by number of cells.
!! The linked list was also sorted such that the cost of grids
!! decreases from list head to list tail

      subroutine arrangeGrids(level, numg)
c
      use amr_module
      implicit double precision (a-h,o-z)
      integer listgrids(numg), cost(numg), index(numg), prevptr
c
c   slow sort for now, putting most expensive grids first on lstart list
c   measure cost by number of cells
c
       mptr = lstart(level)
       do i = 1, numg
         listgrids(i) = mptr
         cost(i) =  (node(ndihi,mptr)-node(ndilo,mptr)+3) *
     1              (node(ndjhi,mptr)-node(ndjlo,mptr)+2)
         index(i) = i
         mptr = node(levelptr, mptr)
       end do
c
c        write(*,*)" before sorting"
c       write(*,*) index
c
       call  qsorti(index, numg, cost)

c       write(*,*)"after sorting"
c       write(*,*) index

c qsort returns in ascending order, repack in descending order
c grids can stay in place, just their levelptrs need to change
       lstart(level) = listgrids(index(numg))  ! last grid is most expensive
       prevptr = listgrids(index(numg))
       do i = 1, numg-1             
          node(levelptr, prevptr) = listgrids(index(numg-i))
          prevptr = listgrids(index(numg-i))
       end do
       node(levelptr,prevptr) = null  !signal the last grid

       return
       end
