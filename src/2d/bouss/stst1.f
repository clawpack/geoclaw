c
c --------------------------------------------------------------
c
      subroutine stst1
c
      use amr_module
      implicit double precision (a-h,o-z)


c
c :::::::::::::::::::::::::::::: STST1 ::::::::::::::::::::::::::::::::
c    intialize a few variables needed before calling user set up
c    routine domain.
c    the spatial and temporal stepsizes are set. the node array
c    is kept as a linked list of free nodes.  "ndfree" points to the
c    head of the list, i.e. first free node.  use first row of each
c    col to hold this pointer, set by the macro "nextfree".
c    the free space list, managed in lfree, will have first and
c    last positions filled with an allocation of zero words,
c    to avoid boundary cases.
c ::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::::
c
      ndfree = 1

      !! now rnode and node are allocatable to allow resizing
      call init_nodes()

      ! node space allocated, now thread it
      do 10 i   = 1, maxgr
         node(nextfree,i) = i+1
 10   continue
c
c the last free node will have a null pointer
 
      node(nextfree, maxgr) = null
c

c     Initialize dynamic memory
      call init_alloc()

      lfine = 1
C       do 20 i  = 1, memsize
C         alloc(i) = WEIRD
C  20   continue
c
c  initialize linked list of alloc storage as well.
c  first and last locations are dummy placeholders of zero words
c  of allocation each, to avoid boundary cases.
c
      do  40 i  = 1, lfdim
        lfree(i,1) = 0
        lfree(i,2) = 0
 40   continue
      lfree(3,1) =memsize + 2
      lfree(2,1) = 1
      lfree(2,2) =memsize
      lenf       = 3


c  need to manage the boundary List too
c     do i = 1, bndListSize
c        bndList(i,nextfree) = i+1
c     end do
c     bndList(bndListSize,nextfree) = null
c     ndfree_bnd = 1
      call initBndryList()
c
c after kcheck integrations of parent grid, move its refinements.
c finest level grid never needs to have its finer subgrids moved.
c
      call initTimers()   ! used to be done here, but needs to be called from restarting too when stst1 not called

      do 60 i   = 1, maxlv
         tvoll(i) = 0.d0
         iregridcount(i) = 0
         avenumgrids(i) = 0
         numgrids(i) = 0
         numcells(i) = 0
         lstart(i) = 0
 60      icheck(i) = 0
c
c finish initializing spatial and counting arrays
c
      level      = 2
 70   if (level .gt. mxnest) go to 80
          hxposs(level) = hxposs(level-1) / dble(intratx(level-1))
          hyposs(level) = hyposs(level-1) / dble(intraty(level-1))
          level         = level + 1
      go to 70
 80   continue


      return
      end
c
c -------------------------------------------------------------------------
c
      subroutine initTimers()

      use amr_module
      !implicit double precision (a-h,o-z)

      timeBufnst         = 0
      timeStepgrid       = 0
      timeStepgridCPU    = 0.d0
      timeBound          = 0
      timeBoundCPU       = 0.d0
      timeGrdfitAll      = 0
      timeFlagger        = 0
      timeFlglvl         = 0
      timeFlglvlTot      = 0
      timeGrdfit2        = 0
      timeLinSolve       = 0
      timeLinSolveCPU    = 0.d0
      timePrepBuild      = 0
      timePrepBuildCPU    = 0.d0
      timeRegridding     = 0
      timeRegriddingCPU  = 0.d0
      timeTick           = 0
      timeTickCPU        = 0.d0
      timeUpdating       = 0
      timeValout         = 0
      timeValoutCPU      = 0.d0
      avgIterRef         = 0.d0
      maxIterRef         = 0.d0

      return
      end
