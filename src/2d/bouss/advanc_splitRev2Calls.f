c
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      use amr_module
      use bouss_module
      !use fixedgrids_module
      use topo_module, only: topo_finalized

      implicit double precision (a-h,o-z)


      logical vtime, doUpdate, implicitSolveDone
      logical debug
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(level))
      integer(kind=8) :: clock_start, clock_finish, clock_rate
      integer(kind=8) :: clock_startStepgrid, clock_finishStepgrid
      integer(kind=8) :: clock_startBound, clock_finishBound
      real(kind=8) cpu_start, cpu_finish
      real(kind=8) cpu_startBound,cpu_finishBound
      real(kind=8) cpu_startStepgrid, cpu_finishStepgrid

      type(matrix_levInfo),  pointer :: minfo


c     maxgr is maximum number of grids  many things are
c     dimensioned at, so this is overall. only 1d array
c     though so should suffice. problem is
c     not being able to dimension at maxthreads


c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells 
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
#ifdef WHERE_AM_I
      write(*,*)"starting advance with level ",level
#endif
      flush(66)

      !debug = .true.
      debug = .false.

      call system_clock(clock_start,clock_rate)
      call cpu_time(cpu_start)
      hx   = hxposs(level)
      hy   = hyposs(level)
      delt = possk(level)
      time = rnode(timemult,lstart(level))
      implicitSolveDone = .false.
      ! need to start on level 1 once enough time passed
      ! include .not. startWithBouss so only print once
      if (level .eq. 1 .and. startBoussTime .le. time 
     &                 .and. .not. startWithBouss) then
            write(*,*)"Using/Switching to Boussinesq equations"
            write(outunit,*)"Using/Switching to Boussinesq equations",
     &                      "at time ",time
            startWithBouss = .true.
      endif
c 
      if (debug) then
         write(*,*)"at start of  advanc level ",level
         call valout(level,level,time,nvar,naux)
      endif
c
      !! initial conditions were put in new storage
      !! save to old. implicit_update counts on this
      if (level .lt. lfine) call copyNew2Old(level,nvar)
c
c     set up for Bouss step
      if (ibouss .gt. 0) then
       if (level .ge. minLevelBouss .and. level .le. maxLevelBouss .and.
     &     startWithBouss) then
         minfo =>  matrix_info_allLevs(level) 
         doUpdate = .true. ! update momenta after implicit solve
         call implicit_update(nvar,naux,level,minfo%numBoussCells,
     &                        doUpdate,time)
         implicitSolveDone = .true. ! also indicates that bound was called
       endif
      endif

      if (debug) then
         write(*,*)"in advanc level ",level," after 1st implicit call"
         call valout(level,level,time,nvar,naux)
      endif


c    dont call bound here, since called for implicit update above
c    and momentum updated in ghost cells in solnUpdate
c    dont want to overwrite
c
cc     maxthreads initialized to 1 above in case no openmp
c!$    maxthreads = omp_get_max_threads()
c
c We want to do this regardless of the threading type

      if (.not. implicitSolveDone) then ! this level needs to call bound
        call system_clock(clock_startBound,clock_rate)
        call cpu_time(cpu_startBound)
!$OMP PARALLEL DO PRIVATE(j,locnew,locaux,mptr,nx,ny,mitot,
!$OMP&                    mjtot,time,levSt),
!$OMP&            SHARED(level,nvar,naux,alloc,intrat,delt,nghost,
!$OMP&                   node,rnode,numgrids,listStart,listOfGrids),
!$OMP&            SCHEDULE (dynamic,1)
!$OMP&            DEFAULT(none)
        do  j = 1, numgrids(level)
            levSt  = listStart(level)
            mptr   = listOfGrids(levSt+j-1)
            nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            mitot  = nx + 2*nghost
            mjtot  = ny + 2*nghost
            locnew = node(store1,mptr)
            locaux = node(storeaux,mptr)
            time   = rnode(timemult,mptr)
c
            call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1                 alloc(locaux),naux)

        end do
!$OMP END PARALLEL DO
        call system_clock(clock_finishBound,clock_rate)
        call cpu_time(cpu_finishBound)
        timeBound = timeBound + clock_finishBound - clock_startBound
        timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound
      endif
c
c
c save coarse level values if there is a finer level for wave fixup
      if (level+1 .le. mxnest) then
         if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
         endif
      endif
c
      time = rnode(timemult,lstart(level))
c      call fgrid_advance(time,delt)
      
      dtlevnew = rinfinity
      cfl_level = 0.d0    !# to keep track of max cfl seen on each level

      if (.not. topo_finalized) then
         call topo_update(time)
          if (debug) then
             write(*,*)"after topoupdate level ",level
             call valout(level,level,time,nvar,naux)
          endif
         endif
c 
      call system_clock(clock_startStepgrid,clock_rate)
      call cpu_time(cpu_startStepgrid)
        
!$OMP PARALLEL DO 
!$OMP&            PRIVATE(j,mptr,nx,ny,mitot,mjtot)  
!$OMP&            PRIVATE(mythread,dtnew,levSt)
!$OMP&            SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat)
!$OMP&            SHARED(nghost,intratx,intraty,hx,hy,naux,listsp)
!$OMP&            SHARED(node,rnode,dtlevnew,numgrids)
!$OMP&            SHARED(listStart,listOfGrids)
!$OMP&            SCHEDULE (DYNAMIC,1)
!$OMP&            DEFAULT(none)
      do  j = 1, numgrids(level)
          levSt  = listStart(level)
          mptr   = listOfGrids(levSt+j-1)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
c
c         call to set up for finite volume update
          call par_advanc(mptr,mitot,mjtot,nvar,naux,dtnew)
!$OMP CRITICAL (newdt)
          !write(*,555) mptr, dtnew
555       format(" grid ",i5," has new dt ",e15.7)
          dtlevnew = dmin1(dtlevnew,dtnew)
!$OMP END CRITICAL (newdt)    

      end do
!$OMP END PARALLEL DO
      call system_clock(clock_finishStepgrid,clock_rate)
      call cpu_time(cpu_finishStepgrid)
      timeStepgrid = timeStepgrid + 
     &               clock_finishStepgrid-clock_startStepgrid
      timeStepgridCPU = timeStepgridCPU + 
     &               cpu_finishStepgrid-cpu_startStepgrid

      if (debug) then
       write(*,*)"after bound/sw step before 2nd impl call level ",level
       call valout(level,level,time,nvar,naux)
      endif
c
c
c     set up for post SWE  Bouss step, for updated corrections
c     needed for finer grid interp, but actually dont update momentum
c     don't need psi on finest level after stepping either 
      if (ibouss .gt. 0) then
       if (level .ge. minLevelBouss .and. level .le. maxLevelBouss .and.
     &   level.lt.mxnest .and. startWithBouss) then
         minfo =>  matrix_info_allLevs(level) 
         time = rnode(timemult,lstart(level))
         doUpdate = .false. ! only need psi, for
         call implicit_update(nvar,naux,level,minfo%numBoussCells,
     &                        doUpdate,time)
       endif
      endif
c
      call system_clock(clock_finish,clock_rate)
      call cpu_time(cpu_finish)
      tvoll(level) = tvoll(level) + clock_finish - clock_start
      tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start

      if (debug) then
         write(*,*)"in advanc level ",level," after 2nd implicit call"
         call valout(level,level,time,nvar,naux)
      endif

c
#ifdef WHERE_AM_I
      write(*,*)"ending   advance with level ",level
#endif
      return
      end
c
c -------------------------------------------------------------
c
       subroutine prepgrids(listgrids,num, level)

       use amr_module
       implicit double precision (a-h,o-z)
       integer listgrids(num)

       mptr = lstart(level)
       do j = 1, num
          listgrids(j) = mptr
          mptr = node(levelptr, mptr)
       end do

      if (mptr .ne. 0) then
         write(*,*)" Error in routine setting up grid array "
         stop
      endif

      return
      end

c
c --------------------------------------------------------------
c
      subroutine par_advanc (mptr,mitot,mjtot,nvar,naux,dtnew)
c
      use amr_module
      use gauges_module, only: update_gauges, num_gauges
      implicit double precision (a-h,o-z)


      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/

      double precision fp(nvar,mitot,mjtot),fm(nvar,mitot,mjtot)
      double precision gp(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)


c
c  :::::::::::::: PAR_ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate this grid. grids are done in parallel.
c  extra subr. used to allow for stack based allocation of
c  flux arrays. They are only needed temporarily. If used alloc
c  array for them it has too long a lendim, makes too big
c  a checkpoint file, and is a big critical section.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      delt  =  possk(level)
      nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      time  = rnode(timemult,mptr)

!$         mythread = omp_get_thread_num()

      locold = node(store2, mptr)
      locnew = node(store1, mptr)

c
c  copy old soln. values into  next time steps soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
!         if (level .lt. mxnest) then
!            ntot   = mitot * mjtot * nvar
!cdir$ ivdep
!            do i = 1, ntot
!              alloc(locold + i - 1) = alloc(locnew + i - 1)
!            enddo
!         endif
c
      xlow = rnode(cornxlo,mptr) - nghost*hx
      ylow = rnode(cornylo,mptr) - nghost*hy

!$OMP CRITICAL(rv)
      rvol = rvol + nx * ny
      rvoll(level) = rvoll(level) + nx * ny
!$OMP END CRITICAL(rv)

c     Call b4step2 here so that time dependent arrays can be filled properly
      locaux = node(storeaux,mptr)
      call b4step2(nghost, nx, ny, nvar, alloc(locnew), 
     &             rnode(cornxlo,mptr), rnode(cornylo,mptr), hx, hy, 
     &             time, dt, naux, alloc(locaux))
c
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         locsvq = locsvf + nvar*lenbc
         locx1d = locsvq + nvar*lenbc
         call qad(alloc(locnew),mitot,mjtot,nvar,
     1            alloc(locsvf),alloc(locsvq),lenbc,
     2            intratx(level-1),intraty(level-1),hx,hy,
     3            naux,alloc(locaux),alloc(locx1d),delt,mptr)
      endif

c        # See if the grid about to be advanced has gauge data to output.
c        # This corresponds to previous time step, but output done
c        # now to make linear interpolation easier, since grid
c        # now has boundary conditions filled in.

c     should change the way print_gauges does io - right now is critical section
c     NOW changed, mjb 2/6/2015.
c     NOTE that gauge subr called before stepgrid, so never get
c     the very last gauge time at end of run.

      if (num_gauges > 0) then
           call update_gauges(alloc(locnew:locnew+nvar*mitot*mjtot),
     .                       alloc(locaux:locaux+naux*mitot*mjtot),
     .                       xlow,ylow,nvar,mitot,mjtot,naux,mptr)
           endif

c
      call stepgrid(alloc(locnew),fm,fp,gm,gp,
     2            mitot,mjtot,nghost,
     3            delt,dtnew,hx,hy,nvar,
     4            xlow,ylow,time,mptr,naux,alloc(locaux),.true.)

      if (node(cfluxptr,mptr) .ne. 0)
     2   call fluxsv(mptr,fm,fp,gm,gp,
     3               alloc(node(cfluxptr,mptr)),mitot,mjtot,
     4               nvar,listsp(level),delt,hx,hy)
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
         locsvf = node(ffluxptr,mptr)
         call fluxad(fm,fp,gm,gp,
     2               alloc(locsvf),mptr,mitot,mjtot,nvar,
     4               lenbc,intratx(level-1),intraty(level-1),
     5               nghost,delt,hx,hy)
      endif

          rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
c
      return
      end
c
c -----------------------------------------------------------
c
      subroutine copyNew2Old(level,nvar)

      use amr_module
      implicit none

      integer, intent(in) :: level, nvar
      integer:: i,j, levSt, mptr,nx,ny,mitot,mjtot,locold, locnew,ntot

      do  j = 1, numgrids(level)
          levSt  = listStart(level)
          mptr   = listOfGrids(levSt+j-1)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locold = node(store2, mptr)
          locnew = node(store1, mptr)

          ntot   = mitot * mjtot * nvar
!cdir$ ivdep
          do i = 1, ntot
            alloc(locold + i - 1) = alloc(locnew + i - 1)
          enddo
      enddo

      return
      end
