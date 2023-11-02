c :::::::::::::::::::: FLGLVL :::::::::::::::::::::::::::::::::
c
!> \callgraph
!! \callergraph
!! Controls the error estimation/flagging bad pts. for
!! an entire level of grids.  returns pointer into alloc
!! where the (x,y) coordinations of the flagged pts. are.
!!
!! \param nvar number of equations for the system
!! \param naux number of auxiliary variables
!! \param lcheck level to be flagged
!! \param nxypts number of flagged points in total
!! \param index starting index (memory address) in alloc of the flagged points (which occupy 2*nxypts locations)
!! \param lbase  base AMR level for current refinement, which stays
!! fixed. Note that **lbase** is always less or equal to **lcheck**
!! \param npts Number of unique flagged cells (after removing
!! duplicates)
!! \param[in] start_time start time of current simulation
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c -----------------------------------------------------------
c
      subroutine flglvl2(nvar,naux,lcheck,nxypts,index,lbase,npts,
     &                   start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)
      integer(kind=8) :: clock_start, clock_finish, clock_rate
c
c
c
#ifdef WHERE_AM_I
      write(*,*)"starting flglvl2"
#endif
      nxypts = 0
      numbad = 0


c     flag arrays- based on either spatial gradients (and/or user defined 
c                  criteria),  or Richardson error estimation,
c                  or user-defined regions
    
      call system_clock(clock_start,clock_rate)
      call flagger(nvar,naux,lcheck,start_time)
      call system_clock(clock_finish,clock_rate)
      timeFlagger = timeFlagger + clock_finish - clock_start


c     buffer the flagged cells (done for each grid patch of flags)
c     also project flags from finer levels onto this level to ensure
c     proper nesting. Finally compute proper domain for each patch
      call system_clock(clock_start,clock_rate)
      call bufnst2(nvar,naux,numbad,lcheck,lbase) 
      call system_clock(clock_finish,clock_rate)
      timeBufnst = timeBufnst + clock_finish - clock_start

      nxypts = nxypts + numbad
c
c  colate flagged pts into flagged points array
c  new version needs to check for proper nesting at this point
c  also needs to sort,  so can remove duplicates.
c
      if (nxypts .gt. 0) then  
c        build domain flags for each grid at level lcheck, instead of
c        previous approach using domain flags over entire  domain
c         call domgrid(lbase,lcheck)   ! will need since there are flagged pts  NOW IN BUFNST2
c
c in new version, there are bad cells but nxypts isnt true count any longer
c since there are duplicates, and proper nesting not yet checked
           index = igetsp(2*nxypts)
           call colate2(alloc(index),nxypts,lcheck,npts,lbase)
      else 
         npts = 0  !npts is number of unique flagged points after removing duplicates
         call freeFlags(lcheck)   ! otherwise storage freed in colate2. perhaps always do it here
      endif

#ifdef WHERE_AM_I
      write(*,*)"ending   flglvl2"
#endif
      return
      end
c
c ---------------------------------------------------------------------------------
c
       subroutine freeFlags(lcheck)

       use amr_module
       implicit double precision (a-h, o-z)

       mptr = lstart(lcheck)
 10          continue
             locamrflags = node(storeflags,mptr)
             locdomflags = node(domflags_base,mptr)
             locdom2 = node(domflags2,mptr)

             nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
             ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
             mbuff = max(nghost,ibuff+1)
             mibuff = nx + 2*mbuff
             mjbuff = ny + 2*mbuff

             ibytesPerDP = 8
             nwords = (mibuff*mjbuff)/ibytesPerDP+1
             call reclam(locdomflags, nwords)
             call reclam(locdom2, nwords)
             call reclam(locamrflags,mibuff*mjbuff)

        mptr = node(levelptr, mptr)
        if (mptr .ne. 0) go to 10

       return
       end
