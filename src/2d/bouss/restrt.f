
c
c ---------------------------------------------------------
c
      subroutine restrt(nsteps,time,nvar,naux)
c
      use amr_module
      use fgmax_module
      use bouss_module
      use gauges_module, only: gauges, num_gauges
      use topo_module, only: aux_finalized
      implicit double precision (a-h,o-z)
      logical   ee
 
      logical foundFile
      dimension intrtx(maxlv),intrty(maxlv),intrtt(maxlv)
      type(fgrid), pointer :: fg

      integer :: num_gauges_previous, i, ii, previous_gauge_num
      integer matlabu0
c
c :::::::::::::::::::::::::::: RESTRT ::::::::::::::::::::::::::::::::
c read back in the check point files written by subr. check.
c
c some input variables might have changed, and also the
c alloc array could have been written with a smaller size at checkpoint
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     !! Now allow user-specified file name !!
c     rstfile  = 'restart.data'

      ! If checkpt_style < 0 then alternating between two checkpoint files.
      ! Set which one to use for the first checkpoint after this restart.
      ! Set check_a to .true. unless fort.chkaaaaa is the file being read.
      ! When alternating checkpoint files used, this keeps proper sequence going
      ! otherwise (checkpt_style > 0) check_a is not used elsewhere.
      check_a = .not. (rstfile == 'fort.chkaaaaa')

      write(6,*) 'Attempting to restart computation using '
      write(6,*) '  checkpoint file: ',trim(rstfile)
      inquire(file=trim(rstfile),exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did not find checkpoint file!"
        stop
      endif
      open(rstunit,file=trim(rstfile),status='old',form='unformatted')
      rewind rstunit

      matlabu0 = matlabu ! value set by amr2, 0 if output_t0, else 1

      !read(rstunit) lenmax,lendim,isize
      !! new version has flexible node size, so need to read current size maxgr
      read(rstunit) lenmax,lendim,isize,maxgr

c     # need to allocate for dynamic memory:
      call restrt_alloc(isize)
      call restrt_nodes(maxgr)

      read(rstunit) (alloc(i),i=1,lendim)
      read(rstunit) hxposs,hyposs,possk,icheck
      read(rstunit) lfree,lenf
      read(rstunit) rnode,node,lstart,newstl,listspStart,listsp,tl,
     1       ibuf,mstart,ndfree,ndfree_bnd,lfine,iorder,mxnold,
     2       intrtx,intrty,intrtt,iregsz,jregsz,
     2       iregst,jregst,iregend,jregend,
     3       numgrids,kcheck1,nsteps,time,
     3       matlabu,aux_finalized,
     4       minLevelBouss_orig,maxLevelBouss_orig,
     4       boussMinDepth_orig,startBoussTime_orig
      read(rstunit) avenumgrids, iregridcount,
     1              evol,rvol,rvoll,lentot,tmass0,cflmax,
     2              tvoll,tvollCPU,timeTick,timeTickCPU,
     3              timeStepgrid,timeStepgridCPU,
     3              timeLinSolve,timeLinSolveCPU,
     3              timePrepBuild,timePrepBuildCpu,
     4              timeBound,timeBoundCPU,
     5              timeRegridding,timeRegriddingCPU,
     6              timeValout,timeValoutCPU

!     WRITE(*,*)"REMEMBER:  RESETTING Solver TIMERS For special tests"
!     timeLinSolve = 0
!     timeLinSolveCPU = 0.d0
!     timePrepBuild = 0
!     timePrepBuildCPU = 0.d0

      if (mxnest .gt. mxnold) then
        do ifg = 1, FG_num_fgrids
          fg => FG_fgrids(ifg)
          read(rstunit) fg%levelmax
          read(rstunit) fg%auxdone(1:mxnold) 
          read(rstunit) fg%x,fg%y,fg%valuemax,fg%tmax,
     &       fg%arrival_time,fg%aux(1:mxnold,:,:),
     &       fg%t_last_updated(1:mxnold)
        end do
      else
        do ifg = 1, FG_num_fgrids
          fg => FG_fgrids(ifg)
          read(rstunit) fg%levelmax
          read(rstunit) fg%auxdone(1:mxnest)
          read(rstunit) fg%x,fg%y,fg%valuemax,fg%tmax,
     &         fg%arrival_time,fg%aux,fg%t_last_updated
        end do
      endif

c     Check for any lagrangian gauges, and if present reset x,y location:

      read(rstunit) num_gauges_previous

      ! Note that it now works (?) to add new gauges for a restart that
      ! were not previously there, but the old ones should still be there
      ! with the same gauge numbers and order as before.
      ! Also note that set_gauges is now called before restrt, so 
      ! we have access to old gauge data, and can overwrite location
      ! read from gauges.data for lagrangian gauges that have moved.

      if (num_gauges_previous > num_gauges) then
         write(6,*) '*** Number of gauges reduced, expect problems!'
      else
         do ii = 1, num_gauges
            previous_gauge_num = gauges(ii)%gauge_num
            read(rstunit) gauges(ii)%gauge_num,
     &                    gauges(ii)%t_last_written,
     &                    gauges(ii)%x_last_written,
     &                    gauges(ii)%y_last_written
            if (gauges(ii)%gauge_num .ne. previous_gauge_num) then
                write(6,*) '*** Gauge number has changed!'
                write(6,*) ii,gauges(ii)%gauge_num,previous_gauge_num
                ! stop?
            endif
            if (gauges(ii)%gtype == 2) then
                ! lagrangian gauge, reset x,y location:
                gauges(ii)%x = gauges(ii)%x_last_written
                gauges(ii)%y = gauges(ii)%y_last_written
                if ((gauges(ii)%t_last_written .ne. NEEDS_TO_BE_SET)
     &             .and. (gauges(ii)%t_last_written .ne. time)) then
                    ! Note this always happens now since gauge is
                    ! not written at final time.  Remove this warning
                    ! until we fixt that?
                    write(6,602) time, gauges(ii)%gauge_num, 
     &                       gauges(ii)%t_last_written
 602                format('Restart at time ', e13.6,
     &                 '  Gauge ', i8, ' is reset from time ', e13.6)
                endif
            endif
         end do
      endif

c
      close(rstunit) 

      matlabu = matlabu - 1 + matlabu0 ! redo previous frame if output_t0

      write(outunit,100) nsteps,time
      write(6,100) nsteps,time
 100  format(/,' RESTARTING the calculation after ',i5,' steps',
     1        /,'  (time = ',e15.7,')')
c
c     error checking that refinement ratios have not changed
c     ### new feature: when using variable refinement in time
c     ### (varRefTime = T) the time ratios are allowed to be different
c     ###  (since they are ignored and calc. on the fly)
c     ### This is not checked for here, since the same amr2.f is now
c         used for geoclaw too, and varRefTime is only available in geoclaw.
c
      do i = 1, min(mxnold-1,mxnest-1)
        if ( (intratx(i) .ne. intrtx(i)) .or.
     .       (intraty(i) .ne. intrty(i)) ) then
c    .       (kratio(i) .ne.  intrtt(i) .and. .not. varRefTime) ) then
        write(outunit,*) 
     .  " not allowed to change existing refinement ratios on Restart"
        write(outunit,*)" Old ratios:"
        write(*,*)      " Old ratios:"
        write(outunit,903)(intrtx(j),j=1,mxnold-1)
        write(*,903)      (intrtx(j),j=1,mxnold-1)
        write(outunit,903)(intrty(j),j=1,mxnold-1)
        write(*,903)      (intrty(j),j=1,mxnold-1)
c       write(outunit,903)(intrtt(j),j=1,mxnold-1)
c       write(*,903)      (intrtt(j),j=1,mxnold-1)
 903    format(6i3)
        stop
       endif
      end do

c     if (varRefTime) then  ! reset intrat to previously saved ratios, not input ratios
      if (.true.) then  ! reset intrat to previously saved ratios, not input ratios
        do i = 1, mxnold-1
            kratio(i) = intrtt(i)
        end do
      endif

c
c adjust free list of storage in case size has changed.
c
      idif = memsize - isize
      if (idif .gt. 0) then
          lfree(lenf,1) = isize + 2
          call reclam(isize+1,idif)
      else if (idif .lt. 0) then
            write(outunit,900) isize, memsize
            write(*,900)       isize, memsize
 900        format(' size of alloc not allowed to shrink with ',/,
     .             ' restart old size ',i7,' current size  ',i7)
            stop
      endif
c
c adjust storage in case mxnest has changed - only allow it to increase,
c
       if (mxnest .eq. mxnold) go to 99

       if (mxnest .lt. mxnold) then
         if (lfine .lt. mxnest) then
             go to 99
         else
             write(outunit,901) mxnold, mxnest
             write(*,      901) mxnold, mxnest
901          format('  mxnest reduced on restart: ',/,
     &            '  old mxnest ',i4, ' new mxnest ',i4)
             write(outunit,*)" reclaiming finer levels from",
     .                mxnest+1," to ",mxnold
             write(*,*)" reclaiming finer levels from",
     .                mxnest+1," to ",mxnold
             do 95 lev = mxnest,mxnold
                mptr = lstart(lev)
                if (lev .gt. mxnest) lstart(lev) = 0   
 85             if (mptr .eq. 0) go to 95
                   if (lev .lt. mxnold) then
                    call reclam(node(cfluxptr,mptr), 5*listsp(lev))
                    node(cfluxptr,mptr) = 0
                   endif
                   nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
                   ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
                   if (lev .gt. 1) then
                      ikeep = nx/intrtx(lev-1)
                      jkeep = ny/intrty(lev-1)
                      lenbc = 2*(ikeep+jkeep)
                   endif
                   if (lev .gt. mxnest) then
                       call reclam
     .                  (node(ffluxptr,mptr),2*nvar*lenbc+naux*lenbc)
                       node(ffluxptr,mptr) = 0
                   endif
                   mitot = nx + 2*nghost
                   mjtot = ny + 2*nghost
                   if (lev .gt. mxnest) then ! if level going away take away first storage
                      call reclam(node(store1,mptr),mitot*mjtot*nvar)
                      node(store1,mptr) = 0
                      if (naux .gt. 0) then ! and aux arrays
                       call reclam(node(storeaux,mptr),mitot*mjtot*naux)
                       node(storeaux,mptr) = 0
                      endif
                   endif
                   if (lev .ge. mxnest .and. lev .lt. mxnold) then  !reclam 2nd level storage too
                      call reclam(node(store2, mptr), mitot*mjtot*nvar)
                      node(store2,mptr) = 0
                   endif
                mold   = mptr
                mptr   = node(levelptr,mptr)
                if (lev .gt. mxnest) call putnod(mold)
                go to 85
 95        end do
c              stop
c  reset lfine to new finest level
             do lev = mxnest,1,-1
                if (lstart(lev) .gt. 0) then
                  lfine = lev
                  write(*,*)" resetting finest level to ",lfine
                  go to 199
                endif 
             end do
         endif
       endif

c if mxnest has increased, add second storage loc to grids at previous mxnest 
        ee = .false.
        do 10 level = 1, mxnold
           if (icheck(level) .ge. kcheck) then
              ee = .true.
           endif
10      continue

           write(*,*)" increasing max num levels from ",mxnold,
     .                 ' to',mxnest
           write(outunit,*)" increasing max num levels from ",mxnold,
     .                 ' to',mxnest

        if (ee .and. flag_richardson) then
c           ## if Richardson used, will delay error est. 1 step til have old soln. vals
             write(*,*)" first Richardson error estimation step"
             write(*,*)" will estimate mostly spatial error "
             write(outunit,*)" first Richardson error estimation step"
             write(outunit,*)" will estimate mostly spatial error  "
         endif

c          #  add second storage location to previous mxnest level
         mptr = lstart(mxnold)
15       if (mptr .eq. 0) go to 25
            mitot = node(ndihi,mptr)-node(ndilo,mptr)+1+2*nghost
            mjtot = node(ndjhi,mptr)-node(ndjlo,mptr)+1+2*nghost
            node(store2,mptr) = igetsp(mitot*mjtot*nvar)
            mptr = node(levelptr,mptr)
            go to 15
25       continue
c
c          # add new info. to spatial and counting arrays
 99        level = lfine + 1
           rrk = dble(kratio(lfine))
35         if (level .gt. mxnest) go to 45
             hxposs(level) = hxposs(level-1) / dble(intratx(level-1))
             hyposs(level) = hyposs(level-1) / dble(intraty(level-1))
             possk (level) = possk (level-1) / rrk
             iregsz(level) = iregsz(level-1) * intratx(level-1)
             jregsz(level) = jregsz(level-1) * intraty(level-1)
             rrk           = kratio(level)
             level         = level + 1
             go to 35
45         continue
c
c
 199   continue

c
c     save array of grids to avoid copying each advanc or update step
c     lbase is 1 here, to start building from level 1
c     only for level 1 is listStart set outside of makeGridList
c     call with lbase 0 to make level 1
      listStart(1) = 1
      call makeGridList(0)
c
c     bndry list for faster ghost cell filling
      call initBndryList()
      do level = 1, lfine
         call makeBndryList(level)
      end do
c

      return
      end
