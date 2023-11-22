c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,cut,nstart,vtime,time,naux,start_time,
     &                rest,dt_max)
c
      use geoclaw_module
      use refinement_module, only: varRefTime
      use amr_module
      use topo_module, only: dt_max_dtopo, num_dtopo, topo_finalized,
     &                       aux_finalized, topo0work
      use gauges_module, only: setbestsrc, num_gauges
      use gauges_module, only: print_gauges_and_reset_nextLoc

      use storm_module, only: landfall, display_landfall_time
      use fgmax_module, only: FG_num_fgrids, FG_fgrids, fgrid
      use fgout_module, only: FGOUT_num_grids, FGOUT_fgrids,
     &                         fgout_write,fgout_grid, FGOUT_ttol

      implicit double precision (a-h,o-z)

      logical vtime,dumpout/.false./,dumpchk/.false./,rest,dump_final
      logical stopFound
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      integer(kind=8) :: clock_start, clock_finish, clock_rate
      integer(kind=8) :: tick_clock_finish, tick_clock_rate
      integer ifg
      character(len=128) :: time_format
      real(kind=8) cpu_start,cpu_finish
      type(fgrid), pointer :: fg
      type(fgout_grid), pointer :: fgout
      logical :: debug

c
c :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
c  main driver routine.  controls:
c        integration  of all grids.
c        error estimation / regridding
c        output counting
c        updating of fine to coarse grids

c  parameters:
c     nstop   = # of coarse grid time steps to be taken
c     iout    = output interval every 'iout' coarse time steps
c               (if 0, not used - set to inf.)
c     vtime   = true for variable timestep, calculated each coarse step
c
c  integration strategy is to advance a fine grid until it catches
c  up to the coarse grid. this strategy is applied recursively.
c  coarse grid goes first.
c
c  nsteps: used to count how number steps left for a level to be
c          integrated before it catches up with the next coarser level.
c  ncycle: counts number of coarse grid steps = # cycles.
c
c  icheck: counts the number of steps (incrementing by 1
c          each step) to keep track of when that level should
c          have its error estimated and finer levels should be regridded.
c ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::
c


      ncycle         = nstart
      call setbestsrc()     ! need at very start of run, including restart
      if (iout .eq. 0) then
c        # output_style 1 or 2
         iout  = iinfinity
         nextout = 0
         if (nout .gt. 0) then
            nextout = 1
            if (nstart .gt. 0) then
c              # restart: make sure output times start after restart time
               do ii = 1, nout
                 if (tout(ii) .gt. time) then
                   nextout = ii
                   go to 2
                 endif
               end do
  2         continue
            endif
         endif
      endif

      nextchk = 1
      if ((nstart .gt. 0) .and. (abs(checkpt_style).eq.2)) then
c        if this is a restart, make sure chkpt times start after restart time
         do ii = 1, nchkpt
           if (tchk(ii) .gt. time) then
              nextchk = ii
              go to 3
              endif
           enddo
  3      continue
         endif

      tlevel(1)      = time

      do 5 i       = 2, mxnest
       tlevel(i) = tlevel(1)
 5     continue

c
c  ------ start of coarse grid integration loop. ------------------
c
 20   if (ncycle .ge. nstop .or. time .ge. tfinal) goto 999

      if (nout .gt. 0) then
          if (nextout  .le. nout) then
             outtime       = tout(nextout)
          else
             outtime       = rinfinity
          endif
      else
          outtime = tfinal
      endif

      if (nextchk  .le. nchkpt) then
         chktime       = tchk(nextchk)
      else
         chktime       = rinfinity
      endif

      dumpout = .false.  !# may be reset below

      if (time.lt.outtime .and. time+1.001*possk(1) .ge. outtime) then
c        ## adjust time step  to hit outtime exactly, and make output
c        #  apr 2010 mjb: modified so allow slightly larger timestep to
c        #  hit output time exactly, instead of taking minuscule timestep
c        #  should still be stable since increase dt in only 3rd digit.
         oldposs = possk(1)
         possk(1) = outtime - time
c        write(*,*)" old possk is ", possk(1)
         diffdt = oldposs - possk(1)  ! if positive new step is smaller


         if (.false.) then  
            write(*,122) diffdt,outtime  ! notify of change
 122        format(" Adjusting timestep by ",e10.3,
     .             " to hit output time of ",e13.6)
c           write(*,*)" new possk is ", possk(1)
            if (diffdt .lt. 0.) then ! new step is slightly larger
              pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
              write(*,123) pctIncrease
 123          format(" New step is ",e9.2," % larger.",
     .               "  Should still be stable")
              endif
            endif


         do i = 2, mxnest
            possk(i) = possk(i-1) / kratio(i-1)
            enddo
         if (nout .gt. 0) then
            nextout = nextout + 1
            dumpout = .true.
            endif
      endif


      if (time.lt.chktime .and. time + possk(1) .ge. chktime) then
c        ## adjust time step  to hit chktime exactly, and do checkpointing
         possk(1) = chktime - time
         do 13 i = 2, mxnest
 13         possk(i) = possk(i-1) / kratio(i-1)
         nextchk = nextchk + 1
        dumpchk = .true.
      else
        dumpchk = .false.
      endif

c
      level        = 1
      ntogo(level) = 1
      dtnew(1:maxlv) = rinfinity
C       do i = 1, maxlv
C          dtnew(i)  = rinfinity
C       enddo

c     We should take at least one step on all levels after any
c     moving topography (dtopo) has been finalized to insure that
c     all aux arrays are consistent with the final topography.
c     The variable aux_finalized is incremented so that we can check
c     if this is true by checking if aux_finalized == 2 elsewhere in code.

      if (aux_finalized .eq. 1 .and. num_dtopo > 0) then
c         # this is only true once, and only if there was moving topo
          deallocate(topo0work)
          endif 
      if (topo_finalized .and. (aux_finalized .lt. 2)) then
          aux_finalized = aux_finalized + 1
          endif

    
c
c     ------------- regridding  time?  ---------
c
c check if either
c   (i)  this level should have its error estimated before being advanced
c   (ii) this level needs to provide boundary values for either of
c        next 2 finer levels to have their error estimated.
c        this only affects two grid levels higher, occurs because
c        previous time step needs boundary vals for giant step.
c  no error estimation on finest possible grid level
c
 60       continue
          if (icheck(level) .ge. kcheck) then
               lbase = level
          else if (level+1 .ge. mxnest) then
               go to 90
          else if (icheck(level+1) .ge. kcheck) then
               lbase = level+1
          else if (level+2 .ge. mxnest) then
               go to 90
          else if (icheck(level+2) .ge. kcheck) then
               lbase = level+2
          else
               go to 90
          endif
          if (lbase .eq. mxnest .or. lbase .gt. lfine) go to 70
c
c regrid level 'lbase+1' up to finest level.
c level 'lbase' stays fixed.
c
          if (rprint) write(outunit,101) lbase
101       format(8h  level ,i5,32h  stays fixed during regridding )

          call system_clock(clock_start,clock_rate)
          call cpu_time(cpu_start)
          debug = .false.
          if (debug) then
              write(*,*)" before regrid lbase ",lbase," time ",
     &                    tlevel(lbase)
              !call valout(lbase+1,lfine,time,nvar,naux)
              call valout(lbase,lfine,time,nvar,naux)
          endif
          call regrid(nvar,lbase,cut,naux,start_time)

          ! output new grids for debugging:
          if (debug)
     .    call valout(lbase+1,lfine,tlevel(lbase+1),nvar,naux)

          if (debug) then
              write(*,*)" after regrid at time ",tlevel(lbase)
              !call valout(lbase+1,lfine,time,nvar,naux)
              call valout(lbase,lfine,time,nvar,naux)
          endif
          call system_clock(clock_finish,clock_rate)
          call cpu_time(cpu_finish)
          timeRegridding = timeRegridding + clock_finish - clock_start
          timeRegriddingCPU=timeRegriddingCPU+cpu_finish-cpu_start

          call setbestsrc()     ! need at every grid change
c         call conck(1,nvar,naux,time,rest)
c         call outtre(lstart(lbase+1),.true.,nvar,naux)
c note negative time to signal regridding output in plots
c         call valout(lbase,lfine,-tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),printout,nvar,naux)
          endif
 70       continue
          do 80  i  = lbase, lfine
 80          icheck(i) = 0
          do 81  i  = lbase+1, lfine
 81          tlevel(i) = tlevel(lbase)
c
c          MJB: modified to check level where new grids start, which is lbase+1
          !if (verbosity_regrid.ge.lbase+1) then
          if (.false.) then  ! don't need to print these every time
                 do levnew = lbase+1,lfine
                     write(6,1006) intratx(levnew-1),intraty(levnew-1),
     &                             kratio(levnew-1),levnew
 1006                format('   Refinement ratios...  in x:', i3, 
     &                 '  in y:',i3,'  in t:',i3,' for level ',i4)
                 end do

              endif

c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue

          ! Only compute the mininum of fg%levelmax if necessary,
          ! and only once per time step on each level, not on each patch:
          do ifg=1,FG_num_fgrids
              fg => FG_fgrids(ifg)
              if (fg%min_level_check == mxnest) then
                  fg%min_levelmax_checked = mxnest
                else
                  fg%min_levelmax_checked = minval(fg%levelmax)
                endif
                enddo
c 
c         ! time after step:
          timenew = tlevel(level)+possk(level)

          ! is it time to update fgmax grids on any level?
          ! if so, set fg%update_now(level) to .true.
          do ifg=1,FG_num_fgrids
              fg => FG_fgrids(ifg)
              if ((timenew >= fg%tstart_max) .and. 
     &            (timenew <= fg%tend_max) .and.
     &            (timenew >= fg%t_last_updated(level)  
     &                         + fg%dt_check) .and.
     &            (level >= fg%min_level_check) .and.
     &            (level >= fg%min_levelmax_checked)) then 

                      fg%update_now(level) = .true.
                      fg%t_last_updated(level) = timenew
                      !write(6,660) level, timenew
 660                  format('+++ update fgmax on level ',i1,
     &                       ' at t = ',f12.2)
              else
                      fg%update_now(level) = .false.
                      endif
              enddo


          call advanc(level,nvar,dtlevnew,vtime,naux)


c Output time info
          timenew = tlevel(level)+possk(level)
          time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ',e12.6)"
          if (display_landfall_time) then
c           Convert time to days            
            timenew = timenew / (3.6d3 * 24d0)
            time_format = "(' AMRCLAW: level ',i2,'  CFL = ',e8.3," //
     &                  "'  dt = ',e10.4,  '  final t = ', f5.2)"
          end if
          if (tprint) then
              write(outunit, time_format) level, cfl_level, 
     &                                    possk(level), timenew
          endif
          if (method(4).ge.level) then
              print time_format, level, cfl_level, possk(level), timenew
          endif

c        # to debug individual grid updates...
c        call valout(level,level,time,nvar,naux)
c
c done with a level of integration. update counts, decide who next.
c
          ntogo(level)  = ntogo(level) - 1
          dtnew(level)  = dmin1(dtnew(level),dtlevnew)
          tlevel(level) = tlevel(level) + possk(level)
          icheck(level) = icheck(level) + 1
c
          if (level .lt. lfine) then
             level = level + 1
c            #  check if should adjust finer grid time step to start wtih
             if (((possk(level-1) - dtnew(level-1))/dtnew(level-1)) .gt.
     .            .05) then
                dttemp = dtnew(level-1)/kratio(level-1)
                ntogo(level) = (tlevel(level-1)-tlevel(level))/dttemp+.9
              else
                ntogo(level) = kratio(level-1)
              endif
             possk(level) = possk(level-1)/ntogo(level)
             go to 60
          endif

c         write(6,*) '+++ in tick, done with level',level,
c    &               ' tlevel = ',tlevel(1:level)

c         When we reach here, we are done with grid patches at 
c         the finest level with grids present at this time

c     ---- fgout output ----
c     This used to be done in stepgrid.f, but only needs to be done
c     after all patches at finest level have been advanced.

      tc0 = tlevel(level)  ! current time on finest level present
      !write(6,*) '+++ tick: tc0 = ',tc0
      
      do ng=1,FGOUT_num_grids
        fgout => FGOUT_fgrids(ng)
        
        do ioutfg=1,fgout%num_output
            if (fgout%output_frames(ioutfg) == -1) then
                ! this time not yet written out
                if (fgout%output_times(ioutfg) < 
     &                 tc0+FGOUT_ttol) then
                     toutfg = fgout%output_times(ioutfg)
c                     write(6,*) '+++ tick call fgrid_out, frame, t: ',
c     &                          ioutfg,toutfg
                     call fgout_write(fgout,toutfg,ioutfg)
                     fgout%output_frames(ioutfg) = ioutfg
                     fgout%next_output_index = ioutfg+1
                endif
            endif
        enddo
      enddo
c     ---- end fgout output ----
c
 105      if (level .eq. 1) go to 110
              if (ntogo(level) .gt. 0) then
c                same level goes again. check for ok time step
 106             if ((possk(level)-dtnew(level))/dtnew(level)
     .                .gt. .05)  then

                    write(6,601) level, time
 601                format(" ***adjusting timestep for level ", i3,
     &                     " at t = ",d16.6)
                    print *,"    old ntogo dt",ntogo(level),possk(level)

c                   adjust time steps for this and finer levels

                    ! try computing ntogo properly (worked better in BoussDev)
                    ! (old way was to repeatedly increment by 1)
                    new_ntogo = ceiling(((tlevel(level-1)
     &                              -tlevel(level)) / dtnew(level)))
                    new_ntogo = min(new_ntogo, ntogo(level)+10)
                    ntogo(level) = new_ntogo
                    possk(level) = (tlevel(level-1)-tlevel(level))
     &                             / ntogo(level)
                    write(*,*) "    NEW ntogo dt ",ntogo(level),
     &                         possk(level)
                    if (varRefTime) then
                      kratio(level-1) = ceiling(possk(level-1)
     &                                  / possk(level))
                    endif

                    go to 106
                 endif

                 if (ntogo(level) .gt. 100) then
                     write(6,*) "**** Too many dt reductions ****"
                     write(6,*) "**** Stopping calculation   ****"
                     write(6,*) "**** ntogo = ",ntogo(level)
                     write(6,1006) intratx(level-1),intraty(level-1),
     &                             kratio(level-1),level
 603                 format("**** Writing extra output frames for ",
     &                      "debugging at level-1 =",i3,
     &                      " and level =",i3)
                     write(6,603) level-1, level
                     if (num_gauges .gt. 0) then
                        do ii = 1, num_gauges
                           call print_gauges_and_reset_nextLoc(ii)
                        end do
                     endif
                     call valout(level-1,level-1,tlevel(level-1),
     &                           nvar,naux)                       
                     call valout(level,level,tlevel(level),nvar,naux)
                     !call outtre(lstart(level),.true.,nvar,naux)
                     stop
                 endif

                 go to 60
              else
                 level = level - 1
                 call system_clock(clock_start,clock_rate)
                 call update(level,nvar,naux)
                 call system_clock(clock_finish,clock_rate)
                 timeUpdating=timeUpdating+clock_finish-clock_start
              endif
          go to 105
c
c  --------------one complete coarse grid integration cycle done. -----
c
c      time for output?  done with the whole thing?
c
 110      continue
          time    = time   + possk(1)
          ncycle  = ncycle + 1
          call conck(1,nvar,naux,time,rest)


      if ( .not.vtime) goto 201

        ! Adjust time steps if variable time step and/or variable
        ! refinement ratios in time
        if (.not. varRefTime) then
          ! find new dt for next cycle (passed back from integration routine).
           do 115 i = 2, lfine
             ii = lfine+1-i
             dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*kratio(ii))
 115       continue
           possk(1) = dtnew(1)
           do 120 i = 2, mxnest
 120         possk(i) = possk(i-1) / kratio(i-1)
        else  ! since refinement ratio in time can change need to set new timesteps in different order
c             ! use same alg. as when setting refinement when first make new fine grids
          dtnew(1) = min(dtnew(1),dt_max)
          if ((num_dtopo>0).and.(topo_finalized.eqv..false.)) then
              dtnew(1) = min(dtnew(1),dt_max_dtopo)
          endif

          possk(1) = dtnew(1)
          do 125 i = 2, lfine
             if (dtnew(i)  .gt. possk(i-1)) then
               kratio(i-1) = 1  ! cant have larger timestep than parent level
               possk(i)    = possk(i-1)
            else
               kratio(i-1) = ceiling(possk(i-1)/dtnew(i))  ! round up for stable integer ratio
               possk(i)    = possk(i-1)/kratio(i-1)        ! set exact timestep on this level
           endif
 125    continue


      endif

 201  if ((abs(checkpt_style).eq.3 .and. 
     &      mod(ncycle,checkpt_interval).eq.0) .or. dumpchk) then
               if (num_gauges .gt. 0) then
                  do ii = 1, num_gauges
                     call print_gauges_and_reset_nextLoc(ii)
                  end do
               endif
               call check(ncycle,time,nvar,naux)
               dumpchk = .true.
       endif

       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
         call valout(1,lfine,time,nvar,naux)
         if (abs(checkpt_style).eq.4) then
            call check(ncycle,time,nvar,naux)
            dumpchk = .true.
         endif
         if (printout) call outtre(mstart,.true.,nvar,naux)
         if (num_gauges .gt. 0) then
            do ii = 1, num_gauges
               call print_gauges_and_reset_nextLoc(ii)
            end do
         endif
       endif

       ! new STOP feature to do immediate checkpt and exit
       inquire(FILE="STOP",exist=stopFound) 
          if (stopFound) then
          write(*,*)"STOP file found. Checkpointing and Stopping"
          write(*,*)"REMEMBER to remove file before restarting"
          write(outunit,*)"STOP file found. Checkpointing and Stopping"
          write(outunit,*)"REMEMBER to remove file before restarting"
          call check(ncycle,time,nvar,naux)
          stop
       endif

      go to 20
c
999   continue

c
c  # computation is complete to final time or requested number of steps
c
       if (ncycle .ge. nstop .and. tfinal .lt. rinfinity) then
c         # warn the user that calculation finished prematurely
          write(outunit,102) nstop
          write(6,102) nstop
  102     format('*** Computation halted after nv(1) = ',i8,
     &           '  steps on coarse grid')
          endif
c
c  # final output (unless we just did it above)
c
      dump_final = ((iout.lt.iinfinity) .and. (mod(ncycle,iout).ne.0))
      if (.not. dumpout) then
          if (nout > 0) then
              dump_final = (tout(nout).eq.tfinal)
              endif
          endif
      
      if (dump_final) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
           if (num_gauges .gt. 0) then
              do ii = 1, num_gauges
                 call print_gauges_and_reset_nextLoc(ii)
              end do
           endif
           dumpout = .true.
      endif

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c
      call system_clock(tick_clock_finish,tick_clock_rate)
      call cpu_time(tick_cpu_finish)
      timeTick = timeTick + tick_clock_finish - tick_clock_start 
      timeTickCPU = timeTickCPU + tick_cpu_finish - tick_cpu_start 


c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c
      ! write gauges first in case lagrangian gauge x,y needed by check         
      if (num_gauges .gt. 0) then
         if (.not. dumpout) then
             ! normally this was done already, unless nout = 0
             do ii = 1, num_gauges
                call print_gauges_and_reset_nextLoc(ii)
             end do
         endif
      endif
      
      if (checkpt_style .ne. 0) then  ! want a chckpt
         ! check if just did it so dont do it twice
         if (.not. dumpchk) call check(ncycle,time,nvar,naux)
      endif

      write(6,*) "Done integrating to time ",time
      return
      end
