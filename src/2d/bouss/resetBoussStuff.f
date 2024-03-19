c
c ---------------------------------------------------
c
       subroutine resetBoussStuff(naux)

       use amr_module
       use bouss_module
       implicit double precision (a-h,o-z)
       type(matrix_levInfo),  pointer :: minfo

       do lev = max(1,minLevelBouss), min(maxLevelBouss,lfine)
          numBoussGrids = 0
          minfo => matrix_info_allLevs(lev)
          do ng = 1, numgrids(lev)
             levSt = listStart(lev)
             mptr  = listOfGrids(levSt+ng-1)
             nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
             ny = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
             nb = node(bpatchNum,mptr)
             numBoussGrids = numBoussGrids + 1
          end do
          if (numBoussGrids .gt. 0) then
             allocate(minfo%matrix_indices(numBoussGrids))
             minfo%numBoussGrids = numBoussGrids
             call setMatrixIndex(lev)
             call setBoussFlag(lev,naux)
          endif
       end do

!      initialize for solver
      if (isolver .eq. 2) then
        newGrids = .true.
#ifdef HAVE_PARDISO
        do lev = 1, mxnest
           if (lev.le.maxLevelBouss .and. lev.ge.minLevelBouss) then
                call pardisoinit(pt(1,lev),mtype,solver,iparm,dparm,
     &                           error)
                write(*,*)"restrt initialized pardiso for lev ",lev
           endif
        end do
#endif
      endif

      if (isolver .eq. 3) then
         newGrids = .true.
      endif

      return
      end
