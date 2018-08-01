!> Integrate all grids at the input **level** by one step of its delta(t)
!!
!! this includes:  
!! - setting the ghost cells 
!! - advancing the solution on the grid
!! - adjusting fluxes for flux conservation step later
! --------------------------------------------------------------
!
#include "amr_macros.H"

subroutine advanc(level,nvar,dtlevnew,vtime,naux)
    use amr_module 
    use parallel_advanc_module
    use fixedgrids_module
    use topo_module, only: topo_finalized
    use geoclaw_module
#ifdef CUDA
    use gauges_module, only: update_gauges, num_gauges
    use memory_module, only: cpu_allocate_pinned, cpu_deallocate_pinned, &
        gpu_allocate, gpu_deallocate
    use cuda_module, only: device_id, id_copy_cflux, toString, id_copy_fflux
    use cuda_module, only: wait_for_all_gpu_tasks, wait_for_stream
    use cuda_module, only: aos_to_soa_r2, soa_to_aos_r2, get_cuda_stream
    use cuda_module, only: compute_kernel_size
    use timer_module
    use cudafor
    use reflux_module, only: qad_cpu2, qad_gpu, &
        fluxad_fused_gpu, fluxsv_fused_gpu
    use cuda_module, only: grid_type
    ! use problem_para_module, only: cc, zz
#ifdef PROFILE
    use profiling_module
#endif
#endif
    implicit none


    logical vtime
    integer omp_get_thread_num, omp_get_max_threads
    integer mythread/0/, maxthreads/1/
    integer listgrids(numgrids(level))
    integer clock_start, clock_finish, clock_rate
    integer clock_startStepgrid,clock_startBound,clock_finishBound
    integer :: level, levst, mptr, locsvq
    integer :: nx, ny, mitot, mjtot, nvar, naux, lenbc, ntot
    real(CLAW_REAL) hx,hy,delt, dtnew, dtlevnew, time
    real(CLAW_REAL) cpu_start, cpu_finish
    real(CLAW_REAL) cpu_startBound, cpu_finishBound
    real(CLAW_REAL) cpu_startStepgrid, cpu_finishStepgrid

#ifdef CUDA
    integer :: locold, locnew, locaux
    integer :: i,j, id
    integer :: cudaResult
    real(CLAW_REAL) :: xlow, ylow
    real(CLAW_REAL) :: cfl_local
    type(dim3) :: numBlocks, numThreads
    real(CLAW_REAL), dimension(:,:), pointer, contiguous :: cfls
    real(CLAW_REAL), dimension(:,:), pointer, contiguous, device :: cfls_d

    ! type(grid_type), allocatable         :: grids(:)
    ! type(grid_type), allocatable, device :: grids_d(:)
    integer :: max_lenbc

#endif

    !     maxgr is maximum number of grids  many things are
    !     dimensioned at, so this is overall. only 1d array
    !     though so should suffice. problem is
    !     not being able to dimension at maxthreads


    !
    !  ::::::::::::::: ADVANC :::::::::::::::::::::::::::::::::::::::::::
    !  integrate all grids at the input  'level' by one step of its delta(t)
    !  this includes:  setting the ghost cells 
    !                  advancing the solution on the grid
    !                  adjusting fluxes for flux conservation step later
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !

    ! get start time for more detailed timing by level
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)

#ifdef PROFILE
    call take_cpu_timer("advanc", timer_advanc)
    call cpu_timer_start(timer_advanc)
#endif

!$OMP PARALLEL PRIVATE(hx, hy, delt, &
!$OMP                  j, levSt, mptr, nx, ny, mitot, mjtot, locold, locnew, locaux, time, &
!$OMP                  id, xlow, ylow, cudaResult, lenbc, locsvq, numBlocks, numThreads, &
!$OMP                  ntot, mythread, dtnew) &
!$OMP          SHARED(node, rnode, alloc) &
!$OMP          DEFAULT(SHARED)
    cudaResult = cudaSetDevice(device_id)

#ifdef PROFILE
    call startCudaProfiler("advanc level "//toString(level),level)
#endif
  
  

    hx   = hxposs(level)
    hy   = hyposs(level)
    delt = possk(level)
    !     this is linear alg.
    !     call prepgrids(listgrids,numgrids(level),level)
    !

    !$OMP MASTER 
    call system_clock(clock_startBound,clock_rate)
    call cpu_time(cpu_startBound)
    !$OMP END MASTER 


    !     maxthreads initialized to 1 above in case no openmp
    ! !$    maxthreads = omp_get_max_threads()




!! ################################################################## 
!! Filling ghost cells 
!! ##################################################################
#ifdef PROFILE
    call startCudaProfiler("bound", 24)
    call take_cpu_timer("Filling ghost cells", timer_bound)
    call cpu_timer_start(timer_bound)
#endif
    ! We want to do this regardless of the threading type
    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr   = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1,mptr)
        locaux = node(storeaux,mptr)
        time   = rnode(timemult,mptr)
        !     
        call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr, alloc(locaux),naux)

    end do
    !$OMP END DO
#ifdef PROFILE
    call cpu_timer_stop(timer_bound)
    call endCudaProfiler() ! bound
#endif

    !$OMP MASTER 
    call system_clock(clock_finishBound,clock_rate)
    call cpu_time(cpu_finishBound)
    timeBound = timeBound + clock_finishBound - clock_startBound
    timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound
    !$OMP END MASTER 



!! ##################################################################
!! Save coarse level values if there is a finer level for wave fixup
!! ##################################################################

#ifdef PROFILE
    call startCudaProfiler("saveqc", 24)
    call take_cpu_timer('saveqc', timer_saveqc)
    call cpu_timer_start(timer_saveqc)
#endif

    !$OMP MASTER 
    ! if (level+1 .le. mxnest) then
    !     if (lstart(level+1) .ne. clawpack_null) then
    !         call saveqc(level+1,nvar,naux)
    !     endif
    ! endif
    call system_clock(clock_startStepgrid,clock_rate)
    call cpu_time(cpu_startStepgrid)
    if (.not. topo_finalized) then
#ifndef CUDA
        time = rnode(timemult,lstart(level))
        call topo_update(time)
#else
        print *, "dtopo is not supported in CUDA version."
        stop
#endif
    endif
    !$OMP END MASTER 

#ifdef PROFILE
    call cpu_timer_stop(timer_saveqc)
    call endCudaProfiler() ! saveqc
#endif


    !$OMP SINGLE
    dtlevnew = rinfinity
    cfl_level = 0.d0    !# to keep track of max cfl seen on each level
    !$OMP END SINGLE

#ifdef PROFILE
    call take_cpu_timer('stepgrid', timer_stepgrid)
    call cpu_timer_start(timer_stepgrid)
    call startCudaProfiler("Stepgrid", 24)
#endif



#ifdef CUDA


    !$OMP SINGLE
    call cpu_allocate_pinned(cfls,1,SPACEDIM,1,numgrids(level))
    call gpu_allocate(cfls_d,device_id,1,SPACEDIM,1,numgrids(level))

    ! TODO: merge this with something else
    cfls_d = 0.d0
    !$OMP END SINGLE 

    ! !$OMP SINGLE 
    ! allocate(grids(numgrids(level)))
    ! allocate(grids_d(numgrids(level)))
    ! max_lenbc = 0
    ! !$OMP END SINGLE 
    !$OMP BARRIER

#ifdef PROFILE
    call take_cpu_timer('qad, advance sol. and copy old sol.', timer_gpu_loop)
    call cpu_timer_start(timer_gpu_loop)

    call take_cpu_timer('Launch qad and stepgrid_soa', timer_launch_compute_kernels)
    call cpu_timer_start(timer_launch_compute_kernels)

    call startCudaProfiler("qad, advance sol. and copy old sol.",74)
    call startCudaProfiler("Launch qad and stepgrid_soa",74)
#endif
!! ##################################################################
!! qad and stepgrid
!! ##################################################################
    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1, mptr)
        id = j

        !  copy old soln. values into  next time step's soln. values
        !  since integrator will overwrite it. only for grids not at
        !  the finest level. finest level grids do not maintain copies
        !  of old and new time solution values.

        locold = node(store2, mptr)
        locnew = node(store1, mptr)
        locaux = node(storeaux,mptr)
    
        if (level .lt. mxnest) then
            ntot   = mitot * mjtot * nvar
            do i = 1, ntot
                alloc(locold + i - 1) = alloc(locnew + i - 1)
            enddo
        endif

        !        # See if the grid about to be advanced has gauge data to output.
        !        # This corresponds to previous time step, but output done
        !        # now to make linear interpolation easier, since grid
        !        # now has boundary conditions filled in.

        !     should change the way print_gauges does io - right now is critical section
        !     no more,  each gauge has own array.

        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy
        if (num_gauges > 0) then
            call update_gauges(alloc(locnew:locnew+nvar*mitot*mjtot), &
                               alloc(locaux:locaux+nvar*mitot*mjtot), &
                               xlow,ylow,nvar,mitot,mjtot,naux,mptr)
        endif


        ! call cpu_allocate_pinned(grid_data(mptr)%ptr, &
        !         1,mitot,1,mjtot,1,nvar)
#ifdef PROFILE
        call take_cpu_timer('aos_to_soa', timer_aos_to_soa)
        call cpu_timer_start(timer_aos_to_soa)
        call startCudaProfiler("aos_to_soa", 14)
#endif
#ifdef PROFILE
        call endCudaProfiler() ! aos_to_soa
        call cpu_timer_stop(timer_aos_to_soa)
#endif


        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy
        locaux = node(storeaux,mptr)

#ifdef PROFILE
        call take_cpu_timer('memory operation', timer_memory)
        call cpu_timer_start(timer_memory)
#endif
        ! call gpu_allocate(grid_data_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(grid_data_d_copy2(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(aux_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,naux)
        ! call gpu_allocate(fms_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(fps_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(gms_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(gps_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        ! call gpu_allocate(sx_d(mptr)%ptr,device_id,1,mitot-1,1,mjtot-2,1,NWAVES)
        ! call gpu_allocate(sy_d(mptr)%ptr,device_id,1,mitot-2,1,mjtot-1,1,NWAVES)
        ! call gpu_allocate(wave_x_d(mptr)%ptr,device_id,1,mitot-1,1,mjtot-2,1,NEQNS,1,NWAVES)
        ! call gpu_allocate(wave_y_d(mptr)%ptr,device_id,1,mitot-2,1,mjtot-1,1,NEQNS,1,NWAVES)
#ifdef PROFILE
        call cpu_timer_stop(timer_memory)
#endif

        ! copy q to GPU
        !$OMP CRITICAL(launch)
        cudaResult = cudaMemcpyAsync(grid_data_d(mptr)%ptr, alloc(locnew), nvar*mitot*mjtot, cudaMemcpyHostToDevice, get_cuda_stream(id,device_id))
        cudaResult = cudaMemcpyAsync(aux_d(mptr)%ptr, alloc(locaux), naux*mitot*mjtot, cudaMemcpyHostToDevice, get_cuda_stream(id,device_id))
        
!         if (associated(fflux_hh(mptr)%ptr)) then
!             lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
!             locsvq = 1 + nvar*lenbc
! 
!             ! CPU version
!             ! cudaResult = cudaMemcpy(fflux_hh(mptr)%ptr, fflux_hd(mptr)%ptr, nvar*lenbc*2+naux*lenbc)
!             ! call qad_cpu2(grid_data(mptr)%ptr,mitot,mjtot,nghost,nvar, &
!             !        fflux_hh(mptr)%ptr,fflux_hh(mptr)%ptr(locsvq),lenbc, &
!             !        intratx(level-1),intraty(level-1),hx,hy, &
!             !        delt,mptr,cc,zz)
!             ! cudaResult = cudaMemcpy(fflux_hd(mptr)%ptr, fflux_hh(mptr)%ptr, nvar*lenbc*2+naux*lenbc)
! 
!             call compute_kernel_size(numBlocks, numThreads, &
!                 1,2*(nx+ny))
!             call qad_gpu<<<numBlocks,numThreads,0,get_cuda_stream(id,device_id)>>>( &
!                    grid_data_d(mptr)%ptr,mitot,mjtot,nghost,nvar, &
!                    fflux_hd(mptr)%ptr,fflux_hd(mptr)%ptr(locsvq),lenbc, &
!                    intratx(level-1),intraty(level-1),hx,hy, &
!                    delt,mptr,max1d,cc,zz)
!         endif


!         if (dimensional_split .eq. 0) then
! !           # Unsplit method
!             call stepgrid_soa( &
!                     grid_data_d(mptr)%ptr,fms_d(mptr)%ptr,fps_d(mptr)%ptr,gms_d(mptr)%ptr,gps_d(mptr)%ptr, &
!                     mitot,mjtot,nghost, &
!                     delt,hx,hy,nvar, &
!                     xlow,ylow,time,mptr,naux,alloc(locaux),& 
!                     numgrids(level),id,cfls_d)
!         else if (dimensional_split .eq. 1) then
! !           # Godunov splitting
!             call stepgrid_dimsplit_soa( &
!                     grid_data_d(mptr)%ptr,fms_d(mptr)%ptr,fps_d(mptr)%ptr,gms_d(mptr)%ptr,gps_d(mptr)%ptr, &
!                     mitot,mjtot,nghost, &
!                     delt,hx,hy,nvar, &
!                     xlow,ylow,time,mptr,naux,alloc(locaux),& 
!                     numgrids(level),id,cfls_d)
!         else 
! !           # should never get here due to check in amr2
!             write(6,*) '*** Strang splitting not supported'
!         endif
!         cudaResult = cudaMemcpyAsync(grid_data(mptr)%ptr, grid_data_d(mptr)%ptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))


! test the new cudaclaw function here
        call stepgrid_cudaclaw(mitot,mjtot,nghost, &
            xlow, xlow+hx*mitot, ylow, ylow+hy*mjtot, delt, &
            grid_data_d_copy2(mptr)%ptr, grid_data_d(mptr)%ptr, &
            aux_d(mptr)%ptr, &
            nvar, naux, &
            cfls_d, numgrids(level), &
            id, device_id) 

        cudaResult = cudaMemcpyAsync(alloc(locnew), grid_data_d(mptr)%ptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))


        !$OMP END CRITICAL(launch)

    enddo
    !$OMP END DO
    
#ifdef PROFILE
    call endCudaProfiler() ! Launch qad and stepgrid_soa
    call cpu_timer_stop(timer_launch_compute_kernels)
#endif

!! ##################################################################
!! Copying old solution
!! ##################################################################
#ifdef PROFILE
    call take_cpu_timer('Copy q to old storage and update gauges', timer_copy_old_solution)
    call cpu_timer_start(timer_copy_old_solution)
    call startCudaProfiler("copy q to old storage and update gauges", 33)
#endif

    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost

        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy

        !$OMP CRITICAL(rv)
        rvol = rvol + nx * ny
        rvoll(level) = rvoll(level) + nx * ny
        !$OMP END CRITICAL(rv)



        ! grids(j)%fm => fms_d(mptr)%ptr
        ! grids(j)%fp => fps_d(mptr)%ptr
        ! grids(j)%gm => gms_d(mptr)%ptr
        ! grids(j)%gp => gps_d(mptr)%ptr
        ! grids(j)%mptr = mptr
        ! grids(j)%nx   = nx 
        ! grids(j)%ny   = ny 
        ! if (level > 1) then
        !     !$OMP CRITICAL (max_lenbc)
        !     max_lenbc = max(max_lenbc, 2*(nx/intratx(level-1)+ny/intraty(level-1)))
        !     !$OMP END CRITICAL (max_lenbc)
        ! endif

    enddo
    !$OMP END DO

#ifdef PROFILE
    call endCudaProfiler()
    call cpu_timer_stop(timer_copy_old_solution)
#endif

    !$OMP MASTER 
    call wait_for_all_gpu_tasks(device_id)
    ! grids_d = grids
    !$OMP END MASTER 
    !$OMP BARRIER
#ifdef PROFILE
    call endCudaProfiler()
    call cpu_timer_stop(timer_gpu_loop)
#endif

    
!! ##################################################################
!! fluxad and fluxsv
!! ##################################################################
#ifdef PROFILE
    call take_cpu_timer('fluxsv, fluxad and soa-aos conversion', timer_fluxsv_fluxad)
    call cpu_timer_start(timer_fluxsv_fluxad)
    call startCudaProfiler("fluxsv, fluxad and soa-aos conversion",12)
#endif

    ! one kernel launch to do fluxad for all grids at this level
    ! we don't do this for the coarsest level
    !$OMP MASTER 
#ifdef PROFILE
    call startCudaProfiler("Launch fluxad and fluxsv",12)
#endif
!     if (level > 1) then
!         call compute_kernel_size(numBlocks, numThreads, &
!             1,max_lenbc,1,numgrids(level))
!         call fluxad_fused_gpu<<<numBlocks,numThreads,0,get_cuda_stream(1,device_id)>>>( &
!                 grids_d, fflux_dd,&
!                 nghost, numgrids(level), intratx(level-1), intraty(level-1), &
!                 delt, hx, hy)
!     endif
!     ! one kernel launch to do fluxsv for all grids at this level
!     ! we don't do this for the fineset level
!     if (level < lfine) then
!         call compute_kernel_size(numBlocks, numThreads, &
!             1,listsp(level),1,numgrids(level))
!         call fluxsv_fused_gpu<<<numBlocks,numThreads,0,get_cuda_stream(2,device_id)>>>( &
!                  grids_d, cflux_dd, fflux_dd, &
!                  nghost, numgrids(level), nvar,listsp(level),delt,hx,hy)
!     endif
#ifdef PROFILE
    call endCudaProfiler() ! Launch fluxad and fluxsv
#endif
    !$OMP END MASTER 
    !$OMP BARRIER

#ifdef PROFILE
    call startCudaProfiler('soa_to_aos',99)
    call take_cpu_timer('soa_to_aos', timer_soa_to_aos)
    call cpu_timer_start(timer_soa_to_aos)
#endif

    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1, mptr)


        rnode(timemult,mptr)  = rnode(timemult,mptr)+delt

    enddo
    !$OMP END DO
#ifdef PROFILE
    call endCudaProfiler() ! soa_to_aos
    call cpu_timer_stop(timer_soa_to_aos)
#endif

    !$OMP MASTER 
    call wait_for_all_gpu_tasks(device_id)
    !$OMP END MASTER 
    !$OMP BARRIER

#ifdef PROFILE
    call take_cpu_timer('memory operation', timer_memory)
    call cpu_timer_start(timer_memory)
#endif
    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)


        ! call gpu_deallocate(grid_data_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(grid_data_d_copy2(mptr)%ptr,device_id)
        ! call gpu_deallocate(aux_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(fms_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(fps_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(gms_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(gps_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(sx_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(sy_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(wave_x_d(mptr)%ptr,device_id)
        ! call gpu_deallocate(wave_y_d(mptr)%ptr,device_id)

    enddo
    !$OMP END DO
#ifdef PROFILE
    call cpu_timer_stop(timer_memory)
#endif

#ifdef PROFILE
    call startCudaProfiler('Transfer fflux to CPU',99)
#endif

    ! !$OMP SINGLE
    ! if (level > 1) then
    !     do j = 1, numgrids(level)
    !         levSt = listStart(level)
    !         mptr = listOfGrids(levSt+j-1)
    !         nx   = node(ndihi,mptr) - node(ndilo,mptr) + 1
    !         ny   = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
    !         lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
    !         cudaResult = cudaMemcpyAsync(fflux_hh(mptr)%ptr, fflux_hd(mptr)%ptr, &
    !             nvar*lenbc*2+naux*lenbc,get_cuda_stream(id_copy_fflux,device_id))
    !     enddo
    ! endif
    ! !$OMP END SINGLE

    ! !$OMP SINGLE
    ! deallocate(grids)
    ! deallocate(grids_d)
    ! !$OMP END SINGLE
    !$OMP BARRIER

#ifdef PROFILE
    call endCudaProfiler() ! Transfer fflux to CPU
    call endCudaProfiler() ! fluxsv, fluxad and soa-aos conversion
    call cpu_timer_stop(timer_fluxsv_fluxad)
#endif

#else
    !$OMP DO SCHEDULE (DYNAMIC,1)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        !
        call par_advanc(mptr,mitot,mjtot,nvar,naux,dtnew)
        !$OMP CRITICAL (newdt)
        dtlevnew = min(dtlevnew,dtnew)
        !$OMP END CRITICAL (newdt)    
    end do
    !$OMP END DO

#endif




#ifdef CUDA


!! ##################################################################
!! CFL reduction
!! ##################################################################
#ifdef PROFILE
    call startCudaProfiler('CFL reduction', 34)
    call take_cpu_timer('CFL reduction', timer_cfl)
    call cpu_timer_start(timer_cfl)
#endif

    ! reduction to get cflmax and dtlevnew
    !$OMP MASTER
    cudaResult = cudaMemcpy(cfls, cfls_d, numgrids(level)*2)
    do j = 1,numgrids(level)
        cfl_local = max(cfls(1,j),cfls(2,j))
        if (cfl_local .gt. cflv1) then
            write(*,810) cfl_local
            write(outunit,810) cfl_local, cflv1
      810   format('*** WARNING *** Courant number  =', d12.4, &
          '  is larger than input cfl_max = ', d12.4)
        endif
        cfl_level = max(cfl_level, cfl_local)
    enddo
    dtlevnew = delt*cfl/cfl_level
    cflmax = max(cflmax, cfl_level)
    call cpu_deallocate_pinned(cfls)
    call gpu_deallocate(cfls_d,device_id)
    !$OMP END MASTER

#ifdef PROFILE
    call cpu_timer_stop(timer_cfl)
    call endCudaProfiler() ! CFL reduction
#endif

#else
    cflmax = max(cflmax, cfl_level)
#endif

#ifdef PROFILE
    call cpu_timer_stop(timer_stepgrid)
    call endCudaProfiler() ! Stepgrid
    call endCudaProfiler() ! advanc level 
#endif
    !$OMP MASTER
    timeStepgrid = timeStepgrid +clock_finish-clock_startStepgrid
    timeStepgridCPU=timeStepgridCPU+cpu_finish-cpu_startStepgrid      
    !$OMP END MASTER

!$OMP END PARALLEL 

#ifdef PROFILE
    call cpu_timer_stop(timer_advanc)
#endif

    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    tvoll(level) = tvoll(level) + clock_finish - clock_start
    tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start

    return
end subroutine advanc
    !
    ! -------------------------------------------------------------
    !
subroutine prepgrids(listgrids, num, level)

    use amr_module
    implicit real(CLAW_REAL) (a-h,o-z)
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
end subroutine prepgrids

