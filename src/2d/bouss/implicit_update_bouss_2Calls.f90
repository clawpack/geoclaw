subroutine implicit_update(nvar,naux,levelBouss,numBoussCells,doUpdate,time)

    ! For source terms that require an implicit solve, in the case when there
    ! is a single levelBouss at which the source terms should be applied.
    
    ! Set up and solve the linear system using sparse matrix and gmres.
    ! Then use the resulting soln to update q by looping over all 
    ! Bouss level patches.
    
    ! Boussinesq version.
    
    use amr_module
    use topo_module, only: aux_finalized
    use bouss_module
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    logical, intent(in) :: doUpdate
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8) :: soln(0:2*numBoussCells)
    real(kind=8) :: rhs(0:2*numBoussCells)
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: i,j,nb,nn,levSt,nx,ny,mitot,mjtot,loc,locaux,mptr
    integer :: locnew, locold, locOther
    real(kind=8) :: dt
    
    integer :: nD, nst, ncc, k
    integer(kind=8) :: clock_startBound,clock_finishBound,clock_rate
    integer(kind=8) :: clock_startLinSolve,clock_finishLinSolve
    real(kind=8) cpu_startBound,cpu_finishBound, time 
    real(kind=8) cpu_startLinSolve,cpu_finishLinSolve
    logical :: debug

#ifdef WHERE_AM_I
    write(*,*) "starting implicit_update for level ",levelBouss
#endif
    dt = possk(levelBouss)
    nD = numBoussCells ! shorthand for size of matrix blocks D1 to D4
    debug = .false.
    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.
    
    minfo => matrix_info_allLevs(levelBouss)

    
! =======================  step 1: fill ghost cells ==============================

! need ghost cells after finite vol update, before next phase of differencingc 
! get start time for more detailed timing by level
       call system_clock(clock_startBound,clock_rate)
       call cpu_time(cpu_startBound)

    !!write(*,*) "before bound loop"
    !!call timestamp()
    !! reset isBouss flag array if topo/aux still changing

    if (aux_finalized .lt. 2) then
      call setBoussFlag(levelBouss,naux)
    endif
!$OMP PARALLEL DO PRIVATE(j,locnew,locaux,mptr,nx,ny,mitot,                   &
!$OMP                     mjtot,levSt),                                       &
!$OMP             SHARED(levelBouss,nvar,naux,alloc,nghost,time,              &
!$OMP                    node,rnode,numgrids,listStart,listOfGrids),          &
!$OMP             SCHEDULE (dynamic,1)                                        &
!$OMP             DEFAULT(none)
      
      do  nn = 1, numgrids(levelBouss)
          levSt  = listStart(levelBouss)
          mptr   = listOfGrids(levSt+nn-1)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
 
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,     &
                     alloc(locaux),naux)
        end do
!$OMP END PARALLEL DO

      call system_clock(clock_finishBound,clock_rate)
      call cpu_time(cpu_finishBound)
      timeBound = timeBound + clock_finishBound - clock_startBound
      timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound

    ! ================  Step 2: build matrix and rhs  ======================

    if (ibouss .eq. 1) then
      ! Madsen-Sorensen, always uses COO format:
      call buildSparseMatrixMScoo(rhs,nvar,naux,levelBouss,numBoussCells)
    else ! ibouss 2 means SGN
      if (.not. crs) then ! COO triplet format
        call prepBuildSparseMatrixSGNcoo(soln,rhs,nvar,naux,levelBouss,numBoussCells)
      else ! CRS format
         minfo%vals = 0
         minfo%cols = -1  ! to recognize if not overwritten with real col indices
        call prepBuildSparseMatrixSGNcrs(soln,rhs,nvar,naux,levelBouss,numBoussCells)
      endif
    endif
    
    !================   Step 4 Solve matrix system =======================

    if (isolver .eq.1) then  ! use gmres
       write(*,*)" No longer supporting gmres option"
       write(outunit,*)" No longer supporting gmres option"
       stop
            
    else if (isolver .eq.2) then  ! use pardiso
    
#ifdef HAVE_PARDISO
       write(*,*)" No longer supporting pardiso option"
       write(outunit,*)" No longer supporting pardiso option"
       stop
       !if (minfo%matrix_nelt .gt. 0) then
       !  call system_clock(clock_startLinSolve,clock_rate)
       !  call cpu_time(cpu_startLinSolve)
       !  ! convert to compressed row format
       !  nst = minfo%matrix_nelt
       !  call st_to_cc_size(nst,minfo%matrix_ia,minfo%matrix_ja,ncc)
       !  call pardiso_driver(soln,rhs,levelBouss,numBoussCells,nst,ncc,aux_finalized)
       !  call system_clock(clock_finishLinSolve,clock_rate)
       !  call cpu_time(cpu_finishLinSolve)
       !  timeLinSolve = timeLinSolve + clock_finishLinSolve - clock_startLinSolve
       !else
       !  go to 99
       !endif
#endif

    else if (isolver .eq.3) then  ! use petsc
#ifdef HAVE_PETSC
       ! now add Petsc. Not all levels are bouss grids so test first
       ! for now test for either CRS or COO option for now
       if (minfo%numColsTot .gt. 0 .or. minfo%matrix_nelt .gt. 0) then
          call system_clock(clock_startLinSolve,clock_rate)
          call cpu_time(cpu_startLinSolve)
          call petsc_driver(soln,rhs,levelBouss,numBoussCells,time,     &
                           aux_finalized,minfo%matrix_nelt)
          call system_clock(clock_finishLinSolve,clock_rate)
          call cpu_time(cpu_finishLinSolve)
          !write(89,*)" level ",levelBouss,"  rhs   soln"
          !do k = 1, 2*numBoussCells
          !   write(89,899) k,rhs(k), soln(k)
  899         format(i5,2e15.7)
          !end do
          timeLinSolve = timeLinSolve + clock_finishLinSolve - clock_startLinSolve
          timeLinSolveCPU = timeLinSolveCPU+cpu_finishLinSolve - cpu_startLinSolve
       else
         go to 99 ! IS THERE ANYTHING TO UPDATE
       endif
#endif
    
    endif

    !  DEBUGGING: have soln, and matrix. Multiply and check diff with rhs
    if (.not. crs .and. debug) then ! matvec not yet implemented for crs
       call testSoln(2*numBoussCells,soln,rhs,minfo%matrix_nelt,minfo%matrix_ia,   &
                  minfo%matrix_ja,minfo%matrix_sa)
    endif

    !================   Step 5  Update solution  =======================
            
    ! Use solution as updated q by looping over grids at levelBouss
    ! This can be done in parallel on patches
    
    ! using forward Euler at this stage to update momenta
    ! we are solving P_t = -S  where S is the solution of the linear system
    ! so update P to P - dt*S

!$OMP  PARALLEL DO PRIVATE(mptr,locnew,locOther,locaux,nx,ny,mitot,mjtot,nb,i,j,levSt),     &
!$OMP  SHARED(listOfGrids,listStart,alloc,soln,levelBouss,ibouss,crs,         &
!$OMP         doUpdate,nvar,naux,numgrids,node,dt,possk,nghost,nD,rhs,mxnest),    &  
!$OMP              SCHEDULE (dynamic,1),                                          &  
!$OMP              DEFAULT(none)
    do  nn = 1, numgrids(levelBouss)
        levSt  = listStart(levelBouss)
        mptr   = listOfGrids(levSt+nn-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1,mptr)  
        locaux = node(storeaux,mptr)  
        if (levelBouss .lt. mxnest) then
           ! old storage location exists
           locOther = node(store2,mptr)  
        else
           ! want to save in same place, beware
           locOther = locnew
        endif
    
        nb = node(bpatchNum, mptr)
        !if (nb > 0) then  ! is a bouss patch
           if (ibouss .eq. 1) then
              call solnUpdate(alloc(locnew),alloc(locOther),nb,mitot,mjtot,nghost,nvar, &
                           dt,soln,nD,mptr,levelBouss,rhs,doUpdate)
           else ! ibouss 2 SGN, saving different stuff, different update
                call solnUpdateSGN(alloc(locnew),alloc(locOther),nb,mitot,mjtot,nghost,nvar, &
                     dt,soln,nD,mptr,levelBouss,rhs,doUpdate,alloc(locaux),naux)
           endif
        !endif
    end do        
!$OMP END PARALLEL DO

    !! next loop to copy from either other grids at same level and to
    !! update ghost cell momenta with newly copied psi
!$OMP PARALLEL DO PRIVATE(mptr,locnew,locaux,locold,nx,ny,mitot,mjtot,nb, &
!$OMP                      i,j,levSt),   &
!$OMP              SHARED(possk,listOfGrids,listStart,alloc,soln,levelBouss, &
!$OMP                     nvar,naux,numgrids,node,dt,nghost,nD,rhs,time),   &  
!$OMP              SCHEDULE (dynamic,1),                                               &  
!$OMP              DEFAULT(none)
    do  nn = 1, numgrids(levelBouss)
        levSt  = listStart(levelBouss)
        mptr   = listOfGrids(levSt+nn-1)
        dt     = possk(levelBouss)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1,mptr)
        locaux = node(storeaux,mptr)
        nb = node(bpatchNum, mptr)
        !if (nb > 0) then  ! is a bouss patch
             ! might want to add this to previous boussUpdate routine and     
             ! and check for ghost cell overlap on fine patch at same
             ! time as update fine -> coarse
             ! still need ghost cell filling. simpleBound only does
             ! nboring grid filling at same level
             call simpleBound(time,nvar,alloc(locnew),mitot,mjtot,mptr,alloc(locaux),naux)
        !endif
    end do      
!$OMP END PARALLEL DO

 99     continue

#ifdef WHERE_AM_I
    write(*,*) "ending   implicit_update for level ",levelBouss
#endif
        

contains

    ! Index into q array
    integer pure function iadd(ivar,i,j)
        integer, intent(in) :: i, j, ivar
        !! ghost cells added to iadd definition, NOT to call
        iadd = loc + ivar-1 + nvar*((j-1+nghost)*mitot+i+nghost-1)
    end function iadd

    ! Index into aux array
    integer pure function iaddaux(iaux,i,j)
        integer, intent(in) :: i, j, iaux  
        !! ghost cells added to iaddaux definition, NOT to call
        iaddaux = locaux + iaux-1 + naux*((j-1+nghost)*mitot+i+nghost-1)
    end function iaddaux



end subroutine implicit_update

subroutine solnUpdate(valnew,valOther,nb,mitot,mjtot,nghost,nvar,dt, &
                      soln,nD,mptr,levelBouss,rhs,doUpdate)

    ! if doing the momentum update, put it in valnew.
    ! but save  psi in valold, since goes with previous
    ! momenta and h (in new bc approach)

    ! if not doing the update, then is the new solve added for
    ! a level AFTER the SWE update, for use in interpolating for
    ! fine psi.  So it goes in the valnew storage.
    
    ! Currently only doing first-order forward Euler update.
    
    use amr_module, only: outunit,mxnest
    use bouss_module, only: matrix_info_allLevs
    use bouss_module, only: matrix_levInfo, matrix_patchIndex
    implicit none

    integer, intent(in) :: nvar,nghost,mitot,mjtot,nD,mptr,nb,levelBouss
    logical, intent(in) :: doUpdate
    real(kind=8), intent(in) :: soln(0:2*nD), dt, rhs(0:2*nD)

    real(kind=8), intent(inout) :: valnew(nvar,mitot,mjtot)
    real(kind=8), intent(inout) :: valOther(nvar,mitot,mjtot)
    ! local variable
    integer :: i,j,k_ij
    integer :: imax1,jmax1,imax2,jmax2
    logical :: debug
    real(kind=8) :: solmax1,solmax2
    real(kind=8) :: valset

    !type(matrix_index), pointer :: mi
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo

    minfo =>matrix_info_allLevs(levelBouss)
    mi => minfo%matrix_indices(nb)

    debug = .false.
    solmax1 = 0.d0
    solmax2 = 0.d0
    imax1 = -1
    jmax1 = -1
    imax2 = -1
    jmax2 = -1

    !! if not updating momentum, just save psi, translating from
    !! soln vector back to patch variables
    !!if (.not. doUpdate .and. nvar .eq. 3) then
    !!    write(*,*) " this option doesn't make sense"
    !!    stop
    !!endif

    !! if finest level grid, dont have 2 levels of storage
    !! no need to save in old since wont need for interpolation

    if (doUpdate) then
       do j=nghost+1, mjtot-nghost
       do i=nghost+1, mitot-nghost
            k_ij = mi%mindex(i-nghost,j-nghost)  ! hu element index for (i,j)
            valnew(2,i,j) = valnew(2,i,j) + dt * soln(k_ij)
            valnew(3,i,j) = valnew(3,i,j) + dt * soln(k_ij + nD)
            if (.false.) write(67,*) '+++ i,j,soln: ', &
                        i-nghost,j-nghost,soln(k_ij),soln(k_ij + nD)
            if (abs(soln(k_ij)) .gt. solmax1) then
               solmax1 = abs(soln(k_ij))
               imax1 = i
               jmax1 = j
            endif
            if (abs(soln(k_ij+nD)) .gt. solmax2) then
               solmax2 = abs(soln(k_ij+nD))
               imax2 = i
               jmax2 = j
            endif

            if (nvar == 5) then
                ! save update in components 4,5 for Dirichlet BC next iteration
                ! for coarser levels, valOther is valOld
                ! on the finest grid we only have one storage loc
                ! and valOther is valNew, since dont need to save old
                ! and new for interpolation
                    valOther(4,i,j) = soln(k_ij)
                    valOther(5,i,j) = soln(k_ij+nD)
            endif
       enddo
       enddo

       ! next update momenta in ghost cells with interpolated psi. Will
       ! be overwritten with neighboring grids at same level in
       ! simpleBound so just do all ghost cells now
       valnew(2:3,1:mitot,1:nghost) = valnew(2:3,1:mitot,1:nghost) +  &
                                  dt* valnew(4:5,1:mitot,1:nghost)
       valnew(2:3,1:mitot,mjtot-nghost+1:mjtot) = valnew(2:3,1:mitot,mjtot-nghost+1:mjtot) + &
                                              dt* valnew(4:5,1:mitot,mjtot-nghost+1:mjtot)
       valnew(2:3,1:nghost,nghost+1:mjtot-nghost) = valnew(2:3,1:nghost,nghost+1:mjtot-nghost) +  &
                                                 dt* valnew(4:5,1:nghost,nghost+1:mjtot-nghost)
       valnew(2:3,mitot-nghost+1:mitot,nghost+1:mjtot-nghost) =                &
                      valnew(2:3,mitot-nghost+1:mitot,nghost+1:mjtot-nghost) +  &
                  dt* valnew(4:5,mitot-nghost+1:mitot,nghost+1:mjtot-nghost)

    else !! just save psi
       do j=nghost+1, mjtot-nghost
       do i=nghost+1, mitot-nghost
            k_ij = mi%mindex(i-nghost,j-nghost)  ! hu element index for (i,j)
            valnew(4,i,j) = soln(k_ij)
            valnew(5,i,j) = soln(k_ij+nD)
            if (abs(soln(k_ij)) .gt. solmax1) then
               solmax1 = abs(soln(k_ij))
               imax1 = i
               jmax1 = j
            endif
            if (abs(soln(k_ij+nD)) .gt. solmax2) then
               solmax2 = abs(soln(k_ij+nD))
               imax2 = i
               jmax2 = j
            endif
        end do
        end do
    endif

    if (debug) then
       write(*,*)"When doUpate = ",doUpdate,":"
       write(*,100) solmax1,imax1,jmax1
 100   format("max soln update for hu ",e15.7, " at ",i5,i5)
       write(*,101) solmax2,imax2,jmax2
 101   format("max soln update for hv ",e15.7, " at ",i5,i5)
    endif


end subroutine solnUpdate

! =========================================================
subroutine solnUpdateSGN(valnew,valOther,nb,mitot,mjtot,nghost,nvar,dt, &
                      soln,nD,mptr,levelBouss,rhs,doUpdate,aux,naux)

    ! if doing the momentum update, put it in valnew.
    ! but save  D from SGN solve  in valold, since goes with previous
    ! momenta and h (in new bc approach)

    ! if not doing the update, then is the new solve added for
    ! a level AFTER the SWE update, for use in interpolating for
    ! fine psi.  So it goes in the valnew storage.
    
    ! Currently only doing first-order forward Euler update.
    
    use amr_module, only: outunit,mxnest
    use bouss_module, only: matrix_info_allLevs,boussMindepth,crs
    use bouss_module, only: matrix_levInfo, matrix_patchIndex
    use bouss_module, only: alpha, origCooFormat
    use amr_module, only: rnode,cornxlo,cornylo,hxposs,hyposs
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: dry_tolerance
    implicit none

    integer, intent(in) :: nvar,naux,nghost,mitot,mjtot,nD,mptr,nb,levelBouss
    logical, intent(in) :: doUpdate
    real(kind=8), intent(in) :: soln(0:2*nD), dt, rhs(0:2*nD)

    real(kind=8), intent(inout) :: valnew(nvar,mitot,mjtot)
    real(kind=8), intent(inout) :: valOther(nvar,mitot,mjtot)
    real(kind=8), intent(in) :: aux(naux,mitot,mjtot)
    ! local variable
    integer :: ng,i,j,k_ij,ivar,iunit,k,imax,jmax
    integer :: ist,iend,jst,jend
    logical :: debug
    real(kind=8) :: eta(mitot,mjtot)
    real(kind=8) :: etax(mitot,mjtot), etay(mitot,mjtot)
    real(kind=8) :: speed(mitot,mjtot), speedmax 
    real(kind=8) :: dx,dy,dxm,dym
    real(kind=8) :: xlowWithGhost,ylowWithGhost,x,y

    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo

    logical IS_INNER_GHOST

    IS_INNER_GHOST(i,j) = (i .le. nghost .or. i .gt. mitot-nghost .or.   &
                     j .le. nghost .or. j .gt. mjtot-nghost)

    if (nvar .ne. 5) then
      write(*,*)" this SGN routine assumes 5 equations"
      stop
    endif

    minfo =>matrix_info_allLevs(levelBouss)
    mi => minfo%matrix_indices(nb)

    dy = hyposs(levelBouss)
    dx = hxposs(levelBouss)

    if (coordinate_system .eq. 2) then
      dym = dy *earth_radius*deg2rad
    else
      dym = dy
    endif

    !debug = .true.
    debug = .false.
    !if (levelBouss .eq. 1) debug = .true.

    !! if not updating momentum, just save psi, translating from
    !! soln vector back to patch variables
    if (debug) then
      write(*,*)"start of solnUpdateSGN, level ",levelBouss," doUpate = ",doUpdate
    endif

    !! if finest level grid, dont have 2 levels of storage
    !! no need to save in old since wont need for interpolation

    ! compute eta used below in update
    eta(:,:) = valnew(1,:,:) + aux(1,:,:)
       
    ! grad eta needed for updates. do one-sided at cell edges
    do j = 1, mjtot
         if (coordinate_system .eq. 2) then
            dxm = dx*earth_radius*deg2rad * cos(aux(3,1,j))
         else
            dxm = dx
         endif
       etax(1,j) = (eta(2,j) - eta(1,j))/dxm
       etax(mitot,j) = (eta(mitot,j) - eta(mitot-1,j))/dxm
    do i = 2, mitot-1
       etax(i,j) = (eta(i+1,j) - eta(i-1,j))/(2.d0*dxm)
    end do
    end do

    do i = 1, mitot
       etay(i,1) = (eta(i,2) - eta(i,1))/dym
       etay(i,mjtot) = (eta(i,mjtot) - eta(i,mjtot-1))/dym
    do j = 2, mjtot-1
       etay(i,j) = (eta(i,j+1) - eta(i,j-1))/(2.d0*dym)
    end do
    end do

    if (doUpdate) then
       do j=nghost+1, mjtot-nghost
       do i=nghost+1, mitot-nghost
            k_ij = mi%mindex(i-nghost,j-nghost)  ! hu element index for (i,j)
            if (mi%isBouss(i-nghost,j-nghost)) then
              ! check if soln == 0 as proxy if cell is more shallow than deep bouss
              ! and no update
                if (.not. crs) then
                  ! COO triplet format
                  if (origCooFormat) then
                   valnew(2,i,j) = valnew(2,i,j) + dt * valnew(1,i,j)*     &
                          (grav/alpha * etax(i,j) - soln(k_ij))
                   valnew(3,i,j) = valnew(3,i,j) + dt * valnew(1,i,j)*     &
                          (grav/alpha * etay(i,j) - soln(k_ij+nD))
                   if (debug) then
                      write(*,901) i,j,k_ij,soln(k_ij),soln(k_ij+nD),rhs(k_ij),rhs(k_ij+nD)
                   endif
                  else ! block format, 1-based
  901              format(3i9,4e18.8)
                   valnew(2,i,j) = valnew(2,i,j) + dt * valnew(1,i,j)*     &
                          (grav/alpha * etax(i,j) - soln(2*k_ij-1))
                   valnew(3,i,j) = valnew(3,i,j) + dt * valnew(1,i,j)*     &
                          (grav/alpha * etay(i,j) - soln(2*k_ij))
                   if (debug) then
                      write(*,901) i,j,k_ij,soln(2*k_ij-1),soln(2*k_ij),rhs(2*k_ij-1),rhs(2*k_ij)
                   endif
                 endif
                else  ! CRS format,  different mapping, 0 based
                  valnew(2,i,j) = valnew(2,i,j) + dt * valnew(1,i,j)*     &
                         (grav/alpha * etax(i,j) - soln(2*(k_ij-1)))
                  valnew(3,i,j) = valnew(3,i,j) + dt * valnew(1,i,j)*     &
                         (grav/alpha * etay(i,j) - soln(2*k_ij-1))
                  if (debug) then
                     write(*,901)i,j,k_ij,soln(2*(k_ij-1)),soln(2*k_ij-1),rhs(2*(k_ij-1)),rhs(2*k_ij-1)
                  endif
                endif

                ! save update in components 4,5 for Dirichlet BC next iteration
                ! for coarser levels, valOther is valOld
                ! on the finest grid we only have one storage loc
                ! and valOther is valNew, since dont need to save old
                ! and new for interpolation
                if (.not. crs) then
                    if (origCooFormat) then
                       valOther(4,i,j) = soln(k_ij)
                       valOther(5,i,j) = soln(k_ij+nD)
                    else ! block format
                       valOther(4,i,j) = soln(2*k_ij-1)
                       valOther(5,i,j) = soln(2*k_ij)
                    endif
                else ! crs format
                    !valOther(4,i,j) = soln(2*k_ij-1)
                    !valOther(5,i,j) = soln(2*k_ij)
                    valOther(4,i,j) = soln(2*(k_ij-1))
                    valOther(5,i,j) = soln(2*k_ij-1)
                 endif
            else  ! not a bouss cell. set to zero
               valOther(4,i,j) = 0.d0
               valOther(5,i,j) = 0.d0
            endif
       enddo
       enddo

       ! next update momenta in ghost cells with interpolated psi. Will
       ! be overwritten with neighboring grids at same level in
       ! simpleBound so just do all ghost cells now
       ! bottom bndry. Use first interior cell to determine if
       ! ghost cells should be updated. CORNERS?

       ! remember mindex only from 1:nx, 1:ny so subtract nghost as needed
       ! below code left in to remember why it doesnt work
       go to 444
       ! left side
       do j = nghost+1, mjtot-nghost
          !k_ij = mi%mindex(1,j-nghost)  ! 1st interior cell for row i
          !if (doBouss(k_ij)) then
          !if (mi%isBouss(1,j-nghost)) then
          do k = 1, nghost
          !if (aux(1,k,j) .lt. 0.d0) then  ! water, maybe check depth for bouss
          if (-aux(1,k,j) .gt. boussMinDepth) then  ! water, maybe check depth for bouss
             valnew(2,k,j) = valnew(2,k,j) + dt*valnew(1,k,j)   &
                          * (grav/alpha*etax(k,j) - valnew(4,k,j))
             valnew(3,k,j) = valnew(3,k,j) + dt*valnew(1,k,j)   &
                          * (grav/alpha*etay(k,j) - valnew(5,k,j))
          endif
          end do
        end do
       ! right side
       do j = nghost+1, mjtot-nghost
          !k_ij = mi%mindex(mitot-2*nghost,j-nghost)  ! last interior cell row i
          !if (doBouss(k_ij)) then
          !if (mi%isBouss(mitot-2*nghost,j-nghost)) then
          do k = mitot-nghost+1, mitot
          !if (aux(1,k,j) .lt. 0.d0) then
          if (-aux(1,k,j) .gt. boussMinDepth) then
             valnew(2,k,j) = valnew(2,k,j) + dt*valnew(1,k,j)         &
                          * (grav/alpha*etax(k,j) - valnew(4,k,j))
             valnew(3,k,j) = valnew(3,k,j) + dt*valnew(1,k,j)         &
                          * (grav/alpha*etay(k,j) - valnew(5,k,j))
          endif
          end do
        end do
       ! bottom side
       do i = nghost+1, mitot-nghost
          !k_ij = mi%mindex(i-nghost,1)  ! 1st interior cell for row i
          !if (doBouss(k_ij)) then
          !if (mi%isBouss(i-nghost,1)) then
          do k = 1, nghost
          !if (aux(1,i,k) .lt. 0.d0) then
          if (-aux(1,i,k) .gt. boussMinDepth) then
             valnew(2,i,k) = valnew(2,i,k) + dt*valnew(1,i,k)   &
                          * (grav/alpha*etax(i,k) - valnew(4,i,k))
             valnew(3,i,k) = valnew(3,i,k) + dt*valnew(1,i,k)   &
                          * (grav/alpha*etay(i,k) - valnew(5,i,k))
          endif
          end do
        end do
       ! top side
       do i = nghost+1, mitot-nghost
          !k_ij = mi%mindex(i-nghost,mjtot-2*nghost)  ! last interior cell row i
          !if (doBouss(k_ij)) then
          !if (mi%isBouss(i-nghost,mjtot-2*nghost)) then
          do k = mjtot-nghost+1,mjtot
          !if (aux(1,i,k) .lt. 0.d0) then
          if (-aux(1,i,k) .gt. boussMinDepth) then
             valnew(2,i,k) = valnew(2,i,k) + dt*valnew(1,i,k)         &
                          * (grav/alpha*etax(i,k) - valnew(4,i,k))
             valnew(3,i,k) = valnew(3,i,k) + dt*valnew(1,i,k)         &
                          * (grav/alpha*etay(i,k) - valnew(5,i,k))
          endif
          end do
        end do
   444 continue

    else !! just save psi. TEST FOR doBOUSS - SET TO 0?
       do j=nghost+1, mjtot-nghost
       do i=nghost+1, mitot-nghost
            k_ij = mi%mindex(i-nghost,j-nghost)  ! hu element index for (i,j)
            if (mi%isBouss(i-nghost,j-nghost)) then
              if (.not. crs) then
                if (origCooFormat) then
                  valnew(4,i,j) = soln(k_ij)
                  valnew(5,i,j) = soln(k_ij+nD)
                  if (debug) then
                     write(*,901)i,j,k_ij,soln(k_ij),soln(k_ij+nD),rhs(k_ij),rhs(k_ij+nD)
                  endif
                else ! block format
                  valnew(4,i,j) = soln(2*k_ij-1)
                  valnew(5,i,j) = soln(2*k_ij)
                  if (debug) then
                     write(*,901)i,j,k_ij,soln(2*k_ij-1),soln(2*k_ij),rhs(2*k_ij-1),rhs(2*k_ij)
                  endif
                endif
              else ! using crs
                 valnew(4,i,j) = soln(2*(k_ij-1))
                 valnew(5,i,j) = soln(2*k_ij-1)
                 if (debug) then
                    write(*,901)i,j,k_ij,soln(2*(k_ij-1)),soln(2*k_ij-1),rhs(2*(k_ij-1)),rhs(2*k_ij-1)
                 endif
              endif
            else ! not bouss cell
              valnew(4,i,j) = 0.d0
              valnew(5,i,j) = 0.d0
            endif
       end do
       end do
    endif

    if (debug) then
       xlowWithGhost = rnode(cornxlo,mptr) - nghost*dx
       ylowWithGhost = rnode(cornylo,mptr) - nghost*dy
       write(34,*)" end of solnUpdateSGN for grid ",mptr," with doUpdate ",doUpdate
       ng = 2 ! indicates to print ghost cells too
       iunit = 34
       call lookAtGrid(valnew,valOther,dx,dy,xlowWithGhost,ylowWithGhost,nvar,   &
                       mitot,mjtot,ng,iunit,doUpdate,levelBouss)
    endif

end subroutine solnUpdateSGN

subroutine symcheck(val,mitot,mjtot)

    ! only test for vertical symmetry here.

    real(kind=8), intent(in) :: val(5,mitot,mjtot)
    integer mitot, mjtot, i,j
    double precision eps,hucDif,hvcDif

    eps = 1.d-10
    do i = 1, mitot
    do j = 1, mjtot/2
       hucDif =  abs(val(4,i,j) - val(4,i,mjtot+1-j))
       hvcDif =  abs(val(5,i,j) + val(5,i,mjtot+1-j))
       if (hucDif .gt. eps .or. hvcDif .gt. eps) then
          write(31,100) i,j,hucDif,hvcDif
 100      format("Cell ",2i5," difs 4/5 = ",2e12.7)
       endif
    end do
    end do

    return
end subroutine symcheck

! =========================================================

      subroutine symcheck2(val,mitot,mjtot)

      ! only test for vertical symmetry here.

      real(kind=8), intent(in) :: val(5,mitot,mjtot)
      integer mitot, mjtot, i,j
      double precision eps,hucDif,hvcDif

      eps = 1.d-10

      hDifmax = 0.d0
      huDifmax = 0.d0
      hvDifmax = 0.d0
      hucDifmax = 0.d0
      hvcDifmax = 0.d0

      write(31,*)"checking vertical symmetry"
      do i = 1, mitot
      do j = 1, mjtot/2
         hDif =  abs(val(1,i,j) - val(1,i,mjtot+1-j))
         huDif =  abs(val(2,i,j) - val(2,i,mjtot+1-j))
         hvDif =  abs(val(3,i,j) + val(3,i,mjtot+1-j))
         hucDif =  abs(val(4,i,j) - val(4,i,mjtot+1-j))
         hvcDif =  abs(val(5,i,j) + val(5,i,mjtot+1-j))
         if (hucDif .gt. eps .or. hvcDif .gt. eps .or. hDif .gt. eps   &
            .or. huDif .gt. eps .or. hvDif .gt. eps) then
            write(31,100) i,j,hDif,huDif,hvDif,hucDif,hvcDif
         endif
         hDifmax = max(hDifmax,hDif)
         huDifmax = max(huDifmax,huDif)
         hvDifmax = max(hvDifmax,hvDif)
         hucDifmax = max(hucDifmax,hucDif)
         hvcDifmax = max(hvcDifmax,hvcDif)
      end do
      end do

      write(31,101)hDifmax,huDifmax,hvDifmax,hucDifmax,hvcDifmax  
      write(*,101)hDifmax,huDifmax,hvDifmax,hucDifmax,hvcDifmax  
 101  format("max vertical diff found: ",5e13.6)
  
      hDifmax = 0.d0
      huDifmax = 0.d0
      hvDifmax = 0.d0
      hucDifmax = 0.d0
      hvcDifmax = 0.d0

      write(31,*)"checking horizntl symmetry"
      do j = 1, mjtot
      do i = 1, mitot/2
         hDif =  abs(val(1,i,j) - val(1,mitot+1-i,j))
         huDif =  abs(val(2,i,j) + val(2,mitot+1-i,j))
         hvDif =  abs(val(3,i,j) - val(3,mitot+1-i,j))
         hucDif =  abs(val(4,i,j) + val(4,mitot+1-i,j))
         hvcDif =  abs(val(5,i,j) - val(5,mitot+1-i,j))
         if (hucDif .gt. eps .or. hvcDif .gt. eps .or. hDif .gt. eps   &
            .or. huDif .gt. eps .or. hvDif .gt. eps) then
            write(31,100) i,j,hDif,huDif,hvDif,hucDif,hvcDif
 100        format("Cell ",2i5," difs = ",5e13.6)
         endif
         hDifmax = max(hDifmax,hDif)
         huDifmax = max(huDifmax,huDif)
         hvDifmax = max(hvDifmax,hvDif)
         hucDifmax = max(hucDifmax,hucDif)
         hvcDifmax = max(hvcDifmax,hvcDif)
      end do
      end do

      write(31,102)hDifmax,huDifmax,hvDifmax,hucDifmax,hvcDifmax  
      write(*,102)hDifmax,huDifmax,hvDifmax,hucDifmax,hvcDifmax  
 102  format("max horizntl diff found: ",5e13.6)
  
      return
      end subroutine symcheck2
