subroutine prepBuildSparseMatrixSGNcoo(soln,rhs,nvar,naux,levelBouss,numBoussCells)

    
    use amr_module
    use bouss_module
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: sea_level
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8) :: soln(0:2*numBoussCells), rhs(0:2*numBoussCells) 
    real(kind=8) fms
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: i,j,nng,levSt,nx,ny,mitot,mjtot,loc,locold,locaux,mptr

    integer(kind=8) :: clock_start, clock_finish, clock_rate
    real(kind=8) :: cpu_start, cpu_finish
    
    real(kind=8), allocatable, dimension(:) :: uv,Auv,rowsum
    real(kind=8), allocatable, dimension(:,:) :: hx, hy, phi
    real(kind=8), allocatable, dimension(:,:) :: Bx, By, Bxy, Bxx, Byy
    real(kind=8), allocatable, dimension(:,:) :: etax, etay
    
    integer :: nD, k
    real(kind=8) :: matrix_rhs1_max, matrix_rhs2_max
    integer :: minfo_matrix_ia(24*numBoussCells)
    integer :: minfo_matrix_ja(24*numBoussCells)
    real(kind=8) :: minfo_matrix_sa(24*numBoussCells)

    external :: MATVEC
    logical :: debug
    logical :: yhi_db, ylo_db, xlo_db, xhi_db

    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.
    
#ifdef WHERE_AM_I
  write(*,*) "starting prepBuildSparseMatrixSGNcoo"
#endif

    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)
    
    debug = .false.
    !debug = .true.
    
    minfo => matrix_info_allLevs(levelBouss)

!
! ================== Step 2  set up for matrix elements  ======================
    

    nD = numBoussCells   ! shorthand for size of matrix blocks D1 to D4
    
    !allocate(uv(2*nD), Auv(2*nD)) ! for testing matvec
    
    ! For Boussinesq, matrix is 2*nD by 2*nD,
    !     A = [I - D1   -D2  ]
    !         [ -D3    I - D4]
    ! and the solution vector consists of 
    !     updates to hu in all Bouss cells, followed by
    !     updates to hv in all Bouss cells
    
    ! Set up matrix A in sparse triad form: 
    !    matrix_sa(k) is the nonzero element if A in 
    !    row matrix_ia(k), and column matrix_ja(a)
    
    ! count number of matrix elements as added, 
    minfo%matrix_nelt  = 0

    rhs(:) = 0.d0
    
    matrix_rhs1_max = 0.d0
    matrix_rhs2_max = 0.d0
    
    ! loop over Bouss patches, must be done sequentially:
    
    levSt  = listStart(levelBouss)
    
    nng_loop: do  nng = 1, numgrids(levelBouss)
    
        mptr   = listOfGrids(levSt+nng-1)
        
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        loc = node(store1,mptr)  ! needed for iadd
        locold = node(store2,mptr)  ! needed for iadd
        if (locold .eq. 0) locold = loc  ! finest level has no old time soln
        locaux = node(storeaux,mptr)  ! needed for iaddaux

        if (origCooFormat) then
           call buildSparseMatrixSGNcoo(alloc(loc),alloc(locold),          &
                        alloc(locaux), soln,rhs,                           &
                        minfo%matrix_ia,minfo%matrix_ja,minfo%matrix_sa,   &
                        numBoussCells,levelBouss,mptr,nx,ny,nvar,naux,     &
                        matrix_rhs1_max,matrix_rhs2_max)
         else
           call buildSparseMatrixSGNcoo_blocks(alloc(loc),alloc(locold),   &
                        alloc(locaux), soln,rhs,                           &
                        minfo%matrix_ia,minfo%matrix_ja,minfo%matrix_sa,   &
                        numBoussCells,levelBouss,mptr,nx,ny,nvar,naux,     &
                        matrix_rhs1_max,matrix_rhs2_max)
         endif
            
    
    enddo  nng_loop


    if (minfo%matrix_nelt .gt. 0 .and. debug) then
    ! output for testing in matlab 
      !minfo%matrix_ia = minfo_matrix_ia
      !minfo%matrix_ja = minfo_matrix_ja
      !minfo%matrix_sa = minfo_matrix_sa
      write(25,*)" coo matrix level ",levelBouss," num rows ",minfo%matrix_nelt
      write(25,909) (j,minfo%matrix_ia(j),minfo%matrix_ja(j),minfo%matrix_sa(j),j=1,minfo%matrix_nelt)
 909  format(3i8,e25.15)
      write(26,*)" rhs size",2*numBoussCells," u and v equations"
      write(26,908) (j,rhs(j),j=1,numBoussCells)
      write(26,908) (j+nD,rhs(j+nD),j=1,numBoussCells)
 908  format(i8,e25.15)
    endif
    
    !============ Debug:
            
    if (.false.) then
        ! check row sums of matrix
        allocate(rowsum(2*nD))
        rowsum = 0.d0
        do k=1,minfo%matrix_nelt
            i = minfo_matrix_ia(k)
            if (minfo%matrix_ja(k) <= nD) then
                ! sums over left half of matrix:
                rowsum(i) = rowsum(i) + minfo%matrix_sa(k)
            endif
            if ((minfo_matrix_ia(k) <= nD) .and. (minfo%matrix_ja(k) > nD)) then
                write(6,*) '+++ k,ia,ja,sa: ',k,minfo_matrix_ia(k),minfo%matrix_ja(k),minfo%matrix_sa(k)
            endif
        enddo
        do i=1,2*nD
            write(6,*) '+++ row ',i,'  rowsum = ',rowsum(i)
        enddo
        
        ! test multiplication:
        ! construct uv vector that varies only in x, uniform in y
        do k = 1,nD
            uv(k) = mod(k,nx)
            uv(k+nD) = 2.d0*mod(k,nx)
            write(6,*) '+++ k,uv(k),uv(k+nD): ',k,uv(k),uv(k+nD)
            enddo
            
        call MATVEC(2*nD, uv, Auv, minfo%matrix_nelt, minfo_matrix_ia, minfo%matrix_ja, minfo%matrix_sa, 0)

        do k=1,nD
            write(6,*) '+++ k,Auv(k),Auv(k+nD): ',k,Auv(k),Auv(k+nD)
            enddo
                
        write(6,*) '+++ Set up matrix of dimension ',2*nD
        write(6,*) '+++ Number of nonzeros: ', minfo%matrix_nelt
        write(6,*) '+++ Size of minfo%matrix_sa: ', size(minfo%matrix_sa)
    endif
    
    ! =============

    if (.false.) then
        allocate(uv(2*nD), Auv(2*nD)) ! for testing matvec

        ! right hand side:
        write(63,*) '========= rhs:'
        do j=1,nD
            write(63,630) j,rhs(j),j+nD,rhs(j+nD)
        enddo
        
 630    format(i6,d26.16, i6,d26.16)
 631    format(4d26.16)
 632    format('Column',i3,'   Row',i3, d26.16, '   Row',i3,d26.16)
        
        ! print columns of matrix:
        do j=1,2*nD
            uv(1:2*nD) = 0.d0
            uv(j) = 1.d0
            
            call MATVEC(2*nD, uv, Auv, minfo%matrix_nelt, minfo_matrix_ia, &
                        minfo%matrix_ja, minfo%matrix_sa, 0)
            !write(63,*) 'Column ',j
            !write(63,631) (Auv(k), k=1,nD)
            !write(63,*) ' '
            !write(63,631) (Auv(k), k=nD+1,2*nD)
            k = 11
            write(63,632) j,k,Auv(k),k+nD,Auv(k+nD)
        enddo
        deallocate(uv,Auv)
    endif
    

    ! =============
    
    
    if (minfo%matrix_nelt > size(minfo%matrix_sa)) then
        write(6,*) '*** Error nelt > size(minfo%matrix_sa), need to allocate larger'
        stop
    endif

    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    timePrepBuild = timePrepBuild+ clock_finish-clock_start
    timePrepBuildCPU = timePrepBuildCPU + cpu_finish-cpu_start
!! debugging
    if (0 .eq. 1) then
    fms = real(clock_finish-clock_start,kind=8)
    write(*,*)" level ",levelBouss
    write(*,400)real(fms,kind=8)/real(clock_rate,kind=8),real(timePrepBuild,kind=8)/real(clock_rate,kind=8)
    write(*,401)cpu_finish-cpu_start,timePrepBuildCPU
 400 format("clock delta ", e15.7," timePrepBuild now ",e15.7)
 401 format("cpu   delta ", e15.7," timePrepBuildCPU  ",e15.7)
    endif
        
#ifdef WHERE_AM_I
  write(*,*) "ending   prepBuildSparseMatrixSGNcoo"
#endif

end subroutine prepBuildSparseMatrixSGNcoo


