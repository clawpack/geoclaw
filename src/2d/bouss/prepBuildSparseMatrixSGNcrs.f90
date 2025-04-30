subroutine prepBuildSparseMatrixSGNcrs(soln,rhs,nvar,naux,levelBouss,numBoussCells)

    
    use amr_module
    use bouss_module
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: sea_level
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8) :: soln(0:2*numBoussCells), rhs(0:2*numBoussCells) 
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: i,j,nng,levSt,nx,ny,mitot,mjtot,loc,locold,locaux,mptr
    integer :: numColsTot

    integer(kind=8) :: clock_start, clock_finish, clock_rate
    real(kind=8) :: cpu_start, cpu_finish
    
    real(kind=8), allocatable, dimension(:) :: uv,Auv,rowsum
    real(kind=8), allocatable, dimension(:,:) :: hx, hy, phi
    real(kind=8), allocatable, dimension(:,:) :: Bx, By, Bxy, Bxx, Byy
    real(kind=8), allocatable, dimension(:,:) :: etax, etay
    
    integer :: nD, k

    external :: MATVEC
    logical :: yhi_db, ylo_db, xlo_db, xhi_db

    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.
    
#ifdef WHERE_AM_I
  write(*,*) "starting prepBuildSparseMatrixSGNcrs"
#endif

    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)
    
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
    
    rhs(:) = 0.d0
    
    ! loop over Bouss patches, must be done sequentially:
    
    levSt  = listStart(levelBouss)
    
!$OMP PARALLEL DO PRIVATE(mptr,nx,ny,mitot,mjtot,loc,locold,locaux),   &
!$OMP             SHARED(levelBouss,numgrids,levSt,alloc,node,nghost,  &
!$OMP                    numBoussCells,nvar,naux,minfo,listOfGrids,    &
!$OMP                    soln,rhs),    &
!$OMP             SCHEDULE(dynamic,1),                                 &
!$OMP             DEFAULT(none)
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

        call buildSparseMatrixSGNcrs(alloc(loc),alloc(locold),         &
                     alloc(locaux),soln,rhs,minfo%rowPtr,minfo%cols,   &
                     minfo%vals, numBoussCells,levelBouss,mptr,nx,ny,nvar,naux)
            
    
    enddo  nng_loop

    ! put last line in CRS row pointers to signify end
    minfo%rowPtr(2*numBoussCells) = minfo%rowPtr(2*numBoussCells-1)+12

    call compressOut(minfo%vals,minfo%rowPtr,minfo%cols,numBoussCells,numColsTot)
    minfo%numColsTot = numColsTot
    minfo%rowPtr(2*numBoussCells) = numColsTot

    ! =============

    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    timePrepBuild = timePrepBuild + clock_finish-clock_start
    timePrepBuildCPU = timePrepBuildCPU + cpu_finish-cpu_start
        
#ifdef WHERE_AM_I
  write(*,*) "ending   prepBuildSparseMatrixSGNcrs"
#endif

end subroutine prepBuildSparseMatrixSGNcrs


