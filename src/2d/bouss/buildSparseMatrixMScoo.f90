subroutine buildSparseMatrixMScoo(rhs,nvar,naux,levelBouss,numBoussCells)

    
    use amr_module
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: sea_level
    use bouss_module
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8) :: rhs(0:2*numBoussCells) 
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: i,j,nng,nb,nelt,levSt,nx,ny,mitot,mjtot,loc,locaux,mptr
    
    real(kind=8), allocatable, dimension(:,:) :: s1, s2, h0etax, h0etay
    real(kind=8), allocatable, dimension(:) :: uv,Auv,rowsum
    
    integer :: nD
    integer :: k, k_ij, k_imj, k_ipj, k_ijm, k_ijp 
    integer :: k_imjp, k_imjm, k_ipjm, k_ipjp
    integer :: kh_ij, kh_imj, kh_ipj, kh_ijm, kh_ijp
    integer :: kh_imjp, kh_imjm, kh_ipjm, kh_ipjp, khu_ij, khv_ij
    integer(kind=8) :: clock_start, clock_finish, clock_rate
    real(kind=8) cpu_start,cpu_finish, time 

    real(kind=8) :: h_ij, h_imj, h_ipj, h_ijm, h_ijp, h2_ij, h3_ij
    real(kind=8) :: h_imjm, h_imjp, h_ipjm, h_ipjp
    real(kind=8) :: hu_ij, hu_imj, hu_ipj, hu_ijm, hu_ijp
    real(kind=8) :: hu_imjm, hu_ipjm, hu_imjp, hu_ipjp
    real(kind=8) :: hv_imjm, hv_ipjm, hv_imjp, hv_ipjp
    real(kind=8) :: b_ij, b_imj, b_ipj, b_ijm, b_ijp
    real(kind=8) :: b_imjm, b_imjp, b_ipjm, b_ipjp
    real(kind=8) :: h0_ij, h0_imj, h0_ipj, h0_ijm, h0_ijp
    real(kind=8) :: h0_imjm, h0_imjp, h0_ipjm, h0_ipjp

    real(kind=8) :: hv_ij, hv_imj, hv_ipj, hv_ijm, hv_ijp
    real(kind=8) :: eta_ij, eta_imj, eta_ipj, eta_ijm, eta_ijp
    real(kind=8) :: dx, dy, dxm, dym, dxdx, dydy, dxdy4
    
    integer :: nelt_ij, nelt_ij_nD, nelt_imj, nelt_ipj, nelt_ijm_nD, nelt_ijp_nD
    real(kind=8) :: sa, huc, hvc, matrix_rhs1_max, matrix_rhs2_max
    integer :: ja, bc, c
    real(kind=8) :: b_min, b_max

    external :: MATVEC
    logical :: debug, noxycorner, noxy
    logical :: yhi_db, ylo_db, xlo_db, xhi_db

    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.

#ifdef WHERE_AM_I
  write(*,*) "starting buildSparseMatrixMScoo"
#endif
    
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)

    debug = .false.
    !debug = .true.
    
    noxy = .false.  ! omit all cross derivative terms (test stability)
    noxycorner = .false.  ! omit cross derivative terms only at patch corners
    
    minfo => matrix_info_allLevs(levelBouss)

!
! ================== Step 2  set up for matrix elements  ======================

    ! Set up matrix A and RHS by looping over all Bouss patches:
    
    dx = hxposs(levelBouss)
    dy = hyposs(levelBouss)
    
     
    if (coordinate_system .eq. 2) then
      ! convert from dy in degrees to meters:
      dym = dy * earth_radius * deg2rad
      ! dxm depends on latitude, done in loops below instead of here
    else
      dym = dy
    endif
    
    dydy = dym**2
    

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
    
    nelt = 0  ! count number of matrix elements as added, 
              ! should be less than 9*(2*nD)**2

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
        locaux = node(storeaux,mptr)  ! needed for iaddaux
        
        ! set booleans indicating whether each patch edge is a Domain Boundary:
        !call setDomainBndries(mptr,yhi_db,ylo_db,xlo_db,xhi_db)
        xlo_db = (node(ndilo,mptr) .eq. 0)
        xhi_db = (node(ndihi,mptr) .eq. (iregsz(levelBouss)-1))
        ylo_db = (node(ndjlo,mptr) .eq. 0)
        yhi_db = (node(ndjhi,mptr) .eq. (jregsz(levelBouss)-1))
    
        nb = node(bpatchNum, mptr)
        
        ! Get integer array providing matrix indices mapping 
        ! cell coordinates i,j to location k in vector of all unknowns: 
           
        mi => minfo%matrix_indices(nb)  ! integer array for this patch

        ! allocate temp arrays needed for this patch:
        ! provide one layer of ghost cells needed when computing
        ! second differences of these quantities, but values in ghost
        ! cells are never actually needed since second diff's won't be
        ! needed in first row of interior cells
        
         allocate(s1(0:nx+1,0:ny+1), s2(0:nx+1,0:ny+1), &
                  h0etax(0:nx+1,0:ny+1), h0etay(0:nx+1,0:ny+1))


        s1(0:nx+1,0:ny+1) = 0
        s2(0:nx+1,0:ny+1) = 0
                                  
        ! First loop over all cells (i,j), one row of ghost cells,
        ! computing s1, s2, h0etax, h0etay,
        ! quantities that involve first derivatives.  In the second loop
        ! over all cells below we will compute matrix_rhs using second
        ! differences of these quantities.
        
        if (debug) then
            !write(*,*) "in buildsparsematrix: calling uvsym for grid ",mptr
            !call uvsymcheck5(alloc(loc),0,mitot,mjtot,nvar)
            do j=-1,ny+2
                do i=-1,nx+2
                    kh_ij = iadd(1,i,j)
                    write(*,900) i,j,alloc(kh_ij),          &
                           alloc(kh_ij+1),alloc(kh_ij+2),   &
                           alloc(kh_ij+3),alloc(kh_ij+4)
 900   format('+++ i,j,h,hu,hv,huc,hvc: ',2i4,5e20.12)
                    enddo
                enddo
        endif
            
         do j=0,ny+1 
            if (coordinate_system .eq. 2) then
              ! aux(3,i,j) is latitude, doesn't depend on i so use i=1
              dxm = dx * earth_radius * deg2rad *              &
                    cos(alloc(iaddaux(3,1,j)))
            else
              dxm = dx
            endif
            
            do i=0,nx+1
                ! alloc indices for h component q(1,i,j) and neighbors:
                kh_ij = iadd(1,i,j)
                kh_imj = iadd(1,i-1,j)
                kh_ipj = iadd(1,i+1,j)
                kh_ijm = iadd(1,i,j-1)
                kh_ijp = iadd(1,i,j+1)
                
                ! use the fact that h,hu,hv are contiguous in alloc, and
                ! also compute eta = h + topo where topo is in aux(1,i,j)

                h_ij = alloc(kh_ij)
                hu_ij = alloc(kh_ij+1)
                hv_ij = alloc(kh_ij+2)
                b_ij = alloc(iaddaux(1,i,j))
                if (b_ij < sea_level) then
                    eta_ij = h_ij + b_ij
                    h0_ij = sea_level - b_ij
                else
                    eta_ij = 0.d0
                    h0_ij = 0.d0
                endif
                
                h_ipj = alloc(kh_ipj)
                hu_ipj = alloc(kh_ipj+1)
                hv_ipj = alloc(kh_ipj+2)
                b_ipj = alloc(iaddaux(1,i+1,j))
                if (b_ipj < sea_level) then
                    eta_ipj = h_ipj + b_ipj
                    h0_ipj = sea_level - b_ipj
                else
                    eta_ipj = 0.d0
                    h0_ipj = 0.d0
                endif

                h_imj = alloc(kh_imj)
                hu_imj = alloc(kh_imj+1)
                hv_imj = alloc(kh_imj+2)
                b_imj = alloc(iaddaux(1,i-1,j))
                if (b_imj < sea_level) then
                    eta_imj = h_imj + b_imj
                    h0_imj = sea_level - b_imj
                else
                    eta_imj = 0.d0
                    h0_imj = 0.d0
                endif

                h_ijp = alloc(kh_ijp)
                hu_ijp = alloc(kh_ijp+1)
                hv_ijp = alloc(kh_ijp+2)
                b_ijp = alloc(iaddaux(1,i,j+1))
                if (b_ijp < sea_level) then
                    eta_ijp = h_ijp + b_ijp
                    h0_ijp = sea_level - b_ijp
                else
                    eta_ijp = 0.d0
                    h0_ijp = 0.d0
                endif
                
                h_ijm = alloc(kh_ijm)
                hu_ijm = alloc(kh_ijm+1)
                hv_ijm = alloc(kh_ijm+2)
                b_ijm = alloc(iaddaux(1,i,j-1))
                if (b_ijm < sea_level) then
                    eta_ijm = h_ijm + b_ijm
                    h0_ijm = sea_level - b_ijm
                else
                    eta_ijm = 0.d0
                    h0_ijm = 0.d0
                endif

               ! only modify s1,s2 from initialized values of 0 if cell i,j
               ! and all four neighbors are wet:
                                                            
                if (.not.(h_ij  .eq. 0. .or.                                  &
                          h_ipj .eq. 0. .or. h_ijm .eq. 0.d0 .or.             &
                          h_ijp .eq. 0. .or. h_imj .eq. 0.d0))     then
                
                   s1(i,j) = (hu_ipj**2/h_ipj - hu_imj**2/h_imj        &
                         + grav*h_ij*(eta_ipj - eta_imj)) / (2.d0*dxm) &
                         + (hu_ijp*hv_ijp/h_ijp - hu_ijm*hv_ijm/h_ijm) &
                           / (2.d0*dym)
                           
                   ! Debug: remove nonlinear terms:
                   !s1(i,j) = grav*h_ij*(eta_ipj - eta_imj) / (2.d0*dxm)
                   !s1(i,j) = (hu_ipj**2/h_ipj - hu_imj**2/h_imj) / (2.d0*dxm)
                           
                   s2(i,j) = (hv_ijp**2/h_ijp - hv_ijm**2/h_ijm        &
                         + grav*h_ij*(eta_ijp - eta_ijm)) / (2.d0*dym) &
                         + (hu_ipj*hv_ipj/h_ipj - hu_imj*hv_imj/h_imj) &
                           / (2.d0*dxm)
                           
                   ! Debug: remove nonlinear terms:
                   !s2(i,j) = 0.d0 !grav*h_ij*(eta_ijp - eta_ijm) / (2.d0*dym)
                   
                   !write(67,*) '+++ i,j,s1,s2: ',i,j,s1(i,j),s2(i,j)
                           
                endif
                           
                h0etax(i,j) = h0_ij * (eta_ipj - eta_imj) / (2.d0*dxm)
                h0etay(i,j) = h0_ij * (eta_ijp - eta_ijm) / (2.d0*dym)
                !write(67,*) '+++ i,j,h0etax,h0etay: ',i,j,h0etax(i,j),h0etay(i,j)
                
                if (debug) then
                    write(67,*) '  h_ij, eta_ij: ',  h_ij,eta_ij
                    write(67,*) '  h0etax:       ',  h0etax(i,j)
                endif
                
            enddo
        enddo

! ================== build matrix and rhs
        ! Second loop over all cells.  ONLY INTERIOR OK NOW
        
        do j=1,ny 
            if (coordinate_system .eq. 2) then
              ! aux(3,i,j) is latitude, doesn't depend on i so use i=1
              dxm = dx * earth_radius * deg2rad *              &
                    cos(alloc(iaddaux(3,1,j)))
            else
              dxm = dx
            endif
            
            dxdx = dxm**2
            dxdy4 = 4.d0*dxm*dym
            
            do i=1,nx

                ! sparse matrix row indices for cell (i,j) are:
                !    k_ij for the row of I-D11 and -D12 in block matrix,
                !    k_ij + nD for the row of -D21 and I-D22
                k_ij = mi%mindex(i,j)
                
                ! in this inner loop we set all nonzero matrix elements
                ! in these two rows of the matrix (and rhs) corresponding
                ! to momenumtum updates to hu and hv in the (i,j) cell
                
                ! sparse matrix column indices for neighboring cells
                ! are the values k computed next (for I-D1 and -D3)
                ! and the corresponding k + nD for the corresponding column
                ! of the -D2 and I-D4.  These indices are needed to
                ! set matrix elements in these columns, corresponding to
                ! momentum components hu and hv in the neighboring cells
                
                k_imj = mi%mindex(i-1,j)
                k_ipj = mi%mindex(i+1,j)
                k_ijm = mi%mindex(i,j-1)
                k_ijp = mi%mindex(i,j+1)
                k_imjm = mi%mindex(i-1,j-1)
                k_imjp = mi%mindex(i-1,j+1)
                k_ipjm = mi%mindex(i+1,j-1)
                k_ipjp = mi%mindex(i+1,j+1)
                
                ! alloc indices for h component q(1,i,j) and neighbors:
                kh_ij = iadd(1,i,j)
                kh_imj = iadd(1,i-1,j)
                kh_ipj = iadd(1,i+1,j)
                kh_ijm = iadd(1,i,j-1)
                kh_ijp = iadd(1,i,j+1)
                kh_imjm = iadd(1,i-1,j-1)
                kh_imjp = iadd(1,i-1,j+1)
                kh_ipjm = iadd(1,i+1,j-1)
                kh_ipjp = iadd(1,i+1,j+1)
                
                ! use the fact that h,hu,hv are contiguous in alloc:
                h_ij = alloc(kh_ij)
                hu_ij = alloc(kh_ij+1)
                hv_ij = alloc(kh_ij+2)
                b_ij = alloc(iaddaux(1,i,j))
                h0_ij = sea_level - b_ij
                
                h_ipj = alloc(kh_ipj)
                hu_ipj = alloc(kh_ipj+1)
                hv_ipj = alloc(kh_ipj+2)
                b_ipj = alloc(iaddaux(1,i+1,j))
                h0_ipj = sea_level - b_ipj

                h_imj = alloc(kh_imj)
                hu_imj = alloc(kh_imj+1)
                hv_imj = alloc(kh_imj+2)
                b_imj = alloc(iaddaux(1,i-1,j))
                h0_imj = sea_level - b_imj

                h_ijp = alloc(kh_ijp)
                hu_ijp = alloc(kh_ijp+1)
                hv_ijp = alloc(kh_ijp+2)
                b_ijp = alloc(iaddaux(1,i,j+1))
                h0_ijp = sea_level - b_ijp
                
                h_ijm = alloc(kh_ijm)
                hu_ijm = alloc(kh_ijm+1)
                hv_ijm = alloc(kh_ijm+2)
                b_ijm = alloc(iaddaux(1,i,j-1))
                h0_ijm = sea_level - b_ijm

                h_ipjp = alloc(kh_ipjp)
                hu_ipjp = alloc(kh_ipjp+1)
                hv_ipjp = alloc(kh_ipjp+2)
                b_ipjp = alloc(iaddaux(1,i+1,j+1))
                h0_ipjp = sea_level - b_ipjp
                
                h_ipjm = alloc(kh_ipjm)
                hu_ipjm = alloc(kh_ipjm+1)
                hv_ipjm = alloc(kh_ipjm+2)
                b_ipjm = alloc(iaddaux(1,i+1,j-1))
                h0_ipjm = sea_level - b_ipjm
                
                h_imjp = alloc(kh_imjp)
                hu_imjp = alloc(kh_imjp+1)
                hv_imjp = alloc(kh_imjp+2)
                b_imjp = alloc(iaddaux(1,i-1,j+1))
                h0_imjp = sea_level - b_imjp
                
                h_imjm = alloc(kh_imjm)
                hu_imjm = alloc(kh_imjm+1)
                hv_imjm = alloc(kh_imjm+2)
                b_imjm = alloc(iaddaux(1,i-1,j-1))
                h0_imjm = sea_level - b_imjm

                ! quantities repeatedly used:
                h2_ij = h0_ij**2 * (Bparam + 0.5d0)
                h3_ij = h0_ij**3 / 6.d0
                
                
                if (k_ij == -1) then
                    write(6,*) '+++ unexpected k_ij=-1 at i,j: ',i,j
                endif
                
                b_min = min(b_ij,b_imj,b_ipj,b_ijm,b_ijp,              &
                        b_imjm, b_imjp, b_ipjm, b_ipjp)
                b_max = max(b_ij,b_imj,b_ipj,b_ijm,b_ijp,              &
                        b_imjm, b_imjp, b_ipjm, b_ipjp)                        
                if ((-b_min < boussMinDepth) .or. (b_max >= 0.d0)) then
                        
                    ! The deepest water in block of 9 cells is shallower 
                    ! than boussMinDepth, or else one of the cells is onshore.

                    ! Revert to SWE in this cell, i.e. two rows
                    ! of matrix should be rows of identity with rhs = 0,
                    ! so correction terms huc,hvc calculated will be 0.
                    
                    ! Replace diagonal block I-D11 by I matrix:
                    nelt = nelt + 1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = k_ij
                    minfo%matrix_sa(nelt) = 1.d0
                    !rhs(k_ij) = 0.d0  ! Not needed due to initialization
                    
                    ! Replace diagonal block I-D22 by I matrix:
                    nelt = nelt + 1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = k_ij + nD
                    minfo%matrix_sa(nelt) = 1.d0
                    !rhs(k_ij + nD) = 0.d0  ! Not needed due to initialization
                    
                    ! -D12 and -D21 are the 0 matrix in this case (no nozeros)
                    
                    cycle ! go on to next i
                endif

                
                ! If we get here, we are dealing with interior point
                ! Add rows k_ij and (k_ij + nD) to matrix:
                ! Row k_ij    corresponds to correction to hu in cell (i,j)
                ! Row k_ij+nD corresponds to correction to hv in cell (i,j)
                ! The nonzero matrix elements will be in the same columns
                ! and also in columns 
                !     k_ipj, k_ipj+nD for the coeff of the (i+1,j) stencil pt
                !     k_imjp, k_imjp+nD for the coeff of (i-1,j+1) point,
                !     etc.
                ! Note that if e.g. k_ipj = -1 then cell (i+1,j) is not in the
                ! union of patches at this level and so there is no 
                ! corresponding column.  
                ! In this case:
                ! For Dirichlet BC with value 0, the stencil value is discarded,
                ! For Neumann BC, the stencil coefficient is instead applied
                ! to the neighboring cell that is in the patch, 
                ! assuming constant extrapolation at each edge.
                
                ! Increment nelt each time we add another nonzero matrix
                ! element.  
                ! Store index nelt_ij for index to diagonal element in row k_ij
                ! Store index nelt_ij_nD for diagonal element in row k_ij+nD

                ! initialize to -1 for debugging in case not set properly:
                nelt_ij = -1
                nelt_imj = -1
                nelt_ipj = -1
                nelt_ij_nD = -1
                nelt_ijm_nD = -1
                nelt_ijp_nD = -1
                
  601           format('U(',i2,',',i2,') --> U(',i2,',',i2,') c=',i1,' nelt=',i8,' sa = ',e25.15)
  602           format('U(',i2,',',i2,') --> V(',i2,',',i2,') c=',i1,' nelt=',i8,' sa = ',e25.15)
  603           format('V(',i2,',',i2,') --> U(',i2,',',i2,') c=',i1,' nelt=',i8,' sa = ',e25.15)
  604           format('V(',i2,',',i2,') --> V(',i2,',',i2,') c=',i1,' nelt=',i8,' sa = ',e25.15)
    
                !-----------------------------------------------
                ! I-D11 block (xx derivatives of hu components):
                
                ! Matrix element row U(i,j) column U(i,j):
                sa = 1.d0 + 2.d0*(h2_ij - h3_ij/h0_ij) / dxdx
                nelt = nelt+1
                minfo%matrix_ia(nelt) = k_ij
                minfo%matrix_ja(nelt) = k_ij
                minfo%matrix_sa(nelt) = sa
                if (debug) write(67,601) i,j,i,j,1,nelt,minfo%matrix_sa(nelt)
                
                nelt_ij = nelt  ! index of diagonal matrix element
                    
                ! Matrix element row U(i,j) column U(i+1,j):
                sa = -(h2_ij - h3_ij/h0_ipj) / dxdx
                if (k_ipj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = k_ipj
                    minfo%matrix_sa(nelt) = sa
                    nelt_ipj = nelt
                    if (debug) write(67,601) i,j,i+1,j,1,nelt,sa
                else
                    ! right boundary of patch
                    if (bc_xhi==0) then
                        ! for 0 Dirchlet BC don't any element
                    else if (bc_xhi==1 .and. xhi_db) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij) = minfo%matrix_sa(nelt_ij) + sa
                        if (debug) write(67,601) i,j,i+1,j,2,nelt_ij,sa
                    else if (bc_xhi==3 .and. xhi_db) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij) = minfo%matrix_sa(nelt_ij) - sa
                        if (debug) write(67,601) i,j,i+1,j,2,nelt_ij,-sa
                    else if (bc_xhi==2 .or. .not. xhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = alloc(iadd(4,i+1,j))
                        rhs(k_ij) = rhs(k_ij) - sa*huc
                        if (debug) write(67,601) i,j,i+1,j,5,-1,sa*huc
                    endif
                    
                    
                endif
                
                ! Matrix element row U(i,j) column U(i-1,j):
                sa = -(h2_ij - h3_ij/h0_imj) / dxdx
                if (k_imj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = k_imj
                    minfo%matrix_sa(nelt) = sa
                    nelt_imj = nelt
                    if (debug) write(67,601) i,j,i-1,j,1,nelt,sa
                else
                    ! left boundary of patch
                    if (bc_xlo==0) then
                        ! for 0 Dirchlet BC don't add any element
                    else if (bc_xlo==1 .and. xlo_db) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij) = minfo%matrix_sa(nelt_ij) + sa
                        if (debug) write(67,601) i,j,i-1,j,2,nelt_ij,sa
                    else if (bc_xlo==3 .and. xlo_db) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij) = minfo%matrix_sa(nelt_ij) - sa
                        if (debug) write(67,601) i,j,i-1,j,2,nelt_ij,-sa
                    else if (bc_xlo==2 .or. .not. xlo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = alloc(iadd(4,i-1,j))
                        rhs(k_ij) = rhs(k_ij) - sa*huc
                        if (debug) write(67,601) i,j,i-1,j,5,-1,sa*huc
                    endif
                endif

                !-----------------------------------------------
                ! I-D22 block (yy derivatives of hv components):
                
                ! Matrix element row V(i,j) column V(i,j):
                sa = 1.d0 + 2.d0*(h2_ij - h3_ij/h0_ij) / dydy
                nelt = nelt+1
                minfo%matrix_ia(nelt) = k_ij + nD
                minfo%matrix_ja(nelt) = k_ij + nD
                minfo%matrix_sa(nelt) = sa
                if (debug) write(67,604) i,j,i,j,1,nelt,minfo%matrix_sa(nelt)
                
                nelt_ij_nD = nelt  ! index of diagonal matrix element
                    
                ! Matrix element row V(i,j) column V(i,j+1):
                sa = -(h2_ij - h3_ij/h0_ijp) / dydy
                if (k_ijp > -1) then
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = k_ijp + nD
                    minfo%matrix_sa(nelt) = sa
                    nelt_ijp_nD = nelt
                    if (debug) write(67,604) i,j,i,j+1,1,nelt,sa
                else
                    ! top boundary of patch
                    if (bc_yhi==0) then
                        ! for 0 Dirchlet BC don't add any element
                    else if (bc_yhi==1 .and. yhi_db) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij_nD) = &
                                minfo%matrix_sa(nelt_ij_nD) + sa
                        if (debug) write(67,604) i,j,i,j+1,2,nelt_ij_nD,sa
                    else if (bc_yhi==3 .and. yhi_db) then
                        ! For wall BC add -sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij_nD) = &
                                minfo%matrix_sa(nelt_ij_nD) - sa
                        if (debug) write(67,604) i,j,i,j+1,2,nelt_ij_nD,-sa
                    else if (bc_yhi==2 .or. .not. yhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = alloc(iadd(5,i,j+1))
                        rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                        if (debug) write(67,604) i,j,i,j+1,5,-1,sa*hvc
                    endif
                endif
                
                ! Matrix element row V(i,j) column V(i,j-1):
                sa = -(h2_ij - h3_ij/h0_ijm) / dydy
                if (k_ijm > -1) then
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = k_ijm + nD
                    minfo%matrix_sa(nelt) = sa
                    nelt_ijm_nD = nelt
                    if (debug) write(67,604) i,j,i,j-1,1,nelt,sa
                else
                    ! bottom boundary of patch
                    if (bc_ylo==0) then
                        ! for 0 Dirchlet BC don't add any element
                    else if (bc_ylo==1 .and. ylo_db) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij_nD) = &
                                minfo%matrix_sa(nelt_ij_nD) + sa
                        if (debug) write(67,604) i,j,i,j-1,2,nelt_ij_nD,sa
                    else if (bc_ylo==3 .and. ylo_db) then
                        ! For wall BC add -sa to diagonal matrix element:
                        minfo%matrix_sa(nelt_ij_nD) = &
                                minfo%matrix_sa(nelt_ij_nD) - sa
                        if (debug) write(67,604) i,j,i,j-1,2,nelt_ij_nD,-sa
                    else if (bc_ylo==2 .or. .not. ylo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = alloc(iadd(5,i,j-1))
                        rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                        if (debug) write(67,604) i,j,i,j-1,5,-1,sa*hvc
                    endif   
                endif                                        

                if (noxy) go to 111  ! omit all cross derivatives
            
                !-----------------------------------------------
                ! -D12 block (cross derivatives):
                ! row ia = k_ij for hu corrections in cell (i,j)
                ! columns k_ipjp + nD and other neighbors for hv corrections                
                
                ! Matrix element row U(i,j) column V(i+1,j+1):
                sa = -(h2_ij - h3_ij/h0_ipjp) / dxdy4
                ja = -1
                c = 0
                if (k_ipjp > -1) then
                    ja = k_ipjp + nD
                    c = 1
                else if (noxycorner .and. (k_ijp == -1) .and. (k_ipj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijp == -1) .and. ((bc_yhi==2).or. .not. yhi_db)) .or. &
                         ((k_ipj == -1) .and. ((bc_xhi==2).or. .not. xhi_db))) then
                    hvc = alloc(iadd(5,i+1,j+1))
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                    if (debug) write(67,602) i,j,i+1,j+1,5,-1,sa*hvc
                else if (k_ijp > -1) then
                    ! at right boundary, use cell above instead
                    if (bc_xhi /= 2) ja = k_ijp + nD
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead,
                    ja = k_ipj + nD
                    if (bc_yhi==3) sa = -sa
                else
                    ! top-right corner
                    ja = k_ij + nD
                    if (bc_yhi==3) sa = -sa  ! v negated at wall
                endif

                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,602) i,j,i+1,j+1,c,nelt,sa
                endif
                
                    
                ! Matrix element row U(i,j) column V(i+1,j-1):
                sa = (h2_ij - h3_ij/h0_ipjm) / dxdy4
                ja = -1
                c = 0
                if (k_ipjm > -1) then
                    ja = k_ipjm + nD
                    c = 1
                else if (noxycorner .and. (k_ijm == -1) .and. (k_ipj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs                
                else if (((k_ijm == -1) .and. ((bc_ylo==2).or. .not. ylo_db)) .or. &
                         ((k_ipj == -1) .and. ((bc_xhi==2).or. .not. xhi_db))) then
                    hvc = alloc(iadd(5,i+1,j-1))
                    rhs(k_ij) = rhs(k_ij) - sa*hvc   
                    if (debug) write(67,602) i,j,i+1,j-1,5,-1,sa*hvc             
                else if (k_ijm > -1) then
                    ! at right boundary, use cell below instead
                    if (bc_xhi /= 2) ja = k_ijm + nD
                else if (k_ipj > -1) then
                    ! at bottom boundary, use cell to right instead,
                    ja = k_ipj + nD
                    if (bc_ylo==3) sa = -sa
                else
                    ! bottom-right corner
                    ja = k_ij + nD
                    if (bc_ylo==3) sa = -sa  ! v negated at wall
                endif   
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,602) i,j,i+1,j-1,c,nelt,sa
                endif
                
                ! Matrix element row U(i,j) column V(i-1,j+1):    
                sa = (h2_ij - h3_ij/h0_imjp) / dxdy4
                ja = -1
                c = 0
                if (k_imjp > -1) then
                    ja = k_imjp + nD
                    c = 1
                else if (noxycorner .and. (k_ijp == -1) .and. (k_imj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijp == -1) .and. ((bc_yhi==2).or. .not. yhi_db)) .or. &
                         ((k_imj == -1) .and. ((bc_xlo==2).or. .not. xlo_db))) then
                    hvc = alloc(iadd(5,i-1,j+1))
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                    if (debug) write(67,602) i,j,i-1,j+1,5,-1,sa*hvc
                else if (k_ijp > -1) then
                    ! at left boundary, use cell above instead
                    if (bc_xlo /= 2) ja = k_ijp + nD
                else if (k_imj > -1) then
                    ! at top boundary, use cell to left instead,
                    ja = k_imj + nD
                    if (bc_yhi==3) sa = -sa
                else
                    ! top-left corner
                    ja = k_ij + nD
                    if (bc_yhi==3) sa = -sa  ! v negated at wall
                endif
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,602) i,j,i-1,j+1,c,nelt,sa
                endif
                
                ! Matrix element row U(i,j) column V(i-1,j-1):
                sa = -(h2_ij - h3_ij/h0_imjm) / dxdy4
                ja = -1
                c = 0
                if (k_imjm > -1) then
                    ja = k_imjm + nD
                    c = 1
                else if (noxycorner .and. (k_ijm == -1) .and. (k_imj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijm == -1) .and. ((bc_ylo==2).or. .not. ylo_db)) .or. &
                         ((k_imj == -1) .and. ((bc_xlo==2).or. .not. xlo_db))) then
                    hvc = alloc(iadd(5,i-1,j-1))
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                    if (debug) write(67,602) i,j,i-1,j-1,5,-1,sa*hvc
                else if (k_ijm > -1) then
                    ! at left boundary, use cell below instead
                    if (bc_xlo /= 2) ja = k_ijm + nD
                else if (k_imj > -1) then
                    ! at bottom boundary, use cell to left instead,
                    ja = k_imj + nD
                    if (bc_ylo==3) sa = -sa
                else
                    ! bottom-left corner
                    ja = k_ij + nD
                    if (bc_ylo==3) sa = -sa  ! v negated at wall
                endif   
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,602) i,j,i-1,j-1,c,nelt,sa
                endif


                !-----------------------------------------------
                ! -D21 block (cross derivatives):
                ! row ia = k_ij + nD  for hv corrections in cell (i,j)
                ! columns k_ipjp and other neighbors for hu corrections
                
                ! Matrix element row V(i,j) column U(i+1,j+1):
                sa = -(h2_ij - h3_ij/h0_ipjp) / dxdy4
                ja = -1
                c = 0
                if (k_ipjp > -1) then
                    ja = k_ipjp
                    c = 1
                else if (noxycorner .and. (k_ijp == -1) .and. (k_ipj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijp == -1) .and. ((bc_yhi==2).or. .not. yhi_db)) .or. &
                         ((k_ipj == -1) .and. ((bc_xhi==2).or. .not. xhi_db))) then
                    huc = alloc(iadd(4,i+1,j+1))
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                    if (debug) write(67,603) i,j,i+1,j+1,5,-1,sa*huc
                else if (k_ijp > -1) then
                    ! at right boundary, use cell above instead,
                    ja = k_ijp
                    if (bc_xhi==3) sa = -sa
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead,
                    if (bc_yhi /= 2) ja = k_ipj
                else
                    ! top-right corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! u negated at wall
                endif
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,603) i,j,i+1,j+1,c,nelt,sa
                endif
                
                ! Matrix element row V(i,j) column U(i+1,j-1):
                sa = (h2_ij - h3_ij/h0_ipjm) / dxdy4
                ja = -1
                c = 0
                if (k_ipjm > -1) then
                    ja = k_ipjm
                    c = 1
                else if (noxycorner .and. (k_ijm == -1) .and. (k_ipj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijm == -1) .and. ((bc_ylo==2).or. .not. ylo_db)) .or. &
                         ((k_ipj == -1) .and. ((bc_xhi==2).or. .not. xhi_db))) then
                    huc = alloc(iadd(4,i+1,j-1))
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                    if (debug) write(67,603) i,j,i+1,j-1,5,-1,sa*huc
                else if (k_ijm > -1) then
                    ! at right boundary, use cell below instead,
                    ja = k_ijm
                    if (bc_xhi==3) sa = -sa
                else if (k_ipj > -1) then
                    ! at bottom boundary, use cell to right instead,
                    if (bc_ylo /= 2) ja = k_ipj
                else
                    ! bottom-right corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! u negated at wall
                endif
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,603) i,j,i+1,j-1,c,nelt,sa
                endif
                
                ! Matrix element row V(i,j) column U(i-1,j+1):
                sa = (h2_ij - h3_ij/h0_imjp) / dxdy4
                ja = -1
                c = 0
                if (k_imjp > -1) then
                    ja = k_imjp
                    c = 1
                else if (noxycorner .and. (k_ijp == -1) .and. (k_imj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijp == -1) .and. ((bc_yhi==2).or. .not. yhi_db)) .or. &
                         ((k_imj == -1) .and. ((bc_xlo==2).or. .not. xlo_db))) then
                    huc = alloc(iadd(4,i-1,j+1))
                    rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                    if (debug) write(67,603) i,j,i-1,j+1,5,-1,sa*huc
                else if (k_ijp > -1) then
                    ! at left boundary, use cell above instead,
                    ja = k_ijp
                    if (bc_xlo==3) sa = -sa
                else if (k_imj > -1) then
                    ! at top boundary, use cell to left instead,
                    if (bc_yhi /= 2) ja = k_imj
                else
                    ! top-left corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! u negated at wall
                endif
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,603) i,j,i-1,j+1,c,nelt,sa
                endif

                ! Matrix element row V(i,j) column U(i-1,j-1):
                sa = -(h2_ij - h3_ij/h0_imjm) / dxdy4
                ja = -1
                c = 0
                if (k_imjm > -1) then
                    ja = k_imjm
                    c = 1
                else if (noxycorner .and. (k_ijm == -1) .and. (k_imj == -1)) then
                    ! test omitting corner cross derivatives for stability?
                    ! in this case don't add sa to matrix or rhs
                else if (((k_ijm == -1) .and. ((bc_ylo==2).or. .not. ylo_db)) .or. &
                         ((k_imj == -1) .and. ((bc_xlo==2).or. .not. xlo_db))) then
                    huc = alloc(iadd(4,i-1,j-1))
                    rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                    if (debug) write(67,603) i,j,i-1,j-1,5,-1,sa*huc
                else if (k_ijm > -1) then
                    ! at left boundary, use cell below instead,
                    ja = k_ijm
                    if (bc_xlo==3) sa = -sa
                else if (k_imj > -1) then
                    ! at bottom boundary, use cell to left instead,
                    ja = k_imj
                else
                    ! bottom-left corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! u negated at wall
                endif
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo%matrix_ia(nelt) = k_ij + nD
                    minfo%matrix_ja(nelt) = ja
                    minfo%matrix_sa(nelt) = sa
                    if (debug) write(67,603) i,j,i-1,j-1,c,nelt,sa
                endif
                
 111            continue
 
                !-------------------------------------------------------
                ! Right hand sides for vector elements k_ij and k_ij+nD
                ! corresponding to b1 and b2 in notes, where
                !    b = -D*s + g*B1*h0**2 * \grad(\grad\cdot(h0*\grad\eta))
                ! with
                !    s = z + g*h*\grad\eta
                !    z = u\grad\cdot(h*u) + (h*u\cdot\grad)u
                !        are the nonlinear terms of SWE
                ! Recall that:
                !    h2_ij = h0_ij**2 * (Bparam + 0.5d0)
                !    h3_ij = h0_ij**3 / 6.d0
                !
                ! s1,s2 are components of   s = z + g*h*\grad\eta
                ! h0etax, h0etay approximate (h0 * eta)_x, (h0 * eta)_y
                !
                ! s1,s2,h0etax,h0etay were all computed in previous loop
                ! over i,j, so we can use neighboring values to now compute
                ! second differences of these.
                
                ! DEBUG:
                !if (abs(rhs(k_ij)) > 0.d0) then
                !    write(6,*) '+++ rhs huc nonzero: ',i,j,rhs(k_ij)
                !endif
                !
                !if (abs(rhs(k_ij+nD)) > 0.d0) then
                !    write(6,*) '+++ rhs hvc nonzero: ',i,j,rhs(k_ij+nD)
                !endif
                

                
              
              if (noxy) then
                ! leave out cross derivative terms
                ! b1, rhs for hu updates:
                rhs(k_ij) = rhs(k_ij) &
                  - h2_ij * (s1(i-1,j) - 2.d0*s1(i,j) + s1(i+1,j)) / dxdx &
                  + h3_ij * (s1(i-1,j)/h0_imj - 2.d0*s1(i,j)/h0_ij &
                          + s1(i+1,j)/h0_ipj) / dxdx &
                  + Bparam * grav * h0_ij**2 &
                      * (h0etax(i-1,j) - 2.d0*h0etax(i,j) + h0etax(i+1,j)) / dxdx 
                          
                ! b2, rhs for hv updates:      
                rhs(k_ij+nD) = rhs(k_ij+nD) &
                  - h2_ij * (s2(i,j-1) - 2.d0*s2(i,j) + s2(i,j+1)) / dydy &
                  + h3_ij * (s2(i,j-1)/h0_ijm - 2.d0*s2(i,j)/h0_ij &
                          + s2(i,j+1)/h0_ijp) / dydy &
                  + Bparam * grav * h0_ij**2 &
                      * (h0etay(i,j-1) - 2.d0*h0etay(i,j) + h0etay(i,j+1)) / dydy
                      
              else
                ! full thing with cross derivatives:             
                ! b1, rhs for hu updates:
                rhs(k_ij) = rhs(k_ij) &
                    -h2_ij * ((s1(i-1,j) - 2.d0*s1(i,j) + s1(i+1,j)) / dxdx &
                          + (s2(i+1,j+1) - s2(i-1,j+1) &
                           - s2(i+1,j-1) + s2(i-1,j-1)) / dxdy4) &
                  + h3_ij * ((s1(i-1,j)/h0_imj - 2.d0*s1(i,j)/h0_ij &
                          + s1(i+1,j)/h0_ipj) / dxdx &
                          + (s2(i+1,j+1)/h0_ipjp - s2(i-1,j+1)/h0_imjp &
                            -s2(i+1,j-1)/h0_ipjm + s2(i-1,j-1)/h0_imjm) / dxdy4) &
                  + Bparam * grav * h0_ij**2 &
                      * ((h0etax(i-1,j) - 2.d0*h0etax(i,j) + h0etax(i+1,j)) / dxdx &
                        + (h0etay(i+1,j+1) - h0etay(i-1,j+1) &
                          -h0etay(i+1,j-1) + h0etay(i-1,j-1)) / dxdy4)
                          
                            
                ! b2, rhs for hv updates:      
                rhs(k_ij+nD) = rhs(k_ij+nD) &
                    -h2_ij * ((s2(i,j-1) - 2.d0*s2(i,j) + s2(i,j+1)) / dydy &
                          + (s1(i+1,j+1) - s1(i-1,j+1) &
                           - s1(i+1,j-1) + s1(i-1,j-1)) / dxdy4) &
                  + h3_ij * ((s2(i,j-1)/h0_ijm - 2.d0*s2(i,j)/h0_ij &
                          + s2(i,j+1)/h0_ijp) / dydy &
                          + (s1(i+1,j+1)/h0_ipjp - s1(i-1,j+1)/h0_imjp &
                            -s1(i+1,j-1)/h0_ipjm + s1(i-1,j-1)/h0_imjm) / dxdy4) &
                  + Bparam * grav * h0_ij**2 &
                      * ((h0etay(i,j-1) - 2.d0*h0etay(i,j) + h0etay(i,j+1)) / dydy &
                        + (h0etax(i+1,j+1) - h0etax(i-1,j+1) &
                          -h0etax(i+1,j-1) + h0etax(i-1,j-1)) / dxdy4)
                
                endif  ! with cross derivatives

                matrix_rhs1_max = max(matrix_rhs1_max, abs(rhs(k_ij)))
                matrix_rhs2_max = max(matrix_rhs2_max, abs(rhs(k_ij+nD)))

              if (debug) then
                  write(67,*) '++b i,j,rhs = ',i,j,rhs(k_ij),rhs(k_ij + nD)
              endif 
                                
                  
                !if (j==3) then ! DEBUG
                if (.false.) then
                    write(67,*) '+++ i,j,rhs = ',i,j,rhs(k_ij),rhs(k_ij + nD)
                    write(67,*) h0_ij
                    write(67,*) h2_ij, ((s1(i-1,j) - 2.d0*s1(i,j) + s1(i+1,j)) / dxdx)
                    write(67,*) h3_ij, ((s1(i-1,j)/h0_imj - 2.d0*s1(i,j)/h0_ij &
                            + s1(i+1,j)/h0_ipj) / dxdx )
                    write(67,*) Bparam * grav * h0_ij**2 &
                        * ((h0etax(i-1,j) - 2.d0*h0etax(i,j) + h0etax(i+1,j)) / dxdx)
                    write(67,*) Bparam, grav, h0_ij**2, &
                      ((h0etax(i-1,j) - 2.d0*h0etax(i,j) + h0etax(i+1,j)) / dxdx)
                    write(67,*) ((h0etax(i-1,j) - 2.d0*h0etax(i,j) + h0etax(i+1,j)) / dxdx)
                    write(67,*) s1(i,j)
                    write(67,*) h0etax(i,j)
                    write(67,*) h_ij + b_ij
                endif

            enddo
        enddo
        
    deallocate(s1, s2, h0etax, h0etay)
    
    if (debug) then
        write(6,*) '+++ levelBouss = ',levelBouss,'  nng = ',nng
        write(6,*) '+++ matrix_rhs_max = ',matrix_rhs1_max, matrix_rhs2_max
    endif
    
    enddo  nng_loop

        
    minfo%matrix_nelt = nelt
    
    !============ Debug:
            
    if (.false.) then
        ! check row sums of matrix
        allocate(rowsum(2*nD))
        rowsum = 0.d0
        do k=1,nelt
            i = minfo%matrix_ia(k)
            if (minfo%matrix_ja(k) <= nD) then
                ! sums over left half of matrix:
                rowsum(i) = rowsum(i) + minfo%matrix_sa(k)
            endif
            if ((minfo%matrix_ia(k) <= nD) .and. (minfo%matrix_ja(k) > nD)) then
                write(6,*) '+++ k,ia,ja,sa: ',k,minfo%matrix_ia(k),minfo%matrix_ja(k),minfo%matrix_sa(k)
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
            
        call MATVEC(2*nD, uv, Auv, nelt, minfo%matrix_ia, minfo%matrix_ja, minfo%matrix_sa, 0)

        do k=1,nD
            write(6,*) '+++ k,Auv(k),Auv(k+nD): ',k,Auv(k),Auv(k+nD)
            enddo
                
        write(6,*) '+++ Set up matrix of dimension ',2*nD
        write(6,*) '+++ Number of nonzeros: ', nelt
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
            
            call MATVEC(2*nD, uv, Auv, nelt, minfo%matrix_ia, &
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
    
    
    if (nelt > size(minfo%matrix_sa)) then
        write(6,*) '*** Error nelt > size(minfo%matrix_sa), need to allocate larger'
        stop
    endif

    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    timePrepBuild = timePrepBuild + clock_finish-clock_start
    timePrepBuildCPU = timePrepBuildCPU + cpu_finish-cpu_start

#ifdef WHERE_AM_I
  write(*,*) "ending   buildSparseMatrixMScoo"
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

end subroutine buildSparseMatrixMScoo


