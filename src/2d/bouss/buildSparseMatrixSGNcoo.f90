
subroutine buildSparseMatrixSGNcoo(q,qold,aux,soln,rhs,                    &
                                     minfo_matrix_ia,minfo_matrix_ja,minfo_matrix_sa, &
                                     numBoussCells,levelBouss,                        &
                                     mptr,nx,ny,nvar,naux,                    &
                                     matrix_rhs1_max,matrix_rhs2_max)

    
    use amr_module
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: sea_level, dry_tolerance
    use bouss_module
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    
    !! pass in numBoussCells for this level so can dimension these arrays
    real(kind=8) :: soln(0:2*numBoussCells), rhs(0:2*numBoussCells) 
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: i,j,nng,nb,nelt,levSt,nx,ny,loc,locaux,mptr
    integer :: mitot, mjtot, neltSt, neltBegin
    integer :: jst, jend, ivar
    
    real(kind=8) :: q(nvar,1-nghost:nx+nghost,1-nghost:ny+nghost)
    real(kind=8) :: qold(nvar,1-nghost:nx+nghost,1-nghost:ny+nghost)
    real(kind=8) :: aux(naux,1-nghost:nx+nghost,1-nghost:ny+nghost)

    real(kind=8), dimension(0:nx+1,0:ny+1) :: Bxx, Byy, Bx, By, Bxy
    real(kind=8), dimension(0:nx+1,0:ny+1) :: phi, eta, w
    real(kind=8) :: etax, etay, hx,hy, wx,wy
    real(kind=8) :: xlow, ylow, x, y
    integer :: ii,jj,kd

    integer :: nD
    integer :: k, k_ij, k_imj, k_ipj, k_ijm, k_ijp 
    integer :: m,n
    integer :: k_imjp, k_imjm, k_ipjm, k_ipjp
    integer :: kh_ij, kh_imj, kh_ipj, kh_ijm, kh_ijp
    integer :: kh_imjp, kh_imjm, kh_ipjm, kh_ipjp, khu_ij, khv_ij

    real(kind=8) :: h_ij
    real(kind=8) :: dx, dy, dxm, dym, dxdx, dydy, dxdy4
    
    integer :: nelt_ij, nelt_ij_nD
    real(kind=8) :: sa, huc, hvc, matrix_rhs1_max, matrix_rhs2_max
    integer :: ja, bc, c
    real(kind=8) :: b_min, b_max
    integer :: minfo_matrix_ia(24*numBoussCells)
    integer :: minfo_matrix_ja(24*numBoussCells)
    real(kind=8) :: minfo_matrix_sa(24*numBoussCells)
    real(kind=8) :: u(1-nghost:nx+nghost,1-nghost:ny+nghost), v(1-nghost:nx+nghost,1-nghost:ny+nghost) 
    real(kind=8) :: phix, phiy

    logical :: debug
    logical :: yhi_db, ylo_db, xlo_db, xhi_db
    logical :: allzero

    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.
    

#ifdef WHERE_AM_I
  write(*,*) "starting buildSparseMatrixSGNcoo"
#endif
    
    debug = .false.
    !debug = .true.
    
    minfo => matrix_info_allLevs(levelBouss)

    !write(31,*)"Start of buildSparseMatrixSGNcoo for grid ",mptr
    !write(*,*)"Start of buildSparseMatrixSGNcoo for grid ",mptr
    mitot = nx + 2*nghost
    mjtot = ny + 2*nghost
    xlow = rnode(cornxlo,mptr)
    ylow = rnode(cornylo,mptr)
    !call symcheck2(q,mitot,mjtot)


!
! ================== Step 2  set up for matrix elements  ======================
    
    ! Set up matrix A and RHS by looping over all Bouss patches:

       if (debug) then
         write(49,*)"Level ",levelBouss," grid ",mptr," at start of buildSparseMatrixSGNcoo"
         do j = 1-nghost, ny+nghost
         do i = 1-nghost, nx+nghost
           write(49,900)i+nghost,j+nghost,(q(m,i,j),m=1,nvar)
 900       format(2i5,5e15.7)
         end do
         end do
       endif
    
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
    nelt = minfo%matrix_nelt  ! some number of matrix elements already added, 
    neltBegin = nelt  ! this is where this grid started
              
        
        ! set booleans indicating whether each patch edge is a Domain Boundary:
        !call setDomainBndries(mptr,yhi_db,ylo_db,xlo_db,xhi_db)
        xlo_db = (node(ndilo,mptr) .eq. 0)
        xhi_db = (node(ndihi,mptr) .eq. (iregsz(levelBouss)-1))
        ylo_db = (node(ndjlo,mptr) .eq. 0)
        yhi_db = (node(ndjhi,mptr) .eq. (jregsz(levelBouss)-1))
    

        
        ! Get integer array providing matrix indices mapping 
        ! cell coordinates i,j to location k in vector of all unknowns: 
           
        nb = node(bpatchNum,mptr)
        mi => minfo%matrix_indices(nb)  ! integer array for this patch

        ! temp arrays needed for this patch:
        ! provide one layer of ghost cells needed when computing
        ! second differences of these quantities, but values in ghost
        ! cells are never actually needed since second diff's won't be
        ! needed in first row of interior cells
        
        ! First loop over all cells (i,j), one row of ghost cells,
        ! computing s1, s2, h0etax, h0etay,
        ! quantities that involve first derivatives.  In the second loop
        ! over all cells below we will compute matrix_rhs using second
        ! differences of these quantities.
        
        ! compute eta used below in update
        ! Remember that ETA not same size as Q
         eta(:,:) = q(1,0:nx+1,0:ny+1) +   &
                  aux(1,0:nx+1,0:ny+1)

         ! will have to FIX for h = 0? Hopefully not that close to 0
         !u(:,:) = q(2,:,:) / q(1,:,:)
         !v(:,:) = q(3,:,:) / q(1,:,:)
         !where (q(1,:,:) /= 0.d0)
         where (abs(q(1,:,:)) > dry_tolerance)
           u(:,:) = q(2,:,:) / q(1,:,:)
           v(:,:) = q(3,:,:) / q(1,:,:)
         elsewhere
           u(:,:) = 0.d0
           v(:,:) = 0.d0
         end where


         do j=0,ny+1 
            if (coordinate_system .eq. 2) then
              ! aux(3,i,j) is latitude, doesn't depend on i so use i=1
              dxm = dx * earth_radius * deg2rad * cos(aux(3,1,j))
            else
              dxm = dx
            endif
            dxdx = dxm**2
            dxdy4 = 4.d0*dxm*dym
            
            ! compute phi. rest will now be computed ! in second loop
            do i=0,nx+1
              
               phi(i,j) = ((v(i+1,j)-v(i-1,j))*(u(i,j+1)-u(i,j-1)) -      &
                          (u(i+1,j)-u(i-1,j))*(v(i,j+1)-v(i,j-1)))/dxdy4  &
                          + ((u(i+1,j)-u(i-1,j))/(2.d0*dxm)               &
                          +  (v(i,j+1)-v(i,j-1))/(2.d0*dym))**2
               Bx(i,j) = (aux(1,i+1,j) - aux(1,i-1,j))/(2.d0*dxm) 
               By(i,j) = (aux(1,i,j+1) - aux(1,i,j-1))/(2.d0*dym)
               Bxy(i,j) = (aux(1,i+1,j+1)-aux(1,i+1,j-1) &
                          -aux(1,i-1,j+1)+aux(1,i-1,j-1))/dxdy4
               Bxx(i,j) = (aux(1,i+1,j)-2.d0*aux(1,i,j)+aux(1,i-1,j))/dxm**2
               Byy(i,j) = (aux(1,i,j+1)-2.d0*aux(1,i,j)+aux(1,i,j-1))/dym**2                          
               w(i,j) = u(i,j)**2 * Bxx(i,j) + 2.d0*u(i,j)*v(i,j)*Bxy(i,j) &
                       +v(i,j)**2 * Byy(i,j)
              
            enddo
        enddo

! ================== build matrix and rhs
        ! Second loop over all cells.  ONLY INTERIOR OK NOW
        
        do j=1,ny 
            if (coordinate_system .eq. 2) then
              ! aux(3,i,j) is latitude, doesn't depend on i so use i=1
              dxm = dx * earth_radius * deg2rad *              &
                    cos(aux(3,1,j))
            else
              dxm = dx
            endif
            dxdx = dxm**2
            dxdy4 = 4.d0*dxm*dym
            
            
            do i=1,nx

                ! sparse matrix row indices for cell (i,j) are:
                !    k_ij for the row of I-D11 and -D12 in block matrix,
                !    k_ij + nD for the row of -D21 and I-D22
                
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
                
                ! use the fact that h,hu,hv are contiguous in alloc:
                h_ij = q(1,i,j)  ! to save typing keep this one
                
                k_ij = mi%mindex(i,j)
                if (k_ij == -1) then
                    write(6,*) '+++ unexpected k_ij=-1 at i,j: ',i,j
                endif
                neltSt = nelt
                if (debug) then
                  write(44,999)i,j,k_ij,k_ij,nelt
 999              format("cell ",2i4," k_ij = ",i5," before shift was ",i5,"   starting with nelt ",i6)
                endif
                
                !b_min = min(aux(1,i,j), aux(1,i-1,j),aux(1,i+1,j),     &
                !            aux(1,i,j-1),aux(1,i,j+1),                 &
                !            aux(1,i-1,j-1),aux(1,i-1,j+1),             &
                !            aux(1,i+1,j-1),aux(1,i+1,j+1))
                !b_max = max(aux(1,i,j), aux(1,i-1,j),aux(1,i+1,j),     &
                !            aux(1,i,j-1),aux(1,i,j+1),                 &
                !            aux(1,i-1,j-1),aux(1,i-1,j+1),             &
                !            aux(1,i+1,j-1),aux(1,i+1,j+1))
                if (.not. mi%isBouss(i,j)) then
                        
                    ! The deepest water in block of 9 cells is shallower 
                    ! than boussMinDepth, or else one of the cells is onshore.

                    ! Revert to SWE in this cell, i.e. two rows
                    ! of matrix should be rows of identity with rhs = 0,
                    ! so correction terms huc,hvc calculated will be 0.
                    
                    !QUESTION: I should be h_ij*I. FIXED
                    ! BUT: then discovered the leading h in matrix is a typo
                    ! fixed for now by multiplying rhs by h_ij too.
                    
                    ! Replace diagonal block I-D11 by I matrix:
                    nelt = nelt + 1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = k_ij

                    ! changed to +1*max on sept 2 2022 to test gamg
                    minfo_matrix_sa(nelt) =  1.d0*max(abs(h_ij),1000.d0)
                    
                    ! Replace diagonal block I-D22 by I matrix:
                    nelt = nelt + 1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = k_ij + nD

                    ! changed to +1*max on sept 2 2022 to test gamg
                    minfo_matrix_sa(nelt) =  1.d0*max(abs(h_ij),1000.d0)
                    
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
                nelt_ij_nD = -1
                
                hx = (q(1,i+1,j) - q(1,i-1,j))/(2.d0*dxm)
                hy = (q(1,i,j+1) - q(1,i,j-1))/(2.d0*dym)
                etax = (eta(i+1,j) - eta(i-1,j))/(2.d0*dxm)
                etay = (eta(i,j+1) - eta(i,j-1))/(2.d0*dym)

                !-----------------------------------------------
                ! I-D11 block (x and xx derivatives of hu components):
                
                ! Matrix element row U(i,j) column U(i,j):
                sa = (1.d0 + alpha*(2.d0*h_ij**2/(3.d0*dxdx) &
                     + 0.5d0*h_ij*Bxx(i,j) + Bx(i,j)*etax))
                nelt = nelt+1
                minfo_matrix_ia(nelt) = k_ij
                minfo_matrix_ja(nelt) = k_ij
                minfo_matrix_sa(nelt) = sa
                
                nelt_ij = nelt  ! index of diagonal matrix element
                
                
                ! Matrix element row U(i,j) column U(i+1,j):
                sa = alpha*(-h_ij*hx/(2.d0*dxm) - h_ij**2/(3.d0*dxdx))
                
                if (k_ipj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = k_ipj
                    minfo_matrix_sa(nelt) = sa
                else
                    ! right boundary of patch
                    if (.not. xhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = q(4,i+1,j)
                        rhs(k_ij) = rhs(k_ij) - sa*huc
                    else if (bc_xhi==1) then 
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij) = minfo_matrix_sa(nelt_ij) + sa
                    else if (bc_xhi==3) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij) = minfo_matrix_sa(nelt_ij) - sa
                    endif                    
                endif
                
                
                ! Matrix element row U(i,j) column U(i-1,j):
                sa = alpha*(h_ij*hx/(2.d0*dxm) - h_ij**2/(3.d0*dxdx))
                if (k_imj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = k_imj
                    minfo_matrix_sa(nelt) = sa
                else
                    ! left boundary of patch
                    if (.not. xlo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = q(4,i-1,j)
                        rhs(k_ij) = rhs(k_ij) - sa*huc
                    else if (bc_xlo==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij) = minfo_matrix_sa(nelt_ij) + sa
                    else if (bc_xlo==3) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij) = minfo_matrix_sa(nelt_ij) - sa
                    endif
                endif


                !-----------------------------------------------
                ! I-D22 block (yy derivatives of hv components):
                
                ! Matrix element row V(i,j) column V(i,j):
                sa = (1.d0 + alpha*(2.d0*h_ij**2/(3.d0*dydy) &
                     + 0.5d0*h_ij*Byy(i,j) + By(i,j)*etay))
                nelt = nelt+1
                minfo_matrix_ia(nelt) = k_ij + nD
                minfo_matrix_ja(nelt) = k_ij + nD
                minfo_matrix_sa(nelt) = sa
                
                nelt_ij_nD = nelt  ! index of diagonal matrix element


                ! Matrix element row V(i,j) column V(i,j+1):
                sa = alpha*(-h_ij*hy/(2.d0*dym) - h_ij**2/(3.d0*dydy))
                if (k_ijp > -1) then
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = k_ijp + nD
                    minfo_matrix_sa(nelt) = sa
                else
                    ! top boundary of patch
                    if (.not. yhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = q(5,i,j+1)
                        rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                    else if (bc_yhi==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij_nD) = &
                                minfo_matrix_sa(nelt_ij_nD) + sa
                    else if (bc_yhi==3) then
                        ! For wall BC add -sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij_nD) = &
                                minfo_matrix_sa(nelt_ij_nD) - sa
                    endif
                endif
                
                ! Matrix element row V(i,j) column V(i,j-1):
                sa = alpha*(h_ij*hy/(2.d0*dym) - h_ij**2/(3.d0*dydy))
                if (k_ijm > -1) then
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = k_ijm + nD
                    minfo_matrix_sa(nelt) = sa
                else
                    ! bottom boundary of patch
                    if (.not. ylo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = q(5,i,j-1)
                        rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                    else if (bc_ylo==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij_nD) = &
                                minfo_matrix_sa(nelt_ij_nD) + sa
                    else if (bc_ylo==3) then
                        ! For wall BC add -sa to diagonal matrix element:
                        minfo_matrix_sa(nelt_ij_nD) = &
                                minfo_matrix_sa(nelt_ij_nD) - sa
                    endif   
                endif                                        

            
                !-----------------------------------------------
                ! -D12 block (cross derivatives):
                ! row ia = k_ij for hu corrections in cell (i,j)
                ! columns k_ipjp + nD and other neighbors for hv corrections                

                ! Matrix element row U(i,j) column V(i,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH VARYING TOPO
                sa = alpha*(0.5d0*h_ij*Bxy(i,j) + etax*By(i,j))
                nelt = nelt+1
                minfo_matrix_ia(nelt) = k_ij
                minfo_matrix_ja(nelt) = k_ij + nD
                minfo_matrix_sa(nelt) = sa
                    
            
                ! Matrix element row U(i,j) column V(i,j+1):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = -alpha*h_ij*(hx + 0.5d0*Bx(i,j))/(2.d0*dym)
                ja = -1
                if (k_ijp > -1) then
                    ja = k_ijp + nD
                else if (.not. yhi_db) then
                    hvc = q(5,i,j+1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else
                    ! at top boundary, use cell ij instead,
                    ja = k_ij + nD
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i,j-1):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = alpha*(h_ij*(hx + 0.5d0*Bx(i,j))/(2.d0*dym))
                ja = -1
                if (k_ijm > -1) then
                    ja = k_ijm + nD
                else if (.not. ylo_db) then
                    hvc = q(5,i,j-1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else
                    ! at bottom boundary, use cell ij instead,
                    ja = k_ij + nD
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif

                
                ! Matrix element row U(i,j) column V(i+1,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = alpha*0.5d0*h_ij*By(i,j) / (2.d0*dxm)
                ja = -1
                if (k_ipj > -1) then
                    ja = k_ipj + nD
                else if (.not. xhi_db) then
                    hvc = q(5,i+1,j) 
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else
                    ! at right boundary, use cell ij instead,
                    ja = k_ij + nD
                    ! no reflection of V at right wall
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif

                
                ! Matrix element row U(i,j) column V(i-1,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = -alpha*0.5d0*h_ij*By(i,j) / (2.d0*dxm)
                ja = -1
                if (k_imj > -1) then
                    ja = k_imj + nD
                else if (.not. xlo_db) then
                    hvc = q(5,i-1,j)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else
                    ! at left boundary, use cell ij instead,
                    ja = k_ij + nD
                    ! no reflection of V at left wall
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                                                
                
                ! Matrix element row U(i,j) column V(i+1,j+1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjp > -1) then
                    ja = k_ipjp + nD
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                 else if ((.not. yhi_db .or. j .lt. ny) .and.    &
                          (.not. xhi_db .or. i .lt. nx)) then
                    hvc = q(5,i+1,j+1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else if (k_ijp > -1) then
                    ! at right boundary, use cell above instead
                    ja = k_ijp + nD
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead, 
                    ja = k_ipj + nD
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                else if (xhi_db .and. i .eq. nx .and. .not. yhi_db) then
                    ja = k_ij + nD
                else
                    ! top-right corner
                    ja = k_ij + nD
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                    
                ! Matrix element row U(i,j) column V(i+1,j-1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjm > -1) then
                    ja = k_ipjm + nD            
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                else if ((.not. ylo_db .or. j>1) .and.            &
                         (.not. xhi_db .or. i .lt. nx)) then
                    hvc = q(5,i+1,j-1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc   
                else if (k_ijm > -1) then
                    ! at right boundary, use cell below instead
                    ja = k_ijm + nD
                else if (k_ipj > -1) then
                    ! at bottom boundary, use cell to right instead,
                    ja = k_ipj + nD
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                else if (xhi_db .and. i.eq.nx .and. .not. ylo_db) then
                    ja = k_ij + nD
                else
                    ! bottom-right corner
                    ja = k_ij + nD
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif

                
                ! Matrix element row U(i,j) column V(i-1,j+1):    
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjp > -1) then
                    ja = k_imjp + nD
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                else if ((.not. yhi_db .or. j .lt. ny) .and. &
                         (.not. xlo_db .or. i > 1)) then
                    hvc = q(5,i-1,j+1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else if (k_ijp > -1) then
                    ! at left domain boundary, use cell above instead
                    ja = k_ijp + nD
                else if (k_imj > -1) then
                    ! at top domain boundary, use cell to left instead,
                    ja = k_imj + nD
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                else if (xlo_db .and. i .eq. 1 .and. .not. yhi_db) then
                    ja = k_ij + nD
                else
                    ! top-left corner
                    ja = k_ij + nD
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i-1,j-1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjm > -1) then
                    ja = k_imjm + nD
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                else if ((.not. ylo_db .or. j>1) .and.      &
                         (.not. xlo_db .or. i>1)) then
                    hvc = q(5,i-1,j-1)
                    rhs(k_ij) = rhs(k_ij) - sa*hvc
                else if (k_ijm > -1) then
                    ! at left boundary, use cell below instead
                   ja = k_ijm + nD
                else if (k_imj > -1) then
                    ! at bottom boundary, use cell to left instead,
                    ja = k_imj + nD
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                else if (xlo_db .and. i .eq. 1 .and. .not. ylo_db) then
                   ja = k_ij + nD
                else
                    ! bottom-left corner
                    ja = k_ij + nD
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif


                !-----------------------------------------------
                ! -D21 block (cross derivatives
                ! row ia = k_ij + nD  for hv corrections in cell (i,j)
                ! columns k_ipjp and other neighbors for hu corrections

                ! Matrix element row V(i,j) column U(i,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH VARYING TOPO
                sa = alpha*(0.5d0*h_ij*Bxy(i,j) + etay*Bx(i,j))
                nelt = nelt+1
                minfo_matrix_ia(nelt) = k_ij + nD
                minfo_matrix_ja(nelt) = k_ij
                minfo_matrix_sa(nelt) = sa
                
                ! Matrix element row V(i,j) column U(i+1,j):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = -alpha*h_ij*(hy + 0.5d0*By(i,j))/(2.d0*dxm)
                ja = -1
                if (k_ipj > -1) then
                    ja = k_ipj
                else if (.not. xhi_db) then
                    huc = q(4,i+1,j)
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                else
                    ! at right boundary, use cell ij instead,
                    ja = k_ij
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i-1,j):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = alpha*h_ij*(hy + 0.5d0*By(i,j))/(2.d0*dxm)
                ja = -1
                if (k_imj > -1) then
                    ja = k_imj
                else if (.not. xlo_db) then
                    huc = q(4,i-1,j)
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                else
                    ! at left boundary, use cell ij instead,
                    ja = k_ij
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i,j+1):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = alpha*0.5d0*h_ij*Bx(i,j) / (2.d0*dym)
                ja = -1
                if (k_ijp > -1) then
                    ja = k_ijp ! FIXED BUG2
                else if (.not. yhi_db) then
                    huc = q(4,i,j+1)
                    rhs(k_ij) = rhs(k_ij) - sa*huc
                else
                    ! at top boundary, use cell ij instead,
                    ja = k_ij
                    ! no reflection of U at top wall
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif

                ! Matrix element row V(i,j) column U(i,j-1): ! FIXED BUG2
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = -alpha*0.5d0*h_ij*Bx(i,j) / (2.d0*dym)
                ja = -1
                if (k_ijm > -1) then
                    ja = k_ijm ! FIXED BUG2
                else if (.not. ylo_db) then
                    huc = q(4,i,j-1)
                    rhs(k_ij) = rhs(k_ij) - sa*huc
                else
                    ! at bottom boundary, use cell ij instead,
                    ja = k_ij ! FIXED BUG2
                    ! no reflection of U at bottom wall
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif

                
                ! Matrix element row V(i,j) column U(i+1,j+1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjp > -1) then
                    ja = k_ipjp
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                else if ((.not. yhi_db .or. j.lt.ny) .and. &
                         (.not. xhi_db .or. i.lt.nx)) then
                    huc = q(4,i+1,j+1)
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                else if (k_ijp > -1) then
                    ! at right boundary, use cell above instead,
                    ja = k_ijp
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead,
                    ja = k_ipj
                else if (yhi_db .and. j .eq. ny .and. .not. xhi_db) then
                    ja = k_ij
                else
                    ! top-right corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i+1,j-1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjm > -1) then
                    ja = k_ipjm
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                else if ((.not. ylo_db .or. j>1) .and. &
                         (.not. xhi_db .or. i .lt. nx)) then
                    huc = q(4,i+1,j-1)
                    rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                else if (k_ijm > -1) then
                    ! at right boundary, use cell below instead,
                    ja = k_ijm
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                else if (k_ipj > -1) then
                    ! at bottom boundary, use cell to right instead,
                    ja = k_ipj
                else if (ylo_db .and. j .eq. 1 .and. .not. xhi_db) then
                    ja = k_ij
                else
                    ! bottom-right corner
                    ja = k_ij ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                    if (debug) write(67,*) i,j,i+1,j-1,c,nelt,sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i-1,j+1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjp > -1) then
                    ja = k_imjp
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                else if ((.not. yhi_db .or. j .lt. ny) .and. &
                         (.not. xlo_db .or. i>1)) then
                    huc = q(4,i-1,j+1)
                    rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                else if (k_ijp > -1) then
                    ! at left boundary, use cell above instead,
                    ja = k_ijp
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                else if (k_imj > -1) then
                    ! at top boundary, use cell to left instead,
                    ja = k_imj
                else if (yhi_db .and. j.eq.ny .and. .not. xlo_db) then
                    ja = k_ij
                else
                    ! top-left corner
                    ja = k_ij !FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif


                ! Matrix element row V(i,j) column U(i-1,j-1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)  
                ja = -1
                if (k_imjm > -1) then
                    ja = k_imjm
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                else if ((.not. ylo_db .or. j>1) .and. &
                         (.not. xlo_db .or. i>1)) then
                    huc = q(4,i-1,j-1)
                    rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                else if (k_ijm > -1) then
                    ! at left boundary, use cell below instead,
                    ja = k_ijm
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                else if (k_imj > -1) then
                    ! at bottom boundary, use cell to left instead,
                    ja = k_imj
                else if (ylo_db .and. j .eq. 1 .and. .not. xlo_db) then
                    ja = k_ij
                else
                    ! bottom-left corner
                    ja = k_ij !FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                endif
                
                if (ja > -1) then
                    ! add a new matrix element for interior, wall, Neumann case
                    nelt = nelt+1
                    minfo_matrix_ia(nelt) = k_ij + nD
                    minfo_matrix_ja(nelt) = ja
                    minfo_matrix_sa(nelt) = sa
                endif
                
 
                !-------------------------------------------------------
                ! Right hand sides for vector elements k_ij and k_ij+nD
                ! Need to add comments on SGN RHS.
                
              if (debug) then
                  write(67,*) '++before  i,j,rhs = ',i,j,rhs(k_ij),rhs(k_ij + nD)
              endif 
              
              !RHS for SGN 
              phix = (phi(i+1,j) - phi(i-1,j))/(2.d0*dxm)
              phiy = (phi(i,j+1) - phi(i,j-1))/(2.d0*dym)
                            
              wx = (w(i+1,j) - w(i-1,j))/(2.d0*dxm)
              wy = (w(i,j+1) - w(i,j-1))/(2.d0*dym)
              
              etax = (eta(i+1,j) - eta(i-1,j))/(2.d0*dxm)
              etay = (eta(i,j+1) - eta(i,j-1))/(2.d0*dym)
                                 
              rhs(k_ij) = rhs(k_ij)  +  (grav/alpha * etax      &
                    + 2.d0*h_ij*(h_ij/3.d0 * phix &
                    + phi(i,j) * (hx + 0.5d0*Bx(i,j)))   &
                    + 0.5d0*h_ij*wx + w(i,j)*etax)

              rhs(k_ij+nD) = rhs(k_ij+nD) +  (grav/alpha * etay &
                    + 2.d0*h_ij*(h_ij/3.d0 * phiy &
                    + phi(i,j) * (hy + 0.5d0*By(i,j)))   &
                    + 0.5d0*h_ij*wy + w(i,j)*etay)


              if (debug) then
                  write(67,*) '++after   i,j,rhs = ',i,j,rhs(k_ij),rhs(k_ij + nD)
              endif 

              if (debug) then
                write(44,998) nelt,nelt-neltSt
 998                format("                               ending with nelt ",i6,"   used ",i5)
              endif
                     

        if (debug) then ! dump matrix to look for singularity 
           write(88,*)" triplets for i,j ",i,j," level  ",levelBouss
           do k = neltBegin+1, nelt
              !write(88,103) minfo_matrix_ia(k),minfo_matrix_ja(k),minfo_matrix_sa(k)
              write(88,103) minfo_matrix_sa(k)
 !103          format(2i8,e16.7)
 103          format(e16.7)
           end do
           write(88,*)

           !close(88)
           !write(89,*)" level ",levelBouss
           !do k = 1,2*numBoussCells
           !   write(89,104) k,rhs(k)
 104          format(i5,e16.7)
           !end do
           !!close(89)
           !!stop
        endif

            enddo
        enddo

        if (debug) then
         if (mptr .eq. 1) then
           write(88,*)" grid ",mptr," level ",levelBouss,"x,y,ia,ja,sa"
           !do ii = 0, nx+1
           !do jj = 0, ny+1
           do ii = 1, nx
           do jj = 1, ny
                k = mi%mindex(ii,jj)
                kd = k + nD
                x = xlow + (ii-.5d0)*dx
                y = ylow + (jj-.5d0)*dy
                write(88,105) x,y,k,minfo_matrix_ia(k),minfo_matrix_ja(k),minfo_matrix_sa(k)
                write(88,105) x,y,k+nD,minfo_matrix_ia(k+nD),minfo_matrix_ja(k+nD),minfo_matrix_sa(k+nD)
 105          format(2e15.7,3i8,e16.7)
 799          format(2i5,4e15.7)
 798          format(6e15.7,3i5)
           end do
           end do
           endif
        endif

        
    
    ! put back in struct. May have more grids following    
    minfo%matrix_nelt = nelt
    
#ifdef WHERE_AM_I
  write(*,*) "ending   buildSparseMatrixSGNcoo"
#endif


    return
end subroutine buildSparseMatrixSGNcoo
