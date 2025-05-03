
subroutine buildSparseMatrixSGNcrs(q,qold,aux,soln,rhs,rowPtr,cols,vals,   &
                                     numBoussCells,levelBouss,                        &
                                     mptr,nx,ny,nvar,naux)

    
    use amr_module
    use geoclaw_module, only: earth_radius, deg2rad, coordinate_system, grav
    use geoclaw_module, only: sea_level, dry_tolerance
    use bouss_module
        
    implicit none
    
    integer, intent(in) :: nvar, naux, levelBouss, numBoussCells
    integer, intent(inout) :: rowPtr(0:2*numBoussCells)
    real(kind=8), intent(inout) :: vals(0:24*numBoussCells)
    integer, intent(inout) :: cols(0:24*numBoussCells)
    
    !! pass in numBoussCells for this level so can dimension these arrays
    real(kind=8) :: soln(0:2*numBoussCells), rhs(0:2*numBoussCells) 
    
    type(matrix_patchIndex), pointer :: mi
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: kc,i,j,nng,nb,levSt,nx,ny,loc,locaux,mptr
    integer :: mitot, mjtot
    integer :: jst, jend, ivar
    integer :: colPtr
    integer :: ii,jj
    
    real(kind=8) :: q(nvar,1-nghost:nx+nghost,1-nghost:ny+nghost)
    real(kind=8) :: qold(nvar,1-nghost:nx+nghost,1-nghost:ny+nghost)
    real(kind=8) :: aux(naux,1-nghost:nx+nghost,1-nghost:ny+nghost)

    real(kind=8) ::etaxM(0:nx+1,0:ny+1)
    real(kind=8) ::etayM(0:nx+1,0:ny+1)
    real(kind=8) :: bu, bv,x, xlow, y, ylow
    real(kind=8) ::brhs(2,1-nghost:nx+nghost,1-nghost:ny+nghost)

    real(kind=8), dimension(0:nx+1,0:ny+1) :: Bxx, Byy, Bx, By, Bxy
    real(kind=8), dimension(0:nx+1,0:ny+1) :: phi, eta, w
    real(kind=8) :: etax, etay, hx,hy, wx,wy

    integer :: nD
    integer :: k, k_ij, k_imj, k_ipj, k_ijm, k_ijp 
    integer :: m,n
    integer :: k_imjp, k_imjm, k_ipjm, k_ipjp
    integer :: kh_ij, kh_imj, kh_ipj, kh_ijm, kh_ijp
    integer :: kh_imjp, kh_imjm, kh_ipjm, kh_ipjp, khu_ij, khv_ij

    real(kind=8) :: h_ij
    real(kind=8) :: dx, dy, dxm, dym, dxdx, dydy, dxdy4
    
    real(kind=8) :: sa, huc, hvc
    integer :: ja, bc, c, indexSt
    real(kind=8) :: b_min, b_max
    real(kind=8) :: u(1-nghost:nx+nghost,1-nghost:ny+nghost), v(1-nghost:nx+nghost,1-nghost:ny+nghost) 
    real(kind=8) :: phix, phiy

    logical :: debug
    logical :: yhi_db, ylo_db, xlo_db, xhi_db
    logical :: allzero, IS_GHOST

    
    ! bc_xlo, bc_ylo, bc_xhi, bc_yhi are set in bouss_module.
    !     0 means Dirchlet value 0 (no correction at each of union of patches)
    !     1 means Neumann (correction in first interior cell = in ghost cell)
    !     2 means use ghost cell values in q(4:5,:,:) as Dirichlet BCs
    !             Multiply by matrix element and move to RHS
    !             Requires nvar==5.

    ! for CRS format, order of unknowns for each row is 
    ! 12 for u equations, 12 for v equation
    ! the unknown are still order in 0 based order: cell k_ij is u eq<> 2*(kij-1) v<>2*k_ij-1
    ! to get columns to be sorted, indices assigned as follows interior
    ! cell puts columns in order that comes out right. This doesnot hold
    ! for ghost cells, so those 12 entires need to be sorted.
    ! The order for u starting from rowPtr index:
    !   0         1          2         3        4        5     6       7 
    ! v_i-1,j-1  v_i,j-1  v_i+1,j-1  u_i-1,j  v_i-1,j  u_ij  v_ij   u_i+1,j
    !  8          9         10          11
    ! v_i+1,j  v_i-1,j+1, v_i,j+1, v_i+1,j+1
    ! The order for v starting from rowPtr+12 index:
    !   12        13         14        15      16        17    18      19
    ! u_i-1,j-1  u_i,j-1  v_i,j-1  u_i+1,j-1  u_i-1,j  u_ij   v_ij   u_i+1,j
    !  20         21        22          23
    ! u_i-1,j+1  u_i,j+1, v_i,j+1, u_i+1,j+1


   IS_GHOST(i,j) = (i.eq.1 .or. i.eq.nx .or. j.eq.1 .or. j.eq.ny) 

#ifdef WHERE_AM_I
  write(*,*) "starting buildSparseMatrixSGNcrs"
#endif
    
    debug = .false.
    !debug = .true.
    etaxM = 0.d0
    etayM = 0.d0
    
    minfo => matrix_info_allLevs(levelBouss)

    mitot = nx + 2*nghost
    mjtot = ny + 2*nghost
    xlow = rnode(cornxlo,mptr)
    ylow = rnode(cornylo,mptr)


!
! ================== Step 2  set up for matrix elements  ======================
    
    ! Set up matrix A and RHS by looping over all Bouss patches:

       if (debug) then
         write(49,*)"Level ",levelBouss," grid ",mptr," at start of cleanBuild"
         do i = 1-nghost, nx+nghost
         do j = 1-nghost, ny+nghost
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
              
    !REORG to put u and v together in blocks

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
                ! for CRS set row using 24 entries per cell i,j: 12 for u and 12 for v.
                rowPtr(2*(k_ij-1)) = 12*2*(k_ij-1)        ! start of u_ij eq
                rowPtr(2*k_ij-1) = rowPtr(2*(k_ij-1))+12  ! start of v_ij eq
                colPtr = rowPtr(2*(k_ij-1))

                ! initialize all cols for these equations to start, even if it's col*unknown = 0
                
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
                    cols(colPtr+5) = 2*(k_ij-1)

                    ! changed to +1*max on sept 2 2022 to test gamg
                    vals(colPtr+5) = vals(colPtr+5)+1.d0*max(abs(h_ij),1000.d0)
                    
                    ! Replace diagonal block I-D22 by I matrix:
                    cols(colPtr+18) = 2*k_ij-1

                    ! changed to +1*max on sept 2 2022 to test gamg
                    vals(colPtr+18) = vals(colPtr+18)+1.d0*max(abs(h_ij),1000.d0)
                    
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
                
                hx = (q(1,i+1,j) - q(1,i-1,j))/(2.d0*dxm)
                hy = (q(1,i,j+1) - q(1,i,j-1))/(2.d0*dym)
                etax = (eta(i+1,j) - eta(i-1,j))/(2.d0*dxm)
                etay = (eta(i,j+1) - eta(i,j-1))/(2.d0*dym)
                etaxM(i,j) = etax
                etayM(i,j) = etay

                !-----------------------------------------------
                ! I-D11 block (x and xx derivatives of hu components):
                
                ! Matrix element row U(i,j) column U(i,j):
                sa = (1.d0 + alpha*(2.d0*h_ij**2/(3.d0*dxdx) &
                     + 0.5d0*h_ij*Bxx(i,j) + Bx(i,j)*etax))
                cols(colPtr+5) = 2*(k_ij-1)
                vals(colPtr+5) = vals(colPtr+5)+sa
                
                
                
                ! Matrix element row U(i,j) column U(i+1,j):
                sa = alpha*(-h_ij*hx/(2.d0*dxm) - h_ij**2/(3.d0*dxdx))
                
                if (k_ipj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    cols(colPtr+7) = 2*(k_ipj-1)
                    vals(colPtr+7) = vals(colPtr+7)+sa
                else
                    ! right boundary of patch
                    if (.not. xhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = q(4,i+1,j)
                        rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*huc
                    else if (bc_xhi==1) then 
                        ! For Neumann BC add sa to diagonal matrix element:
                        vals(colPtr+5) = vals(colPtr+5) + sa
                    else if (bc_xhi==3) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        vals(colPtr+5) = vals(colPtr+5) - sa
                    endif                    
                endif
                
                
                ! Matrix element row U(i,j) column U(i-1,j):
                sa = alpha*(h_ij*hx/(2.d0*dxm) - h_ij**2/(3.d0*dxdx))
                if (k_imj > -1) then
                    ! adjacent cell is interior pt, set matrix element to sa:
                    cols(colPtr+3) = 2*(k_imj-1)
                    vals(colPtr+3) = vals(colPtr+3)+sa
                else
                    ! left boundary of patch
                    if (.not. xlo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        huc = q(4,i-1,j)
                        rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*huc
                    else if (bc_xlo==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        vals(colPtr+5) = vals(colPtr+5) + sa
                    else if (bc_xlo==3) then
                        ! For wall, negate u: add -sa to diagonal matrix element:
                        vals(colPtr+5) = vals(colPtr+5) - sa
                    endif
                endif


                !-----------------------------------------------
                ! I-D22 block (yy derivatives of hv components):
                
                ! Matrix element row V(i,j) column V(i,j):
                sa = (1.d0 + alpha*(2.d0*h_ij**2/(3.d0*dydy) &
                     + 0.5d0*h_ij*Byy(i,j) + By(i,j)*etay))
                cols(colPtr+18) = 2*k_ij-1
                vals(colPtr+18) = vals(colPtr+18)+sa
                
                ! Matrix element row V(i,j) column V(i,j+1):
                sa = alpha*(-h_ij*hy/(2.d0*dym) - h_ij**2/(3.d0*dydy))
                if (k_ijp > -1) then
                    cols(colPtr+22) = 2*k_ijp-1
                    vals(colPtr+22) = vals(colPtr+22)+sa
                else
                    ! top boundary of patch
                    if (.not. yhi_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = q(5,i,j+1)
                        !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                        !rhs(2*k_ij) = rhs(2*k_ij) - sa*hvc
                        rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*hvc
                    else if (bc_yhi==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        vals(colPtr+18) = vals(colPtr+18) + sa
                    else if (bc_yhi==3) then
                        ! For wall BC add -sa to diagonal matrix element:
                        vals(colPtr+18) = vals(colPtr+18) - sa
                    endif
                endif
                
                ! Matrix element row V(i,j) column V(i,j-1):
                sa = alpha*(h_ij*hy/(2.d0*dym) - h_ij**2/(3.d0*dydy))
                if (k_ijm > -1) then
                    cols(colPtr+14) = 2*k_ijm-1
                    vals(colPtr+14) = vals(colPtr+14)+sa
                else
                    ! bottom boundary of patch
                    if (.not. ylo_db) then
                        ! Dirichlet BC using Bouss correction from ghost cell
                        hvc = q(5,i,j-1)
                        !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*hvc
                        !rhs(2*k_ij) = rhs(2*k_ij) - sa*hvc
                        rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*hvc
                    else if (bc_ylo==1) then
                        ! For Neumann BC add sa to diagonal matrix element:
                        vals(colPtr+18) = vals(colPtr+18) + sa
                    else if (bc_ylo==3) then
                        ! For wall BC add -sa to diagonal matrix element:
                        vals(colPtr+18) = vals(colPtr+18) - sa
                    endif   
                endif                                        

            
                !-----------------------------------------------
                ! -D12 block (cross derivatives):
                ! row ia = k_ij for hu corrections in cell (i,j)
                ! columns k_ipjp + nD and other neighbors for hv corrections                

                ! Matrix element row U(i,j) column V(i,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH VARYING TOPO
                sa = alpha*(0.5d0*h_ij*Bxy(i,j) + etax*By(i,j))
                cols(colPtr+6) = 2*k_ij-1 
                vals(colPtr+6) = vals(colPtr+6)+sa
                    
            
                ! Matrix element row U(i,j) column V(i,j+1):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = -alpha*h_ij*(hx + 0.5d0*Bx(i,j))/(2.d0*dym)
                ja = -1
                if (k_ijp > -1) then
                    !ja = k_ijp + nD
                    !ja = 2*k_ijp
                    ja = 2*k_ijp-1
                    cols(colPtr+10) = ja
                    vals(colPtr+10) = vals(colPtr+10) + sa
                else if (.not. yhi_db) then
                    hvc = q(5,i,j+1)
                    !rhs(k_ij) = rhs(k_ij) - sa*hvc
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else
                    ! at top boundary, use cell ij instead,
                    !ja = k_ij + nD
                    !ja = 2*k_ij 
                    ja = 2*k_ij-1 
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i,j-1):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = alpha*(h_ij*(hx + 0.5d0*Bx(i,j))/(2.d0*dym))
                ja = -1
                if (k_ijm > -1) then
                    ja = 2*k_ijm-1
                    cols(colPtr+1) = ja
                    vals(colPtr+1) = vals(colPtr+1) + sa
                else if (.not. ylo_db) then
                    hvc = q(5,i,j-1)
                    !rhs(k_ij) = rhs(k_ij) - sa*hvc
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else
                    ! at bottom boundary, use cell ij instead,
                    ja = 2*k_ij-1
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i+1,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = alpha*0.5d0*h_ij*By(i,j) / (2.d0*dxm)
                ja = -1
                if (k_ipj > -1) then
                    !ja = k_ipj + nD
                    !ja = 2*k_ipj
                    ja = 2*k_ipj-1
                    cols(colPtr+8) = ja
                    vals(colPtr+8) = vals(colPtr+8) + sa
                else if (.not. xhi_db) then
                    hvc = q(5,i+1,j) 
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else
                    ! at right boundary, use cell ij instead,
                    !ja = k_ij + nD
                    !ja = 2*k_ij
                    ja = 2*k_ij-1
                    ! no reflection of V at right wall
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i-1,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = -alpha*0.5d0*h_ij*By(i,j) / (2.d0*dxm)
                ja = -1
                if (k_imj > -1) then
                    !ja = k_imj + nD
                    !ja = 2*k_imj
                    ja = 2*k_imj-1
                    cols(colPtr+4) = ja
                    vals(colPtr+4) = vals(colPtr+4) + sa
                else if (.not. xlo_db) then
                    hvc = q(5,i-1,j)
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else
                    ! at left boundary, use cell ij instead,
                    !ja = k_ij + nD
                    !ja = 2*k_ij
                    ja = 2*k_ij-1
                    ! no reflection of V at left wall
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i+1,j+1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjp > -1) then
                    !ja = k_ipjp + nD
                    !ja = 2*k_ipjp 
                    ja = 2*k_ipjp-1 
                    cols(colPtr+11) = ja
                    vals(colPtr+11) = vals(colPtr+11) + sa
                else if ((.not. yhi_db .or. j.lt. ny)  .and.    &
                         (.not. xhi_db .or. i.lt. nx)) then  ! not last cell
                    ! interior ghost cell not on domain boundary
                    hvc = q(5,i+1,j+1)
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                !    hvc = q(5,i+1,j+1)
                !    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else if (k_ijp > -1) then
                    ! at right boundary, use cell above instead
                    !ja = k_ijp + nD
                    !ja = 2*k_ijp 
                    ja = 2*k_ijp-1 
                    ! v at right bndry doesnt need negating
                    vals(colPtr+10) = vals(colPtr+10) + sa
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead, 
                    !ja = k_ipj + nD
                    !ja = 2*k_ipj
                    ja = 2*k_ipj-1
                    if (bc_yhi==3) sa = -sa  ! reflect V in y for wall bc
                    vals(colPtr+8) = vals(colPtr+8) + sa
                else if (xhi_db .and. i .eq. nx .and. .not. yhi_db) then
                      ! at top right corner of patch, which doesnt touch
                      ! top domain bndry only right
                      ! use cell ij instead with bc
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! no reflection of V at right wall
                      vals(colPtr+6) = vals(colPtr+6) + sa
                else
                    ! top-right corner
                    !ja = k_ij + nD
                    ja = 2*k_ij-1
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                    
                ! Matrix element row U(i,j) column V(i+1,j-1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjm > -1) then
                    !ja = k_ipjm + nD  
                    !ja = 2*k_ipjm 
                    ja = 2*k_ipjm-1 
                    cols(colPtr+2) = ja
                    vals(colPtr+2) = vals(colPtr+2) + sa
                else if ((.not. ylo_db .or. j>1)  .and.   &
                         (.not. xhi_db .or. i .lt. nx)) then
                    hvc = q(5,i+1,j-1)
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc   
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                !    hvc = q(5,i+1,j-1)
                !    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc   
                else if (k_ijm > -1) then
                    ! at right boundary but not bottom
                    !  use cell below instead
                    !ja = k_ijm + nD
                    !ja = 2*k_ijm 
                    ja = 2*k_ijm-1 
                    ! equations for v so doesnt need negating
                    vals(colPtr+1) = vals(colPtr+1) + sa
                else if (k_ipj > -1) then
                    ! at bottom boundary but not right bndry
                    ! use cell to right instead,
                    !ja = k_ipj + nD
                    !ja = 2*k_ipj
                    ja = 2*k_ipj-1
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                    vals(colPtr+8) = vals(colPtr+8) + sa
                else if (xhi_db .and. i .eq. nx .and. .not. ylo_db) then
                      ! at bottom right corner of patch that
                      ! touches right domain boundary but not
                      ! touch bottom
                      ! use cell ij instead with bc
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! no reflection of V at right wall
                      vals(colPtr+6) = vals(colPtr+6) + sa
                else
                    ! bottom-right corner
                    !ja = k_ij + nD
                    !ja = 2*k_ij 
                    ja = 2*k_ij-1 
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i-1,j+1):    
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjp > -1) then
                    !ja = k_imjp + nD
                    !ja = 2*k_imjp 
                    ja = 2*k_imjp-1 
                    cols(colPtr+9) = ja
                    vals(colPtr+9) = vals(colPtr+9) + sa
                else if ((.not. yhi_db .or. j.lt.ny) .and.  &
                         (.not. xlo_db .or. i .gt. 1)) then
                    ! missing cell not outside phys. domain
                    hvc = q(5,i-1,j+1)
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                !    hvc = q(5,i-1,j+1)
                !    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else if (k_ijp > -1) then
                    ! at left domain boundary, use cell above instead
                    ! to pw const extrap V 
                    !ja = k_ijp + nD
                    !ja = 2*k_ijp 
                    ja = 2*k_ijp-1 
                    vals(colPtr+10) = vals(colPtr+10) + sa
                else if (k_imj > -1) then
                    ! at top domain boundary, use cell to left instead,
                    ! since at V bndry impose V bc and reflect
                    !ja = k_imj + nD
                    !ja = 2*k_imj 
                    ja = 2*k_imj-1 
                    if (bc_yhi==3) sa = -sa  ! reflect V in y
                    vals(colPtr+4) = vals(colPtr+4) + sa
                else if (xlo_db .and. i .eq. 1 .and. .not. yhi_db) then
                      ! at top left boundary of patch, left bndry of domain
                      ! not top bndry of domain or would have to use y bcs, thats next clause
                      !use cell ij instead with bc
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! no reflection of V at left wall
                      vals(colPtr+6) = vals(colPtr+6) + sa
                else
                    ! top-left corner of domain
                    ! put in V_ij stencil
                    !ja = k_ij + nD
                    !ja = 2*k_ij 
                    ja = 2*k_ij-1 
                    if (bc_yhi==3) sa = -sa  ! reflect V in y but not x
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                
                
                ! Matrix element row U(i,j) column V(i-1,j-1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjm > -1) then
                    !ja = k_imjm + nD
                    !ja = 2*k_imjm 
                    ja = 2*k_imjm-1 
                    cols(colPtr+0) = ja
                    vals(colPtr+0) = vals(colPtr+0) + sa
                else if ((.not. xlo_db .or. i>1) .and. (.not. ylo_db .or. j>1)) then
                    ! use ghost cells, outside phys. domain
                    hvc = q(5,i-1,j-1)
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*hvc
                else if (k_imj > -1) then 
                    ! at bottom domain boundary, use cell to left instead with bc,
                    !ja = k_imj + nD
                    !ja = 2*k_imj 
                    ja = 2*k_imj-1 
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                    vals(colPtr+4) = vals(colPtr+4) + sa
                else if (k_ijm > -1) then
                    ! at left domain boundary, use cell below instead
                   !ja = k_ijm + nD
                   !ja = 2*k_ijm
                   !write(*,*)"shouldnt be in this clause"
                   ja = 2*k_ijm-1
                   vals(colPtr+1) = vals(colPtr+1) + sa
                else if (xlo_db .and. i .eq. 1 .and. .not. ylo_db) then
                      ! at left boundary, use cell ij instead with bc
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! no reflection of V at left wall
                      vals(colPtr+6) = vals(colPtr+6) + sa
                else
                    ! bottom-left corner of domain
                    !ja = k_ij + nD
                    !ja = 2*k_ij
                    ja = 2*k_ij-1
                    if (bc_ylo==3) sa = -sa  ! reflect V in y
                    vals(colPtr+6) = vals(colPtr+6) + sa
                endif
                

                !-----------------------------------------------
                ! -D21 block (cross derivatives
                ! row ia = k_ij + nD  for hv corrections in cell (i,j)
                ! columns k_ipjp and other neighbors for hu corrections

                ! Matrix element row V(i,j) column U(i,j):
                ! NEW MATRIX ELEMENT FOR SGN WITH VARYING TOPO
                sa = alpha*(0.5d0*h_ij*Bxy(i,j) + etay*Bx(i,j))
                cols(colPtr+17) = 2*(k_ij-1)
                vals(colPtr+17) = vals(colPtr+17)+sa
                
                ! Matrix element row V(i,j) column U(i+1,j):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = -alpha*h_ij*(hy + 0.5d0*By(i,j))/(2.d0*dxm)
                ja = -1
                if (k_ipj > -1) then
                    !ja = k_ipj
                    !ja = 2*k_ipj-1
                    ja = 2*(k_ipj-1)
                    cols(colPtr+19) = ja
                    vals(colPtr+19) = vals(colPtr+19) + sa
                else if (.not. xhi_db) then
                    huc = q(4,i+1,j)
                    !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                    !rhs(2*k_ij) = rhs(2*k_ij) - sa*huc
                    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                else
                    ! at right boundary, use cell ij instead,
                    !ja = k_ij
                    !ja = 2*k_ij-1
                    ja = 2*(k_ij-1)
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i-1,j):
                ! NEW MATRIX ELEMENT FOR SGN
                sa = alpha*h_ij*(hy + 0.5d0*By(i,j))/(2.d0*dxm)
                ja = -1
                if (k_imj > -1) then
                    !ja = k_imj
                    !ja = 2*k_imj-1
                    ja = 2*(k_imj-1)
                    cols(colPtr+16) = ja
                    vals(colPtr+16) = vals(colPtr+16) + sa
                else if (.not. xlo_db) then
                    huc = q(4,i-1,j)
                    !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                    !rhs(2*k_ij) = rhs(2*k_ij) - sa*huc
                    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                else
                    ! at left boundary, use cell ij instead,
                    !ja = k_ij
                    !ja = 2*k_ij-1
                    ja = 2*(k_ij-1)
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i,j+1):
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = alpha*0.5d0*h_ij*Bx(i,j) / (2.d0*dym)
                ja = -1
                if (k_ijp > -1) then
                    !ja = k_ijp ! FIXED BUG2
                    !ja = 2*k_ijp-1 ! FIXED BUG2
                    ja = 2*(k_ijp-1) ! FIXED BUG2
                    cols(colPtr+21) = ja
                    vals(colPtr+21) = vals(colPtr+21) + sa
                else if (.not. yhi_db) then
                    huc = q(4,i,j+1)
                    !rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*huc
                else
                    ! at top boundary, use cell ij instead,
                    !ja = k_ij
                    !ja = 2*k_ij-1
                    ja = 2*(k_ij-1)
                    ! no reflection of U at top wall
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                ! Matrix element row V(i,j) column U(i,j-1): ! FIXED BUG2
                ! NEW MATRIX ELEMENT FOR SGN WITH TOPO
                sa = -alpha*0.5d0*h_ij*Bx(i,j) / (2.d0*dym)
                ja = -1
                if (k_ijm > -1) then
                    !ja = k_ijm ! FIXED BUG2
                    !ja = 2*k_ijm-1 ! FIXED BUG2
                    ja = 2*(k_ijm-1) ! FIXED BUG2
                    cols(colPtr+13) = ja
                    vals(colPtr+13) = vals(colPtr+13) + sa
                else if (.not. ylo_db) then
                    huc = q(4,i,j-1)
                    !rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                    rhs(2*(k_ij-1)) = rhs(2*(k_ij-1)) - sa*huc
                else
                    ! at bottom boundary, use cell ij instead,
                    !ja = k_ij ! FIXED BUG2
                    !ja = 2*k_ij-1 ! FIXED BUG2
                    ja = 2*(k_ij-1) ! FIXED BUG2
                    ! no reflection of U at bottom wall
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i+1,j+1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjp > -1) then
                    !ja = k_ipjp
                    !ja = 2*k_ipjp-1
                    ja = 2*(k_ipjp-1)
                    cols(colPtr+23) = ja
                    vals(colPtr+23) = vals(colPtr+23) + sa
                else if ((.not. yhi_db .or. j.lt.ny)  .and.   &
                         (.not. xhi_db .or. i.lt.nx)) then
                    huc = q(4,i+1,j+1)
                    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                !    huc = q(4,i+1,j+1)
                !    !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                !    !rhs(2*k_ij) = rhs(2*k_ij) - sa*huc
                !    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                else if (k_ijp > -1) then
                    ! at right boundary but not corner
                    ! use cell above instead,
                    ! need u bcs.
                    !ja = k_ijp
                    !ja = 2*k_ijp-1
                    ja = 2*(k_ijp-1)
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                    vals(colPtr+21) = vals(colPtr+21) + sa
                else if (k_ipj > -1) then
                    ! at top boundary, use cell to right instead,
                    !ja = k_ipj
                    !ja = 2*k_ipj-1
                    ja = 2*(k_ipj-1)
                    vals(colPtr+19) = vals(colPtr+19) + sa
                else if (yhi_db .and. j .eq. ny .and. .not. xhi_db) then
                      ! at top right boundary of patch which touches top
                      ! of domain. Just use u value, no bc for u needed
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      vals(colPtr+17) = vals(colPtr+17) + sa
                else
                    ! top-right corner
                    !ja = k_ij ! FIXED BUG
                    !ja = 2*k_ij-1 ! FIXED BUG
                    ja = 2*(k_ij-1) ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i+1,j-1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_ipjm > -1) then
                    !ja = k_ipjm
                    !ja = 2*k_ipjm-1
                    ja = 2*(k_ipjm-1)
                    cols(colPtr+15) = ja
                    vals(colPtr+15) = vals(colPtr+15) + sa
                else if ((.not. ylo_db .or. j>1) .and.       &
                         (.not. xhi_db .or. i .lt. nx)) then
                    ! interior ghost cell, not touching domain bndry
                    huc = q(4,i+1,j-1)
                    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_ipj == -1) .and. (.not. xhi_db))) then
                !    huc = q(4,i+1,j-1)
                !    !rhs(k_ij+nD) = rhs(k_ij+nD) - sa*huc
                !    !rhs(2*k_ij) = rhs(2*k_ij) - sa*huc
                !    rhs(2*k_ij-1) = rhs(2*k_ij-1) - sa*huc
                else if (k_ijm > -1) then
                    ! at right boundary but not phys domain bottom
                    ! use cell below instead,
                    !ja = k_ijm
                    !ja = 2*k_ijm-1
                    ja = 2*(k_ijm-1)
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                    vals(colPtr+13) = vals(colPtr+13) + sa
                else if (k_ipj > -1) then
                    ! at bottom boundary, use cell to right instead,
                    !ja = k_ipj
                    !ja = 2*k_ipj-1
                    ja = 2*(k_ipj-1)
                    vals(colPtr+19) = vals(colPtr+19) + sa
                else if (ylo_db .and. j .eq. 1 .and. .not. xhi_db) then
                    ! in corner of patch, which touches bottom of domain
                    ! but not right bndry
                    ! no reflection of u at bottom wall bndry
                    vals(colPtr+17) = vals(colptr+17) + sa
                else
                    ! bottom-right corner
                    !ja = k_ij ! FIXED BUG
                    !ja = 2*k_ij-1 ! FIXED BUG
                    ja = 2*(k_ij-1) ! FIXED BUG
                    if (bc_xhi==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
                
                ! Matrix element row V(i,j) column U(i-1,j+1):
                sa = alpha*h_ij**2 / (3.d0*dxdy4)
                ja = -1
                if (k_imjp > -1) then
                    !ja = k_imjp
                    !ja = 2*k_imjp-1
                    ja = 2*(k_imjp-1)
                    cols(colPtr+20) = ja
                    vals(colPtr+20) = vals(colPtr+20) + sa
                else if ((.not. yhi_db .or. j<ny)  &
                   .and. (.not. xlo_db .or. i>1)) then
                    huc = q(4,i-1,j+1)
                    rhs(2*k_ij-1) = rhs(2*k_ij-1)  - sa*huc
                !else if (((k_ijp == -1) .and. (.not. yhi_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                !    huc = q(4,i-1,j+1)
                !    !rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                !    !rhs(2*k_ij) = rhs(2*k_ij)  - sa*huc
                !    rhs(2*k_ij-1) = rhs(2*k_ij-1)  - sa*huc
                else if (k_ijp > -1) then
                    ! at left domain boundary but not top bndry
                    !  use cell above but for u eq apply bcs
                    !ja = k_ijp
                    !ja = 2*k_ijp-1
                    ja = 2*(k_ijp-1)
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                    vals(colPtr+21) = vals(colPtr+21) + sa
                else if (k_imj > -1) then
                    ! at top boundary, use cell to left instead,
                    !ja = k_imj
                    !ja = 2*k_imj-1
                    ja = 2*(k_imj-1)
                    vals(colPtr+16) = vals(colPtr+16) + sa
                else if (yhi_db .and. j .eq. ny .and. .not. xlo_db) then
                      ! at left boundary, use cell ij instead with bc
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! no reflection of U needed since not at left wall
                      vals(colPtr+17) = vals(colPtr+17) + sa
                else
                    ! top-left corner
                    !ja = k_ij !FIXED BUG
                    !ja = 2*k_ij-1 !FIXED BUG
                    ja = 2*(k_ij-1) !FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                

                ! Matrix element row V(i,j) column U(i-1,j-1):
                sa = -alpha*h_ij**2 / (3.d0*dxdy4)  
                ja = -1
                if (k_imjm > -1) then
                    !ja = k_imjm
                    !ja = 2*k_imjm-1
                    ja = 2*(k_imjm-1)
                    cols(colPtr+12) = ja
                    vals(colPtr+12) = vals(colPtr+12) + sa
                else if ((.not. ylo_db .or. j>1)  .and. (.not. xlo_db .or. i>1)) then
                   ! interior ghost cell. 0 is first interior cell touching bndry
                    huc = q(4,i-1,j-1)
                    rhs(2*k_ij-1) = rhs(2*k_ij-1)  - sa*huc
                !else if (((k_ijm == -1) .and. (.not. ylo_db)) .or. &
                !         ((k_imj == -1) .and. (.not. xlo_db))) then
                !    huc = q(4,i-1,j-1)
                !    !rhs(k_ij+nD) = rhs(k_ij+nD)  - sa*huc
                !    !rhs(2*k_ij) = rhs(2*k_ij)  - sa*huc
                !    rhs(2*k_ij-1) = rhs(2*k_ij-1)  - sa*huc
                else if (k_ijm > -1) then
                    ! at left boundary, use cell below instead,
                    !ja = k_ijm
                    !ja = 2*k_ijm-1
                    ja = 2*(k_ijm-1)
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                    vals(colPtr+13) = vals(colPtr+13) + sa
                else if (k_imj > -1) then
                    ! at bottom boundary, use cell to left instead,
                    !ja = k_imj
                    !ja = 2*k_imj-1
                    ja = 2*(k_imj-1)
                    vals(colPtr+16) = vals(colPtr+16) + sa
                else if (ylo_db .and. j .eq. 1 .and. .not. xlo_db) then
                      ! at bottom phys domain  boundary
                      ! in corner of patch, not left phys domain
                      !ja = k_ij + nD
                      !ja = 2*k_ij
                      ja = 2*k_ij-1
                      ! just use u from cell i-1,j above for v equation
                      ! dont need u bc
                      vals(colPtr+17) = vals(colPtr+17) + sa
                else
                    ! bottom-left corner
                    !ja = k_ij !FIXED BUG
                    !ja = 2*k_ij-1 !FIXED BUG
                    ja = 2*(k_ij-1) !FIXED BUG
                    if (bc_xlo==3) sa = -sa  ! reflect U in x
                    vals(colPtr+17) = vals(colPtr+17) + sa
                endif
                
 
                !-------------------------------------------------------
                ! Right hand sides for vector elements k_ij and k_ij+nD
                ! Need to add comments on SGN RHS.
                
              !RHS for SGN 
              phix = (phi(i+1,j) - phi(i-1,j))/(2.d0*dxm)
              phiy = (phi(i,j+1) - phi(i,j-1))/(2.d0*dym)
                            
              wx = (w(i+1,j) - w(i-1,j))/(2.d0*dxm)
              wy = (w(i,j+1) - w(i,j-1))/(2.d0*dym)
              
              etax = (eta(i+1,j) - eta(i-1,j))/(2.d0*dxm)
              etay = (eta(i,j+1) - eta(i,j-1))/(2.d0*dym)
                                 
              rhs(2*(k_ij-1)) = rhs(2*(k_ij-1))  +  (grav/alpha * etax      &
                    + 2.d0*h_ij*(h_ij/3.d0 * phix &
                    + phi(i,j) * (hx + 0.5d0*Bx(i,j)))   &
                    + 0.5d0*h_ij*wx + w(i,j)*etax)

              rhs(2*k_ij-1) = rhs(2*k_ij-1) +  (grav/alpha * etay &
                    + 2.d0*h_ij*(h_ij/3.d0 * phiy &
                    + phi(i,j) * (hy + 0.5d0*By(i,j)))   &
                    + 0.5d0*h_ij*wy + w(i,j)*etay)

                    
              ! if ghost cells in stencil might not be sorted so fix
              ! interior cells numbered so CSR entries are sorted
              ! this row's entries start at rowPtr(colPtr+0) to colPtr+11
              ! for u, and colPtr+12 to colPtr+23 for v
              if (IS_GHOST(i,j)) then 
                 indexSt = rowPtr(2*(k_ij-1))
                 call insertionSort(cols(indexSt),vals(indexSt))
                 indexSt = rowPtr(2*k_ij-1)
                 call insertionSort(cols(indexSt),vals(indexSt))
              endif


              if (0 .eq. 1 .or. debug) then
                 allzero = .true.
                 do n = 0, 23
                    if (vals(colPtr+n) .ne. 0.d0) then
                       allzero = .false.
                       exit
                    endif
                 end do
                 if (allzero) then
                    write(*,505) i,j,k_ij
 505                format(" u,v, rows for cell ",2i5," index ",i5," all zeroes")
                 endif

              endif
                     
              ! put initial guess for k_ij and k_ij+nD respective in soln array
             !if (levelBouss .lt. mxnest) then
             !  ! initial guess when doUpdate is true is previous soln
             !  ! when doUpdate was false, but values needed for finer
             !  ! grid interface conditions
             !  soln(k_ij) = qold(4,i,j) 
             !  soln(k_ij+nD) = qold(5,i,j)
             !else
             !  soln(k_ij) = 0.d0
             !  soln(k_ij+nD) = 0.d0 
             !endif

            enddo
        enddo

        if (debug) then ! dump matrix to look for singularity 
           !write(89,*)" level ",levelBouss
           !do k = 1,2*numBoussCells
           !   write(89,104) k,rhs(k)
 104          format(i5,e16.7)
           !end do
           !!close(89)
           !!stop
           
           write(69,*)" grid ",mptr," level ",levelBouss,"i,j,etax,etay,rhsu,rhsv"
           !do ii = 0, nx+1
           !do jj = 0, ny+1
           do ii = 1, nx
           do jj = 1, ny
                k_ij = mi%mindex(ii,jj)
                x = xlow + (ii-.5d0)*dx
                y = ylow + (jj-.5d0)*dy
                if (k_ij .ne. -1) then
                  bu = rhs(2*(k_ij-1))
                  bv = rhs(2*k_ij-1)
                  write(69,798) x,y,etaxM(ii,jj),etayM(ii,jj),bu,bv,ii,jj,k_ij
                else
                  bu = 0.d0
                  bv = 0.d0
                  write(69,798) x,y,etaxM(ii,jj),etayM(ii,jj),bu,bv,ii,jj,k_ij
                endif
              !write(69,799) ii,jj,etaxM(ii,jj),etayM(ii,jj),bu,bv
              !write(69,798) x,y,etaxM(ii,jj),etayM(ii,jj),bu,bv
 799          format(2i5,4e15.7)
 798          format(6e15.7,3i5)
           end do
           end do
        endif
        
    
#ifdef WHERE_AM_I
  write(*,*) "ending   buildSparseMatrixSGNcrs"
#endif


    return
end subroutine buildSparseMatrixSGNcrs

subroutine insertionSort(cols,vals)

    integer, intent(inout) :: cols(0:11)
    real(kind=8), intent(inout) :: vals(0:11)

    integer :: i,j, colSave
    real(kind=8) :: valSave

! in place insertion sort for 12 element, starting index is indexSt
    do i = 1, 11
      colSave = cols(i)
      valSave = vals(i)

      do j = i-1, 0, -1
        if (colSave < cols(j)) then
           cols(j+1) = cols(j)
           vals(j+1) = vals(j)
        else
            exit
        endif
      end do
      cols(j+1) = colSave 
      vals(j+1) = valSave
    end do

   return
end subroutine insertionSort
