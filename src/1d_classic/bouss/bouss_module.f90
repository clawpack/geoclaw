module bouss_module
    
    implicit none

    ! parameters set in setrun:
    integer :: boussEquations       ! Which equations? 0=SWE, 1=Madsen, 2=SGN
    real(kind=8) :: boussMinDepth   ! depth below which to switch to SWE

    integer :: bc_xlo, bc_xhi
    real(kind=8) :: B_param, alpha

    ! for higher-order derivatives:
    real(kind=8), allocatable, dimension(:) :: rinv, dxm, dxc, cm2r, cp2r, c02r
    real(kind=8), allocatable, dimension(:) :: cm2, cp2, c02, cm1, cp1, c01

    ! for initial water depth and switching to SWE:
    real(kind=8) :: hmin
    real(kind=8), allocatable, dimension(:) :: h0
    logical, allocatable, dimension(:) :: useBouss
    
    ! for tridiagonal solver:
    integer, allocatable, dimension(:) :: IPIV
    real(kind=8), allocatable, dimension(:) :: D, DL, DU, DU2    

    save


contains

!======================================================================

    subroutine set_bouss(mx,mbc,mthbc)
    
    ! Set things up for Boussinesq solver, in particular
    ! create and factor tridiagonal matrix for implicit solves.

    use geoclaw_module, only: earth_radius,deg2rad,coordinate_system,sea_level
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2)
    real(kind=8) :: DLi, r, denom
    integer :: i, iunit
    character*25 fname
    real(kind=8), dimension(1-mbc:mx+mbc) :: h02, h03
    integer(kind=4) :: INFO

    iunit = 7
    fname = 'bouss.data'
 !  # open the unit with new routine from Clawpack 4.4 to skip over
 !  # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) boussEquations  ! which equations
    read(7,*) boussMinDepth   ! depth below which to switch to SWE
    
    if (boussEquations==1) B_param = 1.d0/15.d0  ! parameter in MS equations
    if (boussEquations==2) alpha = 1.153  ! parameter in SGN equations
    
    allocate(rinv(0:mx+1),dxm(0:mx+2),dxc(0:mx+1))
    allocate(cm1(0:mx+1),cp1(0:mx+1),c01(0:mx+1))
    allocate(cm2(0:mx+1),cp2(0:mx+1),c02(0:mx+1))
    allocate(cm2r(0:mx+1),cp2r(0:mx+1),c02r(0:mx+1))
    allocate(D(mx+2), DL(mx+1), DU(mx+2), DU2(mx+2), IPIV(mx+2))
    allocate(h0(-1:mx+2), useBouss(-1:mx+2))


    ! Boundary conditions to impose in computing Boussinesq update:
    
    if ((mthbc(1)==2) .or. (mthbc(2)==2)) then
        write(6,*) '*** Periodic BCs not supported in bouss_module'
        stop
    endif

    ! Use Dirichlet BC with value 0 for correction in ghost cells by default:
    bc_xlo = 0
    bc_xhi = 0
    
    ! Check if wall BCs specified (also be used for radial symmetry at r=0):
    
    if (mthbc(1)==3) then
        ! For wall boundary conditions at left boundary:
        bc_xlo = 3
    endif

    if (mthbc(2)==3) then
        ! For wall boundary conditions at right boundary:
        bc_xhi = 3
    endif

    if (mthbc(1)==0) then
        ! For wavemaker BC at left
        bc_xlo = 3
    endif

    ! To try out Neumann BCs:
    ! This doesn't seem to work well, so not a general option
    !bc_xlo = 1
    !bc_xhi = 1
    
    
    !------------------------------------
    ! coefficients needed for second-order derivatives:

    do i=0,mx+2
        if (coordinate_system .eq. 2) then
            ! x is latitude, convert dx to meters:
            dxm(i) = (xcell(i) - xcell(i-1)) * earth_radius * deg2rad
        else
            ! x is meters
            dxm(i) = xcell(i) - xcell(i-1)
        endif
    enddo

    do i=0,mx+1
        if (grid_type == 0) then
            ! uniform grid
            dxc(i) = 2*dxm(i)
            cm2(i) = 1.d0/dxm(i)**2
            c02(i) = -2.d0/dxm(i)**2
            cp2(i) = 1.d0/dxm(i)**2
            ! first derivative coefficients:
            cm1(i) = -1.d0/dxc(i)
            c01(i) = 0.d0
            cp1(i) = 1.d0/dxc(i)
        else
            dxc(i) = (xcell(i+1) - xcell(i-1))  ! = 2*dx for uniform grid!
            cm2(i) = 2.d0 / (dxm(i)*dxc(i))
            cp2(i) = 2.d0 / (dxm(i+1)*dxc(i))
            c02(i) = -(cm2(i) + cp2(i))
            ! first derivative coefficients:
            denom = (dxm(i)*dxm(i+1)*(dxm(i)+dxm(i+1)))
            ! CORRECT??
            cm1(i) = -dxm(i+1)**2 / denom
            c01(i) = (dxm(i+1)**2 - dxm(i)**2) / denom
            cp1(i) = +dxm(i)**2 / denom
        endif
        
        if (coordinate_system == 1) then
            cm2r(i) = cm2(i)
            cp2r(i) = cp2(i)
            c02r(i) = c02(i)
            rinv(i) = 0.d0
        else if (coordinate_system == -1) then
            ! x = radial coordinate r >= 0
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = xcell(i)
            rinv(i) = 1.d0/r
            cm2r(i) = cm2(i) - 1.d0/(r*dxc(i))
            cp2r(i) = cp2(i) + 1.d0/(r*dxc(i))
            c02r(i) = c02(i) - 1.d0/(r**2)
            ! first derivative coefficients for d/dr are unchanged 
            ! note: *not* including 1/r term in this derivative
        else if (coordinate_system == 2) then
            write(6,*) '*** latitude coordinates correctly implemented in Bouss??'
            stop
            
            ! PROBABLY NOT CORRECT -- r derivative??
            ! x = latitude coordinate -90 < x < 90 degrees
            ! r = sin(x) goes to zero at both poles
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = sin(deg2rad*xcell(i))
            rinv(i) = 1.d0/r
            cm2r(i) = cm2(i) - 1.d0/(r*dxc(i))
            cp2r(i) = cp2(i) + 1.d0/(r*dxc(i))
            c02r(i) = c02(i) - 1.d0/(r**2)
            ! first derivative coefficients for d/dr are unchanged 
            ! note: *not* including 1/r term in this derivative        
        endif
    enddo

    ! initial depth and where to switch to SWE:
    
    do i=1-mbc,mx+mbc
        h0(i) = max(sea_level - zcell(i), 0.d0)
    enddo
    
    do i=2-mbc,mx+mbc-1
        hmin = min(h0(i-1), h0(i), h0(i+1))
        useBouss(i) = (hmin > boussMinDepth)
        !write(6,*) '+++ useBouss: ',i,xcell(i),h0(i),useBouss(i)
    enddo        
    !------------------------------------
    
    if (boussEquations==1) then
    
        ! Madsen-Sorensen
        ! Form tridiagonal matrix and factor
        ! Done once at start and then used in each time step
        ! (for SGN, we need to build in each time step)

        call build_tridiag_ms(mx,mbc,mthbc)

    endif

    end subroutine set_bouss

!======================================================================
! Routines for Madsen-Sorensen equations
!======================================================================

subroutine build_tridiag_ms(mx,mbc,mthbc)

    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2)
    
    integer :: i
    real(kind=8) :: DLi, r
    real(kind=8), dimension(1-mbc:mx+mbc) :: h02, h03
    integer(kind=4) :: INFO
    
    !h0 = sea_level - aux(1,:)  # resting depth

    do i=1-mbc,mx+mbc
        !h0(i) = sea_level - zcell(i) ! now set in set_bouss
  666   format(i3,3e16.6)
        h02(i)=(max(0.,h0(i)))**2
        h03(i)=(max(0.,h0(i)))**3
    enddo

    ! initialize to identity matrix:
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell) to form (I-d^2) operator,
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        if ((h0(i) > boussMinDepth)) then
            
          ! Replace this row of identity matrix with dispersion terms
          ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
          ! approximation to d^2 w = w_{xx} in plane wave case, or
          !                  d^2 w = (1/r * (r*w)_r)_r  in radial case
          ! Also h02 is h0**2 and h03 is h0**3:

          D(i+1) = 1.d0 - c02r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i))

          DU(i+1)= -cp2r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i+1))

          DL(i)= -cm2r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i-1))

        !else
        !  write(66,*) 'i, h0(i): ',i,h0(i)
        endif
        
    enddo

    ! left boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_0 = 0
    ! For other BCs, change superdiagonal in row 1:

    if (bc_xlo==1) then
        ! for Neumann BC at left: impose Q_0 = Q_1
        DU(1) = -1.d0     
    endif
    if (bc_xlo==3) then
        ! for wall-reflecting BC at left: impose Q_0 = -Q_1
        DU(1) = 1.d0     
    endif


    ! right boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_{mx+1} = 0
    ! For other BCs, change subdiagonal in row mx+2:

    if (bc_xhi==1) then
        ! for Neumann at right: impose Q_mx = Q_{mx+1}
        DL(mx+1) = -1.d0
    endif
    if (bc_xhi==3) then
        ! for wall at right: impose Q_mx = - Q_{mx+1}
        DL(mx+1) = 1.d0
    endif

    
    if (.false.) then
        ! DEBUG:
        write(66,*) 'D matrix:'
        do i=1,mx+2
            if (i>1) then
                DLi = DL(i-1)
            else
                DLi = -9999.  ! since no subdiag in first row
            endif
            write(66,661) DLi, D(i), DU(i)
  661       format(3e16.6)
        enddo
        write(66,*) '===== end of D'
    endif

    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )

    end subroutine build_tridiag_ms

!------------------------------------------------------------------

    subroutine solve_tridiag_ms(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in set_bouss, so the
    ! arrays D, DL, DU, DU2, IPIV are already set up for the solver.
    
    use geoclaw_module, only: g => grav, sea_level, dry_tolerance
    use geoclaw_module, only: coordinate_system
    use grid_module, only: xcell

    implicit none

    integer, intent(in) :: mx,meqn,mbc,maux
    real(kind=8), intent(in) :: dx

    real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(out) :: psi(mx+2)

    real(kind=8)  depth
    real(kind=8), dimension(1-mbc:mx+mbc) :: h0, h02, h03, eta, hu2
    real(kind=8), dimension(1-mbc:mx+mbc) :: h_eta_x,h0_eta_x,s1,s1_h0
    real(kind=8)  h0_eta_xxx,s1_xx,s1_h0_xx
    integer :: i,j,k,kk,iL,iR
    integer(kind=4) :: INFO
    integer :: LDB

    h0 = sea_level - aux(1,:)

    do i=1-mbc,mx+mbc
        ! compute hu2 = h*u**2
        if (q(1,i).gt.dry_tolerance) then
            hu2(i)= q(2,i)**2/q(1,i)
        else
            hu2(i)= 0.d0
        endif
        
        !hu2(i) = 0.d0  ! Debug: turning off nonlinear terms
        
        ! h0 is resting depth, h02 = h0**2 and h03 = h0**3:
        eta(i)= q(1,i) - h0(i)
        h02(i)=(max(0.d0,h0(i)))**2
        h03(i)=(max(0.d0,h0(i)))**3
        
    enddo

    ! compute these quantities, extending to one set of ghost cells:
    ! h_eta_x  = h * eta_x
    ! h0_eta_x = h0 * eta_x
    
    ! defaults to 0 if h0==0 at one of the stencil points:
    h_eta_x = 0.d0
    h0_eta_x = 0.d0
        
    do i= 0,mx+1
        if (minval(h0(i-1:i+1)) > 0.d0) then
            ! if h0 > 0 at all stencil points.  Note dxc(i) = 2*dx if uniform:
            h_eta_x(i)= q(1,i)*(eta(i+1)-eta(i-1))/dxc(i)
            h0_eta_x(i)= h0(i)*(eta(i+1)-eta(i-1))/dxc(i)
        endif
    enddo
     

    
    do i=0,mx+1
    
        ! compute approximation to
        !   s1 = (h*u**2)_x + g*h*eta_x   (plane wave) or
        !   s1 = (1/r)*(r*h*u**2)_r + r*g*h*eta_r   (radial)
        !      = (h*u**2)_r + (1/r)*h*u**2 + r*g*h*eta_r
    
        s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) + g*h_eta_x(i)
        !s1(i) =  g*h_eta_x(i)  ! debug
        !s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) 
                 
        if (coordinate_system == -1) then
            ! radial: add term for (1/r)*h*u**2:
            s1(i) =  s1(i) + hu2(i) / xcell(i) 
        else if (coordinate_system == 2) then
            ! latitude
            write(6,*) '*** latitude not yet implemented in bouss_module'
            stop
        endif
        
        ! compute s1_h0 = s1 / h0:
        
        if (h0(i) > dry_tolerance) then
            s1_h0(i)=s1(i)/h0(i)
        else
            s1_h0(i)=0.d0
        endif
    enddo
     

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS

       if (h0(i) <= boussMinDepth) then
            ! no dispersive term:
            psi(k) = 0.d0

        else
            ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
            ! approximation to d^2 w = w_{xx} 
            ! or in radial case to d^2 w = (1/r * (rw)_r)_r

            ! compute s1_xx, approximates d^2 s1:

            s1_xx = cp2r(i)*s1(i+1) + c02r(i)*s1(i) + cm2r(i)*s1(i-1)

            ! compute s1_h0_xx, approximates d^2 s1_h0 = d^2(s1 / h0):
            
            s1_h0_xx = cp2r(i)*s1_h0(i+1) + c02r(i)*s1_h0(i) + cm2r(i)*s1_h0(i-1)
        
            ! compute h0_eta_xxx, approximates d^2 h0_eta_x:
            
            h0_eta_xxx = cp2r(i)*h0_eta_x(i+1) + c02r(i)*h0_eta_x(i) + &
                    cm2r(i)*h0_eta_x(i-1)
                    

            ! Right-hand side:
            
            psi(k) = -(B_param+.5d0) * h02(i) * s1_xx &
                        + h03(i)/6.d0 * s1_h0_xx &
                        + B_param * g * h02(i) * h0_eta_xxx
            if (.false.) then
                write(66,*) '+++i,rhs(i+1): ',i, psi(k)
                write(66,*) (B_param+.5d0) * h02(i), s1_xx
                write(66,*) h03(i)/6.d0,  s1_h0_xx 
                write(66,*) B_param * g * h02(i) * h0_eta_xxx
                write(66,*) B_param,  g,  h02(i),  h0_eta_xxx
                write(66,*) s1(i)
                write(66,*) h0_eta_x(i)
                write(66,*) eta(i)
            endif
              
                               
        endif


        ! OLD NOTATION, NOT USED NOW:
        !if ((h0(i) > sw_depth0) .and. (h0(i) < sw_depth1)) then
        !    ! reduce psi linearly between sw_depth0 and sw_depth1:
        !    psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
        !endif
            
    enddo
    
    if (.false.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) '+++i,rhs: ',i, psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, IPIV, psi, LDB, INFO )
                                
    if (.false.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) '+++ i,psi: ',i,psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag_ms



!======================================================================
! Routines for SGN equations
!======================================================================

subroutine build_tridiag_sgn(meqn,mbc,mx,xlower,dx,q,maux,aux)

    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)
    
    ! Locals
    integer :: i
    real(kind=8) :: DLi, r, Brr, denom, hmin
    real(kind=8), dimension(1-mbc:mx+mbc) :: h, B, Br, hr
    real(kind=8), dimension(1-mbc:mx+mbc) :: etar, eta, phi, w

    integer(kind=4) :: INFO

    do i=1-mbc,mx+mbc
        h(i) = q(1,i)  ! depth
        B(i) = aux(1,i)  ! topo
        eta(i) = B(i) + h(i)
    enddo
    
    ! reset useBouss(:) if it should be based on h rather than h0 or slope, etc:
    if (.false.) then
        do i=2-mbc,mx+mbc-1
            hmin = min(h(i-1), h(i), h(i+1))
            useBouss(i) = (h(i) > boussMinDepth)
        enddo
    endif
    
    do i=2-mbc,mx+mbc-1
        Br(i) = cm1(i)*B(i-1) + c01(i)*B(i) + cp1(i)*B(i+1)
        hr(i) = cm1(i)*h(i-1) + c01(i)*h(i) + cp1(i)*h(i+1)
        etar(i) = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
    enddo

    ! initialize to identity matrix:
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell) to form (I-d^2) operator,
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        !if ((h(i) > boussMinDepth)) then
        if (useBouss(i)) then
            
          ! Replace this row of identity matrix with dispersion terms
          ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
          ! approximation to d^2 w = w_{xx} in plane wave case, or
          !                  d^2 w = (1/r * (r*w)_r)_r  in radial case
          ! while cm2, c02, cp2 are coefficients for pure second derivative
          ! even in radial case.
          
          Brr = cm2(i)*B(i-1) + c02(i)*B(i) + cp2(i)*B(i+1)
          
          D(i+1) = 1.d0 + alpha*(-h(i)**2/3.d0 *c02r(i) &
                    -h(i)*hr(i)*(c01(i) + rinv(i)) &
                    + 0.5d0*h(i)*(Brr - rinv(i)*Br(i)) + etar(i)*Br(i)) !CORRECT??
          
          DU(i+1) = alpha*(-h(i)**2/3.d0 *cp2r(i) - h(i)*hr(i)*cp1(i))

          DL(i) = alpha*(-h(i)**2/3.d0 *cm2r(i) - h(i)*hr(i)*cm1(i))
        endif
    enddo

    ! left boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_0 = 0
    ! For other BCs, change superdiagonal in row 1:

    if (bc_xlo==1) then
        ! for Neumann BC at left: impose Q_0 = Q_1
        DU(1) = -1.d0     
    endif
    if (bc_xlo==3) then
        ! for wall-reflecting BC at left: impose Q_0 = -Q_1
        DU(1) = 1.d0     
    endif


    ! right boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_{mx+1} = 0
    ! For other BCs, change subdiagonal in row mx+2:

    if (bc_xhi==1) then
        ! for Neumann at right: impose Q_mx = Q_{mx+1}
        DL(mx+1) = -1.d0
    endif
    if (bc_xhi==3) then
        ! for wall at right: impose Q_mx = - Q_{mx+1}
        DL(mx+1) = 1.d0
    endif

    
    if (.false.) then
        ! DEBUG:
        write(66,*) 'D matrix:'
        do i=1,mx+2
            if (i>1) then
                DLi = DL(i-1)
            else
                DLi = -9999.  ! since no subdiag in first row
            endif
            write(66,661) DLi, D(i), DU(i)
  661       format(3e16.6)
        enddo
        write(66,*) '===== end of D'
    endif

    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )

    end subroutine build_tridiag_sgn

!------------------------------------------------------------------    

    subroutine solve_tridiag_sgn(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! For the SGN equations
    
    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in src1, so the
    ! arrays D, DL, DU, DU2, IPIV are already set up for the solver.
    
    use geoclaw_module, only: g => grav, sea_level, dry_tolerance
    use geoclaw_module, only: coordinate_system
    use grid_module, only: xcell

    implicit none

    integer, intent(in) :: mx,meqn,mbc,maux
    real(kind=8), intent(in) :: dx

    real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(out) :: psi(mx+2)

    real(kind=8), dimension(1-mbc:mx+mbc) :: h, eta, u, B, phi, w
    real(kind=8) :: Br, hr, etar, Brr, wr, phir, ur
    integer :: i,j,k,kk,iL,iR
    integer(kind=4) :: INFO
    integer :: LDB

    h = q(1,:) ! depth

    do i=1-mbc,mx+mbc
        if (q(1,i).gt.dry_tolerance) then
            u(i)= q(2,i)/q(1,i)
        else
            u(i)= 0.d0
        endif            

        B(i) = aux(1,i)  ! topo
        eta(i) = B(i) + h(i)
    enddo

        
    do i=0,mx+1
        Brr = cm2(i)*B(i-1) + c02(i)*B(i) + cp2(i)*B(i+1)  ! pure 2nd deriv
        w(i) = u(i)**2 * Brr
        ur = cm1(i)*u(i-1) + c01(i)*u(i) + cp1(i)*u(i+1)
        phi(i) = ur**2 + rinv(i)*ur*u(i) + (rinv(i)*ur)**2  ! CORRECT??
    enddo        
    

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS

       !if (h(i) <= boussMinDepth) then
       if (.not. useBouss(i)) then
            ! no dispersive term:
            psi(k) = 0.d0

        else
            Br = cm1(i)*B(i-1) + c01(i)*B(i) + cp1(i)*B(i+1)   ! pure 1st deriv
            hr = cm1(i)*h(i-1) + c01(i)*h(i) + cp1(i)*h(i+1)
            etar = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
            wr = cm1(i)*w(i-1) + c01(i)*w(i) + cp1(i)*w(i+1)
            phir = cm1(i)*phi(i-1) + c01(i)*phi(i) + cp1(i)*phi(i+1)
                    

            ! Right-hand side:
            
            psi(k) = g/alpha * etar &
                    + 2*h(i)*(h(i)/3.d0 * phir + phi(i)*(hr + 0.5d0*Br)) &
                    + 0.5d0*h(i)*wr + w(i)*etar
            
            if (.false.) then
                write(66,*) '+++i,rhs(i+1): ',i, psi(k)
            endif                       
        endif


        ! OLD NOTATION, NOT USED NOW:
        !if ((h(i) > sw_depth0) .and. (h(i) < sw_depth1)) then
        !    ! reduce psi linearly between sw_depth0 and sw_depth1:
        !    psi(k) = psi(k) * (h(i)-sw_depth0)/(sw_depth1-sw_depth0)
        !endif

    enddo
    
    if (.false.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) '+++i,rhs: ',i, psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, IPIV, psi, LDB, INFO )
                                
    if (.false.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) '+++ i,psi: ',i,psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag_sgn


end module bouss_module
