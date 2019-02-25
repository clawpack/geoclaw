!> Compute all fluxes at cell edges 
!! \param qold[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2y(maxm,meqn,maux,mbc,mx,my, &
                 qold,aux,dy,dt,cflgrid, &
                 gm,gp,rpn2)
!     ==========================================================
!     This version of step2 relimits the fluxes in order to
!     maintain positivity.
!     to do so set relimit=.true.
!

    use geoclaw_module, only: dry_tolerance
    use amr_module, only: mwaves, mcapa

    implicit none
    
    external rpn2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)
    
    
    ! Looping scalar storage
    integer :: i,m
    real(kind=8) :: dtdy,cfl1d,p,phi,cm,dtdxij,dtdyij
    

    ! Parameters
    ! Relimit fluxes to maintain positivity
    logical, parameter :: relimit = .false.

    cflgrid = 0.d0
    dtdy = dt/dy
    
    gm = 0.d0
    gp = 0.d0


    ! ============================================================================
    !  y-sweeps   (for Godunov split method, called from stepgrid_dimSplit)
    !  loop indices assume that x sweep has gone first, this is second sweep.
    do i = 1, mx
        
        ! Copy data along a slice into 1d arrays:
        q1d(:,1-mbc:my+mbc) = qold(:,i,1-mbc:my+mbc)

        ! Set dt/dy ratio in slice
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! Copy aux slices
        if (maux .gt. 0)  then
            aux2(:,1-mbc:my+mbc) = aux(:,i,1-mbc:my+mbc)
        endif
        
        ! Compute modifications fadd and gadd to fluxes along this slice
        call flux2_dimSplit(2,maxm,meqn,maux,mbc,my,q1d,dtdy1d,aux2, &
                   faddm,faddp,cfl1d,rpn2)

        cflgrid = max(cflgrid,cfl1d)
        ! write(53,*) 'y-sweep: ',cfl1d,cflgrid

        ! Update fluxes
        gm(:,i,1:my+1) = gm(:,i,1:my+1) + faddm(:,1:my+1)
        gp(:,i,1:my+1) = gp(:,i,1:my+1) + faddp(:,1:my+1)

    enddo

    ! Below is not updated for use when dimensional splitting is used
    ! Relimit correction fluxes if they drive a cell negative
    ! if (relimit) then
    !     dtdxij = dtdx
    !     dtdyij = dtdy
    !     do i=1,mx
    !         do j=1,my
    !             if (mcapa > 0) then
    !                 dtdxij = dtdx / aux(mcapa,i,j)
    !                 dtdyij = dtdy / aux(mcapa,i,j)
    !             endif
    !             p = max(0.d0,dtdxij*fm(1,i+1,j)) + max(0.d0,dtdyij*gm(1,i,j+1)) &
    !               - min(0.d0,dtdxij*fp(1,i,j)) - min(0.d0,dtdyij*gp(1,i,j))
    !             phi = min(1.d0,abs(qold(1,i,j) / (p+dry_tolerance)))

    !             if (phi < 1.d0) then
    !                 do m=1,meqn
    !                     if (fp(1,i,j) < 0.d0) then
    !                         cm = fp(m,i,j) - fm(m,i,j)
    !                         fm(m,i,j) = phi * fm(m,i,j)
    !                         fp(m,i,j) = fm(m,i,j) + cm
    !                     endif
    !                     if (gp(1,i,j) < 0.d0) then
    !                         cm = gp(m,i,j) - gm(m,i,j)
    !                         gm(m,i,j) = phi * gm(m,i,j)
    !                         gp(m,i,j) = gm(m,i,j) + cm
    !                     endif
    !                     if (fm(1,i+1,j) > 0.d0) then
    !                         cm = fp(m,i+1,j) - fm(m,i+1,j)
    !                         fp(m,i+1,j) = phi * fp(m,i+1,j)
    !                         fm(m,i+1,j) = fp(m,i+1,j) - cm
    !                     endif
    !                     if (gm(1,i,j+1) > 0.d0) then
    !                         cm = gp(m,i,j+1) - gm(m,i,j+1)
    !                         gp(m,i,j+1) = phi * gp(m,i,j+1)
    !                         gm(m,i,j+1) = gp(m,i,j+1) - cm
    !                     endif
    !                 end do
    !             endif
    !         enddo
    !     enddo
    ! endif

end subroutine step2y
