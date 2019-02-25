!> Compute all fluxes at cell edges 
!! \param qold[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2x(maxm,meqn,maux,mbc,mx,my, &
                 qold,aux,dx,dt,cflgrid, &
                 fm,fp,rpn2)
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
    real(kind=8), intent(in) :: dx,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)
    
    ! Looping scalar storage
    integer :: j,m
    real(kind=8) :: dtdx,cfl1d,p,phi,cm,dtdxij,dtdyij
    
    ! Parameters
    ! Relimit fluxes to maintain positivity
    logical, parameter :: relimit = .false.

    cflgrid = 0.d0
    dtdx = dt/dx
    
    fm = 0.d0
    fp = 0.d0

    ! ==========================================================================
    ! Perform X-Sweeps
    do j = 1-mbc,my+mbc

        ! Copy old q into 1d slice
        q1d(:,1-mbc:mx+mbc) = qold(:,1-mbc:mx+mbc,j)
        
        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif
        
        ! Copy aux array into slices
        if (maux > 0) then
            aux2(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j  )
        endif
        
        ! Store value of j along the slice into common block
        ! *** WARNING *** This may not working with threading

        ! Compute modifications fadd and gadd to fluxes along this slice:
        call flux2_dimSplit(1,maxm,meqn,maux,mbc,mx,q1d,dtdx1d,aux2, &
                   faddm,faddp,cfl1d,rpn2) 

        cflgrid = max(cflgrid,cfl1d)
        ! write(53,*) 'x-sweep: ',cfl1d,cflgrid

        ! Update fluxes
        fm(:,1:mx+1,j) = fm(:,1:mx+1,j) + faddm(:,1:mx+1)
        fp(:,1:mx+1,j) = fp(:,1:mx+1,j) + faddp(:,1:mx+1)

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

end subroutine step2x
