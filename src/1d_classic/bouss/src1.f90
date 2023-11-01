
subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)

    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version integrates manning friction or other friction terms if present
    
    ! Version from 1D branch of GeoClaw originally from Dave George

    ! Also handles radial source term, if desired.
    ! Note: assumes radial about lower boundary, wall bc should be imposed

    
    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient
    use geoclaw_module, only: earth_radius, coordinate_system
    use grid_module, only: xcell
    use bouss_module, only: boussEquations, alpha, cm1, c01, cp1, &
                      solve_tridiag_ms, build_tridiag_sgn, &
                      solve_tridiag_sgn, useBouss


    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    ! Locals
    real(kind=8) :: eta(0:mx+1)
    real(kind=8) :: gamma, rcell, u, tanxR, etax
    real(kind=8) :: rk_stage(1:mx,4), delt
    integer ::  i,k,ii,rk_order
    real(kind=8)  q0(meqn,1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)

    

      if (frictioncoeff.gt.0.d0 .and. friction_forcing) then
          ! integrate source term based on Manning formula
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
                  q(2,i)= q(2,i)/(1.d0 + dt*gamma)
              endif
            enddo
        endif

!      ----------------------------------------------------------------

    if (coordinate_system == -1) then
        ! radial source term for SWE:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                q(1,i) = q(1,i) - dt/xcell(i) * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) - dt/xcell(i) * q(1,i)*u**2
            endif
         enddo
     endif

!      ----------------------------------------------------------------

    if (coordinate_system == 2) then
        ! source term for x = latitude in degrees -90 <= x <= 90,
        ! u = velocity in latitude direction (m/s) on sphere:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                tanxR = tan(xcell(i)*DEG2RAD) / earth_radius
                q(1,i) = q(1,i) + dt * tanxR * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) + dt * tanxR * q(1,i)*u**2
            endif
         enddo
     endif

    ! -------------------------------------------------
    if (boussEquations > 0) then
        
        ! Boussinesq terms 

        rk_order = 1  ! 1 for Forward Euler, 2 for second-order RK
        delt = dt / rk_order

        if (boussEquations == 2) then
            ! for SGN, need to factor matrix each step
            call build_tridiag_sgn(meqn,mbc,mx,xlower,dx,q,maux,aux)
        endif
                
        q0  = q
          
        do i=0,mx+1
           if (q(1,i) > dry_tolerance) then
               eta(i) = q(1,i)+aux(1,i)
           else
               eta(i) = 0.d0
           endif
        enddo

        !-----------------------
        ! First stage (only stage for rk_order == 1):
        
        if (boussEquations == 1) then
          call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! returns solution psi = source term for Madsen
        else if (boussEquations == 2) then
          call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! modify solution psi for source term of SGN:
          do i=1,mx
              if (useBouss(i)) then
                  etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                  psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
              else
                  psi(i+1) = 0.d0
              endif
           enddo
        endif
              
        ! Forward Euler update to momentum q(2,:):
        ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
        q0(2,1:mx) = q0(2,1:mx) + delt*psi(2:mx+1)
          

        !-----------------------
        if (rk_order == 1) then
            q = q0 ! and we are done
              
        else if (rk_order == 2) then

            ! Second stage for 2-stage R-K, solve for psi at midpoint
            ! in time based on q0 computed in first stage:
            
            if (boussEquations == 1) then
                call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! returns solution psi = source term for Madsen
            else if (boussEquations == 2) then
                call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! modify solution psi for source term of SGN:
                do i=1,mx
                    ! note that h,eta were not changed by first stage
                    if (useBouss(i)) then
                        etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                        psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
                    else
                        psi(i+1) = 0.d0
                    endif
                enddo
            endif
              
            ! Second stage is midpoint method dt=delt*2:
            ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
            q(2,1:mx) = q(2,1:mx) + 2.d0*delt*psi(2:mx+1)

        endif ! rk_order == 2

    endif ! end of Bouss terms


end subroutine src1
