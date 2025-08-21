subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: spherical_distance, coordinate_system
    use geoclaw_module, only: RAD2DEG, pi, dry_tolerance, DEG2RAD
    use geoclaw_module, only: rho_air, rho
    use geoclaw_module, only: earth_radius, sphere_source

    use geoclaw_module, only: speed_limit  ! used if no friction
      
    use storm_module, only: wind_forcing, pressure_forcing, wind_drag
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: variable_friction, friction_index

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: u, v, hu0, hv0
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2)
    real(kind=8) :: Ddt, sloc(2)
    real(kind=8) :: tanyR, huv, huu, hvv
    real(kind=8) :: speed, sratio, xs, ys

    ! Algorithm parameters

    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! ----------------------------------------------------------------
    ! Spherical geometry source term(s)
    !
    ! These should be included for shallow water on the sphere, 
    ! at least in the mass term, but were only added in v5.9.0 as an option.
    ! rundata.geo_data.sphere_source can now be set in setrun.py
    ! Set sphere_source = 0 to omit source terms for backward compatibility
    ! sphere_source = 1 should become the default?

    if ((coordinate_system == 2) .and. (sphere_source > 0)) then
        ! add in spherical source term in mass equation 
        ! if sphere_source in [1,2],
        ! and also in momentum equations if sphere_source == 2
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            tanyR = tan(y*DEG2RAD) / earth_radius
            do i=1,mx
                if (q(1,i,j) > dry_tolerance) then
                    ! source term in mass equation:
                    q(1,i,j) = q(1,i,j) + dt * tanyR * q(3,i,j)

                    if (sphere_source == 2) then
                        ! Momentum source terms that drop out if linearized:
                        ! These seem to have very little effect for
                        ! practical problems
                        huv = q(2,i,j)*q(3,i,j)/q(1,i,j)
                        huu = q(2,i,j)*q(2,i,j)/q(1,i,j)
                        hvv = q(3,i,j)*q(3,i,j)/q(1,i,j)
                        q(2,i,j) = q(2,i,j) + dt * tanyR * 2.d0*huv
                        q(3,i,j) = q(3,i,j) + dt * tanyR * (hvv - huu)
                    endif
                endif
            enddo
        enddo
    endif
                
    ! ----------------------------------------------------------------                
    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        if (.not.variable_friction) then
                            do nman = num_manning, 1, -1
                                if (aux(1,i,j) .lt. manning_break(nman)) then
                                    coeff = manning_coefficient(nman)
                                endif
                            enddo
                        else
                            coeff = aux(friction_index,i,j)
                        endif
                        
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                              * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                        dgamma = 1.d0 + dt * gamma
                        q(2, i, j) = q(2, i, j) / dgamma
                        q(3, i, j) = q(3, i, j) / dgamma
                    endif
                endif
            enddo
        enddo
    ! End of friction source term

    else  

      ! if not friction_forcing: 

      ! Use speed_limit to scale back unphysical speeds that might
      ! arise from computing hu / h when h is very small.  This should
      ! only be needed in this routine if friction is not being used or
      ! the on-shore manning_coefficient is very small. 
      ! (In the latter case, the user might need to tweak this by hand.)

      do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (q(1,i,j) > 0.d0) then
                speed = sqrt((q(2,i,j)**2 + q(3,i,j)**2)) / q(1,i,j)
                if (speed > speed_limit) then
                    sratio = speed_limit / speed
                    q(2,i,j) = q(2,i,j) * sratio
                    q(3,i,j) = q(3,i,j) * sratio

                    if (.false.) then
                        ! write out info useful for investigating topo:
                        ! might give LOTS of output if no friction!
                        xs = xlower + (i-0.5d0)*dx
                        ys = ylower + (j-0.5d0)*dy
                        write(6,604) t, i,j, mx,my, dx
                        write(6,603) speed,q(1,i,j),aux(1,i,j),xs,ys

 604                    format('src2 at t =',f10.2, '  i,j,mx,my:',4i4, &
                               '  dx = ',f10.7)
 603                    format('     reset s =',f9.2,'  h=',e11.3, ' B=',f8.2,&
                               '  x,y = ', f11.6,',',f10.6)
                    endif
               endif
            endif
        enddo
      enddo

    ! end of speed limiting
    endif  ! end of else (i.e., if not friction_forcing)


    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Need storm location and direction for sector based wind drag
        sloc = storm_location(t)
        theta = storm_direction(t)
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            do i=1,mx
                xc = xlower + (i - 0.5d0) * dx
                if (q(1,i,j) > dry_tolerance) then
                    psi = atan2(yc - sloc(2), xc - sloc(1))
                    if (theta > psi) then
                        phi = (2.d0 * pi - theta + psi) * RAD2DEG
                    else
                        phi = (psi - theta) * RAD2DEG 
                    endif
                    wind_speed = sqrt(aux(wind_index,i,j)**2        &
                                    + aux(wind_index+1,i,j)**2)
                    tau = wind_drag(wind_speed, phi) * rho_air * wind_speed / rho(1)
                    q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                    q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
    ! Handled in Riemann solver
    ! if (pressure_forcing) then
    !     do j=1,my  
    !         ym = ylower + (j - 1.d0) * dy
    !         yc = ylower + (j - 0.5d0) * dy
    !         yp = ylower + j * dy
    !         do i=1,mx  
    !             xm = xlower + (i - 1.d0) * dx
    !             xc = xlower + (i - 0.5d0) * dx
    !             xp = xlower + i * dx
                
    !             if (coordinate_system == 2) then
    !                 ! Convert distance in lat-long to meters
    !                 dx_meters = spherical_distance(xp,yc,xm,yc)
    !                 dy_meters = spherical_distance(xc,yp,xc,ym)
    !             else
    !                 dx_meters = dx
    !                 dy_meters = dy
    !             endif

    !             ! Extract depths
    !             h = q(1,i,j)

    !             ! Calculate gradient of Pressure
    !             P_gradient(1) = (aux(pressure_index,i+1,j) &
    !                            - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
    !             P_gradient(2) = (aux(pressure_index,i,j+1) &
    !                            - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

    !                 ! Modify momentum in each layer
    !             if (h > dry_tolerance) then
    !                 q(2, i, j) = q(2, i, j) - dt * h * P_gradient(1) / rho(1)
    !                 q(3, i, j) = q(3, i, j) - dt * h * P_gradient(2) / rho(1)
    !             end if
    !         enddo
    !     enddo
    ! endif

end subroutine src2
