! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.

! This routine is not called in GeoClaw currently due to conservation problems
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
      
    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: rho_air, wind_drag
!     use storm_module, only: pressure_tolerance, wind_tolerance
    use storm_module, only: wind_index, pressure_index

    use friction_module, only: friction_index
      
    use geoclaw_module, only: g => grav, dry_tolerance
    use geoclaw_module, only: coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: omega, coordinate_system

    implicit none

    ! Input
    integer, intent(in) :: meqn, mbc, mx1d, maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(inout) :: q1d(meqn, mx1d), aux1d(maux, mx1d)

    ! Local storage
    integer :: i, m, bottom_index, bottom_layer, layer_index
    logical :: found
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), tau
    real(kind=8) :: wind_speed, theta, P_atmos_x, P_atmos_y

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    real(kind=8), parameter :: rho = 1025.d0

    ! Friction forcing
    if (friction_forcing) then

        do i=1,mx1d

            ! Extract depths
            h = q1d(1,i)
            hu = q1d(2,i)
            hv = q1d(3,i)

            ! If depth is near-zero, set momentum to zero
            if (h < depth_tolerance) then
                q1d(2:3,i) = 0.d0
                cycle 
            endif
            
            ! Apply friction source term only if in shallower water
            if (h <= friction_depth) then
                ! Calculate source term
                gamma = sqrt(hu**2 + hv**2) * (g * aux1d(friction_index,i)**2) &
                                            / h**(7.d0/3.d0)
                dgamma = 1.d0 + dt * gamma
                q1d(2, i) = q1d(2, i) / dgamma
                q1d(3, i) = q1d(3, i) / dgamma
            endif
        enddo
    endif
    
    ! Only lat-long coordinate system supported here right now
    if (coriolis_forcing .and. coordinate_system == 2) then

        do i=1,mx1d
            ! aux(3,:,:) stores the y coordinates multiplied by deg2rad
            fdt = 2.d0 * omega * sin(aux1d(3,i)) * dt

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)
    
            q1d(2,i) = q1d(2,i) * a(1,1) + q1d(3,i) * a(1,2)
            q1d(3,i) = q1d(2,i) * a(2,1) + q1d(3,i) * a(2,2)
        enddo
    endif
    
    ! = Wind Forcing =========================================================
    if (wind_forcing) then
        ! Force only the top layer of water
        ! Assumes that the top-most layer goes dry last
        do i=1,mx1d
            if (q1d(1,i) > dry_tolerance) then
                theta = 0.d0
                wind_speed = sqrt(aux1d(wind_index,i)**2 &
                                + aux1d(wind_index+1,i)**2)
!                 if (wind_speed > wind_tolerance) then
                    tau = wind_drag(wind_speed,theta) * rho_air * wind_speed / rho
                    q1d(2,i) = q1d(2,i) + dt * tau * aux1d(wind_index,i)
                    q1d(3,i) = q1d(3,i) + dt * tau * aux1d(wind_index+1,i)
!                 endif
            endif
        enddo
    endif
    ! ========================================================================
    
    ! == Pressure Forcing ====================================================
    ! Need to add dx and dy to calling signature from qad in order for this to 
    ! work, can probably get it from amr_module but need to know which grid we
    ! are working on
    if (.false.) then
        stop "Not sure how to proceed, need direction and the right dx or dy."
    endif

end subroutine src1d
