subroutine rpn2(ixy, maxm, num_eqn, num_waves, num_aux, num_ghost, num_cells, &
                ql, qr, auxl, auxr, fwave, s, amdq, apdq)

    use, intrinsic :: iso_fortran_env, only: real64

    use geoclaw_module, only: g => grav, dry_tolerance, rho
    use storm_module, only: pressure_index

    implicit none 
    
    integer, parameter :: D = real64

    ! Arguments
    integer, intent(in) :: ixy, maxm, num_eqn, num_waves, num_ghost, num_aux, num_cells
    real(kind=D), intent(in) :: ql(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(in) :: qr(num_eqn, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(in) :: auxl(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(in) :: auxr(num_aux, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(out) :: s(num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(out) :: amdq(num_eqn,1-num_ghost:maxm+num_ghost)
    real(kind=D), intent(out) :: apdq(num_eqn,1-num_ghost:maxm+num_ghost)
    
    ! Locals
    integer :: i, k, normal_index, transverse_index
    real(kind=D) :: hl, ul, vl, hr, ur, vr, hbar, uhat, chat, db, dp
    real(kind=D) :: phil, phir, dry_state_l, dry_state_r
    real(kind=D) :: R(3,3)
    real(kind=D) :: delta(3), beta(3)
    
    
    ! Determine normal and tangential directions
    if (ixy == 1) then
        normal_index = 2
        transverse_index = 3
    else
        normal_index = 3
        transverse_index = 2
    end if

    amdq = 0.0_D
    apdq = 0.0_D

    ! Primary loop over each cell
    do i = 2 - num_ghost, num_cells + num_ghost
        
        ! Check for dry states - need merge here to convert to float
        dry_state_l = merge(0.0_D, 1.0_D, qr(1, i - 1) < dry_tolerance)
        dry_state_r = merge(0.0_D, 1.0_D, ql(1, i) < dry_tolerance)

        ! Note that for the states below u is always the normal velocity and
        ! v is always the tangential velocity

        ! Left states
        hl = qr(1, i - 1) * dry_state_l
        ul = qr(normal_index, i - 1) / qr(1, i - 1) * dry_state_l
        vl = qr(transverse_index, i - 1) / qr(1, i - 1) * dry_state_l
        phil = (0.5_D * g * hl**2 + hl * ul**2) * dry_state_l

        ! Forcing
        db = (auxl(1, i) - auxr(1, i - 1))
        dp = 0.d0
        ! dp = (auxl(pressure_index, i) - auxr(pressure_index, i - 1))
    
        ! Right states
        hr = ql(1, i) * dry_state_r
        ur = ql(normal_index, i) / ql(1, i) * dry_state_r
        vr = ql(transverse_index, i) / ql(1, i) * dry_state_r
        phir = (0.5_D * g * hr**2 + hr * ur**2) * dry_state_r

        ! Roe average states (Roe's linearization)
        hbar = 0.5_D * (hr + hl)
        uhat = (sqrt(hr) * ur + sqrt(hl) * ul) / (sqrt(hr) + sqrt(hl))
        chat = sqrt(g * hbar)
    
        ! Flux differences
        delta(1) = hr * ur - hl * ul
        delta(2) = phir - phil + g * hbar * db + hbar * dp / rho(1)
        delta(3) = hr * ur * vr - hl * ul * vl
    
        ! Wave speeds
        s(1, i) = min(uhat - chat, ul - sqrt(g * hl))
        s(3, i) = max(uhat + chat, ur + sqrt(g * hr))
        s(2, i) = 0.5_D * (s(1, i) + s(3, i))
        
        ! Right eigenvectors (columns)
        ! could possibly use vhat instead of vl and vr
        R(1, 1) = 1.0_D
        R(normal_index, 1) = s(1, i)
        R(transverse_index, 1) = vl
        
        R(1, 2) = 0.0_D
        R(normal_index, 2) = 0.0
        R(transverse_index, 2) = 1.0
        
        R(1, 3) = 1.0_D
        R(normal_index, 3) = s(3, i)
        R(transverse_index, 3) = vr
        
        ! Wave strengths
        beta(1) = (s(3, i) * delta(1) - delta(2)) / (s(3, i) - s(1, i))
        beta(3) = (delta(2) - s(1, i) * delta(1)) / (s(3, i) - s(1, i))
        beta(2) = delta(3) - beta(1) * vl - beta(3) * vr

        ! f-waves
        do k = 1, num_waves
            fwave(:, k, i) = beta(k) * R(:, k)
        enddo
    
        ! Fluctuations
        do k=1, num_waves
            amdq(:, i) = amdq(:, i) + merge(fwave(:, k, i), 0.0_D,            &
                                                        s(k, i) < -1e-14)
            apdq(:, i) = apdq(:, i) + merge(fwave(:, k, i), 0.0_D,            &
                                                        s(k, i) >  1e-14)

            amdq(:, i) = amdq(:, i) + merge(0.5_D * fwave(:, k, i), 0.0_D,    &
                                                          -1e-14_D < s(k, i)  & 
                                                    .and. s(k, i) < 1e-14)
            apdq(:, i) = apdq(:, i) + merge(0.5_D * fwave(:, k, i), 0.0_D,    &
                                                          -1e-14_D < s(k, i)  &
                                                    .and. s(k, i) < 1e-14)
        enddo

    enddo ! End of main loop

end subroutine rpn2
