subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.

    ! This version checks for very small or negative depths, resets to 0

    ! It also does various other things:
    ! If monitor_fgmax:
    !    keeps track of maximum value of h,s,hss over entire computation,
    !    and arrival time of wave in each cell
    ! If monitor_runup:
    !    keeps track of first and last wet cell each time to track 
    !    inundation and runup as a function of time
    ! If monitor_total_zeta:
    !    keeps track of total mass of water, or actually mass of water
    !    relative to ocean at rest, using
    !       zeta = h in initially dry cells where B < sea_level
    !       zeta = eta = h+B in initially wet cells where B > sea_level

    use geoclaw_module, only: dry_tolerance, sea_level
    use geoclaw_module, only: DEG2RAD, earth_radius, coordinate_system, pi
    use grid_module, only: xp_edge
    use grid_module, only: monitor_fgmax, hmax, smax, hssmax, arrival_time
    use grid_module, only: fgmax_arrival_tol, fgmax_tstart
    use grid_module, only: monitor_total_zeta, iunit_total_zeta_mass
    use grid_module, only: total_zeta_mass_t0, NOTSET
    use grid_module, only: monitor_runup, iunit_runup, runup_tolerance

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)

    !local variables
    integer :: i,ig,j,m,mvars
    real(kind=8) :: x, capa_latitude, zeta, eta
    real(kind=8) :: speed(1-mbc:mx+mbc)
    real(kind=8) :: total_zeta_mass, cell_zeta_mass, dzeta
    real(kind=8) :: x_first_wet,z_first_wet,x_last_wet,z_last_wet
    integer :: i_first_wet, i_last_wet

    ! accumulate total mass in excess of initial sea level
    ! zeta = eta = h+B offshore and zeta = h onshore.  Should be conserved.

    if (coordinate_system == 2) then
        capa_latitude = 2.d0*pi*earth_radius**2 ! to scale to surface area
    endif

    do i=1-mbc,mx+mbc
        ! if cell nearly dry, reset h and hu to zero:
        if (q(1,i)<=dry_tolerance) then
            q(1,i) = max(q(1,i),0.0)
            do m=2,meqn
               q(m,i)=0.d0
            enddo
            speed(i) = 0.d0
        else
            speed(i) = abs(q(2,i)/q(1,i))
        endif
    enddo
         
    if (monitor_fgmax .and. (t >= fgmax_tstart)) then
        do i=1,mx
            ! check arrival time:
            ! using the same criterion as the default in 2d fgmax, 
            ! but for some cases may be better to also look at decrease
            ! in eta and/or water speed
            eta = q(1,i)+aux(1,i)
            if ((q(1,i) > dry_tolerance) &
                    .and. (eta > sea_level + fgmax_arrival_tol) &
                    .and. (arrival_time(i) == NOTSET))  then
                arrival_time(i) = t
            endif

            ! keep track of max depth, speed, momentum flux:
            hmax(i) = dmax1(hmax(i), q(1,i))
            smax(i) = dmax1(smax(i), speed(i))
            hssmax(i) = dmax1(hssmax(i), q(2,i)*speed(i))
        enddo
    endif
    
    if (monitor_total_zeta) then
        total_zeta_mass = 0.d0 
        do i=1,mx
            ! total mass deviation based on zeta:
            if (aux(1,i) < sea_level) then
                zeta = q(1,i) + aux(1,i)  ! = eta = h+B
            else
                zeta = q(1,i)             ! = h
            endif

            cell_zeta_mass = zeta * (xp_edge(i+1) - xp_edge(i))

            if (coordinate_system == -1) then
                ! for total_mass on circle when using radial symmetry:
                cell_zeta_mass = cell_zeta_mass * pi * &
                         (xp_edge(i+1) + xp_edge(i))
            else if (coordinate_system == 2) then
                ! for total mass on full sphere:
                cell_zeta_mass = cell_zeta_mass * (sin(xp_edge(i+1)*DEG2RAD) - &
                         sin(xp_edge(i)*DEG2RAD))/dx * capa_latitude
            endif

            total_zeta_mass = total_zeta_mass + cell_zeta_mass
        enddo
    
        if (total_zeta_mass_t0 == NOTSET) then
            ! this must be first step
            total_zeta_mass_t0 = total_zeta_mass
        endif
        dzeta = total_zeta_mass - total_zeta_mass_t0
 600    format(f16.2, 2e16.6)
        write(iunit_total_zeta_mass,600) t,total_zeta_mass, dzeta
    endif


    if (monitor_runup) then
        i_first_wet = 1
        do while (q(1,i_first_wet) < runup_tolerance)
            i_first_wet = i_first_wet + 1
            if (i_first_wet == mx) exit
        enddo
        x_first_wet = xp_edge(i_first_wet)
        z_first_wet = aux(1,i_first_wet) + q(1,i_first_wet)
    
        i_last_wet = mx
        do while (q(1,i_last_wet) < runup_tolerance)
            i_last_wet = i_last_wet - 1
            if (i_last_wet == 0) exit
        enddo
        x_last_wet = xp_edge(i_last_wet+1)
        z_last_wet = aux(1,i_last_wet) + q(1,i_last_wet)
    
 601    format(5e15.6)
        write(iunit_runup, 601) t,x_first_wet,z_first_wet,x_last_wet,z_last_wet
    endif
        


end subroutine b4step1

