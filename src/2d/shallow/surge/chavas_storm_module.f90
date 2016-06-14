module chavas_storm_module
  implicit none
  save
  
  type chavas_storm_type
    integer :: num_casts
    real(kind = 8), allocatable :: track(:,:)

    real(kind = 8) :: ambient_pressure = 101.3d3
    real(kind = 8), allocatable :: max_wind_radius(:)
    real(kind = 8), allocatable :: max_wind_speed(:)
    real(kind = 8), allocatable :: central_pressure(:)
    real(kind = 8), allocatable :: velocity(:,:)
  end type chavas_storm_type

contains

  real(kind=8) function get_chi(v)
    !chi is a free parameter given by the ratio of C_d, a drag coefficient and
    ! W_cool, which is related to the ekman transfer rate.  
    ! W_cool is approximately 2 s/mm, while C_d varies as a function of velocity
    real(kind=8), intent(in) :: v
    if (v < 6) then
      get_chi = .308
    else if (v < 35.4) then
      get_chi = .0295*v + .1307
    else
      get_chi = 1.2
    end if
  end function

  real(kind=8) function get_alpha(v_max)
    real(kind=8), intent(in) :: v_max
    !Uses the quadratic fit given in Chavas et al. to predict C_d/C_k.
    get_alpha = 0.00055*v_max**2 - 0.0259*v_max + 0.763
  end function

  real(kind=8) function coriolis_param(latitude)
    real(kind=8), intent(in)    :: latitude
    real(kind=8)                :: omega
    !This might be redundant with another function in geoclaw.  Determines
    ! the value of the coriolis parameter f.
    omega = 7.2921e-5
    coriolis_param = 2.*omega*sin(2*3.14159*latitude/360.)
    return
  end function coriolis_param

  subroutine solve_r_max(f, r_ref, v_ref, v_max, r_max)
    real(kind=8), intent(in)    :: f, r_ref, v_ref, v_max
    real(kind=8)                :: r_max, alpha
    real(kind=8)                :: lower=0, upper, mid, check
    real, parameter     :: RMAX_TOLERANCE=.01
    ! This function takes a value of f, a reference velocity and the radius at
    ! which the velocity occurs, and the max velocity, and outputs the radius
    ! of maximum velocity.
    alpha = get_alpha(v_max)
    upper = r_ref
    do
       !This is a golden section method to find the zero of 'check'.
       mid = (upper + lower)/2
       !check is the difference between the left hand and right hand
       !side of the equation given for the momentum of the inner region of
       ! a hurricane given in Chavas et al.
       check = ((r_ref*v_ref + 0.5*f*r_ref**2)/(mid*v_max + 0.5*f*mid**2))**(2 - alpha)
       check = check - (2*((r_ref/mid)**2)/(2-alpha+alpha*(r_ref/mid)**2))**(1/(2-alpha))
       if (check < 0) then
          upper = mid
       else
          lower = mid
       end if
       if( abs(upper - lower) < RMAX_TOLERANCE ) then
          exit
       end if
    end do
    r_max = mid
  end subroutine solve_r_max

  subroutine integ_m_out(f, r_0, r_a, res, m_out)
    real(kind=8), intent(in)    :: f, r_0, r_a
    integer, intent(in) :: res
    real(kind=8), dimension(1:res)  :: m_out
    integer             :: i
    real(kind=8)                :: a,b,c,r_p,r_m,chi,v
    !Given outer model parameters, computes a radial wind profile for the outer
    !region of a hurricane.
    !The equation being integrated is \partial_r M = (RV)**2/(r_0**2 - r**2)
    !Starting from r_0 and going down in r.
    m_out(res) = .5*f*r_0**2
    r_m = r_0 - (r_0 - r_a)/(res - 1)
    m_out(res - 1) = .5*f*r_m**2 + r_m*f*(r_0 - r_m)
    do i = res - 2, 1, -1
       r_p = r_m
       r_m = r_p - (r_0 - r_a)/(res-1)
       v = (m_out(i+1) - .5*f*r_p**2)/r_p
       chi = get_chi(v)
       m_out(i) = m_out(i+1) + chi*(r_p - r_m)*(r_p*v)**2/(r_p**2 - r_0**2)
    end do
  end subroutine integ_m_out

  subroutine solve_r_0(f, r_a, v_a, res, r_0, r_guess, m_out)
    real(kind=8), intent(in)    :: f, r_a, v_a
    integer, intent(in) :: res
    real(kind=8)                :: r_0, m_diff, r_guess, r_step
    real(kind=8), dimension(1:res)  :: m_out
    real, PARAMETER     :: R_TOLERANCE = .1
    integer :: STEP_FLAG = 1
    !Uses a golden section method to find the value of r_0, the radius
    !at which the hurricane winds drop to 0.

    r_0 = r_guess
    r_step = r_guess
    do
       call integ_m_out(f, r_0, r_a, res, m_out)
       !Check for continuity with the pre-determined inner model.
       if(m_out(1) - r_a*v_a - .5*f*r_a**2 < 0) then
          r_step = r_step/STEP_FLAG
          r_0 = r_0 + r_step
       else
          r_step = r_step/2
          r_0 = r_0 - r_step
          STEP_FLAG = 2
       end if
       if(r_step < R_TOLERANCE) then
          r_guess = r_0
          exit
       end if
    end do
  end subroutine solve_r_0

  real(kind=8) function inner_derivative(f, r_a, r_max, v_max)
    real(kind=8), intent(in)   :: f, r_a, r_max, v_max
    real(kind=8)               :: deriv, alpha
    !differentiates the radial velocity profile of the inner model.
    !very messy, but analytically solvable.
    alpha = get_alpha(v_max)
    deriv = 4.*r_a*(2. - alpha)/(r_max**2.)
    deriv = deriv/(2. - alpha + alpha*(r_a/r_max)**2.)**2.
    deriv = deriv*(r_max*v_max + .5*f*r_max**2.)/(2. - alpha)
    deriv = deriv*( (2.*(r_a/r_max)**2.)**((alpha-1.)/(2.-alpha)) )
    deriv = deriv/( (2. - alpha + alpha*(r_a/r_max)**2.)**((alpha - 1.)/(2. - alpha)) )

    inner_derivative = deriv
    return
  end function inner_derivative

  real(kind=8) function eval_v_in(f, r_max, v_max, r_a)
    real(kind=8), intent(in)   :: f, r_max, v_max, r_a
    real(kind=8)               :: v_a, alpha
    !Evaluate the wind velocity for the inner model of the hurricane.
    !See Chavas et al. for information on the model.
    alpha = get_alpha(v_max)
    v_a = (2.*(r_a/r_max)**2/(2. - alpha + alpha*(r_a/r_max)**2))**(1.0/(2-alpha))
    v_a = v_a*(.5*f*r_max**2 + r_max*v_max)
    v_a = v_a - .5*f*r_a**2
    v_a = v_a/r_a
    eval_v_in = v_a
  end function eval_v_in

  subroutine solve_hurricane(f, r_max, v_max, res, r_0, r_a)
    real(kind=8), intent(in)   :: f, r_max, v_max
    integer, intent(in):: res
    real(kind=8)       :: r_0, r_a, v_a, r_step, r_guess, chi
    real(kind=8)       :: slope_diff
    real(kind=8), dimension(1:res) :: m_out
    real, PARAMETER    :: R_TOLERANCE = 1
    integer :: STEP_FLAG=1
    
    !Determines the appropriate model parameters r_a, v_a, and r_0 for the
    !hurricane wind profile.  Uses a golden section method to refine guesses
    !of r_a, which are used to determine r_0 in subroutine solve_r_0.

    r_a = 2*r_max
    r_guess = 5*r_a
    r_step = r_a
    do
       v_a = eval_v_in(f, r_max, v_max, r_a)
       if (v_a < 0) then
         r_step = r_step/2
         r_a = r_a - r_step
         cycle
       end if
       chi = get_chi(v_a)
       call solve_r_0(f, r_a, v_a, res, r_0, r_guess, m_out)
       !checks for continuity in the first derivative of the momentum.
       slope_diff = inner_derivative(f, r_a, r_max, v_max) &
            - chi*((r_a*v_a)**2)/(r_0**2 - r_a**2)
       if (slope_diff < 0) then
          r_step = r_step/2
          r_a = r_a - r_step
          STEP_FLAG = 2
       else
          r_step = r_step/STEP_FLAG
          r_a = r_a + r_step
       end if
       if (r_step < R_TOLERANCE) then
          exit
       end if
    end do
  end subroutine solve_hurricane

  subroutine set_chavas_storm(storm_data_path, storm, log_unit)
    ! Mostly identical to set_holland_storm, but adapted for use with 
    ! hurdat2 format and changed slightly to account for the difference
    ! in relevant parameters between the two models.

    use amr_module, only: t0, rinfinity
    use holland_storm_module, only: date_to_seconds
    use geoclaw_module, only: spherical_distance, coordinate_system

    character(len = *), optional :: storm_data_path
    type(chavas_storm_type), intent(in out) :: storm
    integer, intent(in) :: log_unit

    integer :: i,j,k

    integer, PARAMETER :: data_file = 701
    integer :: io_status
    real(kind = 8) :: forecast_time,last_time,x(2),y(2),ds,dt,dx,dy,theta
    real(kind = 8) :: lat, lon
    integer :: year, month, day, hour, forecast, max_wind_speed
    integer :: central_wind_pressure, last_storm_index, num_casts
    integer, dimension(1:12) :: ref_wind_extent
    real(kind=8) :: ref_wind_radius, ref_wind_speed
    character(len = 1) :: direction(2)

    character(len = 4) :: file_format = "NOAA"
    !NOAA, in this case, is hurdat2 format.  This is not the same as the
    !format used by the holland_storm module.
    character(len = *), parameter :: NOAA_FORMAT = "(i4, i2, i2,2x, i2, 11x,"//&
                           "f4.0, a1, 2x, f5.0, a1, 2x, i3, 2x, i4, 12(3x, i4))"

    if (coordinate_system /= 2) then
            stop "Chavas storm type only works on lat-long coordinates."
        endif
  
    
      ! Open data file
    if (present(storm_data_path)) then
      print *,'Reading storm date file ',storm_data_path
      open(unit=data_file,file=storm_data_path,status='old', &
           action='read',iostat=io_status)
    else
      print *,'Reading storm date file ./storm.data'
      open(unit=data_file,file="./storm.data",status='old', &
           action='read',iostat=io_status)
    endif
    if (io_status /= 0) then
      print "(a,i2)", "Error opening storm data file. status = ", io_status
      stop 
    endif    
    ! Count number of data lines
    num_casts = 0
    last_time = -rinfinity
    do
      if (file_format == "NOAA") then
      read(data_file,fmt=NOAA_FORMAT,iostat=io_status) year,month,day, &
              hour,lat,direction(2),lon,direction(1), &
              max_wind_speed,central_wind_pressure,ref_wind_extent(1:12)
      else 
        print *, "Error: Unrecognised storm data file format."
        stop
      end if
      if (io_status /= 0) exit
    
      forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
      if (abs(forecast_time - last_time) >= 1.8d3) then
        num_casts = num_casts + 1
      end if
      last_time = forecast_time
    end do
    rewind(data_file)
    
    write(log_unit, "('Forecasts = ', i3)") num_casts
    
    allocate(storm%track(3,num_casts))
    allocate(storm%max_wind_speed(num_casts))
    allocate(storm%max_wind_radius(num_casts))
    allocate(storm%central_pressure(num_casts))
    i = 1
    do while(i < num_casts)
      if (file_format == "NOAA") then
      read(data_file,fmt=NOAA_FORMAT,iostat=io_status) year,month,day, &
              hour,lat,direction(2),lon,direction(1), &
              max_wind_speed,central_wind_pressure,ref_wind_extent(1:12)
      else 
        print *, "Error: Unrecognised storm data file format."
        stop
      end if
      
      forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
      if (abs(forecast_time - last_time) < 1.8d3) then
        cycle
      end if
      i = i+1
      last_time = forecast_time
      storm%track(1,i) = forecast_time
      if (direction(1) == "E") then
        storm%track(2,i) = real(lon,kind=8) 
      else
        storm%track(2,i) = -real(lon,kind=8) 
      end if
      if (direction(2) == "N") then
        storm%track(3,i) = real(lat,kind=8) 
      else
        storm%track(3,i) = -real(lat,kind=8) 
      end if

      storm%max_wind_speed(i) = real(max_wind_speed, kind=8)*0.51444444
      storm%central_pressure(i) = real(central_wind_pressure, kind=8)*100
      j = 0
      !Check to find the highest velocity reference wind given.
      do while (j < 12)
        if (ref_wind_extent(j + 1) == 0) exit
        j = j+1
      end do
      j = j/4
      if (j /= 0) j = j - 1
      ref_wind_radius = 0
      !Find the root mean square of the radius, which should approximate a
      !radial average.
      do k = 1,4
        ref_wind_radius = ref_wind_radius + 0.25*ref_wind_extent(4*j + k)**2
      end do
      ref_wind_radius = 1852*sqrt(ref_wind_radius)
      !Hurdat2 format gives reference wind radii with velocities of 34, 50, and
      !64 knots.  This case statement checks which of those I'm using.
      select case(j)
        case(0) 
          ref_wind_speed = 34*.51444444
        case(1)
          ref_wind_speed = 50*.51444444
        case(2)
          ref_wind_speed = 64*.51444444
      end select
      !Determine radius of max velocity from the reference measurements.
      if (ref_wind_radius /= 0) then
        call solve_r_max(coriolis_param(storm%track(3,i)), ref_wind_radius, &
                         ref_wind_speed, storm%max_wind_speed(i), &
                         storm%max_wind_radius(i))
      else
        !Use regression suggested in Kossin et al
        storm%max_wind_radius(i) = -0.4813*max_wind_speed &
                                   + 75.5154 + 0.6992*storm%track(3,i)
        storm%max_wind_radius(i) = storm%max_wind_radius(i)*1000
      end if
    end do  
    ! Calculate storm speed 
    allocate(storm%velocity(2,num_casts))
    do i=1,num_casts - 1
      ! Calculate velocity based on great circle distance between

      ! locations of storm
      x = storm%track(2:3,i)
      y = storm%track(2:3,i+1)

      dt = storm%track(1,i + 1) - storm%track(1,i)

      ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
                              y(1), 0.5d0 * (x(2) + y(2)))
      storm%velocity(1,i) = sign(ds / dt,y(1) - x(1))

            
      ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
                              0.5d0 * (x(1) + y(1)), y(2))
      storm%velocity(2,i) = sign(ds / dt,y(2) - x(2))
    end do
    ! Use last approximation for velocity point going forward
    storm%velocity(:,num_casts) = storm%velocity(:,num_casts - 1)

    ! Record number of casts
    storm%num_casts = num_casts

    if (t0 < storm%track(1,1)) then
      print *,t0,storm%track(1,1)
      stop "Start time is before first forecast time."
    endif

    ! This is used to speed up searching for correct storm data
    last_storm_index = 2
    last_storm_index = get_index(t0,storm)
    if (last_storm_index == -1) then
      print *,"Forecast not found for time ",t0,'.'
      stop
    endif
    ! Log everything to the surge log file
    write(log_unit,*) ""
    write(log_unit,*) "Storm Track and Strength"
    write(log_unit,*) ""
    do i=1,storm%num_casts
      write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
                                  (storm%velocity(k,i),k=1,2), &
                                   storm%max_wind_speed(i),  &
                                   storm%central_pressure(i)
    enddo
  end subroutine set_chavas_storm

  function get_index(t, storm) result(index)
    
    real(kind=8), intent(in) :: t
    type(chavas_storm_type), intent(in) :: storm
    integer :: index
    integer :: i
    do i = storm%num_casts, 1, -1
       if (storm%track(1, i) < t) then
          index = i
          exit
       else if (i == 1) then
         index = 0 !indicates that t is before the earliest forecast
         exit
       end if
    end do
    return
  end function

  function chavas_storm_location(t, storm) result(location)

    real(kind = 8), intent(in) :: t
    type(chavas_storm_type), intent(in) :: storm
    integer :: index
    real(kind=8) :: t_factor
    real(kind=8) :: location(2)
    !find the location of the eye of the storm by interpolating track data.
    index = get_index(t, storm)
    if (index == 0) then
       !If extrapolating to before the first track, use the first given position
       location(1) = storm%track(2, 1)
       location(2) = storm%track(3, 1)
    else if (index == storm%num_casts) then
      !If extrapolating after the end of the track, use end velocity to guess
      !location
      t_factor = t - storm%track(1, index)
      location(1) = storm%track(2,index) + t_factor*storm%velocity(1, index - 1)
      location(2) = storm%track(3,index) + t_factor*storm%velocity(2, index - 1)
    else
      !otherwise, linearly interpolate between two time adjacent tracks.
      t_factor = t - storm%track(1, index)
      t_factor = t_factor/(storm%track(1, index + 1) - storm%track(1, index))
      location(1) = storm%track(2,index)*(1 - t_factor) + &
                    t_factor*storm%track(2,index + 1)
      location(2) = storm%track(3, index)*(1 - t_factor) + &
                    t_factor * storm%track(3, index+1)
    end if
  end function chavas_storm_location

  real(kind=8) function chavas_storm_direction(t, storm) result(theta)
    
    implicit none
    real(kind=8), intent(in) :: t
    type(chavas_storm_type), intent(in) :: storm
    integer :: index
    !Determine the bearing of the storm movement.  0 is east.
    index = get_index(t, storm)
    theta = atan2(storm%velocity(2, index), storm%velocity(1, index))
  end function chavas_storm_direction

  subroutine set_chavas_storm_fields(maux,mbc,mx,my, &
                                    xlower,ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, storm)
    
    use geoclaw_module, only: coriolis, deg2rad, spherical_distance
    ! Time of the wind field requested
    integer, intent(in) :: maux,mbc,mx,my
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

    ! Storm description, need in out here since we may update the storm
    ! if at next time point
    type(chavas_storm_type), intent(in out) :: storm

    ! Array storing wind and presure field
    integer, intent(in) :: wind_index, pressure_index
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer, PARAMETER :: res=2000
    real(kind=8), save :: v_vec(1:res), p_vec(1:res)
    real(kind=8), save :: m_out_vec(1:res)
    real(kind=8), save :: last_time=-1
    real(kind=8), save :: dr, r_a, r_0
    integer, save :: in_res, out_res
    ! Local storage
    real(kind=8) :: x, y, r, theta, sloc(2), v, next_loc(2)
    real(kind=8) :: f, r_max, v_max, Pc, Pa, dp, wind, tv(2)
    integer :: i,j, r_index
    real(kind=8) :: r_factor, t_factor
    !output the pressure and velocity fields.  Determines the radial wind
    !profile only if the called at a new time.

    !First, interpolate all of the relevant tracked storm parameters.
    i = get_index(t, storm)
    tv = storm%velocity(:, i)
    Pa = storm%ambient_pressure
    sloc = chavas_storm_location(t, storm)
    f = coriolis_param(sloc(2))
    if(i < storm%num_casts) then
      t_factor = (t - storm%track(1,i))/(storm%track(1,i+1) - storm%track(1,i))
      Pc = (1. - t_factor)*storm%central_pressure(i) + &
           t_factor*storm%central_pressure(i+1)
      v_max = (1. - t_factor)*storm%max_wind_speed(i) + &
              t_factor*storm%max_wind_speed(i+1)
      r_max = (1. - t_factor)*storm%max_wind_radius(i) + &
           t_factor*storm%max_wind_radius(i+1)
    else
      Pc = storm%central_pressure(i)
      v_max = storm%max_wind_speed(i)
      r_max = storm%max_wind_radius(i)
    end if
    if (last_time /= t) then
       !Determine the wind profile, but only if the profile isn't already
       !saved. In v_vec and p_vec
       last_time = t
       call solve_hurricane(f, r_max, v_max, res, r_0, r_a)
       dr = real(r_0)/(res - 1)
       !out_res = number of points in outer model/in_res = points in inner model
       out_res = ceiling( real((res - 1)*(r_0 - r_a))/r_0 )
       in_res = res - out_res
       !Determine angular momentum profile for outer model
       call integ_m_out(f, r_0, r_0 - out_res*dr, out_res, m_out_vec)
       do j= 1, out_res
          !convert angular momentum to azimuthal velocity.
          r = r_0 - (j-1)*r_0/res
          m_out_vec(out_res + 1 - j) = (m_out_vec(out_res + 1 - j)- .5*f*r**2)/r
       end do
       p_vec(1) = 0
       v_vec(1) = 0
       do j=2, res
          !combine inner and outer model into single radial wind profile.
          if(j <= in_res) then
             r = (j - 1)*dr
             v_vec(j) = eval_v_in(f, r_max, v_max, r)
          else
             v_vec(j) = m_out_vec(j - in_res)
          end if
          !pressure is determined using the gradient wind equation.
          p_vec(j) = p_vec(j - 1) + (v_vec(j)**2)/r + f*v_vec(j)
       end do
       do j=1, res
          !normalize pressure to match measurements.
          p_vec(j) = (Pa - Pc)*p_vec(j)/p_vec(res) + Pc
       end do
    end if
    do j=1 - mbc,  my + mbc
      !Assigns a clockwise, radially symmetric wind field by piece-wise linear
      !interpolation of v_vec.  Likewise with p_vec.
      y = ylower + (j - 0.5d0) * dy
      do i = 1 - mbc, mx + mbc
        x = xlower + (i - 0.5d0) * dx
        r = spherical_distance(x, y, sloc(1), sloc(2))
        r_index = floor(r/dr)
        r_factor = (r - r_index*dr)/dr
        r_index = r_index + 1
        theta = atan2((y - sloc(2)) * DEG2RAD,(x - sloc(1)) * DEG2RAD)
        if(r_index < res) then
          aux(pressure_index, i, j) = (1-r_factor)*p_vec(r_index) &
                                    + r_factor*p_vec(r_index + 1)
          aux(wind_index,i,j) = -0.9*( (1-r_factor)*v_vec(r_index) + &
                                   r_factor*v_vec(r_index + 1) )*sin(theta)
          aux(wind_index+1,i,j) = 0.9*( (1-r_factor)*v_vec(r_index)  +&
                                   r_factor*v_vec(r_index + 1) )*cos(theta)
        else
          aux(pressure_index,i,j) = Pa
          aux(wind_index,i,j) = 0
          aux(wind_index+1,i,j) = 0
        end if
      end do
    end do
  end subroutine set_chavas_storm_fields

end module
