! Contains a number of useful functions
module utility_module

    implicit none

    ! ISO Time Format String
    character(len=*), parameter :: ISO_time_format = "(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)"

    ! String manipulation
    character( * ), private, parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character( * ), private, parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

contains

   ! ==========================================================================
   ! Check for netcdf file errors when loading data, only active if the
   ! NETCDF FFLAGS are in the Makefile
   ! ==========================================================================
    subroutine check_netcdf_error(ios)
#ifdef NETCDF
        use netcdf
#endif
        implicit none

        integer, intent(in) :: ios
#ifdef NETCDF
        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if
#else
        print *, "GeoClaw was not compiled with NetCDF support."
        stop
#endif
    end subroutine check_netcdf_error


    ! Returns number of arguments in list, assumes arbitrary number of white 
    ! space as delimiter
    integer function get_value_count(line) result(count)

        implicit none

        character(len=*), intent(in) :: line

        character(len=200) :: search_buffer
        integer :: i
        logical :: found

        count = 1

        search_buffer = trim(line(1:index(line,"=:") - 2))

        do i=1,len_trim(search_buffer)
            if (search_buffer(i:i) == " ") then
                if (.not. found) then
                    found = .true.
                    count = count + 1
                endif
            else
                found = .false.
            endif
        enddo

    end function get_value_count

    ! Converts seconds to days truncating at the second decimal place
    real(kind=8) pure function convert2days(seconds) result(days)
        
        implicit none
        real(kind=8), intent(in) :: seconds

        days = real(int(seconds * 1.d2 / 8.64d4) / 1.d2, kind=8)

    end function convert2days


    !======================================
    subroutine parse_values(str, n, values)
    !======================================

    ! Take the input string `str` and parse it to extract any numerical values,
    ! ignoring any character strings. 
    ! Returns `n` values in array `values`. 
    ! Assumes n <= 10.
    ! If you expect value(i) to be an integer, evaluate as nint(value(i)).

    implicit none

    character(len=*), intent(in) :: str
    integer, intent(out) :: n
    real(kind=8), intent(out) :: values(16)

    integer :: pos2,nw,i,e
    character(len=80) :: word(16), str2
    real(kind=8) :: x

    ! First break into words / tokens based on white space.  
    ! Each might be character or numerical:

    nw = 0
    str2 = trim(adjustl(str))
    do while (len(trim(adjustl(str2))) > 0) 
        pos2 = index(str2, " ")
        
        if (pos2 == 0) then
           nw = nw + 1
           word(nw) = trim(adjustl(str2))
           exit
           endif

        nw = nw + 1
        if (nw == 17) then
            write(6,*) '*** too many words on line, str = '
            write(6,*) str
            stop
            endif
        word(nw) = trim(adjustl(str2(1:pos2-1)))
        str2 = trim(adjustl(str2(pos2+1:)))

        enddo

    ! now extract numerical values:
    n = 0
    do i=1,nw
        read(word(i),*,IOSTAT=e) x
        if (e == 0) then
            ! this token is numerical
            n = n+1
            values(n) = x
            endif
        enddo


    end subroutine parse_values

    ! === Geometry Functions ===================================================
    !
    
    !  point_in_rectangle(p, rect) - Returns whether p is in rect.
    !    real(kind=8) :: p(2), rect(4)
    !    rect is in the format (lower x, lower y, upper x, upper y)
    pure logical function point_in_rectangle(p, rect)
        implicit none
        real(kind=8), intent(in) :: p(2), rect(4)
        point_in_rectangle = p(1) >= rect(1) .and. p(1) <= rect(3) .and. &
                             p(2) >= rect(2) .and. p(2) <= rect(4)
    end function point_in_rectangle

    !  rects_intersect(rect1, rect2) - Returns whether rectangles intersect
    !    real(kind=8) :: rect1(4), rect2(4)
    !    rects are in the format (lower x, lower y, upper x, upper y)
    pure logical function rects_intersect(rect1, rect2)
        implicit none
        real(kind=8), intent(in) :: rect1(4), rect2(4)
        rects_intersect = .not. (rect1(3) < rect2(1) .or.  &
                                 rect1(4) < rect2(2) .or.  &
                                 rect1(1) > rect2(3) .or.  &
                                 rect1(2) > rect2(4))
    end function rects_intersect

    ! ==========================================================================
    ! seconds_from_epoch() Calculates seconds from UNIX epoch (1970) from a
    ! datetime. Returns the total seconds from the epoch includes leap years and
    ! days
    ! ==========================================================================
    pure integer function seconds_from_epoch(time) result(seconds_since_epoch)
        implicit none

        ! year, month, day, hour, [minutes], [seconds] (optional)
        integer, intent(in) :: time(:)
        integer :: days_since_epoch
        integer :: years, months, days, hours, minutes, seconds

        ! Handle possibly missing values
        years = time(1)
        months = time(2)
        days = time(3)
        hours = time(4)
        minutes = 0
        seconds = 0
        if (size(time) >= 5) then
            minutes = time(5)
            if (size(time) == 6) then
                seconds = time(6)
            end if 
        end if

        ! Compute days since epoch (1970-01-01)
        days_since_epoch = date_to_days(years, months, days)        &
                                - date_to_days(1970, 1, 1)

        ! Compute total seconds
        seconds_since_epoch = days_since_epoch * 86400      &
                                + hours * 3600              &
                                + minutes * 60              &
                                + seconds

    end function seconds_from_epoch

    ! Helper function: convert date to days since a fixed point
    pure integer function date_to_days(y, m, d)
        integer, intent(in) :: y, m, d
        integer :: a, y_adj, m_adj

        a = (14 - m) / 12
        y_adj = y + 4800 - a
        m_adj = m + 12 * a - 3

        date_to_days = d + (153 * m_adj + 2) / 5    &
                         + 365 * y_adj              &
                         + y_adj / 4                &
                         - y_adj / 100              &
                         + y_adj / 400              &
                         - 32045                
    end function date_to_days


    ! ==========================================================================
    !  function to_upper(input_string)
    !    Converts *input_string* to upper case and returns that string.
    !    Note that this makes a copy of the given string so does not modify the
    !    original string.
    !
    !   Based on code originally by Paul van Delst, CIMSS/SSEC 18-Oct-1999
    !                               paul.vandelst@ssec.wisc.edu
    !   and originally presented in
    !       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    !       1995 Springer-Verlag, New York.
    ! ==========================================================================
    function to_upper (input_string) result(output_string)

        implicit none
    
        ! Input
        character(len=*), intent(in) :: input_string

        ! Output
        character(len(input_string)) :: output_string

        ! Local
        integer :: i, n

        ! -- copy input string
        output_string = input_string

        ! -- loop over string elements
        do i = 1, len( output_string )

            ! -- find location of letter in lower case constant string
            n = index( LOWER_CASE, output_string( i:i ) )

            ! -- if current substring is a lower case letter, make it upper case
            if ( n /= 0 ) output_string( i:i ) = UPPER_CASE( n:n )

        end do

    end function to_upper

    ! ==========================================================================
    !  function to_lower(input_string)
    !    Converts *input_string* to lower case and returns that string.
    !    Note that this makes a copy of the given string so does not modify the
    !    original string.
    !
    !   Based on code originally by: 
    !       Paul van Delst, CIMSS/SSEC 18-Oct-1999
    !       paul.vandelst@ssec.wisc.edu
    !   and originally presented in
    !       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    !       1995 Springer-Verlag, New York.
    ! ==========================================================================
    function to_lower (input_string) result(output_string)

        implicit none

        ! argument and result
        character(len=*), intent(in) :: input_string
        character(len(input_string)) :: output_string

        ! local variables
        integer :: i, n


        ! copy input string
        output_string = input_string

        ! loop over string elements
        do i = 1, len( output_string )

          ! find location of letter in upper case constant string
          n = index( UPPER_CASE, output_string( i:i ) )

          ! if current substring is an upper case letter, make it lower case
          if ( n /= 0 ) output_string( i:i ) = LOWER_CASE( n:n )

        end do

    end function to_lower


end module utility_module
