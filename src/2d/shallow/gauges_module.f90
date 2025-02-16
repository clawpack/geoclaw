! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine setbestsrc
!     Called each time regridding is done to determine which patch to 
!     use for interpolating to each gauge location.
!   subroutine print_gauges
!     Called each time step for each grid patch.
!     Refactored dumpgauge routine to interpolate for all gauges on patch.
!
!     Note: by default all components of q are printed at each gauge.
!     To print something different or a different precision, modify 
!     format statement 100 and/or the write statement that uses it.
!   
! Note: Updated for Clawpack 5.3.0:
!   - the dumpgauge and setbestsrc subroutines have been moved to this module 
!     and the dumpgauge subroutine has been refactored and renamed print_gauges.
!   - dumpgauge.f must be removed from Makefiles.
!   - setbestsrc uses quicksort to sort gauge numbers and
!     then figures out which gauges will be updated by grid, and stores this
!     information in new module variables mbestg1, mbestg2.
!   - print_gauges no longer uses binary search to locate first gauge handled
!     by a grid.  Instead loop over gauges specified by mbestg1, mbestg2.
!
! Note: Updated for Clawpack 5.4.0
!   - refactor so each gauge writes to its own file, and batches the writes instead of 
!     writing one at a time. This will remove the critical section and should speed up gauges a lot
!   - When array is filled, that gauge will write to file and start over. 
!   - Need to save index so know position in array where left off
!   - At checkpoint times, dump all gauges
!
! Note: Updated for Clawpack 5.4.x
!  - Add gauge formatting capabilities
!
! Changed mjb, May, 2018 hackathon, to remove need for mbestg1 and mbestg2.
! looked at since they depended on maxgr, which is now a variable, not constant.
! also, algorithms didn't scale for O(10^5) grids, and O(100) gauges..

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges

    integer, parameter :: MAX_BUFFER = 1000

    ! Gauge data types
    type gauge_type
        ! Gauge number
        integer :: gauge_num
        integer :: gdata_bytes

        character(len=24) :: file_name      ! for header (and data if 'ascii')
        character(len=24) :: file_name_bin  ! used if file_format='binary'

        ! Location in time and space
        real(kind=8) :: x, y, t_start, t_end

        ! Last time recorded
        real(kind=8) :: last_time

        ! Last time and final (x,y) written to file 
        ! (only needed for lagrangian gauges, for checkpointing)
        real(kind=8) :: t_last_written, x_last_written, y_last_written


        ! Output settings
        integer :: file_format, gtype
        real(kind=8) :: min_time_increment
        character(len=10) :: display_format
        logical, allocatable :: q_out_vars(:)
        logical, allocatable :: aux_out_vars(:)
        integer :: num_out_vars

        ! Data buffers - data holds output and time
        real(kind=8), allocatable :: data(:, :)
        integer :: level(MAX_BUFFER)

        ! Where we are in the buffer
        integer :: buffer_index
    end type gauge_type

    ! Gague array
    type(gauge_type), allocatable :: gauges(:)

    ! Gauge source info
    integer, allocatable, dimension(:) ::  mbestsrc, igauge

    logical, parameter :: INTERPOLATE = .true.

contains

    subroutine set_gauges(restart, num_eqn, num_aux, fname)

        use utility_module, only: get_value_count
        use amr_module, only: NEEDS_TO_BE_SET

        implicit none

        ! Input
        logical, intent(in) :: restart
        integer :: num_eqn, num_aux
        character(len=*), intent(in), optional :: fname

        ! Locals
        integer :: i, n, index
        integer :: num, pos, digit
        integer, parameter :: UNIT = 7
        character(len=128) :: header_1
        character(len=40) :: q_column, aux_column
        character(len=15) :: numstr

        if (.not.module_setup) then

            ! Open file
            if (present(fname)) then
                call opendatafile(UNIT,fname)
            else
                call opendatafile(UNIT,'gauges.data')
            endif

            read(UNIT,*) num_gauges
            allocate(gauges(num_gauges))
            
            ! Initialize gauge source data
            allocate(mbestsrc(num_gauges))
            mbestsrc = 0

            ! Original gauge information
            do i=1,num_gauges
                read(UNIT, *) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, &
                              gauges(i)%t_start, gauges(i)%t_end
                ! note that for lagrangian gauges, the x,y values read here 
                ! might be overwritten if this is a restart
                gauges(i)%buffer_index = 1
                ! keep track of last position for lagrangian gauges,
                ! initialize here in case checkpoint happens before 
                ! ever writing this gauge:
                gauges(i)%t_last_written = NEEDS_TO_BE_SET
                gauges(i)%x_last_written = gauges(i)%x
                gauges(i)%y_last_written = gauges(i)%y
            enddo

            ! Read in output formats
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%file_format, i=1, num_gauges)
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%display_format, i=1, num_gauges)
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%min_time_increment, i=1, num_gauges)
            read(UNIT, *)
            read(UNIT, *)
            read(UNIT, *) (gauges(i)%gtype, i=1, num_gauges)


            ! Read in q fields
            read(UNIT, *)
            read(UNIT, *)
            do i = 1, num_gauges

                ! initialize last_time so that first gauge output will be
                ! at time gauges(i)%t_start regardless of min_time_increment:
                gauges(i)%last_time = gauges(i)%t_start - 1.d0 &
                                      - gauges(i)%min_time_increment

                allocate(gauges(i)%q_out_vars(num_eqn))
                read(UNIT, *) gauges(i)%q_out_vars

                ! Count number of vars to be output
                gauges(i)%num_out_vars = 0
                do n = 1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if
                end do
            end do

            ! Read in aux fields (num_aux > 0 for geoclaw)
            read(UNIT, *)
            read(UNIT, *)
            do i = 1, num_gauges
                allocate(gauges(i)%aux_out_vars(num_aux))
                read(UNIT, *) gauges(i)%aux_out_vars

                ! Count number of vars to be output
                do n = 1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        gauges(i)%num_out_vars = gauges(i)%num_out_vars + 1
                    end if
                end do
            end do

            ! Count eta as one of the out vars
            gauges(:)%num_out_vars = gauges(:)%num_out_vars + 1

            close(UNIT)
            ! Done reading =====================================================

            ! Allocate data buffer - Note extra var out due to eta
            do i = 1, num_gauges
                allocate(gauges(i)%data(gauges(i)%num_out_vars + 2, MAX_BUFFER))
            end do

            ! Create gauge output files
            do i = 1, num_gauges
                num = gauges(i)%gauge_num

                ! convert num to string numstr with zero padding if <5 digits
                ! since we want format gauge00012.txt or gauge1234567.txt:
                write (numstr,'(I0.5)') num
                gauges(i)%file_name = 'gauge'//trim(numstr)//'.txt'

                if (gauges(i)%file_format >= 2) then
                    gauges(i)%file_name_bin = 'gauge'//trim(numstr)//'.bin'
                endif


                ! for restart, do not need to know if gauge file already exists,
                ! so some code removed
                
                if (.not. restart) then

                    if (gauges(i)%file_format >= 2) then
                        ! remove old binary file if it exists:
                        !write(6,*) 'Removing old file ',gauges(i)%file_name_bin
                        open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name_bin, &
                             status='unknown', access='stream')
                        close(unit=OUTGAUGEUNIT, status='delete')
                    endif

                    ! open and rewind old ascii file if it exists:
                    open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name,   &
                         status='unknown', position='rewind', form='formatted')
                    ! will be closed after writing header below
                endif
                     

                if (.not. restart) then
                    ! Write header to .txt file:
                    header_1 = "('# gauge_id= ',i0,' " //                 &
                               "location=( ',1e17.10,' ',1e17.10,' ) " // &
                               "num_var= ',i2)"
                    write(OUTGAUGEUNIT, header_1) gauges(i)%gauge_num,    &
                                                  gauges(i)%x,            &
                                                  gauges(i)%y,            &
                                                  gauges(i)%num_out_vars

                    if (gauges(i)%gtype == 1) then
                        ! Standard gauge
                        write(OUTGAUGEUNIT, '(a)') "# Stationary gauge"
                      else if (gauges(i)%gtype == 2) then
                        ! Lagrangian gauge
                        write(OUTGAUGEUNIT, '(a)') "# Lagrangian particle,"// &
                                    " q[2,3] replaced by (x(t),y(t))"
                      endif

                    ! Construct column labels
                    index = 0
                    q_column = "["
                    do n=1, size(gauges(i)%q_out_vars, 1)
                        if (gauges(i)%q_out_vars(n)) then
                            write(q_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if  
                    end do
                    q_column(3 * index + 2:4 + 3 * index) = "],"

                    aux_column = "["
                    index = 0
                    do n=1, size(gauges(i)%aux_out_vars, 1)
                        if (gauges(i)%aux_out_vars(n)) then
                            write(aux_column(3 * index + 2:4 + 3 * index), "(i3)") n
                            index = index + 1
                        end if  
                    end do
                    aux_column(3 * index + 2:4 + 3 * index) = "]"

                    write(OUTGAUGEUNIT, "(a,a,a,a)") "# level, time, q",       &
                                           trim(q_column), " eta, aux",        &
                                           trim(aux_column)
                    if (gauges(i)%file_format == 1) then
                        write(OUTGAUGEUNIT, '(a)') &
                          "# file format ascii, time series follow in this file"
                    else if (gauges(i)%file_format == 2) then
                        write(OUTGAUGEUNIT, '(a)') &
                            "# file format binary64, time series in .bin file"
                    else if (gauges(i)%file_format == 3) then
                        write(OUTGAUGEUNIT, '(a)') &
                            "# file format binary32, time series in .bin file"
                    endif
                    
                   close(OUTGAUGEUNIT)

               endif  ! end of ascii header file


            end do

            module_setup = .true.
        end if

    end subroutine set_gauges


!
! --------------------------------------------------------------------
!
    subroutine setbestsrc(igauge)
!
!     Called every time grids change, to set the best source grid patch
!     for each gauge, i.e. the finest level patch that includes the gauge.
!
!     lbase is grid level that didn't change, but since fine
!     grid may have disappeared, we still have to look starting
!     at coarsest level 1.
!
        use amr_module
        implicit none

        integer, intent(in), optional :: igauge
        integer :: lev, mptr, i, i1, i2

!
! ##  set source grid for each loc from finest level to coarsest
! ##  and goes on to the next gauge once a patch is found.

! ##  This modified version allows an optional igauge argument
! ##  to only update one gauge, for Lagrangian gauges (particle tracking)

        if (present(igauge)) then
            ! only loop over one gauge
            i1 = igauge
            i2 = igauge
          else
            ! normal case of setting for all gauges
            i1 = 1
            i2 = num_gauges
          endif

        ! for debugging, initialize sources to 0 then check that all set
        mbestsrc(i1:i2) = 0

        !! reorder loop for better performance with O(10^5) grids
        !! for each gauge find best source grid for its data
        do 40 i = i1,i2

           do 30 lev = lfine, 1, -1
            mptr = lstart(lev)
 20              if ((gauges(i)%x >= rnode(cornxlo,mptr)) .and. &
                        (gauges(i)%x <= rnode(cornxhi,mptr)) .and. &  
                        (gauges(i)%y >= rnode(cornylo,mptr)) .and. &
                        (gauges(i)%y <= rnode(cornyhi,mptr)) ) then
                        mbestsrc(i) = mptr
                      !! best source found for this gauge, go to next one
                      !! we know its the best because we started at finest level
                      go to 40  ! on to next gauge 
                 else 
                mptr = node(levelptr, mptr)
                    if (mptr .ne. 0) then
                       go to 20  ! try another grid
                    else
                       go to 30  ! try next coarser level grids
            end if
                 end if
 30        continue 

          if (mbestsrc(i) .eq. 0) then
              if ((gauges(i)%x < xlower) .or. (gauges(i)%x > xupper) .or. &
                  (gauges(i)%y < ylower) .or. (gauges(i)%y > yupper)) then
                  print *, "Gauge number ",gauges(i)%gauge_num, &
                      " is outside domain"
                else
                  print *, "ERROR in setting grid src for gauge number ",&
                      gauges(i)%gauge_num
                  print *, "    x,y = ",gauges(i)%x, gauges(i)%y
                endif
              endif
 40     continue 

!!!  NO MORE qsort and mbestg arrays. 
!!! Each grid now loops over mbestsrc array to see which gauges it owns.

    end subroutine setbestsrc

!
! -------------------------------------------------------------------------
!
    subroutine update_gauges(q, aux, xlow, ylow, num_eqn, mitot, mjtot, num_aux, &
                             mptr)
!
!     This routine is called each time step for each grid patch, to output
!     gauge values for all gauges for which this patch is the best one to 
!     use (i.e. at the finest refinement level).  

!     It is called after ghost cells have been filled from adjacent grids
!     at the same level, so bilinear interpolation can be used to 
!     to compute values at any gauge location that is covered by this grid.  

!     The grid patch is designated by mptr.
!     We only want to set gauges i for which mbestsrc(i) == mptr.
!     The array mbestsrc is reset after each regridding to indicate which
!     grid patch is best to use for each gauge.

!     This is a refactoring of dumpgauge.f from Clawpack 5.2 
!     Loops over only the gauges to be handled by this grid, as specified
!     by indices from mbestg1(mptr) to mbestg2(mptr)
!     NO MORE mbestg1 and 2.  Loop through all gauges. No sorting too.

        use amr_module, only: nestlevel, nghost, timemult, rnode, node, maxvar
        use amr_module, only: maxaux, hxposs, hyposs, possk
        use geoclaw_module, only: dry_tolerance
        use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
        use geoclaw_module, only: ambient_pressure
        use storm_module, only: pressure_index

        implicit none

        ! Input
        integer, intent(in) ::  num_eqn, mitot, mjtot, num_aux, mptr
        real(kind=8), intent(in) :: q(num_eqn, mitot, mjtot)
        real(kind=8), intent(in) :: aux(num_aux, mitot, mjtot)
        real(kind=8), intent(in) :: xlow, ylow

        ! Locals
        real(kind=8) :: var(maxvar + maxaux)
        real(kind=8) :: xcent, ycent, xoff, yoff, tgrid, hx, hy
        integer :: level, icell, jcell, iindex, jindex
        integer :: i, j, n, var_index, eta_index
        real(kind=8) :: h(4), mod_dry_tolerance, topo, h_interp
        real(kind=8) :: dt, xg, yg, ug, vg  ! for lagrangian gauges
        logical :: need_to_reset_patch

        ! No gauges to record, exit
        if (num_gauges == 0) then
            return
        endif

        ! Grid info
        tgrid = rnode(timemult, mptr)
        level = node(nestlevel, mptr)
        hx = hxposs(level)
        hy = hyposs(level)
        dt = possk(level)

        ! Main Gauge Loop ======================================================
        do i = 1, num_gauges
            if (mptr .ne. mbestsrc(i)) cycle 
            need_to_reset_patch = .false.  ! set to true if gauge in ghost cell
            if (tgrid < gauges(i)%t_start .or. tgrid > gauges(i)%t_end) then
               cycle
            endif
            ! Minimum increment
            ! TODO Maybe always allow last time output recording?
            if (tgrid - gauges(i)%last_time < gauges(i)%min_time_increment) then
                cycle
            end if

            ! Compute indexing and bilinear interpolant weights
            ! and check if gauge is outside patch
            iindex =  int(.5d0 + (gauges(i)%x - xlow) / hx)
            jindex =  int(.5d0 + (gauges(i)%y - ylow) / hy)
            if ((iindex < 1 .or. iindex > mitot) .or. &
                (jindex < 1 .or. jindex > mjtot)) then
                    if (gauges(i)%gtype == 2) then
                        ! Lagrangian gauge outside domain, ignore
                        ! (maybe should flag so never checked again??)
                        cycle
                    else
                        ! standard gauge so should be in domain
                        write(6,601) gauges(i)%gauge_num
601                     format("ERROR in output of Gauge Data: Gauge ", &
                             i6, " outside patch -- Skipping")
                        !write(6,*) '+++ gauges(i)%gauge_num', gauges(i)%gauge_num
                        !print *, "iindex, mitot, jindex, mjtot: ", &
                        !        iindex, mitot, jindex, mjtot
                        !print *, "gauges(i)%x = ",gauges(i)%x
                        cycle
                    endif
            else if ((iindex < nghost .or. iindex > mitot-nghost) .or. &
                (jindex < nghost .or. jindex > mjtot-nghost)) then
                    ! might want to suppress this warning...
                    if (gauges(i)%gtype == 2) then
                        write(6,602) gauges(i)%gauge_num
602                     format("WARNING in output of Gauge Data: Gauge ", &
                             i6, " in ghost cells, will reset patch")
                        ! Lagrangian gauge moved into ghost cells, 
                        ! so reset this gauge after updating
                        need_to_reset_patch = .true.
                        endif
            end if

            if (INTERPOLATE) then

                xcent  = xlow + (iindex - 0.5d0) * hx
                ycent  = ylow + (jindex - 0.5d0) * hy
                xoff   = (gauges(i)%x - xcent) / hx
                yoff   = (gauges(i)%y - ycent) / hy
    
                ! Gauge interpolation seems to work, so error test is commented out.
                ! For debugging, use the code below...
                !   Note: we expect 0 <= xoff, yoff <= 1 but if gauge is exactly 
                !   at center of cell these might be off by rounding error
    
                !if (xoff .lt. -1.d-4 .or. xoff .gt. 1.0001d0 .or. &
                !    yoff .lt. -1.d-4 .or. yoff .gt. 1.0001d0) then
                !   write(6,*) "*** print_gauges: Interpolation problem at gauge ",&
                !               igauge(i)
                !   write(6,*) "    xoff,yoff: ", xoff,yoff
                !endif
    
    
                ! Modified below from amrclaw/src/2d/gauges_module.f90 
                ! to interpolate only where all four cells are
                ! wet, otherwise just take this cell value:
    
                ! Check for dry cells by comparing h to mod_dry_tolerance, which 
                ! should be smaller than drytolerance to avoid oscillations since  
                ! when h < drytolerance the velocities are zeroed out which can then 
                ! lead to increase in h again.
    
                mod_dry_tolerance = 0.1d0 * dry_tolerance
    
                h(1) = q(1, iindex, jindex) 
                h(2) = q(1, iindex + 1, jindex) 
                h(3) = q(1, iindex, jindex + 1)
                h(4) = q(1, iindex + 1,jindex + 1)
                
                endif
        
            if ((.not. INTERPOLATE) .or.         &
                (h(1) < mod_dry_tolerance) .or.  &
                (h(2) < mod_dry_tolerance) .or.  &
                (h(3) < mod_dry_tolerance) .or.  &
                (h(4) < mod_dry_tolerance)) then

                ! If never interpolating, or if 
                ! one of the cells is dry, just use value from grid cell
                ! that contains gauge rather than interpolating
            
                icell = int(1.d0 + (gauges(i)%x - xlow) / hx)
                jcell = int(1.d0 + (gauges(i)%y - ylow) / hy)

                h_interp = q(1, icell, jcell)
                var_index = 0
                do n = 1, size(gauges(i)%q_out_vars, 1)
                    if (gauges(i)%q_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) = q(n, icell, jcell) 
                    end if
                enddo
                
                ! Note here that we skip one of the var indices to accomodate 
                ! eta here.
                var_index = var_index + 1
                eta_index = var_index

                do n=1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) = aux(n, icell , jcell)
                    end if
                end do

                ! This is the bottom layer and we should figure out the
                ! topography
                topo = aux(1, icell, jcell)

            else
                ! Linear interpolation between four cells
                ! Count for number of variables written to var    
                var_index = 0
            do n = 1, size(gauges(i)%q_out_vars, 1)
                if (gauges(i)%q_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) =                                       &
                           (1.d0 - xoff) * (1.d0 - yoff) * q(n,iindex,jindex)  &
                          + xoff*(1.d0 - yoff) * q(n,iindex+1,jindex)          &
                          + (1.d0 - xoff) * yoff * q(n,iindex,jindex+1)        &
                          + xoff * yoff * q(n,iindex+1,jindex+1)
                    end if
                end do
                
                ! Note here that we skip one of the var indices to accomodate 
                ! eta here.
                var_index = var_index + 1
                eta_index = var_index
                
                do n=1, size(gauges(i)%aux_out_vars, 1)
                    if (gauges(i)%aux_out_vars(n)) then
                        var_index = var_index + 1
                        var(var_index) =                                       &
                         (1.d0 - xoff) * (1.d0 - yoff) * aux(n,iindex,jindex)  &
                        + xoff*(1.d0 - yoff) * aux(n,iindex+1,jindex)  &
                        + (1.d0 - xoff) * yoff * aux(n,iindex,jindex+1)  &
                        + xoff * yoff * aux(n,iindex+1,jindex+1)
                    end if
                end do
                topo = (1.d0 - xoff) * (1.d0 - yoff) * aux(1,iindex,jindex)  &
                        + xoff * (1.d0 - yoff) * aux(1,iindex+1,jindex)  &
                        + (1.d0 - xoff) * yoff * aux(1,iindex,jindex+1)  &
                        + xoff * yoff * aux(1,iindex+1,jindex+1)

                ! Explicitly do depth in case the depth is not computed above
                if (.not. gauges(i)%q_out_vars(1)) then
                    h_interp = (1.d0 - xoff) * (1.d0 - yoff) &
                               * q(1,iindex,jindex)  &
                         + xoff*(1.d0 - yoff) * q(1,iindex+1,jindex)  &
                         + (1.d0 - xoff) * yoff * q(1,iindex,jindex+1)  &
                         + xoff * yoff * q(1,iindex+1,jindex+1)
                else
                    h_interp = var(1)
                end if

            end if

            ! Check to make sure we grabbed all the values
            if (gauges(i)%num_out_vars /= var_index) then
                print *, gauges(i)%num_out_vars, var_index
                print *, gauges(i)%q_out_vars
                print *, gauges(i)%aux_out_vars
                stop "Somehow we did not grab all the values we wanted..."
            end if

            ! Extract surfaces
            var(eta_index) = h_interp + topo

            ! Zero out tiny values to prevent later problems reading data,
            ! as done in valout.f
            do j = 1, gauges(i)%num_out_vars
                if (abs(var(j)) < 1d-90) var(j) = 0.d0
            end do
       
            if (gauges(i)%gtype == 2) then
                ! Lagrangian gauge, update location
                if (var(1) < dry_tolerance) then
                    ug = 0.d0
                    vg = 0.d0
                  else
                    ug = var(2)/var(1)
                    vg = var(3)/var(1)
                  endif
                if (coordinate_system == 1) then
                    ! x,y in meters, ug,vg in m/s
                    xg = gauges(i)%x + dt*ug
                    yg = gauges(i)%y + dt*vg
                  else
                    ! x,y in degrees, ug,vg in m/s
                    xg = gauges(i)%x + dt*ug / (earth_radius*deg2rad &
                                             * cos(deg2rad*gauges(i)%y))
                    yg = gauges(i)%y + dt*vg / (earth_radius*deg2rad)
                  endif
                
                ! Update location and store new xg,yg in place of hu,hv 
                gauges(i)%x = xg
                gauges(i)%y = yg
                var(2) = xg
                var(3) = yg
                endif
                
            ! save info for this time 
            n = gauges(i)%buffer_index
     
            gauges(i)%level(n) = level
            gauges(i)%data(1,n) = tgrid
            do j = 1, gauges(i)%num_out_vars
                gauges(i)%data(1 + j, n) = var(j)
            end do
            
            gauges(i)%buffer_index = n + 1
            if (gauges(i)%buffer_index > MAX_BUFFER) then
                call print_gauges_and_reset_nextLoc(i)
            endif

            gauges(i)%last_time = tgrid

            if (need_to_reset_patch) then
                call setbestsrc(i)
                endif

        end do ! End of gauge loop =============================================

    end subroutine update_gauges
!
! -------------------------------------------------------------------------
! Write out gauge data for the gauge specified
!
    subroutine print_gauges_and_reset_nextLoc(gauge_num)

        implicit none

        ! Input
        integer, intent(in) :: gauge_num

        ! Locals
        integer :: j, k, myunit, nvals, ntimes
        integer :: omp_get_thread_num, mythread
        character(len=32) :: out_format
        real(kind=4), allocatable :: gdata4(:,:)
        real(kind=8), allocatable :: gdata8(:,:)

        ! Loop through gauge's buffer writing out all available data.  Also
        ! reset buffer_index back to beginning of buffer since we are emptying
        ! the buffer here

        ! Rearranged for v5.9.0 to also allow writing binary format

        j = gauges(gauge_num)%buffer_index - 1
        if (j > 0) then
            gauges(gauge_num)%t_last_written = gauges(gauge_num)%data(1, j)
            gauges(gauge_num)%x_last_written = gauges(gauge_num)%data(3, j)
            gauges(gauge_num)%y_last_written = gauges(gauge_num)%data(4, j)
        endif
        
        nvals = gauges(gauge_num)%num_out_vars + 1
        ntimes = gauges(gauge_num)%buffer_index - 1

        ! Open unit dependent on thread number
        mythread = 0
!$      mythread = omp_get_thread_num()
        myunit = OUTGAUGEUNIT + mythread

        if (gauges(gauge_num)%file_format == 1) then

            ! ascii output

            ! Construct output format based on number of output variables and
            ! request format
            write(out_format, "(A7, i2, A6, A1)") "(i5.2,",                    &
                                        gauges(gauge_num)%num_out_vars + 1,    &
                                        gauges(gauge_num)%display_format, ")"

            ! Open gauge file:
            open(unit=myunit, file=gauges(gauge_num)%file_name,       &
                 status='unknown', position='append', form='formatted')

            do j = 1, ntimes
                write(myunit, out_format) gauges(gauge_num)%level(j), &
                      (gauges(gauge_num)%data(k, j), k=1,nvals)
            end do

        else if (gauges(gauge_num)%file_format >= 2) then

            ! binary output

            open(unit=myunit, file=gauges(gauge_num)%file_name_bin, &
                 status='unknown', position='append',access='stream')
            
            if (gauges(gauge_num)%file_format == 3) then
                allocate(gdata4(nvals+1, ntimes))
                gdata4(1, 1:ntimes) = real(gauges(gauge_num)%level(1:ntimes), kind=4)
                gdata4(2:nvals+1, 1:ntimes) = &
                        real(gauges(gauge_num)%data(1:nvals,1:ntimes), kind=4)
                write(myunit) gdata4
                deallocate(gdata4)  
            else
                allocate(gdata8(nvals+1, ntimes))
                gdata8(1, 1:ntimes) = real(gauges(gauge_num)%level(1:ntimes), kind=8)
                gdata8(2:nvals+1, 1:ntimes) = &
                        real(gauges(gauge_num)%data(1:nvals,1:ntimes), kind=8)
                write(myunit) gdata8
                deallocate(gdata8)
            endif


        else
            print *, "Unhandled file format ", gauges(gauge_num)%file_format
            stop
        end if

        ! close file
        close(myunit)

        gauges(gauge_num)%buffer_index = 1                        

      end subroutine print_gauges_and_reset_nextLoc

end module gauges_module
