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

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT = 89
    integer :: num_gauges

    integer, parameter :: MAX_BUFFER = 10

    ! Gauge data types
    type gauge_type
        ! Gauge number
        integer :: gauge_num

        character(len=14) :: file_name

        ! Location in time and space
        real(kind=8) :: x, y, t_start, t_end

        ! Output settings
        character(len=10) :: format
        integer, allocatable :: q_out_vars(:)
        integer, allocatable :: aux_out_vars(:)
        integer :: num_out_vars

        ! Data buffers - data holds output and time
        real(kind=8), allocatable :: data(:, :)
        integer :: level(MAX_BUFFER)

        ! Where we are in the buffer
        integer :: buffer_index
    end type gauge_type

    ! Gague array
    type(gauge_type), allocatable :: gauges(:)

    ! real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable, dimension(:) ::  mbestsrc, mbestorder, &
                          igauge, mbestg1, mbestg2

contains

    subroutine set_gauges(restart, nvar, fname)

        use amr_module, only: maxgr

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname
        logical, intent(in) :: restart
        integer, intent(in) :: nvar

        ! Locals
        integer :: i, ipos, idigit, inum
        integer, parameter :: iunit = 7

        if (.not.module_setup) then
            ! Open file
            if (present(fname)) then
                call opendatafile(iunit,fname)
            else
                call opendatafile(iunit,'gauges.data')
            endif

            read(iunit,*) num_gauges

            allocate(gauges(num_gauges))

            allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
            mbestsrc = 0
            allocate(mbestg1(maxgr), mbestg2(maxgr))
            
            do i=1,num_gauges
                read(iunit,*) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, &
                              gauges(i)%t_start,gauges(i)%t_end
                gauges(i)%buffer_index = 1

                ! Temporary setting here, need input method
                allocate(gauges(i)%q_out_vars(4))
                gauges(i)%q_out_vars = [1, 2, 3]
                allocate(gauges(i)%aux_out_vars(0))
                ! gauges(i)%aux_out_vars = []
                gauges(i)%num_out_vars = size(gauges(i)%q_out_vars, 1) +   &
                                         size(gauges(i)%aux_out_vars, 1)
                allocate(gauges(i)%data(gauges(i)%num_out_vars + 2, MAX_BUFFER))

                gauges(i)%format = "e15.6"
            enddo

            close(iunit)
            
            ! initialize for starters

            do i = 1, num_gauges
               gauges(i)%file_name = 'gaugexxxxx.txt'
               inum = gauges(i)%gauge_num
               do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
                  idigit = mod(inum,10)
                  gauges(i)%file_name(ipos:ipos) = char(ichar('0') + idigit)
                  inum = inum / 10
               end do

    !          status unknown since might be a restart run. maybe need to test and rewind?
               if (restart) then
                  open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name, status='old',        &
                       position='append', form='formatted')
               else
                  open(unit=OUTGAUGEUNIT, file=gauges(i)%file_name, status='unknown',        &
                       position='append', form='formatted')
                  rewind OUTGAUGEUNIT
                  write(OUTGAUGEUNIT, 100) gauges(i)%gauge_num, gauges(i)%x, gauges(i)%y, gauges(i)%num_out_vars + 1
 100              format("# gauge_id= ",i5," location=( ",1e15.7," ",1e15.7," ) num_eqn= ",i2)
                  write(OUTGAUGEUNIT, 101)
 101              format("# Columns: level time q_vars eta aux_vars")
               endif

               close(OUTGAUGEUNIT)

            end do

            module_setup = .true.
        end if

    end subroutine set_gauges


!
! --------------------------------------------------------------------
!
      subroutine setbestsrc()
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

      integer :: lev, mptr, i, k1, ki

!
! ##  set source grid for each loc from coarsest level to finest.
! ##  that way finest src grid left and old ones overwritten
! ##  this code uses fact that grids do not overlap

! # for debugging, initialize sources to 0 then check that all set
        do i = 1, num_gauges
            mbestsrc(i) = 0
        end do

        do lev = 1, lfine  
            mptr = lstart(lev)
            do
                do i = 1, num_gauges
                    if ((gauges(i)%x >= rnode(cornxlo,mptr)) .and. &
                        (gauges(i)%x <= rnode(cornxhi,mptr)) .and. &  
                        (gauges(i)%y >= rnode(cornylo,mptr)) .and. &
                        (gauges(i)%y <= rnode(cornyhi,mptr)) ) then
                        mbestsrc(i) = mptr
                    end if
                end do
                mptr = node(levelptr, mptr)
                if (mptr == 0) exit
            end do 
        end do

        do i = 1, num_gauges
            if (mbestsrc(i) .eq. 0) &
               print *, "ERROR in setting grid src for gauge data", i
        end do

        ! Sort the source arrays for easy testing during integration
        call qsorti(mbestorder, num_gauges, mbestsrc)

!     After sorting,  
!           mbestsrc(mbestorder(i)) = grid index to be used for gauge i
!     and mbestsrc(mbestorder(i)) is non-decreasing as i=1,2,..., num_gauges

!     write(6,*) '+++ mbestorder: ',mbestorder
!     write(6,*) '+++ mbestsrc: ',mbestsrc

!     Figure out the set of gauges that should be handled on each grid:  
!     after loop below, grid k should handle gauges numbered
!          mbestorder(i) for i = mbestg1(k), mbestg1(k)+1, ..., mbestg2(k)
!     This will be used for looping in print_gauges subroutine.

      ! initialize arrays to default indicating grids that contain no gauges:
        mbestg1 = 0
        mbestg2 = 0

        k1 = 0
        do i=1,num_gauges
            ki = mbestsrc(mbestorder(i))
            if (ki > k1) then
                ! new grid number seen for first time in list
                if (k1 > 0) then
                    ! mark end of gauges seen by previous grid
                    mbestg2(k1) = i-1
!                   write(6,*) '+++ k1, mbestg2(k1): ',k1,mbestg2(k1)
                endif
                mbestg1(ki) = i
!               write(6,*) '+++ ki, mbestg1(ki): ',ki,mbestg1(ki)
            endif
            k1 = ki
        enddo
        if (num_gauges > 0) then
            ! finalize 
            mbestg2(ki) = num_gauges
        endif

    end subroutine setbestsrc

!
! -------------------------------------------------------------------------
!
    subroutine update_gauges(q,aux,xlow,ylow,nvar,mitot,mjtot,naux,mptr)
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

        use amr_module
        use geoclaw_module, only: dry_tolerance
  
        implicit none
  
        real(kind=8), intent(in) ::  q(nvar,mitot,mjtot)
        real(kind=8), intent(in) ::  aux(naux,mitot,mjtot)
        real(kind=8), intent(in) ::  xlow,ylow
        integer, intent(in) ::  nvar,mitot,mjtot,naux,mptr
  
        ! local variables:
        real(kind=8) :: var(maxvar)
        real(kind=8) :: xcent,ycent,xoff,yoff,tgrid,hx,hy
        integer :: level,i,j,ioff,joff,iindex,jindex,ivar, ii,i1,i2
        real(kind=8) :: h(4),drytol2,topo,eta
        integer :: icell,jcell, index
        integer :: q_var_start, q_var_end, aux_var_start, aux_var_end
  
!       write(*,*) '+++ in print_gauges with num_gauges, mptr = ',num_gauges,mptr
  
        if (num_gauges == 0) then
           return
        endif
  
        i1 = mbestg1(mptr)
        i2 = mbestg2(mptr)
  
        if (i1 == 0) then
           ! no gauges to be handled by this grid
           return
        endif
  
!       write(6,*) '+++ mbestg1(mptr) = ',mbestg1(mptr)
!       write(6,*) '+++ mbestg2(mptr) = ',mbestg2(mptr)
  
!       # this stuff the same for all gauges on this grid
        tgrid = rnode(timemult,mptr)
        level = node(nestlevel,mptr)
        hx    =  hxposs(level)
        hy    =  hyposs(level)

!     write(*,*) 'tgrid = ',tgrid

        do i = i1,i2
            ii = mbestorder(i)
            if (mptr /= mbestsrc(ii)) then
                print *, '*** should not happen... i, ii, mbestsrc(ii), mptr:'
                print *, i, ii, mbestsrc(ii), mptr
                stop
            endif
            if (tgrid < gauges(ii)%t_start .or. tgrid > gauges(ii)%t_end) then
                ! This gauge should not be output at this time
                cycle
            endif

            ! Bilinear interpolation at gauge location
            ! Note: changed 0.5 to  0.5d0 etc.
            iindex =  int(.5d0 + (gauges(ii)%x - xlow) / hx)
            jindex =  int(.5d0 + (gauges(ii)%y - ylow) / hy)
            if ((iindex < nghost .or. iindex > mitot-nghost) .or. &
                (jindex < nghost .or. jindex > mjtot-nghost)) then
                    print *, "ERROR in output of Gauge Data "
            end if
            xcent  = xlow + (iindex - 0.5d0) * hx
            ycent  = ylow + (jindex - 0.5d0) * hy
            xoff   = (gauges(ii)%x - xcent) / hx
            yoff   = (gauges(ii)%y - ycent) / hy
!  IF WANT TO USE, MODIFY TO TEST FOR ROUNDOFF LEVEL DIFF
!       if (xoff .lt. 0.d0 .or. xoff .gt. 1.d0 .or. &
!           yoff .lt. 0.d0 .or. yoff .gt. 1.d0) then
!          write(6,*)" BIG PROBLEM in DUMPGAUGE", i
!       endif

            ! ## Modified below from amrclaw/src/2d/gauges_module.f90 
            ! ## to interpolate only where all four cells are
            ! ## wet, otherwise just take this cell value:

            ! Check for dry cells by comparing h to drytol2, which should be
            ! smaller than drytolerance to avoid oscillations since when 
            ! h < drytolerance the velocities are zeroed out which can then lead
            ! to increase in h again.

            drytol2 = 0.1d0 * dry_tolerance

            h(1) = q(1, iindex,     jindex) 
            h(2) = q(1, iindex + 1, jindex) 
            h(3) = q(1, iindex,     jindex + 1)
            h(4) = q(1, iindex + 1, jindex + 1) 
            
            q_var_start = 1
            q_var_end = size(gauges(ii)%q_out_vars, 1)
            aux_var_start = q_var_end + 1
            aux_var_end = gauges(ii)%num_out_vars

            if ((h(1) < drytol2) .or.  &
                (h(2) < drytol2) .or.  &
                (h(3) < drytol2) .or.  &
                (h(4) < drytol2)) then
                
                ! One of the cells is dry, so just use value from grid cell
                ! that contains gauge rather than interpolating
                
                icell = int(1.d0 + (gauges(ii)%x - xlow) / hx)
                jcell = int(1.d0 + (gauges(ii)%y - ylow) / hy)
                ! Loop through the q output first
                do ivar=q_var_start, q_var_end
                    var(ivar) = q(gauges(ii)%q_out_vars(ivar),icell,jcell) 
                enddo
                ! Now do the aux output
                do ivar=aux_var_start, aux_var_end
                    var(ivar) = aux(gauges(ii)%aux_out_vars(ivar - q_var_end), icell , jcell)
                end do

                ! This is the bottom layer and we should figure out the
                ! topography
                topo = aux(1, icell, jcell)
            else
                ! Linear interpolation between four cells
                do ivar=q_var_start, q_var_end
                    var(ivar) = (1.d0 - xoff) * (1.d0 - yoff) &
                               * q(gauges(ii)%q_out_vars(ivar),iindex,jindex)  &
                    + xoff*(1.d0 - yoff) * q(gauges(ii)%q_out_vars(ivar),iindex+1,jindex)  &
                    + (1.d0 - xoff) * yoff * q(gauges(ii)%q_out_vars(ivar),iindex,jindex+1)  &
                    + xoff * yoff * q(gauges(ii)%q_out_vars(ivar),iindex+1,jindex+1)
                end do
                do ivar=aux_var_start, aux_var_end
                    var(ivar) = (1.d0 - xoff) * (1.d0 - yoff) &
                               * aux(gauges(ii)%aux_out_vars(ivar - q_var_end),iindex,jindex)  &
                    + xoff*(1.d0 - yoff) * aux(gauges(ii)%aux_out_vars(ivar - q_var_end),iindex+1,jindex)  &
                    + (1.d0 - xoff) * yoff * aux(gauges(ii)%aux_out_vars(ivar - q_var_end),iindex,jindex+1)  &
                    + xoff * yoff * aux(gauges(ii)%aux_out_vars(ivar - q_var_end),iindex+1,jindex+1)
                end do
                topo = (1.d0 - xoff) * (1.d0 - yoff)  &
                        * aux(1,iindex,jindex)  &
                     + xoff * (1.d0 - yoff) * aux(1,iindex+1,jindex)  &
                     + (1.d0 - xoff) * yoff * aux(1,iindex,jindex+1)  &
                     + xoff * yoff * aux(1,iindex+1,jindex+1)
            endif

            ! Extract surfaces
            eta = var(1) + topo

            ! Zero out tiny values to prevent later problems reading data,
            ! as done in valout.f
            do j = 1, gauges(ii)%num_out_vars
                if (abs(var(j)) < 1d-90) var(j) = 0.d0
            end do
            if (abs(eta) < 1d-90) eta = 0.d0
       
            ! save info for this time 
            index = gauges(ii)%buffer_index
     
            gauges(ii)%level(index) = level
            gauges(ii)%data(1,index) = tgrid
            do j = q_var_start, q_var_end
                gauges(ii)%data(1 + j, index) = var(j)
            end do
            gauges(ii)%data(gauges(ii)%num_out_vars + 2, index) = eta
            do j = aux_var_start + 1, aux_var_end + 1
                gauges(ii)%data(1 + j, index) = var(j - 1)
            end do
            
            gauges(ii)%buffer_index = index + 1
            if (gauges(ii)%buffer_index > MAX_BUFFER) then
                call print_gauges_and_reset_nextLoc(ii)  
            endif

        end do! end of loop over all gauges
 
    end subroutine update_gauges
!
! -------------------------------------------------------------------------
!
    subroutine print_gauges_and_reset_nextLoc(gauge_num)
        ! Write out gauge data for the gauge specified

        implicit none

        integer, intent(in) :: gauge_num

        ! Locals
        integer :: j, k, myunit
        integer :: omp_get_thread_num, mythread
        character(len=32) :: out_format

        ! Construct output format based on number of output variables and request format
        write(out_format, "(A7, i2, A6, A1)") "(i5.2,",         &
               gauges(gauge_num)%num_out_vars + 2, gauges(gauge_num)%format, ")"

        ! Open unit dependent on thread number
        mythread = 0
!$      mythread = omp_get_thread_num()
        myunit = OUTGAUGEUNIT + mythread
        open(unit=myunit, file=gauges(gauge_num)%file_name, status='old', &
                          position='append', form='formatted')
      
        ! Loop through gauge's buffer writing out all available data.  Also
        ! reset buffer_index back to beginning of buffer since we are emptying
        ! the buffer here
        do j = 1, gauges(gauge_num)%buffer_index - 1
            write(myunit, out_format) gauges(gauge_num)%level(j),    &
                (gauges(gauge_num)%data(k, j), k=1, gauges(gauge_num)%num_out_vars + 2)
        end do
        gauges(gauge_num)%buffer_index = 1                        

        ! close file
        close(myunit)

    end subroutine print_gauges_and_reset_nextLoc

end module gauges_module
