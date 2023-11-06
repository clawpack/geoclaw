! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Modified for geoclaw_1d/src/1d_nonuniform code based on $CLAW/classic

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine print_gauges
!     Called each time step 
!
!     Note: by default all components of q are printed at each gauge.
!     To print something different or a different precision, modify 
!     format statement 100 and/or the write statement that uses it.
!   

module gauges_module

    implicit none
    save

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges, inum
    real(kind=8), allocatable :: xgauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable, dimension(:) ::  igauge, nextLoc
    
    integer, parameter :: MAXDATA=1000
    real(kind=8), pointer :: gaugeArray(:,:,:)

contains

    subroutine set_gauges(restart, nvar, fname)

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname
        logical, intent(in) :: restart
        integer, intent(in) :: nvar

        ! Locals
        integer :: i, ipos, idigit
        integer, parameter :: iunit = 7
        character*14 ::  fileName

        if (.not. module_setup) then

            ! Open file
            if (present(fname)) then
                call opendatafile(iunit,fname)
            else
                call opendatafile(iunit,'gauges.data')
            endif

            read(iunit,*) num_gauges

            allocate(xgauge(num_gauges))
            allocate(t1gauge(num_gauges), t2gauge(num_gauges))
            allocate(igauge(num_gauges))

            allocate(nextLoc(num_gauges))
            allocate(gaugeArray(nvar+2,MAXDATA,num_gauges)) !+2 for time,eta
            
            do i=1,num_gauges
                read(iunit,*) igauge(i),xgauge(i),t1gauge(i),t2gauge(i)
            enddo

            close(iunit)
            
            ! initialize for starters
            nextLoc  = 1  ! next location to be filled with gauge info

            do i = 1, num_gauges
               fileName = 'gaugexxxxx.txt'    ! NB different name convention too
               inum = igauge(i)
               do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
                  idigit = mod(inum,10)
                  fileName(ipos:ipos) = char(ichar('0') + idigit)
                  inum = inum / 10
               end do

    !          status unknown since might be a restart run. might need to rewind
               if (restart) then
                  open(unit=OUTGAUGEUNIT, file=fileName, status='old',        &
                       position='append', form='formatted')
               else
                  open(unit=OUTGAUGEUNIT, file=fileName, status='unknown',        &
                       position='append', form='formatted')
                  rewind OUTGAUGEUNIT
                  ! write out nvar+1 variables h, hu, eta=h+B
                  write(OUTGAUGEUNIT,100) igauge(i), xgauge(i), nvar+1
  100             format("# gauge_id= ",i5," location=( ",1e15.7," ) num_eqn= ",i2)
                  write(OUTGAUGEUNIT,101)
  101             format("# Columns: level time h hu eta")
               endif

               close(OUTGAUGEUNIT)

            end do

            module_setup = .true.
        end if

    end subroutine set_gauges


!
! -------------------------------------------------------------------------
!
      subroutine update_gauges(q,aux,xlow,nvar,mx,mbc,naux,dx,tgrid)
!
!     This routine is called each time step.

      implicit none

      real(kind=8), intent(in) ::  xlow,dx,tgrid
      integer, intent(in) ::  nvar,mx,mbc,naux
      real(kind=8), intent(in) ::  q(nvar,1-mbc:mx+mbc)
      real(kind=8), intent(in) ::  aux(naux,1-mbc:mx+mbc)

      ! local variables:
      real(kind=8) :: var(nvar)
      real(kind=8) :: xcent,xoff,topo
      integer :: i,iindex,ivar, ii,i1,i2
      integer :: nindex

!     write(*,*) '+++ in print_gauges with num_gauges = ',num_gauges

      if (num_gauges == 0) then
         return
      endif

!     write(*,*) 'tgrid = ',tgrid

      do 10 ii=1,num_gauges
!       write(6,*) '+++ interploting for gauge ', ii
        iindex =  int(.5d0 + (xgauge(ii)-xlow)/dx)
        if (iindex .lt. 1-mbc .or. iindex .gt. mx+mbc-1) &
          write(*,*)"ERROR in output of Gauge Data "
        xcent  = xlow + (iindex-.5d0)*dx
        xoff   = (xgauge(ii)-xcent)/dx
!   IF WANT TO INCLUDE THIS TEST< MODIFY FOR ROUNDOFF LEVEL DIFF
!       if (xoff .lt. 0.d0 .or. xoff .gt. 1.d0) then
!          write(6,*)" BIG PROBLEM in DUMPGAUGE", i
!       endif

!       ## linear interpolation
        do ivar = 1, nvar
           var(ivar) = (1.d0-xoff)*q(ivar,iindex) &
                   + xoff*q(ivar,iindex+1)
!          # for printing without underflow of exponent:
           if (abs(var(ivar)) .lt. 1.d-90) var(ivar) = 0.d0
        end do


       ! save info for this time
        nindex = nextLoc(ii)
 
        gaugeArray(1,nindex,ii) = tgrid
        do ivar = 1, nvar
           gaugeArray(1+ivar,nindex,ii) = var(ivar)
        end do

!       # also store surface eta computed using interpolated topography:
        topo = (1.d0-xoff)*aux(1,iindex) + xoff*aux(1,iindex+1)
        gaugeArray(nvar+2,nindex,ii) = var(1) + topo
        
        nextLoc(ii) = nextLoc(ii) + 1
        if (nextLoc(ii) .gt. MAXDATA) then
          call print_gauges_and_reset_nextLoc(ii,nvar)  
        endif
 10     continue  ! end of loop over all gauges
 
      end subroutine update_gauges
!
! -------------------------------------------------------------------------
!
      subroutine print_gauges_and_reset_nextLoc(gaugeNum,nvar)
!
!    Array of gauge data for this gauge reached max capacity
!    print to file.

      implicit none
      integer :: gaugeNum,nvar,j,inum,k,idigit,ipos,myunit
      integer :: omp_get_thread_num, mythread
      character*14 :: fileName

      ! open file for gauge gaugeNum, figure out name
      ! not writing gauge number since it is in file name now
      ! status is old, since file name and info written for
      ! each file in in set_gauges.
      !
      ! NB: output written in different order, losing backward compatibility


      fileName = 'gaugexxxxx.txt'    ! NB different name convention too
      inum = igauge(gaugeNum)
      do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
         idigit = mod(inum,10)
         fileName(ipos:ipos) = char(ichar('0') + idigit)
         inum = inum / 10
      end do


      open(unit=OUTGAUGEUNIT, file=fileName, status='old',    &
           position='append', form='formatted')
      
      ! called either because array is full (write MAXDATA amount of gauge data)
      ! or checkpoint time, so write whatever is in array and reset.
      ! nextLoc has already been increment before this subr. called
      do j = 1, nextLoc(gaugeNum)-1
        write(OUTGAUGEUNIT,100) 1, (gaugeArray(k,j,gaugeNum),k=1,nvar+2) 
            ! writes 1 for level and gaugeArray includes time in first col
      end do
      nextLoc(gaugeNum) = 1                        

      ! if you want to modify number of digits printed, modify this...
100     format(i5.2, 10e15.7)

      ! close file
      close(OUTGAUGEUNIT)

      end subroutine print_gauges_and_reset_nextLoc

end module gauges_module
