! ============================================================================
!  File:        gauges_module (for 1D)
! ============================================================================
!    Copyright (C) 2010-04-21 Clawpack Developers http://www.clawpack.org
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module gauges_module

      implicit none


      !module for outputing time-series solution at every time-step

      !variable to specify gauges
      integer :: mgauges
      integer, parameter :: OUTGAUGEUNIT = 277
      double precision, allocatable :: xgauge(:),t0gauge(:),tFgauge(:)
      double precision, allocatable :: sol(:)
      integer, allocatable :: igauge(:)

contains
      !================================================================
      subroutine set_gauges(fname)
      !
      ! initialize the gauge locations from spec file, fname
      !================================================================
         implicit none

         ! i/o
         character*25, optional, intent(in) :: fname

         !locals
         character*25 :: file_name
         real(kind=8) :: ygarbage
         integer :: iunit = 127
         integer :: mvars,l
         logical :: found_file

         ! Read user parameters from setgauges.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'setgauges.data'
         endif

         inquire(file=file_name,exist=found_file)
         if (found_file) then
            call opendatafile(iunit, file_name)
            read(iunit,*) mgauges
            allocate(igauge(mgauges),xgauge(mgauges))
            allocate(t0gauge(mgauges),tFgauge(mgauges))
            do l = 1, mgauges
               read(iunit,*) igauge(l),xgauge(l),ygarbage,t0gauge(l),tFgauge(l)
            enddo
            close(iunit)
            open(unit=OUTGAUGEUNIT,file='fort.gauge',status='unknown',form='formatted')
         endif


      end subroutine set_gauges

      !=================================================================
      subroutine return_gauge(meqn,maux,mvars,mx,dx,xlower,xgauge,sol,q,aux)
      !
      ! return values sol(vars) at xlower
      !=================================================================

      !i/o
      integer, intent(in) :: meqn,maux,mvars,mx
      double precision, intent(in) :: dx,xlower,xgauge
      double precision, intent(in) :: q(meqn,mx), aux(maux,mx)
      double precision, intent(out) :: sol(mvars)

      !locals
      integer :: i,m

      i = floor((xgauge - xlower)/dx) + 1

      do m = 1, meqn
         sol(m) = q(m,i)
      enddo
      do m = 1, maux
         sol(m+meqn) = aux(m,i)
      enddo

      return
      end subroutine

end module gauges_module



