subroutine b4step1(mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    ! Called before each call to step1.
    ! Use to set time-dependent aux arrays or perform other tasks.

    ! this version checks for negative depths and outputs gauge information

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)

    !local variables
    integer, i
    real(kind=8) mvars

      do i=1-mbc,maxmx+mbc
         if (q(1,i)<=drytolerance) then
            q(1,i) = max(q(1,i),drytolerance)
            do m=2,meqn
               q(m,i)=0.d0
            enddo
         endif
      enddo

      if (allocated(igauge)) then
         mvars = meqn+maux
         if (.not.allocated(sol)) then
            allocate(sol(1:mvars))
         endif

         do ig = 1,mgauges
            if (t.ge.t0gauge(ig).and.t.le.tFgauge(ig)) then
               call return_gauge(meqn,maux,mvars,mx,dx,xlower, &
     &               xgauge(ig),sol(1:mvars),q(1:meqn,1:mx),aux(1:maux,1:mx))

               write(OUTGAUGEUNIT,*) igauge(ig),1,t,(sol(j),j=1,mvars)
            endif
         enddo
      endif

end subroutine b4step1

