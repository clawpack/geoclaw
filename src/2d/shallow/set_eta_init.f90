
subroutine set_eta_init(mbc,mx,my,xlow,ylow,dx,dy,t,veta)

    ! This routine is called only if variable_eta_init = .true.
    ! The input is a single grid patch at time t, as specified by
    !     mbc,mx,my,xlow,ylow,dx,dy,t
    ! The output is an array 
    !     veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    ! with the desired initial surface elevation at all points on this patch.

    ! This is called by qinit and also by filpatch and filval when refining.

    ! This default version sets veta to the value sea_level everywhere and
    ! then adjusts this based on any dtopo deformation that has occured
    ! up to this time t, as determined by interpolation from the 
    ! (possibly time-dependent) ! dtopo files.

    ! There is also a commented-out section below indicating how you might set 
    ! a higher value in some region that contains an onshore lake, for example.


    use topo_module
    use geoclaw_module, only: sea_level

    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,t
    real(kind=8), intent(inout) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Local
    integer :: i,j,i1,i2,j1,j2,idtopo,jdtopo,kdtopo,m
    integer(kind=8) :: index0_dtopowork, ij
    real(kind=8) :: x,y
    real(kind=8) :: x1_lake,x2_lake,y1_lake,y2_lake, lake_level

    veta = sea_level  ! initialize to sea_level, update below at some (i,j)

    if (.false.) then
        ! This illustrates how you might set a higher value than sea_level
        ! in a rectangle containing an onshore lake, as an example.
        ! Depending on the geometry, a simple rectangle might not suffice.
        ! What's below has been used for Lake Ozette on the Washington coast.
        lake_level = 11.d0  ! meters relative to 0 datum of topofiles
        x1_lake = -124.67d0
        x2_lake = -124.58d0
        y1_lake = 48.03d0
        y2_lake = 48.16d0
        i1 = max(int(floor((x1_lake - xlow)/dx)), 1)
        i2 = min(int(floor((x2_lake - xlow)/dx)), mx)
        j1 = max(int(floor((y1_lake - ylow)/dy)), 1)
        j2 = min(int(floor((y2_lake - ylow)/dy)), my)

        forall(i=i1:i2, j=j1:j2)
            veta(i,j) = lake_level
        end forall

    endif

    ! The code below adjusts veta values set above for any uplift or
    ! subsidence found in this region at the current time t, which assumes
    ! that the water surface moves with the grount initially.

    if (num_dtopo == 0) return

    do m=1,num_dtopo
        if (t < t0dtopo(m)) cycle  ! this dtopo isn't moving yet
        
        ! compute indices in patch that overlap with this dtopo,
        ! cycling to the next dtopo array if we discover there's no overlap:
        
        i1 = max(int(floor((xlowdtopo(m) - xlow)/dx)), 1)
        if (i1 >= mx) cycle
        
        i2 = min(int(floor((xhidtopo(m) - xlow)/dx)), mx)
        if (i2 < 1) cycle

        j1 = max(int(floor((ylowdtopo(m) - ylow)/dy)), 1)
        if (j1 >= my) cycle
        
        j2 = min(int(floor((yhidtopo(m) - ylow)/dy)), my)
        if (j2 < 1) cycle
    
        ! There is some overlap of dtopo with this patch
        ! Next figure out index into time-dependent dtopo based on t:

        if (mtdtopo(m) == 1) then
            ! Special case: instantaneous displacement at one instant in time
            kdtopo = 1
          else
            kdtopo = int(floor((t-t0dtopo(m))/dtdtopo(m)))+1
            kdtopo = min(kdtopo,mtdtopo(m))
            kdtopo = max(kdtopo,1)
          endif

        index0_dtopowork = i0dtopo(m) + int(kdtopo-1, 8)*int(mxdtopo(m), 8)*int(mydtopo(m), 8)

        ! Adjust eta_init by dtopo on part of patch that overlaps dtopo.
        ! Code below assumes dtopo is smooth enough that we can just 
        ! evaluate at one point in space and time, not doing interpolation.
        ! This could be improved as in topo_update.f, but probably not
        ! necessary (??).

        do i=i1,i2
            x = xlow + (i-0.5d0)*dx
            do j=j1,j2
                y = ylow + (j-0.5d0)*dy
                idtopo = int(floor((x-xlowdtopo(m))/dxdtopo(m))) + 1
                idtopo = max(1, min(mxdtopo(m)-1, idtopo))
                jdtopo = int(floor((yhidtopo(m)-y)/dydtopo(m))) + 1
                jdtopo = max(1, min(mydtopo(m)-1, jdtopo))
                ij = index0_dtopowork + int(jdtopo-1, 8)*int(mxdtopo(m),8) &
                                      + int(idtopo-1, 8)

                veta(i,j) = veta(i,j) + dtopowork(ij)

                enddo
            enddo

        enddo
    return

end subroutine set_eta_init
