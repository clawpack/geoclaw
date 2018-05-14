
subroutine set_sea_level(mbc,mx,my,xlow,ylow,dx,dy,t,vsea)

    ! Set variable sea level at time t in the vsea array,
    ! for all grid points on the patch defined by mx,my,xlow,ylow,dx,dy.

    ! These values are used when creating a new refined grid patch near
    ! the shore and it is necessary to determine the depth of water to
    ! initialize grid cells for which the underlying coarse cell was dry.

    ! Set vsea(i,j) for i=1,mx and j=1,my,  to the value of sea level
    ! to be used at x = xlow + (i-1)*dx,  y = ylow + (j-1)*dy.

    ! This version uses the dtopo files to determine the modification needed
    ! to the original scalar sea_level value from the geoclaw_module,
    ! as required at time t, due to subsidence or uplift.

    ! The user can copy and modify this file in order to specify some
    ! other spatial/temporal variation in sea level.

    use topo_module
    use geoclaw_module, only: sea_level

    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,t
    real(kind=8), intent(inout) :: vsea(1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Local
    integer :: i,j,i1,i2,j1,j2,idtopo,jdtopo,kdtopo, &
               index0_dtopowork,ij,m
    real(kind=8) :: x,y

    vsea = sea_level  ! initialize array to sea_level, update elements below.

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

        index0_dtopowork = i0dtopo(m) + (kdtopo-1)*mxdtopo(m)*mydtopo(m)
        !write(6,*) '+++ index0_dtopowork = ',index0_dtopowork

        ! Adjust sea_level by dtopo on part of patch that overlaps dtopo.
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
                ij = index0_dtopowork + (jdtopo-1)*mxdtopo(m) + idtopo-1
                vsea(i,j) = vsea(i,j) + dtopowork(ij)
                enddo
            enddo

        enddo
    return

end subroutine set_sea_level
