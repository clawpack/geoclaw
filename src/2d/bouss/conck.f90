subroutine conck(level, nvar, naux, time, rest)
    !> Conservation check for specified level.
    !! This is mostly a debugging tool and assumes grids don't overlap
    !! Modified for GeoClaw:  sum up zeta = h or h+B
 
    use amr_module, only: node,store1,storeaux,ndilo,ndjlo,ndihi,ndjhi
    use amr_module, only: lstart,nghost,levelptr,outunit,mcapa
    use amr_module, only: alloc,t0,hxposs,hyposs,possk,tmass0
    implicit none
    real(kind=8), intent(in) :: time
    integer, intent(in) :: level,nvar,naux
    logical, intent(in) :: rest

    real(kind=8) :: hx,hy,dt,totmass,zeta
    integer :: mptr,loc,locaux,nx,ny,mitot,mjtot,i,j

    ! grid loop for given level
 
    hx      = hxposs(level)
    hy      = hyposs(level)
    dt      = possk(level)
    totmass = 0.d0

    mptr = lstart(level)
    
    do while (mptr > 0)
        
        loc    = node(store1,mptr)
        locaux = node(storeaux,mptr)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost

        do j  = nghost+1, mjtot-nghost
            do i  = nghost+1, mitot-nghost
                if (alloc(iaddaux(1,i,j,locaux,mitot)) < 0.d0) then
                    ! B<0, zeta = h+B
                    zeta = alloc(iadd(1,i,j,loc,mitot)) &
                           + alloc(iaddaux(1,i,j,locaux,mitot))
                else
                    ! B>0, zeta = h
                    zeta = alloc(iadd(1,i,j,loc,mitot))
                endif
                if (mcapa > 0) then
                    zeta = zeta * alloc(iaddaux(mcapa,i,j,locaux,mitot)) 
                endif
                totmass = totmass + zeta
            enddo
        enddo           
        
        mptr = node(levelptr,mptr)  ! next patch
    enddo
 
    totmass = totmass * hx * hy
    if (time.eq. t0 .and. (level.eq.1) .and. .not. rest) then
        tmass0 = totmass
        write(6,*) 'Total zeta at initial time: ',tmass0
    endif
    write(outunit,77) time, totmass, totmass-tmass0
 77 format('time t = ',e12.5,',  total zeta = ',e22.15, '  diff = ', e11.4)
 
contains

    integer pure function iadd(ivar,i,j,loc,mitot)
        integer, intent(in) :: i, j, ivar, loc, mitot
        iadd = loc + ivar-1 + nvar*((j-1)*mitot+i-1)
    end function iadd

    integer pure function iaddaux(ivar,i,j,locaux,mitot)
        integer, intent(in) :: i, j, ivar, locaux, mitot
        iaddaux = locaux + ivar-1 + naux*((j-1)*mitot+i-1)
    end function iaddaux

end subroutine conck
