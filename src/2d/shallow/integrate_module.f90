! ============================================================================
!  Program:     integrate_module
!  File:        integrate_module.f90
!  Created:     2018-05-27
!  
! ============================================================================
!      
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
!  this module contains integration routines that can fill grids by 
!  integrating multiple raster data sets (ie. sets of DEM file formats)
!  It works for qinit, topo, and aux data types.
! ============================================================================
!  Module contains a few utilities for integrating.
!  Defines a raster set data type. A single raster_set is a set of 
!  rasters/DEMs that define a particular field (eg. topo, component of q or aux)
! ============================================================================
module integrate_module

    
    implicit none

    
    type raster_set
    ! data for a set of rasters/DEMs for initializing grids (topo, qinit, aux)
    ! Work array for all data for a particular field for all t
    real(kind=8), allocatable :: rasterwork(:)

    !raster data 
    character(len=150), allocatable :: rasterfname(:)
    integer :: mrasterfiles,mrastersize
    real(kind=8), allocatable :: xlowraster(:), ylowraster(:), tlowraster(:)
    real(kind=8), allocatable :: xhiraster(:), yhiraster(:), thiraster(:)
    real(kind=8), allocatable :: dxraster(:), dyraster(:)
    real(kind=8), allocatable :: rastertime(:)
    integer, allocatable ::  mxraster(:), myraster(:)

    integer, allocatable :: i0raster(:), mraster(:), mrasterorder(:)
    integer, allocatable :: minlevelraster(:), maxlevelraster(:), irasterfiletype(:)
    integer, allocatable :: rasterID(:),raster0save(:),irasterfield(:)
    logical :: raster_finalized
contains

    
recursive subroutine topoarea(x1,x2,y1,y2,m,area)

    ! Compute the area of overlap of topo with the rectangle (x1,x2) x (y1,y2)
    ! using topo arrays indexed mtopoorder(mtopofiles) through mtopoorder(m) 
    ! (coarse to fine).

    ! The main call to this subroutine has corners of a physical domain for
    ! the rectangle and m = 1 in order to compute the area of overlap of
    ! domain by all topo arrays.  Used to check inputs and insure topo
    ! covers domain.

    ! The recursive strategy is to first compute the area using only topo 
    ! arrays mtopoorder(mtopofiles) to mtopoorder(m+1), 
    ! and then apply corrections due to adding topo array mtopoorder(m).
     
    ! Corrections are needed if the new topo array intersects the grid cell.
    ! Let the intersection be (x1m,x2m) x (y1m,y2m).
    ! Two corrections are needed, first to subtract out the area over
    ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
    ! topo arrays mtopoorder(mtopofiles) to mtopoorder(m+1),
    ! and then adding in the area over this same region using 
    ! topo array mtopoorder(m).

    ! Based on the recursive routine rectintegral that integrates
    ! topo over grid cells using a similar strategy.

    implicit none

    ! arguments
    real (kind=8), intent(in) :: x1,x2,y1,y2
    integer, intent(in) :: m
    real (kind=8), intent(out) :: area

    ! local
    real(kind=8) :: xmlo,xmhi,ymlo,ymhi,x1m,x2m, &
        y1m,y2m, area1,area2,area_m
    integer :: mfid, indicator, i0
    real(kind=8), external :: topointegral  


    mfid = mtopoorder(m)
    i0=i0topo(mfid)

    if (m == mtopofiles) then
         ! innermost step of recursion reaches this point.
         ! only using coarsest topo grid -- compute directly...
         call intersection(indicator,area,xmlo,xmhi, &
             ymlo,ymhi, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

    else
        ! recursive call to compute area using one fewer topo grids:
        call topoarea(x1,x2,y1,y2,m+1,area1)

        ! region of intersection of cell with new topo grid:
        call intersection(indicator,area_m,x1m,x2m, &
             y1m,y2m, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

        
        if (area_m > 0) then
        
            ! correction to subtract out from previous set of topo grids:
            call topoarea(x1m,x2m,y1m,y2m,m+1,area2)
    
            ! adjust integral due to corrections for new topo grid:
            area = area1 - area2 + area_m
        else
            area = area1
        endif
    endif

end subroutine topoarea

    

recursive subroutine rectintegral(x1,x2,y1,y2,m,integral)

    ! Compute the integral of topo over the rectangle (x1,x2) x (y1,y2)
    ! using topo arrays indexed mtopoorder(mtopofiles) through mtopoorder(m) 
    ! (coarse to fine).

    ! The main call to this subroutine has corners of a grid cell for the 
    ! rectangle and m = 1 in order to compute the integral over the cell 
    ! using all topo arrays.

    ! The recursive strategy is to first compute the integral using only topo 
    ! arrays mtopoorder(mtopofiles) to mtopoorder(m+1), 
    ! and then apply corrections due to adding topo array mtopoorder(m).
     
    ! Corrections are needed if the new topo array intersects the grid cell.
    ! Let the intersection be (x1m,x2m) x (y1m,y2m).
    ! Two corrections are needed, first to subtract out the integral over
    ! the rectangle (x1m,x2m) x (y1m,y2m) computed using
    ! topo arrays mtopoorder(mtopofiles) to mtopoorder(m+1),
    ! and then adding in the integral over this same region using 
    ! topo array mtopoorder(m).

    ! Note that the function topointegral returns the integral over the 
    ! rectangle based on a single topo array, and that routine calls
    ! bilinearintegral.


    implicit none

    ! arguments
    real (kind=8), intent(in) :: x1,x2,y1,y2
    integer, intent(in) :: m
    real (kind=8), intent(out) :: integral

    ! local
    real(kind=8) :: xmlo,xmhi,ymlo,ymhi,area,x1m,x2m, &
        y1m,y2m, int1,int2,int3
    integer :: mfid, indicator, mp1fid, i0
    real(kind=8), external :: topointegral  


    mfid = mtopoorder(m)
    i0=i0topo(mfid)

    if (m == mtopofiles) then
         ! innermost step of recursion reaches this point.
         ! only using coarsest topo grid -- compute directly...
         call intersection(indicator,area,xmlo,xmhi, &
             ymlo,ymhi, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

         if (indicator.eq.1) then
            ! cell overlaps the file
            ! integrate surface over intersection of grid and cell
            integral = topointegral( xmlo,xmhi,ymlo, &
                    ymhi,xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                    dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
         else
            integral = 0.d0
         endif

    else
        ! recursive call to compute area using one fewer topo grids:
        call rectintegral(x1,x2,y1,y2,m+1,int1)

        ! region of intersection of cell with new topo grid:
        call intersection(indicator,area,x1m,x2m, &
             y1m,y2m, x1,x2,y1,y2, &
             xlowtopo(mfid),xhitopo(mfid),ylowtopo(mfid),yhitopo(mfid))

        
        if (area > 0) then
        
            ! correction to subtract out from previous set of topo grids:
            call rectintegral(x1m,x2m,y1m,y2m,m+1,int2)
    
            ! correction to add in for new topo grid:
            int3 = topointegral(x1m,x2m, y1m,y2m, &
                        xlowtopo(mfid),ylowtopo(mfid),dxtopo(mfid), &
                        dytopo(mfid),mxtopo(mfid),mytopo(mfid),topowork(i0),1)
    
            ! adjust integral due to corrections for new topo grid:
            integral = int1 - int2 + int3
        else
            integral = int1
        endif
    endif

end subroutine rectintegral


subroutine intersection(indicator,area,xintlo,xinthi, &
           yintlo,yinthi,x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi)

    ! find the intersection of two rectangles, return the intersection
    ! and it's area, and indicator =1
    ! if there is no intersection, indicator =0

      implicit none

      integer, intent(out) :: indicator

      real(kind=8), intent(in) ::  x1lo,x1hi,y1lo,y1hi,x2lo,x2hi,y2lo,y2hi
      real(kind=8), intent(out) :: area,xintlo,xinthi,yintlo,yinthi

      xintlo=dmax1(x1lo,x2lo)
      xinthi=dmin1(x1hi,x2hi)
      yintlo=dmax1(y1lo,y2lo)
      yinthi=dmin1(y1hi,y2hi)


      if (xinthi.gt.xintlo.and.yinthi.gt.yintlo) then
         area = (xinthi-xintlo)*(yinthi-yintlo)
         indicator = 1
      else
         area = 0.d0
         indicator = 0
      endif

end subroutine intersection

end module topo_module
