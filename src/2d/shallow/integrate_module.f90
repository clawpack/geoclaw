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
!  Defines raster data types. A single raster_collection is a set of 
!  rasters/DEMs that define a particular field (eg. topo, component of q or aux)
!  a single raster_type is raster data input (eg. a topo file)
! ============================================================================
module integrate_module

    
    implicit none

    
    type raster_data_type
        ! data for a set of rasters/DEMs for initializing grids (topo, qinit, aux)
        ! Work array for all data for a particular field for all t
        real(kind=8), allocatable :: data(:,:)

        !raster data 
        character(len=150), :: fname
        
        real(kind=8)  :: xlow, ylow, tlow, xhi, yhi, thi
        real(kind=8)  :: dx, dy
        real(kind=8)  :: rastertime

        integer :: mx, my, morder 
        integer :: minlevel, maxlevel 
        integer :: ifiletype, ID, save0data, iq, iaux, ifieldID
        logical :: finalized ! for this raster grid
    end type raster_data_type

    type raster_collection_type
        ! data for a set of rasters/DEMs for initializing grids (topo, qinit, aux)
        ! All arrays for all data for a particular field for all t
        type(raster_data_type), allocatable :: rasters(:)
        integer, allocatable :: raster_IDs(:)
        integer, allocatable :: raster_resolution_order(:)

        integer, :: num_raster_sets
        integer, :: iq, iaux, ifieldID
        logical :: finalized ! for all raster grids in collection
    end type raster_collection_type

contains

    subroutine cellrasterintegrate(cellint,xim,xcell,xip,yjm,ycell,
     &           yjp,raster_collection)

        implicit none

        ! input
        real(kind=8), intent(in) :: xim,xcell,xip,yjm,ycell,yjp
        type(raster_collection_type), intent(in) :: raster_collection
        ! output
        real(kind=8), intent(out) :: cellint

        !local 
        integer :: mrasters,rasternum
        type(raster_data_type) :: raster_data

        !##############################################################################
        ! cellrasterintegrate integrates a unique surface, 
        ! defined from data from a raster collection,
        ! (using the finest data available in any region)
        ! over a cell
        !
        ! The rectangle has coords:
        ! xim <= x <= xip, yjm <= y <= yjp, 
        ! with center (x,y) = (xcell, ycell)

        !The intersection (with one particular raster data) has coords:
        !xintlo <= x <= xinthi, yintlo <= y <= yinthi

        !initialize the integral of the surface
        cellint=0.d0

        !determine the type of integration
        im = 1

        ! first see if the grid cell is entirely in a fine topofile
        ! if so we are done...if not use recursive strategy for grids 
        !as coarse or coarser than first finest intersection 

        mrasters = raster_collection%num_raster_sets
        do m = 1,mrasters
            !look at raster files, from fine to coarse
            rasternum = raster_collection%raster_resolution_order(m))
            raster_data = raster_collection(rasternum)
            !check for intersection of grid cell and this raster set
            cellarea = (xip-xim)*(yjp-yjm)
            call intersection(indicator,area,xmlo,xmhi,ymlo,ymhi,xim,xip,yjm,yjp,raster_data)
            if (indicator.eq.1) then !cell overlaps raster set
                if (area.eq.cellarea) then !cell is entirely in raster set
                ! (should we check if they agree to some tolerance??)
                    !integrate surface and get out of here
                    cellint = cellint + rasterintegral(xmlo,xmhi,ymlo,ymhi,raster_data,im)
                    return
                else
                    go to 222
                endif
            endif
        enddo


 222  continue

      ! this grid cell intersects only raster set m and perhaps coarser:
      call rectintegral(xim,xip,yjm,yjp,m,cellint,raster_collection)

      return
      end

    end subroutine cellgridintegrate

    !-------------------------------------------------------------------------------

    function rasterintegral(xim,xip,yjm,yjp,intmethod,raster_data)
    !===========================================================================

        !######################################################################
        ! rasterintegral integrates a surface over a cell
        ! that is the intersection with a single raster data set
        ! the surface integrated is defined by a piecewise bilinear through the
        ! nodes of the raster data set

        ! The rectangular intersection has coords:
        ! xim <= x <= xip, yjm <= y <= yjp
        !
        ! The Cartesian grid has coords:
        ! xxlow <= x <= xxhi, yylow <= y <= yyhi, with grid cell size dxx by dyy
        ! and mxx by myy cells.
        !
        !                                written by David L. George
        !                                Seattle, WA 7/16/08
        !###########################################################################

        implicit none

        ! input
        real(kind=8), intent(in) :: xim,xip,yjm,yjp
        integer, intent(in) :: intmethod
        type(raster_data_type) :: raster_data
    
        
        ! output
        real(kind=8), intent(out) :: cellint

        ! local
        real(kind=8) :: theintegral,xxhi,xxlow,yyhi,yylow,dxx,dyy,dx,dy
        real(kind=8) :: z11,z22,z12,z21
        integer, :: mxx,myy
        integer, :: jj,jjstart,jjend
        integer, :: ii,iistar,iiend
        integer, :: djjstart,djjend,diistart,diiend
        integer, :: jjz1,jjz2

        xxlow = raster_data%xlow
        yylow = raster_data%ylow
        dxx = raster_data%dx 
        dyy = raster_data%dy
        mxx = raster_data%mx
        myy = raster_data%my

        !initialize:
        theintegral = 0.d0

        xxhi=xxlow+(mxx-1)*dxx
        yyhi=yylow+(myy-1)*dyy

        !TEST FOR SMALL ROUNDING ERROR==========
        if ((xim-xxlow).lt.0.d0.or.(xip-xxhi).gt.0.d0) then
            xim=max(xxlow,xim)
            xip=min(xxhi,xip)
        endif

        if ((yjm-yylow).lt.0.d0.or.(yjp-yyhi).gt.0.d0) then
            yjm=max(yylow,yjm)
            yjp=min(yyhi,yjp)
        endif

        dx=xip-xim
        dy=yjp-yjm

        !INTEGRATE PIECEWISE BILINEAR OVER RECTANGULAR REGION====
        if (intmethod.eq.1) then !use bilinear method

            !find indices that include the rectangular region
            djjstart=(yjm-yylow)/dyy
            jjstart=idint(djjstart)+1

            diistart=(xim-xxlow)/dxx
            iistart=idint(diistart)+1

            diiend=(xip-xxlow)/dxx
            iiend=ceiling(diiend) + 1

            djjend=(yjp-yylow)/dyy
            jjend=ceiling(djjend)+1

            iistart=max(iistart,1)
            jjstart=max(jjstart,1)
            iiend=min(mxx,iiend)
            jjend=min(myy,jjend)


            do jj=jjstart,jjend-1
                y1=yylow + (jj-1.d0)*dyy
                y2=yylow + (jj)*dyy
                !the array zz is indexed from north to south: 
                ! jjz is the actual index of interest in the array zz
                jjz1= myy-jj+1
                jjz2= jjz1-1

                do ii=iistart,iiend-1
                    x1=xxlow + (ii-1.d0)*dxx
                    x2=xxlow + (ii)*dxx

                    z11 = raster_data%data(ii,jjz1)
                    z12 = raster_data%data(ii,jjz2)
                    z21 = raster_data%data(ii+1,jjz1)
                    z22 = raster_data%data(ii+1,jjz2)

                    if (coordinate_system.eq.1) then !cartesian rectangle
                        theintegral = theintegral + &
                        bilinearintegral(xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,z11,z12,z21,z22)
                    elseif (coordinate_system.eq.2) then !integrate on surface of sphere
                        theintegral = theintegral + &
                        bilinearintegral_s(xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy,z11,z12,z21,z22)
                    else
                        write(*,*)  'TOPOINTEGRAL: coordinate_system error'
                    endif
                enddo
            enddo

        else
            write(*,*) 'RASTERINTEGRAL: only intmethod = 1,2 is supported'
        endif

        topointegral= theintegral
        return
        end function rasterintegral

    
recursive subroutine topoarea(x1,x2,y1,y2,m,area,raster_collection)

    ! Compute the area of overlap of top with the rectangle (x1,x2) x (y1,y2)
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

    type(raster_collection_type), intent(in) :: raster_collection

    ! local
    real(kind=8) :: xmlo,xmhi,ymlo,ymhi,x1m,x2m,y1m,y2m,area1,area2,area_m
    integer :: indicator  

    mrasters = raster_collection%num_raster_sets
    rasternum = raster_collection%raster_resolution_order(m)
    raster_data = raster_collection(rasternum)

    if (m == mrasters) then
        ! innermost step of recursion reaches this point.
        ! only using coarsest topo grid -- compute directly...
        call intersection(indicator,area,xmlo,xmhi,ymlo,ymhi,x1,x2,y1,y2,raster_data)
    else
        ! recursive call to compute area using one fewer topo grids:
        call topoarea(x1,x2,y1,y2,m+1,area1,raster_collection)

        ! region of intersection of cell with new topo grid:
        call intersection(indicator,area_m,x1m,x2m,y1m,y2m,x1,x2,y1,y2,raster_data)
             
        if (area_m > 0) then
        
            ! correction to subtract out from previous set of topo grids:
            call topoarea(x1m,x2m,y1m,y2m,m+1,area2,raster_collection)
    
            ! adjust integral due to corrections for new topo grid:
            area = area1 - area2 + area_m
        else
            area = area1
        endif
    endif

end subroutine topoarea

    

recursive subroutine rectintegral(x1,x2,y1,y2,m,integral,raster_collection)

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
    real(kind=8) :: xmlo,xmhi,ymlo,ymhi,x1m,x2m,y1m,y2m,int1,int2,int3,area
    integer :: indicator

    mrasters = raster_collection%num_raster_sets
    rasternum = raster_collection%raster_resolution_order(m)
    raster_data = raster_collection(rasternum)

    if (m == mrasters) then
         ! innermost step of recursion reaches this point.
         ! only using coarsest topo grid -- compute directly...
         call intersection(indicator,area,xmlo,xmhi,ymlo,ymhi,x1,x2,y1,y2,raster_data)

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

end module integrate_module
