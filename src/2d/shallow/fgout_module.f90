module fgout_module

    implicit none
    save

    ! Container for fixed grid data, geometry and output settings
    type fgout_grid
        ! Grid data
        real(kind=8), pointer :: early(:,:,:)
        real(kind=8), pointer :: late(:,:,:)
        
        ! Geometry
        integer :: num_vars(2),mx,my,point_style,fgno,output_format
        real(kind=8) :: dx,dy,x_low,x_hi,y_low,y_hi
        
        ! Time Tracking and output types
        integer :: num_output,last_output_index
        real(kind=8) :: last_output_time,start_time,end_time,dt
    end type fgout_grid    


    logical, private :: module_setup = .false.

    ! Fixed grid arrays and sizes
    integer :: FGOUT_num_grids
    type(fgout_grid), allocatable :: FGOUT_fgrids(:)
    real(kind=8) :: FGOUT_tcfmax

contains
                        
    
    ! Setup routine that reads in the fixed grids data file and sets up the
    ! appropriate data structures
    subroutine set_fgout(fname)

        use amr_module, only: parmunit

        implicit none
        
        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname
        
        ! Local storage
        integer, parameter :: unit = 7
        integer :: i

        if (.not.module_setup) then

            write(parmunit,*) ' '
            write(parmunit,*) '--------------------------------------------'
            write(parmunit,*) 'SETFGOUT:'
            write(parmunit,*) '-----------'

            ! Open data file
            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,'fgout_grids.data')
            endif

            ! Read in data
            read(unit,'(i2)') FGOUT_num_grids
            write(parmunit,*) '  mfgrids = ',FGOUT_num_grids
            if (FGOUT_num_grids == 0) then
                write(parmunit,*) '  No fixed grids specified for output'
                return
            endif
            
            ! Allocate fixed grids (not the data yet though)
            allocate(FGOUT_fgrids(FGOUT_num_grids))

            ! Read in data for each fixed grid
            do i=1,FGOUT_num_grids
                ! Read in this grid's data
                read(unit,*) FGOUT_fgrids(i)%fgno
                read(unit,*) FGOUT_fgrids(i)%start_time
                read(unit,*) FGOUT_fgrids(i)%end_time
                read(unit,*) FGOUT_fgrids(i)%num_output
                read(unit,*) FGOUT_fgrids(i)%point_style
                read(unit,*) FGOUT_fgrids(i)%output_format
                read(unit,*) FGOUT_fgrids(i)%mx, FGOUT_fgrids(i)%my
                read(unit,*) FGOUT_fgrids(i)%x_low, FGOUT_fgrids(i)%y_low
                read(unit,*) FGOUT_fgrids(i)%x_hi, FGOUT_fgrids(i)%y_hi
                
                   
                if (FGOUT_fgrids(i)%point_style .ne. 2) then
                    print *, 'set_fgout: ERROR, unrecognized point_style = ',\
                          FGOUT_fgrids(i)%point_style
                endif
                    
               ! Setup data for this grid
               ! Set dtfg (the timestep length between outputs) for each grid
               if (FGOUT_fgrids(i)%end_time <= FGOUT_fgrids(i)%start_time) then
                   if (FGOUT_fgrids(i)%num_output > 1) then 
                      print *,'set_fgout: ERROR for fixed grid', i
                      print *,'start_time <= end_time yet num_output > 1'
                      print *,'set end_time > start_time or set num_output = 1'
                      stop
                   else
                       FGOUT_fgrids(i)%dt = 0.d0
                   endif
               else
                   if (FGOUT_fgrids(i)%num_output < 2) then
                       print *,'set_fgout: ERROR for fixed grid', i
                       print *,'end_time > start_time, yet num_output = 1'
                       print *,'set num_output > 2'
                       stop
                   else
                       FGOUT_fgrids(i)%dt = (FGOUT_fgrids(i)%end_time &
                                          - FGOUT_fgrids(i)%start_time) &
                                           / (FGOUT_fgrids(i)%num_output - 1)
                   endif
                endif

                ! Initialize last_output_time and index
                FGOUT_fgrids(i)%last_output_time = FGOUT_fgrids(i)%start_time &
                     - FGOUT_fgrids(i)%dt
                FGOUT_fgrids(i)%last_output_index = -1  ! so first output is 0

                ! Set spatial intervals dx and dy on each grid
                if (FGOUT_fgrids(i)%mx > 1) then
                   FGOUT_fgrids(i)%dx = (FGOUT_fgrids(i)%x_hi &
                        - FGOUT_fgrids(i)%x_low) / (FGOUT_fgrids(i)%mx - 1)
                else if (FGOUT_fgrids(i)%mx == 1) then
                   FGOUT_fgrids(i)%dx = 0.d0
                else
                     print *,'set_fgout: ERROR for fixed grid', i
                     print *,'x grid points mx <= 0, set mx >= 1'
                endif

                if (FGOUT_fgrids(i)%my > 1) then
                    FGOUT_fgrids(i)%dy = (FGOUT_fgrids(i)%y_hi &
                         - FGOUT_fgrids(i)%y_low) / (FGOUT_fgrids(i)%my - 1)
                else if (FGOUT_fgrids(i)%my == 1) then
                    FGOUT_fgrids(i)%dy = 0.d0
                else
                    print *,'set_fgout: ERROR for fixed grid', i
                    print *,'y grid points my <= 0, set my >= 1'
                endif 
           
                ! set the number of variables stored for each grid
                ! this should be (the number of variables you want to write out + 1)
                FGOUT_fgrids(i)%num_vars(1) = 6
                
                ! Allocate new fixed grid data array
                allocate(FGOUT_fgrids(i)%early(FGOUT_fgrids(i)%num_vars(1), &
                         FGOUT_fgrids(i)%mx,FGOUT_fgrids(i)%my))
                FGOUT_fgrids(i)%early = nan()
                
                allocate(FGOUT_fgrids(i)%late(FGOUT_fgrids(i)%num_vars(1), &
                         FGOUT_fgrids(i)%mx,FGOUT_fgrids(i)%my))
                FGOUT_fgrids(i)%late = nan()
                
           enddo
           close(unit)
           
           FGOUT_tcfmax=-1.d16

           module_setup = .true.
        end if

    end subroutine set_fgout
    
    
    !=====================FGOUT_INTERP=======================================
    ! This routine interpolates q and aux on a computational grid
    ! to an fgout grid not necessarily aligned with the computational grid
    ! using bilinear interpolation defined on computational grid
    !=======================================================================
    subroutine fgout_interp(fgrid_type,fgrid, &
                            t,q,meqn,mxc,myc,mbc,dxc,dyc,xlowc,ylowc, &
                            maux,aux,maxcheck)
    
        use geoclaw_module, only: dry_tolerance  
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: fgrid_type
        type(fgout_grid), intent(inout) :: fgrid
        integer, intent(in) :: meqn,mxc,myc,mbc,maux,maxcheck
        real(kind=8), intent(in) :: t,dxc,dyc,xlowc,ylowc
        real(kind=8), intent(in) :: q(meqn,1-mbc:mxc+mbc,1-mbc:myc+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mxc+mbc,1-mbc:myc+mbc)
    
        ! Indices
        integer :: ifg,jfg,m,ic1,ic2,jc1,jc2
        integer :: bathy_index,eta_index

        ! Tolerances
        real(kind=8) :: total_depth,depth_indicator,nan_check

        ! Geometry
        real(kind=8) :: xfg,yfg,xc1,xc2,yc1,yc2,xhic,yhic
        real(kind=8) :: geometry(4)
        
        ! Work arrays for eta interpolation
        real(kind=8) :: eta(2,2),h(2,2)
        
        
        ! Alias to data in fixed grid
        integer :: num_vars
        real(kind=8), pointer :: fg_data(:,:,:)
        
        if (maxcheck > 0) then
            write(6,*) '*** Unexpected maxcheck = ',maxcheck
            stop
        endif
        
        ! Setup aliases for specific fixed grid
        if (fgrid_type == 1) then
            num_vars = fgrid%num_vars(1)
            fg_data => fgrid%early
        else if (fgrid_type == 2) then
            num_vars = fgrid%num_vars(1)
            fg_data => fgrid%late
        else
            write(6,*) '*** Unexpected fgrid_type = ', fgrid_type
            stop
            ! fgrid_type==3 is deprecated, use fgmax grids instead
            !num_vars = fgrid%num_vars(2)
            !fg_data => fgrid%often
        endif
            
        xhic = xlowc + dxc*mxc  
        yhic = ylowc + dyc*myc    
        
        ! Find indices of various quantities in the fgrid arrays
        bathy_index = meqn + 1
        eta_index = meqn + 2
    
        !write(59,*) '+++ ifg,jfg,eta,geometry at t = ',t
    
        ! Primary interpolation loops 
        do ifg=1,fgrid%mx
            xfg=fgrid%x_low + (ifg-1)*fgrid%dx
            do jfg=1,fgrid%my
                yfg=fgrid%y_low + (jfg-1)*fgrid%dy
    
                ! Check to see if this coordinate is inside of this grid
                if (.not.((xfg < xlowc.or.xfg > xhic).or.(yfg < ylowc.or.yfg > yhic))) then
    
                    ! find where xfg,yfg is in the computational grid and compute the indices
                    ! and relevant coordinates of each corner
                    ic1 = int((xfg-(xlowc+0.5d0*dxc))/(dxc))+1
                    jc1 = int((yfg-(ylowc+0.5d0*dyc))/(dyc))+1
                    if (ic1.eq.mxc) ic1=mxc-1
                    if (jc1.eq.myc) jc1=myc-1 
                    ic2 = ic1 + 1
                    jc2 = jc1 + 1
                        
                    xc1 = xlowc + dxc * (ic1 - 0.5d0)
                    yc1 = ylowc + dyc * (jc1 - 0.5d0)
                    xc2 = xlowc + dxc * (ic2 - 0.5d0)
                    yc2 = ylowc + dyc * (jc2 - 0.5d0)
         
                    ! Calculate geometry of interpolant
                    ! interpolate bilinear used to interpolate to xfg,yfg
                    ! define constant parts of bilinear
                    geometry = [(xfg - xc1) / dxc, &
                                (yfg - yc1) / dyc, &
                                (xfg - xc1) * (yfg - yc1) / (dxc*dyc), &
                                1.d0]
        
                    ! Interpolate for all conserved quantities and bathymetry
                    forall (m=1:meqn)
                        fg_data(m,ifg,jfg) = interpolate([[q(m,ic1,jc1),q(m,ic1,jc2)], &
                                                          [q(m,ic2,jc1),q(m,ic2,jc2)]], geometry)
                    end forall
                    
                    fg_data(bathy_index,ifg,jfg) = interpolate([[aux(1,ic1,jc1),aux(1,ic1,jc2)], &
                                                                [aux(1,ic2,jc1),aux(1,ic2,jc2)]], geometry)

                    
                    ! Interpolate surface eta, only use wet eta points near the shoreline
                    eta(1,:) = [aux(1,ic1,jc1) + q(1,ic1,jc1), aux(1,ic1,jc2) + q(1,ic1,jc2)]
                    eta(2,:) = [aux(1,ic2,jc1) + q(1,ic2,jc1), aux(1,ic2,jc2) + q(1,ic2,jc2)]
                    h(1,:) = [q(1,ic1,jc1),q(1,ic1,jc2)]
                    h(2,:) = [q(1,ic2,jc1),q(1,ic2,jc2)]
                         
                    depth_indicator= min(h(1,1),h(1,2),h(2,1),h(2,2))
                    total_depth = sum(h)
    
                    ! We are near shoreline
                    if (depth_indicator < dry_tolerance .and. &
                        total_depth > 4.d0 * dry_tolerance) then
                        ! Check to see if each cell around fixed grid point is 
                        ! wet, if not re-balance
                        if (h(1,1) < dry_tolerance) then
                            eta(1,1) =  (h(1,2)*eta(1,2) &
                                       + h(2,1)*eta(2,1) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(1,2) < dry_tolerance) then
                            eta(1,2) =  (h(1,1)*eta(1,1) &
                                       + h(2,1)*eta(2,1) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(2,1) < dry_tolerance) then
                            eta(2,1) =  (h(1,1)*eta(1,1) &
                                       + h(1,2)*eta(1,2) &
                                       + h(2,2)*eta(2,2)) / total_depth
                        endif
                        if (h(2,2) < dry_tolerance) then
                            eta(2,2)=  (h(1,2)*eta(1,2) &
                                      + h(2,1)*eta(2,1) &
                                      + h(1,1)*eta(1,1)) / total_depth
                        endif            
                    endif
                    
                    !if (total_depth <= 4.d0*dry_tolerance) then
                    !    eta(2,2) = nan()
                    !endif
    
                    ! evaluate the interpolant
                    fg_data(eta_index,ifg,jfg) = interpolate(eta,geometry)
                    
                    if (total_depth <= 4.d0*dry_tolerance) then
                        ! surface eta = B topography
                        eta(2,2) = fg_data(bathy_index,ifg,jfg)
                    endif
                    
                    fg_data(num_vars,ifg,jfg) = t
                    
                    !write(59,*) '+++',ifg,jfg
                    !write(59,*) eta
                    !write(59,*) geometry

                    
                endif ! if fgout point is on this grid
            enddo ! fgout y-coordinate loop
        enddo ! fgout x-coordinte loop
    
    end subroutine fgout_interp
    

    !================ fgout_write ==========================================
    ! This routine interpolates in time and then outputs a grid at
    ! time=out_time
    !
    ! files now have the same format as frames produced by outval
    !=======================================================================
    subroutine fgout_write(grid_index,fgrid,out_time,out_index)

        implicit none
        
        ! Subroutine arguments
        type(fgout_grid), intent(inout) :: fgrid
        real(kind=8), intent(in) :: out_time
        integer, intent(in) :: grid_index,out_index
              
        ! I/O
        integer, parameter :: unit = 95
        character(len=15) :: fg_filename
        character(len=4) :: cfgno, cframeno
        integer :: grid_number,ipos,idigit,out_number,columns
        integer :: ifg,ifg1, iframe,iframe1
        
        ! Output format strings 
        ! These are now the same as in outval for frame data, for compatibility
        ! For fgout grids there is only a single grid (ngrids=1)
        ! and we set AMR_level=0, naux=0, nghost=0 (so no extra cells in binary)
        
        character(len=*), parameter :: header_format =                         &
                                    "(i6,'                 grid_number',/," // &
                                     "i6,'                 AMR_level',/,"   // &
                                     "i6,'                 mx',/,"          // &
                                     "i6,'                 my',/"           // &
                                     "e26.16,'    xlow', /, "               // &
                                     "e26.16,'    ylow', /,"                // &
                                     "e26.16,'    dx', /,"                  // &
                                     "e26.16,'    dy',/)"
                                     
        character(len=*), parameter :: t_file_format = "(e18.8,'    time', /," // &
                                           "i6,'                 meqn'/,"   // &
                                           "i6,'                 ngrids'/," // &
                                           "i6,'                 naux'/,"   // &
                                           "i6,'                 ndim'/,"   // &
                                           "i6,'                 nghost'/,/)"
        
        ! Other locals
        integer :: i,j,m
        real(kind=8) :: t0,tf,tau, qaug(6)
        real(kind=8), allocatable :: qeta(:)
    
        !allocate(qeta(4, fgrid%mx, fgrid%my))  ! to store h,hu,hv,eta
        allocate(qeta(4*fgrid%mx*fgrid%my))
        
        
        ! Interpolate the grid in time, to the output time, using 
        ! the solution in fgrid1 and fgrid2, which represent the 
        ! solution on the fixed grid at the two nearest computational times
        do j=1,fgrid%my
            do i=1,fgrid%mx
                ! Fetch times for interpolation, this is done per grid point 
                ! since each grid point may come from a different source
                t0 = fgrid%early(fgrid%num_vars(1),i,j)
                tf = fgrid%late(fgrid%num_vars(1),i,j)
                tau = (out_time - t0) / (tf - t0)
                
                ! Check for small numbers
                forall(m=1:fgrid%num_vars(1)-1,abs(fgrid%early(m,i,j)) < 1d-90)
                    fgrid%early(m,i,j) = 0.d0
                end forall
                forall(m=1:fgrid%num_vars(1)-1,abs(fgrid%late(m,i,j)) < 1d-90)
                    fgrid%late(m,i,j) = 0.d0
                end forall
                
                ! interpolate in time:
                qaug = (1.d0-tau)*fgrid%early(:,i,j) + tau*fgrid%late(:,i,j)
                !write(6,*) '+++ tau, early, late: ',tau,fgrid%early(:,i,j),fgrid%late(:,i,j)
                
                ! Output the conserved quantities and topo value
                !qeta(1,i,j) = qaug(1)  ! h
                !qeta(2,i,j) = qaug(2)  ! hu
                !qeta(3,i,j) = qaug(3)  ! hv
                !qeta(4,i,j) = qaug(5)  ! eta
                qeta(iaddqeta(1,i,j,fgrid%mx)) = qaug(1)
                qeta(iaddqeta(2,i,j,fgrid%mx)) = qaug(2)
                qeta(iaddqeta(3,i,j,fgrid%mx)) = qaug(3)
                qeta(iaddqeta(4,i,j,fgrid%mx)) = qaug(5)

            enddo
        enddo


        ! Make the file names and open output files
        cfgno = '0000'
        ifg = grid_index
        ifg1 = ifg
        do ipos=4,1,-1
            idigit = mod(ifg1,10)
            cfgno(ipos:ipos) = char(ichar('0') + idigit)
            ifg1 = ifg1/10
            enddo

        cframeno = '0000'
        iframe = out_index
        iframe1 = iframe
        do ipos=4,1,-1
            idigit = mod(iframe1,10)
            cframeno(ipos:ipos) = char(ichar('0') + idigit)
            iframe1 = iframe1/10
            enddo
            
        write(6,*) '+++ grid_index, out_index: ',grid_index, out_index
        write(6,*) '+++ cfgno, cframeno: ',cfgno, '    ', cframeno
        !fg_filename = 'fgout' // cfgno // 'frame' // cframeno // '.txt'
        fg_filename = 'fgout' // cfgno // '.q' // cframeno 
        print *, 'Writing to file ', fg_filename

        open(unit,file=fg_filename,status='unknown',form='formatted')

        ! Determine number of columns that will be written out
        columns = fgrid%num_vars(1) - 1
        if (fgrid%num_vars(2) > 1) then
           columns = columns + 2
        endif
        
        !write(6,*) '+++ fgout out_time = ',out_time
        !write(6,*) '+++ fgrid%num_vars: ',fgrid%num_vars(1),fgrid%num_vars(2)
        
        ! Write out header in .q file:
        !write(unit,header_format) out_time,fgrid%mx,fgrid%my, &
        !     fgrid%x_low,fgrid%y_low, fgrid%x_hi,fgrid%y_hi, columns

        write(unit,header_format) fgrid%fgno, 0, fgrid%mx,fgrid%my, &
            fgrid%x_low,fgrid%y_low, fgrid%dx, fgrid%dy
            
        if (fgrid%output_format == 1) then
            ! ascii output added to .q file:
            do j=1,fgrid%my
                do i=1,fgrid%mx
                    write(unit, "(50e26.16)") &
                         (qeta(iaddqeta(m,i,j,fgrid%mx)), m=1,4)
                    !qeta(1,i,j),qeta(2,i,j), &
                    !            qeta(3,i,j),qeta(4,i,j)
                enddo
                write(unit,*) ' '  ! blank line required between rows
            enddo
        endif  
        
        close(unit)
        
        if (fgrid%output_format == 3) then
            ! binary output goes in .b file:
            fg_filename = 'fgout' // cfgno // '.b' // cframeno 
            write(6,*) '+++ fgout filename: ',fg_filename
            open(unit=unit, file=fg_filename, status="unknown",    &
                 access='stream')
            write(unit) qeta
            close(unit)
        endif
        
        deallocate(qeta)

        ! time info .t file:
        fg_filename = 'fgout' // cfgno // '.t' // cframeno 
        write(6,*) '+++ fgout filename: ',fg_filename
        open(unit=unit, file=fg_filename, status='unknown', form='formatted')
        ! time, num_eqn+1, num_grids, num_aux, num_dim, num_ghost:
        write(unit, t_file_format) out_time, 4, 1, 0, 2, 0
        close(unit)
        
        print "(a,i2,a,i2,a,e18.8)",'fgout for grid #',grid_index, &
              '  frame ',out_index,' at time =',out_time
      
        ! Index into qeta for binary output
        ! Note that this implicitly assumes that we are outputting only h, hu, hv
        ! and will not output more (change num_eqn parameter above)
        
    end subroutine fgout_write

    pure integer function iaddqeta(m, i, j, mx)
        implicit none
        integer, intent(in) :: m, i, j, mx
        iaddqeta = 1 + m - 1 + 4 * ((j - 1) * mx + i - 1)
    end function iaddqeta
              
    
    ! =========================================================================
    ! Utility functions for this module
    ! Returns back a NaN
    real(kind=8) function nan()
        real(kind=8) dnan
        integer inan(2)
        equivalence (dnan,inan)
        inan(1)=2147483647
        inan(2)=2147483647
        nan=dnan
    end function nan
    
    ! Interpolation function (in space)
    ! Given 4 points (points) and geometry from x,y,and cross terms
    real(kind=8) pure function interpolate(points,geometry) result(interpolant)
                            
        implicit none
                                
        ! Function signature
        real(kind=8), intent(in) :: points(2,2)
        real(kind=8), intent(in) :: geometry(4)
        
        ! This is set up as a dot product between the approrpriate terms in 
        ! the input data.  This routine could be vectorized or a BLAS routine
        ! used instead of the intrinsics to ensure that the fastest routine
        ! possible is being used
        interpolant = sum([points(2,1)-points(1,1), &
                           points(1,2)-points(1,1), &
                           points(1,1) + points(2,2) - (points(2,1) + points(1,2)), &
                           points(1,1)] * geometry)
                           
    end function interpolate
    
    ! Interpolation function in time
    pure function interpolate_time(num_vars,early,late,tau) result(interpolant)
        
        implicit none
        
        ! Input arguments
        integer, intent(in) :: num_vars
        real(kind=8), intent(in) :: early(num_vars),late(num_vars),tau
        
        ! Return value
        real(kind=8) :: interpolant(num_vars)

        interpolant = (1.d0 - tau) * early(:) + tau * late(:)

    end function interpolate_time
    

end module fgout_module
