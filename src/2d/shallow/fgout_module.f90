module fgout_module

    implicit none
    save

    ! Container for fixed grid data, geometry and output settings
    type fgout_grid
        ! Grid data
        real(kind=8), pointer :: early(:,:,:)
        real(kind=8), pointer :: late(:,:,:)

        ! Geometry
        integer :: mx,my,point_style,fgno,output_format,output_style
        real(kind=8) :: dx,dy,x_low,x_hi,y_low,y_hi

        ! Time Tracking and output types
        integer :: num_output,next_output_index
        real(kind=8) :: start_time,end_time,dt

        integer, allocatable :: output_frames(:)
        real(kind=8), allocatable :: output_times(:)

        integer :: num_vars ! number of variables to interpolate (num_eqn+3)
        integer :: nqout ! number of q components to output (+1 for eta)
        logical, allocatable :: q_out_vars(:) ! which components to print

        !logical, allocatable :: iqout(:)  ! which components to output
        integer :: bathy_index,eta_index,time_index

    end type fgout_grid


    logical, private :: module_setup = .false.

    ! Fixed grid arrays and sizes
    integer :: FGOUT_num_grids
    type(fgout_grid), target, allocatable :: FGOUT_fgrids(:)
    real(kind=8) :: FGOUT_tcfmax
    real(kind=8), parameter :: FGOUT_ttol = 1.d-6 ! tolerance for times


contains


    ! Setup routine that reads in the fixed grids data file and sets up the
    ! appropriate data structures

    subroutine set_fgout(rest,num_eqn,fname)

        use amr_module, only: parmunit, tstart_thisrun
        use utility_module, only: parse_values

        implicit none

        ! Subroutine arguments
        logical, intent(in) :: rest  ! restart?
        integer, intent(in) :: num_eqn
        character(len=*), optional, intent(in) :: fname

        ! Local storage
        integer, parameter :: unit = 7
        integer :: i,k
        type(fgout_grid), pointer :: fg
        real(kind=8) :: ts

        integer :: n, iq
        real(kind=8) :: values(num_eqn+2)
        character(len=80) :: str

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
                fg => FGOUT_fgrids(i)
                ! Read in this grid's data
                read(unit,*) fg%fgno
                read(unit,*) fg%output_style
                read(unit,*) fg%num_output
                allocate(fg%output_times(fg%num_output))
                allocate(fg%output_frames(fg%num_output))

                if (fg%output_style == 1) then
                    read(unit,*) fg%start_time
                    read(unit,*) fg%end_time
                else if (fg%output_style == 2) then
                    read(unit,*) (fg%output_times(k), k=1,fg%num_output)
                    fg%start_time = fg%output_times(1)
                    fg%end_time = fg%output_times(fg%num_output)
                endif

                read(unit,*) fg%point_style
                read(unit,*) fg%output_format
                read(unit,*) fg%mx, fg%my
                read(unit,*) fg%x_low, fg%y_low
                read(unit,*) fg%x_hi, fg%y_hi

                allocate(fg%q_out_vars(num_eqn+2))  ! for q + eta, topo
                do iq=1,num_eqn+2
                    fg%q_out_vars(iq) = .false.
                enddo
                read(unit,'(a)') str
                call parse_values(str, n, values)
                do k=1,n
                    iq = nint(values(k))
                    fg%q_out_vars(iq) = .true.
                enddo
                write(6,*) '+++ q_out_vars:', fg%q_out_vars

                fg%num_vars = num_eqn + 3  ! also interp topo,eta,time
                fg%nqout = 0  ! count how many are to be output
                do k=1,num_eqn+2
                    if (fg%q_out_vars(k)) fg%nqout = fg%nqout + 1
                enddo
                !fg%nqout = fg%nqout + 1  ! for eta


                ! Initialize next_output_index
                ! (might be reset below in case of a restart)
                fg%next_output_index = 1

                if (fg%point_style .ne. 2) then
                    print *, 'set_fgout: ERROR, unrecognized point_style = ',\
                          fg%point_style
                endif

               if (fg%output_style == 1) then
                   ! Setup data for this grid
                   ! Set fg%dt (the timestep length between outputs)
                   if (fg%end_time <= fg%start_time) then
                       if (fg%num_output > 1) then
                          print *,'set_fgout: ERROR for fixed grid', i
                          print *,'start_time <= end_time yet num_output > 1'
                          print *,'set end_time > start_time or set num_output = 1'
                          stop
                       else
                           ! only a single fgout time:
                           fg%dt = 0.d0
                       endif
                   else
                       if (fg%num_output < 2) then
                           print *,'set_fgout: ERROR for fixed grid', i
                           print *,'end_time > start_time, yet num_output = 1'
                           print *,'set num_output > 2'
                           stop
                       else
                           fg%dt = (fg%end_time  - fg%start_time) &
                                               / (fg%num_output - 1)
                           do k=1,fg%num_output
                               fg%output_times(k) = fg%start_time + (k-1)*fg%dt
                           enddo
                       endif
                    endif
                endif

               do k=1,fg%num_output
                   if (rest) then
                       ! don't write initial time or earlier
                       ts = tstart_thisrun+FGOUT_ttol
                   else
                       ! do write initial time
                       ts = tstart_thisrun-FGOUT_ttol
                   endif

                   if (fg%output_times(k) < ts) then
                        ! will not output this time in this run
                        ! (might have already be done when restarting)
                        fg%output_frames(k) = -2
                        fg%next_output_index = k+1
                   else
                        ! will be reset to frameno when this is written
                        fg%output_frames(k) = -1
                   endif
               enddo



                ! Set spatial intervals dx and dy on each grid
                if (fg%mx > 1) then
                   !fg%dx = (fg%x_hi - fg%x_low) / (fg%mx - 1) ! points
                   fg%dx = (fg%x_hi - fg%x_low) / fg%mx   ! cells
                else if (fg%mx == 1) then
                   fg%dx = 0.d0
                else
                     print *,'set_fgout: ERROR for fixed grid', i
                     print *,'x grid points mx <= 0, set mx >= 1'
                endif

                if (fg%my > 1) then
                    !fg%dy = (fg%y_hi - fg%y_low) / (fg%my - 1) ! points
                    fg%dy = (fg%y_hi - fg%y_low) / fg%my  ! cells
                else if (fg%my == 1) then
                    fg%dy = 0.d0
                else
                    print *,'set_fgout: ERROR for fixed grid', i
                    print *,'y grid points my <= 0, set my >= 1'
                endif

                fg%bathy_index = num_eqn + 2
                fg%eta_index = num_eqn + 1
                fg%time_index = num_eqn + 3

                write(6,*) '+++ nqout = ',fg%nqout

                ! Allocate new fixed grid data arrays at early, late time:
                ! dimension (10,:,:) to work for either GeoClaw or D-Claw

                allocate(fg%early(10, fg%mx,fg%my))
                fg%early = nan()

                allocate(fg%late(10, fg%mx,fg%my))
                fg%late = nan()

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
                            maux,aux)

        use geoclaw_module, only: dry_tolerance
        implicit none

        ! Subroutine arguments
        integer, intent(in) :: fgrid_type
        type(fgout_grid), intent(inout) :: fgrid
        integer, intent(in) :: meqn,mxc,myc,mbc,maux
        real(kind=8), intent(in) :: t,dxc,dyc,xlowc,ylowc
        real(kind=8), intent(in) :: q(meqn,1-mbc:mxc+mbc,1-mbc:myc+mbc)
        real(kind=8), intent(in) :: aux(maux,1-mbc:mxc+mbc,1-mbc:myc+mbc)

        integer, parameter :: method = 0 ! interpolate in space?

        ! Indices
        integer :: ifg,jfg,m,ic1,ic2,jc1,jc2

        ! Tolerances
        real(kind=8) :: total_depth,depth_indicator,nan_check

        ! Geometry
        real(kind=8) :: xfg,yfg,xc1,xc2,yc1,yc2,xhic,yhic
        real(kind=8) :: geometry(4)

        real(kind=8) :: points(2,2), eta_tmp

        ! Work arrays for eta interpolation
        real(kind=8) :: eta(2,2),h(2,2)


        ! Alias to data in fixed grid
        integer :: num_vars
        real(kind=8), pointer :: fg_data(:,:,:)


        ! Setup aliases for specific fixed grid
        if (fgrid_type == 1) then
            num_vars = fgrid%num_vars
            fg_data => fgrid%early
        else if (fgrid_type == 2) then
            num_vars = fgrid%num_vars
            fg_data => fgrid%late
        else
            write(6,*) '*** Unexpected fgrid_type = ', fgrid_type
            stop
            ! fgrid_type==3 is deprecated, use fgmax grids instead
        endif

        xhic = xlowc + dxc*mxc
        yhic = ylowc + dyc*myc


        !write(59,*) '+++ ifg,jfg,eta,geometry at t = ',t

        ! Primary interpolation loops
        do ifg=1,fgrid%mx
            xfg=fgrid%x_low + (ifg-0.5d0)*fgrid%dx   ! cell centers
            do jfg=1,fgrid%my
                yfg=fgrid%y_low + (jfg-0.5d0)*fgrid%dy   ! cell centers

                ! Check to see if this coordinate is inside of this grid
                if (.not.((xfg < xlowc.or.xfg > xhic) &
                    .or.(yfg < ylowc.or.yfg > yhic))) then

                    ! find where xfg,yfg is in the computational grid and
                    ! compute the indices
                    !   (Note: may be subject to rounding error if fgout point
                    !    is right on a cell edge!)
                    ic1 = int((xfg - xlowc + dxc)/dxc)
                    jc1 = int((yfg - ylowc + dyc)/dyc)

                    if (method == 0) then

                        ! piecewise constant: take values from cell (ic1,jc1):

                        forall (m=1:meqn)
                            fg_data(m,ifg,jfg) = q(m,ic1,jc1)
                        end forall

                        fg_data(fgrid%bathy_index,ifg,jfg) = aux(1,ic1,jc1)

                        ! for pw constant we take B, h, eta from same cell,
                        ! so setting eta = h+B should be fine even near shore:
                        fg_data(fgrid%eta_index,ifg,jfg) = fg_data(1,ifg,jfg) &
                                + fg_data(fgrid%bathy_index,ifg,jfg)


                    else if (method == 1) then

                        ! bilinear used to interpolate to xfg,yfg
                        ! (not recommended)

                        ! define constant parts of bilinear
                        if (ic1.eq.mxc) ic1=mxc-1
                        if (jc1.eq.myc) jc1=myc-1
                        ic2 = ic1 + 1
                        jc2 = jc1 + 1

                        xc1 = xlowc + dxc * (ic1 - 0.5d0)
                        yc1 = ylowc + dyc * (jc1 - 0.5d0)
                        xc2 = xlowc + dxc * (ic2 - 0.5d0)
                        yc2 = ylowc + dyc * (jc2 - 0.5d0)

                        geometry = [(xfg - xc1) / dxc, &
                                    (yfg - yc1) / dyc, &
                                    (xfg - xc1) * (yfg - yc1) / (dxc*dyc), &
                                    1.d0]


                        ! Interpolate all conserved quantities and bathymetry
                        forall (m=1:meqn)
                            fg_data(m,ifg,jfg) = &
                                interpolate(q(m,ic1:ic2,jc1:jc2), geometry)
                        end forall

                        fg_data(fgrid%bathy_index,ifg,jfg) = &
                                interpolate(aux(1,ic1:ic2,jc1:jc2),geometry)


                        ! surface eta = h + B:

                        ! Note that for pw bilinear interp there may
                        ! be a problem interpolating each separately since
                        ! interpolated h + interpolated B may be much larger
                        ! than eta should be offshore.
                        eta = q(1,ic1:ic2,jc1:jc2) + aux(1,ic1:ic2,jc1:jc2)
                        fg_data(fgrid%eta_index,ifg,jfg) = interpolate(eta,geometry)
                        ! NEED TO FIX
                    endif


                    ! save the time this fgout point was computed:
                    fg_data(fgrid%time_index,ifg,jfg) = t


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
    subroutine fgout_write(fgrid,out_time,out_index)

        use geoclaw_module, only: dry_tolerance
        implicit none

        ! Subroutine arguments
        type(fgout_grid), intent(inout) :: fgrid
        real(kind=8), intent(in) :: out_time
        integer, intent(in) :: out_index

        ! I/O
        integer, parameter :: unit = 87
        character(len=15) :: fg_filename
        character(len=8) :: cfgno, cframeno
        character(len=8) :: file_format
        integer :: grid_number,ipos,idigit,out_number,columns
        integer :: ifg,ifg1, iframe,iframe1

        integer, parameter :: method = 0  ! interpolate in time?

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
                                           "i6,'                 nghost'/," // &
                                           "a10,'             format'/,/)"

        ! Other locals
        integer :: i,j,m,iq,k,jq,num_eqn
        real(kind=8) :: t0,tf,tau, qaug(10)
        real(kind=8), allocatable :: qeta(:,:,:)
        real(kind=4), allocatable :: qeta4(:,:,:)
        real(kind=8) :: h_early,h_late,topo_early,topo_late

        allocate(qeta(fgrid%nqout, fgrid%mx, fgrid%my))  ! to store h,hu,hv,eta
                                                         ! or subset

        ! Interpolate the grid in time, to the output time, using
        ! the solution in fgrid1 and fgrid2, which represent the
        ! solution on the fixed grid at the two nearest computational times
        do j=1,fgrid%my
            do i=1,fgrid%mx

                ! Check for small numbers
                forall(m=1:fgrid%num_vars-1,abs(fgrid%early(m,i,j)) < 1d-90)
                    fgrid%early(m,i,j) = 0.d0
                end forall

                if (method == 0) then

                    ! no interpolation in time, use solution from full step:
                    qaug = fgrid%early(:,i,j)

                    ! note that CFL condition ==> waves can't move more than 1
                    ! cell per time step on each level, so solution from nearest
                    ! full step should be correct to within a cell width
                    ! Better to use early than late since for refinement tracking
                    ! wave moving out into still water.

                else if (method == 1) then

                    ! interpolate in time. May have problems near shore?

                    ! Fetch times for interpolation, this is done per grid point
                    ! since each grid point may come from a different source
                    t0 = fgrid%early(fgrid%num_vars,i,j)
                    tf = fgrid%late(fgrid%num_vars,i,j)
                    tau = (out_time - t0) / (tf - t0)

                    ! check for small values:
                    forall(m=1:fgrid%num_vars-1,abs(fgrid%late(m,i,j)) < 1d-90)
                        fgrid%late(m,i,j) = 0.d0
                    end forall

                    ! linear interpolation:
                    qaug = (1.d0-tau)*fgrid%early(:,i,j) + tau*fgrid%late(:,i,j)

                    ! If resolution changed between early and late time, may be
                    ! problems near shore when interpolating B, h, eta
                    ! separately (at least in case when B changed and point
                    ! was dry at one time and wet the other).
                    ! Switch back to fgrid%early values, only in this case.
                    ! This is implemented below but not extensively tested.

                    if (qaug(1) > 0.d0) then
                        topo_early = fgrid%early(4,i,j)
                        topo_late = fgrid%late(4,i,j)
                        if (topo_early .ne. topo_late) then
                            ! AMR resolution changed between times
                            h_early = fgrid%early(1,i,j)
                            h_late = fgrid%late(1,i,j)
                            if (((h_early < dry_tolerance) &
                                    .and. (h_late >= dry_tolerance)) &
                                .or. ((h_late < dry_tolerance) &
                                    .and. (h_early >= dry_tolerance))) then
                                ! point changed between wet and dry
                                qaug = fgrid%early(:,i,j) ! don't interpolate
                            endif
                        endif
                    endif
                endif

                ! Output the requested quantities and eta value:

                iq = 1
                num_eqn = size(fgrid%q_out_vars) - 2  ! last components eta,B
                do jq=1,num_eqn+2
                    if (fgrid%q_out_vars(jq)) then
                        qeta(iq,i,j) = qaug(jq)
                        iq = iq+1
                    endif
                enddo

                ! now append eta value: (this is now included in loop above)
                !qeta(iq,i,j) = qaug(fgrid%eta_index)
                ! NOTE: qaug(fgrid%bathy_index) contains topo B value

            enddo
        enddo


        ! convert fgrid%fgno to a string of length at least 4 with leading 0's
        ! e.g. '0001' or '12345':
        write(cfgno, '(I0.4)') fgrid%fgno

        ! convert out_index to a string of length at least 4 with leading 0's
        write(cframeno, '(I0.4)') out_index


        fg_filename = 'fgout' // trim(cfgno) // '.q' // trim(cframeno)

        open(unit,file=fg_filename,status='unknown',form='formatted')


        write(unit,header_format) fgrid%fgno, 0, fgrid%mx,fgrid%my, &
            fgrid%x_low,fgrid%y_low, fgrid%dx, fgrid%dy

        if (fgrid%output_format == 1) then
            ! ascii output added to .q file:
            do j=1,fgrid%my
                do i=1,fgrid%mx
                    write(unit, "(50e26.16)") (qeta(k,i,j), k=1,fgrid%nqout)
                enddo
                write(unit,*) ' '  ! blank line required between rows
            enddo
        endif

        close(unit)

        if (fgrid%output_format == 3) then
            ! binary output goes in .b file as full 8-byte (float64):
            fg_filename = 'fgout' // trim(cfgno) // '.b' // trim(cframeno)
            open(unit=unit, file=fg_filename, status="unknown",    &
                 access='stream')
            write(unit) qeta
            close(unit)
        else if (fgrid%output_format == 2) then
            ! binary output goes in .b file as 4-byte (float32):
            fg_filename = 'fgout' // trim(cfgno) // '.b' // trim(cframeno)
            open(unit=unit, file=fg_filename, status="unknown",    &
                 access='stream')
            allocate(qeta4(fgrid%nqout, fgrid%mx, fgrid%my))
            qeta4 = real(qeta, kind=4)
            write(unit) qeta4
            deallocate(qeta4)
            close(unit)
        endif

        deallocate(qeta)

        ! time info .t file:


        if (fgrid%output_format == 1) then
            file_format = 'ascii'
        else if (fgrid%output_format == 2) then
            file_format = 'binary32'
        else if (fgrid%output_format == 3) then
            file_format = 'binary64'
        else
            write(6,*) '*** unrecognized fgrid%output_format = ', &
                        fgrid%output_format
            write(6,*) '*** should be ascii, binary32, or binary64'
        endif

        fg_filename = 'fgout' // trim(cfgno) // '.t' // trim(cframeno)
        open(unit=unit, file=fg_filename, status='unknown', form='formatted')
        ! time, num_eqn+1, num_grids, num_aux, num_dim, num_ghost:
        write(unit, t_file_format) out_time, fgrid%nqout, 1, 0, 2, 0,file_format
        close(unit)

        print "(a,i4,a,i4,a,e18.8)",'Writing fgout grid #',fgrid%fgno, &
              '  frame ',out_index,' at time =',out_time

    end subroutine fgout_write


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
        integer :: icell, jcell

        ! pw bilinear
        ! This is set up as a dot product between the approrpriate terms in
        ! the input data.  This routine could be vectorized or a BLAS routine
        ! used instead of the intrinsics to ensure that the fastest routine
        ! possible is being used
        interpolant = sum([points(2,1)-points(1,1), &
                       points(1,2)-points(1,1), &
                       points(1,1) + points(2,2) - (points(2,1) + points(1,2)), &
                       points(1,1)] * geometry)

    end function interpolate


end module fgout_module
