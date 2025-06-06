! ============================================================================
!  Module for topography data
! ============================================================================
module topo_module

    use amr_module, only: xlower,xupper,ylower,yupper,hxposs,hyposs,nghost
    implicit none

    logical, private :: module_setup = .false.

    ! Work array for topography for all t
    real(kind=8), allocatable :: topowork(:)

    ! Topography file data
    integer :: test_topography
    character(len=150), allocatable :: topofname(:)
    integer :: mtopofiles
    integer(kind=8) :: mtoposize
    real(kind=8), allocatable :: xlowtopo(:), ylowtopo(:), tlowtopo(:)
    real(kind=8), allocatable :: xhitopo(:), yhitopo(:), thitopo(:)
    real(kind=8), allocatable :: dxtopo(:), dytopo(:)
    real(kind=8), allocatable :: topotime(:)
    integer(kind=8), allocatable :: mtopo(:), i0topo(:)

    integer, allocatable :: mxtopo(:), mytopo(:), mtopoorder(:)
    integer, allocatable :: itopotype(:)
    integer, allocatable :: topoID(:),topo0save(:)
    logical :: topo_finalized

    ! Moving topography support
    integer :: imovetopo, aux_finalized

    ! Analytic topography
    real(kind=8), private :: topo_left, topo_right, topo_location
    real(kind=8), private :: topo_x0,topo_x1,topo_x2,topo_basin_depth
    real(kind=8), private :: topo_shelf_depth,topo_shelf_slope,topo_beach_slope

    ! NetCDF4 support
    real(kind=4), parameter :: CONVENTION_REQUIREMENT = 1.0

    ! dtopo variables
    ! Work array
    real(kind=8), allocatable :: dtopowork(:)
    real(kind=8) :: dt_max_dtopo

    ! File data parameters
    character (len=150), allocatable :: dtopofname(:)
    real(kind=8), allocatable :: xlowdtopo(:),ylowdtopo(:),xhidtopo(:)
    real(kind=8), allocatable :: yhidtopo(:),t0dtopo(:),tfdtopo(:)
    real(kind=8), allocatable :: dxdtopo(:),dydtopo(:),dtdtopo(:)
    real(kind=8), allocatable :: tdtopo1(:),tdtopo2(:),taudtopo(:)

    integer, allocatable :: mxdtopo(:), mydtopo(:), mtdtopo(:)
    integer(kind=8), allocatable :: mdtopo(:), i0dtopo(:)
    integer, allocatable :: dtopotype(:), mdtopoorder(:),kdtopo1(:),kdtopo2(:)
    integer(kind=8), allocatable :: index0_dtopowork1(:),index0_dtopowork2(:)

    integer :: num_dtopo
    real(kind=8) dz
    logical, allocatable :: topoaltered(:) !don't think this is needed anymore

    ! Initial topography
    ! Work array for initial topography (only arrays where topo evolves)
    real(kind=8), allocatable :: topo0work(:)
    integer, allocatable :: topo0ID(:)
    integer(kind=8), allocatable :: i0topo0(:)
    integer(kind=8) :: mtopo0size
    integer :: mtopo0files

    real(kind=8) topo_missing

contains

    ! ========================================================================
    ! Read topography files as specified in topography.data
    !
    ! Each topography file has a type stored in topotype(i).
    !   topotype = 1:  standard GIS format: 3 columns: lon,lat,height(m)
    !   topotype = 2:  Header as in DEM file, height(m) one value per line
    !   topotype = 3:  Header as in DEM file, height(m) one row per line
    ! For other formats modify readtopo routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Finest value of topography in a given region will be used for
    ! computation
    ! ========================================================================
    subroutine read_topo_settings(restart,file_name)

        use geoclaw_module

        implicit none

        ! Input arguments
        character(len=*), intent(in), optional :: file_name
        logical, intent(in) :: restart

        ! Locals
        integer, parameter :: iunit = 7
        integer :: i,j,finer_than,rank
        real(kind=8) :: area_i,area_j
        real(kind=8) :: area, area_domain
        if (.not.module_setup) then

            ! Open and begin parameter file output
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SETTOPO:'
            write(GEO_PARM_UNIT,*) '---------'


            if (present(file_name)) then
                call opendatafile(iunit, file_name)
            else
                call opendatafile(iunit, 'topo.data')
            endif
            ! Read in value to use in place of no_data_value in topofile
            read(iunit,*) topo_missing

            ! Read in topography specification type
            read(iunit,"(i1)") test_topography

            ! Primary topography type, read in topography files specified
            if (test_topography == 0) then
                read(iunit,*) mtopofiles

                if (mtopofiles == 0) then
                    write(GEO_PARM_UNIT,*) '   mtopofiles = 0'
                    write(GEO_PARM_UNIT,*) '   No topo files specified, '
                    write(GEO_PARM_UNIT,*) '          will set B(x,y) = 0 in setaux'
                    return
                endif

                mtopofiles = mtopofiles + num_dtopo
                write(GEO_PARM_UNIT,*) '   mtopofiles = ',mtopofiles-num_dtopo
                
                ! Read and allocate data parameters for each file
                allocate(mxtopo(mtopofiles),mytopo(mtopofiles))
                allocate(xlowtopo(mtopofiles),ylowtopo(mtopofiles))
                allocate(tlowtopo(mtopofiles),xhitopo(mtopofiles),yhitopo(mtopofiles))
                allocate(thitopo(mtopofiles),dxtopo(mtopofiles),dytopo(mtopofiles))
                allocate(topofname(mtopofiles),itopotype(mtopofiles))
                allocate(i0topo(mtopofiles),mtopo(mtopofiles),mtopoorder(mtopofiles))
                allocate(topoID(mtopofiles),topotime(mtopofiles),topo0save(mtopofiles))
                allocate(i0topo0(mtopofiles),topo0ID(mtopofiles))

                do i=1,mtopofiles - num_dtopo
                    read(iunit,*) topofname(i)
                    read(iunit,*) itopotype(i)

                    write(GEO_PARM_UNIT,*) '   '
                    write(GEO_PARM_UNIT,*) '   ',topofname(i)
                    write(GEO_PARM_UNIT,*) '  itopotype = ', itopotype(i)
                    if (abs(itopotype(i)) == 1) then
                        print *, 'converting to topotype > 1 might reduce file size'
                        print *, 'python tools for converting files are provided'
                    endif
                    call read_topo_header(topofname(i),itopotype(i),mxtopo(i), &
                        mytopo(i),xlowtopo(i),ylowtopo(i),xhitopo(i),yhitopo(i), &
                        dxtopo(i),dytopo(i))
                    topoID(i) = i
                    mtopo(i) = int(mxtopo(i), 8) * int(mytopo(i), 8)
                enddo

                ! adding extra topo arrays corresponding spatially to dtopo
                ! arrays will be part of time dependent topowork as with all others.
                ! the state of these arrays at t=0 will be stored in topo0work
                topo0save(:)= 0
                do i= mtopofiles - num_dtopo + 1, mtopofiles
                   j = i - mtopofiles + num_dtopo
                   itopotype(i) = dtopotype(j)
                   mxtopo(i) = mxdtopo(j)
                   mytopo(i) = mydtopo(j)
                   xlowtopo(i) = xlowdtopo(j)
                   ylowtopo(i) = ylowdtopo(j)
                   xhitopo(i) = xhidtopo(j)
                   yhitopo(i) = yhidtopo(j)
                   dxtopo(i) = dxdtopo(j)
                   dytopo(i) = dydtopo(j)
                   mtopo(i) = int(mxtopo(i), 8) * int(mytopo(i), 8)
                   topoID(i) = i
                   topotime(i) = -huge(1.0)
                   topo0save(i) = 1
                enddo

                ! Indexing into work array
                i0topo(1)=1
                if (mtopofiles > 1) then
                    do i=2,mtopofiles
                        i0topo(i)=i0topo(i-1) + mtopo(i-1)
                    enddo
                endif

                ! Read topography and allocate space for each file
                mtoposize = sum(mtopo)
                allocate(topowork(mtoposize))

                do i=1,mtopofiles - num_dtopo
                    topoID(i) = i
                    topotime(i) = -huge(1.0)
                    call read_topo_file(mxtopo(i),mytopo(i),itopotype(i),topofname(i), &
                        xlowtopo(i),ylowtopo(i),topowork(i0topo(i):i0topo(i)+mtopo(i)-1))
                    ! set topo0save(i) = 1 if this topo file intersects any
                    ! dtopo file.  This approach to setting topo0save is changed from 
                    ! v5.4.1, where it only checked if some dtopo point lies within the
                    ! topo grid, which might not happen for small scale topo
                    do j=mtopofiles - num_dtopo + 1, mtopofiles
                        if ((xhitopo(i)<xlowtopo(j)) .or. &
                            (xlowtopo(i)>xhitopo(j)) .or. &
                            (yhitopo(i)<ylowtopo(j)) .or. &
                            (ylowtopo(i)>yhitopo(j))) then
                              topo0save(i) = 0
                          else
                              topo0save(i) = 1
                          endif

                    enddo
                enddo

                ! topography order...This determines which order to process topography
                !
                ! The finest topography will be given priority in any region
                ! mtopoorder(rank) = i means that i'th topography file has rank rank,
                ! where the file with rank=1 is the finest and considered first.
                do i=1,mtopofiles
                    finer_than = 0
                    do j=1,mtopofiles
                        if (j /= i) then
                            area_i=dxtopo(i)*dytopo(i)
                            area_j=dxtopo(j)*dytopo(j)
                            if (area_i < area_j) finer_than = finer_than + 1
                            ! if two files have the same resolution, order is
                            ! arbitrarily chosen
                            if ((area_i == area_j).and.(j < i)) then
                                finer_than = finer_than + 1
                            endif
                        endif
                    enddo
                    ! ifinerthan tells how many other files, file i is finer than
                    rank = mtopofiles - finer_than
                    mtopoorder(rank) = i
                enddo

                write(GEO_PARM_UNIT,*) ' '
                write(GEO_PARM_UNIT,*) '  Ranking of topography files', &
                    '  finest to coarsest: ', &
                    (mtopoorder(rank),rank=1,mtopofiles)
                write(GEO_PARM_UNIT,*) ' '

                !set values in topo array for dtopo generated topo
                !this call will also determine which topo arrays to save in topo0work
                do i = mtopofiles - num_dtopo + 1, mtopofiles
                   call set_topo_for_dtopo(mxtopo(i),mytopo(i),dxtopo(i),dytopo(i), &
                        xlowtopo(i),yhitopo(i), topowork(i0topo(i):i0topo(i)+mtopo(i)-1))
                enddo

                !create topo0work array for finest arrays covering dtopo
                !arrays to be saved are indicated in topo0save
                topo_finalized = .true.
                aux_finalized = 2   !# indicates aux arrays properly set with dtopo
                if (num_dtopo>0) then
                   topo_finalized = .false.
                   if (.not. restart) then ! rest is read in
                      aux_finalized = 0  !# will be incremented each time level 1 goes
                   endif

                   i0topo0(1) = 1
                   mtopo0size = dot_product(mtopo,topo0save)
                   allocate(topo0work(mtopo0size))
                   do i = 2,mtopofiles
                      i0topo0(i)= i0topo0(i-1) + mtopo(i-1)*topo0save(i-1)
                   enddo

                   do i = 1,mtopofiles
                      if (topo0save(i)>0) then
                         topo0work(i0topo0(i):i0topo0(i)+mtopo(i)-1) = &
                            topowork(i0topo(i):i0topo(i)+mtopo(i)-1)
                      endif
                   enddo
                endif

                ! Check that topo arrays cover full domain:
                call topoarea(xlower,xupper,ylower,yupper,1,area)
                area_domain = (yupper-ylower)*(xupper-xlower)
                if (abs(area - area_domain) > 1e-2*area_domain) then
                    write(6,*) '**** topo arrays do not cover domain'
                    write(6,*) '**** area of overlap = ', area
                    write(6,*) '**** area of domain  = ', area_domain
                    stop
                else if (abs(area - area_domain) > 1e-12*area_domain) then
                    write(6,*) '**** WARNING'
                    write(6,*) '**** topo arrays do not quite cover domain'
                    write(6,*) '**** area of overlap = ', area
                    write(6,*) '**** area of domain  = ', area_domain
                    write(6,*) '**** error is less than 1% so proceeding...'
                endif

            !---------------tests for analytic bathymetry-------------------
            ! Simple jump discontinuity in bathymetry
            else if (test_topography == 1) then
                topo_finalized = .true.
                read(iunit,"(d16.8)") topo_location
                read(iunit,"(d16.8)") topo_left
                read(iunit,"(d16.8)") topo_right

            ! Idealized ocean shelf
            else if (test_topography == 2 .or. test_topography == 3) then
                topo_finalized = .true.
                read(iunit,"(d16.8)") topo_x0
                read(iunit,"(d16.8)") topo_x1
                read(iunit,"(d16.8)") topo_x2
                read(iunit,"(d16.8)") topo_basin_depth
                read(iunit,"(d16.8)") topo_shelf_depth
                read(iunit,"(d16.8)") topo_beach_slope
                topo_shelf_slope = (topo_basin_depth - topo_shelf_depth) &
                                            / (topo_x0 - topo_x1)
            else
                print *,"Error:  Unknown test topography type ",test_topography
                stop
            endif

            module_setup = .true.
        end if

    end subroutine read_topo_settings

    ! ========================================================================
    !  set_topo_for_dtopo()
    !
    !  Set topography values in new topo arrays that correspond to dtopo spatialy
    !  array values come from the finest topography already in topowork
    ! ========================================================================

    subroutine set_topo_for_dtopo(mx,my,dx,dy,xlow,yhi,newtopo)

        !arguments
        integer, intent(in) :: mx,my
        real(kind=8), intent(in) :: dx,dy,xlow,yhi
        real(kind=8), intent(inout) :: newtopo(1:int(mx, 8) * int(my, 8))

        !locals
        integer :: id,irank,itopo1,itopo2,jtopo1,jtopo2
        integer(kind=8) :: i,j,ij,ijll,ijlr,ijul,ijur
        real(kind=8) :: x,y,xl,xr,yu,yl,zll,zlr,zul,zur,z,dxdy

        do j=1,int(my, 8)
               y = yhi - (j-1)*dy
            do i=1,int(mx, 8)
               x = xlow + (i-1)*dx
               ij = (j-1)*int(mx, 8) + i
               !find intersection starting from finest topo
               !all points must lie in some topo file therefore the
               !finest topo file for all dtopo points will be saved in topo0
               do irank = 1,mtopofiles
                  id = mtopoorder(irank)
                  if (id.gt.mtopofiles-num_dtopo) then
                     !this is another dtopo ==> topo file: skip
                     cycle
                  elseif ( (x>xhitopo(id)).or.(x<xlowtopo(id)).or. &
                          (y>yhitopo(id)).or.(y<ylowtopo(id))) then
                     !no intersection
                     cycle
                  else !lies in this topofile

                     ! Old way of setting topo0save up to v5.4.1.
                     ! This assumed topo file did not
                     ! intersect dtopo file if this point was never reached.
                     ! Not true if topofile has such small extent that it lies between
                     ! dtopo points.  Now instead we set topo0save earlier.
                     !topo0save(id) = 1

                     !find indices for bilinear cell in topo
                     !arrays are in form of DEM...high y values first
                     !note for xy points lying on nodes all indices will be equal
                     itopo1 = int(floor((x-xlowtopo(id))/dxtopo(id)))+1
                     itopo2 = int(ceiling((x-xlowtopo(id))/dxtopo(id)))+1
                     jtopo1 = int(floor((yhitopo(id)-y)/dytopo(id))) + 1
                     jtopo2 = int(ceiling((yhitopo(id)-y)/dytopo(id))) + 1
                     !indices for work array
                     ijll = i0topo(id) + (jtopo2-1)*mxtopo(id) + itopo1 -1
                     ijlr = i0topo(id) + (jtopo2-1)*mxtopo(id) + itopo2 -1
                     ijul = i0topo(id) + (jtopo1-1)*mxtopo(id) + itopo1 -1
                     ijur = i0topo(id) + (jtopo1-1)*mxtopo(id) + itopo2 -1
                     !find x,y,z values for bilinear
                     !z may be from only 1 or 2 nodes for aligned grids
                     !bilinear should still evaluate correctly
                     zll = topowork(ijll)
                     zlr = topowork(ijlr)
                     zul = topowork(ijul)
                     zur = topowork(ijur)
                     xl = xlowtopo(id) + real(itopo1-1,kind=8)*dxtopo(id)
                     xr = xl + dxtopo(id)
                     yu = yhitopo(id) - real(jtopo1-1,kind=8)*dytopo(id)
                     yl = yu - dytopo(id)
                     dxdy = dxtopo(id)*dytopo(id)
                     z = zll*(xr-x)*(yu-y) + zlr*(x-xl)*(yu-y) + zul*(xr-x)*(y-yl) + zur*(x-xl)*(y-yl)
                     newtopo(ij) = z/dxdy
                     !this was the finest topo file, move to next point
                     exit
                  endif
               enddo
            enddo
        enddo

    end subroutine set_topo_for_dtopo


    ! ========================================================================
    !  read_topo_file(mx,my,topo_type,fname,xll,yll,topo)
    !
    !  Read topo file.
    ! ========================================================================

    subroutine read_topo_file(mx,my,topo_type,fname,xll,yll,topo)

#ifdef NETCDF
        use netcdf
#endif

        use geoclaw_module
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,topo_type
        character(len=150), intent(in) :: fname
        real(kind=8), intent(inout) :: topo(1:int(mx, 8)*int(my, 8))
        real(kind=8), intent(in) :: xll,yll

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        logical, parameter :: maketype2 = .false.
        integer :: missing,status,n
        real(kind=8) :: no_data_value,x,y,topo_temp
        real(kind=8) :: values(16)
        character(len=80) :: str
        integer(kind=8) :: i, j, mtot

        ! NetCDF Support
        character(len=64) :: direction, x_dim_name, x_var_name, y_dim_name, &
            y_var_name, z_var_name, var_name
        ! character(len=1) :: axis_string
        real(kind=8), allocatable :: nc_buffer(:, :), xlocs(:), ylocs(:)
        integer(kind=4) :: x_var_id, y_var_id, z_var_id, x_dim_id, y_dim_id
        integer(kind=4) :: xstart(1), ystart(1), mx_tot, my_tot
        integer(kind=4) :: ios, nc_file, dim_ids(2), num_dims, &
            var_type, num_vars, num_dims_tot, z_dim_ids(2)

        mtot = int(mx, 8) * int(my, 8)

        print *, ' '
        print *, 'Reading topography file  ', fname

        select case(abs(topo_type))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                i = 0
                status = 0
                do i=1,mtot
                    read(iunit,fmt=*,iostat=status) x,y,topo_temp
                    if ((i > mtot) .and. (status == 0)) then
                        print *,'*** Error: i > mx*my = ',mtot
                        print *,'*** i, mx, my: ',i,mx,my
                        print *,'*** status = ',status
                        stop
                    endif

                    if (status /= 0) then
                        print *,"Error reading topography file, reached EOF."
                        print *,"  File = ",fname
                        stop
                    else
                        topo(i) = topo_temp
                    endif
                enddo
                
                close(unit=iunit)

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if topo_type=2 or
            ! mx values per line if topo_type=3
            ! ================================================================
            case(2:3)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo
                
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                no_data_value = values(1)

                ! Read in data
                missing = 0
                select case(abs(topo_type))
                    case(2)
                        do i=1,mtot
                            read(iunit,*) topo(i)
                            if (topo(i) == no_data_value) then
                                missing = missing + 1
                                topo(i) = topo_missing
                                ! uncomment next line to print row i
                                ! write(6,600) i
 600                            format('*** missing data, i = ',i6)
                            endif
                        enddo
                    case(3)
                        do j=1,int(my, 8)
                            read(iunit,*) (topo((j-1)*int(mx, 8) + i),i=1,mx)
                            do i=1,int(mx, 8)
                                if (topo((j-1)*int(mx, 8) + i) == no_data_value) then
                                    missing = missing + 1
                                    topo((j-1)*int(mx, 8) + i) = topo_missing
                                    ! uncomment next line to print row j
                                    ! and element i that are missing.
                                    ! write(6,601) i,j
 601                                format('*** missing data, i = ',i6, &
                                           '  j = ',i6)
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    write(6,602) missing
 602                format('WARNING... ',i6, &
                           ' missing data values in this topofile')
                    write(6,603) topo_missing
 603                format('   These values have been set to topo_missing = ',&
                           f13.3, ' in read_topo_file')
                    if (topo_missing == 99999.d0) then
                        print *, 'ERROR... do not use this default value'
                        print *, 'Fix your topofile or set'
                        print *, '  rundata.topo_data.topo_missing in setrun.py'
                        stop
                        endif
                    print *, ' '
                endif

                close(unit=iunit)
            
            ! NetCDF
            case(4)
#ifdef NETCDF
                ! Open file    
                call check_netcdf_error(nf90_open(fname, nf90_nowrite, nc_file))
                
                ! Get number of dimensions and vars
                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, &
                    num_vars))
                
                ! get x and y dimension info
                call get_dim_info(nc_file, num_dims_tot, x_dim_id, x_dim_name, &
                    mx_tot, y_dim_id, y_dim_name, my_tot)
                    
                allocate(xlocs(mx_tot),ylocs(my_tot))
                
                call check_netcdf_error(nf90_get_var(nc_file, x_dim_id, xlocs, start=(/ 1 /), count=(/ mx_tot /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_dim_id, ylocs, start=(/ 1 /), count=(/ my_tot /)))
                xstart = minloc(xlocs, mask=(xlocs.eq.xll))
                ystart = minloc(ylocs, mask=(ylocs.eq.yll))
                deallocate(xlocs,ylocs)
                
                z_var_id = -1
                do n=1, num_vars
                    call check_netcdf_error(nf90_inquire_variable(nc_file, n, var_name, var_type, num_dims, dim_ids))
                    
                    ! Assume dim names are same as var ids
                    if (var_name == x_dim_name) then
                        x_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, x_var_name, x_var_id))
                    else if (var_name == y_dim_name) then
                        y_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, y_var_name, y_var_id))
                    else if (num_dims == 2) then
                        z_var_name = var_name
                        z_dim_ids = dim_ids
                        call check_netcdf_error(nf90_inq_varid(nc_file, z_var_name, z_var_id))
                    end if

                end do
                if (z_var_id == -1) then
                    stop "Unable to find topography data!"
                end if

                ! Load in data
                ! TODO: Provide striding into data if need be
                
                ! only try to access data if theres overlap with the domain
                if ((mx > 0) .and. (my > 0)) then
                    if (z_dim_ids(1) == x_dim_id) then
                        allocate(nc_buffer(mx, my))
                    else if (z_dim_ids(1) == y_dim_id) then
                        allocate(nc_buffer(my, mx))
                    else 
                        stop " NetCDF z variable has dimensions that don't align with x and y"
                    end if

                    call check_netcdf_error(nf90_get_var(nc_file, z_var_id, nc_buffer, &
                                                start = (/ xstart(1), ystart(1) /), &
                                                count = (/ mx, my /)))

                    ! check for lon/lat ordering of z variable
                    if (z_dim_ids(1) == x_dim_id) then
                        do j = 0, int(my, 8) - 1
                            topo(j * int(mx, 8) + 1:(j + 1) * int(mx, 8)) = nc_buffer(:, my - j)
                        end do
                    else if (z_dim_ids(1) == y_dim_id) then
                        do j = 0, int(my, 8) - 1
                            topo(j * int(mx, 8) + 1:(j + 1) * int(mx, 8)) = nc_buffer(my - j, :)
                        end do
                    end if
                    deallocate(nc_buffer)

                    ! do j=1, my
                    !     i = i - 1
                    !     topo((i - 1) * mx + 1:)
                    ! print *, mx, my
                    ! do j=1, my
                    !     row_index = row_index - 1
                    !     print *, row_index
                    !     print *, (row_index-1)*mx + 1, (row_index-1)*mx + mx
                    !     ! call check_netcdf_error(nf90_get_var(nc_file, z_var_id,    &
                    !             ! topo))
                    !     ! call check_netcdf_error(nf90_get_var(nc_file, z_var_id,  &
                    !     !       topo((row_index-1)*mx + 1:(row_index-1)*mx + mx),  &
                    !     !                 start = (/ 1, j /), count=(/ mx, 1 /)))
                    !     call check_netcdf_error(nf90_get_var(nc_file, z_var_id,  &
                    !           nc_buffer, start = (/ 1, j /), count=(/ mx, 1 /)))
                    !     topo((row_index - 1) * mx + 1:(row_index - 1) * mx + mx) = &
                    !             nc_buffer
                    ! end do

                    ! Check if the topography was defined positive down and flip the
                    ! sign if need be.  Just in case this is true but topo_type < 0
                    ! we do not do anything here on this to avoid doing it twice.
                    ios = nf90_get_att(nc_file, z_var_id, 'positive', direction)
                    call check_netcdf_error(nf90_close(nc_file))
                    if (ios == NF90_NOERR) then
                        if (to_lower(direction) == "down") then
                            if (topo_type < 0) then
                                topo = -topo
                            endif
                        end if
                    end if
                end if
#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of GeoClaw."
                stop
#endif
        end select

        ! Handle negative topo types
        if (topo_type < 0) then
            forall(i=1:mtot)
                topo(i) = -topo(i)
            end forall
        endif

        ! ====================================================================
        ! when topo_type=1 and data has x,y,z columns,
        ! set maketype2 to true to create a file new.topo_type2 with only z
        ! values, so next time it will take less time to read in.
        ! only works if dx = dy.
        ! Currently we don't have the right info in this routine to do this
!         if ((topo_type == 1).and.maketype2) then
!             open(unit=29,file='new.tt2',status='unknown',form='formatted')
!             write(29,*) mx, '       mx'
!             write(29,*) my, '       my'
!             write(29,*) xll, '     xllcorner'
!             write(29,*) yll, '     yllcorner'
!             write(29,*) dx, '     cellsize'
!             write(29,*) -9999, '     nodataval'
!             do i=topo_start,itopo
!                 write(29,'(d18.8)') topowork(i)
!             enddo
!             close(unit=29)
!             print *, 'New topo file new.tt2 created'
!             if (dx.ne.dy) then
!                 print *,  ' ** Warning, dx must equal dy for this'
!                 print *,  ' dx = ',dx,'  dy = ',dy
!             endif
!         endif
        ! ====================================================================

    end subroutine read_topo_file

    ! ========================================================================
    ! subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read topo file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - topo_type - (int) Type of topography file (-3 < topo_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy)
#ifdef NETCDF
        use netcdf
#endif

        use geoclaw_module
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Input and Output
        character(len=150), intent(in) :: fname
        integer, intent(in) :: topo_type
        integer, intent(out) :: mx, my
        real(kind=8), intent(out) :: xll, yll, xhi, yhi, dx, dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: topo_size, status, n
        real(kind=8) :: x,y,z,nodata_value
        logical :: found_file
        real(kind=8) :: values(16)
        character(len=80) :: str
        logical :: verbose
        logical :: xll_registered, yll_registered

        ! NetCDF Support
        ! character(len=1) :: axis_string
        ! character(len=6) :: convention_string
        ! integer(kind=4) :: convention_version
        integer(kind=4) :: nc_file
        real(kind=8), allocatable :: xlocs(:),ylocs(:)
        logical, allocatable :: x_in_dom(:),y_in_dom(:)
        integer(kind=4) :: dim_ids(2), num_dims, var_type, num_vars, num_dims_tot
        integer(kind=4), allocatable :: var_ids(:)
        character(len=64) :: var_name, x_var_name, y_var_name, z_var_name
        character(len=64) :: x_dim_name, y_dim_name
        integer(kind=4) :: x_var_id, y_var_id, z_var_id, x_dim_id, y_dim_id

        verbose = .false.

        inquire(file=fname, exist=found_file)
        if (.not. found_file) then
            print *, 'Missing topography file:'
            print *, '   ', fname
            stop
        endif

        select case(abs(topo_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                open(unit=iunit, file=fname, status='unknown',form='formatted')

                ! Initial size variables
                topo_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                topo_size = topo_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    topo_size = topo_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                do
                    read(iunit,fmt=*,iostat=status) x,y,z
                    if (status /= 0) exit
                    topo_size = topo_size + 1
                enddo
                if (status > 0) then
                    print *,"ERROR:  Error reading header of topography file ",fname
                    stop
                endif

                ! Calculate remaining values
                my = topo_size / mx
                xhi = x
                yll = y
                dx = (xhi-xll) / (mx-1)
                dy = (yhi-yll) / (my-1)

            ! ASCII file with header followed by z data
            case(2:3)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                mx = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                my = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                xll = values(1)
                str = to_lower(str)  ! convert to lower case
                xll_registered = (index(str, 'xllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                yll = values(1)
                str = to_lower(str)  ! convert to lower case
                yll_registered = (index(str, 'yllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                dx = values(1)
                if (n == 2) then
                    dy = values(2)
                  else
                    dy = dx
                  endif

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                nodata_value = values(1)

                if (xll_registered) then
                    xll = xll + 0.5d0*dx
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting xllcorner by 0.5*dx to cell center'
                    endif 

                if (yll_registered) then
                    yll = yll + 0.5d0*dy
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting yllcorner by 0.5*dy to cell center'
                    endif 


                xhi = xll + (mx-1)*dx
                yhi = yll + (my-1)*dy
                
            ! NetCDF
            case(4)
#ifdef NETCDF

                ! Open file    
                call check_netcdf_error(nf90_open(fname, nf90_nowrite, nc_file))

                ! NetCDF4 GEBCO topography, should conform to CF metadata
                ! standard
                ! ios = nf90_get_att(nc_file, NF90_GLOBAL, 'Conventions', convention_string)
                ! if (verbose) then
                    ! print *, convention_string
                ! end if
                ! if (ios /= NF90_NOERR .or. convention_string(1:3) /= "CF-") then
                !     print *, "Topography file does not conform to the CF"
                !     print *, "conventions and meta-data.  Please see the"
                !     print *, "information at "
                !     print *, ""
                !     print *, "    cfconventions.org"
                !     print *, ""
                !     print *, "to find out more info."
                !     stop
                ! end if

                ! call parse_values(convention_string(4:6), num_values, convention_version)
                ! if (convention_version(1) < CONVENTION_REQUIREMENT) then
                !     print *, "Topography file conforms to a CF version that"
                !     print *, "is too old (", convention_version, ").  "
                !     print *, "Please refer to"
                !     print *, ""
                !     print *, "    cfconventions.org"
                !     print *, ""
                !     print *, "to find out more info."
                !     stop
                ! end if

                ! get x and y dimension info
                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, &
                    num_vars))
                call get_dim_info(nc_file, num_dims_tot, x_dim_id, x_dim_name, &
                    mx, y_dim_id, y_dim_name, my)
                
                ! allocate vector to hold lon and lat vals and vector for var ids
                allocate(xlocs(mx),ylocs(my),x_in_dom(mx),y_in_dom(my),var_ids(num_vars))

                if (verbose) then
                    print *, "Names = (", x_dim_name, ", ", y_dim_name, ")"
                    print *, "Size = (", mx, ", ", my, ")"
                end if

                ! Read in variables
                call check_netcdf_error(nf90_inq_varids(nc_file, num_vars, var_ids))
                if (verbose) then
                    print *, "n, var_name, var_type, num_dims, dim_ids"
                end if

                do n=1, num_vars
                    call check_netcdf_error(nf90_inquire_variable(nc_file, n, var_name, var_type, num_dims, dim_ids))
                    if (verbose) then
                        print *, n, ": ", var_name, "|", var_type, "|", num_dims, "|", dim_ids
                    end if

                    ! Assume dim names are same as var ids
                    if (var_name == x_dim_name) then
                        x_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, x_var_name, x_var_id))
                        ! x_var_id = n
                        ! x_var_name = var_name
                    else if (var_name == y_dim_name) then
                        y_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, y_var_name, y_var_id))
                        ! y_var_id = n
                        ! y_var_name = var_name
                    ! Assume that this is the topography data
                    else if (num_dims == 2) then
                        z_var_name = var_name
                        call check_netcdf_error(nf90_inq_varid(nc_file, z_var_name, z_var_id))
                        ! z_var_id = n
                        ! z_var_name = var_name
                    else
                        if (verbose) then
                            print *, "Not using var_id ", n
                        end if
                    end if

                end do

                if (verbose) then
                    print *, "x_var_name, x_var_id = ", x_var_name, x_var_id
                    print *, "y_var_name, y_var_id = ", y_var_name, y_var_id
                    print *, "z_var_name, z_var_id = ", z_var_name, z_var_id
                end if

                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xlocs, start=(/ 1 /), count=(/ mx /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, ylocs, start=(/ 1 /), count=(/ my /)))
                
                dx = xlocs(2) - xlocs(1)
                dy = ylocs(2) - ylocs(1)
                
                ! find which locs are within domain (with a dx/dy buffer around domain + ghost cells)
                x_in_dom = (xlocs.gt.(xlower-dx-hxposs(1)*nghost)) .and. (xlocs.lt.(xupper+dx+hxposs(1)*nghost))
                y_in_dom = (ylocs.gt.(ylower-dy-hyposs(1)*nghost)) .and. (ylocs.lt.(yupper+dy+hyposs(1)*nghost))
                
                xll = minval(xlocs, mask=x_in_dom)
                yll = minval(ylocs, mask=y_in_dom)
                xhi = maxval(xlocs, mask=x_in_dom)
                yhi = maxval(ylocs, mask=y_in_dom)
                
                ! adjust mx, my
                mx = count(x_in_dom)
                my = count(y_in_dom)

                call check_netcdf_error(nf90_close(nc_file))
                deallocate(xlocs,ylocs,x_in_dom,y_in_dom)
#else
                print *, "ERROR:  NetCDF library was not included in this build"
                print *, "  of GeoClaw."
                stop
#endif

            case default
                print *, 'ERROR:  Unrecognized topo_type'
                print *, '    topo_type = ',topo_type
                print *, '  for topography file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)
        write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_topo_header

    real(kind=8) pure function test_topo(x) result(topography)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x

        if (test_topography == 1) then
            if (x < topo_location) then
                topography = topo_left
            else
                topography = topo_right
            endif
        else if (test_topography == 2) then
            if (x < topo_x0) then
                topography = topo_basin_depth
            else if (topo_x0 <= x .and. x < topo_x1) then
                topography = topo_shelf_slope * (x-topo_x0) + topo_basin_depth
            else if (topo_x1 <= x .and. x < topo_x2) then
                topography = topo_shelf_depth
            else
                topography = topo_beach_slope * (x-topo_x2) + topo_shelf_depth
            endif
        endif

    end function test_topo


    ! ========================================================================
    !  set_dtopo(fname)
    ! ========================================================================
    ! Read moving topography info from setdtopo.data
    ! Time-dependend topography is used to initiate tsunami, for example.
    !
    ! If num_dtopo = 0, no movement of topography specified.
    !
    ! If num_dtopo = 1, the topo changes dynamically.
    ! the topofile can then be formated in the following ways
    ! dtopotype > 1: file contains a header, with a line for each record
    ! my; mx; mt; xlower; ylower; t0; dx; dy; dt;
    ! dtopotype = 1:
    ! Then the dtopofile should have 4 columns:
    ! time,longitude,latitude,vertical displacement(m) since t0
    !
    ! Longitude and latitude advance in the standard GIS way from
    ! upper left corner across in x and then down in y.
    ! Time column advances most slowly.
    ! ========================================================================
    subroutine read_dtopo_settings(file_name)

        use geoclaw_module

        implicit none

        ! Input arguments
        character (len=25), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 79
        integer :: finer_than,rank
        real(kind=8) :: area_i,area_j
        integer :: i,j

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SETDTOPO:'
        write(GEO_PARM_UNIT,*) '-------------'

        if (present(file_name)) then
            call opendatafile(iunit,file_name)
        else
            call opendatafile(iunit,'dtopo.data')
        endif

        read(iunit,*) num_dtopo
        write(GEO_PARM_UNIT,*) '   num dtopo files = ',num_dtopo
        if (num_dtopo == 0) then
            return
        endif

        ! Allocate and read in dtopo info
        allocate(dtopofname(num_dtopo))
        allocate(mxdtopo(num_dtopo))
        allocate(mydtopo(num_dtopo),mtdtopo(num_dtopo),mdtopo(num_dtopo))
        allocate(xlowdtopo(num_dtopo),ylowdtopo(num_dtopo),t0dtopo(num_dtopo))
        allocate(xhidtopo(num_dtopo),yhidtopo(num_dtopo),tfdtopo(num_dtopo))
        allocate(dtopotype(num_dtopo),i0dtopo(num_dtopo))
        allocate(dxdtopo(num_dtopo),dydtopo(num_dtopo),dtdtopo(num_dtopo))
        allocate(kdtopo1(num_dtopo),kdtopo2(num_dtopo))
        allocate(index0_dtopowork1(num_dtopo),index0_dtopowork2(num_dtopo))
        allocate(tdtopo1(num_dtopo),tdtopo2(num_dtopo),taudtopo(num_dtopo))
        allocate(mdtopoorder(num_dtopo),topoaltered(num_dtopo))

        do i=1,num_dtopo
            read(iunit,*) dtopofname(i)
            read(iunit,*) dtopotype(i)
            write(GEO_PARM_UNIT,*) '   fname:',dtopofname(i)
            write(GEO_PARM_UNIT,*) '   topo type:',dtopotype(i)

            ! Read in header data
            call read_dtopo_header(dtopofname(i),dtopotype(i),mxdtopo(i), &
                mydtopo(i),mtdtopo(i),xlowdtopo(i),ylowdtopo(i),t0dtopo(i), &
                xhidtopo(i),yhidtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i), &
                dtdtopo(i))
            mdtopo(i) = mxdtopo(i) * mydtopo(i) * mtdtopo(i)
        enddo

        read(iunit,*) dt_max_dtopo  
        !  largest allowable dt while dtopo is moving


        ! Indexing into work array
        i0dtopo(1) = 1
        if (num_dtopo > 1) then
            do i=2,num_dtopo
                i0dtopo(i) = i0dtopo(i-1) + mdtopo(i-1)
            enddo
        endif

            ! dtopo order..for updating topo from finest dtopo model
            ! The finest topography will be given priority in any region
            ! mtopoorder(rank) = i means that i'th topography file has rank rank,
            ! where the file with rank=1 is the finest and considered first.
            do i=1,num_dtopo
                finer_than = 0
                do j=1,num_dtopo
                    if (j /= i) then
                        area_i=dxdtopo(i)*dydtopo(i)
                        area_j=dxdtopo(j)*dydtopo(j)
                        if (area_i < area_j) finer_than = finer_than + 1
                        ! if two files have the same resolution, order is
                        ! arbitrarily chosen
                        if ((area_i == area_j).and.(j < i)) then
                            finer_than = finer_than + 1
                        endif
                    endif
                enddo
                ! ifinerthan tells how many other files, file i is finer than
                rank = num_dtopo - finer_than
                mdtopoorder(rank) = i
            enddo

        ! Allocate and read dtopo files
        allocate(dtopowork(sum(mdtopo)))

        do i=1,num_dtopo
            call read_dtopo(mxdtopo(i),mydtopo(i),mtdtopo(i),dtopotype(i), &
                dtopofname(i),dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1))
        enddo

    end subroutine read_dtopo_settings
    ! ========================================================================

    ! ========================================================================
    !  read_dtopo(fname)
    ! ========================================================================
    subroutine read_dtopo(mx,my,mt,dtopo_type,fname,dtopo)

      implicit none

      ! Arguments
      integer, intent(in) :: mx,my,mt,dtopo_type
      character (len=150), intent(in) :: fname
      real(kind=8), intent(inout) :: dtopo(1:int(mx, 8)*int(my, 8)*int(mt, 8))

      ! Local
      integer, parameter :: iunit = 29
      integer :: status
      real(kind=8) :: t,x,y
      integer(kind=8) :: i, j, k, mtot

      mtot = int(mx, 8) * int(my, 8)

      open(unit=iunit, file=fname, status = 'unknown',form='formatted')

      select case(abs(dtopo_type))
         case(1)
            ! ASCII file with 4 columns
            do i = 1,mtot*mt
               read(iunit,fmt=*,iostat=status) t,x,y, dtopo(i)
            enddo

         case(2)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            ! read the data
            do i = 1,mtot*mt
               read(iunit,*) dtopo(i)
            enddo
         case(3)
            ! read header
            do i = 1,9
               read(iunit,*)
            enddo
            do k = 1,int(mt, 8)
               do j = 1,int(my, 8)
                  read(iunit,*) (dtopo((k-1)*mtot + (j-1)*int(mx, 8) + i) , i=1,int(mx, 8))
               enddo
            enddo
      end select

    end subroutine read_dtopo

    ! ========================================================================
    !  subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi,
    !                               yhi,tf,dx,dy,dt)
    ! ========================================================================

    !  Read in dtopo file header and either read or calculate the grid info
    !
    !  :Input:
    !   - fname - (char) Name of the dtopo file
    !   - topo_type - (int) Topography file type (1-3 are valid)
    !
    !  :Output:
    !   - mx,my,mt - (int) Number of grid point in space (mx,my) and time (mt)
    !   - xlow,ylow - (dp) Lower corner spatial coordinate of grid
    !   - xhi,yhi - (dp) Upper corner spatial coodinate of grid
    !   - t0,tf - (dp) Beginning and end times for the file
    !   - dx,dy,dt - (dp) Distance between space (dx,dy) and time (dt) points
    ! ========================================================================
    subroutine read_dtopo_header(fname,topo_type,mx,my,mt,xlow,ylow,t0,xhi, &
        yhi,tf,dx,dy,dt)

        implicit none

        ! Input Arguments
        character (len=150), intent(in) :: fname
        integer, intent(in) :: topo_type

        ! Output Arguments
        integer, intent(out) :: mx,my,mt
        real(kind=8), intent(out) :: xlow,ylow,t0,xhi,yhi,tf,dx,dy,dt

        ! Locals
        integer, parameter :: iunit = 7
        integer :: topo_size,status
        real(kind=8) :: x,y,t,y_old
        logical :: found_file

        ! Open file
        inquire(file=fname,exist=found_file)
        if (.not.found_file) then
            print *, 'Missing dtopo file:'
            print *, '    ', fname
            stop
        endif
        open(unit=iunit,file=fname,status='unknown',form='formatted')

        select case(topo_type)
            ! Old style ASCII dtopo files
            case(1)
                ! Initial size variables
                topo_size = 0
                my = 1
                mt = 1

                ! Read in first values, determines xlow, yhi and t0
                read(iunit,*) t0,xlow,yhi
                topo_size = topo_size + 1
                t = t0
                y_old = yhi
                ! Go through entire file figuring out my, mt and topo_size
                status = 0
                do while (status == 0.and. t.eq.t0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                    if (y /= y_old .and. t.eq.t0 ) then
                        my = my + 1
                        y_old = y
                    endif
                enddo
                mx = (topo_size-1)/my
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) t,x,y
                    topo_size = topo_size + 1
                enddo


                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                mt = (topo_size-1)/ (my*mx)
                xhi = x
                ylow = y
                tf = t
                dx = (xhi-xlow) / (mx-1)
                dy = (yhi-ylow) / (my-1)
                dt = (tf - t0) / (mt-1)

            ! New ASCII headered dtopo files, similar to topography files type
            ! 2 and 3
            case(2:3)
                ! Read in header directly
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) mt
                read(iunit,*) xlow
                read(iunit,*) ylow
                read(iunit,*) t0
                read(iunit,*) dx
                read(iunit,*) dy
                read(iunit,*) dt

                xhi = xlow + dx*(mx-1)
                yhi = ylow + dy*(my-1)
                tf = t0 + dt*(mt-1)
            case default
                print *, 'ERROR:  Unrecognized topography type'
                print *, '    topo_type = ',topo_type
                print *, '  for dtopo file:'
                print *, '   ', fname
                stop
        end select

         close(iunit)
    end subroutine read_dtopo_header


    

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
    integer :: mfid, indicator
    integer(kind=8) :: i0
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
    integer :: mfid, indicator
    integer(kind=8) :: i0
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

#ifdef NETCDF
    subroutine check_netcdf_error(ios)

        use netcdf

        implicit none

        integer, intent(in) :: ios

        if (ios /= NF90_NOERR) then
            print *, "NetCDF IO error: ", ios
            print *, trim(nf90_strerror(ios))
            stop
        end if

    end subroutine check_netcdf_error

    subroutine get_dim_info(nc_file, ndims, x_dim_id, x_dim_name, mx, &
        y_dim_id, y_dim_name, my)
        use netcdf
        implicit none
        integer, intent(in) :: nc_file, ndims
        integer, intent(out) :: x_dim_id, y_dim_id, mx, my
        character (len = *), intent(out) :: x_dim_name, y_dim_name
        integer :: m_tmp, n
        character(20) :: dim_name_tmp

        ! get indices to start at for reading netcdf within domain
        do n=1, ndims
            call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                n, dim_name_tmp, m_tmp))
            if (ANY((/ 'LON      ','LONGITUDE','X        ' /) == Upper(dim_name_tmp))) then
                x_dim_name = dim_name_tmp
                mx = m_tmp
                x_dim_id = n
            else if (ANY((/ 'LAT     ','LATITUDE','Y       ' /) == Upper(dim_name_tmp))) then
                y_dim_name = dim_name_tmp
                my = m_tmp
                y_dim_id = n
            end if
        end do
    end subroutine get_dim_info
#endif

function Upper(s1)  RESULT (s2)
    CHARACTER(*)       :: s1
    CHARACTER(LEN(s1)) :: s2
    CHARACTER          :: ch
    INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
    INTEGER            :: i

    DO i = 1,LEN(s1)
       ch = s1(i:i)
       IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
       s2(i:i) = ch
    END DO
END function Upper

end module topo_module
