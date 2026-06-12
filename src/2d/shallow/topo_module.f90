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
    character(len=512), allocatable :: topofname(:)
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

    ! NetCDF descriptor fields — one slot per topo file, populated from the
    ! key=value block in topo.data by read_netcdf_descriptor (type 4 only).
    character(len=64), allocatable :: nc_var_name(:)     ! topo variable name
    character(len=64), allocatable :: nc_lon_name(:)     ! longitude coord name
    character(len=64), allocatable :: nc_lat_name(:)     ! latitude coord name
    character(len=32), allocatable :: nc_lat_order(:)    ! 'N_to_S' or 'S_to_N'
    character(len=32), allocatable :: nc_dim_order(:)    ! e.g. 'lat,lon'
    character(len=8),  allocatable :: nc_fill_action(:)  ! 'abort' or 'warn'
    real(kind=8),      allocatable :: nc_lon_offset(:)     ! scalar added to file coords: x_domain = x_file + offset
    real(kind=8),      allocatable :: nc_fill_value(:)   ! fill/nodata sentinel
    real(kind=8),      allocatable :: nc_crop_bounds(:,:) ! (4, n): lon0,lon1,lat0,lat1
    logical,           allocatable :: nc_has_fill(:)     ! fill_value was specified
    logical,           allocatable :: nc_has_crop(:)     ! crop_bounds was specified

    ! Per-file preprocessing parameters — read from the 7 new lines in topo.data
    ! for all topo file types (2, 3, 4, 5). Type-1 files are skipped.
    real(kind=8), allocatable :: tp_crop_extent(:,:)  ! (4,n): x1,x2,y1,y2 domain coords; all-zero=no crop
    integer,      allocatable :: tp_coarsen(:)         ! subsampling factor (1=none)
    integer,      allocatable :: tp_buffer(:)          ! buffer cells (grid points) around crop region
    real(kind=8), allocatable :: tp_align(:,:)         ! (2,n): align offset; all-zero=none
    real(kind=8), allocatable :: tp_x_shift(:)         ! coordinate registration offset
    real(kind=8), allocatable :: tp_z_shift(:)         ! bulk datum offset
    logical,      allocatable :: tp_negate_z(:)        ! explicit Z sign flip

    ! Full (un-cropped) file sizes and kept-window start indices for types 2/3,
    ! saved by read_topo_settings so read_topo_file can extract the kept
    ! (cropped/coarsened) lattice: 0-based indices from the file's south-west
    ! corner; kept points are tp_i_start + k*tp_coarsen, k = 0..mxtopo-1.
    integer,      allocatable :: mxtopo_full(:), mytopo_full(:)
    integer,      allocatable :: tp_i_start(:), tp_j_start(:)

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
        integer :: i,j,rank,irank,jrank,krank,idtopo
        real(kind=8) :: area, area_domain, area_maxj, area_tol
        real(kind=8), allocatable :: areatopo(:)
        ! Crop index locals (types 2/3 only)
        integer :: i_start, i_end, j_start, j_end
        real(kind=8) :: x1_c, x2_c, y1_c, y2_c
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

                ! modified topo ordering algorithm for v5.14.0:
                ! don't increment mtopofiles by num_dtopo until
                ! after the initial ordering is done

                write(GEO_PARM_UNIT,*) '   mtopofiles = ',mtopofiles

                ! Allocate arrays for each file,
                ! including space for the topo_for_dtopo array:
                mtopofiles = mtopofiles + num_dtopo  ! increment for allocates
                allocate(mxtopo(mtopofiles),mytopo(mtopofiles))
                allocate(xlowtopo(mtopofiles),ylowtopo(mtopofiles))
                allocate(tlowtopo(mtopofiles),xhitopo(mtopofiles),yhitopo(mtopofiles))
                allocate(thitopo(mtopofiles),dxtopo(mtopofiles),dytopo(mtopofiles))
                allocate(topofname(mtopofiles),itopotype(mtopofiles))
                allocate(i0topo(mtopofiles),mtopo(mtopofiles),mtopoorder(mtopofiles))
                allocate(topoID(mtopofiles),topotime(mtopofiles),topo0save(mtopofiles))
                allocate(i0topo0(mtopofiles),topo0ID(mtopofiles))
                allocate(areatopo(mtopofiles))

                ! NetCDF descriptor arrays — allocated for all files; only
                ! filled in for type 4 entries by read_netcdf_descriptor.
                allocate(nc_var_name(mtopofiles), nc_lon_name(mtopofiles))
                allocate(nc_lat_name(mtopofiles), nc_lat_order(mtopofiles))
                allocate(nc_dim_order(mtopofiles), nc_fill_action(mtopofiles))
                allocate(nc_lon_offset(mtopofiles), nc_fill_value(mtopofiles))
                allocate(nc_crop_bounds(4, mtopofiles))
                allocate(nc_has_fill(mtopofiles), nc_has_crop(mtopofiles))
                nc_var_name      = ''
                nc_lon_name      = ''
                nc_lat_name      = ''
                nc_lat_order     = 'S_to_N'
                nc_dim_order     = 'lat,lon'
                nc_fill_action   = 'abort'
                nc_lon_offset    = 0.0d0
                nc_fill_value    = 0.0d0
                nc_crop_bounds   = 0.0d0
                nc_has_fill      = .false.
                nc_has_crop      = .false.

                allocate(tp_crop_extent(4, mtopofiles), tp_coarsen(mtopofiles))
                allocate(tp_buffer(mtopofiles), tp_align(2, mtopofiles))
                allocate(tp_x_shift(mtopofiles), tp_z_shift(mtopofiles))
                allocate(tp_negate_z(mtopofiles))
                tp_crop_extent = 0.0d0
                tp_coarsen     = 1
                tp_buffer      = 0
                tp_align       = 0.0d0
                tp_x_shift     = 0.0d0
                tp_z_shift     = 0.0d0
                tp_negate_z    = .false.

                allocate(mxtopo_full(mtopofiles), mytopo_full(mtopofiles))
                allocate(tp_i_start(mtopofiles), tp_j_start(mtopofiles))
                mxtopo_full = 0
                mytopo_full = 0
                tp_i_start  = 0
                tp_j_start  = 0

                mtopofiles = mtopofiles - num_dtopo  ! decrement after allocates

                do i=1,mtopofiles
                    read(iunit,*) topofname(i)
                    read(iunit,*) itopotype(i)

                    ! Deprecation warning for topo_type=1 (unstructured).
                    if (abs(itopotype(i)) == 1) then
                        print *, "WARNING: topo_type=1 is deprecated. Convert to topo_type=2, 3,"
                        print *, "  or 4 using topotools.Topography. Unstructured type-1 data is"
                        print *, "  not supported and will produce incorrect results."
                    end if

                    ! Read the 7 new preprocessing parameter lines present in
                    ! topo.data for ALL file types (written by TopographyData.write()).
                    read(iunit,*) tp_crop_extent(1,i), tp_crop_extent(2,i), &
                                  tp_crop_extent(3,i), tp_crop_extent(4,i)
                    read(iunit,*) tp_coarsen(i)
                    read(iunit,*) tp_buffer(i)
                    read(iunit,*) tp_align(1,i), tp_align(2,i)
                    read(iunit,*) tp_x_shift(i)
                    read(iunit,*) tp_z_shift(i)
                    read(iunit,*) tp_negate_z(i)

                    ! Preprocessing is not supported for topo_type=1.
                    if (abs(itopotype(i)) == 1 .and. ( &
                            any(tp_crop_extent(:,i) /= 0.0d0) .or. &
                            tp_coarsen(i) > 1 .or. &
                            tp_buffer(i) /= 0 .or. &
                            any(tp_align(:,i) /= 0.0d0) .or. &
                            tp_negate_z(i) .or. &
                            tp_x_shift(i) /= 0.0d0 .or. &
                            tp_z_shift(i) /= 0.0d0)) then
                        print *, "ERROR: preprocessing attributes (crop, coarsen, shift) are not"
                        print *, "  supported for topo_type=1. Convert to type 2/3/4 first."
                        stop 1
                    end if

                    ! For NetCDF files, read the optional key=value descriptor
                    ! block that follows the preprocessing lines.
                    if (abs(itopotype(i)) == 4) then
                        call read_netcdf_descriptor(iunit, i)
                    end if

                    write(GEO_PARM_UNIT,*) '   '
                    write(GEO_PARM_UNIT,*) '   ',topofname(i)
                    write(GEO_PARM_UNIT,*) '  itopotype = ', itopotype(i)
                    call read_topo_header(topofname(i),itopotype(i),mxtopo(i), &
                        mytopo(i),xlowtopo(i),ylowtopo(i),xhitopo(i),yhitopo(i), &
                        dxtopo(i),dytopo(i), i)

                    ! For types 2/3: apply crop_extent/align/buffer/coarsen by
                    ! computing the kept index window, mirroring Python
                    ! Topography.crop().  The full file sizes and the window
                    ! start are saved so read_topo_file can extract the kept
                    ! lattice.  For type 4 the window is computed inside
                    ! read_topo_header and applied via strided NetCDF reads.
                    if (abs(itopotype(i)) == 2 .or. abs(itopotype(i)) == 3) then
                        mxtopo_full(i) = mxtopo(i)
                        mytopo_full(i) = mytopo(i)

                        if (any(tp_crop_extent(:,i) /= 0.0d0) .or. &
                                tp_coarsen(i) > 1) then

                            ! Base window: indices of file points inside the
                            ! crop extent (first >= lower edge, last <= upper
                            ! edge, with a tolerance so grid-aligned extents
                            ! stay exact), or the full file if only coarsening.
                            ! tp_crop_extent is in domain coordinates; subtract
                            ! tp_x_shift to get file coordinates (the header
                            ! values are in file coords at this point, before
                            ! tp_x_shift is applied below).
                            if (any(tp_crop_extent(:,i) /= 0.0d0)) then
                                x1_c = tp_crop_extent(1,i) - tp_x_shift(i)
                                x2_c = tp_crop_extent(2,i) - tp_x_shift(i)
                                y1_c = tp_crop_extent(3,i)
                                y2_c = tp_crop_extent(4,i)
                                i_start = max(0, ceiling((x1_c - xlowtopo(i)) &
                                                / dxtopo(i) - 1.0d-6))
                                i_end   = min(mxtopo(i)-1, &
                                              floor((x2_c - xlowtopo(i)) &
                                                / dxtopo(i) + 1.0d-6))
                                j_start = max(0, ceiling((y1_c - ylowtopo(i)) &
                                                / dytopo(i) - 1.0d-6))
                                j_end   = min(mytopo(i)-1, &
                                              floor((y2_c - ylowtopo(i)) &
                                                / dytopo(i) + 1.0d-6))
                                if (i_start > i_end .or. j_start > j_end) then
                                    print *, "ERROR: tp_crop_extent does not overlap topo file:"
                                    print *, "  file = ", trim(topofname(i))
                                    print *, "  crop = ", tp_crop_extent(:,i)
                                    stop 1
                                end if
                            else
                                i_start = 0
                                i_end   = mxtopo(i) - 1
                                j_start = 0
                                j_end   = mytopo(i) - 1
                            end if

                            ! Align targets are in domain coordinates; convert
                            ! to file coordinates for x.  All-zero align is the
                            ! "no alignment" sentinel (matches the topo.data
                            ! format written by TopographyData.write()).
                            call apply_align_buffer_coarsen(i_start, i_end, &
                                mxtopo(i), xlowtopo(i), dxtopo(i), &
                                tp_align(1,i) - tp_x_shift(i), &
                                any(tp_align(:,i) /= 0.0d0), &
                                max(1, tp_coarsen(i)), tp_buffer(i))
                            call apply_align_buffer_coarsen(j_start, j_end, &
                                mytopo(i), ylowtopo(i), dytopo(i), &
                                tp_align(2,i), &
                                any(tp_align(:,i) /= 0.0d0), &
                                max(1, tp_coarsen(i)), tp_buffer(i))

                            ! Update stored extents/sizes to the kept lattice
                            ! (still file coords; tp_x_shift is applied below).
                            xlowtopo(i) = xlowtopo(i) + real(i_start,8)*dxtopo(i)
                            xhitopo(i)  = xlowtopo(i) + real(i_end-i_start,8)*dxtopo(i)
                            ylowtopo(i) = ylowtopo(i) + real(j_start,8)*dytopo(i)
                            yhitopo(i)  = ylowtopo(i) + real(j_end-j_start,8)*dytopo(i)
                            mxtopo(i) = (i_end - i_start)/max(1, tp_coarsen(i)) + 1
                            mytopo(i) = (j_end - j_start)/max(1, tp_coarsen(i)) + 1
                            dxtopo(i) = dxtopo(i)*max(1, tp_coarsen(i))
                            dytopo(i) = dytopo(i)*max(1, tp_coarsen(i))
                            tp_i_start(i) = i_start
                            tp_j_start(i) = j_start
                        end if
                    end if

                    ! Apply tp_x_shift to the stored domain extents so that all
                    ! coordinate comparisons (topoarea, rectintegral,
                    ! set_topo_for_dtopo) use the shifted bounds.  For types 2/3
                    ! this is the only place x_shift is applied; for type 4 we
                    ! also shift xlocs inside read_topo_file so the xstart
                    ! lookup against xlowtopo (= xll) still finds the right index.
                    if (tp_x_shift(i) /= 0.0d0) then
                        xlowtopo(i) = xlowtopo(i) + tp_x_shift(i)
                        xhitopo(i)  = xhitopo(i)  + tp_x_shift(i)
                    end if

                    topoID(i) = i
                    mtopo(i) = int(mxtopo(i), 8) * int(mytopo(i), 8)
                    areatopo(i) = dxtopo(i)*dytopo(i)
                enddo

                ! adding extra topo arrays corresponding spatially to dtopo
                ! arrays will be part of time dependent topowork as with all others.
                ! the state of these arrays at t=0 will be stored in topo0work
                topo0save(:)= 0
                do j=1,num_dtopo
                   i = mtopofiles + j
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
                   areatopo(i) = dxtopo(i)*dytopo(i)
                enddo

                ! Indexing into work array
                i0topo(1)=1
                if (mtopofiles+num_dtopo > 1) then
                    do i=2,mtopofiles+num_dtopo
                        i0topo(i)=i0topo(i-1) + mtopo(i-1)
                    enddo
                endif

                ! Read topography and allocate space for each file
                mtoposize = sum(mtopo)
                allocate(topowork(mtoposize))

                do i=1,mtopofiles
                    topoID(i) = i
                    topotime(i) = -huge(1.0)
                    call read_topo_file(mxtopo(i),mytopo(i),itopotype(i),topofname(i), &
                        xlowtopo(i),ylowtopo(i),topowork(i0topo(i):i0topo(i)+mtopo(i)-1), i)
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

                ! Topography priority order.
                ! mtopoorder(rank) = i means the i'th topo file has the given rank.
                ! rank=1 is the finest (highest priority); rank=mtopofiles is coarsest.
                ! rectintegral/topoarea use mtopoorder coarse-to-fine and finer files
                ! always overwrite coarser contributions in overlapping regions.
                !
                ! Convention: first file listed in topo.data = rank 1 = highest priority.
                ! Python's _compute_priority_order writes the finest file first, so the
                ! identity assignment mtopoorder(i)=i preserves that ordering without
                ! any Fortran-side re-sort. dtopo_for_topo files are inserted below.

                area_tol = 1.d-3  ! relative tol used in dtopo insertion below

                do i=1,mtopofiles
                    mtopoorder(i) = i
                end do

                ! now insert topo_for_dtopo arrays into ordering:

                do i=1,num_dtopo
                    idtopo = mtopofiles + i  ! file number of topo_for_dtopo
                    ! initialize topo_for_dtopo at lowest priority:
                    irank = mtopofiles + i
                    mtopoorder(irank) = irank

                    ! now bump up in priority as long as its cell area dx*dy
                    ! remains smaller than the maximum cell area of all
                    ! higher priority (lower rank) topos:

                    do jrank=mtopofiles+i-1, 1, -1
                        ! compute max area dx*dy over the top jrank topofiles
                        area_maxj = areatopo(mtopoorder(1))
                        do krank=2,jrank
                            area_maxj = max(area_maxj, &
                                            areatopo(mtopoorder(krank)))
                        enddo

                        if (areatopo(idtopo) < (1.d0-area_tol)*area_maxj) then
                            ! better resolution than first jrank, bump it up
                            ! (with fudge factor so it must be clearly better)
                            irank = irank - 1
                            mtopoorder(irank+1) = mtopoorder(irank) ! shift
                            mtopoorder(irank) = mtopofiles + i      ! insert
                        endif
                    enddo
                enddo

                ! update mtopofiles to include topo_for_dtopo as well:
                mtopofiles = mtopofiles + num_dtopo

                write(GEO_PARM_UNIT,*) ' '
                write(GEO_PARM_UNIT,*) 'Ranking of topography files', &
                    '  (including topo_for_dtopo) finest to coarsest: '
                write(GEO_PARM_UNIT,*) '   (filenumber is order they appear in'&
                                    // ' setrun.py, topofiles first)'
                write(GEO_PARM_UNIT,*) ' '

  301           format('rank = ',i3,'  filenumber = ',i3,'  dx*dy = ',e16.6)
                do rank=1,mtopofiles
                    write(GEO_PARM_UNIT,301) rank, mtopoorder(rank), &
                            areatopo(mtopoorder(rank))
                enddo
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
    !  apply_align_buffer_coarsen
    !
    !  Adjust a kept index window [i_lo, i_hi] (0-based file indices along one
    !  axis) for the align, buffer, and coarsen preprocessing attributes,
    !  mirroring Python Topography.crop() exactly (topotools.py):
    !
    !    1. align (only when coarsening): among the first c candidate start
    !       offsets, pick the one minimizing the fractional misalignment of
    !       (x - align_val) / (dxf*c); then trim i_hi onto the kept lattice.
    !    2. buffer: expand the window by nbuf*c points each side, clamped to
    !       the available range (clamping re-anchors the lattice, as Python's
    !       slicing does).
    !    3. coarsen: snap i_hi to the last kept lattice point
    !       (i_hi - i_lo becomes a multiple of c).
    !
    !  The kept points are then i_lo + k*c for k = 0..(i_hi-i_lo)/c.
    !
    !  x0 is the coordinate of file index 0 and dxf the signed spacing in
    !  file index order (negative for descending axes, e.g. N-to-S latitude).
    ! ========================================================================
    subroutine apply_align_buffer_coarsen(i_lo, i_hi, n_avail, x0, dxf, &
                                          align_val, has_align, c, nbuf)

        implicit none

        integer, intent(inout) :: i_lo, i_hi
        integer, intent(in) :: n_avail, c, nbuf
        real(kind=8), intent(in) :: x0, dxf, align_val
        logical, intent(in) :: has_align

        integer :: k, ioff
        real(kind=8) :: dx_new, r, frac, best

        if (c > 1 .and. has_align) then
            dx_new = dxf * c
            best = huge(1.0d0)
            ioff = 0
            do k = 0, c - 1
                r = (x0 + real(i_lo + k, 8)*dxf - align_val) / dx_new
                frac = abs(r - nint(r))
                if (frac < best) then
                    best = frac
                    ioff = k
                end if
            end do
            i_lo = i_lo + ioff
            i_hi = i_hi - mod(i_hi - i_lo, c)
        end if

        i_lo = max(0, i_lo - nbuf*c)
        i_hi = min(n_avail - 1, i_hi + nbuf*c)

        ! Snap the upper index onto the kept lattice.
        i_hi = i_lo + ((i_hi - i_lo)/c)*c

    end subroutine apply_align_buffer_coarsen

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

    subroutine read_topo_file(mx,my,topo_type,fname,xll,yll,topo,topo_idx)

#ifdef NETCDF
        use netcdf
#endif

        use geoclaw_module
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,topo_type
        character(len=512), intent(in) :: fname
        real(kind=8), intent(inout) :: topo(1:int(mx, 8)*int(my, 8))
        real(kind=8), intent(in) :: xll,yll
        integer, intent(in), optional :: topo_idx

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        logical, parameter :: maketype2 = .false.
        integer :: missing,status,n
        real(kind=8) :: no_data_value,x,y,topo_temp
        real(kind=8) :: values(16)
        character(len=80) :: str
        integer(kind=8) :: i, j, mtot
        ! For types 2/3 crop/coarsen extraction
        real(kind=8), allocatable :: topo_full(:)
        logical :: crop_active
        integer :: mx_full, my_full, cf
        integer(kind=8) :: mtot_full, icrop, jcrop, i_in_full, j_in_full

        ! NetCDF Support
        character(len=64) :: direction, x_dim_name, x_var_name, y_dim_name, &
            y_var_name, z_var_name, var_name
        real(kind=8), allocatable :: nc_buffer(:, :), xlocs(:), ylocs(:)
        integer(kind=4) :: x_var_id, y_var_id, z_var_id, x_dim_id, y_dim_id
        integer(kind=4) :: xstart(1), ystart(1), mx_tot, my_tot
        integer(kind=4) :: ios, nc_file, dim_ids(2), num_dims, &
            var_type, num_vars, num_dims_tot, z_dim_ids(2)
        logical :: lat_south_to_north

        mtot = int(mx, 8) * int(my, 8)

        print *, ' '
        print *, 'Reading topography file  ', trim(fname)
        print *

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
                ! Skip 5-line header; read no_data_value from line 6.
                do i=1,5
                    read(iunit,*)
                enddo
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                no_data_value = values(1)

                missing = 0

                if (present(topo_idx)) then
                    crop_active = any(tp_crop_extent(:,topo_idx) /= 0.0d0) &
                                  .or. tp_coarsen(topo_idx) > 1
                else
                    crop_active = .false.
                end if

                if (crop_active) then
                    ! ----------------------------------------------------------
                    ! Crop and/or coarsen active: read the full file into a
                    ! temporary array, then extract the kept lattice computed
                    ! by read_topo_settings (tp_i_start/tp_j_start, stride
                    ! tp_coarsen) into the output topo.
                    ! ----------------------------------------------------------
                    mx_full  = mxtopo_full(topo_idx)
                    my_full  = mytopo_full(topo_idx)
                    mtot_full = int(mx_full,8) * int(my_full,8)
                    allocate(topo_full(mtot_full))

                    select case(abs(topo_type))
                        case(2)
                            do i = 1, mtot_full
                                read(iunit,*) topo_full(i)
                                if (topo_full(i) == no_data_value) then
                                    missing = missing + 1
                                    topo_full(i) = topo_missing
                                end if
                            end do
                        case(3)
                            do j = 1, int(my_full,8)
                                read(iunit,*) (topo_full((j-1)*int(mx_full,8)+i), &
                                               i=1,mx_full)
                                do i = 1, int(mx_full,8)
                                    if (topo_full((j-1)*int(mx_full,8)+i) &
                                            == no_data_value) then
                                        missing = missing + 1
                                        topo_full((j-1)*int(mx_full,8)+i) = topo_missing
                                    end if
                                end do
                            end do
                    end select

                    ! Extract the kept mx-by-my lattice from topo_full.  The
                    ! output topo rows run N-to-S while tp_j_start counts from
                    ! the south, so output row jcrop corresponds to the kept
                    ! south-based row tp_j_start + (my - jcrop)*c, which is
                    ! file row my_full - that (file rows are 1-based from the
                    ! north).
                    cf = max(1, tp_coarsen(topo_idx))
                    do jcrop = 1, int(my, 8)
                        j_in_full = int(my_full,8) - (int(tp_j_start(topo_idx),8) &
                                    + (int(my,8) - jcrop)*int(cf,8))
                        do icrop = 1, int(mx, 8)
                            i_in_full = int(tp_i_start(topo_idx),8) &
                                        + (icrop - 1)*int(cf,8) + 1
                            topo((jcrop-1)*int(mx,8) + icrop) = &
                                topo_full((j_in_full-1)*int(mx_full,8) + i_in_full)
                        end do
                    end do

                    deallocate(topo_full)

                else
                    ! ----------------------------------------------------------
                    ! No crop: read directly into the output topo array.
                    ! ----------------------------------------------------------
                    select case(abs(topo_type))
                        case(2)
                            do i=1,mtot
                                read(iunit,*) topo(i)
                                if (topo(i) == no_data_value) then
                                    missing = missing + 1
                                    topo(i) = topo_missing
                                    ! uncomment next line to print row i
                                    ! write(6,600) i
 600                                format('*** missing data, i = ',i6)
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
 601                                    format('*** missing data, i = ',i6, &
                                               '  j = ',i6)
                                    endif
                                enddo
                            enddo
                    end select
                end if

                ! Write a warning if we found any missing values
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

                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, num_vars))

                ! ----------------------------------------------------------------
                ! Dimension discovery: descriptor or fallback heuristic.
                ! ----------------------------------------------------------------
                if (present(topo_idx) .and. &
                        len_trim(nc_lon_name(topo_idx)) > 0) then
                    call check_netcdf_error(nf90_inq_dimid(nc_file, &
                        trim(nc_lon_name(topo_idx)), x_dim_id))
                    call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                        x_dim_id, x_dim_name, mx_tot))
                    call check_netcdf_error(nf90_inq_dimid(nc_file, &
                        trim(nc_lat_name(topo_idx)), y_dim_id))
                    call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                        y_dim_id, y_dim_name, my_tot))
                else
                    call get_dim_info(nc_file, num_dims_tot, x_dim_id, &
                        x_dim_name, mx_tot, y_dim_id, y_dim_name, my_tot)
                end if

                ! Dimension IDs and variable IDs are separate ID spaces in NetCDF.
                ! Look up the coordinate variable IDs by name so nf90_get_var
                ! receives the correct IDs regardless of creation order in the file.
                call check_netcdf_error(nf90_inq_varid(nc_file, trim(x_dim_name), x_var_id))
                call check_netcdf_error(nf90_inq_varid(nc_file, trim(y_dim_name), y_var_id))

                allocate(xlocs(mx_tot), ylocs(my_tot))
                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xlocs, &
                    start=(/ 1 /), count=(/ mx_tot /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, ylocs, &
                    start=(/ 1 /), count=(/ my_tot /)))

                ! Apply lon_offset then tp_x_shift to convert file coordinates
                ! to domain coordinates, so xstart matches the domain-coord xll
                ! that was computed by read_topo_header (which also has both
                ! offsets applied via read_topo_settings).
                if (present(topo_idx)) then
                    if (nc_lon_offset(topo_idx) /= 0.0d0) then
                        xlocs = xlocs + nc_lon_offset(topo_idx)
                    end if
                    if (tp_x_shift(topo_idx) /= 0.0d0) then
                        xlocs = xlocs + tp_x_shift(topo_idx)
                    end if
                end if

                ! Coarsening factor: the kept points are read with a strided
                ! nf90_get_var, so mx/my count kept points and the spacing
                ! between them is cf times the file spacing.
                if (present(topo_idx)) then
                    cf = max(1, tp_coarsen(topo_idx))
                else
                    cf = 1
                end if

                xstart = minloc(xlocs, mask=(xlocs == xll))
                lat_south_to_north = (ylocs(1) < ylocs(my_tot))
                if (lat_south_to_north) then
                    ! S-to-N: southern boundary (yll) is at a low index.
                    ystart = minloc(ylocs, mask=(ylocs == yll))
                else
                    ! N-to-S: the read must start at the northern boundary (yhi),
                    ! which sits at a low index in the descending array.
                    ! Reconstruct yhi from yll, my, and the kept-point spacing.
                    y = yll + (my - 1) * cf * abs(ylocs(2) - ylocs(1))
                    ystart = minloc(ylocs, &
                        mask=(abs(ylocs - y) < 0.5d0 * abs(ylocs(2) - ylocs(1))))
                end if
                deallocate(xlocs, ylocs)

                ! ----------------------------------------------------------------
                ! Find the topo data variable: descriptor or first 2-D variable.
                ! ----------------------------------------------------------------
                if (present(topo_idx) .and. &
                        len_trim(nc_var_name(topo_idx)) > 0) then
                    call check_netcdf_error(nf90_inq_varid(nc_file, &
                        trim(nc_var_name(topo_idx)), z_var_id))
                    call check_netcdf_error(nf90_inquire_variable(nc_file, &
                        z_var_id, var_name, var_type, num_dims, z_dim_ids))
                else
                    z_var_id = -1
                    do n = 1, num_vars
                        call check_netcdf_error(nf90_inquire_variable(nc_file, &
                            n, var_name, var_type, num_dims, dim_ids))
                        if (var_name /= x_dim_name .and. &
                            var_name /= y_dim_name .and. num_dims == 2) then
                            z_var_id = n
                            z_dim_ids = dim_ids
                            exit
                        end if
                    end do
                    if (z_var_id == -1) then
                        print *, "ERROR: Unable to find topo variable in ", trim(fname)
                        stop
                    end if
                end if

                ! ----------------------------------------------------------------
                ! Load data into topo array (only if region overlaps domain).
                ! ----------------------------------------------------------------
                if ((mx > 0) .and. (my > 0)) then
                    if (z_dim_ids(1) == x_dim_id) then
                        allocate(nc_buffer(mx, my))
                    else if (z_dim_ids(1) == y_dim_id) then
                        allocate(nc_buffer(my, mx))
                    else
                        print *, "ERROR: NetCDF z variable dimensions do not align with x/y"
                        stop
                    end if

                    if (z_dim_ids(1) == x_dim_id) then
                        ! Variable stored (lon, lat): start/count already lon-first.
                        call check_netcdf_error(nf90_get_var(nc_file, z_var_id, &
                            nc_buffer, start=(/ xstart(1), ystart(1) /), &
                            count=(/ mx, my /), stride=(/ cf, cf /)))
                    else
                        ! Variable stored (lat, lon): swap start and count to match
                        ! the actual dimension order of the variable.
                        call check_netcdf_error(nf90_get_var(nc_file, z_var_id, &
                            nc_buffer, start=(/ ystart(1), xstart(1) /), &
                            count=(/ my, mx /), stride=(/ cf, cf /)))
                    end if

                    ! Transpose into GeoClaw topo array (N-to-S row order).
                    ! For S-to-N files ystart points to the southern end of the
                    ! domain, so nc_buffer row 1 = south, row my = north; reverse
                    ! with (my - j) to store north first.
                    ! For N-to-S files ystart points to the northern end, so
                    ! nc_buffer row 1 = north already; copy in order with (j + 1).
                    if (z_dim_ids(1) == x_dim_id) then
                        if (lat_south_to_north) then
                            do j = 0, int(my, 8) - 1
                                topo(j*int(mx,8)+1:(j+1)*int(mx,8)) = nc_buffer(:, my - j)
                            end do
                        else
                            do j = 0, int(my, 8) - 1
                                topo(j*int(mx,8)+1:(j+1)*int(mx,8)) = nc_buffer(:, j + 1)
                            end do
                        end if
                    else
                        if (lat_south_to_north) then
                            do j = 0, int(my, 8) - 1
                                topo(j*int(mx,8)+1:(j+1)*int(mx,8)) = nc_buffer(my - j, :)
                            end do
                        else
                            do j = 0, int(my, 8) - 1
                                topo(j*int(mx,8)+1:(j+1)*int(mx,8)) = nc_buffer(j + 1, :)
                            end do
                        end if
                    end if
                    deallocate(nc_buffer)

                    ! ----------------------------------------------------------------
                    ! Fill-value check (descriptor only).
                    ! Python has already verified no fill values in the crop region,
                    ! but we keep this as a safety net for edge-case mismatches.
                    ! ----------------------------------------------------------------
                    if (present(topo_idx) .and. nc_has_fill(topo_idx)) then
                        do i = 1, mtot
                            if (abs(topo(i) - nc_fill_value(topo_idx)) < &
                                    1.0d-6 * max(abs(nc_fill_value(topo_idx)), 1.0d0)) then
                                if (trim(nc_fill_action(topo_idx)) == 'abort') then
                                    print *, "ERROR: Fill value found in topo file ", &
                                        trim(fname), " at index ", i
                                    print *, "  fill_value = ", nc_fill_value(topo_idx)
                                    print *, "  Set fill_action = warn to replace with topo_missing."
                                    stop
                                else
                                    topo(i) = topo_missing
                                end if
                            end if
                        end do
                    end if

                    ! Honour 'positive = down' attribute (flip sign for depth-positive files).
                    ! Only flip when topo_type > 0 to avoid double-flip for negative types.
                    ios = nf90_get_att(nc_file, z_var_id, 'positive', direction)
                    if (ios == NF90_NOERR) then
                        if (to_lower(direction) == "down" .and. topo_type > 0) then
                            topo(1:mtot) = -topo(1:mtot)
                        end if
                    end if
                end if

                call check_netcdf_error(nf90_close(nc_file))
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
        ! Per-file preprocessing applied in-memory after reading (types 2,3,4,5).
        ! Both Python (Topography.read()) and Fortran apply these independently
        ! so the original file on disk is never modified.
        ! ====================================================================
        if (present(topo_idx) .and. abs(topo_type) /= 1) then

            ! 1. negate_z: explicit Z sign flip, independent of topo_type sign.
            if (tp_negate_z(topo_idx)) then
                topo(1:mtot) = -topo(1:mtot)
            end if

            ! 2. z_shift: bulk datum adjustment.  Skip fill/missing cells.
            if (tp_z_shift(topo_idx) /= 0.0d0) then
                do i = 1, mtot
                    if (abs(topo(i) - topo_missing) > &
                            1.0d-6 * max(abs(topo_missing), 1.0d0)) then
                        topo(i) = topo(i) + tp_z_shift(topo_idx)
                    end if
                end do
            end if

            ! 3. x_shift: applied to xlowtopo/xhitopo in read_topo_settings
            !    (after read_topo_header) and to xlocs in the type-4 case above.
            !    No action needed here on the z array.

            ! 4. crop_extent/align/buffer/coarsen: index window computed in
            !    read_topo_settings (types 2/3, extracted from the full read
            !    above) and in read_topo_header (type 4, applied via strided
            !    NetCDF reads).  No action needed here on the z array.

        end if

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
    subroutine read_topo_header(fname,topo_type,mx,my,xll,yll,xhi,yhi,dx,dy,topo_idx)
#ifdef NETCDF
        use netcdf
#endif

        use geoclaw_module
        use utility_module, only: parse_values, to_lower

        implicit none

        ! Input and Output
        character(len=*), intent(in) :: fname
        integer, intent(in) :: topo_type
        integer, intent(out) :: mx, my
        real(kind=8), intent(out) :: xll, yll, xhi, yhi, dx, dy
        integer, intent(in), optional :: topo_idx

        ! Local
        integer, parameter :: iunit = 19
        integer :: topo_size, status, n
        real(kind=8) :: x,y,z,nodata_value
        logical :: found_file
        real(kind=8) :: values(16)
        character(len=80) :: str
        logical :: verbose
        logical :: xll_registered, yll_registered
        ! Crop/align/buffer/coarsen locals (type 4 preprocessing pipeline)
        real(kind=8) :: x1_c, x2_c, y1_c, y2_c, al_x4, al_y4
        integer :: cf4, nbuf4, i_lo4, i_hi4, j_lo4, j_hi4
        logical :: has_al4

        ! NetCDF Support
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

                ! ----------------------------------------------------------------
                ! Discover spatial dimensions.
                ! When a descriptor was provided (topo_idx present, lon_name set)
                ! look up dimension IDs by name directly.  Otherwise fall back to
                ! the name-heuristic in get_dim_info.
                ! ----------------------------------------------------------------
                call check_netcdf_error(nf90_inquire(nc_file, num_dims_tot, num_vars))

                if (present(topo_idx) .and. &
                        len_trim(nc_lon_name(topo_idx)) > 0) then
                    ! Descriptor path: exact name look-up
                    call check_netcdf_error(nf90_inq_dimid(nc_file, &
                        trim(nc_lon_name(topo_idx)), x_dim_id))
                    call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                        x_dim_id, x_dim_name, mx))
                    call check_netcdf_error(nf90_inq_dimid(nc_file, &
                        trim(nc_lat_name(topo_idx)), y_dim_id))
                    call check_netcdf_error(nf90_inquire_dimension(nc_file, &
                        y_dim_id, y_dim_name, my))
                else
                    ! Fallback: auto-discover by name matching
                    call get_dim_info(nc_file, num_dims_tot, x_dim_id, &
                        x_dim_name, mx, y_dim_id, y_dim_name, my)
                end if

                allocate(xlocs(mx), ylocs(my), x_in_dom(mx), y_in_dom(my))

                ! ----------------------------------------------------------------
                ! Look up coordinate variable IDs (CF convention: var name == dim name).
                ! ----------------------------------------------------------------
                call check_netcdf_error(nf90_inq_varid(nc_file, trim(x_dim_name), x_var_id))
                call check_netcdf_error(nf90_inq_varid(nc_file, trim(y_dim_name), y_var_id))

                ! ----------------------------------------------------------------
                ! Find the topo data variable.
                ! Descriptor path: look up by name.  Fallback: first 2-D variable
                ! that is not a coordinate variable.
                ! ----------------------------------------------------------------
                if (present(topo_idx) .and. &
                        len_trim(nc_var_name(topo_idx)) > 0) then
                    call check_netcdf_error(nf90_inq_varid(nc_file, &
                        trim(nc_var_name(topo_idx)), z_var_id))
                else
                    z_var_id = -1
                    do n = 1, num_vars
                        call check_netcdf_error(nf90_inquire_variable(nc_file, &
                            n, var_name, var_type, num_dims, dim_ids))
                        if (var_name /= x_dim_name .and. &
                            var_name /= y_dim_name .and. num_dims == 2) then
                            z_var_id = n
                            exit
                        end if
                    end do
                    if (z_var_id == -1) then
                        print *, "ERROR: Cannot find a 2-D topo variable in ", trim(fname)
                        print *, "  Specify var_name in the topo descriptor block."
                        stop
                    end if
                end if

                ! ----------------------------------------------------------------
                ! Read coordinate arrays.
                ! ----------------------------------------------------------------
                call check_netcdf_error(nf90_get_var(nc_file, x_var_id, xlocs, &
                    start=(/ 1 /), count=(/ mx /)))
                call check_netcdf_error(nf90_get_var(nc_file, y_var_id, ylocs, &
                    start=(/ 1 /), count=(/ my /)))

                dx = abs(xlocs(2) - xlocs(1))
                dy = abs(ylocs(2) - ylocs(1))

                ! ----------------------------------------------------------------
                ! Select the subset of the file to load.
                ! Priority: nc_crop_bounds (file coords) > tp_crop_extent (domain
                ! coords, converted to file coords) > AMR domain fallback.
                ! ----------------------------------------------------------------
                nbuf4 = 0
                if (present(topo_idx) .and. nc_has_crop(topo_idx)) then
                    ! nc_crop_bounds are in FILE coordinates; compare against
                    ! file-coord xlocs before applying lon_offset.
                    x_in_dom = (xlocs >= nc_crop_bounds(1, topo_idx)) .and. &
                               (xlocs <= nc_crop_bounds(2, topo_idx))
                    y_in_dom = (ylocs >= nc_crop_bounds(3, topo_idx)) .and. &
                               (ylocs <= nc_crop_bounds(4, topo_idx))
                else if (present(topo_idx) .and. &
                         any(tp_crop_extent(:,topo_idx) /= 0.0d0)) then
                    ! tp_crop_extent is in domain coordinates (after nc_lon_offset
                    ! and tp_x_shift).  Convert to file coordinates for comparison
                    ! against raw xlocs (which have neither offset applied yet).
                    ! tp_buffer is applied in index space below, mirroring
                    ! Python Topography.crop().
                    x1_c = tp_crop_extent(1,topo_idx) - nc_lon_offset(topo_idx) &
                           - tp_x_shift(topo_idx)
                    x2_c = tp_crop_extent(2,topo_idx) - nc_lon_offset(topo_idx) &
                           - tp_x_shift(topo_idx)
                    y1_c = tp_crop_extent(3,topo_idx)
                    y2_c = tp_crop_extent(4,topo_idx)
                    x_in_dom = (xlocs >= x1_c) .and. (xlocs <= x2_c)
                    y_in_dom = (ylocs >= y1_c) .and. (ylocs <= y2_c)
                    nbuf4 = tp_buffer(topo_idx)
                else
                    x_in_dom = (xlocs > (xlower - dx - hxposs(1)*nghost)) .and. &
                               (xlocs < (xupper + dx + hxposs(1)*nghost))
                    y_in_dom = (ylocs > (ylower - dy - hyposs(1)*nghost)) .and. &
                               (ylocs < (yupper + dy + hyposs(1)*nghost))
                end if

                ! Apply lon_offset AFTER crop_bounds selection: converts file
                ! coordinates to domain coordinates.  crop_bounds above were
                ! compared against file-coord xlocs; xll/xhi below will be in
                ! domain coordinates as GeoClaw expects.
                if (present(topo_idx)) then
                    if (nc_lon_offset(topo_idx) /= 0.0d0) then
                        xlocs = xlocs + nc_lon_offset(topo_idx)
                    end if
                end if

                if (count(x_in_dom) == 0 .or. count(y_in_dom) == 0) then
                    ! No overlap with the requested region: report an empty
                    ! topo; read_topo_file skips loading when mx/my = 0.
                    mx = 0
                    my = 0
                    xll = 0.0d0
                    yll = 0.0d0
                    xhi = 0.0d0
                    yhi = 0.0d0
                else
                    ! Convert the selection masks to index windows (0-based)
                    ! and apply the align/buffer/coarsen pipeline, mirroring
                    ! Python Topography.crop().  The masks select contiguous
                    ! runs, so first/last true bound the window.
                    if (present(topo_idx)) then
                        cf4 = max(1, tp_coarsen(topo_idx))
                        has_al4 = any(tp_align(:,topo_idx) /= 0.0d0)
                        ! Align targets are in domain coordinates; xlocs now
                        ! include lon_offset but not tp_x_shift, so shift the
                        ! x target into the same frame.
                        al_x4 = tp_align(1,topo_idx) - tp_x_shift(topo_idx)
                        al_y4 = tp_align(2,topo_idx)
                    else
                        cf4 = 1
                        has_al4 = .false.
                        al_x4 = 0.0d0
                        al_y4 = 0.0d0
                    end if

                    i_lo4 = 0
                    do while (.not. x_in_dom(i_lo4 + 1))
                        i_lo4 = i_lo4 + 1
                    end do
                    i_hi4 = mx - 1
                    do while (.not. x_in_dom(i_hi4 + 1))
                        i_hi4 = i_hi4 - 1
                    end do
                    j_lo4 = 0
                    do while (.not. y_in_dom(j_lo4 + 1))
                        j_lo4 = j_lo4 + 1
                    end do
                    j_hi4 = my - 1
                    do while (.not. y_in_dom(j_hi4 + 1))
                        j_hi4 = j_hi4 - 1
                    end do

                    call apply_align_buffer_coarsen(i_lo4, i_hi4, mx, &
                        xlocs(1), xlocs(2) - xlocs(1), al_x4, has_al4, &
                        cf4, nbuf4)
                    call apply_align_buffer_coarsen(j_lo4, j_hi4, my, &
                        ylocs(1), ylocs(2) - ylocs(1), al_y4, has_al4, &
                        cf4, nbuf4)

                    ! Kept-lattice extents from the actual coordinate values
                    ! (min/max handles descending axes).  Spacings become the
                    ! kept-point spacings.
                    xll = min(xlocs(i_lo4 + 1), xlocs(i_hi4 + 1))
                    xhi = max(xlocs(i_lo4 + 1), xlocs(i_hi4 + 1))
                    yll = min(ylocs(j_lo4 + 1), ylocs(j_hi4 + 1))
                    yhi = max(ylocs(j_lo4 + 1), ylocs(j_hi4 + 1))
                    mx  = (i_hi4 - i_lo4)/cf4 + 1
                    my  = (j_hi4 - j_lo4)/cf4 + 1
                    dx  = dx*cf4
                    dy  = dy*cf4
                end if

                call check_netcdf_error(nf90_close(nc_file))
                deallocate(xlocs, ylocs, x_in_dom, y_in_dom)
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

! :TODO: Use the utility module's version of these
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
#endif

    ! ========================================================================
    ! subroutine read_netcdf_descriptor(iunit, idx)
    !
    ! Read the optional key=value block that follows the topo_type line in
    ! topo.data for a type-4 (NetCDF) entry.  Lines are consumed until either
    ! a blank line or EOF is reached.  Parsed values are stored in the module-
    ! level nc_* arrays at position idx.
    !
    ! Format (Python DescriptorWriter output):
    !   var_name       = z
    !   lon_name       = longitude
    !   lat_name       = latitude
    !   lon_offset     = 0.0
    !   lat_order      = N_to_S
    !   dim_order      = lat,lon
    !   fill_value     = -9999.0
    !   fill_action    = abort
    !   crop_bounds    = -180.0 180.0 -90.0 90.0
    !
    ! Unknown keys produce a warning and are skipped.
    ! ========================================================================
    subroutine read_netcdf_descriptor(iunit, idx)

        use geoclaw_module, only: GEO_PARM_UNIT

        implicit none

        integer, intent(in) :: iunit, idx

        character(len=512) :: line, key, val
        integer :: eq_pos, comment_pos, ios

        do
            read(iunit, '(a)', iostat=ios) line
            if (ios /= 0) exit          ! EOF

            ! Strip inline comment (! character)
            comment_pos = index(line, '!')
            if (comment_pos > 0) line = line(1:comment_pos-1)
            line = adjustl(trim(line))

            ! Blank line terminates the descriptor block
            if (len_trim(line) == 0) exit

            ! Skip lines without an '=' sign
            eq_pos = index(line, '=')
            if (eq_pos == 0) cycle

            key = adjustl(trim(line(1:eq_pos-1)))
            val = adjustl(trim(line(eq_pos+1:)))

            select case (trim(key))
                case ('var_name')
                    nc_var_name(idx) = trim(val)
                case ('lon_name')
                    nc_lon_name(idx) = trim(val)
                case ('lat_name')
                    nc_lat_name(idx) = trim(val)
                case ('lon_offset')
                    read(val, *) nc_lon_offset(idx)
                case ('lat_order')
                    nc_lat_order(idx) = trim(val)
                case ('dim_order')
                    nc_dim_order(idx) = trim(val)
                case ('fill_value')
                    read(val, *) nc_fill_value(idx)
                    nc_has_fill(idx) = .true.
                case ('fill_action')
                    nc_fill_action(idx) = trim(val)
                case ('crop_bounds')
                    read(val, *) nc_crop_bounds(1, idx), nc_crop_bounds(2, idx), &
                                 nc_crop_bounds(3, idx), nc_crop_bounds(4, idx)
                    nc_has_crop(idx) = .true.
                case default
                    write(GEO_PARM_UNIT, *) 'WARNING: Unknown NetCDF descriptor key: ', &
                        trim(key)
            end select
        end do

        ! Log what was parsed
        write(GEO_PARM_UNIT, *) '  NetCDF descriptor for topo file ', idx, ':'
        write(GEO_PARM_UNIT, *) '    var_name       = ', trim(nc_var_name(idx))
        write(GEO_PARM_UNIT, *) '    lon_name       = ', trim(nc_lon_name(idx))
        write(GEO_PARM_UNIT, *) '    lat_name       = ', trim(nc_lat_name(idx))
        write(GEO_PARM_UNIT, *) '    lon_offset     = ', nc_lon_offset(idx)
        write(GEO_PARM_UNIT, *) '    lat_order      = ', trim(nc_lat_order(idx))
        write(GEO_PARM_UNIT, *) '    fill_action    = ', trim(nc_fill_action(idx))
        if (nc_has_fill(idx)) &
            write(GEO_PARM_UNIT, *) '    fill_value     = ', nc_fill_value(idx)
        if (nc_has_crop(idx)) &
            write(GEO_PARM_UNIT, *) '    crop_bounds    = ', nc_crop_bounds(:, idx)

    end subroutine read_netcdf_descriptor

#ifdef NETCDF
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
