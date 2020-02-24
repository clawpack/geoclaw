module fgmax_module

    use amr_module, only: mxnest
    implicit none
    save


    type fgrid

        ! Derived data type for a single fixed grid.
        ! An array of objects of this type is used to keep track of all
        ! fixed grids in a computation.

        ! identifier number for this fixed grid:
        integer :: fgno ! now used in fgmax_grids.data
        
        integer :: point_style

        ! fixed grid points are (x(k),y(k)), for k=1:npts
        integer :: npts
        real(kind=8), allocatable, dimension(:) :: x,y

        ! arrays valuemax(mv,k), tmax(mv,k), arrival_time(k), aux(level,ma,k)
        ! where mv = 1:num_values  # number of quantities monitored
        !       ma = 1:num_aux
        !        k = 1:npts
        ! New faster version uses klist(1:npts, 0:maxthreads-1)
        !     to select which points overlap each patch, and
        !     t_last_updated(1:mxnest), by level rather than point.
        ! Allocation is done in fgmax_read.
        real(kind=8), allocatable, dimension(:,:,:) :: aux
        real(kind=8), allocatable, dimension(:,:) :: valuemax
        real(kind=8), allocatable, dimension(:,:) :: tmax
        real(kind=8), allocatable, dimension(:) :: t_last_updated
        real(kind=8), allocatable, dimension(:) :: arrival_time
        integer, allocatable, dimension(:) :: levelmax
        integer, allocatable, dimension(:,:) :: klist

        ! time range to monitor:
        real(kind=8) :: tstart_max, tend_max

        ! Desired maximum delta t between updating valuemax:
        ! Will only update at start of time step if end of timestep is
        ! to far beyond last update time
        real(kind=8) :: dt_check

        ! Mininum level to check when updating valuemax or arrival times:
        ! Coarser levels will be ignored to avoid bad data from coarse grids
        integer :: min_level_check
        
        ! Mininum value of levelmax over the grid
        ! Set once each step on each level in advanc to avoid repeated minval's
        integer :: min_levelmax_checked

        ! Tolerance for flagging point for arrival of tsunami.
        ! Flag in fgmax_frompatch if eta tilde > sea_level + arrival_tol
        real(kind=8) :: arrival_tol


        ! Coordinates of corners of bounding box.
        ! This will be useful when generalizing to fgrids not aligned with x-y.
        real(kind=8) :: x1bb,x2bb,y1bb,y2bb

        ! keep track of whether all aux arrays have been computed on a given level:
        logical, allocatable, dimension(:) :: auxdone
        
        ! indicate if it is time to update a level:
        logical, allocatable, dimension(:) :: update_now

        ! new attributes for fgmax masked grids
        real(kind=8) :: xll,yll,xhi,yhi,dx,dy
        integer :: nx,ny
        integer, allocatable, dimension(:,:) :: index
        
        ! interp_method = 0 for pw constant, 1 for bilinear:
        integer :: interp_method
        
    end type

    logical, private :: module_setup = .false.

    ! declare array fgrids of fixed grids, each of type fgrid.
    ! allow them to be targets of pointers for shorthand in code.
    integer, parameter :: FG_MAXNUM_FGRIDS = 50  ! max number of fixed grids
    type(fgrid), target :: FG_fgrids(FG_MAXNUM_FGRIDS)

    ! special value to flag unset portions of arrays:
    real(kind=8), parameter :: FG_NOTSET = -0.99999d99

    ! unit to use for reading input data and writing fgrid results:
    integer, parameter :: FG_UNIT = 45

    ! number of max vals to monitor: Set in set_fgmax
    integer :: FG_NUM_VAL   
    ! number of aux vals to monitor
    integer, parameter :: FG_NUM_AUX = 1

    ! number of fixed grids in use (set by fgmax_read):
    integer :: FG_num_fgrids

    ! turn on debugging output to fort.61:
    logical, parameter :: FG_DEBUG = .false.


contains

    subroutine set_fgmax(fname)

        use amr_module, only: parmunit

        implicit none
        
        ! Subroutine arguments
        character(len=*), intent(in), optional :: fname
        
        ! Local storage
        integer, parameter :: unit = 7
        integer :: ifg
        character(len=150) :: fname_fg
        integer :: num_fgmax_grids, num_fgmax_val

        if (.not.module_setup) then

            write(parmunit,*) ' '
            write(parmunit,*) '--------------------------------------------'
            write(parmunit,*) 'SETFGMAX:'
            write(parmunit,*) '-----------'

            ! Open data file
            if (present(fname)) then
                call opendatafile(unit, fname)
            else
                call opendatafile(unit, 'fgmax_grids.data')
            endif

            ! Read in data
            read(unit,'(i2)') num_fgmax_val   ! name used in setrun.py
            FG_NUM_VAL = num_fgmax_val        ! module variable name
            read(unit,'(i2)') num_fgmax_grids ! name used in setrun.py
            FG_num_fgrids = num_fgmax_grids   ! module variable name
            
            if (FG_num_fgrids > FG_MAXNUM_FGRIDS) then
               write(6,601) FG_num_fgrids
     601       format('*** Too many fixed grids specified: FG_num_fgrids = ',i3,/, &
                   '*** Increase FG_MAXNUM_FGRIDS in fgmax_module.f90')
               stop
               endif

            do ifg=1,FG_num_fgrids
                ! new in v5.7.0: read all data from fgmax_grids.data:
                call fgmax_read(unit, ifg)
                enddo
            write(parmunit,*) ' '
            write(parmunit,*) '--------------------------------------------'
            write(parmunit,*) 'SETFGMAX:'
            write(parmunit,*) '-----------'

            write(parmunit,*) 'FG_NUM_VAL = ', FG_NUM_VAL
            write(parmunit,*) 'FG_num_fgrids = ', FG_num_fgrids

            module_setup = .true.
        end if

    
    end subroutine set_fgmax

end module fgmax_module
