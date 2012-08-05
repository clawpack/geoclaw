! Module containing refinement flagging criteria
module refinement_module

    use geoclaw_module, only: GEO_PARM_UNIT

    implicit none
    save

    ! ========================================================================
    !  Refinement Criteria
    ! ========================================================================
    real(kind=8), allocatable :: wave_tolerance(:)
    real(kind=8), allocatable :: speed_tolerance(:)
    real(kind=8) :: deep_depth
    integer :: max_level_deep
    
    ! ========================================================================
    !  Refinement Regions
    ! ========================================================================
    type region_type
        integer :: min_level,max_level
        real(kind=8) :: x_low,y_low,x_hi,y_hi,t_low,t_hi
    end type region_type
    
    integer :: num_regions
    type(region_type), allocatable :: regions(:)
    
    ! ========================================================================
    !  Flowgrades - Not updated yet, use at your own risk
    ! ========================================================================
    integer :: num_flowgrades
    real(kind=8), allocatable :: flowgradevalue(:)
    integer, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, allocatable :: iflowgrademinlevel(:)

contains
    
    ! =========================================================================
    !  Reads in the refinement control parameters
    ! =========================================================================
    subroutine set_refinement(file_name)
        
        use amr_module, only: mxnest
        use geoclaw_module, only: num_layers
        
        implicit none
        
        ! Arguments
        character(len=*), optional, intent(in) :: file_name
        
        ! Locals
        integer, parameter :: unit = 127
        integer :: i

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'Refinement Control Parameters:'
        write(GEO_PARM_UNIT,*) '------------------------------'

        if (present(file_name)) then
            call opendatafile(unit, file_name)
        else
            call opendatafile(unit, 'refinement.data')
        endif

        ! Basic criteria
        allocate(wave_tolerance(num_layers))
        read(unit,*) wave_tolerance
        allocate(speed_tolerance(mxnest))
        read(unit,*) (speed_tolerance(i),i=1,mxnest)
        read(unit,*) deep_depth
        read(unit,*) max_level_deep
        read(unit,*)
        
        ! Refinement region data
        read(unit,"(i2)") num_regions
        allocate(regions(num_regions))
        do i=1,num_regions
            read(unit,*) regions(i)%min_level, regions(i)%max_level, &
                         regions(i)%t_low, regions(i)%t_hi, &
                         regions(i)%x_low, regions(i)%x_hi, &
                         regions(i)%y_low, regions(i)%y_hi
        enddo
        close(unit)
        
        ! Write out data to parameter file
        write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
        write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
        write(GEO_PARM_UNIT,*) '   maxleveldeep:', max_level_deep
        write(GEO_PARM_UNIT,*) '   depthdeep:', deep_depth
        write(GEO_PARM_UNIT,*) ''
        write(GEO_PARM_UNIT,*) '  num_regions = ',num_regions
        write(GEO_PARM_UNIT,*) '  minlevel, maxlevel, tlow, thi, xlow, xhi, ylow, yhigh values:'
        do i=1,num_regions
            write(GEO_PARM_UNIT,*) regions(i)%min_level, regions(i)%max_level, &
                                   regions(i)%t_low, regions(i)%t_hi, &
                                   regions(i)%x_low, regions(i)%x_hi, &
                                   regions(i)%y_low, regions(i)%y_hi
        enddo
        
    end subroutine set_refinement
    
    
    ! =========================================================================
    ! TODO: This needs to be updated for the new module
    ! =========================================================================
    subroutine set_flow_grades(file_name)

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i

        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET FLOW GRADES:'
        write(GEO_PARM_UNIT,*) '------------'

        ! Read user parameters from setflowgrades.data
        if (present(file_name)) then
            call opendatafile(iunit, file_name)
        else
            call opendatafile(iunit, 'setflowgrades.data')
        endif
        
        read(iunit,*) num_flowgrades

        if (num_flowgrades == 0) then
            write(GEO_PARM_UNIT,*) '  No flow grades specified'
            return
        endif

        ! Allocate arrays
        allocate(flowgradevalue(num_flowgrades),iflowgradevariable(num_flowgrades))
        allocate(iflowgradetype(num_flowgrades),iflowgrademinlevel(num_flowgrades))

        do i=1,num_flowgrades
            read(iunit,*) flowgradevalue(i),iflowgradevariable(i), &
                iflowgradetype(i),iflowgrademinlevel(i)
        enddo

        close(iunit)

        write(GEO_PARM_UNIT,*) '   mflowgrades:',  num_flowgrades

        do i=1,num_flowgrades
            write(GEO_PARM_UNIT,"(d12.3,3i4)") flowgradevalue(i), &
                iflowgradevariable(i),iflowgradetype(i),iflowgrademinlevel(i)

        enddo

    end subroutine set_flow_grades
    
end module refinement_module