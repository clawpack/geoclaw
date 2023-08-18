
subroutine fgmax_finalize()

    ! Print out the maxval and aux arrays and de-allocate storage.
    
    ! New style introduced in v5.7.0:
    !   fgmax0001.txt, etc. instead of fort.FG1.valuemax, fort.FG1.aux
    !   aux(1,:,:) value (topography) on proper level included with values,
    !       instead of separate file giving aux values on all levels.

    use fgmax_module
    use amr_module, only: mxnest

    implicit none
    character(30) :: fname
    character(1) :: cma
    character(4) :: cfgno
    integer :: k,ifg,level,mv,ma,ipos,idigit,ifg1
    type(fgrid), pointer :: fg
    real(kind=8) :: topo

    if (FG_NUM_AUX .ne. 1) then
        write(6,*) '*** Unexpected FG_NUM_AUX = ',FG_NUM_AUX
        write(6,*) '*** Including only aux(1) in fort.FG file'
        endif
        
    cfgno = '0000'
    do ifg=1,FG_num_fgrids

        fg => FG_fgrids(ifg)   

        ifg1 = fg%fgno  ! use the fgno specified for filename, not ifg
        do ipos=4,1,-1
            idigit = mod(ifg1,10)
            cfgno(ipos:ipos) = char(ichar('0') + idigit)
            ifg1 = ifg1/10
            enddo

        !fname = 'fort.FG' // cfgno
        fname = 'fgmax' // cfgno // '.txt'
        print *, 'Writing to file ', fname
        open(unit=FG_UNIT,file=trim(fname),status='unknown',form='formatted')

        do k=1,fg%npts
            do mv=1,FG_NUM_VAL
                if (abs(fg%valuemax(mv,k)) .lt. 1.d-90) then
                    fg%valuemax(mv,k) = 0.d0
                    endif
                enddo
            ! Modified to print topo value after levelmax and before valuemax:
            if (fg%levelmax(k) > 0) then    
                topo = fg%aux(fg%levelmax(k),1,k)
            else
                topo = FG_NOTSET
            endif
            write(FG_UNIT,111) fg%x(k),fg%y(k), fg%levelmax(k), &
                  topo, &
                  (fg%valuemax(mv,k), mv=1,FG_NUM_VAL), &
                  (fg%tmax(mv,k), mv=1,FG_NUM_VAL), fg%arrival_time(k)
 111        format(2e20.11,i4,21e17.8)
            enddo

        close(FG_UNIT)


        ! deallocate(fg%valuemax,fg%levelmax,fg%aux,fg%x,fg%y)
        enddo
    if (FG_DEBUG) then
        write(6,*) '+++ Fixed grid debugging written to fort.61 and fort.65'   
        endif

end subroutine fgmax_finalize
