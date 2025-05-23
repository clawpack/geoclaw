
subroutine fgmax_values(mx,my,meqn,mbc,maux,q,aux,dx,dy, &
                   xlower,ylower,i1,i2,j1,j2,values)

    ! Given a grid q (and aux if needed), set the elements of 
    !   values(mv,i,j)  for mv=1:FG_NUM_VAL
    ! to the desired values that will be output and/or monitored on
    ! the fixed grid(s).
    !
    ! Only the elements with indices in [i1,i2,j1,j2] need be set.

    ! This library routine expects FG_NUM_VAL to be 1, 2, or 5 and sets:
    !   values(1,i,j) = h              (if FG_NUM_VAL >= 1)
    !   values(2,i,j) = speed          (if FG_NUM_VAL >= 2)
    !   values(3,i,j) = momentum       (if FG_NUM_VAL == 5)
    !   values(4,i,j) = momentum flux  (if FG_NUM_VAL == 5)
    !   values(5,i,j) = -depth         (if FG_NUM_VAL == 5)
    ! The max of -depth can be used to determin the minimum depth of water
    ! at a point over the computation, useful in harbors where ships may be
    ! grounded if the depth goes too low.


    use fgmax_module
    use geoclaw_module, only: sea_level, dry_tolerance

    implicit none
    integer, intent(in) :: mx,my,meqn,mbc,maux
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: dx,dy,xlower,ylower
    integer, intent(in) :: i1,i2,j1,j2
    real(kind=8), intent(inout) :: values(FG_NUM_VAL, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    real(kind=8) :: s_dry_tol
    integer :: i,j

    if ((FG_NUM_VAL.ne.1) .and. (FG_NUM_VAL.ne.2) .and. (FG_NUM_VAL.ne.5)) then
        write(6,*) '*** Error -- expecting FG_NUM_VAL = 1, 2, or 5'
        write(6,*) '***   in fgmax_values, found FG_NUM_VAL = ',FG_NUM_VAL
        stop
        endif
        
    s_dry_tol = dry_tolerance
    
    do i=i1,i2
        do j=j1,j2
            values(1,i,j) = q(1,i,j)
            if (FG_NUM_VAL > 1) then
                if (q(1,i,j) > s_dry_tol) then
                    values(2,i,j) = sqrt(q(2,i,j)**2+q(3,i,j)**2) / q(1,i,j)
                  else
                    values(2,i,j) = 0.d0
                  endif
                endif
            if (FG_NUM_VAL > 2) then
                values(3,i,j) = q(1,i,j)*values(2,i,j)       ! hs,   momentum
                values(4,i,j) = values(3,i,j)*values(2,i,j)  ! hs^2, mom flux
                values(5,i,j) = -q(1,i,j)                    ! for min h
                endif
            enddo
        enddo

end subroutine fgmax_values
