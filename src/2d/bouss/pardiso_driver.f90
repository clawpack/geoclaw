
subroutine pardiso_driver(soln,rhs,levelBouss,numBoussCells,nst,ncc,topo_finalized)

    ! use pardiso to solve sparse linear system
    ! called this way so can allocate storage of correct size
    
    use bouss_module
    implicit none
    
    integer, intent(in) :: levelBouss, numBoussCells, nst, ncc
    real(kind=8), intent(in) :: rhs(2*numBoussCells)
    logical, intent(in) :: topo_finalized
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8), intent(out) :: soln(2*numBoussCells)
    
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: k,i,j,mptr
    integer :: numCores, omp_get_max_threads

    integer :: icc(2*numBoussCells+1), ccc(ncc), sys
    integer*8 numeric, symbolic
    real(kind=8) :: acc(ncc), info(90)

    !! paridso vars
    !integer :: mtype,solver,maxfct,mnum,phase,nn,nrhs,error,msglvl
    !real*8 dparm(64)
    !integer*8 pt(64,maxlv)

#ifdef HAVE_PARDISO
    nrhs = 1
    iparm(1) = 0 ! use all default values, except for iparm(3) = # cores?

    numCores = 1 ! initialized to 1 in case not openmp 
!$  numCores = omp_get_max_threads() 
    iparm(3) = numCores

    maxfct = 1
    mnum = 1
    
    mtype = 11
    solver = 0
    
    minfo => matrix_info_allLevs(levelBouss)

    !!write(*,*)"second time"
    !! dump matrix_ia
    !!do k = 1, 10
    !!  write(*,901) k,minfo%matrix_ia(k)
 901  format(2i8)
    !!end do

    !================   Step 4 Solve matrix system =======================
            
    !   call timestamp()

        !! moved to setgrd call pardisoinit(pt(1,levelBouss), mtype, solver, iparm, dparm, error)
        !!if (error .ne. 0) then
         !! if (error.EQ.-10 ) WRITE(*,*) 'No license file found'
         !! if (error.EQ.-11 ) WRITE(*,*) 'License is expired'
         !! if (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
         !! STOP
        !!else
        !!   write(*,*)"license check successful"
        !!endif

        !! convert to compressed format. next call done in implicit_bouss
        ! call st_to_cc_size (nst,minfo%matrix_ia,minfo%matrix_ja,ncc)

        nn = 2*numBoussCells    ! number of rows and cols
        !!! for compressed row, reverse ia and ja, and output arrays too
        !!!call st_to_cc_index (nst,minfo%matrix_ia,minfo%matrix_ja,ncc,nn,icc,ccc)
        call st_to_cc_index (nst,minfo%matrix_ja,minfo%matrix_ia,ncc,nn,ccc,icc)
        !!!call st_to_cc_values(nst,minfo%matrix_ia,minfo%matrix_ja, minfo%matrix_sa,ncc,nn,icc,ccc,acc)
        call st_to_cc_values(nst,minfo%matrix_ja,minfo%matrix_ia, minfo%matrix_sa,ncc,nn,ccc,icc,acc)


        call pardiso_chkmatrix (mtype, nn, acc, icc,ccc, error)
        if (error .ne. 0) then
           write(*,*)"the following error in chk_matrix was found ",error
           stop
        endif

        ! ..  pardiso_printstats(...)
        !     prints information on the matrix to STDOUT.
        !     Use this functionality only for debugging purposes
       !CALL pardiso_printstats (mtype, nn, acc, icc, ccc, nrhs, rhs, error);
       !if (error .NE. 0) then
       !   write(*,*) 'The following ERROR was detected: ', error
       !  stop
       ! endif

        !! reordering and symbolic factorization. Also
        !! allocates all memory needed for factorization
        if (newGrids(levelBouss) .or. .not. topo_finalized) then
           newGrids(levelBouss) = .false.
           phase = 11
           msglvl = 0
           !write(*,*)"before phase 11 call"
           !!call timestamp()
           call pardiso(pt(1,levelBouss), maxfct, mnum, mtype, phase,nn,acc,icc,ccc,  &
                        idum, nrhs, iparm, msglvl, rhs, soln, error, dparm)
                     !  idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)

           IF (error .NE. 0) THEN
              WRITE(*,*) 'The following ERROR was detected: ', error
              STOP
            END IF


           !WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
           !WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

           !  factorization
           phase     = 22  ! only factorization
           !write(*,*)"before phase 22 call"
           !!call timestamp()
           CALL pardiso (pt(1,levelBouss), maxfct, mnum, mtype, phase, nn, acc, icc, ccc,      &
                       idum, nrhs, iparm, msglvl, rhs, soln, error,dparm)
                    !  idum, nrhs, iparm, msglvl, ddum, ddum, error,dparm)

           !WRITE(*,*) 'Factorization completed ...  '
           IF (error .NE. 0) THEN
              WRITE(*,*) 'The following ERROR was detected: ', error
             STOP
           ENDIF 

       endif

 !.. Back substitution and iterative refinement
      phase     = 33  ! only solve
      iparm(8)  = 5   ! max numbers of iterative refinement steps

      !write(*,*)"before phase 33 call"
      !!call timestamp()
      CALL pardiso (pt(1,levelBouss), maxfct, mnum, mtype, phase, nn, acc, icc, ccc,    &
                    idum, nrhs, iparm, msglvl, rhs, soln, error, dparm) 
      !write(*,*)"after phase 33 call"
      !!call timestamp()
      !WRITE(*,*) 'Solve completed ... '
      !write(*,*) "num iterative refinements = iparm(7) =  ", iparm(7)
      avgIterRef = avgIterRef + iparm(7)
      countIterRef = countIterRef + 1
      maxIterRef = max(maxIterRef, iparm(7))
     
!      WRITE(*,*) 'The solution of the system is '
!      DO i = 1, n
!        WRITE(*,*) ' x(',i,') = ', soln(i)
!      END DO

      ! Termination and release of memory
      !phase     = -1           ! release internal memory
      ! moved to amr2 after tick for final release
      !!CALL pardiso (pt(1,levelBouss), maxfct, mnum, mtype, phase, nn, ddum, idum, idum,    &
      !!              idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
#endif

! if no pardiso this is dummy routine that just returns
    return
end subroutine pardiso_driver

