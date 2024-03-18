subroutine petsc_driver(soln,rhs_geo,levelBouss,numBoussCells,time,topo_finalized)

    ! use petsc to solve sparse linear system
    ! called this way so can allocate storage of correct size
    
    use bouss_module
#ifdef HAVE_PETSC
#include <petsc/finclude/petscksp.h>
!!#include "petscmat.h"
    use petscksp
    implicit none
    
    integer, intent(in) :: levelBouss, numBoussCells
    real(kind=8), intent(inout) :: rhs_geo(0:2*numBoussCells)
    real(kind=8), intent(in) :: time
    logical, intent(in) :: topo_finalized
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8), intent(out) :: soln(0:2*numBoussCells)
    !integer :: rowPtr(0:2*numBoussCells), cols(0:24*numBoussCells)
    !real(kind=8)  :: vals(0:24*numBoussCells)
    
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: numCores, omp_get_max_threads,i,nelt
    integer :: j
    logical :: idx

    ! petsc declaration
    integer :: n, nn, ierr
    logical :: flg
    integer*8 :: nz


    !! These PETSC declarations are needed  (moved to amr_module)
    !Mat J
    Vec rhs,solution
    !KSP ksp  ! linear solver ojbect
    PetscInt itnum

!   call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',n,flg,ierr);CHKERRA(ierr)


    minfo => matrix_info_allLevs(levelBouss)

    !================   Step 4 Solve matrix system =======================
            
       nn = 2*numBoussCells    ! number of rows and cols
       n = nn  ! petsc notation below
       nelt = minfo%matrix_nelt
       nz = nelt
       !write(17,*)" level ",levelBouss
       !write(17,777)(j+1,minfo%matrix_ia(j)+1,minfo%matrix_ja(j)+1,minfo%matrix_sa(j),j=0,nelt-1)
 777   format(3i8,e15.6)

       if (newGrids(levelBouss) .or. .not. topo_finalized) then
          ! if first time through nothing to destroy yet.
          ! otherwise destroy old stuff to recreate for new grid setup 
          ! turns out no need to call because petsc initializes and checks
          !if (time .gt. tstart_thisrun) then this was needed since if
          ! time == tstart_thisrun it could be a restart run and nothing to destroy
            call KSPDestroy(ksp(levelBouss),ierr)
            CHKERRA(ierr)
            call MatDestroy(Jr(levelBouss),ierr)
            CHKERRA(ierr)
          !endif

          newGrids(levelBouss) = .false.

          !petsc call to put in compressed row format
          idx = .true.   ! 1 based indexing, so not already decremented

          ! these lines are if reusing matrix with different entries but same structure
          if (.not. crs) then
            ! Using COO triplet sparse matrix storage: 
            call MatCreate(PETSC_COMM_SELF,Jr(levelBouss),ierr)
            CHKERRA(ierr)
            call MatSetSizes(Jr(levelBouss),n,n,n,n,ierr)
            CHKERRA(ierr)
            call MatSetType(Jr(levelBouss),MATSEQAIJ,ierr)
            CHKERRA(ierr)

            ! decrement indices since Petsc zero-based indexing
            do i = 1, nz
               minfo%matrix_ia(i) = minfo%matrix_ia(i) - 1
               minfo%matrix_ja(i) = minfo%matrix_ja(i) - 1
            end do

            call MatSetPreallocationCOO(Jr(levelBouss), nz, &
                                minfo%matrix_ia, minfo%matrix_ja,ierr)
            CHKERRA(ierr)
            ! put 1 indexing back in  
            do i = 1, nz
               minfo%matrix_ia(i) = minfo%matrix_ia(i) + 1
               minfo%matrix_ja(i) = minfo%matrix_ja(i) + 1
            end do
          else 
            ! using CRS sparse matrix format
            call MatCreateSeqAijWithArrays(PETSC_COMM_SELF,2*numBoussCells,  & 
                       2*numBoussCells,minfo%rowPtr,minfo%cols,  &
                       minfo%vals,Jr(levelBouss),ierr)
            CHKERRA(ierr)
            ! per Barry, tell petsc you have blocks of size 2
            call MatSetBlockSize(Jr(levelBouss),2,ierr)
            CHKERRA(ierr)
          endif

          ! both sparse matrix formats do this 
          call KSPCreate(PETSC_COMM_SELF,ksp(levelBouss),ierr)
          CHKERRA(ierr)
          call KSPSetErrorIfNotConverged(ksp(levelBouss),PETSC_TRUE,ierr)
          CHKERRA(ierr)
          call KSPSetFromOptions(ksp(levelBouss),ierr)
          CHKERRA(ierr)

          ! ABOVE CALLS SET OPTIONS FROM ENVIRONMENT VARS
          ! now TURN off any for level 1 you dont want
          call KSPSetReusePreconditioner(ksp(1),PETSC_FALSE,ierr)
          CHKERRA(ierr)
          ! if have previous solution tell PETSC about it
          !call KSPSetInitialGuessNonzero(ksp(levelBouss),PETSC_TRUE,ierr)
          call KSPSetOperators(ksp(levelBouss),Jr(levelBouss),Jr(levelBouss),ierr)
          CHKERRA(ierr)
      else ! when (.not. newGrids(levelBouss)) .and. topo_finalized
         ! just put in new values, reuse same sparse matrix structure
         ! next line notifies matrix has new vals
         ! and has saved pointer to vals in between calls
         call PetscObjectStateIncrease(Jr(levelBouss),ierr)
      endif

      ! this call reuse the matrix values instead of creating new ones
      if (.not. crs) then
         ! next line for COO triplet format
         call MatSetValuesCOO(Jr(levelBouss),minfo%matrix_sa,INSERT_VALUES,ierr)
         CHKERRA(ierr)
      endif

      ! put old soln in soln vec here if have one (WHAT DOES THIS MEAN?)

      if (.not. crs) then
        ! switch to 0-based indexing for COO triplet format
        ! (already 0-based for CRS)
        do i = 1, 2*numBousscells
           rhs_geo(i-1) = rhs_geo(i) 
        end do
      endif 
      
      call VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,rhs_geo,rhs,ierr)
      CHKERRA(ierr)
      call VecCreateSeqWithArray(PETSC_COMM_SELF,1,n,soln,solution,ierr)
      CHKERRA(ierr)

      ! save matrices for debugging. comment out to turn off
      !call MatView(Jr(levelBouss),PETSC_VIEWER_BINARY_SELF,ierr)


      call KSPSolve(ksp(levelBouss),rhs,solution,ierr)
      call KSPGetIterationNumber(ksp(levelBouss), itnum,ierr)
      CHKERRA(ierr)
      itcount(levelBouss) = itcount(levelBouss)+ itnum
      numTimes(levelBouss) = numTimes(levelBouss) + 1

      ! check if need to adjust soln back to 1 indexing
      if (.not. crs) then ! bump both  back up
        do i = 1, 2*numBoussCells
          soln(2*numBoussCells-i+1) = soln(2*numBoussCells-i)
          rhs_geo(2*numBoussCells-i+1) = rhs_geo(2*numBoussCells-i) 
        end do
      endif
      
      !! if itnum > itmax then
      !! call KSPGetResidualNorm(ksp(levelBouss),resmax,ierr)
      CHKERRA(ierr)
      ! petsc already puts solution into my array soln so just return

       ! create and destroy each step since fast
       call VecDestroy(rhs,ierr)
       CHKERRA(ierr)
       call VecDestroy(solution,ierr)
       CHKERRA(ierr)
#endif

end subroutine petsc_driver
