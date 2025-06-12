subroutine petsc_driver(soln,rhs_geo,levelBouss,numBoussCells,time,   &
                        aux_finalized,nelt)

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
    integer, intent(in) :: nelt, aux_finalized
    
    !! pass in numBoussCells for this level so can dimension this array
    real(kind=8), intent(out) :: soln(0:2*numBoussCells)
    
    type(matrix_levInfo),  pointer :: minfo
    
    integer :: numCores, omp_get_max_threads,i
    integer :: j
    integer :: local_ia(nelt),local_ja(nelt)

    ! petsc declaration
    integer :: n, nn, ierr
    logical :: flg
    integer*8 :: nz


    !! These PETSC declarations are needed  
    Vec rhs,solution
    Vec y
    PetscInt itnum
    KSPConvergedReason reason
    !! (moved the following to amr_module)
    !Mat J
    !KSP ksp  ! linear solver ojbect

#ifdef WHERE_AM_I
    write(*,*)" starting petsc_driver"
#endif

    minfo => matrix_info_allLevs(levelBouss)

    !================   Step 4 Solve matrix system =======================
            
       nn = 2*numBoussCells    ! number of rows and cols
       n = nn  ! petsc notation below
       !nelt = minfo%matrix_nelt
       nz = nelt
       !write(17,*)" level ",levelBouss
       !write(17,777)(j+1,minfo%matrix_ia(j)+1,minfo%matrix_ja(j)+1,minfo%matrix_sa(j),j=0,nelt-1)
 777   format(3i8,e15.6)

       if (newGrids(levelBouss) .or.  aux_finalized .lt. 2) then
          ! destroy old solver objects to recreate for new grid setup. No need 
          ! to check if init or restart because petsc checks if null
            call KSPDestroy(ksp(levelBouss),ierr)
            CHKERRA(ierr)
            call MatDestroy(Jr(levelBouss),ierr)
            CHKERRA(ierr)

          newGrids(levelBouss) = .false.

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
            ! and COO format is 1-based
            do i = 1, nz
               local_ia(i) = minfo%matrix_ia(i) - 1
               local_ja(i) = minfo%matrix_ja(i) - 1
            end do

            call MatSetPreallocationCOO(Jr(levelBouss), nz, &
                                local_ia, local_ja,ierr)
            CHKERRA(ierr)
            ! put 1 indexing back in  
            !do i = 1, nz
            !   minfo%matrix_ia(i) = minfo%matrix_ia(i) + 1
            !   minfo%matrix_ja(i) = minfo%matrix_ja(i) + 1
            !end do
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

          ! both sparse matrix formats create solver object
          call KSPCreate(PETSC_COMM_SELF,ksp(levelBouss),ierr)
          CHKERRA(ierr)
          call KSPSetErrorIfNotConverged(ksp(levelBouss),PETSC_TRUE,ierr)
          CHKERRA(ierr)
          ! Next call sets options from environment vars
          call KSPSetFromOptions(ksp(levelBouss),ierr)
          CHKERRA(ierr)

          ! now TURN off any for level 1 you dont want, such as level 1 saving
          ! preconditioner: since it never regrids would never change
          call KSPSetReusePreconditioner(ksp(1),PETSC_FALSE,ierr)
          CHKERRA(ierr)

          ! if want to try initiializing with previous solution this sets it up
          ! if have previous solution tell PETSC about it
          !call KSPSetInitialGuessNonzero(ksp(levelBouss),PETSC_TRUE,ierr)

          call KSPSetOperators(ksp(levelBouss),Jr(levelBouss),Jr(levelBouss),ierr)
          CHKERRA(ierr)
      else ! when (.not. newGrids(levelBouss)) .and. topo_finalized
         ! just put in new values, reuse same sparse matrix structure
         ! next line notifies matrix has new vals
         ! and has saved pointer to vals in between calls
         ! not needed for COO format
         if (crs) then
             call PetscObjectStateIncrease(Jr(levelBouss),ierr)
         endif
      endif

      if (.not. crs) then ! only for COO triplet format
         ! this call reuses the matrix values instead of creating new ones
         ! (for CRS, PetscObjectStsateIncrease is called instead above)
         call MatSetValuesCOO(Jr(levelBouss),minfo%matrix_sa,INSERT_VALUES,ierr)
         CHKERRA(ierr)
      endif

      ! if you want to use previous soln as initial guess,
      ! set soln to old soln here

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
      !call MatView(Jr(levelBouss),PETSC_VIEWER_STDOUT_SELF,ierr)
      !call MatView(Jr(levelBouss),PETSC_VIEWER_BINARY_SELF,ierr)


      call KSPSolve(ksp(levelBouss),rhs,solution,ierr)
      call KSPGetIterationNumber(ksp(levelBouss), itnum,ierr)
      !write(*,*)" took ",itnum," iterations"
      CHKERRA(ierr)
      itcount(levelBouss) = itcount(levelBouss)+ itnum
      numTimes(levelBouss) = numTimes(levelBouss) + 1

      !DEBUG
      !call vecDuplicate(solution,y,ierr)  ! makes new vec y with same size as soln
      !call MatMult(Jr(levelBouss),solution,y,ierr)
      !call VecView(rhs,PETSC_VIEWER_STDOUT_SELF,ierr)
      !call VecView(y,PETSC_VIEWER_STDOUT_SELF,ierr)

      ! still working on this aspect. 
      !call KSPGetConvergedReason(ksp(levelBouss),reason,ierr)
      ! negative reason is bad, leave 
      !if (reason < 0) then
      !write(*,*) 'level ',levelBouss,' reason ',reason
      !endif

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

#ifdef WHERE_AM_I
    write(*,*)" ending   petsc_driver"
#endif

end subroutine petsc_driver
