
module bouss_module

    use amr_module, only: t0, maxlv

#ifdef HAVE_PETSC
#include <petsc/finclude/petscksp.h>
    use petscksp
#endif

    implicit none
    save


    integer :: bc_xlo, bc_xhi, bc_ylo, bc_yhi

    real(kind=8) :: Bparam
    real(kind=8) :: alpha

    logical :: startWithBouss

    logical :: crs, origCooFormat   

    real(kind=8) :: startBoussTime

    !
!   Added for Boussinesq hybrid solver
! ================================
!
    type matrix_patchIndex
       integer, allocatable, dimension(:,:) :: mindex
       logical*1, allocatable, dimension(:,:) :: isBouss
       integer, allocatable, dimension(:,:) :: mindexSolo
       !!! TEST FOR 2 LEVEL CODE, REALLY WANT LISTS
       integer, allocatable, dimension(:,:) :: mfine
       integer, allocatable, dimension(:,:) :: ifine
       integer, allocatable, dimension(:,:) :: jfine
       integer, allocatable, dimension(:,:) :: mcoarse
       integer, allocatable, dimension(:,:) :: icoarse
       integer, allocatable, dimension(:,:) :: jcoarse
    end type

    type matrix_levInfo
    !! matrix_indices is indexed into by patch number, only Bouss Grids at that level
       type (matrix_patchIndex), allocatable, dimension(:) :: matrix_indices
       integer, allocatable, dimension(:) :: matrix_ia, matrix_ja
       real(kind=8), allocatable, dimension(:) :: matrix_sa
       integer :: numBoussGrids, numBoussCells, numBoussCellsSolo, numGhostCount, numUnset
       integer :: numColsTot

       ! for CRS format
       integer, allocatable, dimension(:) :: rowPtr, cols
       real(kind=8), allocatable, dimension(:) :: vals 

       ! intfCountc are the # equations (cells)  added for you as a coarse grid
       ! these are the equations under a fine grid on the first interior cells
       ! that say the coarse cell is the conservative average of the refined cells

       ! intfCountf are the # equations (cells)  if you are the fine grid
       ! these are the equations to interpolate a ghost cell from the
       ! coarser patch
       integer :: intfCount
       integer :: matrix_nelt
    end type

    double precision :: boussMinDepth
    integer :: minLevelBouss, maxLevelBouss
    integer :: isolver
    integer :: ibouss

    type (matrix_levInfo), target :: matrix_info_allLevs(maxlv)

#ifdef HAVE_PARDISO
    ! vars for pardiso solver
    integer :: mtype,solver,maxfct,mnum,phase,nn,nrhs,error,msglvl
    integer iparm(64), idum
    real*8 dparm(64),ddum
    integer*8 pt(64,maxlv)
#endif

    ! if new grids have to redo factorization
    ! either first time in setgrd of in regrid
    logical newGrids(maxlv), newGridsTower(maxlv)
    integer*8  itcount(maxlv), numTimes(maxlv)

#ifdef HAVE_PETSC
    Mat Jr(maxlv)
    KSP ksp(maxlv)
    Mat JrTower(maxlv)
    KSP kspTower(maxlv)
#endif

    ! for restarts:
    integer :: minLevelBouss_orig,maxLevelBouss_orig
    real(kind=8) :: boussMinDepth_orig,startBoussTime_orig

contains

    subroutine set_bouss(rest,time,naux)

    ! Set Bparam and bc choices for Boussinesq implicit solver

    ! Eventually allow more things to be specified in setrun.py
    ! and read them in here.

    use amr_module, only: mthbc, outunit
    implicit none
    logical, intent(in) :: rest
    real(kind=8), intent(in) :: time
    integer, intent(in) :: naux
    
    integer iunit,i
    character(len=25) fname
    real(kind=8) :: eps

#ifdef WHERE_AM_I
    write(*,*) 'starting set_bouss'
#endif

    iunit = 7
    fname = 'bouss.data'
!   # open the unit with new routine from Clawpack 4.4 to skip over
!   # comment lines starting with #:
    call opendatafile(iunit, fname)

    ! set constants for whichever bouss solver being used
    ! MS equation parameter B:
    Bparam = 1.d0/15.d0  

    ! SGN equation parameter alpha
    !alpha = 1.d0
    alpha = 1.153d0


    ! modify write statements to say what value of Bparam or alpha is used:

    read(7,*) ibouss
    if (ibouss .eq. 0) then
       write(*,*)" Using Shallow Water equations"
       write(outunit,*)" Using Shallow Water equations"
    else if (ibouss .eq. 1) then
       write(*,*)" Using Madsen Sorensen equations"
       write(outunit,*)" Using Madsen Sorensen equations"
    else if (ibouss .eq. 2) then
       write(*,*)" Using SGN equations"
       write(outunit,*)" Using SGN equations"
    else
       write(*,*)" Unrecognized option for equation set"
       write(outunit,*)" Unrecognized option for equation set"
       stop
    endif
    
    read(7,*) minLevelBouss
    read(7,*) maxLevelBouss
    read(7,*) boussMinDepth
    read(7,*) isolver
    read(7,*) startBoussTime

    if (ibouss == 1) then
        ! Madsen-Sorensen equations only implemented with COO triplet storage
        ! format (i,j,value) for sparse matrix:
        crs = .false.
    else
        ! SGN implemented for both sparse storage formats, but runs faster
        ! using compressed row storage (CRS, also known as CSR or Yale) format:
        crs = .true.
    endif
    ! if crs is false this will force SGN to use triplet COO form
    ! crs = .false.  

     ! origCoo format ordered unknowns by all u updates then all v updates 
     ! block format is u and v together for a given cell. Better cache
     ! performance, better for debugging and comparing with crs.
     ! also reads into matlab for cond number testing, plotting,
     ! so keeping around
     !origCooFormat = .false.  ! set to false to use block format with coo
     ! setting to avoid uninitialized var, but not used unless crs set to false
     origCooFormat = .true. 

    !------------------------------------------
    if (rest) then
        ! for restart, parameters should not have changed from original run:
        if ((minLevelBouss .ne. minLevelBouss_orig) .or. &
            (maxLevelBouss .ne. maxLevelBouss_orig) .or. &
            (boussMinDepth .ne. boussMinDepth_orig)) then
          write(*,*)"*** Cannot change these bouss params on restart:"
          write(*,*)"*** Old params:"
          write(*,922) minLevelBouss_orig,maxLevelBouss_orig,boussMinDepth_orig
          write(*,*)"*** New params:"
          write(*,922) minLevelBouss,maxLevelBouss,boussMinDepth
  922     format("min/max Bousslevel "2i5," boussMinDepth ",f10.2)
          stop
        endif

        ! check if later startBoussTime than orig, only allowed if not yet using Bouss
        eps = 1.d-13
        ! allowing for roundoff between different compilers in testing
        ! if starting time has  changed
        if (abs(startBoussTime-startBoussTime_orig) .gt. eps)  then
            if (startBoussTime .ge. startBoussTime_orig) then
                if (startBoussTime_orig .lt. time) then
                    write(*,904)
                    write(outunit,904)
  904               format("Cannot postpone Bouss solves if already using it")
                    stop
                else
                    write(*,905) startBoussTime
                    write(outunit,905) startBoussTime
  905               format("Changing Bouss solves to start at time ", e15.7)
                endif
            else 
                write(*,906) time,startBoussTime_orig, startBoussTime
                write(outunit,906) time,startBoussTime_orig,startBoussTime
  906           format(" Already using Bouss at time",d15.7,/, &
                   " new start time ",d15.7," orig start time ",d15.7)
            endif
        endif
        ! allocate the data structures for the bouss info
        call resetBoussStuff(naux)

    endif  ! end checking parameters on restart
    !------------------------------------------

 99   write(*,900) minLevelBouss, maxLevelBouss
      write(outunit,900) minLevelBouss, maxLevelBouss
 900  format("==> Applying Bouss equations to selected grids between levels ",i3," and ",i3)

      write(*,*)"==> Use Bouss. in water deeper than ",boussMinDepth
      write(outunit,*)"==> Use Bouss. in water deeper than ",boussMinDepth

      if (isolver .eq. 1) then
         write(*,*)" No longer supporting GMRES solver"
         stop
      else if (isolver .eq. 2) then
#ifdef HAVE_PARDISO
         !write(*,*)" Using Pardiso solver"
         !write(outunit,*)" Using Pardiso solver"
         write(*,*)"Cannot use expired Pardiso solver"
         write(outunit,*)"Cannot use expired Pardiso solver"
         !stop
#else
         write(*,*)"need to install Pardiso for this option"
         stop
#endif
        else if (isolver .eq. 3) then
#ifdef HAVE_PETSC
         write(*,*)"Using a PETSc solver"
         write(outunit,*)"Using PETSc solver"
#else
         write(*,*)"need to install PETSc for this option"
         stop
#endif
      else
         write(*,*)"Unknown solver",isolver," choose 1,2 or 3"
         write(outunit,*)"Unknown solver",isolver," choose 1,2 or 3"
         stop
      endif

      if (startBoussTime .le. t0) then
         write(*,*)"Using Bouss equations from the start"
         write(outunit,*)"Using Bouss equations from the start"
         startWithBouss = .true.
      else
         write(*,*)"==> Wait until time ",startBoussTime," for before starting Bouss"
         write(*,*)"==> Using SWE until then."
         startWithBouss = .false.
      endif

    close(unit=iunit)



    ! Boundary conditions to impose in computing Boussinesq update:
    if ((mthbc(1)==2) .or. (mthbc(2)==2) &
        .or. (mthbc(3)==2) .or. (mthbc(4)==2)) then
        write(6,*) '*** Periodic BCs newly added to bouss_module'
        !!stop
    endif

    if ((mthbc(1)==4) .or. (mthbc(2)==4) &
        .or. (mthbc(3)==4) .or. (mthbc(4)==4)) then
        write(6,*) '*** Sphere BCs not supported in bouss_module'
        stop
    endif
    
    ! Dirichlet BCs using ghost cell values saved from coarser level:
    ! Requires num_eqn == 5.
    ! This is used in buildSparseMatrix at patch boundaries that are not
    ! domain boundaries, so the values of bc_xlo etc. only affect domain bdries
    ! Default (only if mthbc=0?) unless set otherwise by mthbc tests below

    bc_xlo = 2
    bc_xhi = 2
    bc_ylo = 2
    bc_yhi = 2
    

    ! For extrapolation in x:
    if (mthbc(1)==1) bc_xlo = 1
    if (mthbc(2)==1) bc_xhi = 1

    ! For extrapolation in y:
    if (mthbc(3)==1) bc_ylo = 1
    if (mthbc(4)==1) bc_yhi = 1


    ! For solid wall in x:
    if (mthbc(1)==3) bc_xlo = 3
    if (mthbc(2)==3) bc_xhi = 3

    ! For solid wall in y:
    if (mthbc(3)==3) bc_ylo = 3
    if (mthbc(4)==3) bc_yhi = 3
    
    
    ! FOR A PLANE WAVE CASE, SOLID WALL SHOULD WORK - Don't need what's below
    
    ! For plane wave in x-direction, should instead use Neumann conditions in y:
    ! Can use wall BC in SWE step so assume 
    !  clawdata.bc_lower[1] = clawdata.bc_upper[1] = 3
    !if (mthbc(3)==3) bc_ylo = 1
    !if (mthbc(4)==3) bc_yhi = 1
    !write(6,*) 'For plane wave in x-direction, using Neumann at top,bottom'


    ! For plane wave in y-direction, should instead use Neumann conditions in x:
    !if (mthbc(1)==3) bc_xlo = 1
    !if (mthbc(2)==3) bc_xhi = 1
    !write(6,*) 'For plane wave in y-direction, using Neumann at left,right'

#ifdef WHERE_AM_I
    write(*,*) 'ending   set_bouss'
#endif

    end subroutine set_bouss

end module bouss_module
