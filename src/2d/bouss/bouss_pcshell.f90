! =====================================================================
! bouss_pcshell  (item 4)
!
! PETSc PCSHELL glue for the directional block Gauss-Seidel preconditioner.
! Keeps the PETSc-dependent code out of bouss_tridiag_module (which stays
! standalone/unit-testable).  The actual math is apply_block_gs in that
! module; here we only:
!   - build/refresh the u,v line decompositions for a level (tridiag_pc_setup)
!   - unpack the PETSc residual vector, call apply_block_gs, pack the result
!     (tridiag_pc_apply, the PCShellSetApply callback)
!
! The outer solve is plain PETSc GMRES on PETSC_COMM_SELF (NO
! -mpi_linear_solver_server): GMRES does the A*x matvec on the assembled
! CRS matrix, this shell supplies M^-1.  Parallelism is OpenMP over the
! independent tridiagonal lines inside apply_block_gs.
!
! Decompositions are cached per level in module state; cur_level records
! which level the pending KSPSolve is for (solves are sequential per level,
! no nesting, so a module scalar is safe).
! =====================================================================
#ifdef HAVE_PETSC
module bouss_pcshell

#include <petsc/finclude/petscksp.h>
    use petscksp
    use bouss_tridiag_module, only: line_decomp_t, build_line_decomp,   &
                                    free_line_decomp, apply_block_gs
    use bouss_module, only: matrix_info_allLevs, matrix_levInfo
    use amr_module,   only: maxlv
    implicit none
    private
    public :: tridiag_pc_setup, tridiag_pc_apply

    type(line_decomp_t), save :: gD0(maxlv), gD1(maxlv)
    logical, save :: built(maxlv) = .false.
    integer, save :: cur_level = -1

contains

    ! Rebuild the u (field 0) and v (field 1) line decompositions for `level`
    ! from the currently assembled CRS matrix, and mark it as the current
    ! level for the pending solve.  Call immediately before KSPSolve; the
    ! matrix values change every step so we refresh each solve (O(nnz), cheap).
    subroutine tridiag_pc_setup(level)
        integer, intent(in) :: level
        type(matrix_levInfo), pointer :: minfo
        integer :: ntot, nnz

        minfo => matrix_info_allLevs(level)
        ntot = 2*minfo%numBoussCells
        nnz  = minfo%numColsTot

        if (built(level)) then
            call free_line_decomp(gD0(level))
            call free_line_decomp(gD1(level))
        endif
        call build_line_decomp(minfo%rowPtr, minfo%cols, minfo%vals,    &
                               nnz, ntot, 0, gD0(level))
        call build_line_decomp(minfo%rowPtr, minfo%cols, minfo%vals,    &
                               nnz, ntot, 1, gD1(level))
        built(level) = .true.
        cur_level = level
    end subroutine tridiag_pc_setup

    ! PCShell apply callback:  y = M^-1 x   (multiplicative block G-S).
    subroutine tridiag_pc_apply(pc, x, y, ierr)
        PC             pc
        Vec            x, y
        PetscErrorCode ierr
        PetscScalar, pointer :: xx(:), yy(:)
        type(matrix_levInfo), pointer :: minfo
        integer :: ntot, nnz, d
        real(kind=8), allocatable :: rin(:), yout(:)

        minfo => matrix_info_allLevs(cur_level)
        ntot = 2*minfo%numBoussCells
        nnz  = minfo%numColsTot

        call VecGetArrayRead(x, xx, ierr); if (ierr /= 0) return
        call VecGetArray(y, yy, ierr);     if (ierr /= 0) return

        allocate(rin(0:ntot-1), yout(0:ntot-1))
        do d = 0, ntot-1
            rin(d) = xx(d+1)          ! PETSc arrays are 1-based; dof d -> index d+1
        end do

        call apply_block_gs(gD0(cur_level), gD1(cur_level),             &
                            minfo%rowPtr, minfo%cols, minfo%vals,       &
                            nnz, ntot, rin, yout)

        do d = 0, ntot-1
            yy(d+1) = yout(d)
        end do
        deallocate(rin, yout)

        call VecRestoreArrayRead(x, xx, ierr); if (ierr /= 0) return
        call VecRestoreArray(y, yy, ierr);     if (ierr /= 0) return
        ierr = 0
    end subroutine tridiag_pc_apply

end module bouss_pcshell
#else
module bouss_pcshell
    implicit none
end module bouss_pcshell
#endif
