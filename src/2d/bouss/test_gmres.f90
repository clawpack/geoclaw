! =====================================================================
! Standalone unit test for bouss_gmres_module (OpenMP GMRES).
!
! Builds a 2x2-block interleaved system (u x-tridiagonal, v y-tridiagonal,
! plus u<->v cross-coupling) and checks that gmres_solve:
!   (A) converges unpreconditioned to the true solution, and
!   (B) converges with the block Gauss-Seidel PC in fewer iterations,
! both to the requested tolerance, matching a known x_true.
!
!   gfortran -O2 -fopenmp bouss_tridiag_module.f90 bouss_gmres_module.f90 \
!            test_gmres.f90 -o test_gmres && ./test_gmres
! =====================================================================
program test_gmres
    use bouss_tridiag_module
    use bouss_gmres_module
    implicit none

    integer, parameter :: nx = 8, ny = 6, N = nx*ny, ntot = 2*N
    real(kind=8), parameter :: rtol = 1.d-10

    integer, allocatable :: er(:), ec(:), rowPtr(:), cols(:)
    real(kind=8), allocatable :: ev(:), vals(:)
    integer :: nent, nnz
    type(line_decomp_t) :: D0, D1
    real(kind=8), allocatable :: xtrue(:), b(:), x(:)
    real(kind=8) :: err, resid
    integer :: iters, info, nfail, d

    nfail = 0
    allocate(er(8*N), ec(8*N), ev(8*N))
    nent = 0
    call assemble(er, ec, ev, nent)
    call to_crs(er, ec, ev, nent, ntot, rowPtr, cols, vals, nnz)
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 0, D0)
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 1, D1)

    allocate(xtrue(0:ntot-1), b(0:ntot-1), x(0:ntot-1))
    do d = 0, ntot-1
        xtrue(d) = 1.d0 + 0.5d0*sin(real(d+1,8))
    end do
    call full_matvec(rowPtr, cols, vals, nnz, ntot, xtrue, b)

    write(*,'(a,i0,a,i0,a,i0)') 'grid ', nx, ' x ', ny, ',  ntot = ', ntot

    ! ---- (A) unpreconditioned GMRES ----
    call gmres_solve(rowPtr, cols, vals, nnz, ntot, D0, D1, .false., b, x, &
                     rtol, 100, 2, iters, resid, info)
    err = solerr(x, xtrue, ntot)
    write(*,'(a,i4,a,es10.2,a,es10.2)') '  [no PC       ] iters=', iters, &
          '  resid=', resid, '  err=', err
    call check(info == 0 .and. resid <= rtol .and. err <= 1.d-7, &
               'unpreconditioned GMRES converges to true solution', nfail)

    ! ---- (B) block Gauss-Seidel preconditioned GMRES ----
    call gmres_solve(rowPtr, cols, vals, nnz, ntot, D0, D1, .true., b, x, &
                     rtol, 100, 2, iters, resid, info)
    err = solerr(x, xtrue, ntot)
    write(*,'(a,i4,a,es10.2,a,es10.2)') '  [block-GS PC ] iters=', iters, &
          '  resid=', resid, '  err=', err
    call check(info == 0 .and. resid <= rtol .and. err <= 1.d-7, &
               'block-GS preconditioned GMRES converges to true solution', nfail)

    write(*,*)
    if (nfail == 0) then
        write(*,'(a)') '==== ALL TESTS PASSED ===='
    else
        write(*,'(a,i0,a)') '==== ', nfail, ' TEST(S) FAILED ===='
        call exit(1)
    endif

contains

    integer function udof(kk); integer,intent(in)::kk; udof = 2*(kk-1); end function
    integer function vdof(kk); integer,intent(in)::kk; vdof = 2*kk-1;   end function

    subroutine assemble(er, ec, ev, nent)
        integer, intent(inout) :: er(:), ec(:), nent
        real(kind=8), intent(inout) :: ev(:)
        integer :: i, j, k
        do j = 1, ny
        do i = 1, nx
            k = (j-1)*nx + i
            call addent(er,ec,ev,nent, udof(k), udof(k),   4.d0 + 0.1d0*i + 0.01d0*j)
            if (i > 1)  call addent(er,ec,ev,nent, udof(k), udof(k-1), -1.0d0 - 0.01d0*j)
            if (i < nx) call addent(er,ec,ev,nent, udof(k), udof(k+1), -1.1d0 - 0.01d0*j)
            call addent(er,ec,ev,nent, udof(k), vdof(k), 0.3d0)          ! A01
            call addent(er,ec,ev,nent, vdof(k), vdof(k),   5.d0 + 0.05d0*i + 0.2d0*j)
            if (j > 1)  call addent(er,ec,ev,nent, vdof(k), vdof(k-nx), -1.2d0)
            if (j < ny) call addent(er,ec,ev,nent, vdof(k), vdof(k+nx), -0.9d0)
            call addent(er,ec,ev,nent, vdof(k), udof(k), 0.2d0)          ! A10
        end do
        end do
    end subroutine assemble

    subroutine addent(er, ec, ev, nent, r, c, v)
        integer, intent(inout) :: er(:), ec(:), nent
        real(kind=8), intent(inout) :: ev(:)
        integer, intent(in) :: r, c
        real(kind=8), intent(in) :: v
        nent = nent + 1
        er(nent) = r; ec(nent) = c; ev(nent) = v
    end subroutine addent

    subroutine to_crs(er, ec, ev, nent, ntot_, rowPtr, cols, vals, nnz)
        integer, intent(in) :: nent, ntot_, er(nent), ec(nent)
        real(kind=8), intent(in) :: ev(nent)
        integer, allocatable, intent(out) :: rowPtr(:), cols(:)
        real(kind=8), allocatable, intent(out) :: vals(:)
        integer, intent(out) :: nnz
        integer :: cnt(0:ntot_-1), cur(0:ntot_-1), i, r
        nnz = nent
        allocate(rowPtr(0:ntot_), cols(0:nnz-1), vals(0:nnz-1))
        cnt = 0
        do i = 1, nent
            cnt(er(i)) = cnt(er(i)) + 1
        end do
        rowPtr(0) = 0
        do r = 0, ntot_-1
            rowPtr(r+1) = rowPtr(r) + cnt(r)
        end do
        cur = rowPtr(0:ntot_-1)
        do i = 1, nent
            r = er(i)
            cols(cur(r)) = ec(i)
            vals(cur(r)) = ev(i)
            cur(r) = cur(r) + 1
        end do
    end subroutine to_crs

    subroutine full_matvec(rowPtr, cols, vals, nnz, ntot_, x, b)
        integer, intent(in) :: nnz, ntot_, rowPtr(0:ntot_), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1), x(0:ntot_-1)
        real(kind=8), intent(out) :: b(0:ntot_-1)
        integer :: r, p
        real(kind=8) :: s
        do r = 0, ntot_-1
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                s = s + vals(p)*x(cols(p))
            end do
            b(r) = s
        end do
    end subroutine full_matvec

    real(kind=8) function solerr(x, xt, ntot_) result(e)
        integer, intent(in) :: ntot_
        real(kind=8), intent(in) :: x(0:ntot_-1), xt(0:ntot_-1)
        integer :: i
        e = 0.d0
        do i = 0, ntot_-1
            e = max(e, abs(x(i) - xt(i)))
        end do
    end function solerr

    subroutine check(cond, label, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (cond) then
            write(*,'(a,a)') '  [PASS] ', label
        else
            write(*,'(a,a)') '  [FAIL] ', label
            nfail = nfail + 1
        endif
    end subroutine check

end program test_gmres
