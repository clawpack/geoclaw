! =====================================================================
! Standalone unit test for bouss_tridiag_module (items 1-2).
!
! Builds a small interleaved 2x2-block system on an nx x ny grid:
!   - field 0 (u, even rows): tridiagonal along x-rows (i-neighbours)
!   - field 1 (v, odd  rows): tridiagonal along y-columns (j-neighbours)
!   - plus u<->v cross-coupling entries that extraction MUST ignore
! matching the real SGN matrix structure (compressed CRS, 0-based).
!
! Checks:
!   1. path decomposition covers every field dof exactly once
!   2. #paths(u) = ny (one per row), #paths(v) = nx (one per column)
!   3. Thomas exactly inverts the block that a generic same-parity
!      matvec applies:  solve(A_block * x_true) == x_true
!   4. revert cell: zero a cell's u off-diagonals (leaving its neighbours'
!      coupling to it) -> asymmetric tridiagonal, still solved exactly,
!      line still intact via symmetric-closure adjacency.
!
!   build & run:
!     gfortran -O2 -fopenmp bouss_tridiag_module.f90 test_tridiag.f90 \
!              -o test_tridiag && ./test_tridiag
! =====================================================================
program test_tridiag
    use bouss_tridiag_module
    implicit none

    integer, parameter :: nx = 4, ny = 3, N = nx*ny, ntot = 2*N
    real(kind=8), parameter :: tol = 1.d-10

    integer, allocatable :: er(:), ec(:)
    real(kind=8), allocatable :: ev(:)
    integer, allocatable :: rowPtr(:), cols(:)
    real(kind=8), allocatable :: vals(:)
    integer :: nent, nnz, k22, i
    type(line_decomp_t) :: D0, D1, DR
    real(kind=8), allocatable :: xtrue(:), b(:), x(:)
    real(kind=8) :: err
    integer :: nfail

    nfail = 0
    allocate(er(8*N), ec(8*N), ev(8*N))
    nent = 0
    call assemble(er, ec, ev, nent)
    call to_crs(er, ec, ev, nent, ntot, rowPtr, cols, vals, nnz)
    write(*,'(a,i0,a,i0,a,i0)') 'grid ', nx, ' x ', ny, ',  nnz = ', nnz

    allocate(xtrue(0:ntot-1), b(0:ntot-1), x(0:ntot-1))
    call set_xtrue(xtrue)

    ! ---- field 0 (u) ----
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 0, D0)
    call check_coverage(D0, 'u', nfail)
    call check_int(D0%npaths, ny, '#paths(u) == ny', nfail)
    call check_true(.not. D0%has_cycle, 'u has no cycle', nfail)
    call spmv_field(rowPtr, cols, vals, nnz, ntot, 0, xtrue, b)
    call apply_line_solve(D0, b, x)
    err = fielderr(x, xtrue, ntot, 0)
    call check_err(err, tol, 'u block solve exact', nfail)

    ! ---- field 1 (v) ----
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 1, D1)
    call check_coverage(D1, 'v', nfail)
    call check_int(D1%npaths, nx, '#paths(v) == nx', nfail)
    call check_true(.not. D1%has_cycle, 'v has no cycle', nfail)
    call spmv_field(rowPtr, cols, vals, nnz, ntot, 1, xtrue, b)
    call apply_line_solve(D1, b, x)
    err = fielderr(x, xtrue, ntot, 1)
    call check_err(err, tol, 'v block solve exact', nfail)

    ! ---- revert cell (2,2): zero its u off-diagonals, keep neighbours' ----
    k22 = (2-1)*nx + 2
    call zero_entry(rowPtr, cols, vals, nnz, ntot, udof(k22), udof(k22-1))
    call zero_entry(rowPtr, cols, vals, nnz, ntot, udof(k22), udof(k22+1))
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 0, DR)
    call check_coverage(DR, 'u-revert', nfail)
    call check_int(DR%npaths, ny, '#paths(u-revert) == ny (line intact)', nfail)
    call spmv_field(rowPtr, cols, vals, nnz, ntot, 0, xtrue, b)
    call apply_line_solve(DR, b, x)
    err = fielderr(x, xtrue, ntot, 0)
    call check_err(err, tol, 'u-revert block solve exact (asymmetric)', nfail)

    ! ---- apply_block_gs: exact inverse on a block-LOWER-triangular system ----
    ! With A01=0 the matrix is [[A00,0],[A10,A11]], for which forward block
    ! Gauss-Seidel IS the exact inverse, so block_gs(A*x_true) must recover
    ! x_true.  Exercises the A10 cross-term matvec, sign, and u/v wiring.
    call free_line_decomp(D0); call free_line_decomp(D1); call free_line_decomp(DR)
    nent = 0
    call assemble(er, ec, ev, nent, lower_tri=.true.)
    call to_crs(er, ec, ev, nent, ntot, rowPtr, cols, vals, nnz)
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 0, D0)
    call build_line_decomp(rowPtr, cols, vals, nnz, ntot, 1, D1)
    call set_xtrue(xtrue)
    call full_matvec(rowPtr, cols, vals, nnz, ntot, xtrue, b)
    call apply_block_gs(D0, D1, rowPtr, cols, vals, nnz, ntot, b, x)
    err = 0.d0
    do i = 0, ntot-1
        err = max(err, abs(x(i) - xtrue(i)))
    end do
    call check_err(err, tol, 'block-GS exact on lower-triangular system', nfail)

    ! smoke-test the in-run validation routine on the synthetic matrix
    write(*,*)
    write(*,'(a)') 'validate_line_solve() smoke test:'
    call validate_line_solve(rowPtr, cols, vals, nnz, ntot, 1)

    write(*,*)
    if (nfail == 0) then
        write(*,'(a)') '==== ALL TESTS PASSED ===='
    else
        write(*,'(a,i0,a)') '==== ', nfail, ' TEST(S) FAILED ===='
        call exit(1)
    endif

contains

    integer function udof(kk);  integer,intent(in)::kk; udof = 2*(kk-1); end function
    integer function vdof(kk);  integer,intent(in)::kk; vdof = 2*kk-1;   end function

    subroutine assemble(er, ec, ev, nent, lower_tri)
        integer, intent(inout) :: er(:), ec(:), nent
        real(kind=8), intent(inout) :: ev(:)
        logical, intent(in), optional :: lower_tri
        integer :: i, j, k
        logical :: lt
        lt = .false.
        if (present(lower_tri)) lt = lower_tri
        do j = 1, ny
        do i = 1, nx
            k = (j-1)*nx + i
            ! u-row: x-tridiagonal (i-neighbours) + one cross term to v (A01)
            call addent(er,ec,ev,nent, udof(k), udof(k),   4.d0 + 0.1d0*i + 0.01d0*j)
            if (i > 1)  call addent(er,ec,ev,nent, udof(k), udof(k-1), -1.0d0 - 0.01d0*j)
            if (i < nx) call addent(er,ec,ev,nent, udof(k), udof(k+1), -1.1d0 - 0.01d0*j)
            if (.not. lt) call addent(er,ec,ev,nent, udof(k), vdof(k), 0.3d0)  ! A01 (omit -> lower tri)
            ! v-row: y-tridiagonal (j-neighbours) + one cross term to u (A10)
            call addent(er,ec,ev,nent, vdof(k), vdof(k),   5.d0 + 0.05d0*i + 0.2d0*j)
            if (j > 1)  call addent(er,ec,ev,nent, vdof(k), vdof(k-nx), -1.2d0)
            if (j < ny) call addent(er,ec,ev,nent, vdof(k), vdof(k+nx), -0.9d0)
            call addent(er,ec,ev,nent, vdof(k), udof(k), 0.2d0)          ! A10 (lower block)
        end do
        end do
    end subroutine assemble

    ! full sparse matvec (all columns): b = A x
    subroutine full_matvec(rowPtr, cols, vals, nnz, ntot, x, b)
        integer, intent(in) :: nnz, ntot
        integer, intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1), x(0:ntot-1)
        real(kind=8), intent(out) :: b(0:ntot-1)
        integer :: r, p
        real(kind=8) :: s
        do r = 0, ntot-1
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                s = s + vals(p)*x(cols(p))
            end do
            b(r) = s
        end do
    end subroutine full_matvec

    subroutine addent(er, ec, ev, nent, r, c, v)
        integer, intent(inout) :: er(:), ec(:), nent
        real(kind=8), intent(inout) :: ev(:)
        integer, intent(in) :: r, c
        real(kind=8), intent(in) :: v
        nent = nent + 1
        er(nent) = r; ec(nent) = c; ev(nent) = v
    end subroutine addent

    subroutine to_crs(er, ec, ev, nent, ntot, rowPtr, cols, vals, nnz)
        integer, intent(in) :: nent, ntot, er(nent), ec(nent)
        real(kind=8), intent(in) :: ev(nent)
        integer, allocatable, intent(out) :: rowPtr(:), cols(:)
        real(kind=8), allocatable, intent(out) :: vals(:)
        integer, intent(out) :: nnz
        integer :: cnt(0:ntot-1), cur(0:ntot-1), i, r
        nnz = nent
        allocate(rowPtr(0:ntot), cols(0:nnz-1), vals(0:nnz-1))
        cnt = 0
        do i = 1, nent
            cnt(er(i)) = cnt(er(i)) + 1
        end do
        rowPtr(0) = 0
        do r = 0, ntot-1
            rowPtr(r+1) = rowPtr(r) + cnt(r)
        end do
        cur = rowPtr(0:ntot-1)
        do i = 1, nent
            r = er(i)
            cols(cur(r)) = ec(i)
            vals(cur(r)) = ev(i)
            cur(r) = cur(r) + 1
        end do
    end subroutine to_crs

    ! generic same-parity (block) matvec: b = A_block * x
    subroutine spmv_field(rowPtr, cols, vals, nnz, ntot, parity, x, b)
        integer, intent(in) :: nnz, ntot, parity
        integer, intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in) :: vals(0:nnz-1), x(0:ntot-1)
        real(kind=8), intent(out) :: b(0:ntot-1)
        integer :: r, p, c
        real(kind=8) :: s
        b = 0.d0
        do r = parity, ntot-1, 2
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                c = cols(p)
                if (mod(c,2) == parity) s = s + vals(p)*x(c)
            end do
            b(r) = s
        end do
    end subroutine spmv_field

    subroutine zero_entry(rowPtr, cols, vals, nnz, ntot, r, c)
        integer, intent(in) :: nnz, ntot, r, c
        integer, intent(in) :: rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(inout) :: vals(0:nnz-1)
        integer :: p
        do p = rowPtr(r), rowPtr(r+1)-1
            if (cols(p) == c) vals(p) = 0.d0
        end do
    end subroutine zero_entry

    subroutine set_xtrue(xtrue)
        real(kind=8), intent(out) :: xtrue(0:ntot-1)
        integer :: d
        do d = 0, ntot-1
            xtrue(d) = 1.d0 + 0.5d0*sin(real(d+1,8))
        end do
    end subroutine set_xtrue

    real(kind=8) function fielderr(x, xt, ntot, parity) result(e)
        integer, intent(in) :: ntot, parity
        real(kind=8), intent(in) :: x(0:ntot-1), xt(0:ntot-1)
        integer :: r
        e = 0.d0
        do r = parity, ntot-1, 2
            e = max(e, abs(x(r) - xt(r)))
        end do
    end function fielderr

    subroutine check_coverage(D, label, nfail)
        type(line_decomp_t), intent(in) :: D
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        integer, allocatable :: seen(:)
        integer :: slot, r, cnt
        logical :: ok
        allocate(seen(0:D%ntot-1)); seen = 0
        cnt = 0
        do slot = 1, D%nfield
            r = D%dof(slot)
            seen(r) = seen(r) + 1
            cnt = cnt + 1
        end do
        ok = (cnt == D%nfield) .and. (D%pathStart(D%npaths+1)-1 == D%nfield)
        do r = D%parity, D%ntot-1, 2
            if (seen(r) /= 1) ok = .false.
        end do
        call report('coverage ('//trim(label)//')', ok, nfail)
        deallocate(seen)
    end subroutine check_coverage

    subroutine check_int(got, want, label, nfail)
        integer, intent(in) :: got, want
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        call report(label, got == want, nfail)
        if (got /= want) write(*,'(a,i0,a,i0)') '      got ', got, ' want ', want
    end subroutine check_int

    subroutine check_err(err, tol, label, nfail)
        real(kind=8), intent(in) :: err, tol
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        logical :: ok
        ok = (err <= tol)
        call report(label, ok, nfail)
        write(*,'(a,es10.2)') '      max err = ', err
    end subroutine check_err

    subroutine check_true(cond, label, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail
        call report(label, cond, nfail)
    end subroutine check_true

    subroutine report(label, ok, nfail)
        character(*), intent(in) :: label
        logical, intent(in) :: ok
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,'(a,a)') '  [PASS] ', label
        else
            write(*,'(a,a)') '  [FAIL] ', label
            nfail = nfail + 1
        endif
    end subroutine report

end program test_tridiag
