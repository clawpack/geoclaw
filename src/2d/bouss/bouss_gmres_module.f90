! =====================================================================
! bouss_gmres_module
!
! Self-contained, OpenMP-threaded linear solver to replace the PETSc KSP
! for the SGN Boussinesq system.  Everything in the Krylov loop is threaded:
!   - CSR matvec  (parallel over rows)
!   - vector ops  (dot/axpy/scale as OpenMP reductions/loops)
!   - preconditioner: the block Gauss-Seidel OpenMP Thomas apply
!     (apply_block_gs in bouss_tridiag_module)
!
! Right-preconditioned restarted GMRES(m).  Right preconditioning makes the
! Arnoldi residual equal the TRUE residual ||b - A x||, so the stopping test
! matches PETSc's default (||r_k|| <= rtol*||r_0||, with x_0 = 0 => ||r_0||=||b||).
!
! PETSc-free: depends only on OpenMP and bouss_tridiag_module, so it
! unit-tests with plain gfortran.
! =====================================================================
module bouss_gmres_module

    use bouss_tridiag_module, only: line_decomp_t, apply_block_gs
    implicit none
    private
    public :: csr_matvec, gmres_solve

    ! persistent work buffers for gmres_solve, grown as needed and reused across
    ! solves (avoids re-allocating the ~80 MB Krylov basis V every solve).
    ! gmres_solve is entered on the master thread only, so save is thread-safe.
    real(kind=8), allocatable, save :: V(:,:), H(:,:)
    real(kind=8), allocatable, save :: cs(:), sn(:), g(:), yy(:)
    real(kind=8), allocatable, save :: w(:), z(:), zsum(:), reshist(:)
    integer, save :: wk_ntot = -1, wk_m = -1

contains

    ! y = A x   (compressed CSR, 0-based), threaded over rows
    subroutine csr_matvec(rowPtr, cols, vals, nnz, ntot, x, y)
        integer,      intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1), x(0:ntot-1)
        real(kind=8), intent(out) :: y(0:ntot-1)
        integer :: r, p
        real(kind=8) :: s
        !$omp parallel do schedule(static) default(shared) private(r,p,s)
        do r = 0, ntot-1
            s = 0.d0
            do p = rowPtr(r), rowPtr(r+1)-1
                s = s + vals(p)*x(cols(p))
            end do
            y(r) = s
        end do
        !$omp end parallel do
    end subroutine csr_matvec

    real(kind=8) function omp_dot(n, x, y) result(d)
        integer,      intent(in) :: n
        real(kind=8), intent(in) :: x(0:n-1), y(0:n-1)
        integer :: i
        d = 0.d0
        !$omp parallel do reduction(+:d) default(shared) private(i)
        do i = 0, n-1
            d = d + x(i)*y(i)
        end do
        !$omp end parallel do
    end function omp_dot

    subroutine omp_axpy(n, a, x, y)      ! y = y + a*x
        integer,      intent(in)    :: n
        real(kind=8), intent(in)    :: a, x(0:n-1)
        real(kind=8), intent(inout) :: y(0:n-1)
        integer :: i
        !$omp parallel do default(shared) private(i)
        do i = 0, n-1
            y(i) = y(i) + a*x(i)
        end do
        !$omp end parallel do
    end subroutine omp_axpy

    subroutine omp_scalecopy(n, a, x, y)  ! y = a*x
        integer,      intent(in)  :: n
        real(kind=8), intent(in)  :: a, x(0:n-1)
        real(kind=8), intent(out) :: y(0:n-1)
        integer :: i
        !$omp parallel do default(shared) private(i)
        do i = 0, n-1
            y(i) = a*x(i)
        end do
        !$omp end parallel do
    end subroutine omp_scalecopy

    ! apply the preconditioner: z = M^{-1} r  (block G-S, or identity)
    subroutine apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, ntot, r, z)
        logical,             intent(in)  :: precondition
        type(line_decomp_t), intent(in)  :: D0, D1
        integer,             intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8),        intent(in)  :: vals(0:nnz-1), r(0:ntot-1)
        real(kind=8),        intent(out) :: z(0:ntot-1)
        if (precondition) then
            call apply_block_gs(D0, D1, rowPtr, cols, vals, nnz, ntot, r, z)
        else
            z = r
        endif
    end subroutine apply_pc

    ! -----------------------------------------------------------------
    ! Right-preconditioned restarted GMRES(m), zero initial guess.
    !   A given as compressed CSR (rowPtr,cols,vals); the preconditioner is
    !   block Gauss-Seidel using the line decompositions D0 (u), D1 (v)
    !   (set precondition=.false. for unpreconditioned GMRES).
    !   Converges when ||b - A x||_2 <= rtol * ||b||_2.
    !   iters = total inner iterations, resid = final relative residual,
    !   info  = 0 converged, 1 hit maxcycles*m without converging.
    ! -----------------------------------------------------------------
    subroutine gmres_solve(rowPtr, cols, vals, nnz, ntot, D0, D1, precondition, &
                           b, x, rtol, m, maxcycles, iters, resid, info)
        integer,      intent(in)  :: nnz, ntot, rowPtr(0:ntot), cols(0:nnz-1)
        real(kind=8), intent(in)  :: vals(0:nnz-1)
        type(line_decomp_t), intent(in) :: D0, D1
        logical,      intent(in)  :: precondition
        real(kind=8), intent(in)  :: b(0:ntot-1)
        real(kind=8), intent(out) :: x(0:ntot-1)
        real(kind=8), intent(in)  :: rtol
        integer,      intent(in)  :: m, maxcycles
        integer,      intent(out) :: iters, info
        real(kind=8), intent(out) :: resid

        real(kind=8) :: beta, bnorm, denom, tmp
        logical, save :: hist_printed = .false.
        integer :: i, j, k, cyc, jj

        if (ntot > wk_ntot .or. m > wk_m) then   ! grow persistent buffers
            if (allocated(V)) deallocate(V,H,cs,sn,g,yy,w,z,zsum,reshist)
            allocate(V(0:ntot-1, m+1), H(m+1, m))
            allocate(cs(m), sn(m), g(m+1), yy(m))
            allocate(w(0:ntot-1), z(0:ntot-1), zsum(0:ntot-1))
            allocate(reshist(0:m))
            wk_ntot = ntot; wk_m = m
        endif
        reshist = 0.d0

        x = 0.d0
        iters = 0
        info  = 1
        bnorm = sqrt(omp_dot(ntot, b, b))
        if (bnorm == 0.d0) then
            resid = 0.d0; info = 0
            return
        endif

        do cyc = 1, maxcycles
            ! true residual r = b - A x  (x0=0 first cycle)
            call csr_matvec(rowPtr, cols, vals, nnz, ntot, x, w)
            call omp_axpy(ntot, -1.d0, b, w)      ! w = A x - b
            beta = sqrt(omp_dot(ntot, w, w))      ! ||A x - b|| = ||r||
            resid = beta / bnorm
            reshist(0) = resid
            if (beta <= rtol*bnorm) then; info = 0; exit; endif
            call omp_scalecopy(ntot, -1.d0/beta, w, V(:,1))   ! v1 = (b - A x)/beta
            g = 0.d0; g(1) = beta

            k = 0
            do j = 1, m
                k = j
                call apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, &
                              ntot, V(:,j), z)                        ! z = M^-1 v_j
                call csr_matvec(rowPtr, cols, vals, nnz, ntot, z, w)  ! w = A z
                ! modified Gram-Schmidt
                do i = 1, j
                    H(i,j) = omp_dot(ntot, w, V(:,i))
                    call omp_axpy(ntot, -H(i,j), V(:,i), w)
                end do
                H(j+1,j) = sqrt(omp_dot(ntot, w, w))
                if (H(j+1,j) > 0.d0) then
                    call omp_scalecopy(ntot, 1.d0/H(j+1,j), w, V(:,j+1))
                else
                    V(:,j+1) = 0.d0
                endif
                ! apply previous Givens rotations to new column
                do i = 1, j-1
                    tmp      =  cs(i)*H(i,j) + sn(i)*H(i+1,j)
                    H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j)
                    H(i,j)   =  tmp
                end do
                ! new Givens rotation to eliminate H(j+1,j)
                denom = sqrt(H(j,j)**2 + H(j+1,j)**2)
                cs(j) = H(j,j)/denom
                sn(j) = H(j+1,j)/denom
                H(j,j)   = cs(j)*H(j,j) + sn(j)*H(j+1,j)
                H(j+1,j) = 0.d0
                tmp    =  cs(j)*g(j)
                g(j+1) = -sn(j)*g(j)
                g(j)   =  tmp
                resid  = abs(g(j+1)) / bnorm
                reshist(j) = resid
                iters  = iters + 1
                if (abs(g(j+1)) <= rtol*bnorm) exit
            end do

            ! solve upper-triangular H(1:k,1:k) yy = g(1:k)
            do i = k, 1, -1
                tmp = g(i)
                do j = i+1, k
                    tmp = tmp - H(i,j)*yy(j)
                end do
                yy(i) = tmp / H(i,i)
            end do
            ! x += M^{-1} (V(:,1:k) yy)   [right preconditioning]
            zsum = 0.d0
            do i = 1, k
                call omp_axpy(ntot, yy(i), V(:,i), zsum)
            end do
            call apply_pc(precondition, D0, D1, rowPtr, cols, vals, nnz, ntot, zsum, z)
            call omp_axpy(ntot, 1.d0, z, x)

            if (resid <= rtol) then
                info = 0; exit
            else
                ! only prints when a full cycle finished without converging
                write(*,'(a,i3,a,i5,a,es11.3)') " GMRES restart cycle ",cyc,   &
                      "  total its ",iters,"  rel resid ",resid
                if (.not. hist_printed) then
                    hist_printed = .true.
                    write(*,'(a,es12.4,a,i0)') "  [gmres stall] ||b|| = ",bnorm, &
                          "   ntot = ",ntot
                    write(*,'(a)') "  [gmres stall] intra-cycle relative residual:"
                    do jj = 0, k
                        if (jj <= 6 .or. mod(jj,10) == 0)                        &
                          write(*,'(a,i4,a,es12.4)') "      iter ",jj,":  ",reshist(jj)
                    end do
                endif
            endif
        end do
    end subroutine gmres_solve

end module bouss_gmres_module
