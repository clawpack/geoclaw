! Standalone check: parallel stream-compaction (new compressOut) must produce
! bit-identical cols/vals/rowPtr and numColsTot as the serial version.
!   gfortran -O2 -fopenmp test_compress.f90 -o test_compress && ./test_compress
program test_compress
    implicit none
    integer, parameter :: N = 4000, twoN = 2*N
    integer      :: colsin(0:24*N-1)
    real(kind=8) :: valsin(0:24*N-1)
    integer      :: cs(0:24*N-1), rps(0:twoN)      ! serial result
    real(kind=8) :: vs(0:24*N-1)
    integer      :: cp(0:24*N-1), rpp(0:twoN)      ! parallel result
    real(kind=8) :: vp(0:24*N-1)
    integer :: ntots, ntotp
    integer :: m, s, j, cnt, running, idest, isrc, icount12, k, seed
    integer, allocatable :: cols_tmp(:)
    real(kind=8), allocatable :: vals_tmp(:)
    logical :: ok

    ! synthetic input: each row m occupies slots [12m,12m+11]; slot 0 is a
    ! diagonal (always kept), others randomly a real entry or -1 padding.
    seed = 12345
    do m = 0, twoN-1
        do k = 0, 11
            s = 12*m + k
            seed = mod(seed*1103515245 + 12345, 2147483647)
            if (k == 0) then
                colsin(s) = m;  valsin(s) = 1.d0 + m
            else if (mod(seed,3) == 0) then
                colsin(s) = m*100 + k;  valsin(s) = real(s,8)
            else
                colsin(s) = -1;  valsin(s) = 0.d0
            endif
        end do
    end do

    ! ---- serial compaction (exact copy of original compressOut logic) ----
    cs = colsin;  vs = valsin
    idest = 0;  icount12 = 0
    do isrc = 0, 24*N-1
        if (icount12 == 0) rps(isrc/12) = idest
        if (cs(isrc) /= -1) then
            cs(idest) = cs(isrc);  vs(idest) = vs(isrc);  idest = idest + 1
        endif
        icount12 = icount12 + 1
        if (icount12 == 12) icount12 = 0
    end do
    ntots = idest

    ! ---- parallel compaction (new algorithm) ----
    cp = colsin;  vp = valsin
    ! single parallel region with barriers between phases (mirrors compressOut.f)
    allocate(cols_tmp(0:24*N-1), vals_tmp(0:24*N-1))
    !$omp parallel default(shared) private(m,s,j,cnt)
    !$omp do schedule(static)
    do m = 0, twoN-1
        cnt = 0
        do s = 12*m, 12*m+11
            if (cp(s) /= -1) cnt = cnt + 1
        end do
        rpp(m) = cnt
    end do
    !$omp end do
    !$omp single
    running = 0
    do m = 0, twoN-1
        cnt = rpp(m);  rpp(m) = running;  running = running + cnt
    end do
    ntotp = running
    !$omp end single
    !$omp do schedule(static)
    do m = 0, twoN-1
        j = rpp(m)
        do s = 12*m, 12*m+11
            if (cp(s) /= -1) then
                cols_tmp(j) = cp(s);  vals_tmp(j) = vp(s);  j = j + 1
            endif
        end do
    end do
    !$omp end do
    !$omp do schedule(static)
    do j = 0, ntotp-1
        cp(j) = cols_tmp(j);  vp(j) = vals_tmp(j)
    end do
    !$omp end do nowait
    !$omp end parallel
    deallocate(cols_tmp, vals_tmp)

    ! ---- compare ----
    ok = (ntots == ntotp)
    do m = 0, twoN-1
        if (rps(m) /= rpp(m)) ok = .false.
    end do
    do j = 0, ntots-1
        if (cs(j) /= cp(j))  ok = .false.
        if (vs(j) /= vp(j))  ok = .false.
    end do

    write(*,'(a,i0,a,i0)') 'ntot  serial=', ntots, '  parallel=', ntotp
    if (ok) then
        write(*,'(a)') '==== PASS: parallel compaction identical to serial ===='
    else
        write(*,'(a)') '==== FAIL ===='
        call exit(1)
    endif
end program test_compress
