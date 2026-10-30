c
c -----------------------------------------------------
c
       subroutine compressOut(vals,rowPtr,cols,numBoussCells,numColsTot)

       use bouss_module

       implicit none
       integer numBoussCells, numColsTot
       integer rowPtr(0:2*numBoussCells), cols(0:24*numBoussCells)
       real*8  vals(0:24*numBoussCells)
       integer m, s, j, cnt, running, twoN
c      persistent scratch, reused across solves (avoid re-allocating
c      ~10s of MB every call); grown only when a larger level appears
       integer, allocatable, save :: cols_tmp(:)
       real*8,  allocatable, save :: vals_tmp(:)

c  Compress out the -1 entries in cols (which mark "no entry"); vals
c  follow cols.  Each matrix row occupies 12 uncompressed slots,
c  [12*m, 12*m+11].  Parallel stream compaction, done in a SINGLE parallel
c  region (one fork/join) with barriers between phases instead of three
c  separate PARALLEL DO regions (three fork/joins):
c    (1) count non(-1) entries per row            (!$OMP DO, barrier)
c    (2) prefix-sum counts -> row start offsets,  (!$OMP SINGLE, serial)
c        and grow the scratch buffer               + barrier
c    (3) compact each row into SCRATCH at offset   (!$OMP DO, barrier)
c        (disjoint dests -> race free; separate array -> no in-place hazard)
c    (4) copy compacted data back into cols,vals   (!$OMP DO)
c  Produces exactly the same cols/vals/rowPtr as the serial compaction.

       twoN = 2*numBoussCells

!$OMP PARALLEL DEFAULT(shared) PRIVATE(m,s,j,cnt)

c  (1) count non(-1) entries in each row's 12 slots
!$OMP DO SCHEDULE(static)
       do m = 0, twoN-1
          cnt = 0
          do s = 12*m, 12*m+11
             if (cols(s) .ne. -1) cnt = cnt + 1
          end do
          rowPtr(m) = cnt
       end do
!$OMP END DO
c     (implicit barrier: all counts written before the prefix sum)

c  (2) prefix sum -> starting offsets, and grow scratch.  Serial, so one
c      thread does it; SINGLE's implicit barrier publishes rowPtr/numColsTot
c      and the (re)allocated scratch to all threads before phase (3).
!$OMP SINGLE
       running = 0
       do m = 0, twoN-1
          cnt = rowPtr(m)
          rowPtr(m) = running
          running = running + cnt
       end do
       numColsTot = running
       if (.not. allocated(cols_tmp)) then
          allocate(cols_tmp(0:max(numColsTot,1)-1))
          allocate(vals_tmp(0:max(numColsTot,1)-1))
       else if (size(cols_tmp) .lt. numColsTot) then
          deallocate(cols_tmp, vals_tmp)
          allocate(cols_tmp(0:numColsTot-1))
          allocate(vals_tmp(0:numColsTot-1))
       endif
!$OMP END SINGLE
c     (implicit barrier)

c  (3) compact each row into scratch at its offset
!$OMP DO SCHEDULE(static)
       do m = 0, twoN-1
          j = rowPtr(m)
          do s = 12*m, 12*m+11
             if (cols(s) .ne. -1) then
                cols_tmp(j) = cols(s)
                vals_tmp(j) = vals(s)
                j = j + 1
             endif
          end do
       end do
!$OMP END DO
c     (implicit barrier: scratch fully written before copy-back)

c  (4) copy compacted data back (disjoint indices, race free)
!$OMP DO SCHEDULE(static)
       do j = 0, numColsTot-1
          cols(j) = cols_tmp(j)
          vals(j) = vals_tmp(j)
       end do
!$OMP END DO NOWAIT
c     (nowait: nothing after this in the region; END PARALLEL barrier follows)

!$OMP END PARALLEL

c  scratch kept allocated (persistent) for reuse next call

c  last entry for rowPtr (= numColsTot) is set by the caller on return
       return
       end
