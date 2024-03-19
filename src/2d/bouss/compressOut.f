c
c -----------------------------------------------------
c
       subroutine compressOut(vals,rowPtr,cols,numBoussCells,numColsTot)

       use bouss_module

       implicit none
       integer numBoussCells, numColsTot
       integer rowPtr(0:2*numBoussCells), cols(0:24*numBoussCells)
       real*8  vals(0:24*numBoussCells)
       integer isrc,idest,icount12

c  compress out the -1 indices in cols, which indicates no entries
c  vals follow cols. Adjust rowPtr as needed along the way. 

        idest = 0
        icount12 = 0
        do isrc = 0, 24*numBoussCells-1
           if (icount12 .eq. 0) rowPtr(isrc/12) = idest
           if (cols(isrc) .ne. -1) then
              cols(idest) = cols(isrc)
              vals(idest) = vals(isrc)
              idest = idest + 1
           endif
           icount12 = icount12 + 1
           if (icount12 .eq. 12) icount12 = 0
        end do

        numColsTot = idest
        ! last entry for rowPtr set on return

        return
        end

