
subroutine testSoln( n, soln, rhs, nelt, ia, ja, sa )

    ! matrix-vector multiply, returning y = A*soln to compare with rhs
    ! A is specified in triad form: sa(k) is element at (ia(k),ja(k))
    ! for k = 1:nelt

    implicit none

    integer, intent(in) :: n,nelt,ia(nelt),ja(nelt)
    real(kind=8), intent(in) :: soln(0:n),sa(nelt),rhs(0:n)

    real(kind=8)  :: y(n), phiu, phiv, phiu_max, phiv_max
    integer :: k, kumax, kvmax

    write(*,*)"Testing soln, n, nelt ",n,nelt
    y = 0.d0
    do k=1,nelt
        y(ia(k)) = y(ia(k)) + sa(k) * soln(ja(k))
        write(22,101) k,ia(k),ja(k),sa(k)
 101    format(3i6,e15.7)
    enddo

    phiu_max = 0.d0
    phiv_max = 0.d0
    kumax = -1
    kvmax = -1
    do k = 1, n/2
       phiu = abs(rhs(k)-y(k))
       phiv = abs(rhs(n/2+k)-y(n/2+k))
       if (phiu .gt. phiu_max) then
          phiu_max = phiu
          kumax = k
       endif
       if (phiv .gt. phiv_max) then
          phiv_max = phiv
          kvmax = k+n/2
       endif
    end do
    write(*,*) "phiu_max = ",phiu_max," at row ",kumax
    write(*,*) "phiv_max = ",phiv_max," at row ",kvmax

    return

end subroutine testSoln
