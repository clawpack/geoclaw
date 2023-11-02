
subroutine matvec( n, x, y, nelt, ia, ja, a, isym )

    ! matrix-vector multiply, returning y = A*x
    ! A is specified in triad form: a(k) is element at (ia(k),ja(k))
    ! for k = 1:nelt

    integer, intent(in) :: n,nelt,ia(nelt),ja(nelt),isym
    real(kind=8), intent(in) :: x(n),a(nelt)
    real(kind=8), intent(out) :: y(n)

    integer :: k

    y = 0.d0
    do k=1,nelt
        y(ia(k)) = y(ia(k)) + a(k) * x(ja(k))
    enddo

end subroutine matvec
