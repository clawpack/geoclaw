
subroutine setprob

    ! for problem specific features, copy to your directory and modify

    use grid_module, only: runup_tolerance
    implicit none
    
    runup_tolerance = 1.d-4   ! for wave tank, since dry_tolerance = 1d-6

end subroutine setprob
