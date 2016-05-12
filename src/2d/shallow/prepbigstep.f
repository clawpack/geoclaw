c
c-------------------------------------------------------------------------------------
c  Modified from prepbigstep.f in amrclaw to pass aux array into coarsen.f
c-------------------------------------------------------------------------------------
c
       subroutine prepbigstep(nvar,naux,lcheck,mptr,nx,ny,midub,mjdub,
     .                       valbgc,auxbgc,mi2tot,mj2tot)

       use amr_module, only: rnode,cornylo,cornxlo,evol,hxposs,hyposs
       use amr_module, only: NEEDS_TO_BE_SET,nghost,store2,timemult
       use amr_module, only: timeSetaux,possk,alloc,node,timeSetauxCPU
       use amr_module, only: auxtype
       implicit none

       integer, intent(in) :: nvar, naux, lcheck, mptr, nx, ny
       integer, intent(in) :: midub, mjdub, mi2tot, mj2tot
       real(kind=8), intent(inout) :: valbgc(nvar,mi2tot,mj2tot)
       real(kind=8), intent(inout) :: auxbgc(naux,mi2tot,mj2tot)

c      # Local variables
       real(kind=8) :: valdub(nvar,midub,mjdub)
       real(kind=8) :: auxdub(naux,midub,mjdub)
       real(kind=8) :: fp(nvar,mi2tot,mj2tot),gp(nvar,mi2tot,mj2tot)
       real(kind=8) :: fm(nvar,mi2tot,mj2tot),gm(nvar,mi2tot,mj2tot)
       real(kind=8) :: hx,hy,hx2,hy2,dt,dt2,dtnew2,time,tpre
       integer :: mitot,mjtot,ng2,locold,mx,my
       real(kind=8) :: xlow,ylow,xl,yb

       !for setaux timing
       integer :: clock_start, clock_finish, clock_rate
       real(kind=8) :: cpu_start, cpu_finish
c
          hx  = hxposs(lcheck)
          hy  = hyposs(lcheck)
          hx2 = 2.d0*hx
          hy2 = 2.d0*hy
          dt  = possk(lcheck)
          dt2 = 2. * dt
          time  = rnode(timemult,mptr)
          tpre  = time - dt
c
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          ng2    = 2*nghost
          locold = node(store2,mptr)
          xlow   = rnode(cornxlo,mptr) - nghost*hx2
          ylow   = rnode(cornylo,mptr) - nghost*hy2
c
c         # transfer soln. into grid with twice the ghost cells
          call copysol(valdub,alloc(locold),nvar,mitot,mjtot,
     1              nghost,midub,mjdub,ng2)

c
          if (naux .gt. 0) then
              xl     = rnode(cornxlo, mptr)
              yb     = rnode(cornylo, mptr)
              mx = midub - 4*nghost
              my = mjdub - 4*nghost
              auxdub(1,:,:) = NEEDS_TO_BE_SET  ! signal that needs a val (modified from amrclaw)
c
c             # Generate aux for fine grid
              call system_clock(clock_start, clock_rate)
              call cpu_time(cpu_start)
              call setaux(2*nghost,mx,my,xl,yb,hx,hy,
     &                    naux,auxdub)
              call system_clock(clock_finish, clock_rate)
              call cpu_time(cpu_finish)
              timeSetaux = timeSetaux + clock_finish - clock_start
              timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
c
c             # Generate aux for coarse grid
c              call auxcoarsen(auxdub,midub,mjdub,               # using auxcoarsen should work
c     1                     auxbgc,mi2tot,mj2tot,naux,auxtype)   # but gives strange results currently
              auxbgc(1,:,:) = NEEDS_TO_BE_SET  ! signal that needs a val
              call setaux(nghost,nx/2,ny/2,xl,yb,hx2,hy2,
     1                    naux,auxbgc)
          endif
c
c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,ng2,valdub,midub,mjdub,mptr,
     1               auxdub,naux)
c
c         coarsen by 2 in every direction
          call coarsen(valdub,midub,mjdub,auxdub,
     1                 valbgc,mi2tot,mj2tot,auxbgc,nvar,naux)
c
          call stepgrid(valbgc,fm,fp,gm,gp,
     1                mi2tot,mj2tot,nghost,
     2                dt2,dtnew2,hx2,hy2,nvar,
     3                xlow,ylow,tpre,mptr,naux,auxbgc)
c
c         update counts for error estimation work
          evol = evol + (nx/2)*(ny/2)
c
           return
           end