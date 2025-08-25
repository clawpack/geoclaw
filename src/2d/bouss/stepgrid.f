c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &                  nvar,xlow,ylow,time,mptr,maux,aux,actualstep)
c
c
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q.
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c
c
c      This version of stepgrid, stepgrid_geo.f allows output on
c      fgout grids specified in fgout_grids.data
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use geoclaw_module
      use amr_module
      use fgout_module, only: FGOUT_num_grids, FGOUT_fgrids, 
     &                        FGOUT_tcfmax, fgout_interp, fgout_grid,
     &                        FGOUT_ttol
      implicit double precision (a-h,o-z)

      external rpn2,rpt2



      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
      dimension eta(mitot,mjtot)
c      dimension work(mwork)

      logical :: debug = .false.
      logical :: dump = .false.
      !logical :: dump = .true.
      logical, intent (in) :: actualstep
      type(fgout_grid), pointer :: fgout
      !logical, allocatable :: fgout_interp_needed(:)
      logical :: fgout_interp_needed(FGOUT_num_grids)
      real(kind=8) :: fgout_tnext
c
#ifdef WHERE_AM_I
      write(*,*) "     starting stepgrid grid ",mptr 
#endif
      FGOUT_tcfmax = -rinfinity
      level = node(nestlevel,mptr)


      if (dump) then
         eta = q(1,:,:) + aux(1,:,:)
         write(outunit,*)" at start of stepgrid: dumping grid ",mptr,
     &        "at time ",time
         do i = 1, mitot
         do j = 1, mjtot
!           write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar),
!    .                         (aux(iaux,i,j),iaux=1,maux)
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar),eta(i,j)
 545        format(2i4,6e15.7,/,8x,4e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c
c     # fluxes initialized in step2
c
C       mwork0 = (maxm+2*mbc)*(12*meqn + mwaves + meqn*mwaves + 2)
C c
C       if (mwork .lt. mwork0) then
C          write(outunit,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          write(*      ,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          stop
C       endif
c
c     # partition work array into pieces needed for local storage in
c     # step2 routine. Find starting index of each piece:
c
C       i0faddm = 1
C       i0faddp = i0faddm + (maxm+2*mbc)*meqn
C       i0gaddm = i0faddp + (maxm+2*mbc)*meqn
C       i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
C       i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn
C       i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn
C       i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
C       i0aux1 = i0dtdy1 + (maxm+2*mbc)
C       i0aux2 = i0aux1 + (maxm+2*mbc)*maux
C       i0aux3 = i0aux2 + (maxm+2*mbc)*maux
C c
C c
C       i0next = i0aux3 + (maxm+2*mbc)*maux    !# next free space
C       mused  = i0next - 1                    !# space already used
C       mwork1 = mwork - mused              !# remaining space (passed to step2)

c
      tc0=time !# start of computational step
      tcf=tc0+dt !# end of computational step

c     Check if fgout interpolation needed before and after step:

      !allocate(fgout_interp_needed(FGOUT_num_grids))

!    !$OMP CRITICAL (FixedGrids)

      do ng=1,FGOUT_num_grids
          fgout => FGOUT_fgrids(ng)
          if (fgout%next_output_index > fgout%num_output) then
              fgout_interp_needed(ng) = .false.
          else
              fgout_tnext = fgout%output_times(fgout%next_output_index)
              fgout_interp_needed(ng) =  
     &          ((fgout%x_low < xlowmbc + mx * dx) .and.
     &           (fgout%x_hi  > xlowmbc) .and.
     &           (fgout%y_low < ylowmbc + my * dy) .and.
     &           (fgout%y_hi  > ylowmbc) .and.
     &           (fgout_tnext >= tc0 - FGOUT_ttol) .and.
     &           (fgout_tnext <= tcf + FGOUT_ttol))
          endif
c         write(6,*) '+++ level, before- tc0,tcf,needed: ',
c    &                level,tc0,tcf,fgout_interp_needed(ng)
c         write(6,*) '+++ next index: ',fgout%next_output_index
c         write(6,*) '+++ fgout_tnext, tcf + FGOUT_ttol: ',
c    &                fgout_tnext, tcf + FGOUT_ttol
      enddo
!      !$OMP END CRITICAL (FixedGrids)


c::::::::::::::::::::::::fgout Output:::::::::::::::::::::::::::::::::
c     This has been moved to tick.f, after advancing all patches on 
c     finest level.  No need to check on each patch separately.
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c      This call has been moved out to advanc
c       call b4step2(mbc,mx,my,nvar,q,
c     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux,actualstep)
      
c::::::::::::::::::::::::FGOUT DATA before step:::::::::::::::::::::::
c     # fill in values at fgout points affected at time tc0
      do ng=1,FGOUT_num_grids
          if (fgout_interp_needed(ng)) then
c            # fgout grid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
             fgout => FGOUT_fgrids(ng)
!$OMP        CRITICAL (FixedGrids)
c            write(6,*) '+++ fout_interp(1), tc0, level: ',tc0,level
             call fgout_interp(1,fgout,tc0,q,nvar,mx,my,mbc,
     &                         dx,dy,xlowmbc,ylowmbc,maux,aux)

!$OMP       END CRITICAL (FixedGrids)
      endif
      enddo
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     New fixed grid stuff: Update fixed grid info from this patch...

      call fgmax_frompatch(mx,my,nvar,mbc,maux,q,aux,
     &     dx,dy,xlowmbc,ylowmbc,level,time,time+dt)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     # take one step on the conservation law:
c
      call step2(mbig,nvar,maux,
     &           mbc,mx,my,
     &              q,aux,dx,dy,dt,cflgrid,
     &              fm,fp,gm,gp,rpn2,rpt2)
c
c            
        mptr_level = node(nestlevel,mptr)

c       write(outunit,811) mptr, mptr_level, cflgrid
c811    format(" Courant # of grid ",i5," level",i3," is ",d12.4)
c

!$OMP  CRITICAL (cflm)

        cflmax = dmax1(cflmax,cflgrid)
        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c
!        write(outunit,*)" updating grid ", mptr
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy
         if (mcapa.eq.0) then
          do 50 j=mbc+1,mjtot-mbc
          do 50 i=mbc+1,mitot-mbc
          do 50 m=1,nvar
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j)
     &           - dtdx*(fm(m,i+1,j) - fp(m,i,j))
     &           - dtdy*(gm(m,i,j+1) - gp(m,i,j))
 50       continue
         else
          do 51 j=mbc+1,mjtot-mbc
          do 51 i=mbc+1,mitot-mbc
          do 51 m=1,nvar
c            # with capa array.
           q(m,i,j) = q(m,i,j)
     &           - (dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &           +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
 51       continue
!           write(outunit,543) m,i,j,q(m,i,j),fm(m,i+1,j),fp(m,i,j),
!     .        gm(m,i,j+1), gp(m,i,j)
543       format(3i4,5e25.16)

         endif

c 50      continue
c
c     # Copied here from b4step2 since need to do before saving to qc1d:
      forall(i=1:mitot, j=1:mjtot, q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:meqn,i,j) = 0.d0
      end forall
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
         call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif

c     ::::::::::::::::::::::::fgout data afterstep:::::::::::::::::::::::
c     # fill in values at fgout points affected at time tcf
      do ng=1,FGOUT_num_grids
          if (fgout_interp_needed(ng)) then
c            # fgout grid ng has an output time within [tc0,tcf] interval
c        # and it overlaps this computational grid spatially
!$OMP        CRITICAL (FixedGrids)
             fgout => FGOUT_fgrids(ng)
c            write(6,*) '+++ fout_interp(2), tcf, level: ',tcf,level
             call fgout_interp(2,fgout,tcf,q,nvar,mx,my,mbc,
     &                         dx,dy,xlowmbc,ylowmbc,maux,aux)

!$OMP       END CRITICAL (FixedGrids)
      endif
c        write(6,*) '+++ level,after- tc0,tcf,needed: ',
c    &                level,tc0,tcf,fgout_interp_needed(ng)
c        write(6,*) '+++ next index: ',fgout%next_output_index
c        write(6,*) '+++ fgout_tnext: ',fgout_tnext
      enddo
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
         do 830 i = mbc+1, mitot-1
            do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j),
     .                                 gm(1,i,j),gp(1,i,j)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i,j),fp(m,i,j),
     .            gm(m,i,j),gp(m,i,j)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to
c choose the allowable new time step dtnew.  This will later be
c compared to values seen on other grids.
c
       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
         else
c          # velocities are all zero on this grid so there's no
c          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

c     # give a warning if Courant number too large...
c
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid,cflv1,mptr,mptr_level
            write(outunit,810) cflgrid, cflv1,mptr,mptr_level
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than input cfl_max = ', d12.4,
     &              '  on grid ',i3, ' level ',i3)
            endif
c
      if (dump) then
         eta = q(1,:,:) + aux(1,:,:)
         write(outunit,*)" at end of stepgrid: dumping grid ",mptr,
     &        "at time ",time
         do i = mbc+1, mitot-mbc
         do j = mbc+1, mjtot-mbc
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar),eta(i,j)
c            write(*,545) i,j,(q(i,j,ivar),ivar=1,nvar)
         end do
         end do
      endif

#ifdef WHERE_AM_I
      write(*,*) "     ending   stepgrid grid ",mptr 
#endif
      return
      end


