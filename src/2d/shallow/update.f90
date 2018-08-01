
#define IADD_UP(IVAR,I,J,LOC,NVAR,MITOT) LOC+IVAR-1+NVAR*((J-1)*MITOT+I-1)
#define IADDF_UP(IVAR,I,J,LOCF,NVAR,MI) LOCF+IVAR-1+NVAR*((J-1)*MI+I-1)
#define IADDFAUX_UP(I,J,LOCFAUX,MCAPA,NAUX,MI) LOCFAUX+MCAPA-1+NAUX*((J-1)*MI+(I-1))
#define IADDCAUX_UP(I,J,LOCCAUX,MCAPA,NAUX,MITOT) LOCCAUX+MCAPA-1+NAUX*((J-1)*MITOT+(I-1))
#define IADDFTOPO_UP(I,J,LOCFAUX,NAUX,MI) LOCFAUX+NAUX*((J-1)*MI+(I-1)) 
#define IADDCTOPO_UP(I,J,LOCCAUX,NAUX,MITOT) LOCCAUX+NAUX*((J-1)*MITOT+(I-1)) 

!
! -----------------------------------------------------------
!
      subroutine update (level, nvar, naux)
      
      use geoclaw_module, only: dry_tolerance
!     # modified for shallow water on topography to use surface level eta
!     # rather than depth h = q(i,j,1)
!     # eta = q(i,j,1) + aux(i,j,1)
!
      use amr_module
      implicit double precision (a-h,o-z)


      integer listgrids(numgrids(level))

!
! :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
! update - update all grids at level 'level'.
!          this routine assumes cell centered variables.
!          the update is done from 1 level finer meshes under it.
! input parameter:
!    level  - ptr to the only level to be updated. levels coarser than
!             this will be at a diffeent time.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      lget = level

      if (uprint) write(outunit,100) lget
100   format(19h    updating level ,i5)
!     need to set up data structure for parallel distrib of grids
!     call prepgrids(listgrids,numgrids(level),level)
!
!  grid loop for each level
!
      dt     = possk(lget)

!      mptr = lstart(lget)
! 20   if (mptr .eq. 0) go to 85

#ifdef CUDA
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                     ilo,jlo,ihi,jhi,mkid,iclo,jclo, &
!$OMP                     ichi,jchi,mi,mj,locf,locfaux, &
!$OMP                     iplo,jplo,iphi,jphi,iff,jff,totrat,i,j, &
!$OMP                     capac,bc,etasum,hsum,husum,hvsum,drytol, &
!$OMP                     newt,hf,bf,huf,hvf, &
!$OMP                     etaf,etaav,hav,nwet,hc,huc,hvc, &
!$OMP                     ivar,ico,jco,capa,levSt), &
!$OMP          SHARED(lget,numgrids,listsp,alloc,nvar,naux, &
!$OMP                    intratx,intraty,nghost,uprint,mcapa,node, &
!$OMP                    listOfGrids,listStart,lstart,level,cflux_hh,dry_tolerance), &
!$OMP          DEFAULT(none)
#else
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                     ilo,jlo,ihi,jhi,mkid,iclo,jclo, &
!$OMP                     ichi,jchi,mi,mj,locf,locfaux, &
!$OMP                     iplo,jplo,iphi,jphi,iff,jff,totrat,i,j, &
!$OMP                     capac,bc,etasum,hsum,husum,hvsum,drytol, &
!$OMP                     newt,hf,bf,huf,hvf, &
!$OMP                     etaf,etaav,hav,nwet,hc,huc,hvc, &
!$OMP                     ivar,ico,jco,capa,levSt), &
!$OMP          SHARED(lget,numgrids,listOfGrids,listsp,alloc,nvar,naux, &
!$OMP                    intratx,intraty,nghost,uprint,mcapa,node, &
!$OMP                    listStart,lstart,dry_tolerance,level), &
!$OMP          DEFAULT(none)
#endif


      do ng = 1, numgrids(lget)
         levSt   = listStart(lget)
         mptr    = listOfGrids(levSt+ng-1)
         loc     = node(store1,mptr)
         loccaux = node(storeaux,mptr)
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         ilo     = node(ndilo,mptr)
         jlo     = node(ndjlo,mptr)
         ihi     = node(ndihi,mptr)
         jhi     = node(ndjhi,mptr)

#ifdef CUDA
         if (associated(cflux_hh(mptr)%ptr) .eq. .false.) go to 25
#else
         if (node(cfluxptr,mptr) .eq. 0) go to 25
#endif

! #ifdef CUDA
!          call upbnd(cflux_hh(mptr)%ptr,alloc(loc),nvar, &
!                     naux,mitot,mjtot,listsp(lget),mptr)
! #else
!          call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar, &
!                     naux,mitot,mjtot,listsp(lget),mptr)
! #endif
!
!  loop through all intersecting fine grids as source updaters.
!
 25      mkid = lstart(lget+1)
 30        if (mkid .eq. 0) go to 80
           iclo   = node(ndilo,mkid)/intratx(lget)
           jclo   = node(ndjlo,mkid)/intraty(lget)
           ichi   = node(ndihi,mkid)/intratx(lget)
           jchi   = node(ndjhi,mkid)/intraty(lget)

           mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
           mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
           locf    = node(store1,mkid)
           locfaux = node(storeaux,mkid)
!
!  calculate starting and ending indices for coarse grid update, if overlap
!
         iplo = max(ilo,iclo)
         jplo = max(jlo,jclo)
         iphi = min(ihi,ichi)
         jphi = min(jhi,jchi)

         if (iplo .gt. iphi .or. jplo .gt. jphi) go to 75
!
!  calculate starting index for fine grid source pts.
!
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
         totrat = intratx(lget) * intraty(lget)

         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
         do 70 j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
           if (uprint) then
              write(outunit,101) i,j,mptr,iff,jff,mkid
 101          format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4, &
                    ' of grid ',i4)
              write(outunit,102)(alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)), &
                ivar=1,nvar)
 102          format(' old vals: ',4e25.15)
           endif
!
!
!  update using intrat fine points in each direction
!
!
!
!     # modified by dlg 11/01/05 below:
!
!     Strategy:
!     average eta,h,hu,hv in wet cells only.
!     set hcoarse = ave(eta_wet) - b_coarse, or zero if that is negative

!     set hucoarse = (min(hcoarse,h_average))/h_average)* hu_average
!     note: this conserves mass and momentum as long as all fine cells are wet.
!     Momentum is reduced by the same factor as mass in the coarse cell
!     and is never increased given an increase in mass

      if (mcapa .eq. 0) then
         capac=1.0d0
      else
         capac=alloc(IADDCAUX_UP(i,j,loccaux,mcapa,naux,mitot))
         endif

      bc = alloc(IADDCTOPO_UP(i,j,loccaux,naux,mitot))

      etasum = 0.d0
      hsum = 0.d0
      husum = 0.d0
      hvsum = 0.d0

      nwet=0

      do jco  = 1, intraty(lget)
         do ico  = 1, intratx(lget)
            if (mcapa .eq. 0) then
               capa=1.0d0
            else
               capa=alloc(IADDFAUX_UP(iff+ico-1,jff+jco-1,locfaux,mcapa,naux,mi))
            endif
            hf = alloc(IADDF_UP(1,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 
            bf = alloc(IADDFTOPO_UP(iff+ico-1,jff+jco-1,locfaux,naux,mi))*capa
            huf= alloc(IADDF_UP(2,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 
            hvf= alloc(IADDF_UP(3,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 

            if (hf > dry_tolerance) then
               etaf = hf+bf
               nwet=nwet+1
            else
               etaf = 0.d0
               huf=0.d0
               hvf=0.d0
               endif

               hsum   = hsum + hf
               husum  = husum + huf
               hvsum  = hvsum + hvf
               etasum = etasum + etaf
            enddo
         enddo

      if (nwet.gt.0) then
         etaav=etasum/dble(nwet)
         hav= hsum/dble(nwet)
!         hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
         hc=min(hav,(max(etaav-bc*capac,0.d0)))
         huc=(min(hav,hc)/hsum)*husum
         hvc=(min(hav,hc)/hsum)*hvsum
      else
         hc=0.d0
         huc=0.d0
         hvc=0.d0
         endif

!     # set h on coarse grid based on surface, not conservative near shoreline
      
      alloc(IADD_UP(1,i,j,loc,nvar,mitot)) = hc / capac 
      alloc(IADD_UP(2,i,j,loc,nvar,mitot)) = huc / capac 
      alloc(IADD_UP(3,i,j,loc,nvar,mitot)) = hvc / capac 
!
      if (uprint) write(outunit,103)(alloc(IADD_UP(1,i,j,loc,nvar,mitot)), &
          ivar=1,nvar)
 103  format(' new vals: ',4e25.15)
!
      jff = jff + intraty(lget)
 70   continue
      iff = iff + intratx(lget)
      jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
 71   continue
!
 75   mkid = node(levelptr,mkid)
      go to 30
!
 80   continue
      end do
!$OMP END PARALLEL DO

 99   return
      end
