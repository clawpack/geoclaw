
#define IADD(IVAR,I,J,LOC,NVAR,MITOT) LOC+IVAR-1+NVAR*((J-1)*MITOT+I-1)
#define IADDF(IVAR,I,J,LOCF,NVAR,MI) LOCF+IVAR-1+NVAR*((J-1)*MI+I-1)
#define IADDFAUX(I,J,LOCFAUX,MCAPA,NAUX,MI) LOCFAUX+MCAPA-1+NAUX*((J-1)*MI+(I-1))
#define IADDCAUX(I,J,LOCCAUX,MCAPA,NAUX,MITOT) LOCCAUX+MCAPA-1+NAUX*((J-1)*MITOT+(I-1))
#define IADDFTOPO(I,J,LOCFAUX,NAUX,MI) LOCFAUX+NAUX*((J-1)*MI+(I-1)) 
#define IADDCTOPO(I,J,LOCCAUX,NAUX,MITOT) LOCCAUX+NAUX*((J-1)*MITOT+(I-1)) 

!
! :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
! update - update all grids at level 'level'.
!          this routine assumes cell centered variables.
!          the update is done from 1 level finer meshes under it.
! input parameter:
!    level  - ptr to the only level to be updated. levels coarser than
!             this will be at a diffeent time.
! rewritten in July 2016 from the original update.f
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

subroutine update (level, nvar, naux)

    use geoclaw_module, only: dry_tolerance
    use amr_module
    implicit none
    ! modified for shallow water on topography to use surface level eta
    ! rather than depth h = q(i,j,1) 
    ! eta - q(i,j,1) + aux(i,j,1)

    ! inputs
    integer, intent(in) :: nvar, naux, level
!     real(kind=8), intent(in) :: level

    integer :: ng, levSt, mptr, loc, locaux, nx, ny, mitot, mjtot
    integer :: ilo, jlo, ihi, jhi, mkid, iclo, jclo, ichi, jchi
    integer :: mi, mj, locf, locfaux, iplo, jplo, iphi, jphi
    integer :: iff, jff, nwet, ico, jco, i, j, ivar, loccaux
    integer :: listgrids(numgrids(level)), lget
    real(kind=8) :: dt, totrat, bc, etasum, hsum, husum, hvsum
    real(kind=8) :: hf, bf, huf, hvf, etaf, hav, hc, huc, hvc, capa, etaav
    real(kind=8) :: capac
    character(len=80) :: String

    lget = level

    String = "(19h    updating level ,i5)"
    
    if (uprint) then
        write(outunit, String) lget
    endif

    dt = possk(lget)

#ifdef CUDA
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                     ilo,jlo,ihi,jhi,mkid,iclo,jclo, &
!$OMP                     ichi,jchi,mi,mj,locf,locfaux, &
!$OMP                     iplo,jplo,iphi,jphi,iff,jff,totrat,i,j, &
!$OMP                     capac,bc,etasum,hsum,husum,hvsum, &
!$OMP                     hf,bf,huf,hvf, &
!$OMP                     etaf,etaav,hav,nwet,hc,huc,hvc, &
!$OMP                     ivar,ico,jco,capa,levSt, &
!$OMP                     String), &
!$OMP          SHARED(lget,numgrids,listsp,alloc,nvar,naux, &
!$OMP                    intratx,intraty,nghost,uprint,mcapa,node, &
!$OMP                    listOfGrids,listStart,lstart,level,cflux_hh,dry_tolerance), &
!$OMP          DEFAULT(none)
#else
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                    ilo,jlo,ihi,jhi,mkid,iclo,ichi, &
!$OMP                    jclo,jchi,mi,mj,locf,locfaux,iplo,iphi, &
!$OMP                    jplo,jphi,iff,jff,totrat,i,j,ivar,capac, &
!$OMP                    capa,bc,etasum,hsum,husum,hvsum, &
!$OMP                    levSt,ico,jco,hf,bf,huf,hvf, &
!$OMP                    etaf,etaav,hav,nwet,hc,huc,hvc, String), &
!$OMP             SHARED(numgrids,listOfGrids,level,intratx,intraty, &
!$OMP                   nghost,uprint,nvar,naux,mcapa,node,listsp, &
!$OMP                   alloc,lstart,dry_tolerance,listStart,lget), &
!$OMP            DEFAULT(none)
#endif

    ! need to set up data structure for parallel distribution of grids
    ! call prepgrids(listgrids, numgrids(level), level)

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
        if (associated(cflux_hh(mptr)%ptr)) then
#else
        if (node(cfluxptr, mptr) /= 0) then
#endif

#ifdef CUDA
            call upbnd(cflux_hh(mptr)%ptr,alloc(loc),nvar, &
                    naux,mitot,mjtot,listsp(lget),mptr)
#else
            call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,naux,mitot, &
                    mjtot,listsp(lget),mptr)
#endif
        endif

        mkid = lstart(lget+1)
        do while (mkid /= 0)
            iclo   = node(ndilo,mkid)/intratx(lget)
            jclo   = node(ndjlo,mkid)/intraty(lget)
            ichi   = node(ndihi,mkid)/intratx(lget)
            jchi   = node(ndjhi,mkid)/intraty(lget)
            mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
            mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
            locf    = node(store1,mkid)
            locfaux = node(storeaux,mkid)

            ! calculate the starting adn ending indices for coarse grid update, if overlap
            iplo = max(ilo,iclo)
            jplo = max(jlo,jclo)
            iphi = min(ihi,ichi)
            jphi = min(jhi,jchi)        

            if (iplo <= iphi .and. jplo <= jphi) then
                !calculate the starting index for the fine grid source pts.
                iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
                jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
                totrat = intratx(lget) * intraty(lget)

                do i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
                    do j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
                        if (uprint) then
                            String = "(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,' of grid ',i4)"
                            write(outunit,String) i,j,mptr,iff,jff,mkid                
                            
                            String = "(' old vals: ',4e25.15)"
                            write(outunit,String)(alloc(IADD(ivar,i,j,loc,nvar,mitot)),ivar=1,nvar)
                        endif
                        
                        if (mcapa == 0) then
                            capac = 1.0d0
                        else
                            capac = alloc(IADDCAUX(i,j,loccaux,mcapa,naux,mitot))
                        endif

                        bc = alloc(IADDCTOPO(i,j,loccaux,naux,mitot))

                        etasum = 0.d0
                        hsum = 0.d0
                        husum = 0.d0
                        hvsum = 0.d0

                        nwet = 0

                        do jco = 1, intraty(lget)
                            do ico = 1, intratx(lget)
                                if (mcapa == 0) then
                                    capa = 1.0d0
                                else
                                    capa = alloc(IADDFAUX(iff+ico-1,jff+jco-1,locfaux,mcapa,naux,mi))
                                endif

                                hf = alloc(IADDF(1,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 
                                bf = alloc(IADDFTOPO(iff+ico-1,jff+jco-1,locfaux,naux,mi))*capa
                                huf= alloc(IADDF(2,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 
                                hvf= alloc(IADDF(3,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa 

                                if (alloc(IADDF(1,iff+ico-1,jff+jco-1,locf,nvar,mi)) > dry_tolerance) then
                                    etaf = hf + bf
                                    nwet = nwet + 1
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

                        if (nwet > 0) then
                            etaav = etasum/dble(nwet)
                            hav = hsum/dble(nwet)
                            hc = min(hav, (max(etaav-bc*capac, 0.0d0)))
                            huc = (min(hav, hc) / hsum) * husum
                            hvc = (min(hav, hc) / hsum) * hvsum
                        else
                            hc = 0.0d0
                            huc = 0.0d0
                            hvc = 0.0d0
                        endif

                        alloc(IADD(1,i,j,loc,nvar,mitot)) = hc / capac 
                        alloc(IADD(2,i,j,loc,nvar,mitot)) = huc / capac 
                        alloc(IADD(3,i,j,loc,nvar,mitot)) = hvc / capac 

                        if (uprint) then
                            String = "(' new vals: ',4e25.15)"
                            write(outunit, String)(alloc(IADD(ivar, i, j, loc, nvar,mitot)), ivar=1, nvar)
                        endif

                        jff = jff + intraty(lget)
                    enddo

                    iff = iff + intratx(lget)
                    jff = jplo * intraty(lget) - node(ndjlo, mkid) + nghost + 1
                enddo

                
            endif
            mkid = node(levelptr, mkid)

        enddo
        continue
    enddo
!$OMP END PARALLEL DO

    return


contains



end subroutine
