c
!> print solution and aux. variables to output. 
c =======================================================================
      subroutine outval(val,nvar,mitot,mjtot,mptr,outgrd,naux,aux)
c =======================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension  val(nvar,mitot,mjtot)
      dimension  aux(naux,mitot,mjtot)
      logical    outgrd


c ::::::::::::::::::::::OUTVAL :::::::::::::::::::::::::::::::
c print solution and aux. variables to output. 
c if cell outside domain, don't print soln. value - nothing
c currently in ghost cells.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (.not. outgrd) go to 99
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      cornx = rnode(cornxlo,mptr) -  nghost*hx
      corny = rnode(cornylo,mptr) -  nghost*hy
c
      do 25 j=nghost+1,mjtot-nghost
      do 20 i=nghost+1,mitot-nghost

          x  = cornx + hx*(dble(i)-.5d0)
          y  = corny + hy*(dble(j)-.5d0)
          write(outunit,107) x,y,i,j,(val(ivar,i,j),ivar=1,nvar)
c107      format(2hx=,e9.3,2hy=,e9.3,3h,i=,i3,3h,j=,i3,' a=',
 107      format(e9.3,2x,e9.3,i5,2x,i5,5(e14.7,1x))
c         if (naux.gt.0) write(outunit,108) (aux(iaux,i,j),iaux=1,naux)
 108      format(1x,'aux = ',7(e10.3,1x))

 20   continue
 25   continue

 99   return
      end
