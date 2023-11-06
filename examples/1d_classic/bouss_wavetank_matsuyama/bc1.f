
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # Modified for 1d GeoClaw to extend topo aux(1,:) to ghost cells.
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      use geoclaw_module, only: grav
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)

      pi = 4.d0*atan(1.d0)
      h0 = 4.d0
      
      ampl_eta = 0.03
      ampl_u = ampl_eta * sqrt(grav/h0)
      !write(6,*) '+++ ampl_eta, ampl_u: ',ampl_eta, ampl_u
      !ampl = 0.047  #0.15d0
      tperiod = 20.d0

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,1-ibc)=aux(1,1)
      enddo

      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # wave maker
      if (t < tperiod) then
          s = ampl_u * sin(2.d0*pi*t/tperiod)
      else
          s = 0.d0
      endif
      do ibc=1,mbc
          aux(1,1-ibc) = aux(1,ibc)
          do m=1,meqn
              q(m,1-ibc) = q(m,ibc)
          enddo
          u = q(2,ibc)/q(1,ibc)
          q(2,1-ibc) = (2.d0*s - u) * q(1,1-ibc)
      enddo
      
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
c        # negate the normal velocity:
         q(2,1-ibc) = -q(2,ibc)
      end do
      go to 199

  199 continue

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,mx+ibc)=aux(1,mx)
      enddo

      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
c        # negate the normal velocity:
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end
