c-----------------------------------------------------c
      subroutine stepp(factor)
c-----------------------------------------------------c
c     Momentum step                                   c
c     PI(new) = PI(old) - (1/3)*Beta*dTau*T(Plaq.)    c
c      where T(...) is a projection onto Algebra.     c
c-----------------------------------------------------c
      USE hmc_mod
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u(4)
      TYPE(g_field1) staple, tmp1, tmp2, tmp3
      common/ config/ u

      REAL*8 factor
      TYPE(a_field) c 

      Direction: do mu = 1, 4

         call make_staple(u,staple,mu)   ! Staple for Wilson action
*        call make_istaple(u,staple,mu)  ! Staple for Improved action

         tmp1 = u(mu)         !  tmp1 
c
c                   ---
c                  |   |           
c            tmp2 = ---    +   ---
c                             |   |
c                              ---
c
         tmp2 = tmp1.prodAD.staple     ! u*(staple)_adj
c
c----------------------------------------------------------------------c

c     .....   Projection onto Lie Algebra   .....

         call  projlink(tmp2,tmp3,NV)

         call  algblink(tmp3,c,NV)

c     .....   p(new) = p(old) + fac * c  .....
         p(mu) = p(mu) + ((factor*fac)*c)

      enddo Direction

      return
      end

c-----------------------------------------------------c
      subroutine stepu(factor)
c-----------------------------------------------------c
c     U(new) = exp(dTau*Pi(new)) * U(old)             c
c-----------------------------------------------------c
      USE hmc_mod
      USE field_g
      include '../INCLUDE/para_geometry'

      REAL*8 factor
      TYPE(a_field)  c
      TYPE(g_field0) u(4)
      TYPE(g_field1) tmp1,tmp2,tmp3
      common/ config/ u

      do mu = 1, 4

         tmp1=u(mu)

c     .....   Construct exp(Force)   .....
c      . .    c = dtau * p   . .
         c = (factor*dtau) * p(mu) 

c     . .   tmp2 = exp( i*c(1)*Lambda(1)/2) + i*c(2)*Lambda(2)/2) + ..)
         call  gprojct(c,tmp2)

c     .....   Construct U(new)   .....
         u(mu) = tmp2*tmp1

         call set_wing_g2(u(mu))     

      enddo

      return
      end

c-----------------------------------------------------c
      subroutine steppf(factor)
c-----------------------------------------------------c
c     Molecular Dynamics Step for Pseudo-Fermion      
c     field \phi
c     hop  : Hopping parameter
c     r    : Wilson term
c     eta  : Gaussian distributed Vector
c     X    : W*xi = eta
c-----------------------------------------------------c
      USE debug1    ! for Debug
      USE hmc_mod
      USE field_g
      USE fpara
      include '../INCLUDE/para_geometry'
      REAL*8 factor
*     TYPE(f_field)  tmp1
      TYPE(f_field)  tmp0, tmp1
      TYPE(g_field1) tmp2, tmp3
      TYPE(a_field) c, c2

      TYPE(g_field0) u(4)
      common/ config/ u
      real*8 plaq

c     ...  Solve W^dagger*eta = phi
*     call cg0(eta,phi,eps,imax,2)
      call cg0(eta,phi,2)

c     ...  Fill boundaries
      call set_wing_f1(eta)

c     ...  Solve W*Xi = eta
*     call cg0(xi,eta,eps,imax,1)
      call cg0(xi,eta,1)

c     ...  Fill boundaries
      call set_wing_f1(xi)

      do mu = 1, 4

        !  Construct U(x,mu)*P1
*       tmp1 = mu.fshift.xi 
        tmp0 = mu.fshift.xi 
        tmp1 = hopp(mu)*tmp0
        call vvmat(tmp1,eta,tmp2,1)

        ! ...   Projection onto Lie Algebra   ...
        call  projlink(tmp2,tmp3,NV)
        call  algblink(tmp3,c,NV)
        ! ...  p(new) = p(old) + fac * c  .....
*       p(mu) = p(mu) + ((factor*hopp(mu)*dtau)*c)
        p(mu) = p(mu) + ((factor*dtau)*c)

        !  Construct P2*U_adj(x,mu)
*       tmp1 = (-mu).fshiftB.eta 
        tmp0 = (-mu).fshiftB.eta 
        tmp1 = hopm(mu)*tmp0
        call vvmat(xi,tmp1,tmp2,2)

        ! ...   Projection onto Lie Algebra   ...
        call  projlink(tmp2,tmp3,NV)
        call  algblink(tmp3,c,NV)

        ! ...  p(new) = p(old) + fac * c  .....
*        p(mu) = p(mu) - ((factor*hopm(mu)*dtau)*c)
         p(mu) = p(mu) - ((factor*dtau)*c)

        ! ... Add Clover contribution ...
        IF( ABS(Csw) > 0.00001d0 ) THEN
           CALL dSclover(mu,c2,u)
           p(mu) = p(mu) - ((factor*dtau)*c2)
        ENDIF

      enddo

      return
      end

c-----------------------------------------------------c
      subroutine vvmat(v1,v2,vv,iflag)
c-----------------------------------------------------c
c     iflag = 1
c       vv = |v1><v2| i.e., vv(a,b) = v1(a)*v2_adj(b)
c     iflag = 2
c                           vv(a,b) = v1(a)*v2(b)
c-----------------------------------------------------c
      USE field_g
      USE field_f
      include '../INCLUDE/para_geometry'

      TYPE(f_field), INTENT(IN)::  v1, v2
      TYPE(g_field1),INTENT(OUT):: vv

      if(iflag==1)  then

        do ia = 1, 3
        do ib = 1, 3

        is = 0

          do it = 1, NT
          do iz = 1, NZ
          do iy = 1, NY
          do ix = 1, NX

          is = is + 1

          vv%g(ia,ib,is) 
     &      = v1%f(ia,ix,iy,iz,it,1)*conjg(v2%f(ib,ix,iy,iz,it,1))
     &      + v1%f(ia,ix,iy,iz,it,2)*conjg(v2%f(ib,ix,iy,iz,it,2))
     &      + v1%f(ia,ix,iy,iz,it,3)*conjg(v2%f(ib,ix,iy,iz,it,3))
     &      + v1%f(ia,ix,iy,iz,it,4)*conjg(v2%f(ib,ix,iy,iz,it,4))

          enddo
          enddo
          enddo
          enddo

        enddo
        enddo

      endif

      if(iflag==2)  then

        do ia = 1, 3
        do ib = 1, 3

        is = 0

          do it = 1, NT
          do iz = 1, NZ
          do iy = 1, NY
          do ix = 1, NX

          is = is + 1

          vv%g(ia,ib,is) 
     &      = v1%f(ia,ix,iy,iz,it,1)*v2%f(ib,ix,iy,iz,it,1)
     &      + v1%f(ia,ix,iy,iz,it,2)*v2%f(ib,ix,iy,iz,it,2)
     &      + v1%f(ia,ix,iy,iz,it,3)*v2%f(ib,ix,iy,iz,it,3)
     &      + v1%f(ia,ix,iy,iz,it,4)*v2%f(ib,ix,iy,iz,it,4)

          enddo
          enddo
          enddo
          enddo

        enddo
        enddo

      endif

      return
      end
