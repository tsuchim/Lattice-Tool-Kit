c cg.f
c----------------------------------------------------------------------
*     subroutine  cg0 (x,b,eps,imax,iflag)  ! Changed on Jan.7, 2006
      subroutine  cg0 (x,b,iflag)
c-----------------------------------------------------------------------
c     this routine  solves
c                  W(u) * x = b      for iflag=1
c                  W_adj(u) * x = b  for iflag=2
c     by conjugate gradient 
c
c     eps and imax controle when we stop CG iteration, i.e.,
c     <res|res> < eps or (# of CG iteration)>imax, the iteration is
c     terminated.
c----------------------------------------------------------------------
      USE field_f
      USE fpara             ! Added on Jan.7, 2006
      TYPE(g_field0) U(4)   ! Added on Jan.7, 2006
      common/ config/ U     ! Added on Jan.7, 2006
c      include 'mpif.h'                                                ! MPI
c      include '../INCLUDE/para.h'                                     ! MPI


*     REAL*8 eps
      TYPE(f_field) wxvect, x, b, res, p, q, s
      REAL*8 alpha, beta, c1, c2, c3

      data icheck/ 1/

      IF( abs(Csw) > 0.001 ) THEN
         CALL MKFmunu(U,hop,Csw)
      ENDIF

c     the initial condition  ( for i=0 )
      ! ...  res = b - W*x
      if(iflag==1)  then
         res = b - wxvect(x,2)
      else if(iflag==2)  then
         res = b - wxvect(x,3)
      endif

      rnorm = res * res          ! added on Sept.11, 2012 by AN
      if (rnorm < eps)  then     !
        idone = i                !
        go to 5000               !
      endif                      !

      ! ...  p = W_adj * res
      if(iflag==1)  then
         p = wxvect(res,3)
      else if(iflag==2)  then
         p = wxvect(res,2)
      endif

c
      ! ...  c1 = < p | p >
      c1 = p * p

      ! ...  the iteration starts

      do i = 1, imax

        ! ...  q = W * p
        if(iflag==1)  then
           q = wxvect(p,2)
        else if(iflag==2)  then
           q = wxvect(p,3)
        endif

        ! ...  c2 = < q | q >
        c2 = q * q 

        alpha = c1 / c2

        ! ...  x   = x   + alpha * p   
        x = x + (alpha*p)

        ! ...  res = res - alpha * q 
        res = res - (alpha*q)

c       .....   check of the convergence   ...............
        rnorm = res * res

        if(icheck==1)  then
           WRITE(*,*) "myrank : ", myrank, "rnorm : ", rnorm
        endif

        if (rnorm < eps)  then
          idone = i
          go to 5000
        endif

        ! ...  s = W_adj * res 
        if(iflag==1)  then
           s = wxvect(res,3)
        else if(iflag==2)  then
           s = wxvect(res,2)
        endif

        c3 = s * s
c
        beta = c3 / c1
        c1 = c3

        p = s + (beta*p)

      enddo
      WRITE(*,*) "CG does not converge !!!"

 5000 continue
*     WRITE(*,*) "idone: ", idone

      icheck = 0

      return
      end

c-------------------------------------------------------------------------c
      FUNCTION wxvect(x,iflag)
c-------------------------------------------------------------------------c
c     wxvect = x        for iflag=1                                       c
c              Wx                 2                                       c
c              (W_adj)x           3                                       c
c-------------------------------------------------------------------------c
      USE field_f
      USE fpara
      TYPE(f_field)  wxvect
      TYPE(f_field), INTENT(IN) :: x
      TYPE(f_field) x5
      COMPLEX*16 fac1, fac2  ! added by A.N. on 2008/6/24

      temp  = zero_f

      if (iflag == 1) then           ! y = x

         wxvect = x

      else if (iflag == 2) then      ! y = Wx 

         call set_wing_f1(x) 

         do nu = 1,4

         temp1 =   nu  .fshift. x
         temp2 = (-nu) .fshift. x
         temp = temp + ((hopp(nu)*temp1) + (hopm(nu)*temp2)) 

         enddo

         wxvect =  x - temp

         IF( abs(Csw) > 0.001 ) THEN
           CALL vclover(wxvect,x)
         ENDIF
      
      else if (iflag == 3) then      ! y = (W_adj)x

         call g5xvect(x5,x)

         call set_wing_f1(x5)    ! Added on July 14,2000

         do nu = 1, 4

         temp1 =    nu .fshift. x5
         temp2 = (-nu) .fshift. x5

         ! --------------------------------------------------!
         ! Changed by A.N. 2008/6/24
         ! The following works only when the chemical potential is real.
         ! For general complex chemical potential,
         ! W(\mu)_adj = \gamma_5 W(-CONJG(\mu)) \gamma_5
*        temp  = temp + ((hopm(nu)*temp1)+ (hopp(nu)*temp2))

         if(nu.ne.4) then
            fac1 = hopp(nu)
            fac2 = hopm(nu)
         else
            fac1 = hop * exp(-CONJG(cmu))
            fac2 = hop * exp(+CONJG(cmu))
         endif

         temp  = temp + ((fac1*temp1)+ (fac2*temp2))
         ! --------------------------------------------------!

         enddo

         temp3 =  x5 - temp

         IF( abs(Csw) > 0.001 ) THEN
           CALL vclover(temp3,x5)
         ENDIF

         call g5xvect(wxvect,temp3)

      end if
      
      return
      end function

c--------------------------------------------------------------------------c
      subroutine g5xvect(y,x)
c--------------------------------------------------------------------------c
c     y = gamma_5 * x
c     here
c                  ( -1       )
c        GAMMA5 =  (   -1     )
c                  (     +1   )
c                  (       +1 )
c--------------------------------------------------------------------------c
      USE field_f
      USE fpara

      include '../INCLUDE/para_geometry' 

      INTEGER ncount(4)
      COMPLEX*16 c0(4), cx, cy
      common/ fmultest/ c0, ncount

      TYPE(f_field) x, y

      do ic = 1,NC
      do it = 1,NT
      do iz = 1,NZ
      do iy = 1,NY
      do ix = 1,NX

          y%f(ic,ix,iy,iz,it,1)= -x%f(ic,ix,iy,iz,it,1)
          y%f(ic,ix,iy,iz,it,2)= -x%f(ic,ix,iy,iz,it,2)
          y%f(ic,ix,iy,iz,it,3)= +x%f(ic,ix,iy,iz,it,3)
          y%f(ic,ix,iy,iz,it,4)= +x%f(ic,ix,iy,iz,it,4)

      enddo
      enddo
      enddo
      enddo
      enddo
*     WRITE(*,*) "y (g5xvect):",y%f(1,1,1,1,1,1),y*y

*     do mu = 1, 4
*     cx = 0.d0
*     cy = 0.d0
*     do k = 1, 3
*     do it = 1, 4
*     do iz = 1, 4
*     do iy = 1, 4
*     do ix = 1, 4

*     cx = cx + abs(x%f(k,ix,iy,iz,it,mu))**2
*     cy = cy + abs(y%f(k,ix,iy,iz,it,mu))**2

*     enddo
*     enddo
*     enddo
*     enddo
*     enddo
*     write(*,*) "mu, cx, cy: ", mu, cx, cy
*     enddo

*     write(*,*) "x*x : ", x*x
*     do mu = 1, 4
*     write(*,*) "nounct, c0 : ", ncount(mu), c0(mu)
*     enddo
*     write(*,*) "y*y : ", y*y
*     do mu = 1, 4
*     write(*,*) "nounct, c0 : ", ncount(mu), c0(mu)
*     enddo
   
*     STOP

      return
      end

c----------------------------------------------------------------------c
      subroutine mk_gamma
c----------------------------------------------------------------------c
c     Make gamma matrix
c----------------------------------------------------------------------c
C     THE CONVENTION OF THE GAMMA MATRIX HERE
C     ( EUCLIDEAN CHIRAL REPRESENTATION )
C
C               (       -i )              (       -1 )
C     GAMMA1 =  (     -i   )     GAMMA2 = (     +1   )
C               (   +i     )              (   +1     )
C               ( +i       )              ( -1       )
C
C               (     -i   )              (     -1   )
C     GAMMA3 =  (       +i )     GAMMA4 = (       -1 )
C               ( +i       )              ( -1       )
C               (   -i     )              (   -1     )
C
C               ( -1       )
C     GAMMA5 =  (   -1     )
C               (     +1   )
C               (       +1 )
C
C     ( GAMMA_MU, GAMMA_NU ) = 2*DEL_MU,NU   FOR MU,NU=1,2,3,4   
c----------------------------------------------------------------------c
      USE fpara
      ! fpara contains hop,r,BC,gamma,g5,rpg,rmg

      COMPLEX*16  g0(4,4), g1(4,4), g2(4,4), g3(4,4), g4(4,4)

      COMPLEX*16 ci
      data ci/ (0.0d0,1.0d0)/

      write(*,*) "Now we are in mkgamma"

      g0(1,1)=1.d0; g0(1,2)=0.d0; g0(1,3)=0.d0; g0(1,4)=0.d0
      g0(2,1)=0.d0; g0(2,2)=1.d0; g0(2,3)=0.d0; g0(2,4)=0.d0
      g0(3,1)=0.d0; g0(3,2)=0.d0; g0(3,3)=1.d0; g0(3,4)=0.d0
      g0(4,1)=0.d0; g0(4,2)=0.d0; g0(4,3)=0.d0; g0(4,4)=1.d0

      g1(1,1)=0.d0; g1(1,2)=0.d0; g1(1,3)=0.d0; g1(1,4)=-CI
      g1(2,1)=0.d0; g1(2,2)=0.d0; g1(2,3)=-CI;  g1(2,4)=0.d0
      g1(3,1)=0.d0; g1(3,2)=+CI;  g1(3,3)=0.d0; g1(3,4)=0.d0
      g1(4,1)=+CI;  g1(4,2)=0.d0; g1(4,3)=0.d0; g1(4,4)=0.d0

      g2(1,1)=0.d0; g2(1,2)=0.d0; g2(1,3)=0.d0; g2(1,4)=-1.d0
      g2(2,1)=0.d0; g2(2,2)=0.d0; g2(2,3)=1.d0; g2(2,4)=0.d0
      g2(3,1)=0.d0; g2(3,2)=1.d0; g2(3,3)=0.d0; g2(3,4)=0.d0
      g2(4,1)=-1.d0;g2(4,2)=0.d0; g2(4,3)=0.d0; g2(4,4)=0.d0

      g3(1,1)=0.d0; g3(1,2)=0.d0; g3(1,3)=-CI;  g3(1,4)=0.0
      g3(2,1)=0.d0; g3(2,2)=0.d0; g3(2,3)=0.d0; g3(2,4)=+CI
      g3(3,1)=+CI;  g3(3,2)=0.d0; g3(3,3)=0.d0; g3(3,4)=0.d0
      g3(4,1)=0.d0; g3(4,2)=-CI;  g3(4,3)=0.d0; g3(4,4)=0.d0

      g4(1,1)=0.d0; g4(1,2)=0.d0; g4(1,3)=-1.d0;g4(1,4)=0.d0
      g4(2,1)=0.d0; g4(2,2)=0.d0; g4(2,3)=0.d0; g4(2,4)=-1.d0
      g4(3,1)=-1.d0;g4(3,2)=0.d0; g4(3,3)=0.d0; g4(3,4)=0.d0
      g4(4,1)=0.d0; g4(4,2)=-1.d0;g4(4,3)=0.d0; g4(4,4)=0.d0

      g5(1,1)=-1.d0;g5(1,2)=0.d0; g5(1,3)=0.d0; g5(1,4)=0.d0
      g5(2,1)=0.d0; g5(2,2)=-1.d0;g5(2,3)=0.d0; g5(2,4)=0.d0
      g5(3,1)=0.d0; g5(3,2)=0.d0; g5(3,3)=1.d0; g5(3,4)=0.d0
      g5(4,1)=0.d0; g5(4,2)=0.d0; g5(4,3)=0.d0; g5(4,4)=1.d0

      do j = 1, 4
        do i = 1, 4
 
          gamma(i,j,1) = g1(i,j) 
          gamma(i,j,2) = g2(i,j)
          gamma(i,j,3) = g3(i,j)
          gamma(i,j,4) = g4(i,j)
          gamma(i,j,5) = g5(i,j)

        enddo
      enddo

      do mu = 1, 4
        do j = 1, 4
          do i = 1, 4
          
            rpg(i,j,mu) = r*g0(i,j) + gamma(i,j,mu)
            rmg(i,j,mu) = r*g0(i,j) - gamma(i,j,mu)

          enddo
        enddo
      enddo
 
      return
      end
