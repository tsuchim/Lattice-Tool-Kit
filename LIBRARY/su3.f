c su3.f
c**********************************************************************c
*     SU(3) Library                                                    * 
*     Most of programs here were originally written by Ph.dF           *
c----------------------------------------------------------------------c
      subroutine set_ident(a)
c----------------------------------------------------------------------c
c     Set all input matrices the identity matrix                       
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) a

      do 10  i = 1, NV/NBUSH
     
         a%g(1,1,i) = (1.0d0, 0.0d0)
         a%g(1,2,i) = (0.0d0, 0.0d0)
         a%g(1,3,i) = (0.0d0, 0.0d0)
         a%g(2,1,i) = (0.0d0, 0.0d0)
         a%g(2,2,i) = (1.0d0, 0.0d0)
         a%g(2,3,i) = (0.0d0, 0.0d0)
         a%g(3,1,i) = (0.0d0, 0.0d0)
         a%g(3,2,i) = (0.0d0, 0.0d0)
         a%g(3,3,i) = (1.0d0, 0.0d0)

   10 continue

      return
      end

c----------------------------------------------------------------------c
      subroutine set_zero(a)
c----------------------------------------------------------------------c
c     sets to zero NV/NBUSH matrices                                         c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) a

      do i = 1, NV/NBUSH

        a%g(1, 1, i) = (0.d0,0.d0)
        a%g(1, 2, i) = (0.d0,0.d0)
        a%g(1, 3, i) = (0.d0,0.d0)
        a%g(2, 1, i) = (0.d0,0.d0)
        a%g(2, 2, i) = (0.d0,0.d0)
        a%g(2, 3, i) = (0.d0,0.d0)
        a%g(3, 1, i) = (0.d0,0.d0)
        a%g(3, 2, i) = (0.d0,0.d0)
        a%g(3, 3, i) = (0.d0,0.d0)

      enddo

      return
      end

c----------------------------------------------------------------------c
      subroutine tracev(u,v,s,iflag)
c----------------------------------------------------------------------c
c     iflag = 0                                                        c
c         s(i) =  Re tr u(i) * v_adj(i)                                c
c     iflag = 1                                                        c
c         s(i) =  s(i) + Re tr u(i) * v_adj(i)                         c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) u, v

      REAL*8 s(NV/NBUSH)
      COMPLEX*16 uu(3,3), vv(3,3)

      if (iflag==0)  then

        do i = 1, NV/NBUSH

        uu(1,1) = u%g(1,1,i); uu(1,2) = u%g(1,2,i); uu(1,3) = u%g(1,3,i) 
        uu(2,1) = u%g(2,1,i); uu(2,2) = u%g(2,2,i); uu(2,3) = u%g(2,3,i) 
        uu(3,1) = u%g(3,1,i); uu(3,2) = u%g(3,2,i); uu(3,3) = u%g(3,3,i) 

        vv(1,1) = conjg( v%g(1,1,i) )
        vv(1,2) = conjg( v%g(1,2,i) )
        vv(1,3) = conjg( v%g(1,3,i) )
        vv(2,1) = conjg( v%g(2,1,i) )
        vv(2,2) = conjg( v%g(2,2,i) )
        vv(2,3) = conjg( v%g(2,3,i) )
        vv(3,1) = conjg( v%g(3,1,i) )
        vv(3,2) = conjg( v%g(3,2,i) )
        vv(3,3) = conjg( v%g(3,3,i) )

        s(i) = uu(1,1)*vv(1,1) + uu(1,2)*vv(1,2) + uu(1,3)*vv(1,3)
     &       + uu(2,1)*vv(2,1) + uu(2,2)*vv(2,2) + uu(2,3)*vv(2,3)
     &       + uu(3,1)*vv(3,1) + uu(3,2)*vv(3,2) + uu(3,3)*vv(3,3)

        enddo

      else

        do i = 1, NV/NBUSH
         
        uu(1,1) = u%g(1,1,i); uu(1,2) = u%g(1,2,i); uu(1,3) = u%g(1,3,i) 
        uu(2,1) = u%g(2,1,i); uu(2,2) = u%g(2,2,i); uu(2,3) = u%g(2,3,i) 
        uu(3,1) = u%g(3,1,i); uu(3,2) = u%g(3,2,i); uu(3,3) = u%g(3,3,i) 

        vv(1,1) = conjg( v%g(1,1,i) )
        vv(1,2) = conjg( v%g(1,2,i) )
        vv(1,3) = conjg( v%g(1,3,i) )
        vv(2,1) = conjg( v%g(2,1,i) )
        vv(2,2) = conjg( v%g(2,2,i) )
        vv(2,3) = conjg( v%g(2,3,i) )
        vv(3,1) = conjg( v%g(3,1,i) )
        vv(3,2) = conjg( v%g(3,2,i) )
        vv(3,3) = conjg( v%g(3,3,i) )

        s(i) = s(i)
     &       + uu(1,1)*vv(1,1) + uu(1,2)*vv(1,2) + uu(1,3)*vv(1,3)
     &       + uu(2,1)*vv(2,1) + uu(2,2)*vv(2,2) + uu(2,3)*vv(2,3)
     &       + uu(3,1)*vv(3,1) + uu(3,2)*vv(3,2) + uu(3,3)*vv(3,3)


        enddo

      endif

      return
      end

c----------------------------------------------------------------------c
      subroutine set_rand(a)
c----------------------------------------------------------------------c
c     sets randomly NV/NBUSH matrices
c                                       written by AN
c                                       corrected by SM  9 Mar 2000
c                                       ranf is changed to ran3 by SM   
c
c                                       revised by SC 17 Apr 2000
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) a

      do 100 i = 1, NV/NBUSH

         a%g(1,1,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)
         a%g(1,2,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)
         a%g(1,3,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)
         a%g(2,1,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)
         a%g(2,2,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)
         a%g(2,3,i) = cmplx (ranf(idum)-.5d0, ranf(idum)-.5d0)

 100  continue
 
      call normalize(a)

      return
      end

c----------------------------------------------------------------------c
      subroutine normalize(u)
c----------------------------------------------------------------------c
c     normalizes links                                                 c
c     input  x : arbitrary 3*3 complex matrix                          c
c     output x : SU(3) matrix                                          c
c                                        Feb.29, 2000 written by AN    c
c                                        Mar.07, 2000 checked by SC    c
c                                        Apr.07, 2000 changed by SM    c
c                                        Apr.13, 2000 checked by SC    c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) u
    
      COMPLEX*16  x4, x5, x6, w1, w2
!     w2 is real ??????????

      do 100 i = 1, NV/NBUSH   

        w1 =    u%g(2,1,i)*conjg(u%g(1,1,i) ) 
     &       +  u%g(2,2,i)*conjg(u%g(1,2,i) )
     &       +  u%g(2,3,i)*conjg(u%g(1,3,i) )

        w2 =    u%g(1,1,i)*conjg(u%g(1,1,i) ) 
     &       +  u%g(1,2,i)*conjg(u%g(1,2,i) )
     &       +  u%g(1,3,i)*conjg(u%g(1,3,i) )

        zerock2 = w2
        if(zerock2.eq.0.)  then
           write(*,*) 'w2 is zero  !!  (in normlz)'
           write(*,*) 'u%g(1,1),u%g(1,2),u%g(1,3) : ', 
     &                u%g(1,1,i), u%g(1,2,i), u%g(1,3,i)
           stop 
 1      endif

        w1 = -w1/w2

        x4 = (u%g(2,1,i)) + w1*u%g(1,1,i)
        x5 = (u%g(2,2,i)) + w1*u%g(1,2,i)
        x6 = (u%g(2,3,i)) + w1*u%g(1,3,i)

        w3 = x4*conjg(x4) + x5*conjg(x5) + x6*conjg(x6)

        zerock3 = w3
        if(zerock 3.eq.0.)  then
           write(*,*) 'w3 is zero  !!  (in normlz)'
           write(*,*) 'x4, x5, x6 : ', x4, x5, x6
           stop
        endif

        u%g(2,1,i) = x4
        u%g(2,2,i) = x5
        u%g(2,3,i) = x6

        w3 = 1.d0/sqrt(w3)
        w2 = 1.d0/sqrt(w2)

        u%g(1,1,i) = u%g(1,1,i)*w2
        u%g(1,2,i) = u%g(1,2,i)*w2
        u%g(1,3,i) = u%g(1,3,i)*w2
        u%g(2,1,i) = u%g(2,1,i)*w3
        u%g(2,2,i) = u%g(2,2,i)*w3
        u%g(2,3,i) = u%g(2,3,i)*w3

        if((zerock2*zerock3).eq.0.)  then
           write(*,*) '!! devided by zero !! (in normalize)'
           write(*,*) 'w2 or w3 in normlz is zero !!'
           write(*,*) 'w2, w3 : ', w2, w3   
           stop
        endif 

100   continue

      call m3complv(u)

      return
      end
 
c----------------------------------------------------------------------c
      subroutine m3complv(a)
c----------------------------------------------------------------------c
c     Constructs the third row of su(3) matrix                         c
c                                        Feb.29, 2000 written by AN    c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) a
      REAL*8 aa(18)

      do 100 i = 1, NV/NBUSH

         aa( 1) = REAL( a%g(1,1,i))
         aa( 2) = AIMAG(a%g(1,1,i))
         aa( 3) = REAL( a%g(1,2,i))
         aa( 4) = AIMAG(a%g(1,2,i))
         aa( 5) = REAL( a%g(1,3,i))
         aa( 6) = AIMAG(a%g(1,3,i))
         aa( 7) = REAL( a%g(2,1,i))
         aa( 8) = AIMAG(a%g(2,1,i))
         aa( 9) = REAL( a%g(2,2,i))
         aa(10) = AIMAG(a%g(2,2,i))
         aa(11) = REAL( a%g(2,3,i))
         aa(12) = AIMAG(a%g(2,3,i))

         aa(13) = aa( 3)*aa(11) - aa( 4)*aa(12)
     &          - aa( 5)*aa( 9) + aa( 6)*aa(10)
         aa(14) = aa( 5)*aa(10) + aa( 6)*aa( 9)
     &          - aa( 3)*aa(12) - aa( 4)*aa(11)
         aa(15) = aa( 5)*aa( 7) - aa( 6)*aa( 8)
     &          - aa( 1)*aa(11) + aa( 2)*aa(12)
         aa(16) = aa( 1)*aa(12) + aa( 2)*aa(11)
     &          - aa( 5)*aa( 8) - aa( 6)*aa( 7)
         aa(17) = aa( 1)*aa( 9) - aa( 2)*aa(10)
     &          - aa( 3)*aa( 7) + aa( 4)*aa( 8)
         aa(18) = aa( 3)*aa( 8) + aa( 4)*aa( 7)
     &          - aa( 1)*aa(10) - aa( 2)*aa( 9)

         a%g(3,1,i) = cmplx( aa(13), aa(14) )
         a%g(3,2,i) = cmplx( aa(15), aa(16) )
         a%g(3,3,i) = cmplx( aa(17), aa(18) )

 100  continue

      return
      end

c----------------------------------------------------------------------c
      subroutine  set_rand2 (ur,dia,del)
c----------------------------------------------------------------------c
c     random su(3) matrix around the unit matrix                       c
c     the matrices ur distribute around the unit matrix.               c
c     the shape of this distribution is controled by the               c
c     parameters 'dia' and 'del'                                       c
c     ( when 'dia' is very large or 'del' is very small,               c
c       ur is near to the unit matrix.  )                              c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(g_field1) ur
      REAL*8 dia, del

      COMPLEX*16  v1(3), v2(3), v3(3), v12, v13, v23, det
      parameter( NRAN=8192 )
      REAL*8  rg(2,3,3,NRAN)

      if(mod((NV/NBUSH),2).ne.0)  then
         write(*,*) 
     &   "Sorry NV/NBUSH must be an even number. NV/NBUSH=",
     &    NV/NBUSH
         write(*,*) "(set_rand2)"
         stop
      endif
      nmaxh=(NV/NBUSH)/2

      if(nmaxh>NRAN)  then
         write(*,*) "Sorry, NRAN should be greater than nmaxh"
         write(*,*) "NRAN, nmaxh : ", NRAN, nmaxh
         write(*,*) "(set_rand2)"
         stop
      endif

c     .....  bring random numbers rg.  this is gaussian distrib.  ...
      call  rgauss1(rg,2*3*3*nmaxh)
      do 1000  i = 1, nmaxh

c     .....   generate three ramdom vecors   ...........
      do ic = 1, 3
         v1(ic) = cmplx( del*rg(1,1,ic,i), del*rg(2,1,ic,i) )
         v2(ic) = cmplx( del*rg(1,2,ic,i), del*rg(2,2,ic,i) )
         v3(ic) = cmplx( del*rg(1,3,ic,i), del*rg(2,3,ic,i) )
      enddo

      v1(1) = v1(1) + dia
      v2(2) = v2(2) + dia
      v3(3) = v3(3) + dia

c     .....   schimidt orthoganalization   .............

      vnorm = 0.d0

      do ic = 1, 3
         vnorm = vnorm + conjg(v1(ic))*v1(ic)
      enddo

      vnorm = sqrt(vnorm)

      do ic = 1, 3
         v1(ic) = v1(ic)/vnorm
      enddo

      v12 = 0.0
      do ic = 1, 3
         v12 = v12 + conjg(v1(ic))*v2(ic)
      enddo

      do ic = 1, 3
         v2(ic) = v2(ic) - v1(ic)*v12
      enddo

      vnorm = 0.0
      do ic = 1, 3
         vnorm = vnorm + conjg(v2(ic))*v2(ic)
      enddo

      vnorm = sqrt(vnorm)
      do ic = 1, 3
         v2(ic) = v2(ic)/vnorm
      enddo

      v13 = 0.d0
      do ic = 1, 3
         v13 = v13 + conjg(v1(ic))*v3(ic)
      enddo

      v23 = 0.d0
      do ic = 1, 3
         v23 = v23 + conjg(v2(ic))*v3(ic)
      enddo

      do ic = 1, 3
         v3(ic) = v3(ic) - v1(ic)*v13 - v2(ic)*v23
      enddo

      det = v1(1) * ( v2(2)*v3(3) - v3(2)*v2(3) )
     *    - v1(2) * ( v2(1)*v3(3) - v3(1)*v2(3) )
     *    + v1(3) * ( v2(1)*v3(2) - v3(1)*v2(2) )

      do ic = 1, 3
         v3(ic) = v3(ic)/det
      enddo

      do ic = 1, 3
         ur%g(ic,1,i) = v1(ic)
         ur%g(ic,2,i) = v2(ic)
         ur%g(ic,3,i) = v3(ic)

         ur%g(1,ic,i+nmaxh) = conjg(v1(ic))
         ur%g(2,ic,i+nmaxh) = conjg(v2(ic))
         ur%g(3,ic,i+nmaxh) = conjg(v3(ic))
      enddo

 1000 continue

      return
      end

c-----------------------------------------------------c
      subroutine gaussian(nv,granf,variance)
c-----------------------------------------------------c
c     generates nv random numbers with normal 
c     distribution by Box-Muller  
c-----------------------------------------------------c
      real*8 PI, PI2
      parameter( PI=3.141592653589793238d0, PI2=PI*2.d0) ! PI2=PI*2
      real*8 granf(nv), variance
      real*8 rho, theta

      nvh = nv / 2
      do i = 1, nvh
        rho = sqrt(-2. * alog(ranf()) * variance)
        theta = PI2 * ranf()
        granf(i      ) = rho * cos(theta)
        granf(i + nvh) = rho * sin(theta)
      enddo

      if (2 * nvh .eq. nv) return
      granf(nv) = sqrt(-2. * alog(ranf()) * variance) 
     &          * cos(PI2 * ranf())

      return
      end

c----------------------------------------------------------------------c
      subroutine  grotate (u,g)
c----------------------------------------------------------------------c
c     GAUGE TRANSFORMATION                                             C
c     U : GAUGE FIELDS    G : GAUGE ROTATION MATRIX                    C
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0)  u(4), g
      TYPE(g_field1)  temp1, temp2, temp3

      call set_wing_g2(g)  

      do  mu  = 1, 4

         do ibush = 0, NBUSH-1

*          .....   CALCULATE  U(IS)*G(IS+MU)   .....
           temp1 = mu.gshift.g
           temp2 = u(mu)
           temp3 = temp2 * temp1  ! temp3 = U(is)*G(is+mu)
   
*          .....   CALCULATE  GD(IS) * (U(IS)*G(IS+MU))           .....
*          .....              WHERE GD : COMPLEX CONJUGATE OF G   .....

           temp1 = g
           temp2 = temp1 .ADprod. temp3 
           u(mu) = temp2

         enddo 
      enddo

      do mu = 1, 4
        u(mu)%direction = mu
      enddo

      return
      end

c----------------------------------------------------------------------c
      SUBROUTINE LambdaMul(a,b,k)
c----------------------------------------------------------------------c
c     b = (lambda_k/2)*a
C             lambda_k : GellMann matrices. k=1, 8 
c----------------------------------------------------------------------c
      USE field_g
      implicit none
      include '../INCLUDE/para_geometry'

      INTEGER,        INTENT(IN)  :: k
      TYPE(g_field1), INTENT(IN)  :: a
      TYPE(g_field1), INTENT(OUT) :: b

      INTEGER i
      COMPLEX*16 ci, ci2
      PARAMETER( ci=(0.d0,1.d0), ci2=ci/2.d0)
      REAL*8  sqr3, sqr3inv, sr3ih  !  SQRT(3), 1/SQRT(3), 0.5/SQRT(3)
      PARAMETER( sqr3=1.73205080756887719d0, sqr3inv=1.d0/sqr3 ) 
      PARAMETER( sr3ih = 0.5d0*sqr3inv )
      
      IF( k == 1 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = 0.5d0 * a%g(2,1,i) 
            b%g(1,2,i) = 0.5d0 * a%g(2,2,i)
            b%g(1,3,i) = 0.5d0 * a%g(2,3,i)
            b%g(2,1,i) = 0.5d0 * a%g(1,1,i)
            b%g(2,2,i) = 0.5d0 * a%g(1,2,i)
            b%g(2,3,i) = 0.5d0 * a%g(1,3,i)
            b%g(3,1,i) = (0.0d0, 0.0d0)
            b%g(3,2,i) = (0.0d0, 0.0d0)
            b%g(3,3,i) = (0.0d0, 0.0d0)

         ENDDO

      ENDIF

      IF( k == 2 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = -ci2 * a%g(2,1,i) 
            b%g(1,2,i) = -ci2 * a%g(2,2,i)
            b%g(1,3,i) = -ci2 * a%g(2,3,i)
            b%g(2,1,i) =  ci2 * a%g(1,1,i)
            b%g(2,2,i) =  ci2 * a%g(1,2,i)
            b%g(2,3,i) =  ci2 * a%g(1,3,i)
            b%g(3,1,i) = (0.0d0, 0.0d0)
            b%g(3,2,i) = (0.0d0, 0.0d0)
            b%g(3,3,i) = (0.0d0, 0.0d0)

         ENDDO

      ENDIF

      IF( k == 3 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) =  0.5d0 * a%g(1,1,i) 
            b%g(1,2,i) =  0.5d0 * a%g(1,2,i)
            b%g(1,3,i) =  0.5d0 * a%g(1,3,i)
            b%g(2,1,i) = -0.5d0 * a%g(2,1,i)
            b%g(2,2,i) = -0.5d0 * a%g(2,2,i)
            b%g(2,3,i) = -0.5d0 * a%g(2,3,i)
            b%g(3,1,i) = (0.0d0, 0.0d0)
            b%g(3,2,i) = (0.0d0, 0.0d0)
            b%g(3,3,i) = (0.0d0, 0.0d0)

         ENDDO

      ENDIF

      IF( k == 4 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = 0.5d0 * a%g(3,1,i) 
            b%g(1,2,i) = 0.5d0 * a%g(3,2,i)
            b%g(1,3,i) = 0.5d0 * a%g(3,3,i)
            b%g(2,1,i) = (0.0d0, 0.0d0)
            b%g(2,2,i) = (0.0d0, 0.0d0)
            b%g(2,3,i) = (0.0d0, 0.0d0)
            b%g(3,1,i) = 0.5d0 * a%g(1,1,i)
            b%g(3,2,i) = 0.5d0 * a%g(1,2,i)
            b%g(3,3,i) = 0.5d0 * a%g(1,3,i)

         ENDDO

      ENDIF

      IF( k == 5 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = -ci2 * a%g(3,1,i) 
            b%g(1,2,i) = -ci2 * a%g(3,2,i)
            b%g(1,3,i) = -ci2 * a%g(3,3,i)
            b%g(2,1,i) = (0.0d0, 0.0d0)
            b%g(2,2,i) = (0.0d0, 0.0d0)
            b%g(2,3,i) = (0.0d0, 0.0d0)
            b%g(3,1,i) =  ci2 * a%g(1,1,i)
            b%g(3,2,i) =  ci2 * a%g(1,2,i)
            b%g(3,3,i) =  ci2 * a%g(1,3,i)

         ENDDO

      ENDIF

      IF( k == 6 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = (0.0d0, 0.0d0)
            b%g(1,2,i) = (0.0d0, 0.0d0)
            b%g(1,3,i) = (0.0d0, 0.0d0)
            b%g(2,1,i) = 0.5d0 * a%g(3,1,i) 
            b%g(2,2,i) = 0.5d0 * a%g(3,2,i)
            b%g(2,3,i) = 0.5d0 * a%g(3,3,i)
            b%g(3,1,i) = 0.5d0 * a%g(2,1,i)
            b%g(3,2,i) = 0.5d0 * a%g(2,2,i)
            b%g(3,3,i) = 0.5d0 * a%g(2,3,i)

         ENDDO

      ENDIF

      IF( k == 7 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) = (0.0d0, 0.0d0)
            b%g(1,2,i) = (0.0d0, 0.0d0)
            b%g(1,3,i) = (0.0d0, 0.0d0)
            b%g(2,1,i) = -ci2 * a%g(3,1,i) 
            b%g(2,2,i) = -ci2 * a%g(3,2,i)
            b%g(2,3,i) = -ci2 * a%g(3,3,i)
            b%g(3,1,i) =  ci2 * a%g(2,1,i)
            b%g(3,2,i) =  ci2 * a%g(2,2,i)
            b%g(3,3,i) =  ci2 * a%g(2,3,i)

         ENDDO

      ENDIF

      IF( k == 8 ) THEN

         DO  i = 1, NV/NBUSH
     
            b%g(1,1,i) =  sr3ih * a%g(1,1,i) 
            b%g(1,2,i) =  sr3ih * a%g(1,2,i)
            b%g(1,3,i) =  sr3ih * a%g(1,3,i)
            b%g(2,1,i) =  sr3ih * a%g(2,1,i) 
            b%g(2,2,i) =  sr3ih * a%g(2,2,i)
            b%g(2,3,i) =  sr3ih * a%g(2,3,i)
            b%g(3,1,i) = -sqr3inv * a%g(3,1,i)
            b%g(3,2,i) = -sqr3inv * a%g(3,2,i)
            b%g(3,3,i) = -sqr3inv * a%g(3,3,i)

         ENDDO

      ENDIF

      return
      end
