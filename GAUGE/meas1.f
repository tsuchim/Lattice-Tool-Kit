c-q---------------------------------------------------------------------c
      subroutine  meas1(u,plaq)
c----------------------------------------------------------------------c
c
*     Measure the plaquete energy.                                     *
c
c     This is still temporal version written by SM 11 Apr 2000    
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u(4)
      TYPE(g_field1) staple, temp

      REAL*8 trace(NV/NBUSH) 
      REAL*8 plaq

      plaq = 0.0d0

      do mu = 1, 4
         call set_wing_g2(u(mu))
      enddo

      EvenOdd: do ibush = 0, NBUSH-1

         Direction: do mu = 1, 4

            call make_staple1(u,staple,mu)
            temp = u(mu) 
            call tracev( temp,staple,trace,0)

          Site : do i = 1, NV/NBUSH

            plaq = plaq + trace(i)

          enddo Site

         enddo Direction

       enddo EvenOdd
 
      plaq = plaq / (NV*4*3*3) 
                  !site*mu*nu*color

      return
      end
c----------------------------------------------------------------------c 
*     subroutine  make_staple1(u,staple,mu,ieo) ! changed by AN on 00.12.2
      subroutine  make_staple1(u,staple,mu)
c----------------------------------------------------------------------c
c     Make one staple
c----------------------------------------------------------------------c
      USE field_g
      implicit none

      TYPE(g_field0), INTENT(IN)::   u(4)
      TYPE(g_field1), INTENT(OUT)::  staple
      INTEGER, INTENT(IN)::          mu

      TYPE(g_field1)  temp1, temp2, temp3
      TYPE(ivec2)    idir2
      INTEGER nu
            
      call set_zero(staple)
      
      do nu = 1, 4
      if(nu==mu)  cycle
      
      
c       x+nu temp2
c        .---------.
c        I         I
c  temp1 I         I
c        I         I
c        .         .
c        x        x+mu

         temp1 = u(nu)
         temp2 = nu.gshift.u(mu)
         temp3 = temp1 * temp2
         temp1 = mu.gshift.u(nu)
         staple = staple  + (temp3.prodAD.temp1)

      enddo
c      
      return
      end

c--------------------------------------------------------------------c
      subroutine  make_istaple1(u,staple,mu)
c--------------------------------------------------------------------c
c
c     Make staples for improved action
c          (here 1x1 + 1x2) 
c     c0 : coefficient of 1x1
c     c1 : coefficient of 1x2
c      PARAMETER ( c0=6.1564, c1=(-0.6241) )
c
c     Only the positive direction, i.e., type 1, 3, 5 and 7
c------------------------------------------------------------------c
c
c    +------+                +------+------+  
c    |      |                |             |  
c    |      |                |             |  
c    +      +    +      +    +      +------+   +       ------+  
c                |      |                      |             |  
c                |      |                      |             |  
c                +------+                      +------+------+  
c
c     type 1      type 2         type 3           type 4
c
c                                          +------+ 
c                                          |      | 
c                                          |      | 
c    +------+------+                       +      + 
c    |             |                       |      | 
c    |             |                       |      | 
c    +------+      +   +------+      +     +      +     +      +
c                      |             |                  |      | 
c                      |             |                  |      | 
c                      +------+------+                  +      + 
c                                                       |      |
c        type 5             type 6                      |      |
c                                                       +------+
c                                            type 7      type 8
c
c   This classification is due to Prof.Sakai
c-----------------------------------------------------------------c
      USE field_g
      TYPE(g_field0) u(4)
      TYPE(g_field1)  staple, temp1, temp2, temp3
      TYPE(ivec2)    idir2
      type(ivec3)    idir3
c
      REAL*8  c0, c1
*     PARAMETER ( c0=1.0, c1=0.0 )
      PARAMETER ( c0=3.648, c1=(-0.331) )
c
      call set_zero(staple)
c
c
c-------------------------------------------------------------------c
c     ordinary Wilson part
c-------------------------------------------------------------------c
c

      do nu = 1, 4
      if(nu==mu)  cycle
      
      
c       x+nu temp2
c        .---------.
c        I         I
c  temp1 I         I
c        I         I
c        .         .
c        x        x+mu

         temp1 = u(nu)
         temp2 = nu.gshift.u(mu)
         temp3 = temp1 * temp2
         temp1 = mu.gshift.u(nu)
         temp3 = temp3.prodAD.temp1

         staple = staple + (c0*temp3)

c-------------------------------------------------------------------c
c
c     additional term for improved action
c
c-------------------------------------------------------------------c
c      x+nu
c       +------+------+ 
c       |             | 
c       |             | 
c       +      +------+    clockwise manner 
c       x     x+mu
c    type 3      

         temp1 = u(nu)
         temp2 = nu.gshift.u(mu)
         temp3 = temp1 * temp2

         idir2 = ivec2(mu,nu)
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3*temp1

         idir2 = ivec2(mu,mu)
         temp1 = idir2.gshift.u(nu)
         temp3 = temp3.prodAD.temp1

         temp1 = mu.gshift.u(mu)
         temp3 = temp3.prodAD.temp1

         staple = staple + (c1*temp3)

c-------------------------------------------------------------------c
c
c    x-mu+nu        x+mu-nu                  The node x-mu-nu must be 
c       +------+------+                      replaced by x-mu+nu. The
c       |             |                      replacement was done by
c  temp2|temp1        |                      P.I. on 3 of November.
c       +------+      +     clockwise        
c      x-mu     x     x+mu
c    type 5     

         temp1 = (-mu).gshift.u(mu)
         temp2 = (-mu).gshift.u(nu)
         temp3 = temp1.ADprod.temp2

         idir2 = ivec2(-mu,nu)
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3*temp1

         temp1 = nu.gshift.u(mu)
         temp3 = temp3*temp1

         temp1 = mu.gshift.u(nu)
         temp3 = temp3.prodAD.temp1

         staple = staple + (c1*temp3)

c-------------------------------------------------------------------c
c     x+2nu   x+nu+nu+mu         
c       +------+ 
c temp2 |      | 
c       |      | 
c   x+nu+      + x+mu+nu
c       |      | 
c temp1 |      | 
c       +      + 
c       x      x+mu       clockwise
c
c  type 7  
c
         temp1 = u(nu)
         temp2 = nu.gshift.u(nu)
         temp3 = temp1 * temp2

         idir2 = ivec2(nu,nu)
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3*temp1

         idir2 = ivec2(nu,mu)
         temp1 = idir2.gshift.u(nu)
         temp3 = temp3.prodAD.temp1

         temp1 = mu.gshift.u(nu)
         temp3 = temp3.prodAD.temp1

         staple = staple + (c1*temp3)

      end do

      return
      end

c----------------------------------------------------------------------c
      subroutine  meas1ST(u,plaqS,plaqT)
c----------------------------------------------------------------------c
c
*     Measure the Spatial and Temporal  plaquete energy.               *
c
c     This is still temporal version written by SM 11 Apr 2000    
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      REAL*8, intent(out) :: plaqS, plaqT
      TYPE(g_field0), intent(in) ::  u(4)
      TYPE(g_field1) staple, temp

      REAL*8 trace(NV/NBUSH) 

      plaqS = 0.0d0
      plaqT = 0.0d0

      do mu = 1, 4
         call set_wing_g2(u(mu))
      enddo

      EvenOdd: do ibush = 0, NBUSH-1

         DirectionMu: do mu = 1, 3
         DirectionNu: do nu = mu+1, 4

            call make_staple12(u,staple,mu,nu)
            temp = u(mu) 
            call tracev( temp,staple,trace,0)

            IF( nu==4 ) THEN
              SiteT : do i = 1, NV/NBUSH
                plaqT = plaqT + trace(i)
              enddo SiteT
            ELSE
              SiteS : do i = 1, NV/NBUSH
                plaqS = plaqS + trace(i)
              enddo SiteS
            ENDIF

         enddo DirectionNu
         enddo DirectionMu

       enddo EvenOdd
 
      plaqS = plaqS / (NV*6*3) 
      plaqT = plaqT / (NV*6*3) 
                  !site*mu*nu*color

      return
      end
c----------------------------------------------------------------------c 
      subroutine  make_staple12(u,staple,mu,nu)
c----------------------------------------------------------------------c
c     Make one staple
c----------------------------------------------------------------------c
      USE field_g
      implicit none

      TYPE(g_field0), INTENT(IN)::   u(4)
      TYPE(g_field1), INTENT(OUT)::  staple
      INTEGER, INTENT(IN)::          mu, nu

      TYPE(g_field1)  temp1, temp2, temp3
      TYPE(ivec2)    idir2
            
      call set_zero(staple)
      
      if(nu==mu)  return
      
      
c       x+nu temp2
c        .---------.
c        I         I
c  temp1 I         I
c        I         I
c        .         .
c        x        x+mu

         temp1 = u(nu)
         temp2 = nu.gshift.u(mu)
         temp3 = temp1 * temp2
         temp1 = mu.gshift.u(nu)
         staple = staple  + (temp3.prodAD.temp1)

      return
      end
