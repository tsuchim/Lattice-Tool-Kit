c----------------------------------------------------------------------c 
      subroutine  make_staple(u,staple,mu)
c----------------------------------------------------------------------c
c     Make staple
c----------------------------------------------------------------------c
      USE field_g
      USE aniso
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u(4)
      
      TYPE(g_field1)  staple, temp1, temp2, temp3
      TYPE(ivec2)    idir2
            
      real*8 fac

      call set_zero(staple)
      
      do nu = 1, 4
      if(nu==mu)  cycle
      if( (mu==4) .or. (nu==4) )  then
         fac = gamma_G        ! Temporal 
      else
         fac = 1.d0/gamma_G   ! Spatial
      endif
      
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
         staple = staple + (fac*(temp3.prodAD.temp1))

c        x
c        .         .
c        I         I
c  temp1 I         I
c        I         I
c        .---------.
c      x-nu        x+mu-nu
      
         idir2 = ivec2(mu,-nu)
         temp1 = (-nu).gshift.u(nu)
         temp2 = (-nu).gshift.u(mu)
         temp3 = temp1.ADprod.temp2
         temp1 = idir2.gshift.u(nu)
         staple = staple + (fac*(temp3*temp1))

      end do
      
      return
      end
      
c--------------------------------------------------------------------c
      subroutine  make_istaple(u,staple,mu)
c--------------------------------------------------------------------c
c
c     Make staples for improved action
c          (here 1x1 + 1x2) 
c     c0 : coefficient of 1x1
c     c1 : coefficient of 1x2
c      PARAMETER ( c0=6.1564, c1=(-0.6241) )
c
c     But I do not know how to normalize. 
c                                    by S. Muroya 5 June 2000
c
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

c        x
c        .         .
c        I         I
c  temp1 I         I
c        I         I
c        .---------.
c      x-nu        x+mu-nu
      
         idir2 = ivec2(mu,-nu)
         temp1 = (-nu).gshift.u(nu)
         temp2 = (-nu).gshift.u(mu)
         temp3 = temp1.ADprod.temp2
         temp1 = idir2.gshift.u(nu)
         temp3 = temp3*temp1

         staple = staple + (c0*temp3)
c
c      end do ! switch off improvement here
c
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
c       x      x+mu
c       +      +------+  
c       |             |  
c temp1 |temp2        |  
c       +------+------+   counter clockwise 
c      x+(-nu)
c    type 4

         temp1 = (-nu).gshift.u(nu)
         temp2 = (-nu).gshift.u(mu)
         temp3 = temp1.ADprod.temp2

*        idir2 = ivec2(mu,nu)
         idir2 = ivec2(mu,-nu)      ! corrected on Sept.3,2000
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3*temp1

         idir3 = ivec3(mu,mu,-nu)
         temp1 = idir3.gshift.u(nu)
         temp3 = temp3*temp1

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
c   x-mu temp2 x     x+mu
c       +------+      + 
c  temp1|             | 
c       |             | 
c       +------+------+     counter clockwise
c          
c     x-mu-nu       x+mu-nu
c   type 6     
c
         idir2 = ivec2(-mu,-nu)
         temp1 = idir2.gshift.u(nu)
         temp2 = (-mu).gshift.u(mu)
         temp3 = temp1*temp2

         idir2 = ivec2(-mu,-nu)
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3.ADprod.temp1

         temp1 = (-nu).gshift.u(mu)
         temp3 = temp3*temp1

         idir2 = ivec2(mu,-nu)
         temp1 = idir2.gshift.u(nu)
         temp3 = temp3*temp1

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

c-------------------------------------------------------------------c
c       x     x+mu
c       +      +
c       |      |
c temp2 |      |
c  x-nu +      + x-nu+mu
c       |      |
c temp1 |      |
c       +------+          counter clockwise
c   x-nu-nu   x-nu-nu+mu
c
c   type 8

         idir2 = ivec2(-nu,-nu)
         temp1 = idir2.gshift.u(nu)
         temp2 = (-nu).gshift.u(nu)
         temp3 = temp1*temp2

         idir2 = ivec2(-nu,-nu)
         temp1 = idir2.gshift.u(mu)
         temp3 = temp3.ADprod.temp1      

         idir3 = ivec3(+mu,-nu,-nu)
         temp1 = idir3.gshift.u(nu)
         temp3 = temp3*temp1

         idir2 = ivec2(+mu,-nu)
         temp1 = idir2.gshift.u(nu)
         temp3 = temp3*temp1

         staple = staple + (c1*temp3)

      end do

      return
      end
