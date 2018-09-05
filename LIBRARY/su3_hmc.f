* su3_hmc.f
c**********************************************************************c
*     SU(3) Library  for hmc
c-----------------------------------------------------c
      subroutine  smullink(fac,vin,vout)
c-----------------------------------------------------c
c     vout = fac * vin
c-----------------------------------------------------c
      include '../INCLUDE/para_geometry'
      real*8  fac, vin(8*NV), vout(8*NV)

      do i = 1,8*NV
	vout(i) = fac*vin(i)
      enddo	

      return
      end 

c-----------------------------------------------------c
      subroutine  projlink(vin,vout,nv)
c-----------------------------------------------------c
c     !!!!!   vin and vout should be different vectors
c
c     Projectin of the etraceless antiermite part 
c     vout = x/2 - Tr(x)/6
c     wher   x = vin - Conjg(vin)      
c-----------------------------------------------------c
      USE field_g

      TYPE(g_field1)  vin, vout 
      COMPLEX*16  x12, x13, x23, x21, x31, x32
      COMPLEX*16  v11, v12, v13, v21, v22, v23, v31, v32, v33
      REAL*8  tri, fac13

      fac13 = 1.d0/3.d0

      do i = 1, nv

         v11 = vin%g(1,1,i)
         v22 = vin%g(2,2,i)
         v33 = vin%g(3,3,i)
         tri = fac13 * ( AIMAG(v11) + AIMAG(v22) + AIMAG(v33) )

         vout%g(1,1,i) = CMPLX( 0.d0, AIMAG(v11)-tri ) 
         vout%g(2,2,i) = CMPLX( 0.d0, AIMAG(v22)-tri ) 
         vout%g(3,3,i) = CMPLX( 0.d0, AIMAG(v33)-tri ) 
 
      enddo

      do i = 1, nv

         v12 = vin%g(1,2,i)
         v13 = vin%g(1,3,i)
         v21 = vin%g(2,1,i)
         v23 = vin%g(2,3,i)
         v31 = vin%g(3,1,i)
         v32 = vin%g(3,2,i)

         x12 = v12 - CONJG(v21)
         x13 = v13 - CONJG(v31)
         x23 = v23 - CONJG(v32)
      
         x21 = - CONJG(x12)
         x31 = - CONJG(x13)
         x32 = - CONJG(x23)

         vout%g(1,2,i) = 0.5d0 * x12
         vout%g(1,3,i) = 0.5d0 * x13
         vout%g(2,1,i) = 0.5d0 * x21
         vout%g(2,3,i) = 0.5d0 * x23
         vout%g(3,1,i) = 0.5d0 * x31
         vout%g(3,2,i) = 0.5d0 * x32

      enddo

      return
      end

c-----------------------------------------------------c
      subroutine  algblink(x,c,nv)
c-----------------------------------------------------c
c     Expand 3x3 anti-Hermite Matrix x as
c       x = i * Sum c(l) * Lambda(l) / 2  
c     where Lambda's are Gel-man matrices for SU(3) 
c     algebra  i.e.,
c       i * c(l) = Tr(Lambda(l)*x) 
c-----------------------------------------------------c
      USE field_g
      USE hmc_mod

      REAL*8 sr3, sr3i
      parameter ( sr3=1.73205080756888d0, sr3i=1.d0/sr3 ) ! sqrt(3), 1/sqrt(3)
      TYPE(g_field1)  x 
      TYPE(a_field)   c
      COMPLEX*16 x11, x12, x13, x21, x22, x23, x31, x32, x33

      do i = 1, nv

         x11 = x%g(1,1,i)
         x12 = x%g(1,2,i)
         x13 = x%g(1,3,i)
         x21 = x%g(2,1,i)
         x22 = x%g(2,2,i)
         x23 = x%g(2,3,i)
         x31 = x%g(3,1,i)
         x32 = x%g(3,2,i)
         x33 = x%g(3,3,i)

         c%a(1,i) = ( AIMAG(x12) + AIMAG(x21) )
         c%a(2,i) = ( REAL (x12) - REAL (x21) )
         c%a(3,i) = ( AIMAG(x11) - AIMAG(x22) )
         c%a(4,i) = ( AIMAG(x13) + AIMAG(x31) )
         c%a(5,i) = ( REAL (x13) - REAL (x31) )
         c%a(6,i) = ( AIMAG(x23) + AIMAG(x32) )
         c%a(7,i) = ( REAL (x23) - REAL (x32) )
         c%a(8,i) = sr3i 
     &              * ( AIMAG(x11) + AIMAG(x22)
     &                 -2.d0*AIMAG(x33) )

      enddo

      return
      end

c----------------------------------------------------------------c
      subroutine  gprojct0(x,uout)
c----------------------------------------------------------------c
c     Project Lie algebra x to Lie group element g.
c     g = exp(i*x(1)*lambda(1)/2) * exp(i*x(2)*lambda(2)/2) * ....
c     (approximation)
c     Output is in uout
c----------------------------------------------------------------c
c     !!!!!   This is for SU(3)   !!!!!
c----------------------------------------------------------------c
      USE hmc_mod
      USE field_g
      include '../INCLUDE/para_geometry'
      TYPE(a_field) x
      TYPE(g_field1)  uout, tmp1, tmp2
*     real  x(8,nv),tmp1(18,nv),tmp2(18,nv),uout(18,nv)
    
      call  grplink(x,tmp1,1,NV)
      call  grplink(x,tmp2,2,NV)
*     call  prodlink(tmp1,tmp2,uout,nv,1)
      uout = tmp1*tmp2
      call  grplink(x,tmp2,3,NV)
*     call  prodlink(uout,tmp2,tmp1,nv,1)
      tmp1 = uout*tmp2
      call  grplink(x,tmp2,4,NV)
*     call  prodlink(tmp1,tmp2,uout,nv,1)
      uout = tmp1*tmp2
      call  grplink(x,tmp2,5,NV)
*     call  prodlink(uout,tmp2,tmp1,nv,1)
      tmp1 = uout*tmp2
      call  grplink(x,tmp2,6,NV)
*     call  prodlink(tmp1,tmp2,uout,nv,1)
      uout = tmp1*tmp2
      call  grplink(x,tmp2,7,NV)
*     call  prodlink(uout,tmp2,tmp1,nv,1)
      tmp1 = uout*tmp2
      call  grplink(x,tmp2,8,NV)
*     call  prodlink(tmp1,tmp2,uout,nv,1)
      uout = tmp1*tmp2

      return
      end

c-----------------------------------------------------c
      subroutine  grplink(t,x,iflag,nv)
c-----------------------------------------------------c
c     Projection from SU(3) Lie algebra on its grouop, i.e.,
c     x = exp( i * t(k)*Lambda(k) / 2)
c-----------------------------------------------------c
      USE hmc_mod
      USE field_g
      TYPE(a_field)  t
      TYPE(g_field1) x 
      REAL*8 sr3, sr3i, sr3i2
      parameter( sr3=1.7320508d0, sr3i=1.d0/sr3, sr3i2=2.d0*sr3i )
*     real  t(8,nv), x(18,nv)

      REAL*8  c, s, c1, s1, c2, s2 

      do k1 = 1, 3
      do k2 = 1, 3
        do i = 1, nv
           x%g(k1,k2,i) = 0.d0
        enddo
      enddo
      enddo

      SELECT CASE (iflag) 
      
         CASE(1)
           do i = 1, nv
              c = cos(t%a(1,i) * .5d0)
              s = sin(t%a(1,i) * .5d0)
              x%g(1,1,i) = c 
              x%g(1,2,i) = CMPLX(0.d0,s) 
              x%g(2,1,i) = CMPLX(0.d0,s) 
              x%g(2,2,i) = c 
              x%g(3,3,i) = 1.d0 
           enddo

         CASE(2)
            do i = 1, nv
              c = cos(t%a(2,i) * .5d0)
              s = sin(t%a(2,i) * .5d0)
              x%g(1,1,i) =  c 
              x%g(1,2,i) =  s 
              x%g(2,1,i) = -s 
              x%g(2,2,i) =  c 
              x%g(3,3,i) = 1.d0 
            enddo

         CASE(3)
            do i = 1, nv
              c = cos(t%a(3,i) * .5d0)
              s = sin(t%a(3,i) * .5d0)
              x%g(1,1,i) = CMPLX( c, s) 
              x%g(2,2,i) = CMPLX( c,-s) 
              x%g(3,3,i) = 1.d0 
            enddo

         CASE(4)
            do i = 1, nv
              c = cos(t%a(4,i) * .5d0)
              s = sin(t%a(4,i) * .5d0)
              x%g(1,1,i) =  c 
              x%g(1,3,i) =  CMPLX(0.d0,s) 
              x%g(2,2,i) =  1.d0 
              x%g(3,1,i) =  CMPLX(0.d0,s) 
              x%g(3,3,i) =  c 
            enddo

         CASE(5)
            do i = 1, nv
              c = cos(t%a(5,i) * .5d0)
              s = sin(t%a(5,i) * .5d0)
              x%g(1,1,i) =  c 
              x%g(1,3,i) =  s 
              x%g(2,2,i) =  1.d0 
              x%g(3,1,i) = -s 
              x%g(3,3,i) =  c 
            enddo

         CASE(6)
            do i = 1, nv
              c = cos(t%a(6,i) * .5d0)
              s = sin(t%a(6,i) * .5d0)
              x%g(1,1,i) =  1.d0 
              x%g(2,2,i) =  c 
              x%g(2,3,i) =  CMPLX(0.d0,s) 
              x%g(3,2,i) =  CMPLX(0.d0,s) 
              x%g(3,3,i) =  c 
            enddo

         CASE(7)
            do i = 1, nv
              c = cos(t%a(7,i) * .5d0)
              s = sin(t%a(7,i) * .5d0)
              x%g(1,1,i) =  1.d0 
              x%g(2,2,i) =  c 
              x%g(2,3,i) =  s 
              x%g(3,2,i) = -s 
              x%g(3,3,i) =  c 
            enddo

         CASE(8)
            do i = 1, nv
              c1 = cos( sr3i *t%a(8,i) * .5d0)
              s1 = sin( sr3i *t%a(8,i) * .5d0)
              c2 = cos( sr3i2*t%a(8,i) * .5d0)
              s2 = sin( sr3i2*t%a(8,i) * .5d0)
              x%g(1,1,i) =  CMPLX( c1, s1) 
              x%g(2,2,i) =  CMPLX( c1, s1) 
              x%g(3,3,i) =  CMPLX( c2,-s2) 
            enddo

         CASE DEFAULT
            WRITE(*,*) "The value of iflag is strange ? iflag=", iflag 
            STOP

      END SELECT

      return
      end

c----------------------------------------------------------------c
      subroutine gprojct(u,v)
c----------------------------------------------------------------c
c     Project Lie algebra u to Lie group element g.
c     U = u(1)*lambda(1)/2 + u(2)*lambda(2)/2) + ....
c     g = exp(i*U) 
c     ( by exact diagonalization                           )
c     (  X = Adj(V) * D * V                                )
c     (  where D = diag(d1,d2,d3)i                         )
c     (  g = Adj(V) * diag(exp(i*d1),exp(i*d2),exp(i*d3))) )
c----------------------------------------------------------------c
c for i = 1 to n, v(i) is the SU(3) projection of Lie algebra element u(i)
c     !!!!!   This is for SU(3)   !!!!!
c----------------------------------------------------------------c
      USE hmc_mod
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field1) v, w, ww
      TYPE(a_field)  u

      REAL*8  v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,
     &        v13,v14,v15,v16,v17,v18
      REAL*8  w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,
     &        w13,w14,w15,w16,w17,w18
      REAL*8  ww1, ww2, ww3, ww4, ww5, ww6,
     &        ww7, ww8, ww9, ww10,ww11,ww12,
     &        ww13,ww14,ww15,ww16,ww17,ww18
      REAL*8  trv3, cofac, p3, q, x, arg, theta, e1, e2, e3
      REAL*8  c1, s1, c2, s2, c3, s3
      REAL*8  coeff
      REAL*8  pi, pi23
      parameter( pi=3.141592653589, pi23=2.*pi/3.d0 )

c construct Hermitian matrix v
      call algemat(u,v)

      do i = 1, NV 

        v1  = REAL ( v%g(1,1,i) )
        v2  = AIMAG( v%g(1,1,i) )
        v3  = REAL ( v%g(1,2,i) )
        v4  = AIMAG( v%g(1,2,i) )
        v5  = REAL ( v%g(1,3,i) )
        v6  = AIMAG( v%g(1,3,i) )
        v7  = REAL ( v%g(2,1,i) )
        v8  = AIMAG( v%g(2,1,i) )
        v9  = REAL ( v%g(2,2,i) )
        v10 = AIMAG( v%g(2,2,i) )
        v11 = REAL ( v%g(2,3,i) )
        v12 = AIMAG( v%g(2,3,i) )
        v13 = REAL ( v%g(3,1,i) )
        v14 = AIMAG( v%g(3,1,i) )
        v15 = REAL ( v%g(3,2,i) )
        v16 = AIMAG( v%g(3,2,i) )
        v17 = REAL ( v%g(3,3,i) )
        v18 = AIMAG( v%g(3,3,i) )

c find eigenvalues of v
        trv3 = (v1 + v9 + v17) / 3.d0
        cofac = v1 *  v9  - v3**2  - v4**2
     &        + v1 * v17 -  v5**2  - v6**2
     &        + v9 * v17 - v11**2 - v12**2
        det = v1 * v9 * v17
     &      - v1 * (v11**2 + v12**2)
     &      - v9 * (v5**2 + v6**2)
     &      - v17 * (v3**2 + v4**2)
     &      + (v5 * (v3 * v11 - v4 * v12)
     &      +  v6 * (v3 * v12 + v4 * v11)) * 2.d0
        p3 = cofac / 3.d0 - trv3 ** 2
        q = trv3 * cofac - det - 2.d0 * trv3 ** 3
        x = sqrt(-4.d0 * p3)
        arg = q / (x * p3)
        arg = min(1.d0, max(-1.d0, arg))
        theta = acos(arg) / 3.d0
        e1 = x * cos(theta) + trv3
        theta = theta + pi23
        e2 = x * cos(theta) + trv3
c       theta = theta + pi23
c       e3 = x * cos(theta) + trv3
        e3 = 3.d0 * trv3 - e1 - e2

c solve for eigenvectors

        w1 =   v5 * (v9 - e1) - v3 * v11 + v4 * v12
        w2 = - v6 * (v9 - e1) + v4 * v11 + v3 * v12
        w3 =   (v1 - e1) * v11 - v3 * v5 - v4 * v6
        w4 = - (v1 - e1) * v12 - v4 * v5 + v3 * v6
        w5 = - (v1 - e1) * (v9 - e1) +  v3**2 + v4**2
        w6 = 0.d0

        coeff = 1.d0 / sqrt(w1**2 + w2**2 + w3**2 + w4**2 + w5**2)

        w1 = w1 * coeff
        w2 = w2 * coeff
        w3 = w3 * coeff
        w4 = w4 * coeff
        w5 = w5 * coeff
 
        w7 =   v5 * (v9 - e2) - v3 * v11 + v4 * v12
        w8 = - v6 * (v9 - e2) + v4 * v11 + v3 * v12
        w9 = (v1 - e2) * v11 - v3 * v5 - v4 * v6
        w10 = - (v1 - e2) * v12 - v4 * v5 + v3 * v6
        w11 = - (v1 - e2) * (v9 - e2) +  v3**2 + v4**2
        w12 = 0.d0
 
        coeff = 1.d0 / sqrt(w7**2  + w8**2 + w9**2
     &                    + w10**2 + w11**2)
        w7  = w7  * coeff
        w8  = w8  * coeff
        w9  = w9  * coeff
        w10 = w10 * coeff
        w11 = w11 * coeff
 
        w13 =   v5 * (v9 - e3) - v3 * v11 + v4 * v12
        w14 = - v6 * (v9 - e3) + v4 * v11 + v3 * v12
        w15 =   (v1 - e3) * v11 - v3 * v5 - v4 * v6
        w16 = - (v1 - e3) * v12 - v4 * v5 + v3 * v6
        w17 = - (v1 - e3) * (v9 - e3) +  v3**2 + v4**2
        w18 = 0.d0

        coeff = 1.d0 / sqrt(w13**2 + w14**2 + w15**2
     &                    + w16**2 + w17**2)
        w13 = w13 * coeff
        w14 = w14 * coeff
        w15 = w15 * coeff
        w16 = w16 * coeff
        w17 = w17 * coeff
 
c construct the projection v
        c1 = cos(e1)
        s1 = sin(e1)
        ww1  = w1  * c1 - w2 * s1
        ww2  = w2  * c1 + w1 * s1
        ww3  = w3  * c1 - w4 * s1
        ww4  = w4  * c1 + w3 * s1
        ww5  = w5  * c1 - w6 * s1
        ww6  = w6  * c1 + w5 * s1

        c2 = cos(e2)
        s2 = sin(e2)
        ww7  = w7  * c2 - w8 * s2
        ww8  = w8  * c2 + w7 * s2
        ww9  = w9  * c2 - w10 * s2
        ww10 = w10 * c2 + w9 * s2
        ww11 = w11 * c2 - w12 * s2
        ww12 = w12 * c2 + w11 * s2

        c3 = cos(e3)
        s3 = sin(e3)
        ww13 = w13 * c3 - w14 * s3
        ww14 = w14 * c3 + w13 * s3
        ww15 = w15 * c3 - w16 * s3
        ww16 = w16 * c3 + w15 * s3
        ww17 = w17 * c3 - w18 * s3
        ww18 = w18 * c3 + w17 * s3

        w%g(1,1,i) = CMPLX( w1,  w2  )
        w%g(1,2,i) = CMPLX( w3,  w4  )
        w%g(1,3,i) = CMPLX( w5,  w6  )
        w%g(2,1,i) = CMPLX( w7,  w8  )
        w%g(2,2,i) = CMPLX( w9,  w10 )
        w%g(2,3,i) = CMPLX( w11, w12 )
        w%g(3,1,i) = CMPLX( w13, w14 )
        w%g(3,2,i) = CMPLX( w15, w16 )
        w%g(3,3,i) = CMPLX( w17, w18 )

        ww%g(1,1,i) = CMPLX( ww1,  ww2  )
        ww%g(1,2,i) = CMPLX( ww3,  ww4  )
        ww%g(1,3,i) = CMPLX( ww5,  ww6  )
        ww%g(2,1,i) = CMPLX( ww7,  ww8  )
        ww%g(2,2,i) = CMPLX( ww9,  ww10 )
        ww%g(2,3,i) = CMPLX( ww11, ww12 )
        ww%g(3,1,i) = CMPLX( ww13, ww14 )
        ww%g(3,2,i) = CMPLX( ww15, ww16 )
        ww%g(3,3,i) = CMPLX( ww17, ww18 )

      enddo

      v = w .ADprod. ww      !  v = w_adj * ww

      return
      end

c----------------------------------------------------------------c
      subroutine  algemat(t,x)
c----------------------------------------------------------------c
c     Construct Lie algebra Matri 
c     X = t(1)*lambda(1)/2 + t(2)*lambda(2)/2) + ....
c----------------------------------------------------------------c
c     !!!!!   This is for SU(3)   !!!!!
c----------------------------------------------------------------c
      USE field_g
      USE hmc_mod
      include '../INCLUDE/para_geometry'

      REAL*8 sr3, sr3i, sr3i2
      parameter( sr3=1.7320508075689d0,sr3i=1.d0/sr3,sr3i2=2.d0*sr3i)

*     real  t(8,nv), x(18,nv)
      REAL*8 c1, c2, c3, c4, c5, c6, c7, c8
      TYPE(a_field)   t
      TYPE(g_field1)  x

c     .....   Construct  X = Summ x(k)*Lambda/2  .....
     
      do i = 1, NV 
   
         c1 = t%a(1,i) * .5d0
         c2 = t%a(2,i) * .5d0
         c3 = t%a(3,i) * .5d0
         c4 = t%a(4,i) * .5d0
         c5 = t%a(5,i) * .5d0
         c6 = t%a(6,i) * .5d0
         c7 = t%a(7,i) * .5d0
         c8 = t%a(8,i) * .5d0
     
         x%g(1,1,i) = CMPLX( c3+sr3i*c8 ,  0.d0 )
         x%g(1,2,i) = CMPLX( c1         , -c2   )
         x%g(1,3,i) = CMPLX( c4         , -c5   )
     
         x%g(2,1,i) = CMPLX( c1         ,  c2   )
         x%g(2,2,i) = CMPLX( -c3+sr3i*c8,  0.d0 )
         x%g(2,3,i) = CMPLX( c6         , -c7   )
     
         x%g(3,1,i) = CMPLX( c4         ,  c5   )
         x%g(3,2,i) = CMPLX( c6         ,  c7   )
         x%g(3,3,i) = CMPLX( -sr3i2*c8  ,  0.d0 )
     
      enddo

      return
      end
