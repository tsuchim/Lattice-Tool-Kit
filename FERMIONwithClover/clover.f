!------------------------------------------------------c
      SUBROUTINE  MKFmunu(U,hop,csw)
!------------------------------------------------------c
!     THE CONVENTION OF THE GAMMA MATRIX HERE
!     ( EUCLIDEAN CHIRAL REPRESENTATION )
!
!               (       -i )              (       -1 )
!     GAMMA1 =  (     -i   )     GAMMA2 = (     +1   )
!               (   +i     )              (   +1     )
!               ( +i       )              ( -1       )
!
!               (     -i   )              (     -1   )
!     GAMMA3 =  (       +i )     GAMMA4 = (       -1 )
!               ( +i       )              ( -1       )
!               (   -i     )              (   -1     )
!
!               ( -1       )
!     GAMMA5 =  (   -1     )
!               (     +1   )
!               (       +1 )
!
!     ( GAMMA_MU, GAMMA_NU ) = 2*DEL_MU,NU   FOR MU,NU=1,2,3,4
!
!     SIGMA_MU_NU = i/2 [ GAMMA_MU, GAMMA_NU ]
!
!                 ( -1       )             (   -i     )
!     SIGMA_12 =  (   +1     )  SIGMA_13 = ( +i       )
!                 (     -1   )             (       -i )
!                 (       +1 )             (     +i   )
!
!                 (   -1     )             (   -1     )
!     SIGMA_14 =  ( -1       )  SIGMA_23 = ( -1       )
!                 (       +1 )             (       -1 )
!                 (     +1   )             (     -1   )
!
!                 (   +i     )             ( -1       )
!     SIGMA_24 =  ( -i       )  SIGMA_34 = (   +1     )
!                 (       -i )             (     +1   )
!                 (     +i   )             (       -1 )
!
!------------------------------------------------------------------C
!
!     Total action for fermionic parts
!    
!           S_f = \phi_bar*M_f*\phi
!
!           M_f = M_w + M_c
!
!     Clover Term :
!                              1 
!           M_c = - i*c*kappa*---*sum_{mu,nu}*sigma_{mu,nu}*F_{mu,nu}
!                              2
!
!                            1 
!               F_{mu,nu} = ---( Q_{mu,nu} - Q_{mu,nu}^{\dagger} )
!                            8 
!
!-------------------------------------------------------------------C
!     Make F_mu,nu 
!-------------------------------------------------------------------C
      USE field_g
      implicit real*8(a-h,o-z)
!--   GEOMETRICAL PARAMETERS
      include '../INCLUDE/para_geometry'

      TYPE(g_field0), INTENT(IN) :: U(4)
      REAL*8,         INTENT(IN) :: hop, csw

      COMMON/ cwork/ work1, work2, work3, work4  ! Working space
      TYPE(g_field1) work1, work2, work3, work4  
      COMPLEX*16  coe

      COMPLEX*16 fmunu(3,3,nv,6) ! 1:(mu=1,nu=2), 2:(1,3), 3:(1,4),
                                 ! 4:(2,3),       5:(2,4), 6:(3,4)

      COMMON /clover/ fmunu

*     !  ...  Check
*      WRITE(*,*) "..  Clover term Csw = ", Csw

      CALL  clearC(fmunu,3*3*nv*6)

      ! ... Fill boundary of U
      DO mu = 1, 4
        CALL set_wing_g2(U(mu))
      ENDDO

      ! ... Calculation of 4 leaves under the counter clock order.

      munu = 0
      Direction1 : DO mu = 1, 3
      Direction2 : DO nu = mu+1, 4

      munu = munu + 1
      if (munu>6) stop "munu>6 ?" 

!     1) First leaf, which is located on the right up side.
!                             
!                      (1,3)---------+
!                        |           |
!     nu                 |           | 
!      |                 |    (1)    |
!      |                 |           | 
!      |                 |           |
!      +----> mu       (1,1)-------(1,2)
!
      work1 = U(mu)
      work2 = mu.gshift.U(nu)
      work3 = work1 * work2 
      work1 = U(nu)
      work2 = nu.gshift.U(mu)
      work4 = work3.prodAD.work2 
      work3 = work4.prodAD.work1 
      CALL addfmunu(fmunu(1,1,1,munu),work3,nv)

!     2) Second leaf, which is located on the left up side.
!                             
!                      (1,4)--------(1,2)
!                        |            |
!     nu                 |            | 
!      |                 |    (2)     |
!      |                 |            | 
!      |                 |            |
!      +----> mu       (1,3)--------(1,1)
! 
      work1 = U(nu)
      work2 = ivec2(-mu,nu).gshift.U(mu)
      work3 = work1.prodAD.work2 
      work2 = (-mu).gshift.U(nu)
      work4 = (-mu).gshift.U(mu)
      work1 = work3.prodAD.work2
      work2 = work1 * work4
      CALL addfmunu(fmunu(1,1,1,munu),work2,nv)

!     3) Third leaf, which is located on the left down side.
!                             
!                      (1,2)--------(1,1)
!                        |            |
!     nu                 |            | 
!      |                 |    (3)     |
!      |                 |            | 
!      |                 |            |
!      +----> mu       (1,3)--------(1,4)
!
! 
      work1 = (-nu).gshift.U(nu)
      work2 = ivec2(-mu,-nu).gshift.U(mu)
      work3 = work2 * work1 
      work4 = ivec2(-mu,-nu).gshift.U(nu)
      work1 = work4.ADprod.work3
      work2 = (-mu).gshift.U(mu)
      work3 = work2.ADprod.work1
      CALL addfmunu(fmunu(1,1,1,munu),work3,nv)

!     4) Fourth leaf, which is located on the right down side.
!                             
!                      (1,1)--------(1,4)
!                        |            |
!     nu                 |            | 
!      |                 |    (4)     |
!      |                 |            | 
!      |                 |            |
!      +----> mu       (1,2)--------(1,3)
!
! 
      work1 = (-nu).gshift.U(nu)
      work2 = (-nu).gshift.U(mu)
      work3 = work1.ADprod.work2
      work4 = ivec2(mu,-nu).gshift.U(nu)
      work1 = work3 * work4
      work2 = U(mu)
      work3 = work1.prodAD.work2
      CALL addfmunu(fmunu(1,1,1,munu),work3,nv)

      ENDDO Direction2
      ENDDO Direction1

      !  Check
      CALL trlink(fmunu,s,6*NV)
*     WRITE(*,*) "Check plaquette : ", s/float(NV*6*4*3)
       
      do munu = 1, 6
         CALL cimaglink(fmunu(1,1,1,munu),fmunu(1,1,1,munu),nv)
      enddo 

      coe = (0.0d0,1.0d0)*0.125d0*hop*Csw

      do i = 1, 3
      do j = 1, 3
      do k = 1, nv
      do l = 1, 6
         fmunu(i,j,k,l) = coe*fmunu(i,j,k,l)  
      enddo
      enddo
      enddo
      enddo

      RETURN
      END
 
C-----------------------------------------------------------------------
      SUBROUTINE trlink(u,s,n)
C-----------------------------------------------------------------------
      COMPLEX*16, INTENT(in) :: u(3,3,n)
      REAL*8, INTENT(out):: s

      s = 0.0d0
      do i = 1, n
        s = s + u(1,1,i)+u(2,2,i) + u(3,3,i)
      enddo

      return
      end
C-----------------------------------------------------------------------
      SUBROUTINE  addfmunu (fmunu,v,nv)
C-----------------------------------------------------------------------
      USE field_g
      implicit real*8(a-h,o-z)
      complex*16, INTENT(INOUT)  :: fmunu(3,3,nv)
      TYPE(g_field1), INTENT(IN) :: v

      do i = 1, nv

         fmunu(1,1,i) = fmunu(1,1,i) + v%g( 1,1,i)
         fmunu(1,2,i) = fmunu(1,2,i) + v%g( 1,2,i)
         fmunu(1,3,i) = fmunu(1,3,i) + v%g( 1,3,i)
         fmunu(2,1,i) = fmunu(2,1,i) + v%g( 2,1,i)
         fmunu(2,2,i) = fmunu(2,2,i) + v%g( 2,2,i)
         fmunu(2,3,i) = fmunu(2,3,i) + v%g( 2,3,i)
         fmunu(3,1,i) = fmunu(3,1,i) + v%g( 3,1,i)
         fmunu(3,2,i) = fmunu(3,2,i) + v%g( 3,2,i)
         fmunu(3,3,i) = fmunu(3,3,i) + v%g( 3,3,i)

      enddo

      return
      end


C-----------------------------------------------------------------------
      SUBROUTINE  vclover (vec,X)
C-----------------------------------------------------------------------
C     vec = vec + c*fmunu*(sigma_mu_nu)*X 
C-----------------------------------------------------------------------
      USE field_f
      implicit real*8(a-h,o-z)

      include '../INCLUDE/para_geometry'
      TYPE(f_field), INTENT(IN)    :: X
      TYPE(f_field), INTENT(INOUT) :: vec
      complex*16  fmunu(3,3,NV,6) 
      complex*16  y(3,4,NV) !(Color,Row for dirac indices,Site Num.)
      complex*16  cx(4)!, cy(4)
      complex*16, parameter :: ci=(0.0d0,1.0d0) 
 
      common /clover/ fmunu

      DO i = 1, NV
         y(1,1,i)=0.0d0; y(2,1,i)=0.0d0; y(3,1,i)=0.0d0
         y(1,2,i)=0.0d0; y(2,2,i)=0.0d0; y(3,2,i)=0.0d0
         y(1,3,i)=0.0d0; y(2,3,i)=0.0d0; y(3,3,i)=0.0d0
         y(1,4,i)=0.0d0; y(2,4,i)=0.0d0; y(3,4,i)=0.0d0
      ENDDO


      Color1 : do k1 = 1, 3
      Color2 : do k2 = 1, 3

      i = 0
      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
      do ix = 1, NX
      i = i + 1

            cx(1) = X%f(k2,ix,iy,iz,it,1) 
            cx(2) = X%f(k2,ix,iy,iz,it,2)
            cx(3) = X%f(k2,ix,iy,iz,it,3)
            cx(4) = X%f(k2,ix,iy,iz,it,4)

!     Change the sign:   -fmunu  --> +fmunu
            y(k1,1,i) = y(k1,1,i) 
     &                            + fmunu(k1,k2,i,1)*(-   cx(1))
     &                            + fmunu(k1,k2,i,2)*(-ci*cx(2))
     &                            + fmunu(k1,k2,i,3)*(-   cx(2))
     &                            + fmunu(k1,k2,i,4)*(-   cx(2))
     &                            + fmunu(k1,k2,i,5)*( ci*cx(2))
     &                            + fmunu(k1,k2,i,6)*(-   cx(1))

            y(k1,2,i) = y(k1,2,i) 
     &                            + fmunu(k1,k2,i,1)*(    cx(2))
     &                            + fmunu(k1,k2,i,2)*( ci*cx(1))
     &                            + fmunu(k1,k2,i,3)*(-   cx(1))
     &                            + fmunu(k1,k2,i,4)*(-   cx(1))
     &                            + fmunu(k1,k2,i,5)*(-ci*cx(1))
     &                            + fmunu(k1,k2,i,6)*(    cx(2))

            y(k1,3,i) = y(k1,3,i) 
     &                            + fmunu(k1,k2,i,1)*(-   cx(3))
     &                            + fmunu(k1,k2,i,2)*(-ci*cx(4))
     &                            + fmunu(k1,k2,i,3)*(    cx(4))
     &                            + fmunu(k1,k2,i,4)*(-   cx(4))
     &                            + fmunu(k1,k2,i,5)*(-ci*cx(4))
     &                            + fmunu(k1,k2,i,6)*(    cx(3))

            y(k1,4,i) = y(k1,4,i) 
     &                            + fmunu(k1,k2,i,1)*(    cx(4))
     &                            + fmunu(k1,k2,i,2)*( ci*cx(3))
     &                            + fmunu(k1,k2,i,3)*(    cx(3))
     &                            + fmunu(k1,k2,i,4)*(-   cx(3))
     &                            + fmunu(k1,k2,i,5)*( ci*cx(3))
     &                            + fmunu(k1,k2,i,6)*(-   cx(4))

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      enddo Color2
      enddo Color1

      do mu = 1, 4
      i = 0
      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
      do ix = 1, NX
      i = i + 1
 
         vec%f(1,ix,iy,iz,it,mu) =  
     &             + vec%f(1,ix,iy,iz,it,mu) + y(1,mu,i)
         vec%f(2,ix,iy,iz,it,mu) =  
     &             + vec%f(2,ix,iy,iz,it,mu) + y(2,mu,i)
         vec%f(3,ix,iy,iz,it,mu) =  
     &             + vec%f(3,ix,iy,iz,it,mu) + y(3,mu,i)

      enddo
      enddo
      enddo
      enddo

      enddo

      return
      end

c------------------------------------------------------c
      subroutine cimaglink(x,y,nv)
c------------------------------------------------------c
c     y = x - x_aj
c------------------------------------------------------c
      implicit real*8(a-h,o-z)
      complex*16 x(3,3,nv),y(3,3,nv)
      complex*16 z11, z12, z13, z22, z23, z33

      do i = 1, nv

       z11 = x(1,1,i) - conjg(x(1,1,i)) 
       z12 = x(1,2,i) - conjg(x(2,1,i)) 
       z13 = x(1,3,i) - conjg(x(3,1,i)) 

       z22 = x(2,2,i) - conjg(x(2,2,i)) 
       z23 = x(2,3,i) - conjg(x(3,2,i)) 

       z33 = x(3,3,i) - conjg(x(3,3,i)) 
 
        y(1,1,i) = z11
        y(1,2,i) = z12
        y(1,3,i) = z13

        y(2,1,i) = -conjg(z12)
        y(2,2,i) = z22
        y(2,3,i) = z23 

        y(3,1,i) = -conjg(z13) 
        y(3,2,i) = -conjg(z23) 
        y(3,3,i) = z33

      enddo

      return
      end

c------------------------------------------------------c
      subroutine cimaglink2(x,nv)
c------------------------------------------------------c
c     x <= x - x_aj
c     (Same as cimaglink, but for TYPE(g_field1)
c     The output overwrites the input.
c------------------------------------------------------c
      USE field_g

      INTEGER,        INTENT(IN)  :: nv
      TYPE(g_field1), INTENT(INOUT)  :: x 

      COMPLEX*16 z11, z12, z13, z22, z23, z33

      DO i = 1, nv

        z11 = x%g(1,1,i) - CONJG(x%g(1,1,i)) 
        z12 = x%g(1,2,i) - CONJG(x%g(2,1,i)) 
        z13 = x%g(1,3,i) - CONJG(x%g(3,1,i)) 

        z22 = x%g(2,2,i) - CONJG(x%g(2,2,i)) 
        z23 = x%g(2,3,i) - CONJG(x%g(3,2,i)) 

        z33 = x%g(3,3,i) - CONJG(x%g(3,3,i)) 
 
        x%g(1,1,i) = z11
        x%g(1,2,i) = z12
        x%g(1,3,i) = z13
        x%g(2,1,i) = -CONJG(z12)
        x%g(2,2,i) = z22
        x%g(2,3,i) = z23 
        x%g(3,1,i) = -CONJG(z13) 
        x%g(3,2,i) = -CONJG(z23) 
        x%g(3,3,i) = z33

      ENDDO

      RETURN
      END

!---------------------------------------------------------------!
      SUBROUTINE  clearC(x,n)
!---------------------------------------------------------------!
      COMPLEX*16 x(n)

      DO i = 1, n

        x(i) = (0.d0,0.d0)

      ENDDO

      RETURN
      END

!------------------------------------------------------c
      SUBROUTINE  MKFmunuTest(U,hop,csw)
!------------------------------------------------------c
!     Make Dummy F_mu,nu for Test 
!-------------------------------------------------------------------C
      USE field_g
      implicit real*8(a-h,o-z)
!--   GEOMETRICAL PARAMETERS
      include '../INCLUDE/para_geometry'

      TYPE(g_field0), INTENT(IN) :: U(4)
      REAL*8,         INTENT(IN) :: hop, csw

      COMPLEX*16 fmunu(3,3,nv,6) ! 1:(mu=1,nu=2), 2:(1,3), 3:(1,4),
                                 ! 4:(2,3),       5:(2,4), 6:(3,4)

      COMMON /clover/ fmunu

      CALL  clearC(fmunu,3*3*nv*6)

      write(*,*) "Dummy routine MKFmunuTest is called."

      munu = 6

      DO ic2 = 1, 3
      DO ic1 = 1, 3

        IF(ic1==ic2) THEN
           fmunu(ic1,ic2,1,munu) = 1.d0
        ELSE
           fmunu(ic1,ic2,1,munu) = 0.d0
        ENDIF

      ENDDO
      ENDDO

      RETURN
      END

!---------------------------------------------------------c
      SUBROUTINE CswEstimate1(beta,Csw)
!---------------------------------------------------------c
!     M.Luscher, S.Sint, R.Sommer and P.Weisz, 
!     Nucl. Phys. B491 (1997) 323 
!     dynamical quark 
!---------------------------------------------------------c
      REAL*8, INTENT(IN)  ::  beta
      REAL*8, INTENT(OUT) ::  Csw

      REAL*8  g02
      REAL*8  n, n0, n2, n4, n6
      REAL*8  d, d0, d2

*     beta = 5.4d0

      g02 = 6.d0/beta
      g04 = g02**2
      g06 = g02**3

      n0 = 1.d0
      n2 = -0.656d0
      n4 = -0.152d0
      n6 = -0.054d0
    
      d0 = 1.d0
      d2 = -0.922d0

      n = n0 + n2*g02 + n4*g04 + n6*g06
      d = d0 + d2*g02

      Csw = n/d
*      Csw = 0.d0            ! Wilson fermion case

*      write(*,'(1x,a,e15.7)') "Csw = ", Csw
*      write(*,'(1x,a,e15.7)') "beta = ", beta

      RETURN
      END
!---------------------------------------------------------c
      SUBROUTINE CswEstimate2(beta,Csw)
!---------------------------------------------------------c
!     J.Jansen and R.Sommer, Nucl.Phys. B530 (1998) 185
!     dynamical quark 
!---------------------------------------------------------c
      REAL*8, INTENT(IN)  ::  beta
      REAL*8, INTENT(OUT) ::  Csw

      REAL*8  g02
      REAL*8  n, n0, n2, n4, n6, n8
      REAL*8  d, d0, d2

*     beta = 5.4d0

      g02 = 6.d0/beta
      g04 = g02**2
      g06 = g02**3
      g08 = g02**4

      n0 = 1.d0
      n2 = -0.454d0
      n4 = -0.175d0
      n6 = 0.012d0
      n8 = 0.045d0
    
      d0 = 1.d0
      d2 = -0.720d0

      n = n0 + n2*g02 + n4*g04 + n6*g06 + n8*g08
      d = d0 + d2*g02

      Csw = n/d
*      Csw = 0.d0            ! Wilson fermion case

*      write(*,'(1x,a,e15.7)') "Csw = ", Csw
*      write(*,'(1x,a,e15.7)') "beta = ", beta

      RETURN
      END
       
