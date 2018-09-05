!-----------------------------------------------------!
      SUBROUTINE  dSclover(mu,z,u)
!-----------------------------------------------------!
!     Calculate   dS_clover/dA_mu(x)
!-----------------------------------------------------c
!     Calculate
!        \partial S_c/\partial A^(a)_\mu(x) 
!    
!           S_f = \phi_bar*M_f*\phi = S_w + S_c
!           M_f = M_w + M_c
!-----------------------------------------------------C
!       mu     nu      munu   |   mu     nu      munu
!        1      2         1   |    2      1        -1
!        1      3         2   |    3      1        -2
!        1      4         3   |    4      1        -3
!        2      3         4   |    3      2        -4
!        2      4         5   |    4      2        -5
!        3      4         6   |    4      3        -6
!-----------------------------------------------------!
      USE field_g
      USE fpara
      USE field_f
      USE hmc_mod

      implicit none

      include '../INCLUDE/para_geometry'

      INTEGER,       INTENT(IN)  :: mu
      TYPE(a_field), INTENT(OUT) :: z
      TYPE(g_field0),INTENT(IN)  :: u(4)

      TYPE(g_field1) dF1(8), dF2(8)

      INTEGER  inn(NV,4,2)
      SAVE     inn

      ! ...  Working place
      TYPE(g_field1) work1, work2, work3, work4  
      TYPE(g_field1) gtmp1, gtmp2, gtmp3, gtmp4
      INTEGER  is1, is2
      COMMON/  / work1, work2, work3, work4,
     &           gtmp1, gtmp2, gtmp3, gtmp4, 
     &           is1(NV), is2(NV)

      COMPLEX*16, DIMENSION(3,NV,4) :: v1, v2
      COMPLEX*16, DIMENSION(3,NV,4) :: veta, vxi
      COMPLEX*16, DIMENSION(3,NV,4) :: ftmp1, ftmp2
      real*8 plaq

      INTEGER ix, iy, iz, it, nu, is, ia, ic, ialpha 
      INTEGER iflag, memo
      DATA  memo/ 1/
      SAVE  memo

      ! ...  Make table (only once)  ....................
      IF( memo==1 )  THEN
         CALL fprep3(NX,NY,NZ,NT,inn)
         memo = 0
      ENDIF
      ! .................................................

      ! ... Clear z
      DO ia = 1, 8
      DO is = 1, NV
        z%a(ia,is) = 0.d0
      ENDDO
      ENDDO

      ! ... Set v1=eta, v2=xi 
      DO ialpha = 1, 4
      DO ic = 1, 3
      is = 0
      DO it = 1, NT
      DO iz = 1, NZ
      DO iy = 1, NY
      DO ix = 1, NX
      is = is + 1
    
        veta(ic,is,ialpha) = eta%f(ic,ix,iy,iz,it,ialpha)
        vxi (ic,is,ialpha) =  xi%f(ic,ix,iy,iz,it,ialpha)

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ENDDO
      ENDDO

      DirectionNu : DO nu = 1, 4
      IF(nu==mu) CYCLE

      ! .....................  !
      !    Case 1 and 3        !
      ! .....................  !
        iflag = 1

        CALL  Cal_dFmunu(U,dF1,dF2,
     &        work1, work2, work3, work4,
     &        gtmp1, gtmp2, gtmp3, gtmp4,
     &        mu,nu,iflag)

        ! ... Case 1
        CALL  VxSigxV (veta,vxi,dF1,z,ftmp1,ftmp2,mu,nu)

        ! ... Case 3 
        DO is = 1, NV
          is1(is) = inn(is,mu,1)
        ENDDO

        DO is = 1, NV
          is2(is) = inn(is1(is),nu,1)
          ! For Check 
*         IF( is<5 ) THEN
*            WRITE(*,*) "mu, nu : ", mu, nu
*            WRITE(*,*) "is, is1, is2 : ", is, is1(is),is2(is)
*         ENDIF
        ENDDO
 
        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is2(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is2(is),ialpha)
        ENDDO
        ENDDO
        ENDDO
         
        CALL  VxSigxV (v1,v2,dF2,z,ftmp1,ftmp2,mu,nu)


      ! .....................  !
      !    Case 2 and 4        !
      ! .....................  !

        iflag=2 

        CALL  Cal_dFmunu(U,dF1,dF2,
     &        work1, work2, work3, work4,
     &        gtmp1, gtmp2, gtmp3, gtmp4,
     &        mu,nu,iflag)

        ! ... Case 2
        DO is = 1, NV
          is1(is) = inn(is,mu,1)
        ENDDO

        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is1(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is1(is),ialpha)
        ENDDO
        ENDDO
        ENDDO

        CALL  VxSigxV (v1,v2,dF1,z,ftmp1,ftmp2,mu,nu)

        ! ... Case 4
        DO is = 1, NV
          is1(is) = inn(is,nu,1)
        ENDDO

        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is1(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is1(is),ialpha)
        ENDDO
        ENDDO
        ENDDO
         
        CALL  VxSigxV (v1,v2,dF2,z,ftmp1,ftmp2,mu,nu)

      ! .....................  !
      !    Case 4' and 2'      !
      ! .....................  !
       iflag = 3  

        CALL  Cal_dFmunu(U,dF1,dF2,
     &        work1, work2, work3, work4,
     &        gtmp1, gtmp2, gtmp3, gtmp4,
     &        mu,nu,iflag)

        ! ... Case 4'

        CALL  VxSigxV (veta,vxi,dF1,z,ftmp1,ftmp2,mu,nu)

        ! ... Case 2' 
        DO is = 1, NV
          is1(is) = inn(is,mu,1)
        ENDDO
        DO is = 1, NV
          is2(is) = inn(is1(is),nu,2)
        ENDDO

        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is2(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is2(is),ialpha)
        ENDDO
        ENDDO
        ENDDO
         
        CALL  VxSigxV (v1,v2,dF2,z,ftmp1,ftmp2,mu,nu)

      ! .....................  !
      !    Case 3' and 1'      !
      ! .....................  !
       iflag = 4

        CALL  Cal_dFmunu(U,dF1,dF2,
     &        work1, work2, work3, work4,
     &        gtmp1, gtmp2, gtmp3, gtmp4,
     &        mu,nu,iflag)

        ! ... Case 3'
        DO is = 1, NV
          is1(is) = inn(is,mu,1)
        ENDDO

        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is1(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is1(is),ialpha)
        ENDDO
        ENDDO
        ENDDO

        CALL  VxSigxV (v1,v2,dF1,z,ftmp1,ftmp2,mu,nu)

        ! ... Case 1' 
        DO is = 1, NV
          is1(is) = inn(is,nu,2)
        ENDDO

        DO ic = 1, 3
        DO ialpha = 1, 4
        DO is = 1, NV
          v1(ic,is,ialpha) = veta(ic,is1(is),ialpha)
          v2(ic,is,ialpha) = vxi (ic,is1(is),ialpha)
        ENDDO
        ENDDO
        ENDDO
         
        CALL  VxSigxV (v1,v2,dF2,z,ftmp1,ftmp2,mu,nu)

      ENDDO DirectionNu


      ! z = -2 Re ( .... )
      DO ia = 1, 8
      DO is = 1, NV
        z%a(ia,is) = -2.d0 * z%a(ia,is) 
      ENDDO
      ENDDO

      RETURN
      END

!-----------------------------------------------------!
      SUBROUTINE  Cal_dFmunu(U,dFmunu1,dFmunu2,
     &            work1, work2, work3, work4,
     &            temp1, temp2, temp3, temp4,
     &            mu,nu,iflag)
!-----------------------------------------------------!
!     Calculate dFmunu/dA_\mu(x)
!     iflag=1  --> Case 1 and 3
!     iflag=2  --> Case 2 and 4
!     iflag=3  --> Case 4' and 2'
!     iflag=4  --> Case 3' and 1'
!-----------------------------------------------------!
      USE field_g
      USE fpara

      implicit none

!--   GEOMETRICAL PARAMETERS
      include '../INCLUDE/para_geometry'

      TYPE(g_field0), INTENT(IN)  :: U(4)
      INTEGER,        INTENT(IN)  :: mu, nu, iflag
      TYPE(g_field1), INTENT(OUT) :: dFmunu1(8), dFmunu2(8) 

      ! ...  Working place
      TYPE(g_field1) work1, work2, work3, work4  
      TYPE(g_field1) temp1, temp2, temp3, temp4

      TYPE(ivec2)    idir2 
      COMPLEX*16, parameter :: ci=(0.0d0,1.0d0)
      COMPLEX*16  coe

      real*8 plaq
      INTEGER ia

      coe = (0.0d0,1.0d0)*0.125d0*hop*Csw

!    ...  A) the upper staple  ....................
!                               w3
!       nu                 .-----------+
!                          |           | 
!       /|             w4  |           | w2
!        |                 |           |
!        |                 |           |
!        +----> mu         .-----------.
!                          x    w1

      IF( (iflag==1).or.(iflag==2) )  THEN

        work1 = U(mu)
        work2 = mu.gshift.U(nu)
        work3 = nu.gshift.U(mu)
        work4 = U(nu)

      ENDIF

      ! ...  iflag = 1 (Case 1 and 3)
      IF(iflag==1)  THEN
!                               w3
!                          .-----------+
!  Case 1                  |           |
!       nu                 |           |
!        |              w4 |           | w2
!        |                 |           |
!        +----> mu         o-----------. 
!                          x    w1 
!
!                               w3
!                          .<----------o
!  Case 3                  |           |
!       nu                 |           |
!        |              w4 |           | w2
!        |                 |           |
!        +----> mu         .-----------. 
!                          x    w1 

        temp1 = work1 * work2
        temp2 = work4 * work3

        DO ia = 1, 8

          CALL LambdaMul(temp1,temp3,ia)  ! temp3=(lambda_a/2)*temp1

          temp4 = temp3 .prodAD.temp2
          temp4 = ci * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu1(ia) = coe * temp4

          temp4 = temp2 .ADprod.temp3
          temp4 = ci * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu2(ia) = coe * temp4

        ENDDO

      ENDIF

      ! ...  iflag = 2 (Case 2 and 4)
      IF(iflag==2)  THEN
!                                w3
!                          .-----------+
!  Case 2                  |           |
!       nu                 |           |
!        |             w4  |           | w2
!        |                 |           |
!        +----> mu         .---------->o 
!                          x    w1 
!
!                                w3
!                          o-----------+
!  Case 4                  |           |
!       nu                 |           |
!        |              w4 |           | w2
!        |                 |           |
!        +----> mu         .-----------. 
!                          x    w1 


        temp1 = work2 .prodAD. work3

        DO ia = 1, 8

          CALL LambdaMul(work1,temp2,ia)  ! temp2=(lambda_a/2)*work1

          temp3 = work4 .ADprod. temp2 
          temp4 = temp1 * temp3
          temp4 = ci * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu1(ia) = coe * temp4

          temp4 = temp3 * temp1
          temp4 = ci * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu2(ia) = coe * temp4

        ENDDO

      ENDIF

!    ...  B) the lower staple  ....................
!
!       nu
!       /|
!        |
!        |
!        |                 x    w1
!        +----> mu         .-----------+
!                          |           |
!                      w4  |           | w2
!                          |           |
!                          |           |
!                          .-----------.
!                          x    w3


      IF( (iflag==3).or.(iflag==4) )  THEN

        work1 = U(mu)
        idir2 = ivec2(mu,-nu)
        work2 = idir2.gshift.U(nu)
        work3 = (-nu).gshift.U(mu)
        work4 = (-nu).gshift.U(nu)

      ENDIF
 
      !   ...  iflag = 3 (Case 4' and 2')
      IF(iflag==3)  THEN

!  Case 4'
!       nu
!        |
!        |                 x    w1
!        +----> mu         o<----------. 
!                          |           |
!                      w4  |           | w2
!                          |           |
!                          |           |
!                          .-----------.
!                               w3
!  Case 2'
!       nu
!        |
!        |                 x    w1
!        +----> mu         .<----------. 
!                          |           |
!                      w4  |           | w2
!                          |           |
!                          |           |
!                          .-----------o
!                               w3

        temp1 = work4.ADprod.work3

        DO ia = 1, 8

          CALL LambdaMul(work1,temp2,ia)  ! temp2=(lambda_a/2)*work1
          temp3 = work2 .prodAD. temp2 

          temp4 = temp1 * temp3
          temp4 = (-ci) * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu1(ia) = coe * temp4

          temp4 = temp3 * temp1
          temp4 = (-ci) * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu2(ia) = coe * temp4

        ENDDO

      ENDIF

      ! ...  iflag = 4 (Case 3' and 1')

      IF(iflag==4)  THEN

!  Case 3' 
!       nu
!        |
!        |                 x    w1
!        +----> mu         .<----------o 
!                          |           |
!                      w4  |           | w2
!                          |           |
!                          |           |
!                          .-----------.
!                               w3
!  Case 1' 
!       nu
!        |
!        |                 x    w1
!        +----> mu         .<----------. 
!                          |           |
!                      w4  |           | w2
!                          |           |
!                          |           |
!                          o-----------.
!                               w3

        temp1 = work3 * work2

        DO ia = 1, 8

          CALL LambdaMul(work1,temp2,ia)  ! temp2=(lambda_a/2)*work1
          temp3 = work4 * temp2 

          temp4 = temp3 .ADprod. temp1 
          temp4 = (-ci) * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu1(ia) = coe * temp4

          temp4 = temp1 .prodAD. temp3 
          temp4 = (-ci) * temp4
          CALL cimaglink2(temp4,nv)
          dFmunu2(ia) = coe * temp4

        ENDDO

      ENDIF

      RETURN
      END

!-----------------------------------------------------------------------
      SUBROUTINE  VxSigxV  (v1,v2,u,z,tmp1,tmp2,mu,nu)
!-----------------------------------------------------------------------
!     z(x) = v1^{\dagger}(x) * u(x) * \sigma_{\mu,\nu} * v2(x) 
!     v1, v2 : vectors
!     u : 3x3 matrix
!     sigma : 4x4 matrix 
!     mu, nu : input fixed
!     tmp1, tmp2 : work spaces
!-----------------------------------------------------------------------
      USE field_f
      USE hmc_mod 

      implicit none

      include '../INCLUDE/para_geometry'

      COMPLEX*16, DIMENSION(3,NV,4), INTENT(IN)  :: v1, v2
      TYPE(g_field1),                INTENT(IN)  :: u(8)
      INTEGER,                       INTENT(IN)  :: mu,nu
      TYPE(a_field),                 INTENT(INOUT) :: z
      COMPLEX*16, DIMENSION(3,NV,4), INTENT(INOUT) :: tmp1, tmp2 

      COMPLEX*16, parameter :: ci=(0.0d0,1.0d0)
      COMPLEX*16  s1, s2, s3
      REAL*8 facmunu

      INTEGER  ic, ia, is, nu0, mu0, ialpha

*     WRITE(*,*) "NC: ", Nc

      IF(mu==nu) THEN
         write(*,*)"mu should not be equal to nu (VsSigV)"
         write(*,*) "mu, nu : ", mu, nu
         stop
      ENDIF

      IF(mu<nu) THEN
         facmunu = 1.d0
         mu0 = mu
         nu0 = nu
      ELSE
         facmunu = -1.d0
         mu0 = nu
         nu0 = mu
      ENDIF


      ! ... \sigma_{\mu,\nu} x v2 (in Dirac space)  .... 
      Index_mu : SELECT CASE (mu0)

        CASE (1)  
           SELECT CASE (nu0)
             CASE (2)           ! mu0=1, nu0=2

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = -v2(ic,is,1)
                tmp1(ic,is,2) = +v2(ic,is,2)
                tmp1(ic,is,3) = -v2(ic,is,3)
                tmp1(ic,is,4) = +v2(ic,is,4)
               ENDDO
               ENDDO

             CASE (3)           ! mu0=1, nu0=3

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = -ci * v2(ic,is,2)
                tmp1(ic,is,2) = +ci * v2(ic,is,1)
                tmp1(ic,is,3) = -ci * v2(ic,is,4)
                tmp1(ic,is,4) = +ci * v2(ic,is,3)
               ENDDO
               ENDDO

             CASE (4)           ! mu0=1, nu0=4

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = -v2(ic,is,2)
                tmp1(ic,is,2) = -v2(ic,is,1)
                tmp1(ic,is,3) = +v2(ic,is,4)
                tmp1(ic,is,4) = +v2(ic,is,3)
               ENDDO
               ENDDO

             CASE DEFAULT
                write(*,*) "something is wrong in VsSigV"
                write(*,*) "mu,nu : ", mu, nu
                stop 

           END SELECT

        CASE (2)
           SELECT CASE (nu0)
             CASE (3)           ! mu0=2, nu0=3

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = -v2(ic,is,2)
                tmp1(ic,is,2) = -v2(ic,is,1)
                tmp1(ic,is,3) = -v2(ic,is,4)
                tmp1(ic,is,4) = -v2(ic,is,3)
               ENDDO
               ENDDO

             CASE (4)           ! mu0=2, nu0=4

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = +ci * v2(ic,is,2)
                tmp1(ic,is,2) = -ci * v2(ic,is,1)
                tmp1(ic,is,3) = -ci * v2(ic,is,4)
                tmp1(ic,is,4) = +ci * v2(ic,is,3)
               ENDDO
               ENDDO

             CASE DEFAULT
                write(*,*) "something is wrong in VsSigV"
                write(*,*) "mu,nu : ", mu, nu
                stop 

           END SELECT

        CASE (3)
           SELECT CASE (nu0)
             CASE (4)           ! mu0=3, nu0=4

               DO ic = 1, Nc 
               DO is = 1, NV
                tmp1(ic,is,1) = -v2(ic,is,1)
                tmp1(ic,is,2) = +v2(ic,is,2)
                tmp1(ic,is,3) = +v2(ic,is,3)
                tmp1(ic,is,4) = -v2(ic,is,4)
               ENDDO
               ENDDO

             CASE DEFAULT
                write(*,*) "something is wrong in VsSigV"
                write(*,*) "mu,nu : ", mu, nu
                stop 
           END SELECT

        CASE DEFAULT
           write(*,*) "something is wrong in VsSigV"
           write(*,*) "mu,nu : ", mu, nu
           stop 

      END SELECT Index_mu

      DO ic = 1, Nc
      DO is = 1, NV
        tmp2(ic,is,1) = facmunu * tmp1(ic,is,1)
        tmp2(ic,is,2) = facmunu * tmp1(ic,is,2)
        tmp2(ic,is,3) = facmunu * tmp1(ic,is,3)
        tmp2(ic,is,4) = facmunu * tmp1(ic,is,4)
      ENDDO
      ENDDO

      ! ... tmp1 = u x tmp2 (in color space) 
      Algeb : DO ia = 1, 8

      DO ialpha = 1, 4
      DO is = 1, NV

        s1 = + u(ia)%g(1,1,is)*tmp2(1,is,ialpha) 
     &       + u(ia)%g(1,2,is)*tmp2(2,is,ialpha) 
     &       + u(ia)%g(1,3,is)*tmp2(3,is,ialpha) 
        s2 = + u(ia)%g(2,1,is)*tmp2(1,is,ialpha) 
     &       + u(ia)%g(2,2,is)*tmp2(2,is,ialpha) 
     &       + u(ia)%g(2,3,is)*tmp2(3,is,ialpha) 
        s3 = + u(ia)%g(3,1,is)*tmp2(1,is,ialpha) 
     &       + u(ia)%g(3,2,is)*tmp2(2,is,ialpha) 
     &       + u(ia)%g(3,3,is)*tmp2(3,is,ialpha) 

        tmp1(1,is,ialpha) = s1
        tmp1(2,is,ialpha) = s2
        tmp1(3,is,ialpha) = s3

      ENDDO
      ENDDO

      ! ... v1^{\dagger} * tmp1
      DO is = 1, NV

        s1 = CONJG(v1(1,is,1)) * tmp1(1,is,1)
     &     + CONJG(v1(1,is,2)) * tmp1(1,is,2)
     &     + CONJG(v1(1,is,3)) * tmp1(1,is,3)
     &     + CONJG(v1(1,is,4)) * tmp1(1,is,4)

        s2 = CONJG(v1(2,is,1)) * tmp1(2,is,1)
     &     + CONJG(v1(2,is,2)) * tmp1(2,is,2)
     &     + CONJG(v1(2,is,3)) * tmp1(2,is,3)
     &     + CONJG(v1(2,is,4)) * tmp1(2,is,4)

        s3 = CONJG(v1(3,is,1)) * tmp1(3,is,1)
     &     + CONJG(v1(3,is,2)) * tmp1(3,is,2)
     &     + CONJG(v1(3,is,3)) * tmp1(3,is,3)
     &     + CONJG(v1(3,is,4)) * tmp1(3,is,4)

        z%a(ia,is) = z%a(ia,is) + (s1 + s2 + s3)

      ENDDO

      ENDDO Algeb

      RETURN
      END

C----------------------------------------------------------------------C
      SUBROUTINE  fprep3 (NX,NY,NZ,NT,inn)
C----------------------------------------------------------------------C
C   -  SITES IN THE LATTICE ARE IDENTIFIED BY ONE "ABSOLUTE" INDEX     C
C      (= CUMULATIVE INDEX 'IC'=1,...,NV)                              C
C   -  PREP PREPARES NEAREST NEIGHBOUR LIST 'INN' FOR A HYPERCUBIC     C
C      LATTICE WITH PBC                                                C
C----------------------------------------------------------------------C
      INTEGER,INTENT(IN) :: NX, NY, NZ, NT
      INTEGER,INTENT(OUT):: INN(NX*NY*NZ*NT,4,2)

C//   JF(I,L) = MOD(I+1  ,L)                          ! INDEX FORWARD
C//   JB(I,L) = MOD(I-1+L,L)                          ! INDEX BACKWARD
C                                                     ! CUMULATIVE INDEX
      JF(I,L) = MOD(I+1  ,L) 
      JB(I,L) = MOD(I-1+L,L) 
      IC(IX,IY,IZ,IT,N1,N2,N3) = 1+IX+N1*(IY+N2*(IZ+N3*IT))

      DO  IT = 0, NT-1
      DO  IZ = 0, NZ-1
      DO  IY = 0, NY-1
      DO  IX = 0, NX-1

      I=IC(IX,IY,IZ,IT,NX,NY,NZ)

      inn(I,1,1) = IC (JF(IX,NX),IY       ,IZ       ,IT       ,NX,NY,NZ)
      inn(I,2,1) = IC (IX       ,JF(IY,NY),IZ       ,IT       ,NX,NY,NZ)
      inn(I,3,1) = IC (IX       ,IY       ,JF(IZ,NZ),IT       ,NX,NY,NZ)
      inn(I,4,1) = IC (IX       ,IY       ,IZ       ,JF(IT,NT),NX,NY,NZ)
      inn(I,1,2) = IC (JB(IX,NX),IY       ,IZ       ,IT       ,NX,NY,NZ)
      inn(I,2,2) = IC (IX       ,JB(IY,NY),IZ       ,IT       ,NX,NY,NZ)
      inn(I,3,2) = IC (IX       ,IY       ,JB(IZ,NZ),IT       ,NX,NY,NZ)
      inn(I,4,2) = IC (IX       ,IY       ,IZ       ,JB(IT,NT),NX,NY,NZ)

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
