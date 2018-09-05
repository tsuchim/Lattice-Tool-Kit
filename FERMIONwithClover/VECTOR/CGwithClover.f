c----------------------------------------------------------------------
*     SUBROUTINE  cg0 (x0,b0,eps,imax,iflag)
      SUBROUTINE  cg0 (x0,b0,iflag)
c-----------------------------------------------------------------------
c     this routine  solves
c                  W(U) * X = b       (iflag=1)
c                  W(U)_adj * X = b   (iflag=2)
c     by conjugate gradient
c-----------------------------------------------------------------------
      USE field_f 
      USE fpara
      TYPE(g_field0) U(4)
      common/ config/ U

      TYPE(f_field), INTENT(IN)   :: b0
      TYPE(f_field), INTENT(INOUT):: x0

      include '../INCLUDE/para_geometry' 
      PARAMETER( NDF=2*NC*NC, NDFV=2*NC  ) 

      REAL*8 alpha, beta, c1, c2, c3, rnorm
      
      parameter (NDIM1=NV*NDFV*4, NDIM2=NV*9*NDFV*4 )
 
      REAL*8  x(NV,NDFV,4), b(NV,NDFV,4)
      REAL*8  un(NV,9,NDF)

      INTEGER inn(NV,4,2), msign(NV,4) 

      REAL*8 vtemp(NV*9*NDFV*4), 
     &       res(NV,NDFV,4), p(NV,NDFV,4), q(NV,NDFV,4), s(NV,NDFV,4), 
     &       vtmp1(NV,NDFV,4), vtmp2(NV,NDFV,4)
      common           b, x, vtemp, res, p, q, s, vtmp1, vtmp2

      data  itest/ 1/
      data  memo/ 1/
      SAVE inn, msign, memo

      IF( abs(Csw) > 0.001 ) THEN
         CALL MKFmunu(U,hop,Csw)
      ENDIF

c     ...  Prepare inn and msign.  This should be done only once. 
      IF(memo==1)  THEN
         CALL  fprep2(inn,msign)
         memo = 0
      ENDIF

c     ...  The following routine may be called only once if the
c          configuration is not changed.
      CALL  unear2(un,inn,msign) 
    
      CALL convert_f(b0,b,1)
      CALL convert_f(x0,x,1)

      IF(iflag==1) CALL vxg5(b,b,NV)

      call  clear(vtemp,NDIM2)
      call  clear(res, NDIM1)
      call  clear(p,   NDIM1)
      call  clear(q,   NDIM1)
      call  clear(s,   NDIM1)
      call  clear(vtmp1,  NDIM1)
      call  clear(vtmp2,  NDIM1)

      call  vnear  (NV,b,vtemp,inn)
      call  wxvect2 (un,b,vtemp,vtmp1,vtemp,NV)

      ! ...  res = b - gamma5*W * x    
      call  vnear  (NV,x,vtemp,inn)
      call  wxvect2 (un,x,vtemp,vtmp1,vtemp,NV)
      call  vdiff(res,b,vtmp1,NV)

      ! ...  p = (gamma5*W)_adj) * res =gamma5*W*res
      call  vnear  (NV,res,vtemp,inn)
      call  wxvect2 (un,res,vtemp,p,vtemp,NV)

      ! ...  c1 = < p | p >  
      call  vecnor (p,c1,NV) 

c     ***  Iteration starts  ***
      DO i = 1, imax 

        ! ...  q = gamma5*W * p  
        call  vnear  (NV,p,vtemp,inn)
        call  wxvect2 (un,p,vtemp,q,vtemp,NV)

        ! ...  c2 = < q | q >
        call  vecnor (q,c2,NV)

        alpha = c1 / c2

        ! ...  x   = x   + alpha * p
        call  vsum1 (x,p,+alpha,NV)

        ! ...  res = res - alpha * q     
        call  vsum1 (res,q,-alpha,NV)

c       ***  Check the coNVergence  ***
        call  vecnor(res,rnorm,NV) 
        IF(itest==1) THEN
          write(*,*) "i=", i, "  rnorm=", rnorm
        ENDIF
        if (rnorm.le.eps)  go to 5000

        ! ...  s = (gamma5*W)_adj * res = gamma5*W * res     
        call  vnear  (NV,res,vtemp,inn)
        call  wxvect2 (un,res,vtemp,s,vtemp,NV)

        ! ...  c3 = < s | s >
        call  vecnor (s,c3,NV) 

        beta = c3 / c1
        c1 = c3

        !  ...  p = s + (beta*p)    
        call  vsum2(s,p,beta,NV)

      ENDDO
 5000 continue

      IF(iflag==2) CALL vxg5(x,x,NV)

      CALL convert_f(x0,x,2)

      itest = 0
      return
      END

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

      return
      end
C----------------------------------------------------------------------C
      SUBROUTINE  convert_f(a,b,iflag)
C----------------------------------------------------------------------C
C     a : f_field
C     b : normal array
C     iflag = 1
C       b = a
C     iflag = 2
C       a = b
C----------------------------------------------------------------------C
      USE field_f
      INCLUDE '../INCLUDE/para_geometry' 
      PARAMETER( NDFV=2*3 ) 

      TYPE(f_field) a
      REAL*8 b(NV,NDFV,4)
      INTEGER, INTENT(IN)::iflag

      IF(iflag==1)  THEN

         DO id = 1, 4

           DO it = 1, NT
           DO iz = 1, NZ
           DO iy = 1, NY
           DO ix = 1, NX

              i = 1 + (ix-1) + (iy-1)*NX + (iz-1)*NX*NY
     &              + (it-1)*NX*NY*NZ

              b(i,1,id) = REAL ( a%f(1,ix,iy,iz,it,id) )
              b(i,2,id) = AIMAG( a%f(1,ix,iy,iz,it,id) )
              b(i,3,id) = REAL ( a%f(2,ix,iy,iz,it,id) )
              b(i,4,id) = AIMAG( a%f(2,ix,iy,iz,it,id) )
              b(i,5,id) = REAL ( a%f(3,ix,iy,iz,it,id) )
              b(i,6,id) = AIMAG( a%f(3,ix,iy,iz,it,id) )

           ENDDO
           ENDDO
           ENDDO
           ENDDO

         ENDDO

      ENDIF

      IF(iflag==2)  THEN

         DO id = 1, 4

           DO it = 1, NT
           DO iz = 1, NZ
           DO iy = 1, NY
           DO ix = 1, NX

              i = 1 + (ix-1) + (iy-1)*NX + (iz-1)*NX*NY
     &              + (it-1)*NX*NY*NZ

              a%f(1,ix,iy,iz,it,id) = CMPLX( b(i,1,id), b(i,2,id) )
              a%f(2,ix,iy,iz,it,id) = CMPLX( b(i,3,id), b(i,4,id) )
              a%f(3,ix,iy,iz,it,id) = CMPLX( b(i,5,id), b(i,6,id) )

           ENDDO
           ENDDO
           ENDDO
           ENDDO

         ENDDO

      ENDIF
      RETURN
      END
 
C----------------------------------------------------------------------C
      SUBROUTINE  fprep2 (INN,MSIGN)
C----------------------------------------------------------------------C
C   -  SITES IN THE LATTICE ARE IDENTIFIED BY ONE "ABSOLUTE" INDEX     C
C      (= CUMULATIVE INDEX 'IC'=1,...,NV)                              C
C   -  PREP PREPARES NEAREST NEIGHBOUR LIST 'INN' FOR A HYPERCUBIC     C
C      LATTICE WITH PBC                                                C
C----------------------------------------------------------------------C
      INCLUDE '../INCLUDE/para_geometry' 

      INTEGER,INTENT(OUT):: INN(NV,4,2), MSIGN(NV,4)

      INTEGER IBC(4)
c     DATA IBC/ +1,+1,+1,+1/ 
      DATA IBC/ +1,+1,+1,-1/ 
c     DATA IBC/ -1,-1,-1,-1/ 
C//   JF(I,L) = MOD(I+1  ,L)                          ! INDEX FORWARD
C//   JB(I,L) = MOD(I-1+L,L)                          ! INDEX BACKWARD
C//                                                   ! CUMULATIVE INDEX
      JF(I,L) = MOD(I+1  ,L) 
      JB(I,L) = MOD(I-1+L,L) 
      IC(IX,IY,IZ,IT,N1,N2,N3) = 1+IX+N1*(IY+N2*(IZ+N3*IT))

      WRITE(*,5)  IBC
5     FORMAT(/1X,'#####    BOUNDARY CONDITION FOR FERMION : ',4I3,/)

      DO 10 IT=0,NT-1
      DO 10 IZ=0,NZ-1
      DO 10 IY=0,NY-1
      DO 10 IX=0,NX-1
      I=IC(IX,IY,IZ,IT,NX,NY,NZ)
      INN(I,1,1) = IC (JF(IX,NX),IY       ,IZ       ,IT       ,NX,NY,NZ)
      INN(I,2,1) = IC (IX       ,JF(IY,NY),IZ       ,IT       ,NX,NY,NZ)
      INN(I,3,1) = IC (IX       ,IY       ,JF(IZ,NZ),IT       ,NX,NY,NZ)
      INN(I,4,1) = IC (IX       ,IY       ,IZ       ,JF(IT,NT),NX,NY,NZ)
      INN(I,1,2) = IC (JB(IX,NX),IY       ,IZ       ,IT       ,NX,NY,NZ)
      INN(I,2,2) = IC (IX       ,JB(IY,NY),IZ       ,IT       ,NX,NY,NZ)
      INN(I,3,2) = IC (IX       ,IY       ,JB(IZ,NZ),IT       ,NX,NY,NZ)
      INN(I,4,2) = IC (IX       ,IY       ,IZ       ,JB(IT,NT),NX,NY,NZ)
10    CONTINUE
C  'MSIGN(NV,4)' CONTROLS THE BOUNDARY CONDITION FOR FERMIONS.
C
C       I    I    I    I
C       I    I    I    I <-- IBC(MU)
C       +----+----+----+----           I MU
C       I    I    I    I               I
C       I    I    I    I               I
C       +----+----+----+----           +----->  NU
C       I    I    I    I  ^
C       .    .    .    .  I
C       .    .    .    .  IBC(NU)

      DO 60 MU = 1, 4
      DO 60 IS = 1, NV
   60 MSIGN(IS,MU) = 1

      IX = NX-1
      DO 71 IY = 0, NY-1
      DO 71 IZ = 0, NZ-1
      DO 71 IT = 0, NT-1
      IS = IC(IX,IY,IZ,IT,NX,NY,NZ)
   71 MSIGN(IS,1) = IBC(1)

      IY = NY-1
      DO 72 IZ = 0, NZ-1
      DO 72 IT = 0, NT-1
      DO 72 IX = 0, NX-1
      IS = IC(IX,IY,IZ,IT,NX,NY,NZ)
   72 MSIGN(IS,2) = IBC(2)
      IZ = NZ-1
      DO 73 IT = 0, NT-1
      DO 73 IX = 0, NX-1
      DO 73 IY = 0, NY-1
      IS = IC(IX,IY,IZ,IT,NX,NY,NZ)
   73 MSIGN(IS,3) = IBC(3)
      IT = NT-1
      DO 74 IX = 0, NX-1
      DO 74 IY = 0, NY-1
      DO 74 IZ = 0, NZ-1
      IS = IC(IX,IY,IZ,IT,NX,NY,NZ)
   74 MSIGN(IS,4) = IBC(4)

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE  unear2 (un,inn,msign)
C-----------------------------------------------------------------------
C                                                     +      ^
C UN(NC,1,L) = -K * U (NC,L,mu=1)  UN(NC,5,L) = -K * U (NC,L-1,MU=1) 
C                                                     +      ^
C UN(NC,2,L) = -K * U (NC,L,mu=2)  UN(NC,6,L) = -K * U (NC,L-2,MU=2) 
C                                                     +      ^
C UN(NC,3,L) = -K * U (NC,L,mu=3)  UN(NC,7,L) = -K * U (NC,L-3,MU=3) 
C                                                     +      ^
C UN(NC,4,L) = -K * U (NC,L,mu=4)  UN(NC,8,L) = -K * U (NC,L-4,MU=4) 
C
C                                           HERE 'K' STaNDS FOR 'HOP'
C-----------------------------------------------------------------------
      USE fpara
      USE field_g
      INCLUDE '../INCLUDE/para_geometry' 

      PARAMETER( NDF=18 )

      TYPE(g_field0) u(4)
      common/ config/ u

      INTEGER, INTENT(IN) :: inn(NV,4,2), msign(NV,4) 
      REAL*8,  INTENT(OUT):: un(NV,9,NDF)

      REAL*8 a, b
      REAL*8 u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,
     &       u13,u14,u15,u16,u17,u18

      Direction1 : DO mu = 1, 4

        DO it = 1, NT
        DO iz = 1, NZ
        DO iy = 1, NY
        DO ix = 1, NX

           l = 1 + (ix-1) + (iy-1)*NX + (iz-1)*(NX*NY) 
     &           + (it-1)*(NX*NY*NZ)

           u01 = REAL ( u(mu)%g(1,1,ix,iy,iz,it) )
           u02 = AIMAG( u(mu)%g(1,1,ix,iy,iz,it) )
           u03 = REAL ( u(mu)%g(1,2,ix,iy,iz,it) )
           u04 = AIMAG( u(mu)%g(1,2,ix,iy,iz,it) )
           u05 = REAL ( u(mu)%g(1,3,ix,iy,iz,it) )
           u06 = AIMAG( u(mu)%g(1,3,ix,iy,iz,it) )
           u07 = REAL ( u(mu)%g(2,1,ix,iy,iz,it) )
           u08 = AIMAG( u(mu)%g(2,1,ix,iy,iz,it) )
           u09 = REAL ( u(mu)%g(2,2,ix,iy,iz,it) )
           u10 = AIMAG( u(mu)%g(2,2,ix,iy,iz,it) )
           u11 = REAL ( u(mu)%g(2,3,ix,iy,iz,it) )
           u12 = AIMAG( u(mu)%g(2,3,ix,iy,iz,it) )
           u13 = REAL ( u(mu)%g(3,1,ix,iy,iz,it) )
           u14 = AIMAG( u(mu)%g(3,1,ix,iy,iz,it) )
           u15 = REAL ( u(mu)%g(3,2,ix,iy,iz,it) )
           u16 = AIMAG( u(mu)%g(3,2,ix,iy,iz,it) )
           u17 = REAL ( u(mu)%g(3,3,ix,iy,iz,it) )
           u18 = AIMAG( u(mu)%g(3,3,ix,iy,iz,it) )

           un(l,mu, 1) = MSIGN(L,mu) * (-HOP) * u01
           un(l,mu, 2) = MSIGN(L,mu) * (-HOP) * u02
           un(l,mu, 3) = MSIGN(L,mu) * (-HOP) * u07
           un(l,mu, 4) = MSIGN(L,mu) * (-HOP) * u08
           un(l,mu, 5) = MSIGN(L,mu) * (-HOP) * u13
           un(l,mu, 6) = MSIGN(L,mu) * (-HOP) * u14
           un(l,mu, 7) = MSIGN(L,mu) * (-HOP) * u03
           un(l,mu, 8) = MSIGN(L,mu) * (-HOP) * u04
           un(l,mu, 9) = MSIGN(L,mu) * (-HOP) * u09
           un(l,mu,10) = MSIGN(L,mu) * (-HOP) * u10
           un(l,mu,11) = MSIGN(L,mu) * (-HOP) * u15
           un(l,mu,12) = MSIGN(L,mu) * (-HOP) * u16
           un(l,mu,13) = MSIGN(L,mu) * (-HOP) * u05
           un(l,mu,14) = MSIGN(L,mu) * (-HOP) * u06
           un(l,mu,15) = MSIGN(L,mu) * (-HOP) * u11
           un(l,mu,16) = MSIGN(L,mu) * (-HOP) * u12
           un(l,mu,17) = MSIGN(L,mu) * (-HOP) * u17
           un(l,mu,18) = MSIGN(L,mu) * (-HOP) * u18

        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ENDDO Direction1

      Direction2 : DO mu = 1, 4

        mup4 = mu + 4

        DO it = 1, NT
        DO iz = 1, NZ
        DO iy = 1, NY
        DO ix = 1, NX

           it1 = it
           iz1 = iz
           iy1 = iy
           ix1 = ix
           IF(mu==4) it1 = it-1
           IF(mu==3) iz1 = iz-1
           IF(mu==2) iy1 = iy-1
           IF(mu==1) ix1 = ix-1


           l = 1 + (ix-1) + (iy-1)*NX + (iz-1)*(NX*NY) 
     &           + (it-1)*(NX*NY*NZ)

           u01 = REAL ( u(mu)%g(1,1,ix1,iy1,iz1,it1) )
           u02 = AIMAG( u(mu)%g(1,1,ix1,iy1,iz1,it1) )
           u03 = REAL ( u(mu)%g(1,2,ix1,iy1,iz1,it1) )
           u04 = AIMAG( u(mu)%g(1,2,ix1,iy1,iz1,it1) )
           u05 = REAL ( u(mu)%g(1,3,ix1,iy1,iz1,it1) )
           u06 = AIMAG( u(mu)%g(1,3,ix1,iy1,iz1,it1) )
           u07 = REAL ( u(mu)%g(2,1,ix1,iy1,iz1,it1) )
           u08 = AIMAG( u(mu)%g(2,1,ix1,iy1,iz1,it1) )
           u09 = REAL ( u(mu)%g(2,2,ix1,iy1,iz1,it1) )
           u10 = AIMAG( u(mu)%g(2,2,ix1,iy1,iz1,it1) )
           u11 = REAL ( u(mu)%g(2,3,ix1,iy1,iz1,it1) )
           u12 = AIMAG( u(mu)%g(2,3,ix1,iy1,iz1,it1) )
           u13 = REAL ( u(mu)%g(3,1,ix1,iy1,iz1,it1) )
           u14 = AIMAG( u(mu)%g(3,1,ix1,iy1,iz1,it1) )
           u15 = REAL ( u(mu)%g(3,2,ix1,iy1,iz1,it1) )
           u16 = AIMAG( u(mu)%g(3,2,ix1,iy1,iz1,it1) )
           u17 = REAL ( u(mu)%g(3,3,ix1,iy1,iz1,it1) )
           u18 = AIMAG( u(mu)%g(3,3,ix1,iy1,iz1,it1) )

           un(l,mup4, 1) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u01
           un(l,mup4, 2) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u02
           un(l,mup4, 3) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u07
           un(l,mup4, 4) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u08
           un(l,mup4, 5) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u13
           un(l,mup4, 6) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u14
           un(l,mup4, 7) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u03
           un(l,mup4, 8) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u04
           un(l,mup4, 9) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u09
           un(l,mup4,10) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u10
           un(l,mup4,11) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u15
           un(l,mup4,12) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u16
           un(l,mup4,13) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u05
           un(l,mup4,14) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u06
           un(l,mup4,15) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u11
           un(l,mup4,16) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u12
           un(l,mup4,17) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u17
           un(l,mup4,18) = MSIGN(INN(L,MU,2),MU)*(-HOP)*u18

         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ENDDO Direction2

      DO mu  =  5, 8
        DO l  =  1, NV

          un(l,mu,2 ) = -un(l,mu,2 )
          un(l,mu,10) = -un(l,mu,10)
          un(l,mu,18) = -un(l,mu,18)
          a = un(l,mu,3 )
          b = un(l,mu,4 )
          un(l,mu,3 ) =  un(l,mu,7 )
          un(l,mu,4 ) = -un(l,mu,8 )
          un(l,mu,7 ) = a
          un(l,mu,8 ) = -b
          a = un(l,mu,5 )
          b = un(l,mu,6 )
          un(l,mu,5 ) =  un(l,mu,13)
          un(l,mu,6 ) = -un(l,mu,14)
          un(l,mu,13) =  a
          un(l,mu,14) = -b
          a = un(l,mu,11)
          b = un(l,mu,12)
          un(l,mu,11) =  un(l,mu,15)
          un(l,mu,12) = -un(l,mu,16)
          un(l,mu,15) =  a
          un(l,mu,16) = -b

        ENDDO
      ENDDO

      DO k = 1, 18
        DO l = 1, NV
          un(l,9,k) = 0.0
        ENDDO
      ENDDO
 
      RETURN
      END


C----------------------------------------------------------------------C
      SUBROUTINE  clear  (v,n)
C----------------------------------------------------------------------C
      REAL*8, INTENT(OUT):: v(n)

      DO i = 1, n
        v(i) = 0.0d0
      ENDDO

      END

C----------------------------------------------------------------------
      SUBROUTINE  vdiff (v,v1,v2,nv)
C-----------------------------------------------------------------------
C     v = v1 - v2
C-----------------------------------------------------------------------
      PARAMETER( NDFV=2*3 )
      REAL*8,INTENT(IN) :: v1(nv*NDFV,4), V2(nv*NDFV,4)
      REAL*8,INTENT(OUT):: v(nv*NDFV,4)

      DO L = 1, nv*NDFV
        v(L,1) = v1(L,1) - v2(L,1)
        v(L,2) = v1(L,2) - v2(L,2)
        v(L,3) = v1(L,3) - v2(L,3)
        v(L,4) = v1(L,4) - v2(L,4)
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE  vecnor (vec,c,nv)
C-----------------------------------------------------------------------
C     C = < VEC I VEC >
C-----------------------------------------------------------------------
      PARAMETER( NDFV=2*3 )
      REAL*8,INTENT(IN) :: vec(nv*NDFV,4)
      REAL*8,INTENT(OUT):: c

      c = 0.0d0

      DO l = 1, nv*NDFV 
        c = c + vec(l,1)**2 + vec(l,2)**2
     *        + vec(l,3)**2 + vec(l,4)**2
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE  vnear (nv,x,xn,inn)
C-----------------------------------------------------------------------
C                                  ^                                  ^
C  XN(I,1,NC,MU) = X(L,NC,MU)  L=I+1; XN(I,5,NC,MU) = X(L,NC,MU)  L=I-1
C                                  ^                                  ^
C  XN(I,2,NC,MU) = X(L,NC,MU)  L=I+2; XN(I,6,NC,MU) = X(L,NC,MU)  L=I-2
C                                  ^                                  ^
C  XN(I,3,NC,MU) = X(L,NC,MU)  L=I+3; XN(I,7,NC,MU) = X(L,NC,MU)  L=I-3
C                                  ^                                  ^
C  XN(I,4,NC,MU) = X(L,NC,MU)  L=I+4; XN(I,8.NC,MU) = X(L,NC,MU)  L=I-4
C
C----------------------------------------------------------------------
      PARAMETER( NDFV=2*3 )  

      INTEGER,INTENT(IN):: inn(nv,4,2)
      REAL*8,INTENT(IN) :: x(nv,NDFV*4)
      REAL*8,INTENT(OUT):: xn(nv,9,NDFV*4)

      DO ncmu = 1, ndfv*4
      DO L = 1, NV
      xn(L,1,ncmu) = x(inn(L,1,1), ncmu)
      xn(L,2,ncmu) = x(inn(L,2,1), ncmu)
      xn(L,3,ncmu) = x(inn(L,3,1), ncmu)
      xn(L,4,ncmu) = x(inn(L,4,1), ncmu)
      xn(L,5,ncmu) = x(inn(L,1,2), ncmu)
      xn(L,6,ncmu) = x(inn(L,2,2), ncmu)
      xn(L,7,ncmu) = x(inn(L,3,2), ncmu)
      xn(L,8,ncmu) = x(inn(L,4,2), ncmu)
      ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE  vsum1 (vec1,vec2,ab,nv)
C-----------------------------------------------------------------------
C     VEC1 = VEC1 + AB * VEC2
C-----------------------------------------------------------------------
      PARAMETER( NDFV=2*3 )  
      INTEGER,INTENT(IN):: nv
      REAL*8,INTENT(IN):: ab
      REAL*8,INTENT(INOUT)::  VEC1(nv*NDFV,4)
      REAL*8,INTENT(IN)   ::  VEC2(nv*NDFV,4)

      DO L  = 1, nv*NDFV
        vec1(L,1) = vec1(L,1)  +  ab * vec2(L,1)
        vec1(L,2) = vec1(L,2)  +  ab * vec2(L,2)
        vec1(L,3) = vec1(L,3)  +  ab * vec2(L,3)
        vec1(L,4) = vec1(L,4)  +  ab * vec2(L,4)
      ENDDO

      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE  vsum2 (vec1,vec2,b,nv)
C-----------------------------------------------------------------------
C     vec2 = vec1 + b * vec2 
C-----------------------------------------------------------------------
      PARAMETER( NC=3, NDFV=2*NC  )     
      INTEGER,INTENT(IN):: nv
      REAL*8,INTENT(IN) :: b
      REAL*8,INTENT(IN) :: vec1(nv*NDFV,4)
      REAL*8,INTENT(OUT):: vec2(nv*NDFV,4)

      DO L = 1, nv*NDFV 
        vec2(L,1) = vec1(L,1)  +  b * vec2(L,1)
        vec2(L,2) = vec1(L,2)  +  b * vec2(L,2)
        vec2(L,3) = vec1(L,3)  +  b * vec2(L,3)
        vec2(L,4) = vec1(L,4)  +  b * vec2(L,4)
      ENDDO

      RETURN
      END

c---------------------------------------------------------------------
      SUBROUTINE  vxg5 (vin,vout,nv)
c-----------------------------------------------------------------------
c     vout = gamma_5 * vin
c----------------------------------------------------------------------
      PARAMETER( NC=3, NDFV=2*NC  ) 
      INTEGER,INTENT(IN):: nv
      REAL*8,INTENT(IN) :: vin(nv*NDFV,4)
      REAL*8,INTENT(OUT):: vout(nv*NDFV,4)

      DO ncl = 1, nv*NDFV
        vout(ncl,1) = -vin(ncl,1)
        vout(ncl,2) = -vin(ncl,2)
        vout(ncl,3) =  vin(ncl,3)
        vout(ncl,4) =  vin(ncl,4)
      ENDDO

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE  wxvect2 (un,vec,vecn,vout,uxvec,nv) 
c-----------------------------------------------------------------------
c     vout = ( gamma5*w(u,hop) ) * vec 
c-----------------------------------------------------------------------
      USE fpara

      PARAMETER( NDF=18, NDFV=2*3 )

      INTEGER,INTENT(IN)::nv
      REAL*8  un(nv,9,NDF), vec(nv,NDFV,4), vecn(nv,9,NDFV,4),
     *        vout(nv,NDFV,4), uxvec(nv,9,NDFV,4)

      call  uxvect3 (nv,un, vecn, uxvec)
      call  diraca  (vec,uxvec,vout,nv,r)

      IF( abs(Csw) > 0.001 ) THEN
        CALL vclover2(vout,vec)
      ENDIF

      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE  uxvect3 (nv,un,vv,uxvec)
c-----------------------------------------------------------------------
c     unxvv =  un * vv  in color space 
c-----------------------------------------------------------------------
      PARAMETER( NDF=18, NDFV=2*3 )
      REAL*8  un(nv*9,ndf), vv(nv*9,ndfv,4), uxvec(nv*9,ndfv,4)
      REAL*8  v1, v2, v3, v4, v5, v6

      do 5   nu  = 1, 4
      do 10  lmu = 1, nv*9

      v1 = vv(lmu,1,nu)
      v2 = vv(lmu,2,nu)
      v3 = vv(lmu,3,nu)
      v4 = vv(lmu,4,nu)
      v5 = vv(lmu,5,nu)
      v6 = vv(lmu,6,nu)

      vv(lmu,1 ,nu) =
     *            un(lmu,1 )*v1  - un(lmu,2 )*v2 
     *         +  un(lmu,7 )*v3  - un(lmu,8 )*v4 
     *         +  un(lmu,13)*v5  - un(lmu,14)*v6 

      vv(lmu,2 ,nu) =
     *            un(lmu,1 )*v2  + un(lmu,2 )*v1 
     *         +  un(lmu,7 )*v4  + un(lmu,8 )*v3 
     *         +  un(lmu,13)*v6  + un(lmu,14)*v5 

      vv(lmu,3 ,nu) =
     *            un(lmu,3 )*v1  - un(lmu,4 )*v2 
     *         +  un(lmu,9 )*v3  - un(lmu,10)*v4 
     *         +  un(lmu,15)*v5  - un(lmu,16)*v6 

      vv(lmu,4 ,nu) =
     *            un(lmu,3 )*v2  + un(lmu,4 )*v1 
     *         +  un(lmu,9 )*v4  + un(lmu,10)*v3 
     *         +  un(lmu,15)*v6  + un(lmu,16)*v5 

      vv(lmu,5 ,nu) =
     *            un(lmu,5 )*v1  - un(lmu,6 )*v2 
     *         +  un(lmu,11)*v3  - un(lmu,12)*v4 
     *         +  un(lmu,17)*v5  - un(lmu,18)*v6 

      vv(lmu,6 ,nu) =
     *            un(lmu,5 )*v2  + un(lmu,6 )*v1 
     *         +  un(lmu,11)*v4  + un(lmu,12)*v3 
     *         +  un(lmu,17)*v6  + un(lmu,18)*v5 

   10 continue

    5 continue

      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE  diraca (vec,x,vout,nv,r)
c-----------------------------------------------------------------------
c     vout(n) = gamma5*x(n) + sigma(mu=1,4)( (r - gamma_mu)*x(n+mu)
c                                    + (r + gamma_mu)*x(n-mu) )
c     x is uxvec in sub.wxvect2
c-----------------------------------------------------------------------
      PARAMETER( NDF=18, NC=3, NDFV=2*NC )

      REAL*8  vec(nv,ndfv,4), x(nv,9,ndfv,4), vout(nv,ndfv,4)
      REAL*8  r

      do 100  ic = 1, NC
      ii = 2*ic
      ir = ii - 1
      do 120  l = 1, nv
      vout(l,ir,1) = - vec(l,ir,1)
     *  - r*x(l,1,ir,1) +   x(l,1,ii,4) - r*x(l,5,ir,1) -   x(l,5,ii,4)
     *  - r*x(l,2,ir,1) -   x(l,2,ir,4) - r*x(l,6,ir,1) +   x(l,6,ir,4)
     *  - r*x(l,3,ir,1) +   x(l,3,ii,3) - r*x(l,7,ir,1) -   x(l,7,ii,3)
     *  - r*x(l,4,ir,1) -   x(l,4,ir,3) - r*x(l,8,ir,1) +   x(l,8,ir,3)

      vout(l,ii,1) = - vec(l,ii,1)
     *  - r*x(l,1,ii,1) -   x(l,1,ir,4) - r*x(l,5,ii,1) +   x(l,5,ir,4)
     *  - r*x(l,2,ii,1) -   x(l,2,ii,4) - r*x(l,6,ii,1) +   x(l,6,ii,4)
     *  - r*x(l,3,ii,1) -   x(l,3,ir,3) - r*x(l,7,ii,1) +   x(l,7,ir,3)
     *  - r*x(l,4,ii,1) -   x(l,4,ii,3) - r*x(l,8,ii,1) +   x(l,8,ii,3)
  120 continue

      do 140  l = 1, nv
      vout(l,ir,2) = - vec(l,ir,2)
     *  - r*x(l,1,ir,2) +   x(l,1,ii,3) - r*x(l,5,ir,2) -   x(l,5,ii,3)
     *  - r*x(l,2,ir,2) +   x(l,2,ir,3) - r*x(l,6,ir,2) -   x(l,6,ir,3)
     *  - r*x(l,3,ir,2) -   x(l,3,ii,4) - r*x(l,7,ir,2) +   x(l,7,ii,4)
     *  - r*x(l,4,IR,2) -   X(L,4,IR,4) - R*X(L,8,IR,2) +   X(L,8,IR,4)

      vout(l,ii,2) = - vec(l,ii,2)
     *  - r*x(l,1,ii,2) -   x(l,1,ir,3) - r*x(l,5,ii,2) +   x(l,5,ir,3)
     *  - r*x(l,2,ii,2) +   x(l,2,ii,3) - r*x(l,6,ii,2) -   x(l,6,ii,3)
     *  - r*x(l,3,ii,2) +   x(l,3,ir,4) - r*x(l,7,ii,2) -   x(l,7,ir,4)
     *  - r*x(l,4,ii,2) -   x(l,4,ii,4) - r*x(l,8,ii,2) +   x(l,8,ii,4)
  140 continue

      do 160  l = 1, nv
      vout(l,ir,3) = + vec(l,ir,3)
     *  +   x(l,1,ii,2) + r*x(l,1,ir,3) -   x(l,5,ii,2) + r*x(l,5,ir,3)
     *  -   x(l,2,ir,2) + r*x(l,2,ir,3) +   x(l,6,ir,2) + r*x(l,6,ir,3)
     *  +   x(l,3,ii,1) + r*x(l,3,ir,3) -   x(l,7,ii,1) + r*x(l,7,ir,3)
     *  +   x(l,4,ir,1) + r*x(l,4,ir,3) -   x(l,8,ir,1) + r*x(l,8,ir,3)

      vout(l,ii,3) = + vec(l,ii,3)
     *  -   x(l,1,ir,2) + r*x(l,1,ii,3) +   x(l,5,ir,2) + r*x(l,5,ii,3)
     *  -   x(l,2,ii,2) + r*x(l,2,ii,3) +   x(l,6,ii,2) + r*x(l,6,ii,3)
     *  -   x(l,3,ir,1) + r*x(l,3,ii,3) +   x(l,7,ir,1) + r*x(l,7,ii,3)
     *  +   x(l,4,ii,1) + r*x(l,4,ii,3) -   x(l,8,ii,1) + r*x(l,8,ii,3)
  160 continue

      do 180  l = 1, nv
      vout(l,ir,4) = + vec(l,ir,4)
     *  +   x(l,1,ii,1) + r*x(l,1,ir,4) -   x(l,5,ii,1) + r*x(l,5,ir,4)
     *  +   x(l,2,ir,1) + r*x(l,2,ir,4) -   x(l,6,ir,1) + r*x(l,6,ir,4)
     *  -   x(l,3,ii,2) + r*x(l,3,ir,4) +   x(l,7,ii,2) + r*x(l,7,ir,4)
     *  +   x(l,4,ir,2) + r*x(l,4,ir,4) -   x(l,8,ir,2) + r*x(l,8,ir,4)

      vout(l,ii,4) = + vec(l,ii,4)
     *  -   x(l,1,ir,1) + r*x(l,1,ii,4) +   x(l,5,ir,1) + r*x(l,5,ii,4)
     *  +   x(l,2,ii,1) + r*x(l,2,ii,4) -   x(l,6,ii,1) + r*x(l,6,ii,4)
     *  +   x(l,3,ir,2) + r*x(l,3,ii,4) -   x(l,7,ir,2) + r*x(l,7,ii,4)
     *  +   x(l,4,ii,2) + r*x(l,4,ii,4) -   x(l,8,ii,2) + r*x(l,8,ii,4)
  180 continue

  100 continue

      return
      end

C-----------------------------------------------------------------------
      SUBROUTINE  vclover2 (vec,X)
C-----------------------------------------------------------------------
C     vec = vec + gamma5*c*fmunu*(sigma_mu_nu)*X 
C-----------------------------------------------------------------------
      IMPLICIT REAl*8(a-h,o-z)

      include '../INCLUDE/para_geometry'
      PARAMETER( NDF=18, NDFV=2*3 )

*     TYPE(f_field), INTENT(IN)    :: X
*     TYPE(f_field), INTENT(INOUT) :: vec
      REAL*8, DIMENSION(NV,NDFV,4), INTENT(IN)       :: X
      REAL*8, DIMENSION(NV,NDFV,4), INTENT(INOUT)    :: vec

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

*     do it = 1, NT
*     do iz = 1, NZ
*     do iy = 1, NY
*     do ix = 1, NX
      DO i  = 1, NV

*           cx(1) = X%f(k2,ix,iy,iz,it,1) 
*           cx(2) = X%f(k2,ix,iy,iz,it,2)
*           cx(3) = X%f(k2,ix,iy,iz,it,3)
*           cx(4) = X%f(k2,ix,iy,iz,it,4)
            cx(1) = CMPLX( X(i,2*k2-1,1), X(i,2*k2,1) ) 
            cx(2) = CMPLX( X(i,2*k2-1,2), X(i,2*k2,2) ) 
            cx(3) = CMPLX( X(i,2*k2-1,3), X(i,2*k2,3) ) 
            cx(4) = CMPLX( X(i,2*k2-1,4), X(i,2*k2,4) ) 

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


      enddo Color2

      !  y = gamma5 * y
      DO i  = 1, NV
         y(k1,1,i) = - y(k1,1,i)
         y(k1,2,i) = - y(k1,2,i)
         y(k1,3,i) = + y(k1,3,i)
         y(k1,4,i) = + y(k1,4,i)
      ENDDO

      enddo Color1

      do mu = 1, 4
*     do it = 1, NT
*     do iz = 1, NZ
*     do iy = 1, NY
*     do ix = 1, NX
      DO i  = 1, NV
 
*        vec%f(1,ix,iy,iz,it,mu) =  
*    &             + vec%f(1,ix,iy,iz,it,mu) + y(1,mu,i)
*        vec%f(2,ix,iy,iz,it,mu) =  
*    &             + vec%f(2,ix,iy,iz,it,mu) + y(2,mu,i)
*        vec%f(3,ix,iy,iz,it,mu) =  
*    &             + vec%f(3,ix,iy,iz,it,mu) + y(3,mu,i)

         vec(i,1,mu) = vec(i,1,mu) + REAL(  y(1,mu,i) )
         vec(i,2,mu) = vec(i,2,mu) + AIMAG( y(1,mu,i) )
         vec(i,3,mu) = vec(i,3,mu) + REAL(  y(2,mu,i) )
         vec(i,4,mu) = vec(i,4,mu) + AIMAG( y(2,mu,i) )
         vec(i,5,mu) = vec(i,5,mu) + REAL(  y(3,mu,i) )
         vec(i,6,mu) = vec(i,6,mu) + AIMAG( y(3,mu,i) )

      enddo
      enddo

      return
      end
