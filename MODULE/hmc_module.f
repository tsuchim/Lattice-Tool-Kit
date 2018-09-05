C hmc_module.f
C**********************************************************************C
      MODULE hmc_mod
C**********************************************************************C
C     module for HMC
C----------------------------------------------------------------------C
      USE field_g
      USE field_f
      PRIVATE u0
      PRIVATE temp, temp1, temp2, temp3 

      include '../INCLUDE/para_geometry'
      PRIVATE  NX, NY, NZ, NT, NV, NVH
      PRIVATE  NBUSH
      PRIVATE  NXH, NYH, NZH, NTH
      PRIVATE  NDW

      parameter( NDFALG=8 )


      TYPE a_field   !  Algebra field 
         REAL*8, DIMENSION(NDFALG,NV) :: a
      END TYPE

c     ...  common variables in HMC ...
      TYPE(a_field)  p(4)        ! Momentum PI conjugate to A_mu(x)
      TYPE(g_field0) usave(4)
      TYPE(f_field)  phi, eta, xi 
      REAL*8  beta, betamd, dtau, fac, eold, enew, eymold, eymnew
      INTEGER ntraj, nstep
      LOGICAL fermions

c     ...  Definitions of Operators  ...
      INTERFACE OPERATOR(+)      ! A1 + A2
         MODULE PROCEDURE aadd
      END INTERFACE

      INTERFACE OPERATOR(-)      ! A1 - A2
         MODULE PROCEDURE asub
      END INTERFACE

      INTERFACE OPERATOR(*)      ! (double scalar)*A
         MODULE PROCEDURE asmul
      END INTERFACE

      INTERFACE OPERATOR(*)      ! (double complex)*A
         MODULE PROCEDURE acmul
      END INTERFACE

      INTERFACE OPERATOR(*)      ! A1 * A2 
         MODULE PROCEDURE amul
      END INTERFACE

      CONTAINS
c......................................................................c
c     Definition of operators
c----------------------------------------------------------------------c
      FUNCTION aadd(x,y) RESULT(z)
c----------------------------------------------------------------------c
c     z = x + y
c----------------------------------------------------------------------c
      TYPE(a_field), INTENT(IN):: x, y
      TYPE(a_field) z

      do k = 1, NDFALG
      do i = 1, NV
        z%a(k,i) = x%a(k,i) + y%a(k,i)
      enddo
      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION asub(x,y) RESULT(z)
c----------------------------------------------------------------------c
c     z = x - y
c----------------------------------------------------------------------c
      TYPE(a_field), INTENT(IN):: x, y
      TYPE(a_field) z

      do k = 1, NDFALG
      do i = 1, NV
        z%a(k,i) = x%a(k,i) - y%a(k,i)
      enddo
      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION asmul(s,x) RESULT(y)
c----------------------------------------------------------------------c
c     y = s*x
c----------------------------------------------------------------------c
      REAL*8, INTENT(IN) :: s
      TYPE(a_field), INTENT(IN):: x
      TYPE(a_field) y

      do k = 1, NDFALG
      do i = 1, NV
        y%a(k,i) = s * x%a(k,i)
      enddo
      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION acmul(c,x) RESULT(y)
c----------------------------------------------------------------------c
c     y = c*x
c----------------------------------------------------------------------c
      COMPLEX*16, INTENT(IN) :: c
      TYPE(a_field), INTENT(IN):: x
      TYPE(a_field) y

      do k = 1, NDFALG
      do i = 1, NV
        y%a(k,i) = c * x%a(k,i)
      enddo
      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION amul(x,y) RESULT(s)
c----------------------------------------------------------------------c
c     s = x*y
c----------------------------------------------------------------------c
      TYPE(a_field), INTENT(IN):: x, y
      REAL*8 s,ssum

      s = 0.d0
      do k = 1, NDFALG
      do i = 1, NV
        s = s + x%a(k,i)*y%a(k,i) 
      enddo
      enddo


      END FUNCTION

c----------------------------------------------------------------------c
      END MODULE
