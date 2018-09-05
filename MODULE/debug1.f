C debug1.f
C**********************************************************************C
      MODULE debug1 
C**********************************************************************C
*     Module for gauge data structure for standard plaquette actions   *
*     for which we need eve/odd decomposision.                         * 
C----------------------------------------------------------------------C
      USE field_g
      USE hmc_mod

      include '../INCLUDE/para_geometry'
      PRIVATE  NX, NY, NZ, NT, NV, NVH
      PRIVATE  NBUSH
      PRIVATE  NXH, NYH, NZH, NTH
      PRIVATE  NDW


      TYPE adumy
         REAL*8, DIMENSION(8,NV,4) :: a
      END TYPE

      TYPE fdumy
         REAL*8, DIMENSION(6,4,NV) :: f
      END TYPE

      TYPE(adumy) atemp
      TYPE(fdumy) ftemp

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE dsubst1
      END INTERFACE

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE dsubst2
      END INTERFACE

      CONTAINS 
c......................................................................c
c     Definition of operators
c----------------------------------------------------------------------c
      SUBROUTINE  dsubst1(a,b)
c----------------------------------------------------------------------c
      TYPE(a_field), INTENT(OUT):: a(4)
      TYPE(adumy), INTENT(IN) :: b

      do nu = 1, 4
      do is = 1, NV
         a(nu)%a(1,is) = b%a(1,nu,is)
         a(nu)%a(2,is) = b%a(2,nu,is)
         a(nu)%a(3,is) = b%a(3,nu,is)
         a(nu)%a(4,is) = b%a(4,nu,is)
         a(nu)%a(5,is) = b%a(5,nu,is)
         a(nu)%a(6,is) = b%a(6,nu,is)
         a(nu)%a(7,is) = b%a(7,nu,is)
         a(nu)%a(8,is) = b%a(8,nu,is)
      enddo
      enddo

      END SUBROUTINE

c----------------------------------------------------------------------c
      SUBROUTINE  dsubst2(a,b)
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(OUT):: a
      TYPE(fdumy), INTENT(IN) :: b

      do ia = 1, 4
      is = 0
      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
      do ix = 1, NX
      is = is + 1
         a%f(1,ix,iy,iz,it,ia) = dcmplx( b%f(1,ia,is), b%f(2,ia,is) ) 
         a%f(2,ix,iy,iz,it,ia) = dcmplx( b%f(3,ia,is), b%f(4,ia,is) ) 
         a%f(3,ix,iy,iz,it,ia) = dcmplx( b%f(5,ia,is), b%f(6,ia,is) ) 
      enddo
      enddo
      enddo
      enddo
      enddo

      END SUBROUTINE
c----------------------------------------------------------------------c
      END MODULE


