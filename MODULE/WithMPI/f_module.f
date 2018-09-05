C**********************************************************************C
      MODULE field_f
C**********************************************************************C
c                                    written by AN on May 25, 2000 
c----------------------------------------------------------------------c
      USE field_g

      include '../INCLUDE/para_geometry'
      PRIVATE  NX, NY, NZ, NT, NV, NVH
      PRIVATE  NBUSH
      PRIVATE  NXH, NYH, NZH, NTH
      PRIVATE  NDW

      parameter( NDIM_F=NC*(NX+2)*(NY+2)*(NZ+2)*(NT+2)*4 )
 

      TYPE f_field
         SEQUENCE
         COMPLEX*16, DIMENSION(NC,0:NX+1,0:NY+1,0:NZ+1,0:NT+1,4) :: f   
                       ! (Color, x,     y,     z,     t, Dirac)
      END TYPE

      TYPE(g_field0) u0(4)

      TYPE(f_field)  temp, temp1, temp2, temp3
             
      TYPE(f_field)  zero_f
      data zero_f%f / NDIM_F*(0.d0,0.d0)/

      INTERFACE OPERATOR(+)
         MODULE PROCEDURE fadd
      END INTERFACE

      INTERFACE OPERATOR(-)
         MODULE PROCEDURE fsub
      END INTERFACE

      INTERFACE OPERATOR(*)
         MODULE PROCEDURE fmul   ! Define as the inner product
      END INTERFACE

      INTERFACE OPERATOR(*)
         MODULE PROCEDURE dmul   ! Define as scalar(double) * vector 
      END INTERFACE

      INTERFACE OPERATOR(*)
         MODULE PROCEDURE cmul   ! Define as scalar(complex) * vector 
      END INTERFACE

      INTERFACE OPERATOR(.fshift.)
         MODULE PROCEDURE fshift
      END INTERFACE

      INTERFACE OPERATOR(.fshiftB.)
         MODULE PROCEDURE fshiftB
      END INTERFACE


      CONTAINS 
c......................................................................c
c     Definition of operators
c----------------------------------------------------------------------c
      FUNCTION fadd(a,b) RESULT(c)
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(IN):: a, b
      TYPE(f_field) c
      COMPLEX*16  f1, f2, f3

      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          f1 = a%f(1,ix,iy,iz,it,mu) + b%f(1,ix,iy,iz,it,mu)
          f2 = a%f(2,ix,iy,iz,it,mu) + b%f(2,ix,iy,iz,it,mu)
          f3 = a%f(3,ix,iy,iz,it,mu) + b%f(3,ix,iy,iz,it,mu)

          c%f(1,ix,iy,iz,it,mu) = f1
          c%f(2,ix,iy,iz,it,mu) = f2
          c%f(3,ix,iy,iz,it,mu) = f3

        enddo
        enddo
        enddo
        enddo

      enddo
         
      END FUNCTION

c---------------------------------------------------------------------c
      FUNCTION fsub(a,b) RESULT(c)
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(IN):: a, b
      TYPE(f_field) c
      COMPLEX*16  f1, f2, f3

      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          f1 = a%f(1,ix,iy,iz,it,mu) - b%f(1,ix,iy,iz,it,mu)
          f2 = a%f(2,ix,iy,iz,it,mu) - b%f(2,ix,iy,iz,it,mu)
          f3 = a%f(3,ix,iy,iz,it,mu) - b%f(3,ix,iy,iz,it,mu)

          c%f(1,ix,iy,iz,it,mu) = f1
          c%f(2,ix,iy,iz,it,mu) = f2
          c%f(3,ix,iy,iz,it,mu) = f3

        enddo
        enddo
        enddo
        enddo

      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION fmul(a,b) RESULT(c)
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(IN):: a,b
      COMPLEX*16 c,csum
 
      include 'mpif.h'                                    ! MPI
      include '../INCLUDE/para.h'                         ! MPI
      include '../INCLUDE/para_geometry'
        
      c = (0.d0,0.d0)

      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

           c = c + conjg(a%f(1,ix,iy,iz,it,mu))*b%f(1,ix,iy,iz,it,mu)
     &           + conjg(a%f(2,ix,iy,iz,it,mu))*b%f(2,ix,iy,iz,it,mu)
     &           + conjg(a%f(3,ix,iy,iz,it,mu))*b%f(3,ix,iy,iz,it,mu)

        enddo
        enddo
        enddo
        enddo

      enddo

      call MPI_ALLREDUCE(c,csum,1,MPI_COMPLEX16,MPI_SUM,         ! MPI
     &                   MPI_COMM_WORLD,IERR)                    ! MPI
      c = csum                                                   ! MPI

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION dmul(s,a) RESULT(b)
c----------------------------------------------------------------------c
c     scalar (real*8) * Vector
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(IN) :: a
      TYPE(f_field)  b
      REAL*8, INTENT(IN) :: s
         
      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          b%f(1,ix,iy,iz,it,mu) = s * a%f(1,ix,iy,iz,it,mu)
          b%f(2,ix,iy,iz,it,mu) = s * a%f(2,ix,iy,iz,it,mu)
          b%f(3,ix,iy,iz,it,mu) = s * a%f(3,ix,iy,iz,it,mu)

        enddo
        enddo
        enddo
        enddo

      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION cmul(s,a) RESULT(b)
c----------------------------------------------------------------------c
c     scalar (complex*16) * Vector
c----------------------------------------------------------------------c
      TYPE(f_field), INTENT(IN) :: a
      TYPE(f_field)  b
      COMPLEX*16, INTENT(IN) :: s
         
      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          b%f(1,ix,iy,iz,it,mu) = s * a%f(1,ix,iy,iz,it,mu)
          b%f(2,ix,iy,iz,it,mu) = s * a%f(2,ix,iy,iz,it,mu)
          b%f(3,ix,iy,iz,it,mu) = s * a%f(3,ix,iy,iz,it,mu)

        enddo
        enddo
        enddo
        enddo

      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION fshift(mu,a) RESULT(b)
c----------------------------------------------------------------------c
c     mu = 0  b = a                                                    c
c        > 0  b(x) = (r-gamma_mu) U(x,mu) a(x+\mu)                     c
c        < 0  b(x) = (r+gamma_mu) U(x-(-\mu))_adj a(x-(-\mu))          c      
c----------------------------------------------------------------------c
      USE fpara

      TYPE(f_field), INTENT(IN):: a
      TYPE(f_field)  b
      INTEGER, INTENT(IN):: mu

      TYPE(g_field0) u0(4)
      common/ config/ u0

      COMPLEX*16  u11,u12,u13, u21,u22,u23, u31,u32,u33
      COMPLEX*16  e1, e2, e3, e4, f1, f2, f3, f4

      INTEGER idel(4)

      if( mu == 0 )  then
         b = a
      endif

      if( mu > 0 ) then

         idel(1)=0; idel(2)=0; idel(3)=0; idel(4)=0
         idel(mu) = 1

         ! ... Color multiplication
         do ialpha = 1, 4

         do it = 1, NT 
         it1 = it + idel(4)

         do iz = 1, NZ 
         iz1 = iz + idel(3)

         do iy = 1, NY 
         iy1 = iy + idel(2)

         do ix = 1, NX 
         ix1 = ix + idel(1)

            u11 = u0(mu)%g(1,1,ix,iy,iz,it)
            u12 = u0(mu)%g(1,2,ix,iy,iz,it)
            u13 = u0(mu)%g(1,3,ix,iy,iz,it)
            u21 = u0(mu)%g(2,1,ix,iy,iz,it)
            u22 = u0(mu)%g(2,2,ix,iy,iz,it)
            u23 = u0(mu)%g(2,3,ix,iy,iz,it)
            u31 = u0(mu)%g(3,1,ix,iy,iz,it)
            u32 = u0(mu)%g(3,2,ix,iy,iz,it)
            u33 = u0(mu)%g(3,3,ix,iy,iz,it)

             f1 = u11 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u12 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u13 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             f2 = u21 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u22 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u23 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             f3 = u31 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u32 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u33 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             b%f(1,ix,iy,iz,it,ialpha) = f1
             b%f(2,ix,iy,iz,it,ialpha) = f2
             b%f(3,ix,iy,iz,it,ialpha) = f3

         enddo
         enddo
         enddo
         enddo
         enddo

         ! ... Dirac multiplication
         do ic = 1, NC

           do it = 1, NT
           do iz = 1, NZ
           do iy = 1, NY
           do ix = 1, NX

             e1 = b%f(ic,ix,iy,iz,it,1)
             e2 = b%f(ic,ix,iy,iz,it,2)
             e3 = b%f(ic,ix,iy,iz,it,3)
             e4 = b%f(ic,ix,iy,iz,it,4)

             ! (gamma_mu * vector)
             f1 = rmg(1,1,mu)*e1 + rmg(1,2,mu)*e2 ! Corrected on Sep.14,2000
     &          + rmg(1,3,mu)*e3 + rmg(1,4,mu)*e4
             f2 = rmg(2,1,mu)*e1 + rmg(2,2,mu)*e2 
     &          + rmg(2,3,mu)*e3 + rmg(2,4,mu)*e4
             f3 = rmg(3,1,mu)*e1 + rmg(3,2,mu)*e2 
     &          + rmg(3,3,mu)*e3 + rmg(3,4,mu)*e4
             f4 = rmg(4,1,mu)*e1 + rmg(4,2,mu)*e2 
     &          + rmg(4,3,mu)*e3 + rmg(4,4,mu)*e4

             b%f(ic,ix,iy,iz,it,1) = f1
             b%f(ic,ix,iy,iz,it,2) = f2
             b%f(ic,ix,iy,iz,it,3) = f3
             b%f(ic,ix,iy,iz,it,4) = f4
            
           enddo
           enddo
           enddo
           enddo

         enddo

      endif    

      if( mu < 0 ) then

         idel(1)=0; idel(2)=0; idel(3)=0; idel(4)=0
         idel(-mu) = 1

         ! ... Color multiplication
         do ialpha = 1, 4

         do it = 1, NT
         it1 = it - idel(4)

         do iz = 1, NZ
         iz1 = iz - idel(3)

         do iy = 1, NY
         iy1 = iy - idel(2)

         do ix = 1, NX
         ix1 = ix - idel(1)

            u11 = conjg( u0(-mu)%g(1,1,ix1,iy1,iz1,it1) )
            u12 = conjg( u0(-mu)%g(2,1,ix1,iy1,iz1,it1) )
            u13 = conjg( u0(-mu)%g(3,1,ix1,iy1,iz1,it1) )
            u21 = conjg( u0(-mu)%g(1,2,ix1,iy1,iz1,it1) )
            u22 = conjg( u0(-mu)%g(2,2,ix1,iy1,iz1,it1) )
            u23 = conjg( u0(-mu)%g(3,2,ix1,iy1,iz1,it1) )
            u31 = conjg( u0(-mu)%g(1,3,ix1,iy1,iz1,it1) )
            u32 = conjg( u0(-mu)%g(2,3,ix1,iy1,iz1,it1) )
            u33 = conjg( u0(-mu)%g(3,3,ix1,iy1,iz1,it1) )

             f1 = u11 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u12 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u13 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             f2 = u21 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u22 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u23 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             f3 = u31 * a%f(1,ix1,iy1,iz1,it1,ialpha)
     &          + u32 * a%f(2,ix1,iy1,iz1,it1,ialpha)
     &          + u33 * a%f(3,ix1,iy1,iz1,it1,ialpha)

             b%f(1,ix,iy,iz,it,ialpha) = f1
             b%f(2,ix,iy,iz,it,ialpha) = f2
             b%f(3,ix,iy,iz,it,ialpha) = f3

         enddo
         enddo
         enddo
         enddo
         enddo

         ! (gamma_mu * vector)
         do ic = 1, NC

           do it = 1, NT
           do iz = 1, NZ
           do iy = 1, NY
           do ix = 1, NX

             e1 = b%f(ic,ix,iy,iz,it,1)
             e2 = b%f(ic,ix,iy,iz,it,2)
             e3 = b%f(ic,ix,iy,iz,it,3)
             e4 = b%f(ic,ix,iy,iz,it,4)

             ! (gamma_mu * vector)
*            f1 = rmg(1,1,-mu)*e1 + rmg(1,2,-mu)*e2 
*    &          + rmg(1,3,-mu)*e3 + rmg(1,4,-mu)*e4
*            f2 = rmg(2,1,-mu)*e1 + rmg(2,2,-mu)*e2 
*    &          + rmg(2,3,-mu)*e3 + rmg(2,4,-mu)*e4
*            f3 = rmg(3,1,-mu)*e1 + rmg(3,2,-mu)*e2 
*    &          + rmg(3,3,-mu)*e3 + rmg(3,4,-mu)*e4
*            f4 = rmg(4,1,-mu)*e1 + rmg(4,2,-mu)*e2 
*    &          + rmg(4,3,-mu)*e3 + rmg(4,4,-mu)*e4
             f1 = rpg(1,1,-mu)*e1 + rpg(1,2,-mu)*e2  ! Corrected on Sep.14,2000
     &          + rpg(1,3,-mu)*e3 + rpg(1,4,-mu)*e4
             f2 = rpg(2,1,-mu)*e1 + rpg(2,2,-mu)*e2 
     &          + rpg(2,3,-mu)*e3 + rpg(2,4,-mu)*e4
             f3 = rpg(3,1,-mu)*e1 + rpg(3,2,-mu)*e2 
     &          + rpg(3,3,-mu)*e3 + rpg(3,4,-mu)*e4
             f4 = rpg(4,1,-mu)*e1 + rpg(4,2,-mu)*e2 
     &          + rpg(4,3,-mu)*e3 + rpg(4,4,-mu)*e4

             b%f(ic,ix,iy,iz,it,1) = f1
             b%f(ic,ix,iy,iz,it,2) = f2
             b%f(ic,ix,iy,iz,it,3) = f3
             b%f(ic,ix,iy,iz,it,4) = f4

           enddo
           enddo
           enddo
           enddo

         enddo

      endif    

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION fshiftB(mu,a) RESULT(b)
c----------------------------------------------------------------------c
c     mu = 0  b = a                                                    c
c        > 0  b(x) = a_adj(x-\mu)*(r-gamma_mu) U(x-\mu,mu)
c        < 0  b(x) = a_adj(x+\mu)*(r+gamma_mu) U(x,mu)_adj
c----------------------------------------------------------------------c
      USE fpara

      TYPE(f_field), INTENT(IN):: a
      TYPE(f_field)  b
      INTEGER, INTENT(IN):: mu

      TYPE(g_field0) u0(4)
      common/ config/ u0

      COMPLEX*16  u11,u12,u13, u21,u22,u23, u31,u32,u33
      COMPLEX*16  e1, e2, e3, e4, f1, f2, f3, f4

      INTEGER idel(4)

      if( mu == 0 )  then
         b = a
      endif

      if( mu > 0 ) then

        write(*,*) "Sorry this case if not yet ready"
        write(*,*) "fshiftB  mu = ", mu
        stop

      endif    

      if( mu < 0 ) then

         idel(1)=0; idel(2)=0; idel(3)=0; idel(4)=0
         idel(-mu) = 1

         ! ... Color multiplication
         do ialpha = 1, 4

         do it = 1, NT
         it1 = it + idel(4)

         do iz = 1, NZ
         iz1 = iz + idel(3)

         do iy = 1, NY
         iy1 = iy + idel(2)

         do ix = 1, NX
         ix1 = ix + idel(1)

            u11 = conjg( u0(-mu)%g(1,1,ix,iy,iz,it) )
            u12 = conjg( u0(-mu)%g(2,1,ix,iy,iz,it) )
            u13 = conjg( u0(-mu)%g(3,1,ix,iy,iz,it) )
            u21 = conjg( u0(-mu)%g(1,2,ix,iy,iz,it) )
            u22 = conjg( u0(-mu)%g(2,2,ix,iy,iz,it) )
            u23 = conjg( u0(-mu)%g(3,2,ix,iy,iz,it) )
            u31 = conjg( u0(-mu)%g(1,3,ix,iy,iz,it) )
            u32 = conjg( u0(-mu)%g(2,3,ix,iy,iz,it) )
            u33 = conjg( u0(-mu)%g(3,3,ix,iy,iz,it) )

             f1 = conjg( a%f(1,ix1,iy1,iz1,it1,ialpha) ) * u11
     &          + conjg( a%f(2,ix1,iy1,iz1,it1,ialpha) ) * u21
     &          + conjg( a%f(3,ix1,iy1,iz1,it1,ialpha) ) * u31

             f2 = conjg( a%f(1,ix1,iy1,iz1,it1,ialpha) ) * u12
     &          + conjg( a%f(2,ix1,iy1,iz1,it1,ialpha) ) * u22
     &          + conjg( a%f(3,ix1,iy1,iz1,it1,ialpha) ) * u32

             f3 = conjg( a%f(1,ix1,iy1,iz1,it1,ialpha) ) * u13
     &          + conjg( a%f(2,ix1,iy1,iz1,it1,ialpha) ) * u23
     &          + conjg( a%f(3,ix1,iy1,iz1,it1,ialpha) ) * u33

             b%f(1,ix,iy,iz,it,ialpha) = f1
             b%f(2,ix,iy,iz,it,ialpha) = f2
             b%f(3,ix,iy,iz,it,ialpha) = f3

         enddo
         enddo
         enddo
         enddo
         enddo

         ! (vector * gamma_mu )
         do ic = 1, NC

           do it = 1, NT
           do iz = 1, NZ
           do iy = 1, NY
           do ix = 1, NX

             e1 = b%f(ic,ix,iy,iz,it,1)
             e2 = b%f(ic,ix,iy,iz,it,2)
             e3 = b%f(ic,ix,iy,iz,it,3)
             e4 = b%f(ic,ix,iy,iz,it,4)

             f1 = e1*rpg(1,1,-mu) + e2*rpg(2,1,-mu)
     &          + e3*rpg(3,1,-mu) + e4*rpg(4,1,-mu)
             f2 = e1*rpg(1,2,-mu) + e2*rpg(2,2,-mu)
     &          + e3*rpg(3,2,-mu) + e4*rpg(4,2,-mu)
             f3 = e1*rpg(1,3,-mu) + e2*rpg(2,3,-mu)
     &          + e3*rpg(3,3,-mu) + e4*rpg(4,3,-mu)
             f4 = e1*rpg(1,4,-mu) + e2*rpg(2,4,-mu) 
     &          + e3*rpg(3,4,-mu) + e4*rpg(4,4,-mu)

             b%f(ic,ix,iy,iz,it,1) = f1
             b%f(ic,ix,iy,iz,it,2) = f2
             b%f(ic,ix,iy,iz,it,3) = f3
             b%f(ic,ix,iy,iz,it,4) = f4

           enddo
           enddo
           enddo
           enddo

         enddo

      endif    

      END FUNCTION
c----------------------------------------------------------------------c
      END MODULE
