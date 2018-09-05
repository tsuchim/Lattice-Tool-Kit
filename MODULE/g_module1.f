C g_module1.f
C**********************************************************************C
      MODULE field_g
C**********************************************************************C
*     Module for gauge data structure for standard plaquette actions   *
*     for which we need eve/odd decomposision.                         * 
c----------------------------------------------------------------------C
      include '../INCLUDE/para_geometry'
      PRIVATE  NX, NY, NZ, NT, NV, NVH
      PRIVATE  NBUSH
      PRIVATE  NXH, NYH, NZH, NTH
      PRIVATE  NDW

      parameter( NC=3 )
      INTEGER  ibush
      PUBLIC   ibush


      TYPE g_field0
         SEQUENCE
         COMPLEX*16, DIMENSION(NC,NC,-NDW+1:NX+NDW,-NDW+1:NY+NDW,
     &                               -NDW+1:NZ+NDW,-NDW+1:NT+NDW) :: g
         INTEGER parity, direction
      END TYPE

      TYPE g_field1
         SEQUENCE
         COMPLEX*16, DIMENSION(NC,NC,NV/NBUSH) :: g
         INTEGER parity, direction
      END TYPE

      TYPE ivec2
         INTEGER  x, y
      END TYPE
             
      TYPE ivec3
         INTEGER  x, y, z
      END TYPE

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE gsubst1
      END INTERFACE

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE gsubst2
      END INTERFACE

      INTERFACE OPERATOR(+)
         MODULE PROCEDURE gadd
      END INTERFACE

      INTERFACE OPERATOR(-)
         MODULE PROCEDURE gsub
      END INTERFACE

      INTERFACE OPERATOR(*)          ! A*B
         MODULE PROCEDURE prod0
      END INTERFACE

      INTERFACE OPERATOR(*)          ! (double scalar)*B
         MODULE PROCEDURE smul
      END INTERFACE

      INTERFACE OPERATOR(.prodAD.)    ! A*B_adj
         MODULE PROCEDURE prod1
      END INTERFACE

      INTERFACE OPERATOR(.ADprod.)    ! A_adj*B
         MODULE PROCEDURE prod2
      END INTERFACE

      INTERFACE OPERATOR(.gshift.)
         MODULE PROCEDURE gshift1
      END INTERFACE

      INTERFACE OPERATOR(.gshift.)
         MODULE PROCEDURE gshift2
      END INTERFACE

      INTERFACE OPERATOR(.gshift.)
         MODULE PROCEDURE gshift3
      END INTERFACE

      CONTAINS 
c......................................................................c
c     Definition of operators
c----------------------------------------------------------------------c
      SUBROUTINE  gsubst1(a,b)
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(OUT):: a
      TYPE(g_field0), INTENT(IN) :: b

      icum = 0
      do it = 1, NT
        do iz = 1, NZ
          do iy = 1, NY
            do ix = 1+mod(1+ibush+iy+iz+it,2), NX, 2
       
               icum = icum + 1
               a%g(1,1,icum) = b%g(1,1,ix,iy,iz,it) 
               a%g(2,1,icum) = b%g(2,1,ix,iy,iz,it) 
               a%g(3,1,icum) = b%g(3,1,ix,iy,iz,it) 
               a%g(1,2,icum) = b%g(1,2,ix,iy,iz,it) 
               a%g(2,2,icum) = b%g(2,2,ix,iy,iz,it) 
               a%g(3,2,icum) = b%g(3,2,ix,iy,iz,it) 
               a%g(1,3,icum) = b%g(1,3,ix,iy,iz,it) 
               a%g(2,3,icum) = b%g(2,3,ix,iy,iz,it) 
               a%g(3,3,icum) = b%g(3,3,ix,iy,iz,it) 

            enddo
          enddo
        enddo
      enddo

      a%parity = b%parity
      a%direction = b%direction

      END SUBROUTINE

c----------------------------------------------------------------------c
      SUBROUTINE  gsubst2(a,b)
c----------------------------------------------------------------------c
      TYPE(g_field0), INTENT(OUT):: a
      TYPE(g_field1), INTENT(IN) :: b

      icum = 0
      do it = 1, NT
        do iz = 1, NZ
          do iy = 1, NY
            do ix = 1+mod(1+ibush+iy+iz+it,2), NX, 2
       
               icum = icum + 1
               a%g(1,1,ix,iy,iz,it) = b%g(1,1,icum) 
               a%g(2,1,ix,iy,iz,it) = b%g(2,1,icum) 
               a%g(3,1,ix,iy,iz,it) = b%g(3,1,icum) 
               a%g(1,2,ix,iy,iz,it) = b%g(1,2,icum) 
               a%g(2,2,ix,iy,iz,it) = b%g(2,2,icum) 
               a%g(3,2,ix,iy,iz,it) = b%g(3,2,icum) 
               a%g(1,3,ix,iy,iz,it) = b%g(1,3,icum) 
               a%g(2,3,ix,iy,iz,it) = b%g(2,3,icum) 
               a%g(3,3,ix,iy,iz,it) = b%g(3,3,icum) 

            enddo
          enddo
        enddo
      enddo

      a%parity = b%parity
      a%direction = b%direction

      END SUBROUTINE

c----------------------------------------------------------------------c
      FUNCTION gadd(a,b) RESULT(c)
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(IN):: a, b
      TYPE(g_field1) c

      do i = 1, NV/NBUSH

         c%g( 1, 1, i) = a%g( 1, 1, i) + b%g( 1, 1, i)
         c%g( 1, 2, i) = a%g( 1, 2, i) + b%g( 1, 2, i)
         c%g( 1, 3, i) = a%g( 1, 3, i) + b%g( 1, 3, i)
         c%g( 2, 1, i) = a%g( 2, 1, i) + b%g( 2, 1, i)
         c%g( 2, 2, i) = a%g( 2, 2, i) + b%g( 2, 2, i)
         c%g( 2, 3, i) = a%g( 2, 3, i) + b%g( 2, 3, i)
         c%g( 3, 1, i) = a%g( 3, 1, i) + b%g( 3, 1, i)
         c%g( 3, 2, i) = a%g( 3, 2, i) + b%g( 3, 2, i)
         c%g( 3, 3, i) = a%g( 3, 3, i) + b%g( 3, 3, i)

      enddo
         
      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION gsub(a,b) RESULT(c)
c                                    written by SC 28 Mar 2000
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(IN):: a,b
      TYPE(g_field1) c
         
      do i = 1, NV/NBUSH

         c%g( 1, 1, i) = a%g( 1, 1, i) - b%g( 1, 1, i)
         c%g( 1, 2, i) = a%g( 1, 2, i) - b%g( 1, 2, i)
         c%g( 1, 3, i) = a%g( 1, 3, i) - b%g( 1, 3, i)
         c%g( 2, 1, i) = a%g( 2, 1, i) - b%g( 2, 1, i)
         c%g( 2, 2, i) = a%g( 2, 2, i) - b%g( 2, 2, i)
         c%g( 2, 3, i) = a%g( 2, 3, i) - b%g( 2, 3, i)
         c%g( 3, 1, i) = a%g( 3, 1, i) - b%g( 3, 1, i)
         c%g( 3, 2, i) = a%g( 3, 2, i) - b%g( 3, 2, i)
         c%g( 3, 3, i) = a%g( 3, 3, i) - b%g( 3, 3, i)

      enddo

      END FUNCTION
     
c----------------------------------------------------------------------c
      FUNCTION gshift1(nu,a) RESULT(b)
c----------------------------------------------------------------------c
      INTEGER, INTENT(IN):: nu
      TYPE(g_field0), INTENT(IN):: a
      TYPE(g_field1) b

      integer idel(4)

      if( nu==0 )  then
         b = a
         return
      endif

      idel(1) = 0; idel(2) = 0; idel(3) = 0; idel(4) = 0
      if( nu > 0 ) then
         idel(+nu) = +1
      else
         idel(-nu) = -1
      endif

      ieo = a%parity
      mu = a%direction
                         
*     ibush1 = 1 - ibush   ! Corrected on July 13,2000 by A.N.
      icum = 0
      do it = 1, NT
        it1 = it + idel(4)

        do iz = 1, NZ
          iz1 = iz + idel(3)

          do iy = 1, NY
             iy1 = iy + idel(2)

*            do ix = 1+mod(1+ibush1+iy+iz+it,2), NX, 2 ! Corrected on July
             do ix = 1+mod(1+ibush+iy+iz+it,2), NX, 2  ! 13,2000 by A.N.
               ix1 = ix + idel(1)

               icum = icum +1
               b%g( 1, 1, icum) = a%g( 1, 1, ix1,iy1,iz1,it1 )
               b%g( 1, 2, icum) = a%g( 1, 2, ix1,iy1,iz1,it1 )
               b%g( 1, 3, icum) = a%g( 1, 3, ix1,iy1,iz1,it1 )
               b%g( 2, 1, icum) = a%g( 2, 1, ix1,iy1,iz1,it1 )
               b%g( 2, 2, icum) = a%g( 2, 2, ix1,iy1,iz1,it1 )
               b%g( 2, 3, icum) = a%g( 2, 3, ix1,iy1,iz1,it1 )
               b%g( 3, 1, icum) = a%g( 3, 1, ix1,iy1,iz1,it1 )
               b%g( 3, 2, icum) = a%g( 3, 2, ix1,iy1,iz1,it1 )
               b%g( 3, 3, icum) = a%g( 3, 3, ix1,iy1,iz1,it1 )
         
             end do
          end do
        end do
      end do
             
      b%parity = 1-ieo
      b%direction = mu
                  
      END FUNCTION
       
c----------------------------------------------------------------------c
      FUNCTION gshift2(dir,a) RESULT(b)
c----------------------------------------------------------------------c
      TYPE(ivec2),   INTENT(IN):: dir
      TYPE(g_field0), INTENT(IN):: a
      TYPE(g_field1) b

      integer idel1(4), idel2(4)

c     ....  Stop for the exceptional case   ...
      if( (dir%x.ne.0) .and. (dir%y.ne.0) )  then
         if( dir%x==dir%y )  then
             write(*,*) "Sorry this case is not yet considered!"
             write(*,*) "dir : ", dir
             stop
         endif
      endif

      if((dir%x==0) .and. (dir%y==0))  then
         b = a
         return
      endif
             
      idel1(1) = 0; idel1(2) = 0; idel1(3) = 0; idel1(4) = 0
      idel2(1) = 0; idel2(2) = 0; idel2(3) = 0; idel2(4) = 0
      if( dir%x > 0 )  then
          idel1(+dir%x) = +1
      else if( dir%x < 0 ) then  ! Corrected on July,12,2000
          idel1(-dir%x) = -1
      endif
      if( dir%y > 0 )  then
*         idel1(+dir%y) = +1
          idel2(+dir%y) = +1   ! Corrected on July,11,2000
      else if( dir%y < 0 ) then  ! Corrected on July,12,2000
*         idel1(-dir%y) = -1
          idel2(-dir%y) = -1   ! Corrected on July,11,2000
      endif
       
      ieo = a%parity
      mu =  a%direction
                      
*     if( (dir%x.ne.0) .and. (dir%y.ne.0) )  then ! Corrected on 
*         ibush1 = ibush                          ! July 13,2000 
*     else                                        ! (A.N.)
*         ibush1 = 1 - ibush
*     endif 

      icum = 0
      do it = 1, NT 
        it1 = it  + idel1(4)
        it1 = it1 + idel2(4)

        do iz = 1, NZ
          iz1 = iz  + idel1(3)
          iz1 = iz1 + idel2(3)

          do iy = 1, NY
            iy1 = iy  + idel1(2)
            iy1 = iy1 + idel2(2)

*           do ix = 1+mod(1+ibush1+iy+iz+it,2), NX, 2 ! Corrected on
            do ix = 1+mod(1+ibush+iy+iz+it,2), NX, 2  ! July 13,2000(A.N.)
              ix1 = ix  + idel1(1)
              ix1 = ix1 + idel2(1)

              icum = icum +1
              b%g( 1, 1, icum) = a%g( 1, 1, ix1,iy1,iz1,it1 )
              b%g( 1, 2, icum) = a%g( 1, 2, ix1,iy1,iz1,it1 )
              b%g( 1, 3, icum) = a%g( 1, 3, ix1,iy1,iz1,it1 )
              b%g( 2, 1, icum) = a%g( 2, 1, ix1,iy1,iz1,it1 )
              b%g( 2, 2, icum) = a%g( 2, 2, ix1,iy1,iz1,it1 )
              b%g( 2, 3, icum) = a%g( 2, 3, ix1,iy1,iz1,it1 )
              b%g( 3, 1, icum) = a%g( 3, 1, ix1,iy1,iz1,it1 )
              b%g( 3, 2, icum) = a%g( 3, 2, ix1,iy1,iz1,it1 )
              b%g( 3, 3, icum) = a%g( 3, 3, ix1,iy1,iz1,it1 )
             
            enddo
          enddo
        enddo
      enddo

      if( (dir%x==0) .and. (dir%y==0) ) then
         jeo = ieo
      else if (  (dir%x==0) .or. (dir%y==0) ) then
         jeo = 1-ieo
      else
         jeo = ieo
      endif

      b%parity = jeo
      b%direction = mu
                  
      END FUNCTION
             
c----------------------------------------------------------------------c
      FUNCTION gshift3(dir,a) RESULT(b)
c----------------------------------------------------------------------c
      TYPE(ivec3),   INTENT(IN):: dir
      TYPE(g_field0), INTENT(IN):: a
      TYPE(g_field1) b

      write(*,*) "In g_module1.f, g_shift3 should not be used"
      b=a
      stop

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION prod0(a,b) RESULT(c)
c----------------------------------------------------------------------c
c     c = a*b                             written by AN 28 March 2000
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(IN):: a, b
      TYPE(g_field1) c 
      
      do i = 1, NV/NBUSH

         c%g(1,1,i) = a%g(1,1,i) * b%g(1,1,i)
     &              + a%g(1,2,i) * b%g(2,1,i)
     &              + a%g(1,3,i) * b%g(3,1,i)
         c%g(1,2,i) = a%g(1,1,i) * b%g(1,2,i)
     &              + a%g(1,2,i) * b%g(2,2,i)
     &              + a%g(1,3,i) * b%g(3,2,i)
         c%g(1,3,i) = a%g(1,1,i) * b%g(1,3,i)
     &              + a%g(1,2,i) * b%g(2,3,i)
     &              + a%g(1,3,i) * b%g(3,3,i)

         c%g(2,1,i) = a%g(2,1,i) * b%g(1,1,i)
     &              + a%g(2,2,i) * b%g(2,1,i)
     &              + a%g(2,3,i) * b%g(3,1,i)
         c%g(2,2,i) = a%g(2,1,i) * b%g(1,2,i)
     &              + a%g(2,2,i) * b%g(2,2,i)
     &              + a%g(2,3,i) * b%g(3,2,i)
         c%g(2,3,i) = a%g(2,1,i) * b%g(1,3,i)
     &              + a%g(2,2,i) * b%g(2,3,i)
     &              + a%g(2,3,i) * b%g(3,3,i)

         c%g(3,1,i) = a%g(3,1,i) * b%g(1,1,i)
     &              + a%g(3,2,i) * b%g(2,1,i)
     &              + a%g(3,3,i) * b%g(3,1,i)
         c%g(3,2,i) = a%g(3,1,i) * b%g(1,2,i)
     &              + a%g(3,2,i) * b%g(2,2,i)
     &              + a%g(3,3,i) * b%g(3,2,i)
         c%g(3,3,i) = a%g(3,1,i) * b%g(1,3,i)
     &              + a%g(3,2,i) * b%g(2,3,i)
     &              + a%g(3,3,i) * b%g(3,3,i)

      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION smul(s,a) RESULT(b)
c----------------------------------------------------------------------c
c     b = s*a                             
c----------------------------------------------------------------------c
      REAL*8, INTENT(IN)::s
      TYPE(g_field1), INTENT(IN):: a
      TYPE(g_field1) b 
      
      do i = 1, NV/NBUSH

         b%g(1,1,i) = s * a%g(1,1,i)
         b%g(1,2,i) = s * a%g(1,2,i)
         b%g(1,3,i) = s * a%g(1,3,i)
         b%g(2,1,i) = s * a%g(2,1,i)
         b%g(2,2,i) = s * a%g(2,2,i)
         b%g(2,3,i) = s * a%g(2,3,i)
         b%g(3,1,i) = s * a%g(3,1,i)
         b%g(3,2,i) = s * a%g(3,2,i)
         b%g(3,3,i) = s * a%g(3,3,i)

      enddo

      END FUNCTION

c----------------------------------------------------------------------c
      FUNCTION prod1(a,b) RESULT(c)  ! A*B_adj
c
c                                       written by SM 13 Apr 2000
c                                       checked by SC 13 Apr 2000
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(IN):: a, b
      TYPE(g_field1) c 

      do i = 1,NV/NBUSH
c
c     initialize result matrix c  
c
        c%g( 1, 1, i) = (0.D0,0.D0) 
        c%g( 1, 2, i) = (0.D0,0.D0) 
        c%g( 1, 3, i) = (0.D0,0.D0) 
        c%g( 2, 1, i) = (0.D0,0.D0) 
        c%g( 2, 2, i) = (0.D0,0.D0) 
        c%g( 2, 3, i) = (0.D0,0.D0) 
        c%g( 3, 1, i) = (0.D0,0.D0) 
        c%g( 3, 2, i) = (0.D0,0.D0) 
        c%g( 3, 3, i) = (0.D0,0.D0) 
c
      enddo
c
      do j = 1, 3             ! This is vector type ordering
       do i = 1, NV/NBUSH     !
c
         c%g( 1, 1, i) = c%g( 1, 1, i) + a%g( 1, j, i) 
     &                                 * conjg( b%g( 1, j, i))
         c%g( 1, 2, i) = c%g( 1, 2, i) + a%g( 1, j, i)  
     &                                 * conjg( b%g( 2, j, i))
         c%g( 1, 3, i) = c%g( 1, 3, i) + a%g( 1, j, i)  
     &                                 * conjg( b%g( 3, j, i))
         c%g( 2, 1, i) = c%g( 2, 1, i) + a%g( 2, j, i)  
     &                                 * conjg( b%g( 1, j, i))
         c%g( 2, 2, i) = c%g( 2, 2, i) + a%g( 2, j, i)  
     &                                 * conjg( b%g( 2, j, i))
         c%g( 2, 3, i) = c%g( 2, 3, i) + a%g( 2, j, i)  
     &                                 * conjg( b%g( 3, j, i))
         c%g( 3, 1, i) = c%g( 3, 1, i) + a%g( 3, j, i)  
     &                                 * conjg( b%g( 1, j, i))
         c%g( 3, 2, i) = c%g( 3, 2, i) + a%g( 3, j, i)  
     &                                 * conjg( b%g( 2, j, i))
         c%g( 3, 3, i) = c%g( 3, 3, i) + a%g( 3, j, i)  
     &                                 * conjg( b%g( 3, j, i))
c 
        enddo
      enddo


      END FUNCTION
c----------------------------------------------------------------------c
      FUNCTION prod2(a,b) RESULT(c)  ! A_adj*B
c----------------------------------------------------------------------c
c      c = a_adj * b                      written by AN 28 March 2000  c
c                                         checked by SM 06 April 2000  c
c                                         changed by AN 28 April 2000  c
c----------------------------------------------------------------------c
      TYPE(g_field1), INTENT(IN):: a, b
      TYPE(g_field1) c 

      do i = 1, NV/NBUSH

        c%g(1,1,i) = conjg(a%g(1,1,i)) * b%g(1,1,i)
     &             + conjg(a%g(2,1,i)) * b%g(2,1,i)
     &             + conjg(a%g(3,1,i)) * b%g(3,1,i)
        c%g(1,2,i) = conjg(a%g(1,1,i)) * b%g(1,2,i)
     &             + conjg(a%g(2,1,i)) * b%g(2,2,i)
     &             + conjg(a%g(3,1,i)) * b%g(3,2,i)
        c%g(1,3,i) = conjg(a%g(1,1,i)) * b%g(1,3,i)
     &             + conjg(a%g(2,1,i)) * b%g(2,3,i)
     &             + conjg(a%g(3,1,i)) * b%g(3,3,i)

        c%g(2,1,i) = conjg(a%g(1,2,i)) * b%g(1,1,i)
     &             + conjg(a%g(2,2,i)) * b%g(2,1,i)
     &             + conjg(a%g(3,2,i)) * b%g(3,1,i)
        c%g(2,2,i) = conjg(a%g(1,2,i)) * b%g(1,2,i)
     &             + conjg(a%g(2,2,i)) * b%g(2,2,i)
     &             + conjg(a%g(3,2,i)) * b%g(3,2,i)
        c%g(2,3,i) = conjg(a%g(1,2,i)) * b%g(1,3,i)
     &             + conjg(a%g(2,2,i)) * b%g(2,3,i)
     &             + conjg(a%g(3,2,i)) * b%g(3,3,i)

        c%g(3,1,i) = conjg(a%g(1,3,i)) * b%g(1,1,i)
     &             + conjg(a%g(2,3,i)) * b%g(2,1,i)
     &             + conjg(a%g(3,3,i)) * b%g(3,1,i)
        c%g(3,2,i) = conjg(a%g(1,3,i)) * b%g(1,2,i)
     &             + conjg(a%g(2,3,i)) * b%g(2,2,i)
     &             + conjg(a%g(3,3,i)) * b%g(3,2,i)
        c%g(3,3,i) = conjg(a%g(1,3,i)) * b%g(1,3,i)
     &             + conjg(a%g(2,3,i)) * b%g(2,3,i)
     &             + conjg(a%g(3,3,i)) * b%g(3,3,i)

      enddo
      END FUNCTION

c----------------------------------------------------------------------c
      END MODULE
