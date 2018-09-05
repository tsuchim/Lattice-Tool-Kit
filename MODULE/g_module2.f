C g_module2.f
C**********************************************************************C
      MODULE field_g
C**********************************************************************C
*     Module for gauge data structure for improved actions             *
*     for which we need eve/odd hyper-cube decomposision.              * 
C----------------------------------------------------------------------C
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

      parameter(MASK1=1)


      if( ibush < 16 )  then
         ieo_h = 0              ! Even Hyper-Cube
         ilocal = ibush         ! Local link number of Hyper Cube
      else
         ieo_h = 1              ! Odd Hyper-Cube
         ilocal = ibush - 16    ! Local link number of Hyper Cube
      endif

c     ...  ilocal = (ixl,iyl,izl,itl) = ixl + iyl*2 + izl*2^2 + itl*2^3
      ixl = IAND(ilocal,MASK1)
      iyl = IAND(ISHFT(ilocal,-1),MASK1)
      izl = IAND(ISHFT(ilocal,-2),MASK1)
      itl = IAND(ISHFT(ilocal,-3),MASK1)

      icum = 0

      do ith = 1, NTH

        it = 1 + 2*(ith-1) + itl

        do izh = 1, NZH

          iz = 1 + 2*(izh-1) + izl

          do iyh = 1, NYH

            iy = 1 + 2*(iyh-1) + iyl

            do ixh = 1+mod(1+ieo_h+iyh+izh+ith,2), NXH, 2

            ix = 1 + 2*(ixh-1) + ixl
       
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

      parameter(MASK1=1)

      if( ibush < 16 )  then
         ieo_h = 0              ! Even Hyper-Cube
         ilocal = ibush         ! Local link number of Hyper Cube
      else
         ieo_h = 1              ! Odd Hyper-Cube
         ilocal = ibush - 16    ! Local link number of Hyper Cube
      endif

c     ...  ilocal = (ixl,iyl,izl,itl) = ixl + iyl*2 + izl*2^2 + itl*2^3
      ixl = IAND(ilocal,MASK1)
      iyl = IAND(ISHFT(ilocal,-1),MASK1)
      izl = IAND(ISHFT(ilocal,-2),MASK1)
      itl = IAND(ISHFT(ilocal,-3),MASK1)

      icum = 0

      do ith = 1, NTH

        it = 1 + 2*(ith-1) + itl

        do izh = 1, NZH

          iz = 1 + 2*(izh-1) + izl

          do iyh = 1, NYH

            iy = 1 + 2*(iyh-1) + iyl

            do ixh = 1+mod(1+ieo_h+iyh+izh+ith,2), NXH, 2

            ix = 1 + 2*(ixh-1) + ixl
       
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

      parameter(MASK1=1)
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

      if( ibush < 16 )  then
         ieo_h = 0              ! Even Hyper-Cube
         ilocal = ibush         ! Local link number of Hyper Cube
      else
         ieo_h = 1              ! Odd Hyper-Cube
         ilocal = ibush - 16    ! Local link number of Hyper Cube
      endif

c     ...  ilocal = (ixl,iyl,izl,itl) = ixl + iyl*2 + izl*2^2 + itl*2^3
      ixl = IAND(ilocal,MASK1)
      iyl = IAND(ISHFT(ilocal,-1),MASK1)
      izl = IAND(ISHFT(ilocal,-2),MASK1)
      itl = IAND(ISHFT(ilocal,-3),MASK1)

      icum = 0

      do ith = 1, NTH

        it = 1 + 2*(ith-1) + itl
        it1 = it + idel(4)

        do izh = 1, NZH

          iz = 1 + 2*(izh-1) + izl
          iz1 = iz + idel(3)

          do iyh = 1, NYH

            iy = 1 + 2*(iyh-1) + iyl
            iy1 = iy + idel(2)

            do ixh = 1+mod(1+ieo_h+iyh+izh+ith,2), NXH, 2

               ix = 1 + 2*(ixh-1) + ixl
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
             
      jxl = mod(ixl+idel(1)+2,2)
      jyl = mod(iyl+idel(2)+2,2)
      jzl = mod(izl+idel(3)+2,2)
      jtl = mod(itl+idel(4)+2,2)
      b%parity = jxl + jyl*2 + jzl*(2**2) + jtl*(2**3) 
      b%direction = a%direction 
                  
      END FUNCTION
       
c----------------------------------------------------------------------c
      FUNCTION gshift2(dir,a) RESULT(b)
c----------------------------------------------------------------------c
      TYPE(ivec2),   INTENT(IN):: dir
      TYPE(g_field0), INTENT(IN):: a
      TYPE(g_field1) b

      parameter(MASK1=1)
      integer idel1(4), idel2(4)

      if((dir%x==0) .and. (dir%y==0))  then
         b = a
         return
      endif
             
      idel1(1) = 0; idel1(2) = 0; idel1(3) = 0; idel1(4) = 0
      idel2(1) = 0; idel2(2) = 0; idel2(3) = 0; idel2(4) = 0

      if( dir%x > 0 )  then
          idel1(+dir%x) = +1
      else if( dir%x < 0 ) then  
          idel1(-dir%x) = -1
      endif

      if( dir%y > 0 )  then
          idel2(+dir%y) = +1
      else if( dir%y < 0 ) then
          idel2(-dir%y) = -1
      endif
       
      if( ibush < 16 )  then
         ieo_h = 0              ! Even Hyper-Cube
         ilocal = ibush         ! Local link number of Hyper Cube
      else
         ieo_h = 1              ! Odd Hyper-Cube
         ilocal = ibush - 16    ! Local link number of Hyper Cube
      endif

c     ...  ilocal = (ixl,iyl,izl,itl) = ixl + iyl*2 + izl*2^2 + itl*2^3
      ixl = IAND(ilocal,MASK1)
      iyl = IAND(ISHFT(ilocal,-1),MASK1)
      izl = IAND(ISHFT(ilocal,-2),MASK1)
      itl = IAND(ISHFT(ilocal,-3),MASK1)
                      
      icum = 0
      do ith = 1, NTH                     ! Hyper-cube coordinate

        it = 1 + 2*(ith-1) + itl          ! site coordinate
        it1 = it  + idel1(4) + idel2(4)   ! shifted coordinate

        do izh = 1, NZH

          iz = 1 + 2*(izh-1) + izl
          iz1 = iz  + idel1(3) + idel2(3)

          do iyh = 1, NYH

            iy = 1 + 2*(iyh-1) + iyl
            iy1 = iy  + idel1(2) + idel2(2)

            do ixh = 1+mod(1+ieo_h+iyh+izh+ith,2), NXH, 2

              ix = 1 + 2*(ixh-1) + ixl
              ix1 = ix  + idel1(1) + idel2(1)

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

      jxl = mod(ixl+idel1(1)+idel2(1)+2,2)
      jyl = mod(iyl+idel1(2)+idel2(2)+2,2)
      jzl = mod(izl+idel1(3)+idel2(3)+2,2)
      jtl = mod(itl+idel1(4)+idel2(4)+2,2)

      b%parity = jxl + jyl*2 + jzl*(2**2) + jtl*(2**3) 
      b%direction = a%direction 
                  
      END FUNCTION
             
c----------------------------------------------------------------------c
      FUNCTION gshift3(dir,a) RESULT(b)
c----------------------------------------------------------------------c
      TYPE(ivec3),   INTENT(IN):: dir
      TYPE(g_field0), INTENT(IN):: a
      TYPE(g_field1) b

      parameter(MASK1=1)
      integer idel1(4), idel2(4), idel3(4)

      if((dir%x==0) .and. (dir%y==0) .and. (dir%z==0))  then
         b = a
         return
      endif
             
      idel1(1) = 0; idel1(2) = 0; idel1(3) = 0; idel1(4) = 0
      idel2(1) = 0; idel2(2) = 0; idel2(3) = 0; idel2(4) = 0
      idel3(1) = 0; idel3(2) = 0; idel3(3) = 0; idel3(4) = 0

      if( dir%x > 0 )  then
          idel1(+dir%x) = +1
      else if( dir%x < 0 ) then  
          idel1(-dir%x) = -1
      endif

      if( dir%y > 0 )  then
          idel2(+dir%y) = +1
      else if( dir%y < 0 ) then
          idel2(-dir%y) = -1
      endif
       
      if( dir%z > 0 )  then
          idel3(+dir%z) = +1
      else if( dir%z < 0 ) then
          idel3(-dir%z) = -1
      endif

      if( ibush < 16 )  then
         ieo_h = 0              ! Even Hyper-Cube
         ilocal = ibush         ! Local link number of Hyper Cube
      else
         ieo_h = 1              ! Odd Hyper-Cube
         ilocal = ibush - 16    ! Local link number of Hyper Cube
      endif

c     ...  ilocal = (ixl,iyl,izl,itl) = ixl + iyl*2 + izl*2^2 + itl*2^3
      ixl = IAND(ilocal,MASK1)
      iyl = IAND(ISHFT(ilocal,-1),MASK1)
      izl = IAND(ISHFT(ilocal,-2),MASK1)
      itl = IAND(ISHFT(ilocal,-3),MASK1)
                      
      icum = 0
      do ith = 1, NTH                                ! Hyper-cube coordinate

        it = 1 + 2*(ith-1) + itl                     ! site coordinate
        it1 = it  + idel1(4) + idel2(4) + idel3(4)   ! shifted coordinate

        do izh = 1, NZH

          iz = 1 + 2*(izh-1) + izl
          iz1 = iz  + idel1(3) + idel2(3) + idel3(3)

          do iyh = 1, NYH

            iy = 1 + 2*(iyh-1) + iyl
            iy1 = iy  + idel1(2) + idel2(2) + idel3(2)

            do ixh = 1+mod(1+ieo_h+iyh+izh+ith,2), NXH, 2

              ix = 1 + 2*(ixh-1) + ixl
              ix1 = ix  + idel1(1) + idel2(1) + idel3(1)

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

      jxl = mod(ixl+idel1(1)+idel2(1)+idel3(1)+4,2)
      jyl = mod(iyl+idel1(2)+idel2(2)+idel3(2)+4,2)
      jzl = mod(izl+idel1(3)+idel2(3)+idel3(3)+4,2)
      jtl = mod(itl+idel1(4)+idel2(4)+idel3(4)+4,2)

      b%parity = jxl + jyl*2 + jzl*(2**2) + jtl*(2**3) 
      b%direction = a%direction 
                  
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
