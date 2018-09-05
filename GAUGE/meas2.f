c----------------------------------------------------------------------c
      subroutine  CalPol(u,Pol,avePol)
c----------------------------------------------------------------------c
c
*     Measure Polyakov lines 
c
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0), INTENT(IN) :: u(4)
      COMPLEX*16, INTENT(OUT)    :: Pol(3,3,NX,NY,NZ), avePol

      COMPLEX*16 tmp1(3,3), tmp2(3,3)

      do mu = 1, 4
         call set_wing_g2(u(mu))
      enddo

      DO iz = 1, NZ
      DO iy = 1, NY
      DO ix = 1, NX

         Pol(1,1,ix,iy,iz) = U(4)%g(1,1,ix,iy,iz,1)
         Pol(1,2,ix,iy,iz) = U(4)%g(1,2,ix,iy,iz,1)
         Pol(1,3,ix,iy,iz) = U(4)%g(1,3,ix,iy,iz,1)
         Pol(2,1,ix,iy,iz) = U(4)%g(2,1,ix,iy,iz,1)
         Pol(2,2,ix,iy,iz) = U(4)%g(2,2,ix,iy,iz,1)
         Pol(2,3,ix,iy,iz) = U(4)%g(2,3,ix,iy,iz,1)
         Pol(3,1,ix,iy,iz) = U(4)%g(3,1,ix,iy,iz,1)
         Pol(3,2,ix,iy,iz) = U(4)%g(3,2,ix,iy,iz,1)
         Pol(3,3,ix,iy,iz) = U(4)%g(3,3,ix,iy,iz,1)

      ENDDO
      ENDDO
      ENDDO

      DO it = 2, NT

      DO iz = 1, NZ
      DO iy = 1, NY
      DO ix = 1, NX

         tmp1(1,1) = Pol(1,1,ix,iy,iz)
         tmp1(1,2) = Pol(1,2,ix,iy,iz) 
         tmp1(1,3) = Pol(1,3,ix,iy,iz)
         tmp1(2,1) = Pol(2,1,ix,iy,iz)
         tmp1(2,2) = Pol(2,2,ix,iy,iz)
         tmp1(2,3) = Pol(2,3,ix,iy,iz)
         tmp1(3,1) = Pol(3,1,ix,iy,iz)
         tmp1(3,2) = Pol(3,2,ix,iy,iz)
         tmp1(3,3) = Pol(3,3,ix,iy,iz)

         tmp2(1,1) = U(4)%g(1,1,ix,iy,iz,it)
         tmp2(1,2) = U(4)%g(1,2,ix,iy,iz,it)
         tmp2(1,3) = U(4)%g(1,3,ix,iy,iz,it)
         tmp2(2,1) = U(4)%g(2,1,ix,iy,iz,it)
         tmp2(2,2) = U(4)%g(2,2,ix,iy,iz,it)
         tmp2(2,3) = U(4)%g(2,3,ix,iy,iz,it)
         tmp2(3,1) = U(4)%g(3,1,ix,iy,iz,it)
         tmp2(3,2) = U(4)%g(3,2,ix,iy,iz,it)
         tmp2(3,3) = U(4)%g(3,3,ix,iy,iz,it)

         Pol(1,1,ix,iy,iz) 
     &      = tmp1(1,1)*tmp2(1,1) + tmp1(1,2)*tmp2(2,1) 
     &      + tmp1(1,3)*tmp2(3,1)
         Pol(1,2,ix,iy,iz) 
     &      = tmp1(1,1)*tmp2(1,2) + tmp1(1,2)*tmp2(2,2) 
     &      + tmp1(1,3)*tmp2(3,2)
         Pol(1,3,ix,iy,iz) 
     &      = tmp1(1,1)*tmp2(1,3) + tmp1(1,2)*tmp2(2,3) 
     &      + tmp1(1,3)*tmp2(3,3)

         Pol(2,1,ix,iy,iz) 
     &      = tmp1(2,1)*tmp2(1,1) + tmp1(2,2)*tmp2(2,1) 
     &      + tmp1(2,3)*tmp2(3,1)
         Pol(2,2,ix,iy,iz) 
     &      = tmp1(2,1)*tmp2(1,2) + tmp1(2,2)*tmp2(2,2) 
     &      + tmp1(2,3)*tmp2(3,2)
         Pol(2,3,ix,iy,iz) 
     &      = tmp1(2,1)*tmp2(1,3) + tmp1(2,2)*tmp2(2,3) 
     &      + tmp1(2,3)*tmp2(3,3)

         Pol(3,1,ix,iy,iz) 
     &      = tmp1(3,1)*tmp2(1,1) + tmp1(3,2)*tmp2(2,1) 
     &      + tmp1(3,3)*tmp2(3,1)
         Pol(3,2,ix,iy,iz) 
     &      = tmp1(3,1)*tmp2(1,2) + tmp1(3,2)*tmp2(2,2) 
     &      + tmp1(3,3)*tmp2(3,2)
         Pol(3,3,ix,iy,iz) 
     &      = tmp1(3,1)*tmp2(1,3) + tmp1(3,2)*tmp2(2,3) 
     &      + tmp1(3,3)*tmp2(3,3)


      ENDDO
      ENDDO
      ENDDO

      ENDDO

      avePol = 0.d0

      DO iz = 1, NZ
      DO iy = 1, NY
      DO ix = 1, NX

         avePol = avePol + ( Pol(1,1,ix,iy,iz) 
     &                   +   Pol(2,2,ix,iy,iz) 
     &                   +   Pol(3,3,ix,iy,iz) )

      ENDDO
      ENDDO
      ENDDO

      avePol = avePol/DBLE(NX*NY*NZ)

      RETURN
      END
