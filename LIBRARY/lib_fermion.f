c-------------------------------------------------c
      function z2_ran_f(idumy) result(x)
c-------------------------------------------------c
c     Random number function for Z2 Noise
c-------------------------------------------------c
      USE field_f

      INTEGER, INTENT(IN):: idumy
      TYPE(f_field) x

      include '../INCLUDE/para_geometry'
      parameter( SQR2=1.41421356d0, SQR2I=1.d0/SQR2)

      do mu = 1, 4                ! Dirac
      do ic = 1, 3                ! Color 
         do it = 1, NT            ! Site
         do iz = 1, NZ
         do iy = 1, NY
         do ix = 1, NX

           iran1 = 2.0 * ranf(idum)
           iran2 = 2.0 * ranf(idum)
           kr = 2*iran1 - 1
           ki = 2*iran2 - 1
           xr = SQR2I * kr
           xi = SQR2I * ki
           x%f(ic,ix,iy,iz,it,mu) = CMPLX(xr,xi) 

         enddo
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end


c-------------------------------------------------c
      function gauss_ran_f(idumy) result(x)
c-------------------------------------------------c
c     Random number function for Gaussian  Noise
c-------------------------------------------------c
      USE field_f

      INTEGER, INTENT(IN):: idumy
      TYPE(f_field) x

      include '../INCLUDE/para_geometry'
      REAL*8 PI, PI2
      parameter( PI=3.14159265, PI2=2.d0*PI )

      do mu = 1, 4                ! Dirac
      do ic = 1, 3                ! Color 
         do it = 1, NT            ! Site
         do iz = 1, NZ
         do iy = 1, NY
         do ix = 1, NX

           v1 = sqrt( -log(ranf(idum)+1.e-10) )
           v2 = PI2 * ranf(idum)
c
           xr = v1 * cos(v2)
           xi = v1 * sin(v2)

           x%f(ic,ix,iy,iz,it,mu) = CMPLX(xr,xi) 

         enddo
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end
