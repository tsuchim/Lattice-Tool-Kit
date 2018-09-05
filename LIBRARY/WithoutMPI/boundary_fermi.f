c boundary_fermi.f
c fermion part of "boundary.f"
c----------------------------------------------------------------------c
      subroutine set_wing_f1(a)
c----------------------------------------------------------------------c
c     Fill wings of fermion vectors according to the boudary condition.c
c----------------------------------------------------------------------c
      USE field_f
      USE fpara

      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'
*      parameter( NC=3 )

      TYPE(f_field) a

      !  X-direction
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
           do k = 1, NC
                       a%f(k, 0,  iy, iz, it, ialpha)
     &       = BC(1) * a%f(k, NX, iy, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo

      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
           do k = 1, NC
                       a%f(k, NX+1,iy, iz, it, ialpha) 
     &       = BC(1) * a%f(k, 1,   iy, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
c
      !  Y-direction
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do ix = 1, NX 
           do k = 1, NC
                       a%f(k, ix, 0,  iz, it, ialpha)
     &       = BC(2) * a%f(k, ix, NY, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
c
c
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do ix = 1, NX
           do k = 1, NC
                       a%f(k, ix, NY+1,iz, it, ialpha)
     &       = BC(2) * a%f(k, ix, 1,   iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
 

        do ialpha = 1,4 

        !  Z-direction
        do it = 1, NT
        do iy = 1, NY
        do ix = 1, NX 

          a%f(1,ix,iy,0,it,ialpha) = BC(3) * a%f(1,ix,iy,NZ,it,ialpha)
          a%f(2,ix,iy,0,it,ialpha) = BC(3) * a%f(2,ix,iy,NZ,it,ialpha)
          a%f(3,ix,iy,0,it,ialpha) = BC(3) * a%f(3,ix,iy,NZ,it,ialpha)

          a%f(1,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(1,ix,iy,1,it,ialpha)
          a%f(2,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(2,ix,iy,1,it,ialpha)
          a%f(3,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(3,ix,iy,1,it,ialpha)

        enddo
        enddo
        enddo

        !  T-direction
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          a%f(1,ix,iy,iz,0,ialpha) = BC(4) * a%f(1,ix,iy,iz,NT,ialpha)
          a%f(2,ix,iy,iz,0,ialpha) = BC(4) * a%f(2,ix,iy,iz,NT,ialpha)
          a%f(3,ix,iy,iz,0,ialpha) = BC(4) * a%f(3,ix,iy,iz,NT,ialpha)

          a%f(1,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(1,ix,iy,iz,1,ialpha)
          a%f(2,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(2,ix,iy,iz,1,ialpha)
          a%f(3,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(3,ix,iy,iz,1,ialpha)

        enddo
        enddo
        enddo

      enddo

      return
      end

c----------------------------------------------------------------------c
      subroutine set_wing_f_eo0(a,iflg)
c----------------------------------------------------------------------c
c     Fill wings of fermion vectors according to the boudary condition.c
c----------------------------------------------------------------------c
      USE field_f
      USE fpara

      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'

      TYPE(f_field) a

c     ...  Starting Point   ...
      kld_sta = mod(1+1+iflg,2)
      krd_sta = mod(NX+1+iflg,2)
      klu_sta = mod(1+NY+iflg,2) 

      !  X-direction
c     .... iup ... 
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do iy = 1 + krd_sta, NY, 2
           do k = 1, NC
                       a%f(k, 0,  iy, iz, it, ialpha)
     &       = BC(1) * a%f(k, NX, iy, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo

      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do iy = 1+kld_sta,NY,2
           do k = 1, NC
                      a%f(k, NX+1,iy, iz, it, ialpha)
     &      = BC(1) * a%f(k, 1,   iy, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
 

      !  Y-direction
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do ix =  1+klu_sta, NX, 2
           do k = 1, NC
                       a%f(k, ix, 0,  iz, it, ialpha)
     &       = BC(2) * a%f(k, ix, NY, iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
 
      do ialpha = 1, 4
        do it = 1, NT
        do iz = 1, NZ
        do ix = 1+kld_sta, NX, 2
           do k = 1, NC
                       a%f(k, ix, NY+1,iz, it, ialpha)
     &       = BC(2) * a%f(k, ix, 1,   iz, it, ialpha)
           enddo
        enddo
        enddo
        enddo
      enddo
 

        do ialpha = 1,4 

        !  Z-direction
        do it = 1, NT
        do iy = 1, NY
        do ix = 1, NX 

          a%f(1,ix,iy,0,it,ialpha) = BC(3) * a%f(1,ix,iy,NZ,it,ialpha)
          a%f(2,ix,iy,0,it,ialpha) = BC(3) * a%f(2,ix,iy,NZ,it,ialpha)
          a%f(3,ix,iy,0,it,ialpha) = BC(3) * a%f(3,ix,iy,NZ,it,ialpha)

          a%f(1,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(1,ix,iy,1,it,ialpha)
          a%f(2,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(2,ix,iy,1,it,ialpha)
          a%f(3,ix,iy,NZ+1,it,ialpha) = BC(3) * a%f(3,ix,iy,1,it,ialpha)

        enddo
        enddo
        enddo

        !  T-direction
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          a%f(1,ix,iy,iz,0,ialpha) = BC(4) * a%f(1,ix,iy,iz,NT,ialpha)
          a%f(2,ix,iy,iz,0,ialpha) = BC(4) * a%f(2,ix,iy,iz,NT,ialpha)
          a%f(3,ix,iy,iz,0,ialpha) = BC(4) * a%f(3,ix,iy,iz,NT,ialpha)

          a%f(1,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(1,ix,iy,iz,1,ialpha)
          a%f(2,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(2,ix,iy,iz,1,ialpha)
          a%f(3,ix,iy,iz,NT+1,ialpha) = BC(4) * a%f(3,ix,iy,iz,1,ialpha)

        enddo
        enddo
        enddo

      enddo

      return
      end
