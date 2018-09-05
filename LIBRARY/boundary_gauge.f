c boundary.f
c---------------------------------------------------------------------c
      subroutine set_wing_g2 (u)
c---------------------------------------------------------------------c
c     Fix wings of the gauge field u                                  c
c     2-dimensional decomposition in X-Y plane                        c
c     NDW : width of fringe (1 or 2)                                  c
c---------------------------------------------------------------------c
      USE field_g

*     include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u

c     ! X direction
c     ....  Now we send data
      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
         do id = 1,NDW   
           do k1 = 1, NC
           do k2 = 1, NC
               u%g(k1,k2, -NDW+id,     iy, iz, it) 
     &       = u%g(k1,k2, NX+(id-NDW), iy, iz, it)
           enddo
           enddo
         enddo
      enddo
      enddo
      enddo
c
      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
         do id = 1,NDW
           do k1 = 1, NC
           do k2 = 1, NC
              u%g(k1, k2, NX+id, iy, iz, it) 
     &      = u%g(k1, k2,    id, iy, iz, it)
           enddo
           enddo
         enddo
      enddo
      enddo
      enddo


c     ! Y direction
c     ....  Now we send data
      do it = 1, NT
      do iz = 1, NZ
      do ix = -NDW+1, NX+NDW
        do id = 1,NDW   
           do k1 = 1, NC
           do k2 = 1, NC
               u%g(k1, k2, ix, -NDW+id,     iz, it) 
     &       = u%g(k1, k2, ix, NY+(id-NDW), iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
      enddo
c
      do it = 1, NT
      do iz = 1, NZ
      do ix = -NDW+1, NX+NDW
        do id = 1,NDW 
           do k1 = 1, NC
           do k2 = 1, NC
               u%g(k1, k2, ix, NY+id, iz, it) 
     &       = u%g(k1, k2, ix, id,    iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
      enddo


      ! Z-direction
      do id = 1,NDW
        do it = 1, NT
        do iy = -NDW+1, NY+NDW
        do ix = -NDW+1, NX+NDW
          do k1 = 1,NC
          do k2 = 1, NC  
           u%g(k1,k2,ix,iy,id-NDW,it)=u%g(k1,k2,ix,iy,NZ+(id-NDW),it)
           u%g(k1,k2,ix,iy,NZ+id,it)= u%g(k1,k2,ix,iy,id,it)
          enddo
          enddo 
        enddo
        enddo
        enddo
      enddo

      !  T-direction
      do id = 1,NDW
        do iz = -NDW+1, NZ+NDW
        do iy = -NDW+1, NY+NDW
        do ix = -NDW+1, NX+NDW
  
          do k1 = 1,NC
          do k2 = 1,NC
            u%g(k1,k2,ix,iy,iz,id-NDW)= u%g(k1,k2,ix,iy,iz,NT+(id-NDW))
            u%g(k1,k2,ix,iy,iz,NT+id) = u%g(k1,k2,ix,iy,iz,id)

          enddo
          enddo
        enddo
        enddo
        enddo
      enddo

      return
      end

c---------------------------------------------------------------------c
      subroutine set_wing_g_eo1 (u,iflg)
c---------------------------------------------------------------------c
c     Fix wings of the gauge field u                                  c
c     even/odd case                                                   c
c     odd ; iflg=1, even ; iflg=0                                     c
c     2-dimensional decomposition in X-Y plane.                       c
c---------------------------------------------------------------------c
      USE field_g

*     include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u

c     ...  Check  ...
      if(NDW.ne.1)  then
         write(*,*) "This routine (set_wing_g_eo1) is for even/odd."
         write(*,*) "Therefore NDW should be 1, but NDW=", NDW
         stop
      endif

c     ! X direction
c     ....  Now we send data

c     ...  iup ...
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(NX+1+iz+it+iflg,2)
        do iy = 1+ieo_sta, NY, 2
           do k1 = 1, NC
           do k2 = 1, NC
               u%g(k1,k2, 0,  iy, iz, it) 
     &       = u%g(k1,k2, NX, iy, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo

      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+1+iz+it+iflg,2)
        do iy = 1+ieo_sta, NY, 2
           do k1 = 1, NC
           do k2 = 1, NC
               u%g(k1, k2, NX+1, iy, iz, it)
     &       = u%g(k1, k2, 1,    iy, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo


c     ! Y direction

      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+NY+iz+it+iflg,2)
        do  ix = 1-ieo_sta, NX+1, 2
           do k1 = 1, NC
           do k2 = 1, NC
              u%g(k1, k2, ix, 0, iz, it) = u%g(k1, k2, ix, NY, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
c
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+1+iz+it+iflg,2)
        do ix = 1-ieo_sta, NX+1, 2
           do k1 = 1, NC
           do k2 = 1, NC
              u%g(k1, k2, ix, NY+1, iz, it) = u%g(k1, k2, ix, 1, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo


      ! Z-direction

      do it = 1, NT
      do iy = 0, NY+1

      ieo_sta = mod(0+iy+NZ+it+iflg,2)
      do ix = 0+ieo_sta, NX+1, 2
         do k1 = 1,NC
         do k2 = 1, NC  
            u%g(k1,k2,ix,iy,0,it)    = u%g(k1,k2,ix,iy,NZ,it)
         enddo
         enddo 
      enddo

      ieo_sta = mod(0+iy+1+it+iflg,2)
      do ix = 0+ieo_sta, NX+1, 2
         do k1 = 1,NC
         do k2 = 1, NC  
            u%g(k1,k2,ix,iy,NZ+1,it) = u%g(k1,k2,ix,iy,1,it)
         enddo
         enddo 
      enddo

      enddo
      enddo

      !  T-direction
      do iz = 0, NZ+1
      do iy = 0, NY+1

      ieo_sta = mod(0+iy+iz+NT+iflg,2)
      do ix = 0+ieo_sta, NX+1, 2

         do k1 = 1,NC
         do k2 = 1,NC
            u%g(k1,k2,ix,iy,iz,0)    = u%g(k1,k2,ix,iy,iz,NT)
         enddo
         enddo

      enddo

      ieo_sta = mod(0+iy+iz+1+iflg,2)
      do ix = 0+ieo_sta, NX+1, 2

         do k1 = 1,NC
         do k2 = 1,NC
            u%g(k1,k2,ix,iy,iz,NT+1) = u%g(k1,k2,ix,iy,iz,1)
         enddo
         enddo

      enddo

      enddo
      enddo

      return
      end

