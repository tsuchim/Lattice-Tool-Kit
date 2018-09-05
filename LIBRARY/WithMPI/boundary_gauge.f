c boundary.f
c-----------------------------------------------------------------------c
      subroutine make_table (itable,iprocs,jprocs,myrank,
     &                       myranki,myrankj,mpi_proc_null)
c-----------------------------------------------------------------------c
      integer itable(-1:iprocs, -1:jprocs)

      irank = 0
      do i = 0, iprocs - 1
         do j = 0, jprocs - 1
            itable(i,j) = irank
                   if (myrank == irank) then
                       myranki = i
                       myrankj = j
                   endif
            irank = irank + 1
         enddo
      enddo

      ! periodic boundary condition

      do j = 0, jprocs-1
         itable(-1,j)     = itable(iprocs-1,j)
         itable(iprocs,j) = itable(0,j)
      enddo

      do i = 0, iprocs-1
         itable(i,-1)     = itable(i,jprocs-1)
         itable(i,jprocs) = itable(i,0)
      enddo

      return
      end

c---------------------------------------------------------------------c
      subroutine set_wing_g2 (u)
c---------------------------------------------------------------------c
c     Fix wings of the gauge field u                                  c
c     2-dimensional decomposition in X-Y plane                        c
c     NDW : width of fringe (1 or 2)                                  c
c---------------------------------------------------------------------c
      USE field_g

      include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u
*     COMPLEX*16 u(NC*NC,0:NX+1,0:NY+1,0:NZ+1,0:NT+1)

      COMPLEX*16 works1(NC,NC,NDW*NY*NZ*NT),
     &           works2(NC,NC,NDW*NY*NZ*NT), 
     &           works3(NC,NC,NDW*(NX+2*NDW)*NZ*NT),
     &           works4(NC,NC,NDW*(NX+2*NDW)*NZ*NT), 
     &           workr1(NC,NC,NDW*NY*NZ*NT),
     &           workr2(NC,NC,NDW*NY*NZ*NT),
     &           workr3(NC,NC,NDW*(NX+2*NDW)*NZ*NT),
     &           workr4(NC,NC,NDW*(NX+2*NDW)*NZ*NT)
      integer istatus(mpi_status_size)

c     ...  Neibouring PE's   ...
      iup   = itable(myranki+1, myrankj)
      idown = itable(myranki-1, myrankj)
      jup   = itable(myranki  , myrankj+1)
      jdown = itable(myranki  , myrankj-1)
c

c     ! X direction
c     ....  Now we send data
      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
          do id = 1,NDW   
          icum = icum + 1
          do k1 = 1, NC; do k2 = 1, NC
             works1(k1,k2, icum) = u%g(k1,k2, NX+(id-NDW), iy, iz, it)
           enddo; enddo
           enddo
        enddo
      enddo; enddo

      nsend = NC*NC*NDW*NY*NZ*NT  
      call mpi_isend(works1(1,1,1), nsend, mpi_complex16, 
     &               iup, 1, mpi_comm_world, jsend1, ierr)
c
      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
           do id = 1,NDW
           icum = icum + 1
           do k1 = 1, NC; do k2 = 1, NC
             works2(k1, k2, icum) = u%g(k1, k2, id, iy, iz, it)
           enddo; enddo
           enddo
        enddo
      enddo; enddo
c
      nsend = NC*NC*NDW*NY*NZ*NT
      call mpi_isend(works2(1,1,1), nsend, mpi_complex16, 
     &               idown, 1, mpi_comm_world, jsend2, ierr)
c
c     .....  Now we receive data   ...
      nrecv = NC*NC*NDW*NY*NZ*NT
      call mpi_irecv(workr1(1,1,1), nrecv, mpi_complex16, 
     &               idown, 1, mpi_comm_world, jrecv1, ierr)

      nrecv = NC*NC*NDW*NY*NZ*NT
      call mpi_irecv(workr2(1,1,1), nrecv, mpi_complex16, 
     &               iup, 1, mpi_comm_world, jrecv2, ierr)

      call mpi_wait(jsend1, istatus, ierr)
      call mpi_wait(jsend2, istatus, ierr)
      call mpi_wait(jrecv1, istatus, ierr)
      call mpi_wait(jrecv2, istatus, ierr)

      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
          do id = 1,NDW
          icum = icum + 1   
          do k1 = 1, NC; do k2 = 1, NC
            u%g(k1, k2, NX+id, iy, iz, it) = workr2(k1, k2, icum)
          enddo; enddo
          enddo
        enddo
      enddo; enddo

      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
          do id = 1,NDW
          icum = icum + 1 
          do k1 = 1, NC; do k2 = 1, NC
            u%g(k1, k2, -NDW+id, iy, iz, it) = workr1(k1, k2,  icum)
          enddo; enddo
          enddo
        enddo
      enddo; enddo

c     ! Y direction
c     ....  Now we send data
      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do ix = -NDW+1, NX+NDW
        do id = 1,NDW   
        icum = icum + 1
           do k1 = 1, NC; do k2 = 1, NC
             works3(k1, k2, icum) = u%g(k1, k2, ix, NY+(id-NDW), iz, it)
           enddo; enddo
           enddo
        enddo
      enddo; enddo
c
      nsend = NC*NC*NDW*(NX+2*NDW)*NZ*NT
      call mpi_isend(works3(1,1,1), nsend, mpi_complex16, 
     &               jup, 1, mpi_comm_world, isend1, ierr)
c
      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do ix = -NDW+1, NX+NDW
          do id = 1,NDW 
          icum = icum + 1
           do k1 = 1, NC; do k2 = 1, NC
             works4(k1, k2, icum) = u%g(k1, k2, ix, id, iz, it)
             enddo;enddo
             enddo
          enddo
      enddo; enddo
c
      nsend = NC*NC*NDW*(NX+2*NDW)*NZ*NT
      call mpi_isend(works4(1,1,1), nsend, mpi_complex16, 
     &               jdown, 1, mpi_comm_world, isend2, ierr)

c     .....  Now we receive data   ...
      nrecv = NC*NC*NDW*(NX+2*NDW)*NZ*NT
      call mpi_irecv(workr3(1,1,1), nrecv, mpi_complex16, 
     &               jdown, 1, mpi_comm_world, irecv1, ierr)

      nrecv = NC*NC*NDW*(NX+2*NDW)*NZ*NT
      call mpi_irecv(workr4(1,1,1), nrecv, mpi_complex16, 
     &               jup, 1, mpi_comm_world, irecv2, ierr)

      call mpi_wait(isend1, istatus, ierr)
      call mpi_wait(isend2, istatus, ierr)
      call mpi_wait(irecv1, istatus, ierr)
      call mpi_wait(irecv2, istatus, ierr)

      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do ix = -NDW+1, NX+NDW
         do id = 1,NDW   
         icum = icum + 1
          do k1 = 1, NC; do k2 = 1, NC
            u%g(k1, k2, ix, NY+id, iz, it) = workr4(k1, k2,  icum)
          enddo; enddo
          enddo
        enddo
      enddo; enddo

      icum = 0
      do it = 1, NT; do iz = 1, NZ
        do ix = -NDW+1, NX+NDW
        do id = 1,NDW
        icum = icum + 1
          do k1 = 1, NC; do k2 = 1, NC
            u%g(k1, k2, ix, -NDW+id, iz, it) = workr3(k1, k2, icum)
          enddo; enddo
          enddo
        enddo
      enddo; enddo

      ! Z-direction
      do id = 1,NDW
      do it = 1, NT
      do iy = -NDW+1, NY+NDW
      do ix = -NDW+1, NX+NDW
      do k1 = 1,NC; do k2 = 1, NC  

        u%g(k1,k2,ix,iy,id-NDW,it)=u%g(k1,k2,ix,iy,NZ+(id-NDW),it)
        u%g(k1,k2,ix,iy,NZ+id,it)= u%g(k1,k2,ix,iy,id,it)

      enddo; enddo 
      enddo
      enddo
      enddo
      enddo

      !  T-direction
      do id = 1,NDW
      do iz = -NDW+1, NZ+NDW
      do iy = -NDW+1, NY+NDW
      do ix = -NDW+1, NX+NDW
      do k1 = 1,NC; do k2 = 1,NC

        u%g(k1,k2,ix,iy,iz,id-NDW)= u%g(k1,k2,ix,iy,iz,NT+(id-NDW))
        u%g(k1,k2,ix,iy,iz,NT+id) = u%g(k1,k2,ix,iy,iz,id)

      enddo; enddo
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

      include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u

      COMPLEX*16 bufrs(NC,NC,NY*NZ*NT),bufls(NC,NC,NY*NZ*NT), 
     &           bufrr(NC,NC,NX*NZ*NT),buflr(NC,NC,NX*NZ*NT), 
     &           bufus(NC,NC,NY*NZ*NT),bufds(NC,NC,NY*NZ*NT),
     &           bufur(NC,NC,NX*NZ*NT),bufdr(NC,NC,NX*NZ*NT)
      integer istatus(mpi_status_size)


c     ...  Check  ...
      if(NDW.ne.1)  then
         write(*,*) "This routine (set_wing_g_eo1) is for even/odd."
         write(*,*) "Therefore NDW should be 1, but NDW=", NDW
         stop
      endif

c     ...  Neibouring PE's   ...
      iup   = itable(myranki+1, myrankj)
      idown = itable(myranki-1, myrankj)
      jup   = itable(myranki  , myrankj+1)
      jdown = itable(myranki  , myrankj-1)
c
c     ! X direction
c     ....  Now we send data

c     ...  iup ...
      nrs = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(NX+1+iz+it+iflg,2)
        do iy = 1+ieo_sta, NY, 2
        nrs = nrs + 1
           do k1 = 1, NC
           do k2 = 1, NC
             bufrs(k1,k2, nrs) = u%g(k1,k2, NX, iy, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo

      nsend = NC*NC*nrs
      call mpi_isend(bufrs(1,1,1), nsend, mpi_complex16, 
     &               iup, 1, mpi_comm_world, jsend1, ierr)
c
c     ... idown ..
      nls = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+1+iz+it+iflg,2)
        do iy = 1+ieo_sta, NY, 2
        nls = nls + 1
           do k1 = 1, NC
           do k2 = 1, NC
             bufls(k1, k2, nls) = u%g(k1, k2, 1, iy, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
c
      nsend = NC*NC*nls
      call mpi_isend(bufls(1,1,1), nsend, mpi_complex16, 
     &               idown, 1, mpi_comm_world, jsend2, ierr)
c
c     .....  Now we receive data   ...
      nrecv = NC*NC*NY*NZ*NT/2
      call mpi_irecv(buflr(1,1,1), nrecv, mpi_complex16, 
     &               idown, 1, mpi_comm_world, jrecv1, ierr)

      nrecv = NC*NC*NY*NZ*NT/2
      call mpi_irecv(bufrr(1,1,1), nrecv, mpi_complex16, 
     &               iup, 1, mpi_comm_world, jrecv2, ierr)

      call mpi_wait(jsend1, istatus, ierr)
      call mpi_wait(jsend2, istatus, ierr)
      call mpi_wait(jrecv1, istatus, ierr)
      call mpi_wait(jrecv2, istatus, ierr)
c
c     ... from left ...
      nlr = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(0+1+iz+it+iflg,2)
        do iy =1+ieo_sta, NY, 2
        nlr = nlr + 1   
          do k1 = 1, NC
          do k2 = 1, NC
            u%g(k1, k2, 0, iy, iz, it) = buflr(k1, k2, nlr)
          enddo
          enddo
        enddo
      enddo
      enddo

c     ... from right ...
      nrr = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod((NX+1)+1+iz+it+iflg,2)
        do iy = 1+ieo_sta, NY, 2
          nrr = nrr + 1 
          do k1 = 1, NC
          do k2 = 1, NC
            u%g(k1, k2, NX+1, iy, iz, it) = bufrr(k1, k2,  nrr)
          enddo
          enddo
        enddo
      enddo
      enddo

c     ! Y direction
c     ....  Now we send data
c     ... jup ...
      nus = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+NY+iz+it+iflg,2)
        do  ix = 1-ieo_sta, NX+1, 2
        nus = nus + 1
           do k1 = 1, NC
           do k2 = 1, NC
             bufus(k1, k2, nus) = u%g(k1, k2, ix, NY, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
c
      nsend = NC*NC*nus
      call mpi_isend(bufus(1,1,1), nsend, mpi_complex16, 
     &               jup, 1, mpi_comm_world, isend1, ierr)
c
c     ... jdown ...
      nds = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+1+iz+it+iflg,2)
        do ix = 1-ieo_sta, NX+1, 2
          nds = nds + 1
           do k1 = 1, NC
           do k2 = 1, NC
             bufds(k1, k2, nds) = u%g(k1, k2, ix, 1, iz, it)
           enddo
           enddo
        enddo
      enddo
      enddo
c
      nsend = NC*NC*nds
      call mpi_isend(bufds(1,1,1), nsend, mpi_complex16, 
     &               jdown, 1, mpi_comm_world, isend2, ierr)

c     .....  Now we receive data   ...
      nrecv = NC*NC*(NX+2)*NZ*NT/2
      call mpi_irecv(bufdr(1,1,1), nrecv, mpi_complex16, 
     &               jdown, 1, mpi_comm_world, irecv1, ierr)

      nrecv = NC*NC*(NX+2)*NZ*NT/2
      call mpi_irecv(bufur(1,1,1), nrecv, mpi_complex16, 
     &               jup, 1, mpi_comm_world, irecv2, ierr)

      call mpi_wait(isend1, istatus, ierr)
      call mpi_wait(isend2, istatus, ierr)
      call mpi_wait(irecv1, istatus, ierr)
      call mpi_wait(irecv2, istatus, ierr)

c     ... from down ... 
      ndr = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+0+iz+it+iflg,2)
        do ix = 1-ieo_sta, NX+1, 2
         ndr = ndr + 1
          do k1 = 1, NC
          do k2 = 1, NC
            u%g(k1, k2, ix, 0, iz, it) = bufdr(k1, k2,  ndr)
          enddo
          enddo
        enddo
      enddo
      enddo

c     ... from up ...
      nur = 0
      do it = 1, NT
      do iz = 1, NZ
        ieo_sta = mod(1+(NY+1)+iz+it+iflg,2)
        do ix =  1-ieo_sta, NX+1, 2
        nur = nur + 1
          do k1 = 1, NC
          do k2 = 1, NC
            u%g(k1, k2, ix, NY+1, iz, it) = bufur(k1, k2, nur)
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

