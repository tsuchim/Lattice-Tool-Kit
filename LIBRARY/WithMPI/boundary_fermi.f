c boundary_fermi.f
c fermion part of "boundary.f"
c----------------------------------------------------------------------c
      subroutine set_wing_f1(a)
c----------------------------------------------------------------------c
c     Fill wings of fermion vectors according to the boudary condition.c
c----------------------------------------------------------------------c
      USE field_f
      USE fpara

      include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'
*      parameter( NC=3 )

      TYPE(f_field) a

      COMPLEX*16 works1(NC,NY*NZ*NT,4), works2(NC,NY*NZ*NT,4), 
     &           works3(NC,NX*NZ*NT,4), works4(NC,NX*NZ*NT,4), 
     &           workr1(NC,NY*NZ*NT,4), workr2(NC,NY*NZ*NT,4),
     &           workr3(NC,NX*NZ*NT,4), workr4(NC,NX*NZ*NT,4)
      integer istatus(mpi_status_size)
      REAL*8 fac

c     ...  Neibouring PE's   ...
      iup   = itable(myranki+1, myrankj)
      idown = itable(myranki-1, myrankj)
      jup   = itable(myranki  , myrankj+1)
      jdown = itable(myranki  , myrankj-1)
c
c     ....  Now we send data
      !  X-direction
      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
        icum = icum + 1
           do k = 1, NC
             works1(k, icum, ialpha) = a%f(k, NX, iy, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo

      nsend = NC*NY*NZ*NT*4
      call mpi_isend(works1(1,1,1), nsend, mpi_complex16,
     &               iup, 1, mpi_comm_world, jsend1, ierr)
c
      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
        icum = icum + 1
           do k = 1, NC
             works2(k, icum, ialpha) = a%f(k, 1, iy, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*NY*NZ*NT*4
      call mpi_isend(works2(1,1,1), nsend, mpi_complex16,
     &               idown, 1, mpi_comm_world, jsend2, ierr)
c
c     .....  Now we receive data   ...
      nrecv = NC*NY*NZ*NT*4
      call mpi_irecv(workr1(1,1,1), nrecv, mpi_complex16,
     &               idown, 1, mpi_comm_world, jrecv1, ierr)

      nrecv = NC*NY*NZ*NT*4
      call mpi_irecv(workr2(1,1,1), nrecv, mpi_complex16,
     &               iup, 1, mpi_comm_world, jrecv2, ierr)

      call mpi_wait(jsend1, istatus, ierr)
      call mpi_wait(jsend2, istatus, ierr)
      call mpi_wait(jrecv1, istatus, ierr)
      call mpi_wait(jrecv2, istatus, ierr)

      if (myranki == IPROCS-1) then
         fac = BC(1)
      else
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
        icum = icum + 1
          do k = 1, NC
            a%f(k,NX+1,iy,iz,it,ialpha)=fac*workr2(k,icum,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

      if (myranki == 0) then
         fac = BC(1)
      else 
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1, NY
        icum = icum + 1
          do k = 1, NC
            a%f(k,0,iy,iz,it,ialpha)=fac*workr1(k,icum,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

      !  Y-direction
      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1, NX 
        icum = icum + 1
           do k = 1, NC
             works3(k, icum, ialpha) = a%f(k, ix, NY, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*NX*NZ*NT*4
      call mpi_isend(works3(1,1,1), nsend, mpi_complex16,
     &               jup, 1, mpi_comm_world, isend1, ierr)
c
      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1, NX
        icum = icum + 1
           do k = 1, NC
             works4(k, icum, ialpha) = a%f(k, ix, 1, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*NX*NZ*NT*4
      call mpi_isend(works4(1,1,1), nsend, mpi_complex16,
     &               jdown, 1, mpi_comm_world, isend2, ierr)


c     .....  Now we receive data   ...
      nrecv = NC*NX*NZ*NT*4
      call mpi_irecv(workr3(1,1,1), nrecv, mpi_complex16,
     &               jdown, 1, mpi_comm_world, irecv1, ierr)

      nrecv = NC*NX*NZ*NT*4
      call mpi_irecv(workr4(1,1,1), nrecv, mpi_complex16,
     &               jup, 1, mpi_comm_world, irecv2, ierr)

      call mpi_wait(isend1, istatus, ierr)
      call mpi_wait(isend2, istatus, ierr)
      call mpi_wait(irecv1, istatus, ierr)
      call mpi_wait(irecv2, istatus, ierr)
c
      if(myrankj==JPROCS-1)  then
         fac = BC(2)
      else
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1, NX
        icum = icum + 1
          do k = 1, NC
            a%f(k,ix,NY+1,iz,it,ialpha)=fac*workr4(k,icum,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

      if (myrankj == 0) then
         fac = BC(2)
      else 
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      icum = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1, NX
        icum = icum + 1
          do k = 1, NC
            a%f(k,ix,0,iz,it,ialpha)=fac*workr3(k,icum,ialpha)
          enddo
        enddo
      enddo; enddo
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

      include 'mpif.h'
      include '../INCLUDE/para.h'
      include '../INCLUDE/para_geometry'
*      parameter( NC=3 )

      TYPE(f_field) a

      COMPLEX*16 bufrs(NC,NY*NZ*NT,4), bufls(NC,NY*NZ*NT,4), 
     &           bufrr(NC,NX*NZ*NT,4), buflr(NC,NX*NZ*NT,4), 
     &           bufus(NC,NY*NZ*NT,4), bufds(NC,NY*NZ*NT,4),
     &           bufur(NC,NX*NZ*NT,4), bufdr(NC,NX*NZ*NT,4)
      integer istatus(mpi_status_size)
      REAL*8 fac

c     ...  Neibouring PE's   ...
      iup   = itable(myranki+1, myrankj)
      idown = itable(myranki-1, myrankj)
      jup   = itable(myranki  , myrankj+1)
      jdown = itable(myranki  , myrankj-1)
c
c     ...  Starting Point   ...
      kld_sta = mod(1+1+iflg,2)
      krd_sta = mod(NX+1+iflg,2)
      klu_sta = mod(1+NY+iflg,2) 

c     ....  Now we send data
      !  X-direction
c     .... iup ... 
      do ialpha = 1, 4
      nrs = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1 + krd_sta, NY, 2
        nrs = nrs + 1
           do k = 1, NC
             bufrs(k, nrs, ialpha) = a%f(k, NX, iy, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo

      nsend = NC*nrs*4
      call mpi_isend(bufrs(1,1,1), nsend, mpi_complex16,
     &               iup, 1, mpi_comm_world, jsend1, ierr)
c
c     ... idown ... 
      do ialpha = 1, 4
      nls = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1+kld_sta,NY,2
        nls = nls + 1
           do k = 1, NC
             bufls(k, nls, ialpha) = a%f(k, 1, iy, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*nls*4
      call mpi_isend(bufls(1,1,1), nsend, mpi_complex16,
     &               idown, 1, mpi_comm_world, jsend2, ierr)
c
c     .....  Now we receive data   ...
      nrecv = NC*NY*NZ*NT*4/2
      call mpi_irecv(buflr(1,1,1), nrecv, mpi_complex16,
     &               idown, 1, mpi_comm_world, jrecv1, ierr)

      nrecv = NC*NY*NZ*NT*4/2
      call mpi_irecv(bufrr(1,1,1), nrecv, mpi_complex16,
     &               iup, 1, mpi_comm_world, jrecv2, ierr)

      call mpi_wait(jsend1, istatus, ierr)
      call mpi_wait(jsend2, istatus, ierr)
      call mpi_wait(jrecv1, istatus, ierr)
      call mpi_wait(jrecv2, istatus, ierr)

c     ... starting point ...
      kld_sta2 = 1 - kld_sta
      krd_sta2 = 1 - krd_sta
      klu_sta2 = 1 - klu_sta

c     ... from left ..
      if (myranki == 0) then
         fac = BC(1)
      else
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      nlr = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy = 1+kld_sta2, NY, 2
        nlr = nlr + 1
          do k = 1, NC
            a%f(k,0,iy,iz,it,ialpha)=fac*buflr(k,nlr,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

c     ... from right ...       
      if (myranki == IPROCS-1) then
         fac = BC(1)
      else 
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      nrr = 0   
      do it = 1, NT; do iz = 1, NZ
        do iy =  1+krd_sta2, NY, 2
        nrr = nrr + 1
          do k = 1, NC
            a%f(k,NX+1,iy,iz,it,ialpha)=fac*bufrr(k,nrr,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

      !  Y-direction
c     ... jup ...
      do ialpha = 1, 4
      nus = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix =  1+klu_sta, NX, 2
        nus = nus + 1
           do k = 1, NC
             bufus(k, nus, ialpha) = a%f(k, ix, NY, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*nus*4
      call mpi_isend(bufus(1,1,1), nsend, mpi_complex16,
     &               jup, 1, mpi_comm_world, isend1, ierr)
c
c     ... jdown ...
      do ialpha = 1, 4
      nds = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1+kld_sta, NX, 2
        nds = nds + 1
           do k = 1, NC
             bufds(k, nds, ialpha) = a%f(k, ix, 1, iz, it, ialpha)
           enddo
        enddo
      enddo; enddo
      enddo
c
      nsend = NC*nds*4
      call mpi_isend(bufds(1,1,1), nsend, mpi_complex16,
     &               jdown, 1, mpi_comm_world, isend2, ierr)


c     .....  Now we receive data   ...
      nrecv = NC*NX*NZ*NT*4/2
      call mpi_irecv(bufdr(1,1,1), nrecv, mpi_complex16,
     &               jdown, 1, mpi_comm_world, irecv1, ierr)

      nrecv = NC*NX*NZ*NT*4/2
      call mpi_irecv(bufur(1,1,1), nrecv, mpi_complex16,
     &               jup, 1, mpi_comm_world, irecv2, ierr)

      call mpi_wait(isend1, istatus, ierr)
      call mpi_wait(isend2, istatus, ierr)
      call mpi_wait(irecv1, istatus, ierr)
      call mpi_wait(irecv2, istatus, ierr)
c
c     ... from down ...
      if(myrankj==0)  then
         fac = BC(2)
      else
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      ndr = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix = 1+kld_sta2, NX, 2
        ndr = ndr + 1
          do k = 1, NC
            a%f(k,ix,0,iz,it,ialpha)=fac*bufdr(k,ndr,ialpha)
          enddo
        enddo
      enddo; enddo
      enddo

c     ... from up ... 
      if (myrankj == JPROCS-1) then
         fac = BC(2)
      else 
         fac = 1.0d0
      end if

      do ialpha = 1, 4
      nur = 0   
      do it = 1, NT; do iz = 1, NZ
        do ix =  1+klu_sta2, NX, 2
        nur = nur + 1
          do k = 1, NC
            a%f(k,ix,NY+1,iz,it,ialpha)=fac*bufur(k,nur,ialpha)
          enddo
        enddo
      enddo; enddo
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
















