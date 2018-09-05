c----------------------------------------------------------------------c
      subroutine init(u,istart,ndelay)
c----------------------------------------------------------------------c
*     Initialize the gauge configuration                               * 
*     istart = 1  Cold start                                           *
*              2  Hot  start                                           *
*              3  File start                                           *
*                                           written by AN Feb.15,2000  *
*                                           checked by SC Feb.17,2000  *
c----------------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'
      include '../INCLUDE/para.h'                                     ! MPI

      TYPE(g_field0) u(4)
      TYPE(g_field1) temp 
      COMPLEX*16  temp2(3,3,NV)                                       ! MPI

      REAL*8  beta

c     .....  prepare for the random numbers   .....
      if((istart==1).or.(istart==2))  then
         ndelay = 1 + myrank
         call cinit3(ndelay) 
      endif

      if(istart==1)  then
   
        IF (myrank == 0) THEN                                         ! MPI

         write(*,*) ".....  Cold start"

         do ibush = 0, NBUSH-1 
           do mu = 1, 4
             call set_ident(temp)
             u(mu) = temp
             u(mu)%direction = mu
             u(mu)%parity = ibush
           enddo 
         enddo

        ENDIF                                                         ! MPI

      endif

      if(istart==2)  then
   
        IF (myrank == 0) THEN                                         ! MPI

         write(*,*) ".....  Hot start"

         do ibush = 0, NBUSH-1 
           do mu = 1, 4
             call set_rand(temp)
             u(mu) = temp
             u(mu)%direction = mu
             u(mu)%parity = ibush
           enddo 
         enddo

        ENDIF                                                         ! MPI

      endif

      if(istart==3)  then
   
        IF (myrank == 0) THEN                                         ! MPI

         write(*,*) ".....  File start"

         open(unit=1,file='fort.1',form='unformatted',
     &        action='read',status='old')
         read(1) beta, ndelay
         do mu = 1, 4
           read(1)  (((((( u(mu)%g(ic,jc,ix,iy,iz,it),
     &               ic=1,NC),jc=1,NC),
     &               ix=1,NX),iy=1,NY),iz=1,NZ),it=1,NT)
         enddo
         close(unit=1)

        ENDIF                                                         ! MPI

        CALL MPI_BCAST(ndelay,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)    ! MPI
        call cinit3(ndelay) 

      endif


c     *****   Broadcast Link variables U                            ! MPI 
      do ibush = 0, NBUSH-1                                         ! MPI
      do mu = 1, 4                                                  ! MPI

      IF (myrank == 0) THEN                                         ! MPI
          temp = u(mu)                                              ! MPI
          DO ic1 = 1, 3                                             ! MPI
          DO ic2 = 1, 3                                             ! MPI
          DO is = 1, NV                                             ! MPI
          temp2(ic1,ic2,is) = temp%g(ic1,ic2,is)                    ! MPI
          ENDDO                                                     ! MPI
          ENDDO                                                     ! MPI
          ENDDO                                                     ! MPI
      ENDIF                                                         ! MPI

      CALL MPI_BCAST                                                ! MPI
     &    (temp2,3*3*NV,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)        ! MPI

      IF (myrank .ne. 0) THEN                                       ! MPI
          DO ic1 = 1, 3                                             ! MPI
          DO ic2 = 1, 3                                             ! MPI
          DO is = 1, NV                                             ! MPI
          temp%g(ic1,ic2,is) = temp2(ic1,ic2,is)                    ! MPI
          ENDDO                                                     ! MPI
          ENDDO                                                     ! MPI
          ENDDO                                                     ! MPI
          u(mu) = temp                                              ! MPI
      ENDIF                                                         ! MPI

  
      u(mu)%direction = mu                                          ! MPI
      u(mu)%parity = ibush                                          ! MPI

      enddo                                                         ! MPI
      enddo                                                         ! MPI

c     *****   Fill wings (fringes)   ****
      do mu = 1, 4
        call set_wing_g2(u(mu)) 
      enddo

      return
      end
