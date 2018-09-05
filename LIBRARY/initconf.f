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

      TYPE(g_field0) u(4)
      TYPE(g_field1) temp 

      REAL*8  beta

c     .....  prepare for the random numbers   .....
      if((istart==1).or.(istart==2))  then
         ndelay = 1 + myrank
         call cinit3(ndelay) 
      endif

      if(istart==1)  then
   

         write(*,*) ".....  Cold start"

         do ibush = 0, NBUSH-1 
           do mu = 1, 4
             call set_ident(temp)
             u(mu) = temp
             u(mu)%direction = mu
             u(mu)%parity = ibush
           enddo 
         enddo


      endif

      if(istart==2)  then
   

         write(*,*) ".....  Hot start"

         do ibush = 0, NBUSH-1 
           do mu = 1, 4
             call set_rand(temp)
             u(mu) = temp
             u(mu)%direction = mu
             u(mu)%parity = ibush
           enddo 
         enddo


      endif

      if(istart==3)  then
   

         write(*,*) ".....  File start"

         open(unit=1,file='fort.1',form='unformatted',
     &        action='read',status='old')
         read(1) beta, ndelay
         write(*,*) "beta, ndelay: ", beta, ndelay
         do mu = 1, 4
           read(1)  (((((( u(mu)%g(ic,jc,ix,iy,iz,it),
     &               ic=1,NC),jc=1,NC),
     &               ix=1,NX),iy=1,NY),iz=1,NZ),it=1,NT)
           write(*,*) "... mu = ", mu
           write(*,*) "u: ", u(mu)%g(1,1,1,1,1,1), 
     &                       u(mu)%g(3,3,Nx,Ny,Nz,Nt)
         enddo
         close(unit=1)


        call cinit3(ndelay) 

      endif

c     *****   Fill wings (fringes)   ****
      do mu = 1, 4
        call set_wing_g2(u(mu)) 
      enddo

      return
      end
