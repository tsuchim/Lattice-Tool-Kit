c----------------------------------------------------------------------c
      program qcd0
c----------------------------------------------------------------------c
*     Gauge Update Program.                                            *
c----------------------------------------------------------------------c
      USE field_g
      USE aniso
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u(4)
      common/ config/ u

      REAL*8 beta, plaq
      TYPE(g_field0) g
      TYPE(g_field1) temp
      COMPLEX*16 :: Pol(3,3,NX,NY,NZ), avePol


c     *****   Read input parameters   *****
      call  rdparam(beta,gamma_G,istart,nsweep0,nsweep1)
      write(*,*) 'beta=', beta


c     *****   starting condition   *****
      call  init(u,istart,ndelay)

      call  meas1(u,plaq)
      write(*,*) "Plaq at the begining : ", plaq

      call  CalPol(u,Pol,avePol)
      write(*,*) "Pol at the begining : ", avePol

c     *****   test (Random Gauge Transformation)  ****
      do ibush = 0, NBUSH-1
         call set_rand(temp)
         g = temp 
      enddo
      call grotate (u,g)
      do mu = 1, 4
        call set_wing_g2(u(mu))  ! Added on July 13, 2000 (AN)
      enddo
      call  meas1(u,plaq)
      write(*,*) "Plaq after random gauge trans. : ", plaq

c     *****   Monte Carlo sweeps   *****
      SWEEP0 : do n0 = 1, nsweep0

         SWEEP1 : do n1 = 1, nsweep1

*           call  update1(u,beta)   !  for Metropolis 
            call  update3(u,beta)   !  for the pseudo-heat bath

         enddo SWEEP1

         call  meas1(u,plaq)
         write(*,*) "n0, plaq : ", n0, plaq
         call  CalPol(u,Pol,avePol)
         write(*,*) "n0, Pol : ", n0, avePol

         ! Normalize link variables
         do ibush = 0, NBUSH-1
           do mu = 1, 4
             temp = u(mu)
             call normalize(temp)
             u(mu) = temp
           enddo
         enddo
         call  meas1(u,plaq)
         write(*,*) "Check: PLAQ after normalization= ", plaq

*        ! ... Save configurations
*        ilog = 10 + n0 
*        write(ilog) beta, ndelay
*        do mu = 1, 4
*          write(ilog) (((((( u(mu)%g(ic,jc,ix,iy,iz,it),
*    &               ic=1,NC),jc=1,NC),
*    &               ix=1,NX),iy=1,NY),iz=1,NZ),it=1,NT)
*        enddo
         
      enddo SWEEP0

c     *****   save the configuraton *****
      call  save_conf(u,beta,ndelay)


      end

c----------------------------------------------------------------------c
      subroutine rdparam(beta,gamma_G,istart,nsweep0,nsweep1)
c----------------------------------------------------------------------c
*     Subrouitine to read parameters
c----------------------------------------------------------------------c
      include '../INCLUDE/para_geometry'
      REAL*8 beta, gamma_G

        read(*,*) beta  
        read(*,*) istart
        read(*,*) nsweep0, nsweep1
        read(*,*) gamma_G



        write(*,*) "beta : ", beta
        write(*,*) "Lattice size : ", NX,NY,NZ,NT
        write(*,*) "istart : ", istart
        write(*,*) "nsweep0,nsweep1 : ", nsweep0,nsweep1
        write(*,*) "gamma_G : ", gamma_G

        write(*,*) "NBUSH : ", NBUSH
        if(NBUSH==1)  then
          write(*,*) "(Hybrid MC Program)"
        else if(NBUSH==2)  then
          write(*,*) "(Even/Odd for Gauge)"
        else if(NBUSH==32)  then
          write(*,*) "(Hypercube for Gauge)"
        else
          write(*,*) "Nbush is correct ???"
        endif

   

      RETURN
      END
