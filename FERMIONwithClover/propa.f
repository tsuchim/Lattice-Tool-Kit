c------------------------------------------------------c
      SUBROUTINE propa 
c------------------------------------------------------c
c     Calculate meson and quark propagators
c------------------------------------------------------c
      include '../INCLUDE/para_geometry'

*      REAL*8,INTENT(IN) :: eps   ! Stopping accuracy in CG
*      INTEGER,INTENT(IN):: imax  ! Iteration number  in CG

      COMPLEX*16 g_pi(NT)  ! g_pi  : pion propagator

c     ...  pion propagator
      call q_propa(g_pi)

      write(*,*) 'Pion Propagator'
      DO it = 1, NT
        write(*,*) 'it =', it, ' g_pi = ',g_pi(it)
      ENDDO

      return
      END

c--------------------------------------------------------c
      SUBROUTINE q_propa(g_pi)
c--------------------------------------------------------c
c     Quark and Pion propagator with simple point source c
c                                                        c
c     Quark propagator in Momentum Space with simple     c
c     point source                                       c
c                                                        c
c     INPUT eps,imax : Stopping condition of CG          c
c     OUTPUT g_pi : connected propagators                c
c--------------------------------------------------------c
      USE field_f
      include '../INCLUDE/para_geometry'

      COMPLEX*16,INTENT(OUT):: g_pi(NT)

      TYPE(f_field) x,b

      INTEGER mu,nu                   ! Dirac indices
      INTEGER ix_i, iy_i, iz_i, it_i  ! Source Points
      COMPLEX*16 CI, sum, ftfac 
      PARAMETER( CI=(0.d0,1.d0) )
      REAL*8 fac, PI

      PI = 4.d0*atan(1.d0)

c     ... Soruce Point
      ix_i = 1 
      iy_i = 1
      iz_i = 1
      it_i = 1 

c     ...  Initialize
      DO it = 1,NT
        g_pi(it) = cmplx(0.0d0,0.0d0)
      ENDDO 
      
      SourceColor : DO ic = 1 , 3
      SourceDirac : DO mu = 1 , 4  ! Dirac index for Source

        b = zero_f    ! Initialize

c       ...  Set the right hand of Wx=b
        b%f(ic,ix_i,iy_i,iz_i,it_i,mu)=cmplx(1.d0,0.d0)
 
        iflag = 1  ! (iflag=1 Solve Wx=b, iflag=2 Solve W_adj*x=b)
c        call cg0(x,b,eps,imax,iflag)
        call cg0(x,b,iflag)
 
        DO it = 1 , NT
        DO ix = 1 , NX
        DO iy = 1 , NY
        DO iz = 1 , NZ

          ! ... pion propagator
          DO ic1 = 1 , 3
          DO mu1 = 1 , 4

            g_pi(it) = g_pi(it) +
     &                 x%f(ic1,ix,iy,iz,it,mu1)
     &                 * conjg(x%f(ic1,ix,iy,iz,it,mu1))

          ENDDO
          ENDDO

        ENDDO
        ENDDO
        ENDDO
        ENDDO


      ENDDO SourceDirac
      ENDDO SourceColor


      DO it = 1,NT
        g_pi(it) = g_pi(it)/(NX*NY*NZ)
      ENDDO
c
      return
      END
