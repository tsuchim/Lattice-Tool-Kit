c--------------------------------------------------------c
      SUBROUTINE save_conf(u,beta,ndelay)
c--------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0),INTENT(IN):: u(4)
      REAL*8,INTENT(IN) :: beta
      INTEGER,INTENT(IN):: ndelay

      OPEN(unit=2,file='fort.2',form='unformatted',
     &        action='write',status='unknown')

          nincr = 1
      
      write(2) beta, ndelay+nincr
      do mu = 1, 4
        write(2) (((((( u(mu)%g(ic,jc,ix,iy,iz,it),
     &            ic=1,NC),jc=1,NC),
     &            ix=1,NX),iy=1,NY),iz=1,NZ),it=1,NT)
      enddo

      CLOSE(2)

      END
