c rand_sr.f
*     REAL*8 vran(10000)
*
*     ndelay = 1
*     call cinit3(ndelay)
*
*     call rands(vran,10000)
*
*     x = ranf(idum)
*
c----------------------------------------------------------------------c
      SUBROUTINE cinit3(ndelay)
c----------------------------------------------------------------------c
      write(*,*) "ndelay : ", ndelay
      CALL init3
      CALL delay3(ndelay)
      END

c----------------------------------------------------------------------c
      SUBROUTINE init3
c----------------------------------------------------------------------c
      parameter( IP=521, IQ=32 )
      parameter( MACRM= 40, MACRI=1)
      COMMON /crand/ iw(0:IP-1), jr, kr
      INTEGER ib(0:IP-1)

      ix=MACRI

      DO i = 0, IP-1
        ix=ix*69069
        ib(i)=ishft(ix,-31)
      ENDDO

      jr=0
      kr=IP-IQ
      do 30 j=0,IP-1
        iwork=0
        do 20 i=0,31
          iwork=iwork*2+ib(jr)
          ib(jr)= IEOR(ib(jr),ib(kr))
          jr=jr+1
          if (jr.eq.IP) jr=0
          kr=kr+1
          if (kr.eq.IP) kr=0
   20   continue
        iw(j)=ishft(iwork,-1)
   30 continue

      END

c----------------------------------------c
      SUBROUTINE delay3(lambda)
c----------------------------------------c
      PARAMETER( IP=521, IQ=32 )
      PARAMETER( MACRM=40, MACRI=1)

      INTEGER iwk(0:2*IP-2),c(0:2*IP-1),ib(0:IP+31)
      COMMON /crand/ iw(0:IP-1), jr, kr  

      mu = MACRM
      DO i = 0,IP-1
        iwk(i) = iw(i)
      ENDDO

      DO i = IP, 2*IP-2
        iwk(i) =  IEOR(iwk(i-IP),iwk(i-IQ))
      ENDDO

      DO i = 0, mu-1
        ib(i) = 0
      ENDDO

      m = lambda
      nb = mu-1

      DO WHILE(m>IP-1)  
        nb = nb+1
        ib(nb) = mod(m,2)
        m = m/2
      ENDDO

      DO i = 0, IP-1
        c(i) = 0
      ENDDO
      c(m) = 1

      DO j = nb, 0, -1

        DO i = IP-1, 0, -1
          c(2*i+ib(j)) = c(i)
          c(2*i+1-ib(j)) = 0
        ENDDO

        DO i = 2*IP-1, IP, -1
          c(i-IP) =  IEOR(c(i-IP),c(i))
          c(i-IQ) =  IEOR(c(i-IQ),c(i))
        ENDDO

      ENDDO


      DO j = 0,IP-1

        iwork = 0
        DO i = 0,IP-1
          iwork =  IEOR(iwork,c(i)*iwk(j+i))
        ENDDO

        iw(j) = iwork

      ENDDO

      RETURN
      END

c----------------------------------------------------------------------c
      FUNCTION ranf(idumy)
c----------------------------------------------------------------------c
c     random number a la lagged fibonacci
c     idumy plays no role, i.e., it is dummy
c----------------------------------------------------------------------c
      REAL*4 ranf
      PARAMETER( IP=521, IQ=32 )
      PARAMETER( MACRM= 40, MACRI=1)
      PARAMETER(FNORM=0.465661e-9)
      COMMON /crand/ ir(0:IP-1), j, k

      ir(j)= IEOR(ir(j),ir(k))
      irnd=ir(j)
      frnd=irnd*FNORM
      ranf=frnd

      j=j+1
      if(j.ge.IP) j=j-IP
      k=k+1
      if(k.ge.IP) k=k-IP

      return
      end

c----------------------------------------------------c
      SUBROUTINE rands(rn,n)
c----------------------------------------------------c
c     Same for ranf, but produce array
c----------------------------------------------------c
c     INPUT  n     : No. of random numbers required
c     OUTPUT rn(n) : uniform random numbers 
c                    (Double precision)
c----------------------------------------------------c
      REAL*8,  INTENT(OUT):: rn(n)
      INTEGER, INTENT(IN)::  n

      PARAMETER( IP=521, IQ=32 )
      REAL*8 FNORM
      PARAMETER(FNORM=0.465661d-9)

      COMMON /crand/ iw(0:IP-1), jr, kr

      nhasu = mod(n,IQ)
      nrepeat = (n-nhasu) / IQ

C     ...  Handle hasu  ...

      itot = 1

      IF(nhasu.NE.0) THEN  

        DO i = 1, nhasu

          irnd = ieor(iw(jr),iw(kr))  
          rn(itot)  = FNORM * DBLE(irnd)
          iw(jr)  = irnd 

          jr = jr + 1
          kr = kr + 1
          IF(jr .GE. IP) jr = jr-IP
          IF(kr .GE. IP) kr = kr-IP
          itot = itot + 1

        ENDDO

      ENDIF

      IF(nrepeat==0)  RETURN
      Repeat: DO nn = 1, nrepeat

        DO i = 1, IQ

          irnd = ieor(iw(jr),iw(kr))  
          rn(itot) = FNORM * DBLE(irnd)
          iw(jr) = irnd 

          jr = jr + 1
          kr = kr + 1
          IF(jr .GE. IP) jr = jr-IP
          IF(kr .GE. IP) kr = kr-IP
          itot = itot + 1

        ENDDO

      ENDDO Repeat

      RETURN
      END

c--------------------------------------------------------------c
      SUBROUTINE rgauss1(xr,nup)
c--------------------------------------------------------------c
c     gaussian random number                                   c
c     normalized to unity                                      c
c     notice:  nup should be even.                             c
c--------------------------------------------------------------c
      INTEGER, INTENT(IN) :: nup
      REAL*8,  INTENT(OUT):: xr(nup)

      REAL*8 PI, PI2
      PARAMETER(PI=3.1415926535897932,PI2=2.d0*PI)

      nup2 = nup/2
      IF (nup.ne.(2*nup2)) THEN
         write(*,'(1x,a,a,i5)')
     &   ' !!!!  Attention nup in rgauss should be even',
     *   '  nup = ', nup
         stop
      ENDIF
      
      DO i = 1, nup2

        i2 = 2*i
        i1 = i2 - 1

        v1 = sqrt(-log(ranf(idum)+1.e-10)*2.)
*       v2 = 6.28318*ranf(idum)
        v2 = PI2*ranf(idum)      ! Corrected Feb.10,2005
 
        xr(i1) = v1 * cos(v2)
        xr(i2) = v1 * sin(v2)

      ENDDO

      return
      end
