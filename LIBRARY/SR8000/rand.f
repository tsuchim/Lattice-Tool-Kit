c rand_sr.f
c Random numbers for SR8000
*
*     REAL*8 vran(20480)
*
*     ...  Initialization
*     ndelay = 1
*     call cinit3p(ndelay)
*     call cinit3 (ndelay)
*
*     ....
*
*     call rands(vran,2048)
*
c----------------------------------------------------------------------c
      SUBROUTINE cinit3p(ndelay)
c----------------------------------------------------------------------c
      PARAMETER( IP=521, IQ=32, NPE=512 )
      COMMON /randp/ iw(NPE,0:IP-1), jr, kr
      INTEGER iwtmp(0:IP-1), lambda(NPE)

      CALL init3s( iwtmp, jr, kr ) 

      DO n = 0, IP-1
        DO ipe = 1, NPE
          iw(ipe,n) = iwtmp(n)
        ENDDO
      ENDDO

      DO ipe = 1, NPE
        lambda(ipe) = ndelay
        ndelay = ndelay + 1
      ENDDO

      call delay3p(lambda,iw) 

      RETURN
      END

c----------------------------------------------------------------------c
      SUBROUTINE init3s(iw,jr,kr)
c----------------------------------------------------------------------c
      PARAMETER( IP=521, IQ=32 )
      PARAMETER( MACRM= 40, MACRI=1)

      INTEGER iw(0:IP-1), jr, kr
      INTEGER ib(0:IP-1)

      ix = MACRI
      DO i = 0,IP-1
        ix = ix*69069
        ib(i) = ishft(ix,-31)
      ENDDO

      jr = 0
      kr = IP-IQ

      DO j = 0, IP-1

        iwork = 0

        DO i = 0, 31

          iwork = iwork*2 + ib(jr)
          ib(jr) =  ieor(ib(jr),ib(kr))
          jr = jr + 1
          if (jr.eq.IP) jr = 0
          kr = kr + 1
          if (kr.eq.IP) kr = 0

        ENDDO

        iw(j) = ishft(iwork,-1)

      ENDDO

      RETURN
      END

c----------------------------------------c
      SUBROUTINE delay3s(lambda,iw)
c----------------------------------------c
      PARAMETER( IP=521, IQ=32 )
      PARAMETER( MACRM=40, MACRI=1)

      INTEGER iw(0:IP-1)
      INTEGER iwk(0:2*IP-2),c(0:2*IP-1),ib(0:IP+31)

      mu = MACRM
      DO i = 0,IP-1
        iwk(i) = iw(i)
      ENDDO

      DO i = IP, 2*IP-2
        iwk(i) =  ieor(iwk(i-IP),iwk(i-IQ))
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
          c(i-IP) =  ieor(c(i-IP),c(i))
          c(i-IQ) =  ieor(c(i-IQ),c(i))
        ENDDO

      ENDDO


      DO j = 0,IP-1

        iwork = 0
        DO i = 0,IP-1
          iwork =  ieor(iwork,c(i)*iwk(j+i))
        ENDDO

        iw(j) = iwork

      ENDDO

      RETURN
      END

c----------------------------------------c
      SUBROUTINE delay3p(lambda,iw)
c----------------------------------------c
      PARAMETER( IP=521, IQ=32, NPE=512 )
      PARAMETER( MACRM=40, MACRI=1)

      INTEGER,INTENT(IN):: lambda(NPE)
      INTEGER,INTENT(INOUT):: iw(NPE,0:IP-1)
      INTEGER iwk(0:2*IP-2),c(0:2*IP-1,NPE),ib(0:IP+31,NPE)

      mu = MACRM

*POPTION PARALLEL
      DO ipe = 1, NPE

*VOPTION INDEP
        DO i = 0,IP-1
          iwk(i) = iw(ipe,i)
        ENDDO

*VOPTION INDEP
        DO i = IP, 2*IP-2
          iwk(i) =  ieor(iwk(i-IP),iwk(i-IQ))
        ENDDO

*VOPTION INDEP
        DO i = 0, mu-1
          ib(i,ipe) = 0
        ENDDO

        m = lambda(ipe)
        nb = mu-1
  
        DO WHILE(m>IP-1)  
          nb = nb+1
          ib(nb,ipe) = mod(m,2)
          m = m/2
        ENDDO

*VOPTION INDEP
        DO i = 0, IP-1
          c(i,ipe) = 0
        ENDDO
        c(m,ipe) = 1
  
        DO j = nb, 0, -1
  
          DO i = IP-1, 0, -1
            c(2*i+ib(j,ipe),ipe) = c(i,ipe)
            c(2*i+1-ib(j,ipe),ipe) = 0
          ENDDO
  
          DO i = 2*IP-1, IP, -1
            c(i-IP,ipe) =  ieor(c(i-IP,ipe),c(i,ipe))
            c(i-IQ,ipe) =  ieor(c(i-IQ,ipe),c(i,ipe))
          ENDDO

        ENDDO


        DO j = 0,IP-1

          iwork = 0
          DO i = 0,IP-1
            iwork =  ieor(iwork,c(i,ipe)*iwk(j+i))
          ENDDO

          iw(ipe,j) = iwork

        ENDDO

      ENDDO

      RETURN
      END

c----------------------------------------------------c
      SUBROUTINE rands(rn,n)
c----------------------------------------------------c
c     Same for ranf, but produce array
c----------------------------------------------------c
c     INPUT  n     : No. of random numbers required
c     OUTPUT rn(n) : uniform random numbers 
c                    (Double precision)
c----------------------------------------------------c
      PARAMETER( IP=521, IQ=32, NPE=512 )
      REAL*8 FNORM
      PARAMETER(FNORM=0.465661d-9)

      REAL*8 rn(NPE,n/NPE)
      COMMON /randp/ iw(NPE,0:IP-1), jr, kr

*     IF(MOD(n,NPE).NE.0)  THEN
*        WRITE(*,*) "n is not multiple of NPE"
*        WRITE(*,*) "n, NPE : ", n, NPE
*        STOP
*     ENDIF

      nn = n/NPE
      nhasu = mod(nn,IQ)
      nrepeat = (nn-nhasu) / IQ

C     ...  Handle hasu  ...

      itot = 1

      IF(nhasu.NE.0) THEN  

*POPTION PARALLEL
        DO i = 1, nhasu

*VOPTION INDEP
          DO ipe = 1, NPE

             irnd = ieor(iw(ipe,jr),iw(ipe,kr))  
             rn(ipe,itot)  = FNORM * DBLE(irnd)
             iw(ipe,jr)  = irnd 

          ENDDO

          jr = jr + 1
          kr = kr + 1
          IF(jr .GE. IP) jr = jr-IP
          IF(kr .GE. IP) kr = kr-IP
          itot = itot + 1

        ENDDO

      ENDIF

      IF(nrepeat==0)  RETURN
      Repeat: DO nn = 1, nrepeat

*POPTION PARALLEL
        DO i = 1, IQ

*VOPTION INDEP
          DO ipe = 1, NPE

            irnd = ieor(iw(ipe,jr),iw(ipe,kr))  
            rn(ipe,itot) = FNORM * DBLE(irnd)
            iw(ipe,jr) = irnd 

          ENDDO

          jr = jr + 1
          kr = kr + 1
          IF(jr .GE. IP) jr = jr-IP
          IF(kr .GE. IP) kr = kr-IP
          itot = itot + 1

        ENDDO

      ENDDO Repeat

      RETURN
      END

c----------------------------------------------------------------------c
      FUNCTION ranf(idumy)
c----------------------------------------------------------------------c
c     random number a la lagged fibonacci
c     idumy plays no role, i.e., it is dummy
c----------------------------------------------------------------------c
      REAL*4 ranf
      INTEGER p,q
      PARAMETER( IP= 521, IQ= 32 )
      PARAMETER( macrm= 40, macri=1)
      PARAMETER(p=IP,q=IQ,fnorm=0.465661e-9)
      INTEGER ir(0:p-1)
      COMMON /rand/ ir,j,k
      SAVE rand

      ir(j)= ieor(ir(j),ir(k))
*     ir(j)= ixor(ir(j),ir(k))
      irnd=ir(j)
      frnd=irnd*fnorm
      ranf=frnd

      if (ranf.eq.1.0) then
	 write(*,*) 'ranf = 1.0 !'
         stop 
      endif

      j=j+1
      if(j.ge.p) j=j-p
      k=k+1
      if(k.ge.p) k=k-p

      return
      end

c----------------------------------------------------------------------c
      subroutine cinit3(ndelay)
c----------------------------------------------------------------------c
      call init3
      call delay3(ndelay)
      return
      end

c----------------------------------------------------------------------c
      subroutine init3
c----------------------------------------------------------------------c
      parameter( IP=521, IQ=32 )
      parameter( macrm= 40, macri=1)
      common /rand/ iw(0:IP-1),jr,kr

      integer ib(0:IP-1)
      ix=macri
      do 10 i=0,IP-1
        ix=ix*69069
        ib(i)=ishft(ix,-31)
   10 continue
      jr=0
      kr=IP-IQ
      do 30 j=0,IP-1
        iwork=0
        do 20 i=0,31
          iwork=iwork*2+ib(jr)
          ib(jr)= ieor(ib(jr),ib(kr))
*         ib(jr)= ixor(ib(jr),ib(kr))
          jr=jr+1
          if (jr.eq.IP) jr=0
          kr=kr+1
          if (kr.eq.IP) kr=0
   20   continue
        iw(j)=ishft(iwork,-1)
   30 continue

      end

c-----------------------------------c
      subroutine delay3(lambda)
c-----------------------------------c
      parameter( IP=521, IQ=32 )
      parameter( macrm= 40, macri=1)
      common /rand/ iw(0:IP-1),jr,kr
      integer iwk(0:2*IP-2),c(0:2*IP-1),ib(0:IP+31)

      mu=macrm
      do 110 i=0,IP-1
        iwk(i)=iw(i)
  110 continue
      do 120 i=IP,2*IP-2
        iwk(i)= ieor(iwk(i-IP),iwk(i-IQ))
*       iwk(i)= ixor(iwk(i-IP),iwk(i-IQ))
  120 continue
      do 210 i=0,mu-1
        ib(i)=0
  210 continue
      m=lambda
      nb=mu-1
  220 continue
        if(m.le.IP-1) goto 300
        nb=nb+1
        ib(nb)=mod(m,2)
        m=m/2
        goto 220
  300 do 310 i=0,IP-1
        c(i)=0
  310 continue
      c(m)=1
      do 340 j=nb,0,-1
        do 320 i=IP-1,0,-1
          c(2*i+ib(j))=c(i)
          c(2*i+1-ib(j))=0
  320   continue
        do 330 i=2*IP-1,IP,-1
          c(i-IP)= ieor(c(i-IP),c(i))
*         c(i-IP)= ixor(c(i-IP),c(i))
          c(i-IQ)= ieor(c(i-IQ),c(i))
*         c(i-IQ)= ixor(c(i-IQ),c(i))
  330   continue
  340 continue
      do 420 j=0,IP-1
        iwork=0
        do 410 i=0,IP-1
          iwork= ieor(iwork,c(i)*iwk(j+i))
*         iwork= ixor(iwork,c(i)*iwk(j+i))
  410   continue
        iw(j)=iwork
  420 continue
      end

c--------------------------------------------------------------c
      subroutine  rgauss1(xr,nup)
c--------------------------------------------------------------c
c     gaussian random number                                   c
c     normalized to unity                                      c
c     notice:  nup should be even.                             c
c--------------------------------------------------------------c
      REAL*8  xr(nup)

      nup2 = nup/2
      if (nup.ne.(2*nup2))  write(6,9000)
 9000 format(' !!!!  attention   nup in rgauss should be even',
     *       '  nup = ',i5)
      
      do 10  i = 1, nup2
      i2 = 2*i
      i1 = i2 - 1

      v1 = sqrt(-log(ranf(idum)+1.e-10)*2.)
      v2 = 6.28318*ranf(idum)
c
      xr(i1) = v1 * cos(v2)
      xr(i2) = v1 * sin(v2)

   10 continue
      return
      end

C-----------------------------------------------------C
      SUBROUTINE rgauss2(xr,nup)
C-----------------------------------------------------C
C     Same as rgauss1 but parallelized.    
C-----------------------------------------------------C
C     generates nv random numbers with normal 
C     distribution by Box-Muller
C     pi2 = pi * 2
C-----------------------------------------------------C
C     Input
C       nup     : Dimension 
C     Working space
C       vwork(NV/NBUSH) (Double precision)
C     Output
C       xr(nup)    (Double precision)
C-----------------------------------------------------C
      include '../INCLUDE/para_geometry'

      REAL*8 PI, PI2
      PARAMETER(PI=3.1415926535897932,PI2=2.d0*PI)
      REAL*8 xr(nup)
      REAL*8 vwork(3*3*NV/NBUSH)
      REAL*8 rho, theta, rand1, rand2

      IF((3*3*NV/NBUSH) < nup) THEN
         WRITE(*,*) "working erea vwork has too small dimension"
         WRITE(*,*) "NV/NBUSH, nup : ", NV/NBUSH,nup
         STOP
      ENDIF

      CALL rands(vwork,nup) 

      nuph = nup / 2

*POPTION PARALLEL
*VOPTION INDEP(xr,vwork)
      DO i = 1, nuph
        rho   = SQRT(-2.d0 * LOG(vwork(2*i-1)+1.e-10) )
        theta = PI2 * vwork(2*i)
        xr(2*i-1)  = rho * COS(theta)
        xr(2*i)    = rho * SIN(theta)
      ENDDO

      IF ((2*nuph) .NE.  nup) THEN 
         rand1 = ranf(idum)
         rand2 = ranf(idum)
         xr(nup) = SQRT(-2.d0 * LOG(rand1) )
     &             * COS(PI2 * rand2)              
      ENDIF

      RETURN
      END
