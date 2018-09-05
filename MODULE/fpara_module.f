* fpara_module.f
c--------------------------------------------------------------------c
      MODULE fpara
c--------------------------------------------------------------------c
      REAL*8 hop, r           ! hop: Hopping parameter, r: Wilson term
      REAL*8 Csw              ! Clover term coefficient 
*     REAL*8 cmu              ! chemical potential
*     REAL*8 hopp(4), hopm(4) ! hopp(1)=hopp(2)=hopp(3)=hop
                              ! hopm(1)=hopm(2)=hopm(3)=hop
                              ! hopp(4)=exp(+cmu)*hop, hopm(4)=exp(-cmu)*hop
      COMPLEX*16 cmu              
      COMPLEX*16 hopp(4), hopm(4)
      REAL*8 eps              ! stopping condition for CG
      REAL*8 BC(4)            ! Boundary condition
      INTEGER imax
      INTEGER flagMD          ! flag for Molycular Dynamics

      DATA BC/+1.d0,+1.d0,+1.d0,-1.d0/

      COMPLEX*16 gamma(4,4,5), g5(4,4), rpg(4,4,4), rmg(4,4,4)
      !  gamma(alpha,beta,mu) = gamma_mu(alpha,beta)
      !  gamma(alpha,beta,5) = g5(alpha,beta)
      !  rpg = r + gamma, rmg = r - gamma

      END MODULE  
