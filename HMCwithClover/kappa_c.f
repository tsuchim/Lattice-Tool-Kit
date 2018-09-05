!------------------------------------------------------!
      program CalKappac 
!------------------------------------------------------!
!    J.Jansen and R.Sommer, Nucl.Phys. B530 (1998) 185
!------------------------------------------------------!

      REAL*8 beta
      REAL*8 kappa_c, kappa, mass, a_inv

      REAL*8  g02
      REAL*8  n, n0, n2, n4, n6, n8, n10

      mass = 0.15d0 
      a_inv = 2.09385
      beta = 6.05

      g02 = 6.d0/beta
      g04 = g02**2
      g06 = g02**3
      g08 = g02**4
      g10 = g02**5

      n0 = 0.125d0
      n2 = 0.008439857d0
      n4 = 0.0085d0
      n6 = -0.0272d0
      n8 = 0.0420d0
      n10 = -0.0204d0

      n = n0 + n2*g02 + n4*g04 + n6*g06 + n8*g08 + n10*g10

      kappa_c = n

      kappa = 1.d0/(2.d0*mass/a_inv + 1.d0/kappa_c)

      write(*,'(1x,a,e15.7)') "kappa_c = ", kappa_c
      write(*,'(1x,a,e15.7)') "kappa   = ", kappa
      write(*,'(1x,a,e15.7)') "beta    = ", beta
      write(*,'(1x,a,e15.7)') "g^2     = ", g02
      write(*,'(1x,a,e15.7,a,1x)') "mass    = ", mass, '[GeV]' 
      write(*,'(1x,a,e15.7)') "m*a     = ", mass/a_inv  

      stop
      end
