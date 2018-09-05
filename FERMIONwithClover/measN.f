!-----------------------------------------------------------------------
      SUBROUTINE CalNdens(rho,chirho,Nnz)
!-----------------------------------------------------------------------
!    Calculation of number density <n> with Noize method
!    and susceptibility
!    constructed from Nagata's cal_ndens_nz(rho,chirho) 
!-----------------------------------------------------------------------
      USE field_f
      USE fpara
      include '../INCLUDE/para_geometry'

      INTEGER, INTENT(IN) :: Nnz  !  Number of Noize
      COMPLEX*16  zsum, zsum2, zfac1, zfac2
      COMPLEX*16, INTENT(OUT) :: rho, chirho

      TYPE(f_field) zxi, ztmp1, ztmp2, gauss_ran_f
      TYPE(f_field) x1, x2, x3

      Nf = 2  ! Number of flavors

      zfac1 = - hop * exp(+CONJG(cmu))
      zfac2 = - hop * exp(-CONJG(cmu))

      zsum  = (0.d0,0.d0)
      zsum2 = (0.d0,0.d0)
      
      write(*,*) "No. of Noise : ", Nnz 

      DO inz = 1, Nnz

        zxi = gauss_ran_f(idumy)
        iflag = 1
        CALL cg0(x1,zxi,iflag)   ! X1 = D^{-1}*xi

        nu = 4
        ztmp1 = zfac1 * (nu .fshift. x1)
        ztmp2 = zfac2 * ((-nu) .fshift. x1)
        x2 = ztmp1 - ztmp2       ! X2 = D'*X1

        zsum = zsum + zxi * x2   ! zsum = xi_adj * x2 = <xi|D'*D^{-1}|xi>

        iflag = 1
        CALL cg0(x3,x2,iflag)    ! X3 = D^{-1}*X2

        nu = 4 
        ztmp1 = zfac1 * (nu .fshift. x3)
        ztmp2 = zfac2 * ((-nu) .fshift. x3)
        x1 = ztmp1 - ztmp2       ! X1 = D'*X3

        zsum2 = zsum2 + zxi * x1 ! xsum = <xi|D'*D^{-1}*D'*D^{-1}|xi>
      
      ENDDO

      rho = DBLE(Nf) * Real(zsum) /DBLE(Nnz*NX*NY*NZ*NT)
      chirho = DBLE(Nf) * Real(zsum2) /DBLE(Nnz*NX*NY*NZ*NT)

      RETURN
      END 
