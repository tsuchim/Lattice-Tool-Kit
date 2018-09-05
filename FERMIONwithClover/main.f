c main.f (for fermion) 
c----------------------------------------------------------------------
      program ftest 
c-----------------------------------------------------------------------
      USE field_f
      USE fpara
C     include 'mpif.h'                                                ! MPI
      include '../INCLUDE/para_geometry'
C     include '../INCLUDE/para.h'                                     ! MPI

      TYPE(g_field0) u(4)
      common/ config/ u
      TYPE(g_field1) gtemp

      REAL*8 beta
      TYPE(f_field) x, b

      INTEGER :: Nnz  !  Number of Noize
      COMPLEX*16  :: rho, chirho

      COMPLEX*16 ci  ! For test
      ci = (0.d0,1.d0) ! For test

C     myrank = 0  ! dummy for non-MPI case                              ! MPI
c     ...  MPI Preparations  ........................                   ! MPI
C     call mpi_init(ierr)                                               ! MPI
C     call mpi_comm_size (mpi_comm_world,nprocs,ierr)                   ! MPI
C     call mpi_comm_rank (mpi_comm_world,myrank,ierr)                   ! MPI
c     ...............................................                   ! MPI

c     *****   prepare table   *****                                     ! MPI
C     call make_table                                                   ! MPI
C    &     (itable,iprocs,jprocs,myrank,myranki,myrankj,mpi_proc_null)  ! MPI

      call  rdparam(hop,r,Csw,cmu,eps,imax,istart) 

c     *****  Set hopping parameters with chemical potential
      do mu = 1, 3
        hopp(mu) = hop
        hopm(mu) = hop
      enddo

      hopp(4) = exp(+cmu)*hop
      hopm(4) = exp(-cmu)*hop

c     ...  Check paramaters NBUSH and NDW  ..
      IF(NBUSH.NE.1) THEN
         write(*,*) "This is ferimon program. NBUSH is not 1:",NBUSH
         stop
      ENDIF
      IF(NDW.NE.1)  THEN
         write(*,*) "This is ferimon program. NDW is not 1:",NDW
         stop
      ENDIF

c     ...  Initialize the gauge field
      call  init(u,istart,ndelay)

c     ...  Make Gamma Matrices
      call mk_gamma

c     .....  Calculate meson propagators
      call propa

!     Calculation of number density <n> with Noize method
!     and susceptibility
      Nnz = 200
      call CalNdens(rho,chirho,Nnz)
      write(*,*) "rho, chirho: ", rho, chirho

c     ...  MPI Finalize   ...........................                 ! MPI
C9999 call mpi_finalize(ierr)                                         ! MPI
c     ...............................................                 ! MPI

      end

c----------------------------------------------------------------------c
      subroutine  rdparam(hop,r,Csw,cmu,eps,imax,istart) 
c----------------------------------------------------------------------c
*     Subrouitine to read parameters
c----------------------------------------------------------------------c
C     include 'mpif.h'                                                ! MPI
      include '../INCLUDE/para_geometry'
C     include '../INCLUDE/para.h'                                     ! MPI 
      REAL*8 hop, r, Csw, eps, beta
      COMPLEX*16 cmu
      LOGICAL clover

C     IF (myrank == 0) THEN                                           ! MPI  
        read(*,*) hop, r
        read(*,*) clover 
        read(*,*) beta 
        read(*,*) cmu 
        read(*,*) istart 
        read(*,*) imax 
        read(*,*) eps 

*       Csw = 0.d0
*       IF(clover) call CswEstimate1(beta,Csw)
        Csw = 1.5504d0

C     ENDIF                                                           ! MPI

C     CALL MPI_BCAST(hop,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)           ! MPI
C     CALL MPI_BCAST(r,  1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)           ! MPI
C     CALL MPI_BCAST(Csw,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)           ! MPI
C     CALL MPI_BCAST(cmu,1,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)       ! MPI
C     CALL MPI_BCAST(eps,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)           ! MPI
C     CALL MPI_BCAST(istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      ! MPI
C     CALL MPI_BCAST(imax,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      ! MPI 

C     IF (myrank==0) THEN                                             ! MPI
        write(*,'(a,f7.4,2x,f5.1)') "hop, r : ",  hop, r
        write(*,'(a,f7.4)')          "Csw : ",    Csw
        write(*,'(a,f7.4,1x,f7.4)')          "cmu : ",    cmu
        write(*,'(a,i2)')            "istart : ", istart 
        write(*,'(a,i6)')            "imax : ",   imax 
        write(*,'(a,e12.3)')         "eps : ",    eps 
        write(*,'(a,4i3)')           "Nx,Ny,Nz,Nt: ", NX,NY,NZ,NT 
C     ENDIF                                                            ! MPI

      IF(NBUSH.NE.1)  THEN
        write(*,*) "This is Fermion program, and NBUSH must be 1"
        write(*,*) "But your NBUSH = ", NBUSH
        stop
      ENDIF
      IF(NDW.NE.1)  THEN
        write(*,*) "This is Fermion program, and NDW must be 1"
        write(*,*) "But your NDW = ", NDW
        stop
      ENDIF

      return
      end
