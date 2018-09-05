c--------------------------------------------------c
      PROGRAM mkdata
c--------------------------------------------------c
c     Estimate error bars using Jack-knife         c
c--------------------------------------------------c
c     NMAX : the possible maximum number of data   c
c     nclust : # of clusters                       c
c     ntherm : # of thermalization sweeps          c
c--------------------------------------------------c
      PARAMETER( NMAX=10000 )

      REAL  x(NMAX)

      nclust = 4
      ntherm = 0 

      write(*,*) '# of clusters : ', nclust, 
     &         '  # of thermalization sweeps : ', ntherm

      ndata = 0

    1 continue
      ndata = ndata + 1
      read(*,*,end=1000)  x(ndata)
      goto 1

 1000 continue
      ndata = ndata - 1
      write(*,*) 'ndata : ', ndata
      call  jack(x,ave,errbar,ndata,nclust,ntherm)
      write(*,2100)  ' x = ', ave, ' +/- ', errbar

 2100 format(1x,a,e15.7,a,e15.7)
 2000 continue

      stop
      end 
      
c-----------------------------------------------------------c
      subroutine jack(x,ave, errbar,ndata,nclust,ntherm)
c-----------------------------------------------------------c
c     Data analysis program   c
c     using Jack-Knife        c
c        ndata : # of data
c        nclust : # of clusters
c        ntherm : # of thermalization sweeps
c        nbin   : n_bin
c-----------------------------------------------------------c
      parameter(nclmax=1000)
      real x(ndata)
      real xcl(nclmax), work(nclmax)
      real*8 sum

      if(nclust.gt.nclmax) then
         write(*,*) 'Please increase nclmax to nclust'
         write(*,*) 'nclmax : ', nclmax, ' nclust : ', nclust
         stop
      endif

      nbin  = (ndata-ntherm)/nclust
      ncheck = nbin*nclust + ntherm
      if(ncheck.ne.ndata)  then 
          write(*,*) '(ndata-ntherm) is not a multiple of nclust !'
          write(*,*) 'ndata, ntherm, nclust : ', ndata, ntherm, nclust
          stop
      endif

      m = 0
      mclust = 0
      sum = 0.0d0
      do 1  n = 1, ndata
      if(n.le.ntherm)  goto 1
      m = m + 1
      sum = sum + x(n)
      if(m.eq.nbin) then
         mclust = mclust + 1
         xcl(mclust) = sum/nbin         
         sum = 0.0d0
         m = 0
      endif
    1 continue

      if(nclust.ne.mclust)  then
         write(*,*) 'nclust is not equal to mclust ??'
         write(*,*) 'nclust, mclust : ', nclust, mclust
         stop
      endif

      call step1(xcl,work,nclust)
      call step2(work,ave,sig,nclust)
      
      errbar = sqrt((nclust-1)*sig)

c     write(*,5000)  ave, errbar
c5000 format(1x,'<x> = ',e15.7,2x,'+/-',2x,e15.7)

      return
      end
       
c-------------------------------------------------------------c
      subroutine  step1(xcl,xtild,nclust)
c-------------------------------------------------------------c
c     xtild(1)      = <       xcl(2)+xcl(3)+ ...,+xcl(nclust)>
c     xtild(2)      = <xcl(1)+      + ...,       +xcl(nclust)>
c          .....
c     xtild(nclust) = <xcl(1)+xcl(2)+xcl(3)+ ...,            >
c-------------------------------------------------------------c
      real  xcl(nclust), xtild(nclust)
      real*8  sum

      do 10  i = 1, nclust
      sum = 0.0d0
      do 1000  j = 1, nclust
      if(i.ne.j)  then
         sum = sum + xcl(j)
      endif
 1000 continue
      xtild(i) = sum/(nclust-1)
   10 continue

      return
      end

c--------------------------------------c
      subroutine step2(x,ave,sig,n)
c--------------------------------------c
c     ave : <x>                        c
c     sig : <(x-<x>)**2>               c
c--------------------------------------c
      real x(n)
      real*8 sum1, sum2

      sum1 = 0.0d0
      do 10 i = 1, n
   10 sum1 = sum1 + x(i)
      sum1 = sum1/n

      sum2 = 0.0d0
      do 20  i = 1, n
   20 sum2 = sum2 + (x(i)-sum1)**2
      sum2 = sum2/n

      ave = sum1
      sig = sum2

      return
      end
