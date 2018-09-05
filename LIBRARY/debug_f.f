c debug_f.f
c--------------------------------------------c
c     Debugging routines for fermion
c--------------------------------------------c
      subroutine show_f(x,ia,mu)
c--------------------------------------------c
c     Show fermion field x                   c
c                  ia:Color, mu:Dirac        c
c--------------------------------------------c
      USE field_f
      include '../INCLUDE/para_geometry'

      TYPE(f_field)  x

      write(*,*) "Color : ", ia
      write(*,*) "Dirac : ", mu

      do it = 1, NT
      do iz = 1, NZ
      do iy = 1, NY
        write(*,'(a,3i2,a,i2)') "iy,iz,it:",iy,iz,it," ix=1,",NX
        write(*,1000)  (x%f(ia,ix,iy,iz,it,mu), ix=1,NX)
      enddo
      enddo
      enddo

 1000 format(8("(",f6.3,",",f6.3,")",2x))
      
      return
      end
