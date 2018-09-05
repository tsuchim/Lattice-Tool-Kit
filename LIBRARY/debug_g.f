c debug_g.f
c------------------------------------------------------------c
c     Debugging routines for gauge
c------------------------------------------------------------c
      subroutine  ckunit (u,icheck)
c------------------------------------------------------------c
c     unitarity check                                        c
c     if norm(ui) > eps or i det(u)-1 i  > eps ,             c
c     then icheck = 1  and print some informations, where    c
c         ui =  u * (h.c.of u) .                             c
c     eps is given in data statement.                        c
c------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field1) u
      complex*16  ui(3,3), det
      data  eps/ 1.0e-04/

      icheck = 0

      ncount = 0

      do 10  n = 1, NV/NBUSH

      do 100  i = 1, 3
      do 100  j = 1, 3
      ui(i,j) = u%g(i,1,n) * conjg( u%g(j,1,n) )
     *        + u%g(i,2,n) * conjg( u%g(j,2,n) )
     *        + u%g(i,3,n) * conjg( u%g(j,3,n) )
  100 continue

      unorm = abs(ui(1,1)-1.) + abs(ui(2,2)-1.) + abs(ui(3,3)-1.)
     *      + abs(ui(1,2)) + abs(ui(1,3))
     *      + abs(ui(2,1)) + abs(ui(2,3))
     *      + abs(ui(3,1)) + abs(ui(3,2))

      det = u%g(1,1,n) * ( u%g(2,2,n)*u%g(3,3,n)-u%g(2,3,n)*u%g(3,2,n) )
     *    - u%g(2,1,n) * ( u%g(1,2,n)*u%g(3,3,n)-u%g(1,3,n)*u%g(3,2,n) )
     *    + u%g(3,1,n) * ( u%g(1,2,n)*u%g(2,3,n)-u%g(1,3,n)*u%g(2,2,n) )

      abdet = abs( det - 1.d0 )

      if ( (unorm>eps) .or. (abdet>eps) )  then

        if(ncount<20)  then

         write(6,500)  n, unorm,det, ( (u%g(i,j,n),j=1,3),i=1,3 )
  500    format(/1x,'!!!!   attention   the unitarity is violated at',
     *       i5,'-th matrix.',/1x,'!!!!  unorm = ',e15.5,5x,
     *       'det = (',e15.5,' , ',e15.5,' )',/1x,'!!!!   u = ',
     *        3(/1x,'!!!!',10x,3(e15.5,3x,e15.5,6x)) )
         icheck = 1

        endif

        ncount = ncount + 1

      end if

   10 continue

      return
      end

c------------------------------------------------------------c
c     set test configuration to U for the debug
c------------------------------------------------------------c
      subroutine set_u_test(u)
c------------------------------------------------------------c
      USE field_g
      include '../INCLUDE/para_geometry'

      TYPE(g_field0) u(4)

      do mu = 1, 4

        do it = 1, NT
        do iz = 1, NZ
        do iy = 1, NY
        do ix = 1, NX

          do ic = 1, 3
          do jc = 1, 3
            u(mu)%g(ic,jc,ix,iy,iz,it) = 0.d0
          enddo
          enddo

          u(mu)%g(1,1,ix,iy,iz,it) = 1000*ix + 100*iy + 10*iz + it
     &                             + 0.1*mu
          u(mu)%direction = mu

        enddo
        enddo
        enddo
        enddo

      enddo

      end

