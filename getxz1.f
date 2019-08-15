c   getxz - as defined in dlib.f is used in werami
c   getxz - as defined here (getxz1.f) is used in vertex, meemum, convex
c----------------------------------------------------------------------
      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the zcoor array loaded in resub.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jd, id, ids, icoor

      double precision xt
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18)
c                                  xcoordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,h4,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h4,h9)
c----------------------------------------------------------------------
c      icoor = jcoor(id)

      do i = 1, istg(ids,1)

         xt = 0d0 

         do j = 1, ispg(ids,1,i) - 1
            icoor = icoor + 1
            x(1,i,j) = zco(icoor)
            x3(jd,1,i,j) = zco(icoor)
            xt = xt + zco(icoor)
         end do 

         xt = 1d0 - xt
         x(1,i,j) = xt
         x3(jd,1,i,j) = xt

      end do 

      end