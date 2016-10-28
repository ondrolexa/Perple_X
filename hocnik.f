      subroutine hocnik (fo2)
c----------------------------------------------------------------------
c  program to calculate GCOH fluid properties as a function of XO 
c  using HSMRK/MRK hybrid see Connolly and Cesare (1992) for details.

c  modified to account for diamond stability and to use CORK for
c  h2o and co2, 4/27/04, JADC.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(nsp),nit,ier

      double precision ghh2o,ghco2,ghch4,kh2o,kch4,kco,kco2,kc2h6

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      double precision gmh2o,gmco2,gmch4,vm
      common/ cstchx /gmh2o,gmco2,gmch4,vm(3)

      double precision y,g,v
      common / cstcoh /y(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer ibuf,hu,hv,hw,hx   
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save ins
      data ins/ 5,4,1,3,2,16,10*0/
c----------------------------------------------------------------------

      nit = 0
      oy1 = 0d0

      call setup (ghh2o,ghco2,ghch4,kh2o,kco2,kco,kch4,kc2h6)

      c1 = dexp (kch4) * p
      c2 = dexp (2d0*kc2h6 - 3d0*kch4) *p
      c3 = dexp (kco2 - 2d0*kco) * p 
      c4 = dexp (kh2o - kco) * p

      do 

         t1 = c1 / phi(2) * phi(1) ** 2
         t2 = dsqrt(c2 * (t1 * phi(2)) ** 3) / phi(6)
         t3 = c3 / phi(5) * phi(4) ** 2
         t4 = c4 / phi(3) * phi(1) * phi(4)
c                                 mass balance eq
         m = rat - (t3 * y(4) ** 2 + 2 * t4 * y(1) * y(4) + 2 * y(4)) /(
     #3 * t2 * y(1) ** 3 + 2 * t1 * y(1) ** 2 + t4 * y(1) * y(4) + y(1))
c                                 diff(m,y1)
         dm1 = -2 * t4 * y(4) / (3 * t2 * y(1) ** 3 + 2 * t1 * y(1)**2+ 
     #t4 * y(1) * y(4) + y(1)) + (t3 * y(4) ** 2 + 2 * t4 * y(1) * y(4) 
     #+ 2 * y(4)) / (3 * t2 * y(1) ** 3 + 2 * t1 * y(1) ** 2 + t4 * y(1)
     # * y(4) + y(1)) ** 2 * (9 * t2 * y(1) ** 2 + 4 * t1 * y(1) + t4 * 
     #y(4) + 1)
c                                 diff(m,y4)
         dm4 = -(2 * t3 * y(4) + 2 * t4 * y(1) + 2) / (3 * t2 *y(1)**3+
     # 2 * t1 * y(1) ** 2 + t4 * y(1) * y(4) + y(1)) + (t3 * y(4) ** 2 +
     # 2 * t4 * y(1) * y(4) + 2 * y(4)) / (3 * t2 * y(1) ** 3 + 2 * t1 *
     # y(1) ** 2 + t4 * y(1) * y(4) + y(1)) ** 2 * t4 * y(1)
c                                 closure eq
         c = -t2 * y(1) ** 3 - t1 * y(1) ** 2 - t3 * y(4) ** 2 - t4*y(1)
     # * y(4) - y(1) - y(4) + 1
c                                 diff(c,y1)
         dc1 = -3 * t2 * y(1) ** 2 - 2 * t1 * y(1) - t4 * y(4) - 1
c                                 diff(c,y4)
         dc4 = -2 * t3 * y(4) - t4 * y(1) - 1
c                                 newton-raphson increments
         dy1 = -(c * dm4 - dc4 * m) / (dc1 * dm4 - dc4 * dm1)
         dy2 =  (c * dm1 - dc1 * m) / (dc1 * dm4 - dc4 * dm1)
c                                 increment
         y1 = dinc(y1,dy1)
         y4 = dinc(y4,dy4)
c                                 back calculate y's
         y(1) = y1
         y(4) = y4
         y(2) = t1 * y1**2 
         y(6) = t2 * y1**3
         y(5) = t3 * y4**2
         y(3) = t4 * y1 * y4
c                                 compute ytot
         ytot = 0d0

         do i = 1, 6 
            ytot = ytot + y(i)
         end do 
c                                 renormalize
         do i = 1, 6
            y(i) = y(i) / ytot
         end do 

         nit = nit + 1

         if (dabs(y1-oy1).lt.nopt(5)*y1) exit

         if (nit.gt.iopt(21)) then 
            call warn (176,xh2o,nit,'HOCNIK')
            bad = .true.
            exit          
         end if

         call mrkmix (ins, 5, 1)

         g(3) = ghh2o * g(3)
         g(5) = ghco2 * g(5) 
         g(2) = ghch4 * g(2)

         oy1 = y1 

      end do 

      if (bad) then 

         vol = 0d0 
         fh2o = dlog(p*1d4)
         fco2 = fh2o
         
      else 

         if (hu.eq.0) then 

            fh2o = dlog(g(3)*p*y(3))
            fco2 = dlog(g(5)*p*y(5))
            fo2 = 2d0*(dlog(g(4)*p*y4) - kco)

         else 

            fh2o = dlog(g(1)*p*y1)
            fco2 = 2d0*(dlog(gco*p*y4) - kco)

         end if 

         vol = vol + y(3) * vm(1) + y(5) * vm(2) + y(2) * vm(3)

      end if 

      end

      double precision function dinc (y,dy)
c----------------------------------------------------------------------
c function to increment 0 < y < 1 by dy, if the increment would cause
c y to exceed its bounds, the increment is taken to be half the 
c interval to the bound. 10/2016, JADC.
c----------------------------------------------------------------------
      implicit none

      double precision y, dy
c----------------------------------------------------------------------
      if (y+dy.ge.1d0) then
         dinc = 0.5d0 + y / 2d0
      else if (y+dy.le.0d0) then 
         dinc = y / 2d0
      else 
         dinc = y + dy
      end if 

      end 