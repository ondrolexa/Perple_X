      subroutine gcohx6 (fo2,hybrid)
c----------------------------------------------------------------------
c  program to calculate GCOH fluid properties as a function of XO 
c  see Connolly (1995) and/or coh_speciation_with_ethane.mws for 
c  details. the routine uses a hybrid EoS (pure CO2 and H2O from 
c  Pitzer & Sterner 2004, CH4 from Kerrick & Jacobs 1981, MRK for 
c  all activities and all other fugacities) if hybrid = .true., else
c  uses MRK for all purposes. 

c  modified to account for diamond stability and to use CORK for
c  h2o and co2. 4/27/04, JADC

c  modified to use pseos (pitzer and sterner, 1994) for h2o/co2, 
c  3/5/2015, JADC

c  modified to include ethane and use 2d newton-raphson. 10/28/16, JADC.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(nsp),nit,i

      logical bad, hybrid 

      double precision ghh2o,ghco2,ghch4,kh2o,kch4,kco,kco2,kc2h6,oy5,
     *       c1,c2,c3,c4,t1,t2,t3,t4,m,dm3,dm5,c,dc3,dc5,nh,rat,
     *       dy5,dy3,y5,y3,fo2,ytot,t4y3,t3y3,t2y5,t4y5,det,x

      double precision dinc
      external dinc

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

      double precision units, r13, r23, r43, r59, r1, r2
      common/ cst59 /units, r13, r23, r43, r59, r1, r2

      save ins
      data ins/1,2,3,4,5,16,10*0/
c----------------------------------------------------------------------
      nit = 0
      oy5 = 0d0
      bad = .false.
c                                 setup calculates equilibrium constants
c                                 and computes pure species fugacities
      call setup (ghh2o,ghco2,ghch4,kh2o,kco2,kco,kch4,kc2h6)

      c1 = dexp (kch4) * p
      c2 = dexp (2d0*kc2h6 - 3d0*kch4) *p
      c3 = dexp (kco2 - 2d0*kco) * p 
      c4 = dexp (kh2o - kco) * p
c                                 initial guess, assume near binary
      x = 1d0 + xo 
      rat = xo/(1d0-xo) 

      if (xo.gt.r13) then
         y3 = dsqrt((c3*x*(3d0*xo - 1d0)))/c3/x
         y5 = 2d0 * (1d0 - xo)/c4/y3/x
      else if (xo.eq.r13) then 
         y3 = 1d0/dsqrt(c4)
         y5 = y3 
      else 
         y5 = dsqrt((c1*x*(1d0 - 3d0*xo)))/c1/x
         y3 = 4d0*xo/c4/y5/x
      end if

      do
c                                 iteration loop
         t1 = c1 / g(ins(4)) * g(ins(5)) ** 2
         t2 = dsqrt(c2 * (t1 * g(ins(4))) ** 3) / g(ins(6))
         t3 = c3 / g(ins(2)) * g(ins(3)) ** 2
         t4 = c4 / g(ins(1)) * g(ins(5)) * g(ins(3))

         t4y3 = t4*y3
         t3y3 = t3*y3
         t2y5 = t2*y5
         t4y5 = t4*y5

         nh = ((3d0*t2y5 + 2d0*t1)*y5 + t4y3 + 1d0)*y5
c                                 no/nh
         x = (t3*y3 + (1d0 + t4y5)/2d0)*y3/nh
c                                 mass balance eq
         m = rat - x
c                                 diff(m,y5)
         dm5 = (x*((9d0*t2y5 + 4d0*t1)*y5 + t4y3 + 1d0) - t4y3/2d0)/nh
c                                 diff(m,y3)
         dm3 = (x*t4y5 - 2d0*t3y3 - 0.5d0 - t4y5/2d0)/nh
c                                 closure eq
         c = ((-t2y5 - t1)*y5 - t4y3 - 1d0)*y5 + 1d0 - (t3y3 + 1d0)*y3 
c                                 diff(c,y5)
         dc5 = (-3d0*t2y5 - 2d0*t1)*y5 - t4y3 - 1d0
c                                 diff(c,y3)
         dc3 = -2d0*t3y3 - t4y5 - 1d0
c                                 newton-raphson increments
         det = dc5 * dm3 - dc3 * dm5
         dy5 = -(c * dm3 - dc3 * m) / det
         dy3 =  (c * dm5 - dc5 * m) / det
c                                 add the increment
         y5 = dinc(y5,dy5)
         y3 = dinc(y3,dy3)
c                                 back calculate y's
         y(ins(5)) = y5
         y(ins(3)) = y3
         y(ins(4)) = t1 * y5**2 
         y(ins(6)) = t2 * y5**3
         y(ins(2)) = t3 * y3**2
         y(ins(1)) = t4 * y5 * y3
c                                 compute ytot
         ytot = 0d0

         do i = 1, 6 
            ytot = ytot + y(ins(i))
         end do 
c                                 renormalize
         do i = 1, 6
            y(ins(i)) = y(ins(i)) / ytot
         end do 
c                                 check for convergence, could 
c                                 do better than this
         if (dabs(y5-oy5).lt.nopt(5)*y5) exit
c                                 check if iteration count exceeded
         if (nit.gt.iopt(21)) then 
            call warn (176,y5,nit,'HOCNIK')
            bad = .true.
            exit
         end if
c                                 calculate new fugacity coefficients
         call mrkmix (ins, 6, 1)

         if (hybrid) then 
            g(ins(1)) = ghh2o * g(ins(1))
            g(ins(2)) = ghco2 * g(ins(2)) 
            g(ins(4)) = ghch4 * g(ins(4))
         end if 

         oy5 = y5

         nit = nit + 1

         y5 = y(ins(5))
         y3 = y(ins(3))

      end do 

      if (bad) then 

         vol  = 0d0 
         fh2o = dlog(p*1d4)
         fco2 = fh2o
         fo2  = fh2o
         
      else 

         if (hu.eq.0) then 
c                                 normal fugacities
            fh2o = dlog(g(ins(1))*p*y(ins(1)))
            fco2 = dlog(g(ins(2))*p*y(ins(2)))
            fo2 = 2d0*(dlog(g(ins(3))*p*y3) - kco)

         else 
c                                 projecting through graphite
            fh2o = dlog(g(ins(5))*p*y5)
            fco2 = 2d0*(dlog(g(ins(3))*p*y3) - kco)

         end if

      end if 

      if (hybrid) 
     *   vol = vol + y(ins(1))*vm(1) + y(ins(2))*vm(2) + y(ins(4))*vm(3)

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