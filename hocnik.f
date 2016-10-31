      subroutine gcohx6 (fo2,hybrid)
c----------------------------------------------------------------------
c  program to calculate GCOH fluid properties as a function of XO 
c  see Connolly (1995) and/or coh_speciation_with_ethane.mws for 
c  details. the routine uses a hybrid EoS (pure CO2 and H2O from 
c  Pitzer & Sterner 2004, CH4 from Kerrick & Jacobs 1981, MRK for 
c  all activities and all other fugacities) if hybrid = .true., else
c  uses MRK for all purposes.

c  replaces hocgra and hocmrk.

c  10/28/16, JADC.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ins(6),jns(3),isp,jsp,nit,i

      logical bad, hybrid 

      double precision oy5,fo2,ytot,t4y3,t3y3,t2y5,t4y5,det,x,dy5,dy3,
     *       c1,c2,c3,c4,t1,t2,t3,t4,m,dm3,dm5,c,dc3,dc5,nh,rat,y5,y3

      double precision dinc
      external dinc

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      double precision gh,vh
      common/ csthyb /gh(nsp),vh(nsp)

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

      double precision eqk
      common / csteqk /eqk(nsp)

      save isp, jsp, ins, jns
      data isp, jsp, ins, jns/6,3,1,2,3,4,5,16,1,2,4/
c----------------------------------------------------------------------
      nit = 0
      oy5 = 0d0
      bad = .false.
c                                 check for in bounds composition
      call xochk 
c                                 compute equilibrium constants in csteqk
      call seteqk (ins,isp,elag)
c                                 compute pure mrk fluid properties
      call mrkpur (ins,isp)
c                                 compute hybrid pure fluid props
      if (hybrid) call hybeos (jns,jsp)

      c1 = dexp (eqk(4)) * p
      c2 = dexp (2d0*eqk(16) - 3d0*eqk(4)) * p
      c3 = dexp (eqk(2) - 2d0*eqk(3)) * p 
      c4 = dexp (eqk(1) - eqk(3)) * p
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
            do i = 1, jsp 
               g(jns(i)) = gh(jns(i)) * g(jns(i))
            end do 
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
            fo2 = 2d0*(dlog(g(ins(3))*p*y3) - eqk(3))

         else 
c                                 projecting through graphite
            fh2o = dlog(g(ins(5))*p*y5)
            fco2 = 2d0*(dlog(g(ins(3))*p*y3) - eqk(3))

         end if

      end if 

      if (hybrid) then
         do i = 1, jsp 
            vol = vol + y(jns(i))*vh(jns(i))
         end do 
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

      subroutine seteqk (ins,isp,ac)
c----------------------------------------------------------------------
c program to compute standard (from the elements convention) ln equilibrium 
c constants for 1 mole of the nsp fluid species. 

c ac - ln(a[graphite])
c 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision dg,t2,t3,ac

      integer i,j,isp,ins(*)

      double precision eqk
      common / csteqk /eqk(nsp)

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

c----------------------------------------------------------------------
      t2 = t * t
      t3 = t2 * t
c                                correct activity of graphite
c                                for diamond stability if necessary:
      call dimond (dg)
c                                graphite pressure effect and activity
c                                corrections:
      dg = dg + ac + p*( 1.8042d-06 + (0.058345d0 - 8.42d-08*p)/t ) 
c                                graphite activity effect, note for
c                                certain routines (xoxsrk, ifug = 19/20
c                                the agph term should not be added)!!!
      do i = 1, isp 

         j = ins(i)

         if (j.eq.1) then 
c                                h2o/robie
            eqk(j) = -7.028214449d0 + 30607.34044d0/t  
     *               - 475034.4632d0/t2 + 50879842.55d0/t3
         else if (j.eq.2) then 
c                                co2/robie
            eqk(j) = .04078341613d0 + 47681.676177d0/t 
     *              - 134662.1904d0/t2 + 17015794.31d0/t3 + dg
         else if (j.eq.3) then 
c                                co/robie
            eqk(j) = 10.32730663d0  + 14062.7396777d0/t
     *              - 371237.1571d0/t2 + 53515365.95d0/t3 + dg
         else if (j.eq.4) then 
c                                ch4/robie
            eqk(j) = -13.86241656d0 + 12309.03706d0/t
     *             - 879314.7005d0/t2 + .7754138439d8/t3 + dg
         else if (j.eq.6) then 
c                                h2s/ohmoto and kerrick
            eq(j) = 10115.3d0/t - 0.791d0 * dlog (t) + 0.30164d0
         else if (j.eq.8) then 
c                                so2/ohmoto and kerrick
            eq(j) = 43585.63147d0/t - 8.710679055d0
         else if (j.eq.9) then 
c                                cos/ohmoto and kerrick
            eq(j) = 10893.52964d0/t - 9.988613730d0
         else if (j.eq.16) then 
c                                c2h6/HSC, 10/2016
            eqk(j) =  4.09702552d7/t3 - 8.01186095d5/t2 + 1.39350247d4/t
     *                      - 2.64306669d1 + 2d0 * dg
         end if 

      end do 

      end 

      subroutine dimon1 (agph)
c-----------------------------------------------------------------------
c dimon1 tests if p-t conditions are in the diamond stability field
c if they are it computes the activity of graphite needed to represent
c diamond stability for C-O-H fluid speciation routines. 

c dimon1 should replace dimond which automatically adds a graphite
c activity correction

c  agph - natural log of graphite activity to be added to graphite
c        free energy to equate g(diamond) = g(graphite). 

c polynomial functions fit by S. Poli, 4/27/04
c                                 JADC, 4/27/04
c-----------------------------------------------------------------------
      implicit none 

      double precision ptrans,agph

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

c                                 get transition pressure
c                                 (should really use the G
c                                 function for consistency):

      ptrans = 5284.165053d0 + (33.21515773d0 - .002106330992d0*t)*t

      if (p.gt.ptrans) then 
c                                 compute corrected graphite activity
         agph = 0.008423508384179629d0 
     *        + (4.693008650307614d-11 * p - 3.850380793502567d-5)*p
     *        + (1.4126916053951515d-3 + 1.393226795939807d-8*p 
     *        - 5.887505938975768d-7 * t)*t

      end if 

      end 

      subroutine hybeos (jns, jsp)
c---------------------------------------------------------------------
c set up routine for hybrid fluid EoS calculations. computes the 
c (unecessay?) delta volumes and pure fluid fugacity coefficient 
c rations used to convert the mrk fugacities to hybrid fugacities.

c the choice of the pure fluid eos are hard wired (pseos, pitzer & sterner
c 1994 for h2o/co2; hsmrk, kerrick and jacobs 1981 for methane).

c this routine is to replace hscrkp and cstchx, 10/2015 JADC
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,jns(*),jsp

      double precision hsfch4,fg(nsp)

      external hsfch4
 
      double precision gh,vh
      common/ csthyb /gh(nsp),vh(nsp)

      double precision x,g,v
      common/ cstcoh /x(nsp),g(nsp),v(nsp)

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      do i = 1, jsp

         j = jns(i)

         vh(j) = -v(j)
         gh(j) = g(j)
         
         if (j.le.2) then 
c                                 water/co2, pitzer and sterner 1994:
            call pseos (v(j),fg(j),j)

         else if (j.eq.4) then 
c                                 methane hsmrk kerrick and jacobs 191.
            fg(j) = hsfch4 (v(j))

         else

            stop

         end if 
c                                 the fugacity coefficient of the pure gas
         g(j) = dexp(fg(j))/p
c                                 the hybrid delta volume (hyb-mrk), it's 
c                                 doubtful this thing is really used, if it 
c                                 is it must be in fluids.
         vh(j) = vh(j) + v(j)
c                                 the hybrid/mrk pure fluid fugacity ratio
         gh(j) = g(j)/gh(j)

      end do 

      end

      subroutine xochk
c----------------------------------------------------------------------
c xo check for speciation routines that can't handle xo = 0, 1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision p,t,xo,u1,u2,tr,pr,r,ps
      common / cst5 /p,t,xo,u1,u2,tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                check if xo is <1, >0,
c                                reset if necessary
      if (xo.lt.nopt(5)) then
         xo = nopt(5)
      else if (dabs(xo-1d0).lt.nopt(5)) then
         xo = 1d0 - nopt(5)
      end if 

      end 