      subroutine browop           
c-----------------------------------------------------------------------
c  Brodholt & Wood, 1993 equation of state for pure water and
c  the hsmrk for co2 and the mrk for activities of h2o and co2 in mixtures
c  this and ancillary subroutines were written by 
c  Stefano Poli, Milan, 1996. modified, but untested, July 96, JADC.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision a1,a2,a3,a4,c1,c2,c3,c4,b1,b2,a,
     *                 rt,t12,t15,t25,b,e0,e1,e3,v1,pdv,cor

      integer ins(1),i

      double precision x,g,v
      common/ cstcoh /x(nsp),g(nsp),v(nsp)

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,r,ps

      double precision vol
      common/ cst26 /vol

      save a1,a2,a3,a4,c1,c2,c3,c4,b1,b2,ins

      data a1,a2,a3,a4,c1,c2,c3,c4,b1,b2/-582468d0,
     *     -3038.79d0,-9.24574d-3,3.02674d9,36490.5d0,-1.02451d7,
     *     -1.79681d8,2.18437d9,-3.90463d-2,-0.991078d0/

      data ins/1/
c-----------------------------------------------------------------------
      if (p.lt.1d4) then

         write (*,*) 'Brodholt & Wood EoS invalid at P < 10000 b'
         stop

      end if

      call mrkpur (ins, 1)

      vol = v(1)

      do i = 1, 100

         b = b1 + b2*vol
         e0 = vol - b
         e1 = b + vol
         e3 = 1d0 + b2

         cor = (rt/e0+(c1-a/t12/e1+(c2+(c3+c4/vol)/vol)/vol)/vol - p) 
     *         / 
     *         ((((-4d0*c4/vol-3d0*c3)/vol-2d0*c2)/vol-c1)/vol**2
     *         -rt*e3/e0/e0 
     *         +(e3/e1+1d0/vol)/vol/e1*(a4/t25+a2*t12+a1/t12+a3*t15))

         vol = vol - cor

         if (dabs(cor).lt.1d-3) exit

      end do 

      b = b1 + b2*vol

      fh2o =  (vol*p - v1*1d4 - a/t12/b1*dlog((vol+b)/vol)
     *    - c1*dlog(vol) + (c2 + (c3/2d0 + c4/vol/3d0)/vol)/vol
     *    + pdv ) / rt - dlog(vol-b)/(1d0-b2)

      end


      subroutine hcrkp 
c----------------------------------------------------------------------
c  hcrk calculates the log of water and co2 fugacities using the
c  Halbach Chatterjee 82 equation of state for pure water and
c  the hsmrk for co2 and the mrk for activities of h2o and co2 in 
c  mixtures. this and ancillary subroutines were written by 
c  Stefano Poli, Milan, 1996. modified, but untested, July 96, JADC.
c----------------------------------------------------------------------
      implicit none

      double precision rcm,pdisc,p0,a,b,pincr,volhc,vdp,aincr,c,d,pres

      integer nincr,i

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      double precision p,t,xc,u1,u2,tr,pr,rr,ps
      common/ cst5 /p,t,xc,u1,u2,tr,pr,rr,ps

      save rcm, pdisc, nincr

      data rcm, pdisc, nincr/ 83.144126d0, 1d3, 1000/
c----------------------------------------------------------------------

      if (p.lt.pdisc) then

         write (*,*) 'Halbach & Chatterjee EoS invalid at P < 1000 b'
         stop

      end if
c                        now Halbach & Chatterjee part
      a = 2d0
      b = 4d0
c                        Newton-Raphson numerical integration
      pincr = (p - pdisc) / nincr
      vdp =  pincr * (volhc(pdisc,t) - volhc(p,t))/3d0 
      aincr = pincr 

      do i = 1, nincr - 1
         d = b
         c = a
         b = c
         a = d
         pres = pdisc + aincr
         vdp = vdp + a*pincr*volhc(pres,t)/3d0
         aincr = aincr + pincr
      end do 

         fh2o =  vdp/rcm/t + fh2o

9999  end

      function volhc (p,t)
c--------------------------------------------------------------
      implicit none

      integer iroots, ineg, ipos

      double precision ev(3),rcm,a1,a2,a3,b1,b2,b3,b4,b5,b6,p,t,
     *                 volhc,t12,ahc,bhc,c1,c2,c3,vmin,vmax

      save rcm,a1,a2,a3,b1,b2,b3,b4,b5,b6

      data rcm,a1,a2,a3,b1,b2,b3,b4,b5,b6/83.144126d0,1.616d8,-4.989d4,
     *     -7.358d9,3.4505d-4,3.898d-9,-2.7756d-15,6.3944d-2,2.3776d-5,
     *     4.5717d-10/
c--------------------------------------------------------------
      t12 = dsqrt (t)

      ahc = (a1 + a2*t + a3/t)/ p / t12

      bhc = (1d0 + p*(b1 + p*(b2 + b3*p))) 
     *      / (b4 + p*(b5 + b6*p))

      c1 = -rcm*t / p   
      c2 = ahc - rcm*bhc/p - bhc**2
      c3 = -ahc*bhc

      call roots3 (c1,c2,c3,ev,vmin,vmax,iroots,ineg,ipos)

      if (iroots.eq.3) then
         volhc = vmax
      else
         volhc = ev(1)
      end if

      end



      subroutine delany
c---------------------------------------------------------------------
c  delany calculates ln(fh2o) at P>10 KB using the delaney and
c  helgeson fit (AJS ca 74/78?), routine by S. Poli. this routine
c  is returning delta fh2o!
c-----------------------------------------------------------------------
      implicit none

      double precision p0,b(5,5),c(5,5),gh2o(2),tc

      integer j,l

      double precision p,tk,xc,u1,u2,tr,pr,r,ps
      common / cst5 /p,tk,xc,u1,u2,tr,pr,r,ps

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      save p0, b, c

c     data for delany & helgeson 1974 eos for h2o above 10 kbar

      data b/-5.6130073d0,-1.5285559, -2.6092451d0, 1.7140501d0,
     *  -6.0126987d0,3.8101798d0,1.3752390d0,3.5988857d0,-1.6860893d0,
     *  0d0,-2.1167697d0,-1.5586868d0,-2.7916588d0,0d0,0d0,2.0266445d0,
     *  6.6329577d0, 0d0, 0d0, 0d0,-8.3225572d0,0d0, 0d0, 0d0, 0d0/
c
      data c/4d0,1d0,-2d0,-5d0,-9d0,-1d0,-4d0,-8d0,-11d0,0d0,-6d0,-9d0,
     *    -14d0,0d0,0d0,-11d0,-15d0,0d0,0d0,0d0,-17d0,0d0,0d0,0d0,0d0/
 
      if (p.lt.1d4) then

         write (*,*) 'Delany & Helgeson EoS invalid at P < 10 kbar'
         stop

      end if

      tc = tk - 273.15d0

      gh2o(1) = 0d0
      gh2o(2) = 0d0

      do j = 1, 5
         do l = 1, 6 - j 
           gh2o(1) = gh2o(1)+(b(j,l)*10**c(j,l))*(tc**(j-1))*(p**(l-1))
           gh2o(2) = gh2o(2)+(b(j,l)*10**c(j,l))*(tc**(j-1))*(p0**(l-1))
         end do 
      end do 

      fh2o = ((gh2o(1)-gh2o(2))/(1.9872d0*tk)) 

99    end 



      subroutine saxfei
c---------------------------------------------------------------------
c       subroutine to calculate compressibilities and fugacity coefs
c       according to saxena & fei (1987a, b, 1988) for species in the
c       c-o-h system
 
c       literature: saxena s.k. & fei y. (1987a):
c                   fluids at crustal pressures and temperatures. contr
c                   mineral. petrol., vol. 95, pp 370 - 375.
c
c                   saxena s.k. & fei y. (1987b):
c                   high pressure and hig temperature fluid fugacities.
c                   geochim. cosmochim,. acta, vol. 51, pp 783 - 791.
c
c                   saxena s.k. & fei y. (1988):
c                   the pressure-volume-temperature equation of hydrogen
c                   geochim. cosmochim. acta, vol. 52, pp 1195 - 1196.
c---------------------------------------------------------------------
c       written by peter ulmer, geophysical lab, ciw, june, 1988
c       modified for vertex/frendly by j. connolly, may, 1989.
c       never tested.
c---------------------------------------------------------------------
c          z = compressibility = pv/rt
c             species are: 1 = co2
c                          2 = co
c                          3 = ch4
c                          4 = o2
c                          5 = h2o
c                          6 = h2
c---------------------------------------------------------------------
      implicit none
 
      double precision tcr(6),pcr(6),fcoeff(6),a(3),b(2),c(3),d(2),
     *                 aw(4),bw(3),cw(4),dw(4),at(6),tr,
     *                 bt(6),ct(6),dt(6),pr1(6),pr2(6)
 
      double precision p,t,xco2,u1,u2,trr,prr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,trr,prr,r,ps

      double precision f
      common/ cst11 /f(3)
 
      save tcr, pcr, a, b, c, d, aw, bw, cw, dw
 
      data tcr / 304.2d0, 133.15d0, 190.7d0, 154.8d0, 647.3d0, 33.1d0/
      data pcr / 73.9d0, 35d0, 46.4d0, 50.8d0, 221.3d0, 13d0/
 
      data a / 2.0614, -2.2351, -0.39411/
      data b / 5.5125d-2, 3.9344d-2/
      data c / -1.8935d-6, -1.1092d-5, -2.1892d-5/
      data d / 5.0527d-11, -6.3033d-21/
 
      data aw / 1.4937d0, -1.8626d0, 0.80003d0, -0.3941d0/
      data bw / 4.2410d-2, 2.4097d-2, -8.9634d-3/
      data cw / -9.016d-7, -6.1345d-5, 2.2380d-5, 5.2335d-7/
      data dw / -7.6707d-9, 4.1108d-8, -1.4798d-8, -6.3033d-21/
c---------------------------------------------------------------------
 
c       calculation of the a, b, c, d terms for the different species
 
c       calculation for co2, co, ch4, and o2  ** change to index i
c       the indices on lhs of equalities.
 
        if (p.lt.5d3) then
           write (*,*) 'this version of saxfei is invalid at ',
     *                 'pressures less than 5000 bars.'
           stop
        end if
 
           tr = t/tcr(1)
           at(1) = a(1) + a(2)*tr**(-2) + a(3)*dlog(tr)
           bt(1) = b(1)/tr + b(2)*tr**(-2)
           ct(1) = c(1)/tr + c(2)*tr**(-2) + c(3)*tr**(-3)
           dt(1) = d(1)/tr + d(2)*tr**3

c---------------------------------------------------------------------
c       calculation of a, b, c, d for h2o
 
        tr = t/tcr(5)
        at(5) = aw(1) + aw(2)*tr**(-2) + aw(3)*tr**(-3) + aw(4)*dlog(tr)
        bt(5) = bw(1)/tr + bw(2)*tr**(-2) + bw(3)*tr**(-3)
        ct(5) = cw(1)/tr + cw(2)*tr**(-2) + cw(3)*tr**(-3) +
     .          cw(4)*dlog(tr)
        dt(5) = dw(1)/tr + dw(2)*tr**(-2) + dw(3)*tr**(-3) +
     .          dw(4)*tr**3
 
c       calculation of the fugacity coefficients according to
c       saxena & fei (1987b) for co2
 
        pr1(1) = p/pcr(1)
        pr2(1) = 5d3/pcr(1)
        fcoeff(1) = r*t * (at(1)*dlog(pr1(1)/pr2(1)) +
     .              bt(1)*(pr1(1)-pr2(1)) +
     .    ct(1)*(pr1(1)**2 - pr2(1)**2)/2d0 + dt(1)*
     .    (pr1(1)**3 - pr2(1)**3)/3) - r*t*dlog(p)
     .    + 1d3*(-40.468d0 + 0.06702d0*t + 8.3481d0*dlog(t))
        fcoeff(1) = fcoeff(1)/(r*t)
        fcoeff(1) = dexp(fcoeff(1))
 
      if (xco2.ne.0d0) then
         f(2) = dlog(xco2 * fcoeff(1) * p)
      else
         f(2) = 0d0
      end if
 
c      calculation of the fugacity coefficients according to
c      saxena & fei (1987b) for h2o
 
      pr1(5) = p/pcr(5)
      pr2(5) = 5d3/pcr(5)
      fcoeff(5) = r*t * (at(5)*log(pr1(5)/pr2(5)) +
     .              bt(5)*(pr1(5)-pr2(5)) +
     .    ct(5)*(pr1(5)**2 - pr2(5)**2)/2d0 + dt(5)*
     .    (pr1(5)**3 - pr2(5)**3)/3d0) - r*t*dlog(p)
     .    + 1d3*(-130.517d0 + 0.06497d0*t + 19.4831d0*dlog(t))
      fcoeff(5) = fcoeff(5)/(r*t)
      fcoeff(5) = dexp(fcoeff(5))
 
      if (xco2.ne.1d0) then
         f(1) = dlog( (1d0-xco2) * fcoeff(5) * p)
      else
         f(1) = 0d0
      end if
 
      end

