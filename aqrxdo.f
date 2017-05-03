      subroutine slvnt1 (gsolv)
c-----------------------------------------------------------------------
c computes solvent p-t-composition dependent properties: permittivity (eps), 
c molar mass, debye-hueckel (adh), and gHKF function.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision gsolv, vsolv, cdh, xt, xp, ysum, ysolv(nsp)  

      double precision gfunc, ghybrid, gcpd

      external gfunc, ghybrid, gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision yf,g,v,eps,vf,eps0
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),eps(nsp),vf(nsp),eps0(nsp)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy

      double precision vol
      common/ cst26 /vol

      save cdh 
      data cdh/-42182668.74d0/
c----------------------------------------------------------------------
      gsolv  = 0d0 
      msol   = 0d0
      ysum   = 0d0 

      do i = 1, ns
c                                 solvent mass, kg/mol compound
         msol   = msol + y(i) * fwt(jnd(i))
c                                 g mech mix term for solvent:
         gsolv  = gsolv + aqg(i) * y(i)

         ysum = ysum + y(i)

      end do
c                                 compute normalized fractions and correct solvent 
c                                 gibbs energy for ideal solute concentration
      do i = 1, ns
 
         ysolv(i) = y(i)/ysum
         
      end do
c                                 compute and add in solvent activities
      gsolv = gsolv + ysum*( ghybrid (ysolv) + rt*dlog(ysum) )
c                                 the call to ghybrid loads the partial 
c                                 molar volumes (cstcoh) and total volume (cst26)
c                                 from routine mrkmix
      vsolv = ysum * vol
c                                 get permittivity based on mrk volumetric properties
      call geteps (epsln)
c                                 set p/t to reference condition to get reference 
c                                 permitivity 
      xp = p
      xt = t

      p = pr
      t = tr 

      call mrkmix (ins, ns, 1)

      call geteps (epsln0)

      p = xp
      t = xt                              
c                                 Debye-Hueckel factor, A[cgs] = -q^3*sqrt(NA)/(4*Pi*k^(3/2))
c                                 *(msolg/(10*vh2o))^(1/2)/(epsilon*T)^(3/2) for ln(gamma) = +A*....
c                                 A = cdh*(msolkg/(vh2ojbar))^(1/2)/(epsilon*T)^(3/2)
c                                 this, like epslon, could be improved by using non-ideal 
c                                 volumes
      adh = cdh * dsqrt(1d1*msol/vsolv/(epsln*t)**3)
c                                 shock et al 1992 g function (cgs solvent density),
c                                 used by hkf
      gf = gfunc (msol*1d3/vsolv) 

      end 

      subroutine slvnt2 (gsolv)
c-----------------------------------------------------------------------
c computes solvent p-t-composition dependent properties: permittivity (eps), 
c molar mass (msol, this could be stored in thermo), debye-hueckel (adh).

c assumes pure species volumes in cohhyb have been initialized. and that
c the molar speciation is stored in y.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer l, k

      double precision gsolv, mo(m4), lng0

      double precision gcpd

      external gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy
c---------------------------------------------------------------------- 
      lng0 = 0d0 
c                                 molalities and ionic strength
      do k = sn1, nqs
c                                 ln molality of solutes
         mo(k) = y(k)/msol
         lng0 = lng0 + q2(k) * mo(k)

      end do 
c                                 at this point lng0 is ionic strength
      lng0 = lng0/2d0 
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
c                                 Davies extension.
      lng0 = adh*dsqrt(lng0)/(1d0 + dsqrt(lng0)) + 0.2d0*lng0
c                                 add in the solute gibbs energies
c                                 neutral solutes, ideal
      do l = sn1, sn

         if (mo(l).le.0d0) cycle
         gsolv = gsolv + y(l) * (gcpd(jnd(l),.true.) + rt*dlog(mo(l)))

      end do
c                                 ionic solutes, Davies D-H extension
      do k = l, nqs

         if (y(k).le.0d0) cycle
         gsolv = gsolv + y(k) * (gcpd(jnd(k),.true.)
     *                        + rt*(dlog(mo(k))+ lng0*q2(k)))

      end do 

      end

      subroutine solut0 (id)
c-----------------------------------------------------------------------
c computes solute endmember g's for aqueous solutions, used only for 
c output purposes by routine calpr0.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id

      double precision gunk, ysum

      double precision gfunc, gcpd

      external gfunc, gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision yf,g,v,eps,v0,eps0
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),eps(nsp),v0(nsp),eps0(nsp)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(m4,k5),spct(h9),spnams(m4,h9)
c----------------------------------------------------------------------
      rt = r*t 
      ysum = 0d0 

      do i = 1, ns
 
         ysum = ysum + ysp(i,id)
c                                 solvent endmember gibbs energy
         aqg(i) = gcpd(jnd(i),.true.)

      end do
c                                 call mrkmix to get solvent volumetric
      do i = 1, ns
c                                 load normalized molecular fluid composition
         yf(ins(i)) = ysp(i,id)/ysum

      end do 
         
      call mrkmix (ins,isp,1)
c                                 solvent properties
      call slvnt1 (gunk)
c                                 solute gibbs energies
      do i = sn1, nqs

         aqg(i) = gcpd(jnd(i),.true.)

      end do 

      end

      subroutine geteps (epsln)
c-----------------------------------------------------------------------
c computes dielectric constant for the molecular species of a solvent
c assumes species volumes or partial molar volumes in cohhyb have been 
c initialized.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j

      double precision po(nsp,8), trt, rho, epsln, eps

      double precision epsh2o

      external epsh2o

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision y,g,v,eps1,vf,eps0
      common/ cstcoh /y(nsp),g(nsp),v(nsp),eps1(nsp),vf(nsp),eps0(nsp)
c                                Harvey & Lemmon provide additional data for 
c                                ethylene and long-chain hydrocarbons. A_mu
c                                is zero for all species listed here, therefore
c                                Eq 5 of H&L:

c                                P/rho (cm3/mol) = A + A_mu/T + B*rho + C*rho^D

c                                is simplified by dropping the second term the 
c                                coefficients are a f(T) viz A = a0 + a1*(T/Tr - 1)...

c                                po(i,1:8) - a0, a1, A_mu, b0, b1, c0, c1, D
      data ((po(i,j),j=1,8),i=1,nsp)/ 
c                                1 - H2O
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 
c                                2 - CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                3 - CO approximated by O2
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0, 
     *     1.5d0,
c                                4 - CH4
     *     6.5443d0, 1.33d-2, 0d0,8.4578d0, 3.7196d0, -352.97d0, 
     *     -100.65d0, 2d0,
c                                5 - H2
     *     2.0306d0, 5.6d-3, 0d0, 0.181d0, 0.021d0, -7.4d0, 0d0, 2d0,
c                                6 - H2S approximate by ?
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 
c                                7 - O2
     *     3.9578d0, 6.5d-3, 0d0, 0.575d0, 1.028d0, -8.96d0, -5.15d0, 
     *     1.5d0,
c                                8 - SO2 approximated by CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                9 - COS approximated by CO2
     *     7.3455d0, 3.35d-3, 0d0, 83.93d0, 145.1d0, -578.8d0, -1012d0, 
     *     1.55d0,
c                                10 - N2
     *     4.3872d0, 2.26d-3, 0d0, 2.206d0, 1.135d0, -169d0, -35.83d0, 
     *     2.1d0,
c                                11 - NH3 approximated by ?
     *     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,
c                                12-15 Si-O high T species
     *     32*0d0,
c                                16 - Ethane
     *     11.1552d0, 0.0112d0, 0d0, 36.759d0, 23.639d0, -808.03d0, 
     *     -378.84d0, 1.75d0,
c                                17 - dilutant
     *      8*0d0/

      save po 
c----------------------------------------------------------------------
      trt = t/273.16d0 - 1d0

      epsln = 0d0 

      do i = 1, ns - 1

         j = ins(i)
c                                 rho = 1/vcm3, v(j) is initialized by lnfpur
c                                 in gcpd and is in cm3/mol
         rho = 1d0/v(j)
c                                 Eq 5 of H&L 2005 for  polarization/rho
c                                 for polar species need to add po(j,3)
         eps = po(j,1) + po(j,2)*trt + (po(j,4) + po(j,5)*trt)*rho
     *                + (po(j,6) + po(j,7)*trt)*rho**po(j,8)

         if (po(j,3).eq.0d0) then 
c                                 invert clausius-mosotti relation for dielectric
c                                 constant (Eq 1) non-polar molecules
            eps = (2d0*eps*rho + 1d0) / (1d0 - rho*eps)
         else
c                                 polarized species, currently none, but H2S...
            eps = po(j,3) + (9d0*eps*rho + 1d0 
     *          + 3d0*dsqrt((3d0*eps*rho)**2 + 2d0*eps*rho + 1d0))/4d0

         end if
c                                  Looyenga mixing rule justified by Mountain &
c                                  Harvey 2015, modified here to use volume fraction
c                                  computed from partial molar volumes.
         epsln  = epsln  + vf(j) * eps**(1d0/3d0)

      end do 
c                                  add in water (ns^th species)
      epsln  = epsln + vf(ins(i)) * epsh2o (v(ins(i))/1d1)**(1d0/3d0)

      epsln  = epsln **3

      end 