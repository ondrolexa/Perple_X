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

      double precision yf,g,v,vf
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),vf(nsp)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

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

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh
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

      double precision yf,g,v,vf
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),vf(nsp)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

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

      double precision y,g,v,vf
      common/ cstcoh /y(nsp),g(nsp),v(nsp),vf(nsp)
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

      subroutine aqrxdo (jd,lu)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, l, ichg, jchg(l9), it, ind(l9), id, lu, iexp

      logical bad, output

      character text*200

      double precision c(l9), q(l9), mo(l9), dg(l9), ahy, xis, blk(k5),
     *                 d(l9), q2(l9), lng0, is, gamm0, totm, g0(l9),
     *                 gso(nsp), ysum, ph0, v0(nsp), vf0(nsp), tmass,
     *                 tsmas, tsmol, smol(k5), dn, xdn, errkw

      double precision gcpd, solve, gfunc

      external gcpd, solve, gfunc

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      character cname*5
      common/ csta4  /cname(k5)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision fwt
      common/ cst338 /fwt(k10)

      logical hsccon
      double precision atwt, sel
      common/ cst45 /atwt(k0), sel(k0), hsccon

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision yf,g,v,vf
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),vf(nsp)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      character names*8
      common/ cst8  /names(k1)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq
c----------------------------------------------------------------------
      if (.not.mus) then 
         write (lu,*) ' no mus, cannot back calculate speciation'
         return
      end if 

      output = .true. 

      ichg = 0 
      is   = 0d0
      rt   = r*t
      ysum = 0d0 
c                                 get the solvent permittivities,
c                                 call mrkmix to get solvent volumetric props
      do i = 1, ns
 
         ysum = ysum + x(1,i)

      end do

      do i = 1, ns
c                                 load normalized molecular fluid composition
         yf(ins(i)) = x(1,i)/ysum

      end do 
         
      call mrkmix (ins,isp,1)
c                                  save pmv's and v fractions for output
      do i = 1, ns
         vf0(ins(i)) = vf(ins(i))
         v0(ins(i)) = v(ins(i))
      end do 
c                                  ysum is just a dummy at this point. 
      call slvnt1 (ysum)

      do i = 1, ns 

         gso(i) = 0d0

         do j = 1, icp

            if (isnan(mu(j))) cycle

            gso(i) = gso(i) + mu(j)*cp(j,jnd(i))

         end do 

      end do 
c                                 compute solute properties 
      do i = 1, aqct 

         if (i.eq.ihy) cycle

         k = aqst + i 

         q(i) = thermo(6,k)
         q2(i) = q(i)**2
c                                 dg is the solvent oxide potentials - g
         g0(i) = gcpd(k,.false.)
         dg(i) = -g0(i)

         do j = 1, jbulk 
c                                 if oxide components, but no excess oxygen
c                                 mu(O2) is a nan. 
            if (isnan(mu(j))) cycle

            dg(i) = dg(i) + (aqcp(j,i) - q(i)*aqcp(j,ihy)) * mu(j)

         end do 
c                                 normalize by RT
         dg(i) = dg(i)/rt

         if (q(i).ne.0d0) then 

            ichg = ichg + 1
            jchg(ichg) = i 
c                                  this is now c(i)*a(H+)^(q(i)) = mo(i)*gamma(i)*q(i)
            d(i) = q(i)*dexp(dg(i))
            c(i) = d(i)

         else 
c                                  neutral species assumed to be ideal, molality is
            mo(i) = dexp(dg(i))
      
         end if 

      end do 

      gamm0 = 1d0 
      it = 0 
      mo(ihy) = 1e-6
      xdn = 1d0
      iexp = 0
c                                  iterative loop for ionic strength
      do 
c                                  solve charge balance for H+
         mo(ihy) = solve(c,q,mo(ihy),jchg,ichg,bad)

         ahy = mo(ihy)*gamm0
c                                  back calculate charged species molalities
c                                  and ionic strength
         xis = is
         is = 0d0

         do i = 1, ichg 
            j = jchg(i)
            mo(j) = ahy**(q(j))*c(j)/q(j)
            if (mo(j).gt.10d0) then
c               write (lu,*) aqnam(j), mo(j), q(j)
            end if 
            is = is + q2(j) * mo(j)
         end do

         is = is / 2d0 

         dn = is - xis

         if (dabs(dn).gt.1d0/2d0**iexp) then 
            dn = dn/dabs(dn)/2d0**iexp 
            if (dn*xdn.lt.0d0) iexp = iexp + 1
         end if 

         is = xis + dn

         xdn = dn 

         if (bad) then 
            write (lu,*) 'bombed'
            exit
         end if 
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
         lng0 = adh*dsqrt(is)/(1d0 + dsqrt(is)) + 0.2d0*is
         gamm0 = dexp(lng0)
c                                 check for convergence
         if (dabs(xis-is)/is.lt.nopt(5)) then 
            exit
         else if (it.gt.iopt(21)) then 
            bad = .true.
         end if 

         it = it + 1
c                                 update coefficients
         do i = 1, ichg 
            j = jchg(i)
            c(j) = d(j)*dexp(lng0*(1d0-q2(j)))
         end do

      end do

      if (output) then
c                                 neutral pH
         ph0 = (g0(ihy)+g0(ioh)-gso(ns))/2d0/rt/2.302585d0

         do i = 1, icp
            blk(i) = 0d0
         end do  
c                                 compute mole fractions, total moles first
         do i = 1, ns 
c                                 moles/kg-solvent 
            do j = 1, icp 
               blk(j) = blk(j) + x(1,i)*cp(j,jnd(i))/msol
            end do 
         end do

         do i = 1, aqct
            ind(i) = i 
            do j = 1, icp
               blk(j) = blk(j) + mo(i)*aqcp(j,i)
            end do 
         end do

         totm = 0d0

         do i = 1, icp 
            totm = totm + blk(i)
         end do 

         call rankem (mo,ind,aqct,iopt(32))

         write (lu,1000)

         errkw = ((g0(ioh)-gso(ns))/rt + dlog(mo(ihy)*mo(ioh)*gamm0**2
     *             ))/2.302585d0

         write (text,1050) -dlog10(ahy),dabs(errkw),ph0,is,gamm0
         call deblnk (text)
         write (lu,'(400a)') (chars(j), j = 1, length)
         write (text,1070) epsln,epsln0,msol*1d3
         call deblnk (text)
         write (lu,'(400a)') (chars(j), j = 1, length)

         if (jdaq.eq.20) then

            write (lu,1100)

         else

            write (lu,1040)

         end if 

         do i = 1, iopt(32)

            k = ind(i)

            if (jdaq.eq.20) then 
c                                 check if the species is the solution model
               l = 0

               do j = sn1, nqs
                  if (jnd(j)-aqst.eq.k) then
                     l = j
                     exit
                  end if 
               end do

               if (l.ne.0) then 
                  write (lu,1080) aqnam(k),int(thermo(6,k+aqst)),
     *                            mo(k),mo(k)/totm,x(1,l),
     *                           int(g0(k)+rt*(dlog(mo(k))+lng0*q2(k))),
     *                           int(g0(k))
               else
                  write (lu,1090) aqnam(k),int(thermo(6,k+aqst)),
     *                            mo(k),mo(k)/totm,
     *                           int(g0(k)+rt*(dlog(mo(k))+lng0*q2(k))),
     *                           int(g0(k))
               end if 

            else 

               write (lu,1010) aqnam(k),int(thermo(6,k+aqst)),
     *                         mo(k),mo(k)/totm,
     *                         int(g0(k)+rt*(dlog(mo(k))+lng0*q2(k))),
     *                         int(g0(k))

            end if 
         end do 

         write (lu,1020)

         do i = 1, ns 

            write (lu,1030) names(jnd(i)), x(1,i)/msol, x(1,i),
     *                      vf0(ins(i)), v0(ins(i)), 
     *                      int(gso(i)), int(gcpd(jnd(i),.true.))
         end do 

         write (lu,1060)
c                                 bulk fluid composition 
         tmass = 0d0
         totm = 0d0 

         do i = 1, icp
            totm = totm + blk(i)
            tmass = tmass + atwt(i)*blk(i)
         end do 

         if (jdaq.eq.20) then 

            tsmol = 0d0
            tsmas = 0d0
 
            do i = 1, icp
               smol(i) = 0d0
            end do 

            do i = 1, nqs
               
               k = jnd(i)

               do j = 1, icp

                  if (i.lt.sn1) then 
                     dn = x(1,i) * cp(j,k)
                  else 
                     dn = x(1,i) * aqcp(j,k-aqst)
                  end if 

                  smol(j) = smol(j) + dn

                  tsmol = tsmol + dn
                  tsmas = tsmas + dn*atwt(j)

               end do

            end do 

            write (lu,1130)

            do i = 1, icp
               write (lu,1110) cname(i),blk(i)/totm*1d2,
     *                                  smol(i)/tsmol*1d2,
     *                                  blk(i)*atwt(i)/tmass*1d2,
     *                                  smol(i)*atwt(i)/tsmas*1d2
            end do

         else 

            write (lu,1120)

            do i = 1, icp
               write (lu,1110) cname(i),blk(i)/totm*1d2,
     *                                  blk(i)*atwt(i)/tmass*1d2
            end do

         end if 

      end if

1000  format (/,'Back-calculated solute speciation in the solvent:',/)
1010  format (a8,4x,i2,3x,g12.6,3x,g12.6,5x,i8,5(2x,g12.6))
1020  format (/,'Solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'vol_fraction',
     *        2x,'v,cm3/mol*',2x,'g,J/mol*',3x,'g0,J/mol***')
1030  format (a8,2x,f7.4,5x,f7.5,8x,f7.5,5x,f7.3,3x,i8,4x,i8)
1040  format (/,'Solute endmember properties:',//,10x,'charge',3x,
     *       'molality',5x,'mol_fraction',6x,'g,J/mol*',5x,'g0,J/mol**')
1050  format ('pH = ',f6.3,'+/-',f5.3,
     *        '; neutral_pH = ',f6.3,'; ionic_strength = ',
     *        g10.4,'; gamma/q^2 = ', g10.4)
1060  format (/,'*partial molar, **molal ref. state, ***molar ref. ',
     *        'state.',/)
1070  format ('permittivity = ', g10.4,
     *        '; reference state permittivity = ',g10.4
     *        '; solvent molar mass, g/mol = ',f8.4)
1080  format (a8,4x,i2,3x,g12.6,3x,g12.6,3x,g12.6,5x,i8,5(2x,g12.6))
1090  format (a8,4x,i2,3x,g12.6,3x,g12.6,20x,i8,5(2x,g12.6))
1100  format (/,'Solute endmember properties:',/,
     *        48x,'optimized',/,10x,'charge',3x,
     *       'molality',5x,'mol_fraction',3x,'mol_fraction',6x,
     *       'g,J/mol*',5x,'g0,J/mol**')
1110  format (1x,a8,2x,4(g12.6,3x))
1120  format (/,'Back-calculated fluid bulk composition:',//,
     *        13x,'mol %',11x,'wt %')
1130  format (/,'Back-calculated vs optimized fluid bulk composition:',
     *        //,26x,'optimized',21x,'optimized',/,
     *        13x,'mol %',10x,'mol %',10x,'wt %',11x,'wt %')
c            write (*,1030) names(jnd(i)),x(1,i)/msol,x(1,i),int(gso(i)),
c     *                      v(ins(i)),eps(i),v0(ins(i)),eps0(i)
c1030  format (a8,2x,f7.4,5x,f7.5,5x,i8,2x,f9.4,4x,f8.5,5x,f9.4,5x,f8.5)


      end 

      subroutine rankem (a,ind,r,n)
c-----------------------------------------------------------------------
c rank the n largest values of array a(left:right) in array ind. assumes ind has 
c been initialized (left:right). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, imax, ind(*), r, left, right, n

      double precision a(*), amax
c----------------------------------------------------------------------
      left = 1
      right = r

      do 

         amax = -1d99

         do i = left, right 
            if (a(ind(i)).gt.amax) then
               imax = i
               amax = a(ind(i))
            end if 
         end do 

         j = ind(left)
         ind(left) = ind(imax)
         ind(imax) = j 

         left = left + 1

         if (left.eq.n) return

      end do 

      end 


      double precision function solve (c,q,x,jchg,ichg,bad)
c-----------------------------------------------------------------------
c function solve for hydronium molality (x) from charge balance for aqrxdo.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer jchg(*), ichg, it, i, j

      logical bad

      double precision c(*), q(*), x, y, z, f, df, dx

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------

      it = 0 

      do

         f = x
         df = 1d0
         it = it + 1

         do i = 1, ichg 

            j = jchg(i)

            y = x**(q(j)) * c(j)
            z = y*q(j)/x 

            f = f + y
            df = df + z

         end do  

         dx = -f/df

         x = x + dx

         if (dx/x.lt.nopt(5)) then 
            bad = .false.
            exit
         else if (x.lt.0d0.or.it.gt.iopt(21)) then 
            bad = .true.
            exit
         end if 

      end do

      solve = x

      end 

      subroutine setxyp (ids,id)
c-----------------------------------------------------------------------
c for static pseudocompounds load the compositional coordinates from xco
c into simple compositional arrays for compound id of solution ids. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ids, id, i

      double precision xt
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)
c-----------------------------------------------------------------------
c                                 y's are primsatic 1-d coordinates
c                                 x's are prismatic 2-d coordinates
c                                 p0a's are independent 1-d endmember coordinates 
c                                       of a fully disordered o/d model
c                                 pa's are, ultimately, independent 1-d endmember 
c                                       coordinates of the stably orderded configuration.

c                                 y's are kept distinct form p0's because in the 
c                                 adaptive phase of minimization the all compositions
c                                 are stored in the x version. these are then converted
c                                 to y's and o/d models subsequently convert
c                                 these to p0a's. it is not clear whether the original
c                                 y's are subsequently referenced

c                                 currently the only models that use the x form are 
c                                 nastia's alloy special cases ksmod 29-31. these should be
c                                 checked.
      xt = 0d0 

      if (lorder(ids).and.lrecip(ids)) then 

         do i = 1, nstot(ids) - 1
            p0a(i) = xco(jco(id)+i) 
            pa(i) = p0a(i)
            xt = xt + pa(i)
         end do 

         p0a(i) = 1d0 - xt  
         pa(i) = p0a(i)
  
      else if (lorder(ids)) then

         do i = 1, lstot(ids) - 1
            p0a(i) = xco(jco(id)+i) 
            pa(i) = p0a(i)
            xt = xt + pa(i)
         end do

         p0a(i) = 1d0 - xt
         pa(i) = p0a(i)

      else 

         do i = 1, lstot(ids)- 1
            y(i) = xco(jco(id)+i) 
            xt = xt + y(i)
         end do

         y(i) = 1d0 - xt

       end if

       end