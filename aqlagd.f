      subroutine aqlagd (id,bad)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c configured to be called from resub with output to the (molar normalized)
c arrays g2/cp2
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, it, iexp, id

      logical bad

      double precision c(l9), mo(l9), dg(l9), xis, blk(k5), dn,
     *                 d(l9), lng0, is, gamm0, totm, g0(l9), 
     *                 gso(nsp), lnkw,
     *                 gtot, smo, xdn, slvmo(nsp)

      double precision gcpd, solve, gfunc

      external gcpd, solve, gfunc

      integer ion, ichg, jchg
      double precision q, q2, qr, dcp
      common/ cstaq /q(l9),q2(l9),qr(l9),dcp(k5,l9),jchg(l9),ichg,ion

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

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

      double precision atwt
      common/ cst45 /atwt(k0) 

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision yf,g,v,vf
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),vf(nsp)

      double precision gh,vh,gp
      common/ csthyb /gh(nsp),vh(nsp),gp(nsp) 

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

      integer kd
      double precision x3, caq, ionst, tmolal
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),ionst(k5),tmolal(k5),kd

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c----------------------------------------------------------------------
      if (.not.mus) then 

         bad = .true.
 
         return

      else 

         bad = .false.

      end if 

      is   = 0d0

      call slvnt3 (gso)

      g0(ion) = gcpd(aqst+ion,.false.)
c                                 compute solute properties 
      do i = 1, aqct 

         k = aqst + i
c                                 dg is the solvent oxide potentials - g
         g0(i) = gcpd(k,.false.)
         dg(i) = -g0(i) + qr(i)*g0(ion)

         do j = 1, jbulk 
c                                 if oxide components, but no excess oxygen
c                                 mu(O2) is a nan. 
            if (isnan(mu(j))) cycle

            dg(i) = dg(i) + dcp(j,i) * mu(j)

         end do 
c                                 normalize by RT
         dg(i) = dg(i)/rt

         if (q(i).ne.0d0) then 
c                                  this is now c(i)*a(ion)^(q(i)) = mo(i)*gamma(i)*q(i)
c                                  the rhs q(i) is because the eq to be solved will be
c                                  sum (q(i)*m(i)) = 0. 
            d(i) = q(i)*dexp(dg(i))
            c(i) = d(i)

         else 
c                                  neutral species assumed to be ideal, molality is
            mo(i) = dexp(dg(i))
      
         end if 

      end do 

      gamm0 = 1d0 
      it = 0
      lnkw = (gso(ns)-g0(ioh))/rt
      mo(ion) = dexp(lnkw/2d0)
      xdn = 1d0
      iexp = 0
c                                  iterative loop for ionic strength, this is the achilles 
c                                  heel cause the solute g's are based only on the previous
c                                  chemical potentials.
      do 
c                                  solve charge balance for H+
         mo(ion) = solve(c,qr,mo(ion),jchg,ichg,bad)
c                                  back calculate charged species molalities
c                                  and ionic strength
         xis = is
         is = 0d0

         do i = 1, ichg 
            j = jchg(i)
            mo(j) = c(j) * mo(ion)**qr(j) / q(j)
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

         if (bad) exit
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
            c(j) = d(j)*(gamm0*q2(ion))**qr(j)/(gamm0*q2(j))
         end do

      end do
c                                 back calculated bulk composition
      if (bad) then 

         call warn (99,0d0,0,'AQLAGD did not converge on solute '//
     *                       'speciation')
 
         return

      end if 

      do j = 1, icp
         blk(j) = 0d0
      end do  

      smo = 0d0 
      gtot = 0d0
c                                 everything on a molal basis
c                                 first the solutes
      do i = 1, aqct
c                                 total g
         if (q2(i).ne.0d0) then 
            gtot = gtot + mo(i) * (g0(i) + rt*dlog(mo(i)*gamm0*q2(i)))
         else 
            gtot = gtot + mo(i) * (g0(i) + rt*dlog(mo(i)))
         end if 
c                                  total molality
         smo = smo + mo(i)
c                                 accumulate component moles
         do j = 1, icp
            blk(j) = blk(j) + mo(i)*aqcp(j,i)
         end do 

      end do
c                                 for the solvent mole fractions
c                                 need to accumulate total
c                                 molality first
      do i = 1, ns 
c                                 solvent molality is
         slvmo(i) = yf(ins(i))/msol
c                                 total molality
         smo = smo + slvmo(i)
c                                 moles/kg-solvent 
         do j = 1, icp 
            blk(j) = blk(j) + slvmo(i)*cp(j,jnd(i))
         end do

      end do

      do i = 1, ns
         caq(id,i) = slvmo(i)/smo
         if (caq(id,i).eq.0d0) cycle
         gtot = gtot + slvmo(i) * (gso(i) + rt*dlog(caq(id,i)))
      end do 
c                                 bulk fluid composition 
      totm = 0d0 

      do j = 1, icp
         totm = totm + blk(j)
      end do
c                                 load into molar normalized arrays 
c                                 used by resub
      g2(jphct) = gtot/totm

      do j = 1, icp
         cp2(j,jphct) = blk(j)/totm
      end do

      do i = 1, aqct 
         caq(id,ns+i) = mo(i)
      end do 

      tmolal(id) = smo
      ionst(id) = is 

      end

      subroutine slvnt3 (gso)
c-----------------------------------------------------------------------
c for lagged speciation calculations get solvent properties
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, k

      double precision gso(nsp), dum

      double precision gcpd

      external gcpd

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision yf,g,v,vf
      common/ cstcoh /yf(nsp),g(nsp),v(nsp),vf(nsp)

      double precision gh,vh,gp
      common/ csthyb /gh(nsp),vh(nsp),gp(nsp) 

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c----------------------------------------------------------------------
      rt   = r*t
c                                 doing lagged aqueous speciation
      if (ns.gt.1) then
c                                 the incoming pure solvent fractions from 
c                                 resub are in y(spec) cxt7 and load
c                                 into the fluid species indexed array yf
         do i = 1, ns
            yf(ins(i)) = y(i)
         end do 
c                                 a multi species solvent is present: 
c                                 get solvent dielectric constants,
c                                 call mrkmix to get solvent volumetric props

c                                 for now, don't know if the pure species have
c                                 been computed (this will be the case in production)
         do k = 1, ns
            gso(k) = gcpd(jnd(k),.false.)
         end do
c                                  next get the activity coefficient ratios for the
c                                  hybrid eos (csthyb g(i) = phi0(EoS)/phi0(MRK)
         call hybeos (ins,ns)  
c                                  now get the hybrid fugacity coefficients 
c                                  calculate hybrid fugacity coefficients
         call mrkmix (ins, ns, 1)
c
         do k = 1, ns
c                                  now add in the pure solvent fugacities
c                                  at this point, assuming the fugacity coeff 
c                                  is independent of speciation then the partial
c                                  g of a solvent species will be 
c                                  g(i) = gs0(i) + RT ln x(i)
c                                  where sumsol is the sum of the solute
c                                  mole fractions. 
            i = ins(k)
            gso(k) = gso(k) + rt * dlog(g(i)/gp(i))

         end do 
c                                  ysum is just a dummy.
         call slvnt1 (dum)

      else 
c                                 solvent is pure water 
         call slvnt0 (gso(1),dum)

         yf(ins(1)) = y(1)

      end if 

      end 


      subroutine aqlgot (jd)
c-----------------------------------------------------------------------
c output lagged aqueous speciation for 2d tab output by werami.  
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, jd

      double precision tmass

      double precision gcpd

      external gcpd

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      double precision atwt
      common/ cst45 /atwt(k0) 

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer kd, na1, na2, na3, na4
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,na4,kd
c----------------------------------------------------------------------
      if (jd.eq.0) then 

         do i = 1, iprop
            prop(i) = nopt(7)
         end do 
 
         return

      end if 

      if (lopt(23)) then
c                                 molar bulk composition
         do i = 1, icp
            prop(i) = cp3(i,jd)*1d2
         end do 

      else 
c                                 mass bulk composition
         tmass = 0d0

         do i = 1, icp
            tmass = tmass + atwt(i)*cp3(i,jd)
         end do 

         do i = 1, icp
            prop(i) = cp3(i,jd)*atwt(i)/tmass*1d2
         end do 

      end if 

      k = icp

      do i = 1, nsa 

         k = k + 1

         if (i.lt.sn1) then 
c                                 solvent speciation
            if (lopt(26)) then
c                                 mole fraction
               prop(k) = caq(jd,i)
            else 
c                                 molality
               prop(k) = caq(jd,i)/msol
            end if  

         else 
c                                 solute speciation
            if (lopt(27)) then
c                                 molality
               prop(k) = caq(jd,i)
            else 
c                                 mole fraction
               prop(k) = caq(jd,i)/totm
            end if 

         end if  

      end do 

      call slvnt3 (gso)
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
      gamm0 = dexp(adh*dsqrt(caq(jd,na1))/(1d0 + dsqrt(caq(jd,na1))) 
     *             + 0.2d0*caq(jd,na1))
c                                 other properties:
c                                 pH
      prop(k+2) = -dlog10(caq(jd,ns+ihy)*gamm0)
c                                 pH - pH0
      prop(k+1) = prop(k+2) + (gso(ns)-gcpd(aqst+ioh,.false.))/rt/2d0/2.302585d0
c                                 err_log10(K_w)
      prop(k+3) = dabs(caq(jd,))
c                                 dielectric constant
      prop(k+4) = epsln
c                                 ionic strength
      prop(k+5) = caq(jd,na1)
c                                 solute molality?
      prop(k+6) = caq(jd,na2)

      end 