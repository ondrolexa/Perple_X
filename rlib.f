c-----------------------------------------------------------------------

c RLIB - a library of subprograms called by VERTEX, FRENDLY, ACTCOR
c BUILD, ADIABAT, and ISOCHOR.

c Copyright (c) 1990 by James A. D. Connolly, Institute for Mineralogy
c & Petrography, Swiss Federal Insitute of Technology, CH-8092 Zurich,
c SWITZERLAND. All rights reserved.
 
c Please do not distribute this source.

c-----------------------------------------------------------------------

      character*10 function gname (id)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id
 
      character names*8
      common/ cst8  /names(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c-----------------------------------------------------------------------
      if (id.lt.0) then
         gname = names(-id)
      else if (id.gt.0) then 
         gname = fname(id)
      end if 

      end 

      subroutine gexces (id,dg)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from configurational entropy, excess properties, and dqf corrections
c for special cases (internal EoS's or Van Laar gexces does not return
c the excess energy, see routines fexces, gvlaar, and toop).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps,dg
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision exces
      common/ cst304 /exces(m3,k1)
c-----------------------------------------------------------------------
      
      dg = exces(1,id) + t * exces(2,id) + p * exces(3,id)

      end 

      double precision function gfluid (y)
c-----------------------------------------------------------------------
c gfluid returns the fugacities computed from a binary phase described 
c by an internal EoS referenced by cfluid, the composition of the phase
c is y (the amount of component 2).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision fo2, fs2, y

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)
c-----------------------------------------------------------------------
      xco2 = y

      call cfluid (fo2,fs2)

      gfluid = r*t*((1d0-y)*f(1) + y*f(2))

      end 


      subroutine fexces (id,dg)
c-----------------------------------------------------------------------
c gexces evaluates the contribution to the gibbs energy of a pseudocompound
c arising from dqf corrections, and excess properties computed from an
c internal fluid EoS
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision dg,f,fo2,fs2

      common/ cst11 /f(3)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      double precision exces
      common/ cst304 /exces(m3,k1)
c-----------------------------------------------------------------------

      dg = exces(1,id) + t * exces(2,id) + p * exces(3,id)

      xco2 = cp(iff(2),id)

      call cfluid (fo2,fs2)

      dg = dg + r*t*(cp(iff(1),id)*f(1) + cp(iff(2),id)*f(2))

      end 

      double precision function gzero (id)
c----------------------------------------------------------------------
c gzero computes the 1 bar reference pressure free energy of a compound 
c identified by the argument 'id'. no fugacity terms are added for real
c gas species (cf., gcpd).   
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, j

      double precision g, ndu, vdp

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision mu
      common/ cst39 /mu(i6)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer eos
      common/ cst303 /eos(k10)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)
c----------------------------------------------------------------------
 
      g = thermo(1,id)
c                                 -sdt
     *      + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *      - t * (thermo(5,id) + (thermo(7,id) - thermo(24,id)*t) * t))
     *      - (thermo(6,id) + thermo(10,id) / t) / t
     *      + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)
 
      do j = 1, jmct      
c                                 -ndu   
         g = g - vnumu(j,id) * mu(j)
      end do
c                                 transitions
      vdp = 0d0
      ndu = 0d0 

      if (lct(id).ne.0) call mtrans (g,vdp,ndu,id)

      gzero = g

      end 

      recursive double precision function gcpd (id,proj)
c-----------------------------------------------------------------------
c gcpd computes the gibbs free energy of a compound identified by
c the arguement 'id' from the thermochemical parameters stored
c in the array 'thermo' which is located in common block 'cst1'.
c the parameters are: g(pr,tr) [thermo(1,id)], s(pr,tr) [thermo
c (2,id)], v(pr,tr) [thermo(3,id)], the coefficients of the extended
c maier-kelley heat capacity equation:
 
c           cp(pr,t)=a+b*t+c/(t*t)+d*t*t+e/t**(1/2)+f/t+g/t**3
 
c [thermo(4-10,id)], and the coefficients of the volumetric 
c equations given in the program documentation (Eqs 2.1-2.3):

c        b8 = 0 =>
 
c           v(p,t) = v(pr,tr) + b2*(t-tr) + b4*(p-pr)
c                           + b6*(p-pr)**2 + b7*(t-tr)**2

c        -3d0 < b8 < 0 =>

c           v(p,t) = v(pr,tr) * exp[ b3*(t-tr) + b8*(p-pr) ]

c        b8 < -3d0 =>

c           p = 3 * k * f * (2f + 1)^(5/2) * (1 - 3/2 * f * (4 + b8))
c           f = 1/2 * ((v0/v)^(2/3) - 1)
c           v0 = v(pr,tr) + b1 * (t - tr)
c           k = -v0 / (b2 + b3 * (t - tr))

c        b8 > 0 =>

c           v(p,t) = v(pr,t) * [1 - b8*p/(b8*p+Kt)]^(1/b8)
c               Kt = b6 + b7 * t
c            alpha = b1 + b2*t + b3/t - b4/t^2 + b5/t^(1/2)      

c Ghiorso parameters:

c v(J/bar) dv/dt(J/K/bar) dv/dp(J/bar^2) dv/dp/dt(J/bar^2/K) -K'

c Perple_X parameters (1st index in array thermo):

c    3            11           12             13              18

c N.B. for some reason Ghiorso chooses a reference T = 1673 for tbe
c mechanical EoS parameters, this is hardwired as variable trv.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, iwarn, oldid, j

      logical proj

      double precision ialpha, vt, trv, pth, vdp, ndu, vdpbm3, gsixtr, 
     *                 gstxgi, fs2, fo2, kt, gval, gmake, gkomab, kp,
     *                 a, b, c, gstxlq, glacaz, v1, v2, gmet, gterm2, 
     *                 km, kmk, lnfpur, gaq, ghkf

      external vdpbm3, gsixtr, gstxgi, gmake, gkomab, gstxlq, glacaz, 
     *         gaq,    lnfpur, gmet, gterm2, ghkf

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)   

      character*8 names
      common/ cst8   /names(k1) 

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer make
      common / cst335 /make(k10)

      integer eos
      common/ cst303 /eos(k10)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      double precision f
      common/ cst11 /f(3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      double precision mu
      common/ cst39 /mu(i6)

      save kt,trv,iwarn,oldid 
      data kt,trv,iwarn,oldid/0d0,1673.15d0,0,0/
c---------------------------------------------------------------------

      if (make(id).ne.0) then 
c                                 the phase is a made phase, compute 
c                                 and sum the component g's.
         gval = gmake (id)
         goto 999

      else if (eos(id).eq.5) then
c                                 sixtrude 05 JGR EoS 
         gval = gsixtr (id)
         goto 999

      else if (eos(id).eq.6) then
c                                 stixrude JGI '05 Eos
         gval = gstxgi (id) 
         goto 999
         
      else if (eos(id).eq.11) then
c                                 stixrude EPSL '09 Liquid Eos
         gval = gstxlq (id) 
         goto 999

      else if (eos(id).eq.12) then
c                                read SGTE data and evaluate EOS after Brosh '07,'08
         gval = gmet (id)
         goto 999

      else if (eos(id).eq.14) then
c                                read SGTE data and evaluate EOS after Brosh'07,'08
c                                (modified for liquid carbon)
         gval = gterm2 (id)
         goto 999

      else if (eos(id).eq.15) then 
c                                Anderson density extrapolation aqueous species EoS
         gval = gaq (id)
         goto 999

      else if (eos(id).eq.16) then 
c                                DEW/HKF aqueous species formulation
         gval = ghkf (id)
         goto 999

      end if 
c                                 all other EoS's with Cp function
      gval = thermo(1,id)
c                                 -sdt
     *      + t * (thermo(2,id) - thermo(4,id) * dlog(t)
     *      - t * (thermo(5,id) + (thermo(7,id) - thermo(24,id)*t) * t))
     *      - (thermo(6,id) + thermo(10,id) / t) / t
     *      + thermo(8,id) * dsqrt(t) + thermo(9,id)*dlog(t)
c                                 vdp-ndu term:
      if (eos(id).eq.8) then 
c                                 HP Tait EoS, einstein thermal pressure
         pth = thermo(11,id)*(1d0/(dexp(thermo(15,id)/t)-1d0)
     *         -thermo(19,id))

         v1 = 1d0 + (p -pth)*thermo(17,id)
         v2 = 1d0 + (pr-pth)*thermo(17,id)

         if (v1.lt.0d0.or.v2.lt.0d0) then
c                                 destabilize the phase
            vdp = 1d12*p

            if (iwarn.le.50.and.oldid.ne.id) then 
               call warn (46,t,id,names(id)) 
               iwarn = iwarn + 1
               oldid = id
               if (iwarn.eq.50) call warn (49,t,46,'GCPD_HP_Tait')
            end if 

          else 
c                                 int(vdp,p=Pr..Pf)
             vdp = (thermo(16,id)*(
     *            (v1**thermo(18,id) - v2**thermo(18,id))
     *            /thermo(20,id)-p+pr)+p-pr)*thermo(3,id)

          end if 

      else if (eos(id).eq.9) then 
c                                 True tait used for melt endmembers
c                                 kt = b6 + b5*(T-Tr)
c                                 vt = v0*exp(b1*(T-Tr))
          kt = thermo(16,id) + thermo(15,id) * (t-tr)         
c                                 tait parameters, "c" is 1 - tait c
          a = thermo(19,id)/(thermo(19,id)+kt*thermo(17,id))
          b = thermo(18,id)/kt-thermo(21,id)
          c = 1d0 - (thermo(19,id)+kt*thermo(17,id))
     *             /(thermo(20,id)-kt*thermo(17,id))
c                                 int(vdp,p=Pr..Pf)
          vdp = ((((p*b+1d0)**c-(Pr*b+1d0)**c)/b/c+pr-p)*a-pr+p)*
c                                 vt
     *          thermo(3,id)*dexp(thermo(11,id)*(t-tr))
 
      else if (eos(id).eq.10) then 
c                                 ideal gas EoS
         vdp = r*t*dlog(p/pr)

      else if (eos(id).eq.13) then 
c                                 komabayashi/omori polynomials, murnaghan EoS
         vt = thermo(3,id) * dexp((thermo(11,id) + thermo(12,id)*t)*t +
     *        thermo(13,id)*dlog(t) + thermo(14,id)/t + thermo(23,id))
         kt = 1d0 / (thermo(15,id) + t*(thermo(16,id) + t*(thermo(17,id)
     *                             + thermo(18,id)*t))) 
c                                 komatose form, k'(T)
         kp = thermo(19,id) + thermo(20,id)*t*dlog(t)

         km = kp - 1d0
         kmk = km/kp

         vdp = vt * kt**(1d0/kp)/km 
     *            * ((kt + kp*p )**kmk - (kt + kp*pr)**kmk)

      else if (thermo(18,id).eq.0d0) then 
c                                 normal polynomial:
          vdp =  p * (thermo(3,id) 
     *               + (thermo(17,id) * t + thermo(12,id)) * t
     *               + (thermo(16,id) * p + thermo(14,id)) * p)

      else if (thermo(18,id).gt.0d0) then
c                                 murnaghan EoS:
c                                 int(alpha,T=Tr..T)
         ialpha = (thermo(11,id) + thermo(12,id)*t)*t +
     *            thermo(13,id)*dlog(t) + thermo(14,id)/t + 
     *            thermo(15,id)*dsqrt(t) + thermo(23,id)

         if (lopt(8)) then 
c                                 use hollad & powell's approximate form
            vt = thermo(3,id)*(1d0 + ialpha)

         else        
c                                 v(t,pr), correct form
            vt = thermo(3,id)*dexp(ialpha)

         end if 
c                                
         if (lopt(4)) then 
c                                 compute kt using Anderson-Gruneisen parameter
c                                 and expansivity ala Helffrich & Connolly 2009.
            kt = thermo(16,id)*dexp(-thermo(21,id)*ialpha)

         else 

            kt = thermo(16,id) + thermo(17,id)*t
c                                 a ****wit has entered a ridiculous
c                                 temperature
            if (kt.lt.0d0) then 

               if (iwarn.lt.50.and.id.ne.oldid) then 
                  call warn (46,t,id,names(id)) 
                  iwarn = iwarn + 1
                  oldid = id
                  if (iwarn.eq.50) 
     *               call warn (49,t,46,'GCPD_Murnaghan')
               end if 
c                                 destabalize the phase
               gcpd = 1d12*p

               return 

            end if 

         end if  
c                                 Murnaghan EoS:
c                                 vdp = V(1,T)*KT**(1/K')/(K'-1)
c                                 * [(KT+K'*p)**(1-1/K') -
c                                 (KT+K'*pr)**(1-1/K')]
         vdp = vt * kt**(1d0/thermo(18,id))/thermo(22,id) 
     *            *((kt+thermo(18,id)*p)**thermo(19,id)
     *            -(kt+thermo(20,id))**thermo(19,id))

      else if (thermo(18,id).lt.-3d0) then
c                                 3rd order Birch-Murnaghan
c                                 Ghirso Eos is a special case
c                                 indicated by thermo(16,id) = K0 = 0
         if (thermo(16,id).eq.0d0) then 
c                                 assume ghiorso's expressions for v0 and k (KT)
c                                 vt = Volume at 1 bar, T
            vt = thermo(3,id) + thermo(11,id) * (t - trv)

            kt = -vt / (thermo(12,id) + thermo(13,id) * (t - trv) )

         else 
c                                 int(alpha,T=Tr..Tf)
            ialpha = (thermo(11,id) + thermo(12,id)*t)*t +
     *               thermo(13,id)*dlog(t) + thermo(14,id)/t + 
     *               thermo(15,id)*dsqrt(t) + thermo(23,id)
c                                 v(t,pr)
            vt = thermo(3,id)*dexp(ialpha)
         
            if (lopt(4)) then 
c                                 compute kt using Anderson-Gruneisen parameter
c                                 and expansivity ala Helffrich & Connolly 2009.
               kt = thermo(16,id)*dexp(-thermo(21,id)*ialpha)

            else 

               kt = thermo(16,id) + thermo(17,id)*t 

            end if 

         end if 
c                                 a ****wit has entered a ridiculous
c                                 temperature
         if (kt.lt.0d0.or.vt.lt.0d0) then 

            if (iwarn.lt.50.and.oldid.ne.id) then 
               call warn (46,t,id,names(id)) 
               iwarn = iwarn + 1
               oldid = id
               if (iwarn.eq.50) call warn (49,t,46,'GCPD_BM3')
            end if 
c                                 destabilize the phase
            vdp = 1d12*p

         else

            vdp = vdpbm3 (vt,kt,thermo(18,id))

         end if 

      else 
c                                 gottschalk.
         vdp = thermo(11,id)*dexp(thermo(13,id)*t)*
     *               (1d0 - dexp(thermo(18,id)*(p - pr)))

      end if

      gval = gval + vdp

c                                 check for transitions:
      if (ltyp(id).ne.0) call mtrans (gval,vdp,ndu,id)

c                                 check for temperature dependent
c                                 order/disorder:
      if (idis(id).ne.0) call disord (gval,idis(id))
c                                 fluids in the saturated
c                                 or thermodynamic composition space, these
c                                 are not used for mixture props.
      if (eos(id).gt.100) then

         if (eos(id).eq.201.or.eos(id).eq.202) then
c                                 species has been identified as a special composant
c                                 and eos is set by ifug
            if (eos(id).eq.201) then 
         
               xco2 = 0d0 
               call cfluid (fo2,fs2)
               gval = gval + r*t*f(1)
            
            else
          
               xco2 = 1d0
               call cfluid (fo2,fs2)
               gval = gval + r*t*f(2)

            end if 

         else if (eos(id).lt.117) then  
c                                 call appropriate pure fluid EoS
            gval = gval + r*t * lnfpur(eos(id))
            
         else if (eos(id).ge.600.and.eos(id).le.603) then
c                                 komabayashi & fei (2010) EoS for Fe
            gval = gkomab(eos(id),id,vdp)

         else if (eos(id).eq.605) then 
c                                 Stoichiometic O rk fluid 
            xco2 = 0d0
            call cfluid (fo2,fs2)
c                                 real O fluid (O and O2 species)
            gval = gval + r*t*f(1)
c                                 this is -RT(lnk2+lnk3)/2 (rksi5 k's)
c    *         -0.3213822427D7 / t + 0.6464888248D6 - 0.1403012026D3*t

         else if (eos(id).ge.610.and.eos(id).le.637) then
c                                 lacaze & Sundman (1990) EoS for Fe-Si-C alloys and compounds
c                                 Xiong et al., 2011 for Fe-Cr alloys
            gval = gval + glacaz(eos(id)) + vdp + thermo(1,id)  
         
         end if          

      end if
c                                 kill melt endmembers if T < T_melt 
999   if (ifp(id).lt.0.and.t.lt.nopt(20)) gval = gval + 1d6
c                                 do legendre transform if proj
      if (proj) then 
c                                 mobile components
         do j = 1, jmct      
            gval = gval - vnumu(j,id) * mu(j)
         end do

      end if 

      gcpd = gval

      end

      subroutine zeroi (iarray,index,ivalue,n)

      implicit none

      integer n, iarray(n), index, ivalue, i

      do i = 1, index
         iarray(i) = ivalue
      end do 

      end

      subroutine xchk (xmin, xmax, xinc, tname) 

      implicit none

      double precision xmin, xmax, xinc

      character tname*10

      if (xmax.gt.1d0) then 
         call warn (21,xmax,1,tname)
         xmax = 1d0
      end if 

      if (xmin.lt.0d0) then 
         call warn (22,xmin,1,tname) 
         xmin = 0d0
      end if 

      if (xmax.lt.xmin) then
         call warn (23,xmax,1,tname)
         xmax = 1d0
         xmin = 0d0
      end if 

      if (xinc.le.0d0) then 
         call warn (23,xinc,1,tname)
         xinc = 1d0       
      end if

      end

      logical function badz (z)
c----------------------------------------------------------------------
      implicit none

      double precision z

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      if (z.gt.-zero.and.z.le.r1) then
         badz = .false.
      else 
         badz = .true.
      end if

      end 

      subroutine assptx
c---------------------------------------------------------------------
c assptx assigns the current value of v1 and v2 to the array ptx and
c increments the ordinate counter ipt2.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ipt2
      double precision ptx 
      common/ cst32 /ptx(l5),ipt2

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      ipt2 = ipt2 + 2

      if (ipt2.gt.l5) ipt2 = l5

      ptx(ipt2-1) = v(iv1)
      ptx(ipt2)   = v(iv2)
 
      end
 
      subroutine concrt
c---------------------------------------------------------------------
c concrt determines convergence criteria and limits for reasonable
c solutions for the routine univeq. the array delt is also used in
c the routine slope.
 
c input:  dv - array of default search increments.
c         vmax,vmin - arrays containing the maximum and minimum values
c           of the independent and dependent intensive variables.
c         iv - array indexes variables in vmax, vmin, and output arrays
c output: ulim, blim - arrays with upper and lower limits for reasonable
c           solutions (vmax,vmin+/-4*dv)
c         delt - array containing the finite difference increments
c           (vmax-vmin/10**5), delt also serves as the test value
c           for convergence.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i

      double precision ddv
 
      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)
c---------------------------------------------------------------------
      do i = 1, l2
         if (dv(i).lt.0d0) call error (34,dv(i),i,'CONCRT')
         if (i.eq.3) then 
            ulim(i) = vmax(i)
            blim(i) = vmin(i)
         else if (i.eq.1.or.i.eq.2) then
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
            if (blim(i).lt.0d0) blim(i) = 1d0 
         else 
            ulim(i) = vmax(i)+dv(i)
            blim(i) = vmin(i)-dv(i)
         end if 
         ddv = vmax(i)-vmin(i)
         if (ddv.lt.0d0) call error (35,ddv,i,'CONCRT')
      end do 
 
      end
 
      subroutine conver (g,s,v,a,b,c,d,e,f,gg,c8,
     *                   b1,b2,b3,b4,b5,b6,b7,b8,
     *                   b9,b10,b11,b12,b13,tr,pr,r,ieos)
c---------------------------------------------------------------------
c convert thermodynamic equation of state from a pr tr reference state
c to minimize references to to reference conditions and constants.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ieos

      double precision g,s,v,a,b,c,d,e,f,gg,c8,b1,b2,b3,b4,b5,b6,b7,b8,
     *                 b9,b10,b11,b12,b13,tr,pr,n,v0,k00,k0p, dadt0,
     *                 gamma0,q0,etas0,g0,g0p,r,c1,c2, alpha0, beta0, 
     *                 yr,theta,psi,eta

      double precision emodu
      common/ cst318 /emodu(k15)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                constants for anderson electrolyte extrapolation (ieos = 15)
      save alpha0, beta0, dadt0
      data alpha0, beta0, dadt0 /25.93d-5,45.23d-6,9.5714d-6/
c                                constants for hkf electrolyte formulation (ieos = 16)
      save psi, theta, yr, eta 
      data psi, theta, yr, eta/2600d0, 228d0, -5.79865d-5, 694656.968d0/
c----------------------------------------------------------------------
c                                first conditional reformulates and returns for eos:
c                                      1, 5, 6, 11, 12, 15, 16, 101, 102
c                                reformulates and continues to second conditional for
c                                      eos < 100 and special cases (see final conditional)
      if (ieos.eq.1) then 
c                                G(P,T) polynomial forms, e.g., Helgeson et al 1978 (AJS)
c                                Berman 1988 (J Pet).
         g  = g
c                                add the integral of sdt at tr,pr:
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *       + c8 / 4d0 * tr**4
c                                add the constant components of the
c                                integral -vdp at (t,pr):
     *       - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0
c                                sign on b7 corrected June 16, 2004. PJ Gorman. 
     *       - b6 * pr**3/3d0 - b7 * tr * tr * pr
c                                S* (T)
         s  = a - b2*pr - s + a*dlog(tr) + b*tr - c/tr/tr/2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f/tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 /3d0
c                                b7 term added June 16, 2004. PJ Gorman.
     *        + b7 * 2d0 * pr * tr
c                                V* (P)
         v  = v - b2 * tr - b4 * pr  
c                                b6 term added June 16, 2004. PJ Gorman.
     *          + b6 * pr * pr
c                                sign on b7 corrected June 16, 2004. PJ Gorman.
     *          + b7 * tr * tr
c                                a*  (-T log T)
c        a  = a
c                                b*  (-T*T)
         b  = b7 * pr + b / 2d0
c                                c*  (-1/T)
         c  = c / 2d0
c                                d*  (-T**3)
         d  = d / 6d0
c                                e*  (sqrt T)
         e  = 4d0 * e
c                                f*  (log T)
c        f  = f
c                                gg* (-1/T**2)
         gg = gg/6d0
c                                c8* T**4
         c8 = c8/12d0
c                                b2* (PT)
         b2 = b2
c                                b7 term added June 16, 2004. PJ Gorman.
     *          - 2d0 * b7 * tr 
c                                b3* ((P-Pr) exp (T / c1))
c        b3 = b3
c                                b4* (P**2)
         b4 = b4/2d0
c                                b6 term added June 16, 2004. PJ Gorman.
     *          - b6 * pr 
c                                b5* (exp (P/ c2))
c        b5 = c2 * b5
c                                b6* (P**3)
         b6 = b6/3d0
c                                b7* (P*T**2)
c        b7 = b7

         return

      else if (ieos.eq.5.or.ieos.eq.6) then 
c                              Mie-Gruneisen Solid Models:
         if (ieos.eq.5) then
c                              stixrude & bukowinski JGR '93 +
c                              stixrude & lithgow-bertelloni 2005a (JGR)
            n = s

         else
c                              stixrude & lithgow-bertelloni GJI '05
            n = -s

         end if 

         v0     = -v 
         k00    = a
         k0p    = b
         gamma0 = d
         q0     = e
         etas0  = f
         g0     =  emodu(1)
         g0p    =  emodu(2)            
c                                 nr9
         b1 = 9d0*n*r
c                                 c1
         b2 = 9d0*k00*v0
c                                 c2
         b3 = k0p/2d0-2d0
c                                 c3
         b4 = 3d0*b2*b3
c                                 aii
         b5 = 6d0*gamma0
c                                 aiikk
         b6 = -12d0*gamma0 + 36d0*gamma0**2 - 18d0*q0*gamma0
c                                 as 
         b7 = -(gamma0 + etas0)
c                                 aiikk2 
         b8 = b6/2d0
c                                 aii2   
         b9 = b5/2d0
c                                 nr9t0  
         b10 = b1*tr
         b11 = (3d0*k00*g0p-5d0*g0)
         b12 = ((6d0*g0p-24d0+4.5d0*k0p)*k00-14d0*g0)

         return
         
      else if (ieos.eq.11) then 
c                                Mie-Gruneisen Stixrude liquid Model:
c                                G0 = F0 
c                                S0 = S0 => S0 - Cv  
c                                V0 = V0   
c                                a  = Cv  
c                                b  = K0 => 4.5d0*K0*V0  
c                                c  = K0' => 4.5d0*K0*V0*(K'-4)    
c                                d  = y0 => y0 - y'   
c                                e  = y'    
c                                f  = T0
c                                --- dependent ---
c                                gg = (S0-Cv-Cv*y0)*T0
c                                b1 = Cv*(y0+ln(T0))-S0+Cv
c                                b2 = ln(v0)
         
         gg = (s - a - a*d)*f
         b1 = a*(d+dlog(f)) - s + a
         b2 = dlog(v)

         s = s - a
         b = 4.5d0*b*v
         c = b*(c-4d0)
         d = d - e

         return

      else if (ieos.eq.12.or.ieos.eq.14) then 
c                                 calphad format, don't do anything.
         return 

      else if (ieos.eq.15) then 
c                                 H&P aqueous species, rearrange constants
c                                 and load into thermo(10-14) => (gg,b1,b2,b3,b4)
c                                 the HOH flag is in b3, move to b11
         b11 = b3 

         gg = -s + tr*b + (a - b*tr)/tr/dadt0*alpha0
         b1 =  (a - b*tr)/tr/dadt0  
         b2 = -b/2d0
         b3 = tr*(s - b/2d0*tr) + g - pr*v 
     *        + (a-b*tr)/tr/dadt0*(beta0*pr-alpha0*tr)  
         b4 = v - (a-b*tr)/tr/dadt0*beta0

         return 
 
      else if (ieos.eq.16) then 
c                                 coming in b3 is a flag (HOH) for H+, OH-, returned 
c                                 as b11 (thermo(21))
        b11 = b3 
c                                 DEW/HKF aqueous species
c                                 psi, theta, epsr, yr are generic parameters. 
c                                 coming in HKF species parameters loaded as:
c                                 g, s, v,   a, b, c,  d,   e,  f, gg, b1, b2
c                                 (actually v and cp0 are not loaded)
c                                 and correspond to (HKF notation):
c                                 g, s, v, cp0, w,  q, a1, a2, a3, a4, c1, c2
c                                 the compound constants on output will be 
c        b2 = -s + c1*dlog(tr) + c1 + w*yr + dlog(tr/(tr-theta))*c2/theta**2 => b8 in HKF_G.mws
         b3 = -s + b1*dlog(tr) + b1 + b*yr 
     *            + dlog(tr/(tr-theta))*b2/theta**2
c        b3 = (-w*yr-c1+s)*tr + (-1/epsilonr+1)*w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9
c        may 2, 2017 modified to:
c        b3 = (-w*yr-c1+s)*tr + w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9
         b4 = (-b*yr-b1+s)*tr + b - d*pr 
     *                        - e*dlog(psi+pr) + g + b2/theta
c        b4 = -a3*pr-a4*ln(psi+pr) => b10
         b5 = -f*pr - gg*dlog(psi+pr)
c        b5 = -c2/(tr-theta)/theta => b11
         b6 = -b2/(tr-theta)/theta
c        b6 = c2/theta^2 => b12
         b7 = b2/theta**2
c        b7 = -c1-c2/theta^2
         b8 = -(b1+b7)
c                                 the reference condition born radius (thermo 19) 
         b9 = 5d9 * eta * c**2 / (1.622323167d9 * eta * c + 5d9 * b)

         return 
c                                 remaining standard forms have caloric polynomial
      else if (ieos.le.202.or.ieos.eq.604.or.ieos.eq.605.or.
     *         ieos.eq.606.or.ieos.eq.700.or.ieos.eq.701.or.
     *         ieos.eq.702) then
c                                G(Pr,T) polynomial 
         g  = g
     *       + s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *       - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *       - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *       + c8 / 4d0 * tr**4
         s  = a - s + a * dlog(tr) + b * tr - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 /3d0

         b  = b / 2d0
         c  = c / 2d0
         d  = d / 6d0
         e  = 4d0 * e
         gg = gg/6d0
         c8 = c8/12d0
c                                 fluid special case, this is sloppy
         if (ieos.gt.100.and.ieos.lt.120.or.ieos.eq.201.or.ieos.eq.202)
     *      return

      end if 
c                                 -------------------------------------
c                                 mechanical eos:
      if (ieos.eq.3) then 
c                                 for Ghiorso et al.'s PMELTS formulation
c                                 b6 (K0) is 0 and alpha is a constant to an
c                                 arbitrary reference, therefore parameters are left as is.
      else if (ieos.eq.7) then 
c                                 exponential polynomial on volume, e.g., Haas et al 1986,
c                                 Gottschalk 1997. b3 = alpha, -b8 = beta: 
         b1 = -v / b8 / dexp (b3*tr)

      else if (ieos.eq.8.or.ieos.eq.9) then 
c                                 The HP Tait EoS:
         if (ieos.eq.8) then 
c                                 HP 2010 Full Tait EoS
c                                 on input b1=alpha0, b5=theta, b6=k0, b7=k", b8=k'
c                                 thermal pressure coefficients:
c                                 alpha0/k0/xi0*theta -> b1
            b1 = 1d0/b5*b1*b6*tr**2/dexp(b5/tr)*(dexp(b5/tr)-1d0)**2
            b9 = 1d0/(dexp(b5/tr)-1d0)
c                                 Tait a parameter -> b6
            c1 = (1d0+b8)/(1d0+b8+b6*b7)
c                                 Tait b parameter -> b7
            c2 = b8/b6 - b7/(1d0+b8)
c                                 Tait c parameter -> b8 = 1-c
            b8 = 1d0 - (1d0+b8+b6*b7)/(b8**2+b8-b6*b7)
            b7 = c2
            b6 = c1
            b10 = b7*b8

         else 
c                                 True tait used for melt endmembers
c                                 kt = b6 + b5*(T-Tr)
c                                 vt = v0*exp(b1*(T-Tr))
            b9 = 1d0 + b8
            b10 = b8*b9
            b11 = b7/b9
c                                 tait parameters computed as f(T)
         end if 

      else if (ieos.eq.10) then 
c                                 ideal gas, could make a reference pressure correction here. 
      else if (ieos.eq.13) then 
c                                 1) alpha = b1 + b2*T + b3/T + b4/T^2 
c                                    which is reformulated here to 
c                                    int(alpha,T=Tr..Tf) = b13 + b1*T + b2*T^2 + b3*ln(T) + b4/T 
         b2 = b2/2d0 
         b4 = -b4 
         b13 = -(b1*tr + b2*tr*tr + b3*dlog(tr) + b4/tr)
c                                 2)  K' is a f(T) 
         b9 = b9 - b10*tr*dlog(tr)

      else if (b8.ne.0d0) then 
c                                 All remaining forms (ieos = 2, 4, >100) assume:
c                                 1) alpha = b1 + b2*T + b3/T + b4/T^2 + b5/sqrt(T)
c                                    which is reformulated here to 
c                                    int(alpha,T=Tr..Tf) = b13 + b1*T + b2*T^2 + b3*ln(T) + b4/T + b5*sqrt(T)
         b2 = b2/2d0 
         b4 = -b4 
         b5 = 2d0*b5
         b13 = -(b1*tr + b2*tr*tr + b3*dlog(tr) + b4/tr + b5*dsqrt(tr))
c                                 old parameter test forms:

c                                 2) if lopt(4) not true, then the isothermal bulk modulus is
c                                    K = b6 + b7*(T-Tr)
c                                    which is reformulated to 
c                                    K = b6 + b7*T
         if (.not.lopt(4)) b6 = b6 - b7*tr
c                                 operation solving constants for the Murnaghan (ieos=2), these are not
c                                 used for the BM3.
         b9 = (1d0-1d0/b8) 
         b10 = b8*pr 
         b12 = b8 - 1d0 
c                                 anderson-gruneisen parameter is assumed = K' (abs(b8)) except for
c                                 special EoS forms           
         if (ieos.gt.300) then 
            b11 = -s
         else 
c                                 special EoS, anderson-gruneisen stored in 
c                                 s-position. 
            b11 = dabs(b8) 
         end if

      end if 

      end
 
      subroutine chkphi (ichk,name,good)
c-----------------------------------------------------------------------
c ichk = 0 and 2 -> test for saturated entities
c ichk = 1 and 3 -> test for non-saturated entities
c ichk = 4 -> look for phases that consist entirely of constrained components
c ichk > 1  and not 4 -> do not compare against excluded list (for make definitions). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*8 name
 
      integer ichk,i,j

      logical good

      integer iam
      common/ cst4 /iam
 
      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k0)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      good = .true.
c                               reject the data if excluded in input:
      if (ichk.lt.2.or.ichk.eq.4) then 
         do i = 1, ixct
            if (name.eq.exname(i)) goto 90
         end do
      end if 
c                               reject phases which have a component
c                               not included in the input:
      do i= 1, icmpn
         if (icout(i).eq.0.and.comp(i).ne.0d0) goto 90
      end do
c                               reject phases with negative/zero compositions
      tot = 0d0

      do j = 1, icmpn
         if (comp(j).lt.0d0.and.comp(j).gt.-nopt(5)) then
            comp(j) = 0d0
         else if (comp(j).lt.0d0.and.(ieos.ne.15.and.ieos.ne.16)) then
c                               use ichk to avoid multiple messages
            if (ichk.eq.1.and.iam.eq.1.or.iam.eq.2.or.iam.eq.4.or.
     *          iam.eq.5.or.iam.eq.6.or.iam.eq.15) 
     *                                       call warn (13,tot,j,name)

         end if

         tot = tot + comp(j)

      end do 

      if (tot.eq.0d0) goto 90 
c                               do a check to make sure that the phase does
c                               not consist of just mobile components
      tot = 0d0

      do j = 1, icp + isat + ifct
         tot = tot + comp(ic(j))
      end do  

      if (jmct.gt.0.and.ichk.ne.4.and.tot.eq.0d0) goto 90

      if (lopt(5).and.iam.ne.5.and.iam.ne.6) then 
c                               auto_exclude any phase without a pressure EoS from the 
c                               thermodynamic-saturated component data space, the 
c                               more general exclusion would preclude the use of such data
c                               unless it is explicitly identified as the basis for a fugacity
c                               variable. 
          if (tot.ne.0d0.and.ieos.eq.0) then
c                               got a bad operator, check that it doesn't match a 
c                               special component
             good = .false.

             if (lopt(7)) then

                do i = 1, ispec
                   if (name.eq.cmpnt(idspe(i))) then 
                      good = .true.
                      exit
                   end if 
                end do 

             end if

             if (.not.good) then

                if (ichk.eq.1.and.(iam.eq.1.or.iam.eq.15)
     *                        .or.iam.eq.2.or.iam.eq.4) 
     *                                 call warn (16,tot,j,name)

                goto 90

             end if 

          end if 

       end if 
c                               the following is not executed by build:
c                               if ichk ne 0 reject phases that consist entirely
c                               of saturated components: 
      if (ichk.eq.0.or.ichk.eq.2) return
c                               phases with null composition, saved if ichk = 4
c                               otherwise rejected.
      tot = 0d0

      do j = 1, icp
         tot = tot + comp(ic(j))
      end do 

      if (tot.eq.0d0.and.ichk.ne.4) then
     
         goto 90

      else if (tot.ne.0d0.and.ichk.eq.4) then

         goto 90

      else if (ichk.eq.4.and.tot.eq.0d0) then 
c                               reject a null phase if it contains only 
c                               saturated components, since these phases
c                               are already saved in the sat list.

c                               null phase option over-ride
         if (.not.lopt(37)) goto 90

         tot = 0d0

         do j = icp + 1, icp + ifct
            tot = tot + comp(ic(j))
         end do 

         do j = icp + ifct + isat + 1, icp + ifct + isat + jmct
            tot = tot + comp(ic(j))
         end do

         if (tot.eq.0d0) goto 90

      end if

      return

90    good = .false.
 
      end

      logical function findph (igood)
c-----------------------------------------------------------------------
c findph checks if a phase loaded by getphi consists entirely of 
c component igood, if it does it returns true. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer igood, i

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
c-----------------------------------------------------------------------

      if (comp(igood).ne.0d0) then 
c                                 the phase has the component,
c                                 does it have others?
         do i= 1, icmpn
            if (i.ne.igood.and.comp(i).ne.0d0) then
               findph = .false.
               return
            end if  
         end do

         findph = .true.

      else 

         findph = .false.

      end if 

      end

      subroutine lambw (dg,ld)
c---------------------------------------------------------------------
c calculate the energy of an order-disorder transition using the 
c Bragg-Williams model (Holland and Powell, '96), 0-d speciation.
c    input: ld - pointer to the phase in therlm
c   output: dg - energy change of ordering
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ld

      double precision dg,h,w

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c                                 enthalpy of complete disordering
      h = therlm(1,1,ld) + therlm(2,1,ld)*p
c                                 interaction energy
      w = therlm(3,1,ld) + therlm(4,1,ld)*p
 
      call speci0 (dg,h,w,therlm(5,1,ld),therlm(6,1,ld),
     *                    therlm(7,1,ld),therlm(8,1,ld))
 
      end

      subroutine lamla0 (dg,intvdp,ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using  the
c     Landau model as implemented INCORRECTLY in thermocalc pre-ds6. 

c     The correct implementation of the HP98 Landau model is given
c     by function lamla1.

c     input variables
 
c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm
c                        intvdp = the vdp integral of the phase
 
c     returned - dg - difference in free energy due to the transition

c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ld

      double precision dg,tc,tc0,q2,intvdp
 
      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      tc0 = therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 the hp98 form is 
         q2 = dsqrt(1d0-t/tc)

      else 

         q2 = 0d0

      end if 
c                                 See landau_d55.mws
c                                                JADC Jan 26, 2012.
      dg = therlm(2,1,ld) * 
c                                 This is the ds5 PX version
     *   ( (t-tc)*q2*0.6666667d0 - therlm(8,1,ld)*t + therlm(4,1,ld) )
c                                 This is the ds5 TC version
c    * ((tc0/3d0*q2**2 + (t-tc))*q2 + therlm(7,1,ld) - t*therlm(8,1,ld))
c                                 + int(vt,p)
     *     + therlm(6,1,ld)*intvdp
 
      end

      subroutine lamla1 (dg,intvdp,ld)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using  the
c     Landau model (Holland and Powell '98) but as apparently 
c     implemented for DS6 (i.e., HP 2010) 
 
c     input variables
 
c                        t = temperature in k
c                        p = pressure in bar
c                        ld = pointer to phase in therlm
c                        intvdp = the vdp integral of the phase
 
c     returned - dg - difference in free energy due to the transition
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ld

      double precision dg,tc,tc0,q2,intvdp
 
      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      tc0 =  therlm(1,1,ld)
      tc = tc0 + therlm(3,1,ld)*(p-pr)

      if (t.lt.tc) then
c                                 partially disordered:
         q2 = dsqrt((tc-t)/tc0)

      else 

         q2 = 0d0

      end if 
c                                 This is the hp98 version, differs
c                                 from what's in the TC code such that
c                                 dGPX - dGTC = (tc0-tc)*q2^3/3.
c                                 See landau_d60.mws
c                                                JADC Jan 26, 2012.
c      dg = therlm(2,1,ld) * 
c     *  (therlm(7,1,ld) + t*(q2 - therlm(8,1,ld)) + tc*(q2**3/3d0 - q2))
c        + vdp...
c                                 TC version, according to hans vrijmoed 
c                                 e-mail 9/10/12 this is correct, should 
c                                 check again....
      dg = therlm(2,1,ld) * 
     *   (therlm(7,1,ld) + t*(q2 - therlm(8,1,ld)) 
     *                   - tc*q2 + tc0*q2**3/3d0)
c                                 + int(vt,p)
     *     + therlm(6,1,ld)*intvdp
 
      end

      double precision function gtrans (id,j)
c-----------------------------------------------------------------------
c gtrans computes the reference pressure free energy of a compound 
c aboves its jth transition. 
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j
 
      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
 
      gtrans = therlm(12,j,id)
c                                 -sdt
     *       + t * (therlm(3,j,id) - therlm(5,j,id) * dlog(t)
     *       - t * (therlm(6,j,id) + therlm(8,j,id) * t))
     *       - (therlm(7,j,id) + therlm(11,j,id) / t) / t
     *       + therlm(9,j,id) * dsqrt(t) + therlm(10,j,id)*dlog(t)

      end 

      double precision function gclpht (id,j)
c-----------------------------------------------------------------------
c gclpht computes the reference pressure free energy of a compound 
c aboves its jth transition (SGTE data type)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j
 
      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
 
      gclpht = therlm(5,j,id) + therlm(6,j,id)*t + 
     *         therlm(7,j,id)*t*dlog(t) + therlm(8,j,id)/t + 
     *         therlm(9,j,id)/t**2 + therlm(10,j,id)/t**3 + 
     *         therlm(11,j,id)/t**9 + therlm(12,j,id)*t**2
c                                 added terms for eleanor, for SGTE 1st therlm indices > 3
c                                 are incremented by 1 relative to input array tm index
c                                 parameters are otherwise unmodified. thus therlm(13...)
c                                 is tm(12...) = t12. JADC, 12/3/2017.
     *        + therlm(13,j,id)*t**3 + therlm(14,j,id)*dsqrt(t)
     *        + therlm(15,j,id)*dlog(t) 
     
      end   

      subroutine lamhel (p,t,g,vdp,ld,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS).
 
c     there is something seriously wrong in this routine!!!!
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        lct = number of transitions
c                        g = initial free energy
 
c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ld,lct,i,jtran
      double precision t,g,gtrans,vdp,trtp,p,dt,pstar

      external gtrans
 
      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
      if (t.lt.therlm(1,1,ld)) return

      do i = 1, lct 
         if (t.lt.therlm(1,i,ld)) then 

            if (i.eq.1) then
c                                 T<T lowest transition, ignore
c                                 possibility of < 0 clapeyron slope
c                                 and exit:
                return
             else
c                                 at the i-1 th transition:
                jtran = i - 1
                exit
             end if 
          else 
             jtran = i 
          end if 
      end do 

      g = gtrans(ld,jtran) + vdp 
c                                 add vtrans*dp terms, this is
c                                 only set up for one transition
c                                 with a non-zero clapeyron.
c                                 should write warning. 
      if (therlm(2,1,ld).eq.0d0) return

      trtp = therlm(1,1,ld) + (p-pr)/therlm(2,1,ld)
      dt = t - therlm(1,1,ld)

      if (t .gt. trtp) then
c                                 1 bar polymorph is stable at p,
c                                 p*thermo(3,id) is 0 bar standard
c                                 state pdv polymorph contribution &
c                                 therlm(4,j,ld) is Delta V of trans
c                                 contribution.
         pstar = dt*therlm(2,1,ld) + pr
 
         g = g + (p-pstar)*therlm(4,1,ld)
 
      else
c                                 the 1 bar polymorph isn't stable.
         g = g + dt * therlm(2,1,ld) * therlm(4,1,ld)
 
      end if
 
      end


      subroutine calpht (t,g,ld,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a standard cp transition.
 
c     input variables

c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        lct = number of transitions
c                        g = initial free energy
 
c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ld, lct, i, jtran
      double precision t, g, gclpht 

      external gclpht
 
      double precision therlm,therdi
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
      if (t.lt.therlm(1,1,ld)) return

      do i = 1, lct 
         if (t.lt.therlm(1,i,ld)) then 

            if (i.eq.1) then
c                                 T<T lowest transition, ignore
c                                 possibility of < clapeyron slope
c                                 and exit:
                return
             else
c                                 at the i-1 th transition:
                jtran = i - 1
                exit
             end if 
          else 
             jtran = i 
          end if 
      end do 

c                                    SGTE data format
      g = gclpht (ld,jtran) 
 
      end
 
      subroutine lamqtz (p,t,g,ndu,ld,id)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of helgeson et al 1978 (AJS). eq. (114) (corrected for the
c     misplaced parenthesis and the sign on c sub alpha), this is
C     probably bogus for coesite note that only one transition
c     is allowed.
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        ld = pointer to phase in therlm
c                        g = initial free energy
 
c     constants
c                Tr = 298.15 K
c                Pr = 1. bar
c                S(lope) = 38.5 bar/K
c                ba = 0.066 j/bar
c                aa = 549.82 K
c                ca = -4.973d-6 j/bar/bar
c                trt = 848. K
 
c     returned - g - modified phase free energy
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,ld

      double precision p,t,ps,g,ndu,pdv
 
      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      double precision trt, tr, pr, s, aa, ba, ca , vb
 
      save trt, tr, pr, s, aa, ba, ca , vb
 
      data trt, tr, pr, s, aa, ba, ca, vb / 848., 298.15, 1d0,
     *                  38.5, 549.82, 0.066, -4.973d-6, 2.372 /
 
c      trtp = trt + (p-pr) / s
 
      if (t.gt.trt) then
         ps = (t-trt) * therlm(2,1,ld) + pr
      else
         ps = pr
      end if
c                                 if above the ref P transition T
c                                 do the cp integration:
      if (t.gt.trt) g = therlm(8,1,ld) + (p-ps) * thermo(3,id)
     *                - therlm(3,1,ld) * (t-trt) 
     *     + therlm(5,1,ld) * (t - trt - t * dlog(t/trt))
     *     - (therlm(7,1,ld) + therlm(6,1,ld) * t * trt * trt)
     *     * (t - trt)**2 / 2d0 / t / trt / trt
 
c                                 now add in pdv terms to
c                                 the free energy, note that
c                                 pdv term for the ref polymorph
c                                 is already in:
 
      pdv  =  vb * (ps-pr)
     *      - ca * ( (2d0 * pr * (p-ps) - (p * p - ps * ps) ) / 2d0
     *              + s * (t-tr) * (p-ps) )
     *      + s * (ba + aa*ca*s) * (t-tr)
     *                           * dlog ((aa+p/s) / (aa + ps/s))
 
      g = g + pdv + ndu
 
      end
 
      subroutine lamubc (p,t,gspk,k,lct)
c---------------------------------------------------------------------
c     calculate the extra energy of a lamdba transition using model
c     of berman and brown (1985, contribs. miner. petro.)
 
c     input variables
c                        p = pressure in bars
c                        t = temperature in k
c                        k = pointer to phase in therlm
 
c     returned - gspk - energy because of lambda spike
 
c     modified from perkins et al 1986 (computers and geosciences).
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer j,lct,k

      double precision gspk,ctrans,aspk2,bspk2,ct2,ct3,a1,b1,c1,t92,t93,
     *                 tr92,tr93,dhspk,dsspk,t9,tr,teq,tq1bar,p,tr9,
     *                 dvdtr,dvdp,dstr,abspk,t
 
      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)
c---------------------------------------------------------------------
      gspk=0d0
 
      do j = 1, lct
 
         tq1bar = therlm(3,j,k)
 
         if (tq1bar.eq.0d0) cycle
 
            tr    = therlm(7,j,k)
            teq   = therlm(4,j,k) * (p-1d0) + tq1bar
            ctrans = tq1bar - teq
            tr9    = tr - ctrans
 
            if (t.lt.tr9) cycle 
 
            aspk2 = therlm(1,j,k)
            bspk2 = therlm(2,j,k)
            dvdtr = therlm(5,j,k)
            dvdp  = therlm(6,j,k)
            dstr  = therlm(8,j,k) / tq1bar
            abspk = therlm(9,j,k)
 
            if (t .gt. teq) then
               t9 = teq
            else
               t9 = t
            end if
 
            ct2 = ctrans * ctrans
            ct3 = ct2 * ctrans
 
            a1 = aspk2 * ctrans + 2d0 * abspk * ct2 + bspk2 * ct3
            b1 = aspk2 + 4d0 * abspk * ctrans + 3d0 * bspk2 * ct2
            c1 = 2d0 * abspk + 3d0 * ctrans * bspk2
 
            t92 = t9 * t9
            t93 = t9 *t92
            tr92 = tr9 * tr9
            tr93 = tr92 * tr9
 
            dhspk = a1 * (t9 - tr9)
     *            + b1 * (t92 - tr92) / 2d0
     *            + c1 * (t93 - tr93) / 3d0
     *            + bspk2 * (t9*t93 - tr93*tr9) / 4d0
 
            dsspk = a1 * (dlog(t9) - dlog(tr9))
     *            + b1 * (t9 - tr9)
     *            + c1 * (t92 - tr92) / 2d0
     *            + bspk2 * (t93 - tr93) / 3d0
 
            gspk = gspk - (t9 * dsspk) + dhspk
 
            if (t.gt.teq) gspk = gspk - (dstr + dsspk ) * (t-teq)
 
            gspk = gspk + dvdtr * (p-1d0) * (t9-tr)
     *                  + dvdp * ((p*p-1d0)/2d0 - (p-1d0))
 
      end do 

      end

      subroutine disord (gval,id)
c----------------------------------------------------------------------
c compute t-dependent disorder contribution, g is the
c gibbs energy of disordering, therdi(8,id) is the t of onset of
c disordering, therdi(9,id) is the t of complete disorder.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision gval,dh,tt,ds,trr

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      trr = therdi(8,id)
      if (t .lt. trr) return

      tt  = t
      if (t.gt.therdi(9,id)) tt = therdi(9,id)
 
      dh = therdi(1,id) * (tt - trr)
     *     + 2d0 * therdi(2,id) * (dsqrt(tt) - dsqrt(trr))
     *     - therdi(3,id) * (1d0 / tt - 1d0 / trr)
     *     + therdi(5,id) * dlog(tt/trr)
     *     + therdi(6,id) * (tt*tt - trr*trr) / 2d0
     *     + therdi(7,id) * (tt**3 - trr**3) / 3d0
 
      ds = therdi(1,id) * dlog(tt/trr)
     *     - 2d0 * therdi(2,id) * (tt**(-0.5d0) - trr**(-0.5d0))
     *     - therdi(3,id) * (1d0/tt/tt - 1d0/trr/trr) / 2d0
     *     - therdi(5,id) * (1d0/tt - 1d0 / trr)
     *     + therdi(6,id) * (tt - trr)
     *     + therdi(7,id) * (tt*tt - trr*trr) / 2d0
 
      gval = gval + dh - (t * ds)

      if (therdi(4,id).ne.0d0) gval = gval + dh/therdi(4,id) * (p - pr)

      end
 
      subroutine loadit (id,make,nchk)
c---------------------------------------------------------------------
c loadit loads descriptive data for phases and species (name,comp,
c and therm) into the appropriate arrays (names,comps,thermo,vf,
c and vs).  the arguement 'id' indexes the phase in the arrays.
c note that loadit also computes reference state constants which
c are dependent on the state function being used and on its
c analytical expression.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer id,i,j,k

      logical make, nchk 

      double precision gzero
      external gzero

      double precision z(14),smax,t0,qr2,vmax,dt,g1,g2
 
      character*8 names
      common/ cst8   /names(k1)

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      character*8 name
      common/ csta6 /name

      double precision emodu
      common/ cst318 /emodu(k15)

      integer ic
      common/ cst42 /ic(k0)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ikp
      common/ cst61 /ikp(k1)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer eos
      common/ cst303 /eos(k10)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer junk
      double precision del, rand
      common/ cst321 /del(11,k10),rand(12,k10),junk

      double precision delta
      common/ cst325 /delta(11)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      character*8 eoscmp
      common/ cst98 /eoscmp(2)

      integer iam
      common/ cst4 /iam

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct
      
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision fwt
      common/ cst338 /fwt(k10)
c---------------------------------------------------------------------

      if (id+1.gt.k10) call error (1,0d0,id,'k10')

      ipoint = iphct
c                               check for duplicates
      if (nchk) then 
         do i = jmct + 1, iphct
            if (name.ne.names(i)) cycle
            call error (73,g1,i,name) 
         end do  
      end if       
c                               load name and phase flag
      names(id) = name
c                               moduli flags, indicate availability of
c                               bulk and shear modulus
      if (ieos.eq.5.or.ieos.eq.6) then 
c                               stixrude formulations, both available
         iemod(id) = 2

      else 
c                               if ikind = 0, no explicit moduli 
c                                        = 1, just shear
c                                        = 2, shear and bulk   
         iemod(id) = ikind

      end if 

      ikp(id) = 0
      ifp(id) = 0
      idis(id) = 0
      lmda(id) = 0

      eos(id) = ieos

      if (lopt(7)) then

         do k = 1, ispec

            if (name.ne.cmpnt(idspe(k)).and.name.ne.eoscmp(k)) cycle 

            if (ifct.eq.0) then 
c                                 there is no saturated phase
c                                 assign it the default molecular fluid eos
               eos(id) = 200 + k 

            else if (k.eq.1.and.idfl.ne.2.or. 
     *               k.eq.2.and.idfl.ne.1) then 
c                                 saturated phase, and it's not component(k), ergo
c                                 will be computed by ufluid. this only will work
c                                 for ispec < 3.
               eos(id) = ieos 

            end if
c                                 set fluid flag, this flag is 
c                                 used only to match fluid endmembers
c                                 with fluid pseudocompounds
            ifp(id) = 1  
c                                 gflu used to indicate whether a fluid is 
c                                 in the calculation. not clear why gflu
c                                 is set if saturated phase (formerly it
c                                 was set only if saturated phase was 
c                                 ABSENT, 1/5/2017).
            gflu = .true. 

            exit 

         end do 
          
      end if 
c                                 use ieos flag to signal melt endmembers
c                                 in ifp array, this is only used by gcpd.
      if (eos(id).eq.3.or.eos(id).eq.9.or.eos(id).eq.11) then

         ifp(id) = -1

      else if (eos(id).eq.10.or.eos(id).ge.101.and.eos(id).le.116.or.
     *         eos(id).eq.201.or.eos(id).eq.202) then 
 
         gflu = .true.
         ifp(id) = 1

      end if 
c                               load stoichiometry of components.
      fwt(id) = 0 

      do i = 1, icomp
         cp(i,id) = comp(ic(i))
         fwt(id) = fwt(id) + cp(i,id)*atwt(i)
      end do
c                               convert to kg/mol
      fwt(id) = fwt(id)/1d3
c                               compositional array for frendly
      if (iam.eq.5.and.id.le.k5) then 
         do i = 1, k0
            cp0(i,id) = comp(i)
         end do 
      end if 
c                               and just mobile components
      do i = 1, jmct 
         vnumu(i,id) = comp(ic(i+jprct))
      end do 

      if (make) return
c                               if aqueous solute species store name and 
c                               compositional data in special arrays (in 
c                               principle may need vnumu as well).                  
      if (ieos.eq.15.or.ieos.eq.16) then
c                                aqst is only properly initialized by programs that call input2, 
c                                build sets it to -1 as a flag.
c                                
         if (aqst.eq.-1) aqst = iphct - 1

         aqct = iphct - aqst 

         if (aqct.gt.l9) call error (1,r,aqct,'l9 (max aq species)')

         aqnam(aqct) = name
         aqtot(aqct) = 0d0

         do k = 1, icomp
            aqcp(k,aqct) = comp(ic(k))
            if (k.le.icp) aqtot(aqct) = aqtot(aqct) + comp(ic(k))
         end do
c                               locate H+/OH-, at this point HOH is in thermo(13), after
c                               conver it's in thermo(21)
         if (thermo(13,k10).eq.1d0) then 
            ihy = aqct
         else if (thermo(13,k10).eq.2d0) then 
            ioh = aqct
         end if 

      end if
c                               load elastic props if present
      if (iemod(id).ne.0) then 
c                               kmod initialized 0 in main.
         kmod = 1        
         do i = 1, k15
            emod(i,id) = emodu(i)
         end do

      end if 
c                               compute reference state constants,
c                               etc..
      do i = 1, k4
         thermo(i,id) = thermo(i,k10)
      end do 
c                               load errors for MC calculations
      do i = 1, 11
         del(i,id) = delta(i)
      end do 
 
      call conver (
c                               g0, s0, v0
     *             thermo(1,id),thermo(2,id),thermo(3,id),
c                               c1-c8
     *             thermo(4,id),thermo(5,id),thermo(6,id),
     *             thermo(7,id),thermo(8,id),thermo(9,id),
     *             thermo(10,id),thermo(24,id),
c                               b1-b12                          
     *             thermo(11,id),thermo(12,id),thermo(13,id),
     *             thermo(14,id),thermo(15,id),thermo(16,id),
     *             thermo(17,id),thermo(18,id),thermo(19,id),
     *             thermo(20,id),thermo(21,id),thermo(22,id),
c                               b13 on return
     *             thermo(23,id),
c                               ref stuff
     *             tr,pr,r,eos(id))
 
      if (tr.eq.0d0) then
         thermo(1,id) = thermo(1,k10)
         thermo(2,id) = thermo(2,k10)
      end if
c                              lmda transitions:
      if (ilam.ne.0) then

         lamin = lamin + 1

         if (lamin.gt.k9) call error (1,0d0,lamin,'k9')
 
         if (jlam.eq.5) then 
c                                 holland and powell, bragg-williams model:
c                                 enthalpy change of disordering
            therlm(1,1,lamin) = tm(1,1) - pr*tm(2,1)
c                                 volume change of disordering
            therlm(2,1,lamin) = tm(2,1)
c                                 excess enthalpy
            therlm(3,1,lamin) = tm(3,1) 
c                                 excess volume
            therlm(4,1,lamin) = tm(4,1)
c                                 n
            therlm(5,1,lamin) = tm(5,1)
c                                 fac - unused?
            therlm(6,1,lamin) = tm(6,1)
c                                 n+1
            therlm(7,1,lamin) = tm(5,1) + 1d0
c                                 f 
            therlm(8,1,lamin) = tm(5,1)/(tm(5,1) + 1d0)

         else if (jlam.eq.4) then
c                                 holland and powell, landau model:
            do j = 1, ilam       

               smax = tm(2,j)
               t0 = tm(1,j)
               vmax = tm(3,j)
               qr2 = dsqrt (1d0 - tr/t0)

               therlm(1,j,lamin) = t0
               therlm(2,j,lamin) = smax
c                                 this makes therlm(3) dt/dp
               therlm(3,j,lamin) = vmax/smax
c                                 PX ds5 landau
               therlm(4,j,lamin) = (2d0*t0 + tr)*qr2/3d0
c                                 TC ds6 landau
               therlm(7,j,lamin) = t0*(qr2 - qr2**3/3d0)
               therlm(8,j,lamin) = qr2
c                                 Vdp coefficient
               therlm(6,j,lamin) = vmax*qr2/thermo(3,k10)

            end do 

         else if (jlam.eq.1) then
c                              ubc:
            do j = 1, ilam
 
               therlm(1,j,lamin)=tm(1,j)*tm(1,j)
               therlm(2,j,lamin)=tm(2,j)*tm(2,j)
               therlm(9,j,lamin)=tm(1,j)*tm(2,j)
 
               do k = 3, 8
                  therlm(k,j,lamin)=tm(k,j)
               end do 
            end do 
 
            
         else if (jlam.eq.2.or.jlam.eq.3) then
c                              helgeson:
            p = pr
            lmda(id) = lamin
            ltyp(id) = jlam
c                              set special eos flag to zero to prevent 
c                              gzero from including special terms (currently 
c                              only ieos 605). 
            eos(id) = 0
c                                 now convert paramters:
            do k = 1, ilam
c                                 load into therlm:
               therlm(1,k,lamin) = tm(1,k)
               therlm(2,k,lamin) = tm(2,k)
c                                 c1-c7 cp coefficients
               do j = 4, 10
                  therlm(j+1,k,lamin) = tm(j,k)
               end do 
c                                 the c8 heat capacity coefficient
               therlm(13,k,lamin) = tm(11,k)

               t = tm(1,k)
c                                 temporary counter value
               lct(id)  = k-1
c                                 g at trt:
               therlm(12,k,lamin) = gzero(id)
c                                 delta v trans:
               therlm(4,k,lamin) = 0d0
               if (tm(2,k).ne.0d0) therlm(4,k,lamin) = tm(3,k)/tm(2,k)
c                                 s + st at trt:
               if (t*1d-3.lt.1d0) then 
                  dt = 1d0
               else
                  dt = 1d-3*t
               end if

               t = tm(1,k) + dt 
               g1 = gzero(id)

               t = tm(1,k) - dt 
               g2 = gzero(id) 

               therlm(3,k,lamin) = tm(3,k) - (g1 - g2)/2d0/dt
c                              streamline the eos:
               do j = 1, 13
                  z(j) = 0d0 
               end do 


               call conver (
c                                 g0,s0, dummy
     *                     therlm(12,k,lamin),therlm(3,k,lamin), z(1),
c                                 c1-c8
     *                     therlm(5,k,lamin), therlm(6,k,lamin),
     *                     therlm(7,k,lamin), therlm(8,k,lamin),  
     *                     therlm(9,k,lamin), therlm(10,k,lamin),
     *                     therlm(11,k,lamin),therlm(13,k,lamin),
c                                 dummies (b1-b12)
     *                     z(14),z(2),z(3),z(4),z(5),z(6),z(7),z(8),
     *                     z(9),z(10),z(11),z(12),z(13),
c                                 ref stuff
     *                     tm(1,k),pr,r,0)

            end do 

            eos(id) = ieos

        else if (jlam.eq.6) then 

           do k = 1, ilam
c                              SGTE format
c                              load into therlm:
              therlm(1,k,lamin) = tm(1,k)
              therlm(2,k,lamin) = tm(2,k)

              do j = 4, 14
c                              nastia originally had t12-t15 here, but there
c                              doesn't appear to be any data with more than
c                              the first 11 parameters
c                              ----------------------------------------------
c                              added t12 (t^3), t13 (sqrt(t)), t14 (ln(t)) for eleanor.
c                              note these DO NOT (and DID NOT) correspond to the coefficients
c                              in Eq 4 as currently stated the thermo data file html. 
c                              JADC, 12/3/2017.
                 therlm(j+1,k,lamin) = tm(j,k)

              end do 

           end do

         else 

            write (*,*) 'no such transition model'
            call errpau

         end if
 
         lmda(id) = lamin
         lct(id)  = ilam
         ltyp(id) = jlam
 
      end if
c                              t dependent order: load berman and brown
c                              parameters, (this should be cleaned up)
      if (idiso.ne.0) then
         idsin = idsin + 1
         idis(id) = idsin
         do j = 1, m8
            therdi(j,idsin) = td(j)
         end do 
      end if
 
      end

      subroutine univeq (i,ier)
c---------------------------------------------------------------------
c univeq computes the equilibrium condition of a univariant reaction.
c using a finite difference newton-raphson approximation. the method
c will fail when the derivative of the state function, with respect to
c the state variable v(i), goes to infinity.
c---------------------------------------------------------------------
      implicit none

      double precision vi, delv, del, u, gr, b, xg

      integer i,j,ier
 
      include 'perplex_parameters.h'

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr
 
      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision ddv,vmax,vmin       
      common/ cst9 /vmax(l2),vmin(l2),ddv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
 
      ier = 0
 
      vi = v(i)
      del = delt(i)
      b = blim(i)
      u = ulim(i)
c                                 phase composition special case:
      if (i.eq.3) then 
         if (vi.lt.1d1*del) then 
            del = dabs(vi)/1d1
         else if (1d0-vi.lt.1d1*del) then 
            del = dabs(1d0-vi)/1d1
         end if 
      end if 
 
      if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30
 
      do j = 1, 100

         call grxn (gr)

         v(i) = vi + del
         call incdep (i)

         call grxn (dgr) 
         xg = dgr
         dgr = dgr - gr

         if (dgr.eq.0d0) exit
 
         delv = gr*del/dgr

         if (dabs(delv/ddv(i)).gt.1d0) then
 
c            v(i) = vi - del
c            del = del*ddv(i)/delv/2d0
c            if (dabs(del/delt(i)).lt.1d-6) goto 30 
c            cycle
c                                  changed 7/13/2014
             delv = dabs(delv)/delv * ddv(i)

         end if 

         vi = vi - delv

         if (vi+dabs(del).gt.u.or.vi-dabs(del).lt.b) goto 30

         v(i) = vi

         call incdep (i)
c                                 used to be on the basis of utol,
c                                 but this allows uniform relative error.
         if (dabs(delv).lt.del) return

      end do 
c                                 error code 1: iv and dv must be
c                                 switched
      ier = 1
      return
c                                 error code 2: value too far out of
c                                 range, refine iv increment.
30    ier = 2
      end

      subroutine unver (g,s,v,a,b,c,d,e,f,gg,c8,b1,b2,b3,b4,b5,b6,b7,b8,
     *                  b9,b10,b11,tr,pr,ieos)
c----------------------------------------------------------------------
c convert thermodynamic equation of state from a 0-0 reference state
c to a pr-tr reference state.

c corrections made corresponding to PJ's corrections in conver, 
c June 16, 2004.
c----------------------------------------------------------------------
      implicit none

      integer ieos

      include 'perplex_parameters.h'

      double precision v,gg,a,b,c,d,e,f,g,b1,b2,b4,b5,b6,b7,b8,pr,tr,s,
     *                 c8,b3,b9,b10,b11
c----------------------------------------------------------------------
c                               Stixrude's EoS, Aq, CALPHAD exit without
c                               doing anything
      if (ieos.eq. 5.or.ieos.eq. 6.or.ieos.eq.11.or.ieos.eq.12.or.
     *    ieos.eq.14.or.ieos.eq.15) then 
          
          return

      else if (ieos.eq.16) then 
c                                 HKF electrolyte model. 
          b3 = b11
          return

      end if 

      c8 = 12d0 * c8
      gg = 6d0 * gg
      e  = e / 4d0
      d  = 6d0 * d
      c  = 2d0 * c



      if (b8.eq.0d0) then 
c                                normal vdp term:
         b6 = 3d0 * b6
         b4 = 2d0 * (b4 + b6 * pr)
         b2 = b2 + 2d0 * b7 * tr
         b  = 2d0 * (b - b7*pr)
         v  = v + b2 * tr + b4 * pr - b6 * pr * pr - b7 * tr * tr
         s  = -1d0 * ( s - ( a - b2 * pr  + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 / 3d0
     *        + b7 * 2d0 * pr * tr ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *        - c8 * tr**4 / 4d0
     *        - v * pr + b2 * tr * pr + b4 * pr * pr / 2d0 
     *        - b6 * pr**3/3d0 - b7 * tr * tr * pr )
      else 

         b  = 2d0 * b

         s  = -1d0 * ( s - ( a + a * dlog(tr) + b * tr
     *        - c / tr / tr /2d0
     *        + d * tr * tr / 2d0 - 2d0 * e / dsqrt(tr) - f / tr
     *        - gg / tr**3 / 3d0 + c8 * tr**3 / 3d0 ) )
         g  = g - ( s * tr - a * tr - b * tr * tr / 2d0 + c / tr
     *        - d * tr**3 / 3d0 - 2d0 * e * dsqrt(tr)
     *        - f * dlog(tr) + gg / tr / tr / 2d0 + f
     *        - c8 * tr**4 / 4d0)

         if (b8.gt.0d0.or.(b8.le.-3d0.and.b6.ne.0d0)) then 
c                                 murnaghan or bm3:
            b2 = 2d0 * b2
            b4 = -b4
            b5 = b5 / 2d0  
c                                 convert b6 back to K(Tr)
            b6 = b6 - b7*tr
 
         else if (b8.le.-3d0.and.b6.eq.0d0) then 
c                                 ghirso, do nothing.

         else 
c                                 v = f(exponential(beta))
c                                 zero b1, so users don't get confused.
            b1 = 0d0

         end if  

      end if 
 
      end
 
c-----------------------------------------------------------------------
 
c PLIB -  subroutines common to FRENDLY and VERTEX.
 
c-----------------------------------------------------------------------
 
      subroutine slope (iv1,iv2,s)
c---------------------------------------------------------------------
c slope computes the slope (div(1)/div(2)) of a univariant
c equilibria by approximating the derivatives dg/div(1) and
c dg/div(2) by finite differences.
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer iv(2),iv1,iv2,i

      double precision dg(2),gr,s,gval
 
      double precision delt,dtol,utol,ptol
      common/ cst87 /delt(l2),dtol,utol,ptol

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
 
      iv(1) = iv1
      iv(2) = iv2
      call grxn (gr)

      do i = 1, 2
         v(iv(i)) = v(iv(i)) + delt(iv(i))
         call incdep (iv(i))
         call grxn (gval)
         dg(i) = (gval-gr)/delt(iv(i))
c                                 note the possibility exists here
c                                 for the value of an intensive
c                                 parameter to excede allowed
c                                 limits.  although this is very
c                                 unlikely a limit test like the
c                                 one done in 'univeq' could be
c                                 used.
         v(iv(i)) = v(iv(i)) - delt(iv(i))
         call incdep (iv(i))
      end do 

      s = -dg(2) / dg(1)
 
      end
 
      subroutine switch (div,ivi,ivd,jer)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer ivi,ivd,jer,iovd

      double precision div,vo,s,odiv

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision blim, ulim, dgr
      common/ cxt62 /blim(l2),ulim(l2),dgr
c---------------------------------------------------------------------
c                                 reset intensive variables
      call reptx 
      vo = v(ivi)
c                                 determine the sign for the
c                                 ivd increment
      call slope(ivd,ivi,s)
      jer = 0
c
      if (s.eq.0d0) then
c                                 a zero slope shouldn't occur
c                                 someplace is a bug
         jer = 1
         return
      end if
c
      odiv = div
      div = dv(ivd)
      if (odiv.gt.0d0) goto 10
      if (s.gt.0d0) div = -dv(ivd)
      goto 20
10    if (s.lt.0d0) div = -dv(ivd)
c                                 estimate a new value for v(ivi)
20    v(ivi) = v(ivi)+div/s
c     call incdep (ivi) 
c                                 switch variables
      if ((v(ivi).gt.vmin(ivi)).and.(v(ivi).lt.vmax(ivi))) then
         goto 30
      else if (v(ivi).lt.blim(ivi).or.v(ivi).gt.ulim(ivi)) then 
         jer = 1
         return
      end if         
c                                 call to incdep moved from above 3/2/2011
c                                 to prevent problem with negative T in 
c                                 subinc (calculation with fugacity). 
      call incdep (ivi) 

      div = div/5d0

      if (dabs(div).lt.dv(ivd)/1d6) then
         jer = 1
         return
      end if 

      v(ivi) = vo
      call incdep (ivi) 
      goto 20

30    iovd = ivd
      ivd = ivi
      ivi = iovd

      end

      subroutine reptx
c-----------------------------------------------------------------------
c reset - resets indendent potentials to the last known stable condition 
c along a univariant curve 
c----------------------------------------------------------------------- 
      implicit none

      double precision ptx

      integer ipt2

      include 'perplex_parameters.h'
    
      common/ cst32 /ptx(l5),ipt2

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      v(iv(1)) = ptx(ipt2-1)
      v(iv(2)) = ptx(ipt2)
      call incdp0
   
      end

      subroutine unlam (tm,id)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ilam,id,jd,i,j,k
  
      double precision tm(m7,m6), z(12), g1, g0, s0, gcpd

      external gcpd

      double precision therdi, therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer eos
      common/ cst303 /eos(k10)
c-----------------------------------------------------------------------
      if (ltyp(id).eq.0) return

      jd = lmda(id)

      do i = 1, m7
         do j = 1, m6
            tm(i,j) = 0d0
         end do 
      end do 

      if (ltyp(id).eq.5) then 
c                                 Bragg-Williams model:
          do j = 1, 6
             tm(j,1) = therlm(j,1,jd)
          end do 
        
          tm(1,1) = tm(1,1) + pr*tm(2,1)      

      else if (ltyp(id).eq.4) then
c                                 Landau model        
c                                 HP Landau model
         do j = 1, lct(id)
            tm(1,j) = therlm(1,j,jd)
            tm(2,j) = therlm(2,j,jd)
            tm(3,j) = therlm(3,j,jd) * tm(2,j) 
         end do 

      else if (ltyp(id).eq.1) then
c                                 UBC 
         do j = 1, lct(id)
            tm(1,j) = dsqrt (therlm(1,j,jd))
            tm(2,j) = dsqrt (therlm(2,j,jd)) 
         end do         

      else if (ltyp(id).eq.2.or.ltyp(id).eq.3) then
c                                 Helgeson generic, maybe q/coe too. 
         p = pr
c                                 temporary counter
         ilam = lct(id)

         do i = ilam, 1 , -1
c                                 get s transition:
c                                 load into therlm:
            tm(1,i) = therlm(1,i,jd) 
            tm(2,i) = therlm(2,i,jd) 

            do j = 5, 11
               tm(j-1,i) = therlm(j,i,jd)  
            end do 

            tm(11,i) = therlm(13,i,jd)

            t = tm(1,i)
c                              set transition type to null
c                              for call to gphase 
            lct(id) = i - 1   
c                             -s at trt, this should be 
c                             changed to centered rel diff:
            g1 = gcpd (id,.false.)
            t = t + 1d-3

            tm(3,i) =  (gcpd (id,.false.) - g1)/1d-3

            g0 = therlm(12,i,jd)
            s0 = therlm(3,i,jd)

            do k = 1, 9
               z(k) = 0d0
            end do 

            call unver (
c                                 g0, s0, dummy
     *                  g0,s0,z(1),
c                                 c1-c8 
     *                  tm(4,i),tm(5,i),tm(6,i),tm(7,i),
     *                  tm(8,i),tm(9,i),tm(10,i),tm(13,i),
c                                 dummies 
     *                  z(1),z(2),z(3),z(5),z(6),z(7),z(8),z(9),z(10),
     *                  z(11),z(12),
c                                 ref stuff
     *                  tm(1,i),pr,eos(id)) 

            tm(3,i) = s0 + tm(3,i)

         end do 

         lct(id) = ilam

      end if
      
      end 

      subroutine incdep (ind) 
c-----------------------------------------------------------------------
c either indep or incdp0 are called whenever any primary potential
c variables are changed to reevaluate secondary variables. this 
c cumbersome structure is necessitated by the fact that computational
c variables were mapped directly to the thermodynamic variables in 
c cst5. a more rational strategy, which is used in werami/pssect, is
c to separate the computational and thermodynamic variables.

c incdep:

c  1) if ind = iind, computes the dependent variable if P(T) or T(P)

c  ind  - is the index of the current variable. 
c  iind - is the index of the independent variable that v(idep) is a 
c         function of.
c  idep - is the index of the dependent variable.

c  2) if jmct > 0, assigns the independent chemical potentials
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind
       
      double precision var

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c----------------------------------------------------------------------
      if (ind.eq.iind) then 
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2 
     *                          + var*(c3 + var*c4))) 
c    *                           var*(c3 + var*(c4 + c5*var)))) 
      end if 

      if (jmct.gt.0) call subinc 

      end 

      subroutine incdp1 
c-----------------------------------------------------------------------
c automatically update variables if P(T) or T(P)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
       
      double precision var

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep
c----------------------------------------------------------------------
      if (iind.ne.0) then 

         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2 
     *                          + var*(c3 + var*c4))) 
c    *                           var*(c3 + var*(c4 + c5*var)))) 
      end if 

      end

      subroutine subinc 
c-----------------------------------------------------------------------
c assigns the independent chemical potentials, called by incdep and 
c incdp0
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i
       
      double precision gref, xp, gcpd

      external gcpd

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mu
      common/ cst39 /mu(i6)
c----------------------------------------------------------------------
      do i = 1, jmct

            if (imaf(i).eq.1) then 
c                                 the primary variable is a chemical
c                                 potential.
               mu(i) = v(3+i)

            else 
c                                 the primary variable is a fugacity or
c                                 an activity.
               if (imaf(i).eq.2) then 
c                                 fugacity
                  xp = v(1)
                  v(1) = pr
                  gref = gcpd (idaf(i),.false.)
                  v(1) = xp

               else 
c                                 activity
                  gref = gcpd (idaf(i),.false.)

               end if 

               mu(i) = gref + r*v(2)*v(3+i)*2.302585093d0

             end if 

      end do 
c                                 this line is here to prevent an optimization 
c                                 bug with compaq visual fortran. 
c     xp = xp * xp

      end 

      subroutine incdp0 
c----------------------------------------------------------------------
c incdep 1) conditionally computes the dependent variable if one exists 
c (idep ne 0), 2) assigns independent chemical potentials.

c  idep - is the index of the dependent variable 
c  iind - is the index of the independent variable that v(idep) is a 
c         function of.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
       
      double precision var 

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct
c----------------------------------------------------------------------
      if (idep.ne.0) then 
         var = v(iind)
         v(idep) = c0 + var*(c1 + var*(c2 
     *                          + var*(c3 + var*c4))) 
c    *                            var*(c3 + var*(c4 + c5*var)))) 
      end if 

      if (jmct.gt.0) call subinc

      end 

      double precision function depvar (var) 
c--------------------------------------------------------------------
c depvar computes the dependent variable from the independent variable
c var
c--------------------------------------------------------------------
      implicit none

      double precision var 

      integer iind,idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*c4)))

c      depvar = c0 + var*(c1 + var*(c2 + var*(c3 + var*(c4 + c5*var))))

      end 

      integer function match (idim,ier,name)
c----------------------------------------------------------------------
c find the index of endmember identified by name
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer idim, ier 

      character*8 name

      character mname*8
      common/ cst18a /mname(m4)

      ier = 0

      do match = 1, idim
         if (name.eq.mname(match)) exit
      end do  

      if (match.gt.idim) ier = 1
      
      end 

      subroutine readn (i,idim,tname)
c----------------------------------------------------------------------
c readn - read idim endmember names expected format into mname(i+1..i+idim)

c  name_1 name_2 ... | comment
c  name_n+1....      | comment

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
 
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, idim, ict, i

      character name*8, tname*10 

      character mname*8
      common/ cst18a /mname(m4)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 

      call readcd (n9,len,ier,.false.)
      if (ier.ne.0) goto 90

      ibeg = 1
      ict = i

      do while (ict-i.lt.idim) 
c                                 find the name
         call readnm (ibeg,iend,len,ier,name)
         if (ier.ne.0) goto 90
         ict = ict + 1
         mname(ict) = name

         if (ibeg.ge.len.and.ict.lt.idim) then
            call readcd (n9,len,ier,.false.)
            ibeg = 1
            if (ier.ne.0) goto 90
         end if 

      end do

      return      

90    write (*,1000) tname,(chars(i),i=1,len),name
      
      call errpau
      
1000  format ('**error ver200** READN bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last name read was: ',a,/)

      end 


      subroutine readda (rnums,idim,tname)
c----------------------------------------------------------------------
c readda - reads idim numbers, discarding all data in records 
c          beyond the comment marker "|" or the end of the record.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      integer len, idim, kdim, jdim, i, isnum, ier

      character tname*10, nums*(lchar)

      double precision rnums(*)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)
c----------------------------------------------------------------------
c                                 read card scans for non blank data
c                                 card:
      kdim = 1
      jdim = 0
      isnum = 0 
      len = 0
      ier = 1

      do while (jdim.lt.idim)

         call readcd (n9,len,ier,.true.)
         if (ier.ne.0) exit
c                                 got data, count how many "numbers"
         do i = 1, len
            if (chars(i).ne.' '.and.isnum.eq.0) then
               jdim = jdim + 1
               isnum = 1
            else if (chars(i).eq.' ') then
               isnum = 0
            end if 
         end do
c                                 if there is no comment
c                                 marker there can be more 
c                                 data than anticipated:
         if (jdim.gt.idim) jdim = idim
c                                 write numbers to string
         write (nums,*) (chars(i), i = 1, len),' '
         
         read (nums,*,iostat=ier) (rnums(i), i = kdim, jdim)
         if (ier.ne.0) exit

         kdim = jdim + 1

      end do 

      if (ier.gt.0) then 

         write (*,1000) tname, (chars(i),i = 1, len)
         write (*,1020)
         
         call errpau

      else if (ier.lt.0) then 

         write (*,1010) tname
         write (*,1020)
         
         call errpau

      end if 

1000  format ('**error ver209** READDA bad data, currently',
     *        ' reading solution model: ',/,a,/,'data was:',/,400a)
1010  format ('**error ver210** READDA read to end of file',
     *        ' reading solution model: ',/,a)
1020  format ('READDA was expecting numeric data.',/)

      end 

      subroutine readx (idim,tname)
c----------------------------------------------------------------------
c readx - read excess function for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        W(name-name-...) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, jend, len, ier, iscan, lord, imax, match, idim, 
     *        i, j, iscnlt

      double precision nums(m3) 

      character name*8, begin*5, eod*3, tname*10, values*30, key*22

      logical ok

      external iscnlt, iscan

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      character*2 strgs*3, mstrg, dstrg, tstrg*3, wstrg*3, e16st*3
      common/ cst56 /strgs(32),mstrg(6),dstrg(m8),tstrg(m7),wstrg(m16),
     *               e16st(13)
c----------------------------------------------------------------------

      iterm = 0
      iord = 0 

      call readcd (n9,len,ier,.true.)

      write (begin,'(5a)') (chars(i),i=1, 5)

      if (begin.eq.'ideal') then
         return
      else if (begin.ne.'begin') then 
         goto 90
      end if 

      do i = 1, m1
         do j = 1, m2
            isub(i,j,1) = 0 
         end do 
      end do 
 
      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find expansion type
         if (chars(2).eq.'k'.or.chars(2).eq.'K') then 
            xtyp = 1
         else 
            xtyp = 0
         end if    
c                                 find brackets
         ibeg = iscan (1,len,'(') + 1
         imax = iscan (1,len,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         iterm = iterm + 1
         if (iterm.gt.m1) call error (48,wg(1,1),m1,tname)

         lord = 0

         do while (ibeg.lt.imax)

            call readnm (ibeg,jend,imax,ier,name)
            if (ier.ne.0) goto 90

            lord = lord + 1
            if (lord.gt.m2) call error (49,wg(1,1),m2,tname)

            isub(iterm,lord,1) = match (idim,ier,name)  

            if (ier.ne.0) goto 90         

         end do 

         if (xtyp.eq.0) then 

            ibeg = imax + 2
c                                 read standard form margules pt functions
            if (lord.gt.iord) iord = lord

            call redlpt (nums,ibeg,jend,len,ier) 

            if (ier.ne.0) goto 90 

            do i = 1, m3
               wg(iterm,i) = nums(i)
            end do 

         else
c                                 set "perplex" order 
            rkord(iterm) = 0
            iord = 2 
c                                 rk form, read a new card for each term 
            do j = 1, m17

               do i = 1, m16
                  wk(i,j,iterm) = 0d0
               end do 

            end do

            do  

               ibeg = 1
c                                 check for end of data 
               call readcd (n9,len,ier,.true.)
               write (begin,'(3a)') (chars(i), i = 1, 3)
               if (begin.eq.'end') then
                  return
               else if (begin.eq.'Wk(') then
                  exit
               end if
c                                 we have a data card
               rkord(iterm) = rkord(iterm) + 1

               do 
c                                 locate end of keyword
                  if (ibeg.ge.icom) exit 
                  jend = iscan (ibeg,icom,'=') - 1
                  if (jend.ge.icom) exit
c                                 write keyword
                  write (key,'(22a1)',iostat=ier) (chars(i),i=ibeg,jend)
                  if (ier.ne.0) call error (23,wg(1,1),ier,key) 
c                                 locate data
                  ibeg = iscnlt (jend+2,icom,' ')
                  jend = iscan (ibeg,icom,' ')
c                                 write data 
                  write (values,'(80a1)',iostat=ier) (chars(i),
     *                                                      i=ibeg,jend)
                  if (ier.ne.0) call error (23,wg(1,1),ier,key) 
c                                 shift pointer to next key
                  ibeg = iscnlt(jend,icom,' ')
c                                 assign data
                  ok = .false.

                  do i = 1, m16
                     if (key.eq.wstrg(i)) then 
                        read (values,*,iostat=ier) 
     *                                          wk(i,rkord(iterm),iterm)
                        if (ier.ne.0) call error (23,wg(1,1),ier,key) 
                        ok = .true.
                        exit 
                     end if 
                  end do 

                  if (ok) cycle

                  call error (9,wg(1,1),i,key)

               end do

            end do

         end if 

      end do

      return

90    write (*,1000) tname,(chars(i),i=1,len)
      write (*,1010) name
      
      call errpau
      
1000  format ('**error ver200** READX bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a)
1010  format ('last name read was: ',a,/,
     *        'usually this error is due to a mispelled ',
     *        'endmember name.',/)

      end 

      subroutine readop (idim,jlaar,kstot,reach,stck,norf,tname)
c----------------------------------------------------------------------
c readop - tail of solution model to find optional dqf,
c          van laar size parameters, or reach_increment

c readop - reads data until it finds an     "end_of_model" record

c          van laar data is identified by a "begin_van_la" record
c          dqf data is identified by a      "begin_dqf_co" record
c          or the reach factor is           "reach_increm" record

c readop returns:

c          jlaar = 1 if van laar data found (else 0).
c          idqf  > 0 if dqf data found (else 0).
c          reach, set = 0 if no reach factor found. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical stck, norf

      double precision reach 

      integer ier, idim, jlaar, kstot, i

      character tname*(*), key*22, val*3, nval1*12, nval2*12, nval3*12, 
     *          strg*40, strg1*40

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      jlaar = 0 
      idqf = 0 
      reach = 0d0
      stck = .true.
      norf = .false.

      do 

         call redcd1 (n9,ier,key,val,nval1,nval2,nval3,strg,strg1)

         if (key.eq.     'end_of_model') then 
            
            exit

         else if (key.eq.'begin_model ') then
c                              found new model, current
c                              model does not of end_of_model keyword
            write (*,1000) tname,(chars(i),i=1,length)
            
            call errpau

         else if (key.eq.'begin_van_laar_sizes') then 
c                              read van laar data: 
            jlaar = 1
            call readvl (idim,kstot,tname)
            cycle 

         else if (key.eq.'begin_dqf_corrections') then 
c                              read dqf data:
            call readdq (idim,tname)
            cycle

         else if (key.eq.'reach_increment') then 

            read (val,*) i
            reach = dfloat(i)

         else if (key.eq.'site_check_override') then 

            stck = .false.

         else if (key.eq.'refine_endmembers') then 

            norf = .true.

         else 

            write (*,1010) tname,(chars(i),i=1,length)
            write (*,1020)
            
            call errpau

         end if

      end do  

1000  format (/,'**error ver200** READOP missing "end_of_model"',
     *          ' keyword at end',' of solution model:',a,/)
      
1010  format (/,'**error ver210** READOP bad data, currently',
     *          ' reading solution model: ',a,' data was:',/,400a)

1020  format (/,'This error is most probably due to an out-of-date ',
     *          'solution model file.',//,
     *          'Copy the current version from:',//,
     *          'www.perplex.ethz.ch/perplex/datafiles/',
     *          'solution_model.dat',//)

      end 

      subroutine readvl (idim,kstot,tname)
c----------------------------------------------------------------------
c readvl - read van laar volumes for a solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        alpha(name) number number number

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, len, ier, iscan, imax, match, idim, index

      character name*8, eod*3, tname*10

      double precision nums(m3)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      integer kstot,jend,i,ict,jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod
c----------------------------------------------------------------------

      ict = 0 
 
      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find brackets
         ibeg = iscan (1,len,'(') + 1
         imax = iscan (1,len,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         ict = ict + 1
         if (ict.gt.m4) goto 91

         call readnm (ibeg,jend,imax,ier,name)
         if (ier.ne.0) goto 90

         index = match (idim,ier,name)  
         if (ier.ne.0) goto 90         

         ibeg = imax + 2

         call redlpt (nums,ibeg,jend,len,ier) 

         if (ier.ne.0) goto 90 

         do i = 1, m3
            vlaar(i,index) = nums(i)
         end do 

      end do

      if (ict.lt.kstot) goto 91 

      return

90    write (*,1000) tname,(chars(i),i=1,len),vlaar(i,index)
      write (*,1001)
      call errpau

91    write (*,1010) tname
      call errpau
      
1000  format ('**error ver200** READVL bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

1010  format (' **error ver201** READVL bad data, currently',
     *        ' reading solution model: ',a,/,
     *        ' this model requires 1 size parameter for',
     *        ' each independent endmember, READVL found ',i2,
     *        ' parameters.',/)

      end 

      subroutine readdq (idim,tname)
c----------------------------------------------------------------------
c readvl - read dqf corrections for solution models endmembers, assumes
c data on one line of less than 240 characters, the expected format

c        dqf(name) number number number

c        output:

c          idqf           - the number of endmembers with dqf corrections
c          indq(idqf)     - pointer to corrected endmember in the phase name array
c          dqf(1..3,idqf) - the dqf parameters for the correction

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, len, ier, iscan, imax, match, idim

      character name*8, eod*3, tname*10

      double precision nums(m3)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      integer jend,i,idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf
c----------------------------------------------------------------------

      eod = ' '

      do while (eod.ne.'end')  

         call readcd (n9,len,ier,.true.)
         if (ier.ne.0) goto 90

         write (eod,'(3a)') (chars(i),i=1,3) 
c                                 find brackets
         ibeg = iscan (1,len,'(') + 1
         imax = iscan (1,len,')') - 1

         if (ibeg.gt.len.or.imax.gt.len) cycle
c                                 data found
         idqf = idqf + 1
 
         call readnm (ibeg,jend,imax,ier,name)
         if (ier.ne.0) goto 90

         indq(idqf) = match (idim,ier,name)  
         if (ier.ne.0) goto 90         

         ibeg = imax + 2

         call redlpt (nums,ibeg,jend,len,ier) 
         if (ier.ne.0) goto 90 

         do i = 1, m3
            dqf(i,idqf) = nums(i)
         end do 

      end do

      return

90    write (*,1000) tname,(chars(i),i=1,len),dqf(i,idqf)
      write (*,1001)
      
      call errpau
      
1000  format ('**error ver200** READDQ bad data, currently',
     *        'reading solution model: ',a,' data was:',/,400a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)
1001  format (/,'usually this error is caused by a mispelled ',
     *          'endmember name.',/)

      end 

      subroutine readr (coeffs,enth,inds,idim,nreact,tname)
c----------------------------------------------------------------------
c readr - read stoichiometric reaction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        name "=" (acoef(i), mame(i), i= 2..nreact) Enthalpy_value

c nreact >= 4 + nord for reciprocal reactions
c nreact = -1 on input for ordered/disorder reactions
c enthalpy_value is only read if on input nreact = -1

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, iscan, i, nreact, inds(k7), 
     *        idim, match, iscnlt

      double precision coeffs(k7), enth(3), rnum

      character name*8, tname*10 

      external iscan, iscnlt

      character mname*8
      common/ cst18a /mname(m4)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)
c----------------------------------------------------------------------
      ier = 0 

      call readcd (n9,len,ier,.true.)
      if (ier.ne.0) goto 90

      ibeg = 1
c                                 first name
      call readnm (ibeg,iend,len,ier,name)

      if (ier.ne.0) goto 90

      if (nreact.eq.-1) then 
c                                 if nreact = -1, new name
         idim = idim + 1
         mname(idim) = name
         inds(1) = idim

      else 
c                                 reciprocal sol, get index
         inds(1) = match(idim,ier,name)
         if (ier.ne.0) goto 90
              
      end if 
c                                 find marker '='
      ibeg = iscan (1,len,'=') + 1

      i = 2

      do 
c                                 find a stoichiometric coeff
         call readfr (rnum,ibeg,iend,len,ier)
         if (ier.ne.0) exit 

         coeffs(i) = rnum
c                                 find the name
         call readnm (ibeg,iend,len,ier,name)

         if (ier.ne.0) goto 90

         if (i.gt.k7) call error (1,0d0,i,'k7')

         inds(i) = match(idim,ier,name)
         if (ier.ne.0) goto 90

         if (nreact.gt.0.and.i.eq.nreact) exit 
c                                 increment counter and find 
c                                 the next coeff+name
         i = i + 1
           
      end do

      if (nreact.eq.-1) then      
c                                 ordered compound, read
c                                 enthalpy, find marker '='
         ibeg = iscan (ibeg,len,'=') + 2
         call redlpt (enth,ibeg,iend,len,ier)

         nreact = i - 2

         if (ier.ne.0) goto 90

      else if (i.lt.3) then 
c                                 a reaction with only 2 species 
c                                 is unexpected, write error
         goto 90 

      else 

         nreact = i - 1

      end if 

      return      

90    write (*,1000) tname,(chars(i),i=1,len),name,rnum
      
      call errpau
      
1000  format ('**error ver200** READR bad data, currently ',
     *        'reading solution model: ',a,' data was:',/,400a,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

      end 

      subroutine readz (coeffs,inds,ict,idim,tname,tag)
c----------------------------------------------------------------------
c readz - read site fraction data for a solution model, assumes
c data on one line of less than 240 characters, the expected format

c        comments "=" a0 (acoef(i), mame(i), i= 1...end_of_data)

c end_of_data is either a "|" or the end of the record.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ibeg, iend, len, ier, iscan, ict, inds(k7),
     *        match, idim, i, iscnlt

      double precision rnum, coeffs(k7)

      character name*8, tname*10, tag*3

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)
c----------------------------------------------------------------------
      ict = 0 
      do i = 1, k7
         inds(i) = 0
         coeffs(i) = 0
      end do 

      call readcd (n9,len,ier,.true.)
      if (ier.ne.0) goto 90 
c                                 this first segment is only done for
c                                 the readlm routine:
      ibeg = 1
c                                 read the fist word
      call readnm (ibeg,iend,len,ier,name)
      read (name,'(a)') tag

      if (tag.eq.'end') then 
c                                 if called from readlm tag may be the
c                                 end of data marker
         return

      else 

         i = match(idim,ier,name)

         if (ier.eq.0) then 
c                                 the "comment" is a name for
c                                 readlm
            ict = ict + 1
            inds(ict) = i

         end if
 
      end if 
c                                 find start of data marker '='
      ibeg = iscan (iend,len,'=') + 1
      ict = ibeg
c                                 find a number
      call readfr (rnum,ibeg,iend,len,ier)
      if (ier.ne.0) goto 90 
c                                 find the next non-blank chracter
c                                 if it's text, then the expression
c                                 has no constant, else save the 
c                                 constant.
      if (chars(iscnlt(iend+1,lchar,'/')).lt.'A') then
c                                 assuming ascii collating sequence,
c                                 the next character is numeric, so 
c                                 the number read previously is the a0 constant
         coeffs(1) = rnum

      else 
c                                 no constant, reset character pointer.
         coeffs(1) = 0d0 
         ibeg = ict

      end if 
c                                 the rest of the data should
c                                 consist of coefficients followed
c                                 by names
      ict = 1 

      do while (ibeg.lt.len) 
c                                 find the number
         call readfr (rnum,ibeg,iend,len,ier)

         if (ier.ne.0) then 
c                                 if called from readlm may be the 
c                                 legitimate text string "delta"
            call readnm (ibeg,iend,len,ier,name)

            if (name.eq.'delta') then 

               ibeg = iscan (iend,len,'=') + 1
               call readfr (rnum,ibeg,iend,len,ier)
               if (ier.ne.0) goto 90
               coeffs(ict+1) = rnum
               exit 

            else 
c                                 invalid data 
               goto 90

            end if 

         end if

         if (ier.ne.0) goto 90
c                                 find the name
         call readnm (ibeg,iend,len,ier,name)
c                                 for constant bounds, read delta         
         if (name.eq.'delta') then 
c                                 the lower bound should be the last number read
            coeffs(ict) = rnum
c                                 next find the delta
            ibeg = iscan (iend,len,'=') + 1
            call readfr (rnum,ibeg,iend,len,ier)
            if (ier.ne.0) goto 90
            coeffs(ict+1) = rnum

            exit 

         end if 
         
         
         if (ier.ne.0) goto 90

         ict = ict + 1
         coeffs(ict) = rnum 
         inds(ict) = match(idim,ier,name)

         if (ier.ne.0) then 
            write (*,1010) name,tname,(chars(i),i=1,len)
            
            call errpau
            
         end if 
           
      end do

      return

90    write (*,1000) tname,(chars(i),i=1,len),name,rnum
      
      call errpau
      
1010  format (/,'**error ver201** invalid name: ',a,' in an expression',
     *        ' for solution model: ',a,/,' data was:',/,400a)
1000  format (/,'**error ver200** READZ bad data, currently',
     *        ' reading solution model: ',a,' data was:',/,400a,/,
     *        'last name read was: ',a,/,
     *        'last number (or real equivalent) was: ',g12.6,/)

      end 

      subroutine smodel (tname)
c-----------------------------------------------------------------
c smodel reads models for configurational entropy
c of solutions with dependent site mixing or disordered
c endmembers, see documentation section 1.3.1, equation (8b)
c-----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      integer i,j,knsp,k,l

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer nsub,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub(m10,m11,m0,m12),
     *      nttyp(m10,m11,m0)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf
c-----------------------------------------------------------------

      if (nsite.gt.m10) call error (1,a0(1,1),nsite,'m10')
c
      if (nsite.lt.isite) call error (30,a0(1,1),isite,tname)
c                                 for each site
      do i = 1, nsite 
c                                 read # of species, and site 
c                                 multiplicty.
         read (n9,*,err=90) knsp,smult(i)
         nspm1(i) = knsp - 1
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, nspm1(i)
c                                 read # of terms in the 
c                                 site fraction function and a0.
            read (n9,*,err=90) nterm(i,j), a0(i,j)
            if (nterm(i,j).gt.m0) call error (33,a0(1,1),m0,tname)
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 read term type:
               read (n9,*,err=90) nttyp(i,j,k)
               if (nttyp(i,j,k).lt.m12) then
                  read (n9,*,err=90) acoef(i,j,k),(nsub(i,j,k,l),
     *                               l = 1, nttyp(i,j,k))            
               else
                  call error (29,a0(1,1),nttyp(i,j,k),tname)
               end if
            end do 

         end do 
      end do 

      return

90    call error (20,a0(1,1),i,tname)

      end 

      double precision function gsixtr (id)
c-----------------------------------------------------------------------
c gsixtr computes G from the EoS formulated by Sixtrude & Bukowski '90, 
c JGR 95:19311-19325.

c Sixtrude parameters:

c    F0    n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, izap

      double precision a0, v0, v, df, f, dv, root, nr9t0, nr9t,
     *                 gamma0, k00, plg, c1, c2, c3, f1, gvq, 
     *                 q, vq, v23, theta0, tol, k0p, a, ethv

      double precision nr9, qm1, d2f, tht, tht0, etht, etht0, ltht, df1,
     *                 ltht0, dtht, dtht0, d2tht, d2tht0, g, g0, dg, 
     *                 dg0, d2g, d2g0, dfc, d2fc, dft, d2ft, dft0, d2ft0
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision smu
      common/ cst323 /smu

      character*8 names 
      common/ cst8   /names(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save izap
      data izap /0/
c----------------------------------------------------------------------
c                                 assign local variables:
      a0     =  thermo(1,id)
      v0     = -thermo(3,id)

      theta0 =  thermo(6,id)
      gamma0 =  thermo(7,id)
      q      =  thermo(8,id)

      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      qm1    = q - 1d0
      nr9t   = nr9*t
      nr9t0  = thermo(20,id)
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008. 
c                                 the dfth0 could be loaded as a 
c                                 constant. 
      tht    = theta0/t
      tht0   = theta0/tr
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dft   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3 
     *        - dlog(1d0-exp(-tht)))
      dft0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3 
     *        - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dft+dft0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dft-dft0) + k00)

      if (root.gt.0d0) then 
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0
      else 
         v = v0
      end if 

      f1 = 1d9
      itic = 0 
c                                 change to use relative tolerance
c                                 JADC March 1, 2005. formerly 1d-1 bar.
      tol = 1d-6*p

      do while (dabs(f1).gt.tol)

         itic = itic + 1

         vq = (v/v0)**q
         gvq = gamma0*vq
         v23 = (v0/v)**r23
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = 5d0/9d0*v23/v**2

         tht  =  theta0*dexp(-gamma0*((v/v0)**q-1d0)/q)/t
         if (tht.lt.1d-10) goto 90 
         tht0 =  tht*t/tr
         etht  = dexp(-tht )
         etht0 = dexp(-tht0)
         ltht  = dlog(1d0 - etht )
         ltht0 = dlog(1d0 - etht0)
c                                 diff(theta/T,v)
         dtht  = -gvq/v*tht
         dtht0 = -gvq/v*tht0
c                                 diff(theta/T,v,v)
         d2tht  = gvq*tht /v**2*(gvq - qm1)
         d2tht0 = gvq*tht0/v**2*(gvq - qm1)

         g   = plg(tht )
         g0  = plg(tht0)
         dg  = tht **2*ltht *dtht
         dg0 = tht0**2*ltht0*dtht0
         d2g  = ((2d0*ltht  + tht *etht /(1d0-etht ))*dtht **2 + 
     *          tht *ltht *d2tht )*tht
         d2g0 = ((2d0*ltht0 + tht0*etht0/(1d0-etht0))*dtht0**2 + 
     *          tht0*ltht0*d2tht0)*tht0

         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f

         dft  = nr9t /tht **3*(dg  -3d0/tht *g *dtht )
         dft0 = nr9t0/tht0**3*(dg0 -3d0/tht0*g0*dtht0)

         d2ft =  nr9t /tht **3*(3d0/tht *(dtht *(4d0/tht *g *dtht  
     *                        - 2d0*dg ) - g *d2tht ) + d2g )
         d2ft0 = nr9t0/tht0**3*(3d0/tht0*(dtht0*(4d0/tht0*g0*dtht0 
     *                        - 2d0*dg0) - g0*d2tht0) + d2g0)

         f1  = -dfc - dft + dft0 - p

         df1 = -d2fc - d2ft + d2ft0 

         dv = f1/df1
         v = v - dv

         if (v.le.0d0.or.v/v0.gt.2d1.or.
     *       itic.gt.100.or.dabs(f1).gt.1d40) goto 90
 
      end do
c                                 if everything is ok, now get 
c                                 helmoltz energy:
      goto 10 
c                                 if we get here, failed to converge
90    if (izap.lt.10) then
         write (*,1000) t,p,names(id)
         izap = izap + 1
         if (izap.eq.10) call warn (49,r,369,'GETLOC')
      end if 
c                                 destabilize the phase:
      gsixtr = 1d12*p

      return 
         
10    vq = (v/v0)**q
      f = 0.5d0*(v0/v)**r23 - 0.5d0
      tht  =  theta0*dexp(-gamma0*(vq-1d0)/q)/t
      tht0 =  tht*t/tr 

      a = a0 + c1*f**2*(0.5d0 + c2*f) 
     *  + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

      gsixtr = a + p*v - t*thermo(10,id)
c                                 thermal energy/v
      ethv = (dft0-dft)/gamma0/vq
c                                 etas0 = thermo(9,id)
c                                 g0 = emod(1,id)
c                                 g0p = emod(2,id)
c                                 etas = thermo(9,id)*v/v0
c                                 adiabatic shear modulus
      smu = (1d0 + 2d0*f)**(2.5d0)*
     *      (emod(1,id)*(1d0 - 5d0*f) + f*emod(2,id)*3d0*k00)
     *    -  thermo(9,id)*v/v0*ethv

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

      end 

      double precision function ghkf (id)
c-----------------------------------------------------------------------
c ghkf computes apparent G for aqueous species HKF formulation

c assumes perimittivity (epsln) and HKF g-function (gf) have been computed
c in common cxt37 (by slvnt2 or aqrxdo). 

c HKF parameters are loaded into thermo as:

c thermo(1 ,id) = G0
c thermo(2 ,id) = S0
c thermo(5 ,id) = w (omega0)
c thermo(6 ,id) = q (charge)
c thermo(7 ,id) = a1
c thermo(8 ,id) = a2
c thermo(9 ,id) = a3
c thermo(10,id) = a4
c thermo(11,id) = c1
c thermo(12,id) = c2
c thermo(13,id) = -s + c1*dlog(tr) + c1 + w*yr + dlog(tr/(tr-theta))*c2/theta**2 => b8 in HKF_G.mws
c thermo(14,id) = (-w*yr-c1+s)*tr + w - a1*pr - a2*ln(psi+pr) + g + c2/theta => b9 
c thermo(15,id) = -a3*pr-a4*ln(psi+pr) => b10
c thermo(16,id) = -c2/(tr-theta)/theta => b11
c thermo(17,id) = c2/theta^2 => b12
c thermo(18,id) = -c1-c2/theta^2
c thermo(19,id) = reference born radius 5d10*eta*q^2/(1622323167*eta*q+5d10*omega0)
c thermo(20,id) = q^2
c-----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer id

      double precision ft, fp, omega, psi, theta, z, eta

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer iam
      common/ cst4 /iam

      save psi, theta, eta
      data psi, theta, eta/2600d0, 228d0, 694656.968d0/
c-----------------------------------------------------------------------
      if (id.eq.ihy) then 
c                                 assumes proton is the only species 
c                                 with zero G0, return G_H+(P,T) = 0.
         ghkf = 0
         return 

      else if (iam.eq.5) then 
c                                 needs solvent properties if frendly
        call slvnt0 (fp,ft)

      end if 

      z = thermo(6,id)

      if (z.ne.0d0) then 
c                                 ionic species
         omega = eta * z * (z/(thermo(19,id) + dabs(z)*gf) 
     *                      - 1d0/(3.082d0 + gf))

      else 
c                                 neutral species
         omega = thermo(5,id)

      end if 

      ft = t - theta
      fp = dlog(psi+p)

      ghkf = thermo(14,id) + (thermo(13,id) + thermo(17,id)*dlog(ft) 
     *                                      + thermo(18,id)*dlog(t))*t 
     *     + thermo(16,id)*ft 
     *     + thermo(7,id)*p + thermo(8,id)*fp 
     *     + (thermo(9,id)*p + thermo(10,id)*fp + thermo(15,id))/ft
     *     + omega*(1d0/epsln - 1d0) - thermo(5,id)/epsln0

      end 

      double precision function gfunc (rho)
c----------------------------------------------------------------------
c Shock et al 1992 dielectic g function [HKF_g_function.mws]. 
c g is in angstrom. rho is the CGS solvent density. 

c hacked to use solvent density, this will, at a minimum violate the
c conditions based on p/t
c----------------------------------------------------------------------
      implicit none 

      integer iwarn

      double precision g, tf, psat2, rho

      external psat2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save iwarn
      data iwarn/0/
c---------------------------------------------------------------------
      if (rho.gt.1d0) then
c                                 region III, rho = 1 g/cm3, g = 0
         g = 0d0

      else 
c                                 region I function
         g = ((-6.557892d-6*t + 9.3295764d-3)*t
     *       -4.096745422)*((1d0 - rho)) ** 
     *       ((1.268348e-5*t - 1.767275512e-2)*t + 9.98834792)

         if (t.gt.428.15.and.p.lt.1d3) then
c
            tf = (t/300d0 - 1.427166667d0)
c                                 add region II perturbation term
            g  = g - 
     *           (tf**4.8d0 + 0.366666D-15*tf**16) 
     *         * ((((5.01799d-14*p - 5.0224d-11)*p - 1.504074d-7)*p 
     *               + 2.507672d-4)*p - 0.1003157d0)

         end if 
c                                 check on physical conditions
         if ((rho.lt.0.35d0.and.p.gt.500d0).or.
     *       (t.gt.623.15.and.p.lt.500d0).or.
     *       (t.le.623.15.and.p.lt.psat2(t))) then 
c                                 warn

            if (iwarn.lt.10) then
               write (*,1000) t, p
               iwarn = iwarn + 1
               if (iwarn.eq.10) call warn (49,r,277,'GFUNC')
            end if 
         end if 
      end if 

      gfunc = g 

1000  format (/,'**warning ver277** T= ',f8.2,' K P=',f9.1,' bar ',
     *       'is beyond the limits for',/,'the HKF g function.',/)

      end 



      double precision function epsh2o (v)
c----------------------------------------------------------------------
c Sverjenski 2014 dielectic constant for pure water, v is the molar 
c volume of water in j/bar [HKF_epsilon.mws].
c----------------------------------------------------------------------
      implicit none 

      double precision v

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      epsh2o = dexp(-0.8016651D-4 * t + 0.4769870482D1 - 0.6871618D-1 * 
     *         dsqrt(t - 0.27315D3)) * (0.1801526833D1 / v) ** 
     *         (-0.1576377D-2 * t + 0.1185462878D1 + 0.6810288D-1 * 
     *         dsqrt(t - 0.27315D3))

      end 


      double precision function duah2o (fug)
c----------------------------------------------------------------------
c Duan 2005 pure water volume of water in j/bar [HKF_duan_h2o.mws].
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision prt,b,c,d,e,f,g,expg,gamm,v,vi,veq,dveq,fug,dv

      integer it, iwarn

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5  /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save iwarn
      data iwarn/0/
c---------------------------------------------------------------------
c                                 CORK volume guess and backup fugacity
      call crkh2o (p,t,v,fug)

      prt = p/r/t
      it = 0

      gamm = 1.05999998e-02

      b = 1.957197778 - 6821674.863d0/t**2 + 3047984261d0/t**3
      c = 3.531471196 + 9821873.173d0/t**2 - 7411448875d0/t**3
      d = 16.71639581 - 6007496.747d0/t**2 + .1540316803d11/t**3
      e = -4.611555959 + 11372008.36d0/t**2 - .136192675d11/t**3
      f = -2033.267066d0 / t
      g = -0.002765323035d0 * t
c                                 iteration loop for volume
      do 

         vi = 1d0/v

         expg = dexp(-gamm/v/v)
c                                 p(v)/rt
         veq = -vi - b*vi**2 + (-f*expg-c)*vi**3 + (-g*expg-d)*vi**5 
     *             - e*vi**6
c                                 diff(veq,v)
         dveq = -veq*vi + b*vi**3 + 2d0*(f*expg+c)*vi**4  
     *           + (-2d0*f*expg*gamm + 4d0*g*expg + 4d0*d)*vi**6 
     *           + 5d0*e*vi**7 - 2d0*g*expg*gamm*vi**8

         dv = -(prt + veq)/dveq

         if (dv.lt.0d0.and.v+dv.lt.0d0) then 

            v = v*0.8d0

         else 

            v = v + dv

         end if 

         if (dabs(dv/v).lt.nopt(5)) then

            exit
          
         else if (v.lt.0d0.or.it.gt.iopt(21)) then
c                                 use cork fugacities
            iwarn = iwarn + 1

            if (iwarn.le.50) then 
               write (*,1000) p,t,v
               if (iwarn.eq.50) call warn (49,p,93,'DUAH2O')
            end if 

            exit 

         end if 

         it = it + 1

      end do 

      duah2o = v

1000  format (/,'**warning ver093** DUAH2O did not converge at:',
     *        3(1x,g12.6))

      end 

      double precision function gaq (id)
c-----------------------------------------------------------------------
c gaq computes apparent G for aqueous species with the Anderson et al. 
c (GCA 1991) density model as extended by Holland & Powell (JMG 1998)
c and modified to account for P-Pr integration.  

c parameters, 0 indicates property at Tr,Pr:

c t10 := -s0+Tr*b0+(-Tr*b0+c0)/Tr/dadt0*alpha0;
c t11 := (-Tr*b0+c0)/Tr/dadt0;
c t12 := -1/2*b0;
c t13 := Tr*s0+g0-pr*v0-1/2*b0*Tr^2+(-Tr*b0+c0)/Tr/dadt0*(-Tr*alpha0+beta0*pr);
c t14 := v0-(-Tr*b0+c0)/Tr/dadt0*beta0;

c vh2o = water volume
c vh2o0 = water volume
c alpha0 = water expansivity
c beta0 = water compressibility
c dadT0 = temperature derivative of alpha
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      double precision vh2o, vh2o0, fh2o, tp
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save vh2o0
      data vh2o0/18.7231148995863/
c----------------------------------------------------------------------
      if (thermo(21,id).eq.1d0) then 
c                                 thermo(21) flags hydronium. cause hp 
c                                 use HSC convention, gaq(H+) is not zero. 
         gaq = thermo(1,id)
         return 

      end if 
c                                 compare to last p-t, to save expensive
c                                 calls to the h2o eos.
      call pseos (vh2o,fh2o,1)
c                                 tp instead of t is a h-p innovation
      if (t.lt.500d0) then 
         tp = t 
      else 
         tp = 500d0
      end if 

      gaq = thermo(13,id) 
     *      + t*(thermo(10,id) + thermo(11,id)*dlog(vh2o0/vh2o)/tp 
     *         + thermo(12,id)*t) + thermo(14,id)*p 

      end 

      double precision function gstxlq (id)
c-----------------------------------------------------------------------
c gstxlq computes G from the liquid EoS formulated by Sixtrude et al 2009

c Stxrude parameters:

c   thermo(1)  = F0 
c   thermo(2)  = S0 - Cv  
c   thermo(3)  = V0   
c   thermo(4)  = Cv  
c   thermo(5)  = 4.5d0*K0*V0  
c   thermo(6)  = 4.5d0*K0*V0*(K'-4)    
c   thermo(7)  = y0 - y'   
c   thermo(8)  = y'    
c   thermo(9)  = T0
c   --- dependent ---
c   thermo(10) = (S0-Cv-Cv*y0)*T0
c   thermo(11) = Cv*(y0+ln(T0))-S0+Cv
c   thermo(12) = ln(v0)
  
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, izap

      logical bad

      double precision a5, a6, a7, a8, a9, a10, a11, v0, v, v2, df, f, 
     *                 d2a, v23, tol, d2f, df2, da, a1
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      character*8 names 
      common/ cst8   /names(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save izap
      data izap /0/
c----------------------------------------------------------------------
c                                 assign local variables:
      v0 = thermo(3,id)

      a10 = thermo(4,id)*(thermo(9,id)-t)*thermo(7,id)
      a7  = (thermo(11,id)-thermo(4,id)*dlog(t))*t+thermo(10,id) 
     *      - a10*thermo(12,id)
      a8  = thermo(5,id)
      a9  = thermo(6,id)
      a11 = thermo(4,id)*(thermo(9,id)-t)*thermo(8,id)/v0
      a6  = 3d0*a9
      a5  = 2d0*a8
c                                 initial guess for volume, taylor(diff(a,v),v=v0,3)
      a1 = (p+a11)*v0
      v = v0 + (9d0*(3d0*a8+a9)/(a5+a1*9d0)**2*(a1+a10) - 1d0)
     *         *9d0*v0*(a10+a1)/(a5+a1*9d0) 

      if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0

      itic = 0 
      tol = 1d-6*p

      do 

         itic = itic + 1
c                                 f, and derivatives
         v23 = (v0/v)**r23
         v2  = v**2
         f   = 0.5d0*v23 - 0.5d0
         df  = -v23/v/3d0
         df2 = df*df
         d2f = r59*v23/v2
c                                 a 1st and 2nd derivatives
c                                 da is actually diff(a,v) + p
         da  = (a5 + a6*f)*f*df + a10/v + a11 + p
         d2a = (df2 + f*d2f)*a5 + (2d0*df2 + f*d2f)*a6*f - a10/v2

         v = v - da/d2a

         if (v.le.0d0.or.itic.gt.100.or.dabs(da).gt.1d40) then 
            bad = .true.
            exit 
         else if (dabs(da).lt.tol) then
            bad = .false.
            exit 
         end if 
 
      end do

      if (bad) then  
c                                 if we get here, failed to converge
         if (izap.lt.10) then
            write (*,1000) t,p,names(id)
            izap = izap + 1
            if (izap.eq.10) call warn (49,r,369,'GSTXLQ')
         end if 
c                                 destabilize the phase.
         gstxlq  = 1d12*p

      else 
c                                 everything ok, final f:
         f = 0.5d0*(v0/v)**r23 - 0.5d0 
c                                 g = helmholtz enery + pv
         gstxlq  = (a8 + a9*f)*f**2 + a7 + a10*dlog(v) + a10 + a11*v 
     *             + p*v + thermo(1,id)

      end if  
                   
1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude Liq EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

      end 

      double precision function gstxgi (id)
c-----------------------------------------------------------------------
c gstxgi computes G from the EoS formulated by Sixtrude & Lithgow-B

c Stxrude parameters:

c    F0    -n  -v(J/bar) K0(bar) K0'   theta0 gamma0 q etaS0 Smag

c Perple_X parameters (1st index in array thermo):

c    1     2     3        4     5      6      7    8    9    10

c and for shear modulii

c    G0(bar)   G0'

c are in the Perple_X array emod(i,id) with i

c      1        2
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id, itic, izap

      logical bad

      double precision v0, v, df, f, dv, gamma0, k00, k0p,
     *           plg, c1, c2, c3, f1, aiikk, aiikk2, nr9t,
     *           root, aii, etas, a, ethv, gamma, da, nr9t0,
     *           fpoly, fpoly0, letht, letht0, z, aii2, 
     *           v23, tol, t1, t2, a2f

      double precision nr9, d2f, tht, tht0, etht, etht0, df1,
     *                 dtht, dtht0, d2tht, d2tht0, 
     *                 dfc, d2fc, dfth, d2fth, dfth0, d2fth0
 
      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision smu
      common/ cst323 /smu

      character*8 names 
      common/ cst8   /names(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save izap
      data izap /0/
c----------------------------------------------------------------------
c                                 assign local variables:
      v0     = -thermo(3,id)
      nr9    = thermo(11,id)
      c1     = thermo(12,id)
      c2     = thermo(13,id)
      c3     = thermo(14,id)
      aii    = thermo(15,id)
      aiikk  = thermo(16,id)
      aiikk2 = thermo(18,id)
      aii2   = thermo(19,id)
      nr9t0  = thermo(20,id) 

      t1     = thermo(6,id)/t
      t2     = t/tr
      nr9t   = nr9*t
c                                 initial guess for volume:
c                                 taylor(diff(FTH,v),v=v0,1)
c                                 JADC Feb 26, 2008. 
c                                 the dfth0 could be loaded as a 
c                                 constant. 
      tht    = t1
      tht0   = tht*t2
      gamma0 = thermo(7,id)
      k00    = thermo(4,id)
      k0p    = thermo(5,id)

      dfth   = nr9t*gamma0/v0*(3d0*plg(tht)/tht**3 
     *         - dlog(1d0-exp(-tht)))
      dfth0  = nr9t0*gamma0/v0*(3d0*plg(tht0)/tht0**3 
     *         - dlog(1d0-exp(-tht0)))
c                                 taylor(diff(FC,v),v=v0,2)
c     v       = (k00-dfth+dfth0-p)/k00*v0
c                                 taylor(diff(FC,v),v=v0,3)
      root = k00*((2d0+2d0*k0p)*(p+dfth-dfth0) + k00)

      if (root.gt.0d0) then 
         v = ((2d0+k0p)-dsqrt(root)/k00)*v0/(1d0+k0p)
         if (v.lt.v0/1d1.or.v.gt.v0*1d1) v = v0
      else 
         v = v0
      end if 

      itic = 0 
c                                 change to use relative tolerance
c                                 JADC March 1, 2005. formerly 1d-1 bar.
      tol = 1d-6*p

      do 

         itic = itic + 1
c                                 f, and derivatives
         v23 = (v0/v)**r23
         f = 0.5d0*v23 - 0.5d0
         df = -v23/v/3d0
         d2f = r59*v23/v**2
c                                 cold part derivatives
         dfc = (c3*f+c1)*f*df
         d2fc = (2d0*c3*f+c1)*df**2+(c3*f+c1)*f*d2f
c                                 debye T/T (tht)
         z  = 1d0+(aii+aiikk2*f)*f

         if (z.lt.0d0) then
            bad = .true.
            exit 
         end if 

         root = dsqrt(z)

         tht   = t1*root
         tht0  =  tht*t/tr
c                                 tht derivatives
         a2f   = aii2+aiikk2*f
         da    = a2f/root
         dtht  = t1*da*df
         d2tht = t1*((aiikk2/root-a2f**2/z**1.5d0)*df**2
     *               + da*d2f) 

         dtht0 = dtht*t2
         d2tht0 = d2tht*t2
c                                 polylog functions:
         fpoly   = 3d0*plg(tht )/tht**3
         fpoly0  = 3d0*plg(tht0)/tht0**3
c                                 thermal part derivatives:
         etht  = dexp(-tht )

         if (1d0-etht.lt.0d0) then
            bad = .true.
            exit 
         end if 

         letht = dlog(1d0-etht)

         dfth = (letht-fpoly)*nr9t*dtht/tht
         d2fth = ((4d0*dtht**2/tht-d2tht)*(fpoly-letht)
     *         + dtht**2*etht/(1d0-etht))*nr9t/tht

         etht0 = dexp(-tht0)

         if (1d0-etht0.lt.0d0) then
            bad = .true.
            exit 
         end if 
 
         letht0 = dlog(1d0-etht0)

         dfth0 = (letht0-fpoly0)*nr9t0*dtht0/tht0
         d2fth0 = ((4d0*dtht0**2/tht0-d2tht0)*(fpoly0-letht0)
     *          + dtht0**2*etht0/(1d0-etht0))*nr9t0/tht0

         f1  = -dfc - dfth + dfth0 - p

         df1 = -d2fc - d2fth + d2fth0 

         dv = f1/df1

c         if (dabs(dv).gt.1d-2) dv = 1d-2*dv/dabs(dv)
c                                 the above trap was probably to avoid negative 
c                                 volume, but causes problems at super-earth pressure
c                                 replaced below nov 17, 2017.
         if (v - dv.lt.0d0) dv = v/2d0

         v = v - dv

         if (itic.gt.100.or.dabs(f1).gt.1d40) then 
            bad = .true.
            exit 
         else if (dabs(f1).lt.tol) then
            bad = .false.
            exit 
         end if 
 
      end do

      if (bad) then  
c                                 if we get here, failed to converge
         if (izap.lt.10) then
            write (*,1000) t,p,names(id)
            izap = izap + 1
            if (izap.eq.10) call warn (49,r,369,'GSTX')
         end if 
c                                 destabilize the phase.
         gstxgi  = 1d12*p

      else 

c                                 everything is ok, now get 
c                                 helmoltz energy:
         f = 0.5d0*(v0/v)**r23 - 0.5d0 
         z = 1d0+(aii+aiikk2*f)*f
         root = dsqrt(z)
c                                 final estimate for tht
         tht   = t1*root
         tht0  = tht*t2
c                                 helmholtz enery
         a = thermo(1,id) + c1*f**2*(0.5d0 + c2*f) 
     *     + nr9*(t/tht**3*plg(tht ) -tr/tht0**3*plg(tht0))

         gstxgi = a + p*v - t*thermo(10,id)
c                                 z = (theta/theta0)^2
         gamma = (2d0*f+1d0)*(aii+aiikk*f)/6d0/z
         etas = - gamma - thermo(17,id)/z*(2d0*f + 1d0)**2
c                                 thermal energy/V, based on
c                                 previous v estimate
         if (gamma.ne.0d0) then
            ethv = (dfth0-dfth)/gamma
         else 
            ethv = 0d0
         end if 
c                                 adiabatic shear modulus
         smu = (1d0 + 2d0*f)**(2.5d0)*(
     *         emod(1,id) + f*(thermo(21,id) + thermo(22,id)*f))
     *       - etas*ethv

      end if  
                   
1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Sixtrude GI EoS.',
     *        ' Phase ',a,' will be destabilized.',/)

      end 

      double precision function plg (t)
c-----------------------------------------------------------------------
c evaluates debye integral: int((ln(1-exp(-t))*t^2),t=0..t)
c-----------------------------------------------------------------------
      implicit none

      integer i

      double precision t, p1, p2, p3, p4, dinc

      double precision wmach(9)
      common /ax02za/wmach
c-----------------------------------------------------------------------

      p1 = dexp(-t)
      p2 = t*t
      p3 = 2d0*t

      plg = -2.1646464674223d0

      do i = 1, 100000

         p4 = dfloat(i)
         dinc = (p2 + (p3 + 2d0/p4)/p4)*p1**i/p4/p4  
         plg = plg + dinc

         if (dinc.lt.wmach(3)) exit

      end do 

      end 

      double precision function vdpbm3 (vt,k,kprime)
c-----------------------------------------------------------------------
c vdpbm3 computes the vdp integral of a compound identified by id
c that is described by Birch-Murnaghan 3rd order EoS.
c    vt - is the volume at Pr & T
c    k  - is the bulk modulus at Pr & T
c    kprime - is -K' at Pr and T
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer itic, jerk

      double precision k, vt, rat, rat2, c0, c1, c2, 
     *                 c3, c4, c5, a0, a1, v, df, f, dv, kprime

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save jerk 
      data jerk /0/
c----------------------------------------------------------------------
c                                 constants:
      a0 = 0.375d0 * vt * k
      a1 = -0.125d0 * vt**2 * k
      c0 = (-28d0 -6d0 * kprime) * vt * a0
      c1 = (12d0 + 3d0 * kprime) * vt**2 * a0
      c2 = (16d0 + 3d0 * kprime) * a0
      c3 = a1 * vt * (-196d0 - 42d0 * kprime)
      c4 = a1 * (80d0 + 15d0 * kprime)
      c5 = a1 * vt * (108d0 + 27d0 * kprime)
c                                 use murnaghan guess for volume. GH, 6/23/16
c                                 initial guess for volume:
      dv = 1d0
      v = vt * (1d0 - kprime*p/k)**(dv/kprime)
      itic = 0 

      do while (dabs(dv).gt.1d-5)

         itic = itic + 1
         rat = (vt/v)**r13
         rat2 = rat**2
         f = p  + ((c0*v*rat+c1+c2*v**2*rat2)/v**3)
         df = (c3/rat2+c4*v/rat+c5)/v**4
         dv = f/df
         v = v - dv

         if (v.le.0d0.or.v.gt.1d6.or.itic.gt.20) then

            if (jerk.lt.10) then 
               jerk = jerk + 1
               write (*,1000) t,p
               if (jerk.eq.10) call warn (49,r,369,'VDPBM3')
            end if 
 
            vdpbm3 = 1d12*p

            return

         end if 

      end do
c                                 and the vdp integral is:
      f = 0.5d0*((vt/v)**r23-1d0)
c                                 checked in BM3_integration.mws
      vdpbm3 = p*v - vt*(pr-4.5d0*k*f**2*(1d0-f*(4d0+kprime))) 

1000  format (/,'**warning ver369** failed to converge at T= ',f8.2,' K'
     *       ,' P=',f9.1,' bar',/,'Using Birch-Murnaghan ',
     *        'EoS, probably for Ghiorso et al. MELTS/PMELTS endmember',
     *        ' data.',/,
     *        'The affected phase will be destabilized.',/)

      end 

      subroutine cartes (ksite,ids)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on site ksite of 
c solution ids (or sname). called by subdiv.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jsp, ksite, ids

      double precision ycum

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)
c----------------------------------------------------------------------
      ycum = 0d0
      jsp = ndim(ksite,ids)

      if (jsp.eq.0) then 
c                                 a relict site with only one species
c                                 left over from reform, i have know idea
c                                 if this works. i doubt it does, until
c                                 jan 8, 2016 it was setting an array (y)
c                                 that was not returned to the calling 
c                                 program. 
         prism(1) = xmn(ksite,1)
         npairs = 1
         return

      end if 

      call chopit (ycum,1d0,0,jsp,ksite,ids,0,.false.) 

      end

      subroutine blanko (text,chars,nchar,ilen)
c------------------------------------------------------------------- 
c unblnk - find the last non-blank character in text, 
c     text - character string 
c-------------------------------------------------------------------
      implicit none

      integer nchar, i, ilen

      character text*(*), chars(*)*1
 
      read (text,1000) (chars(i),i=1,ilen)
c                                 scan for blanks:
      do nchar = ilen, 1, -1
         if (chars(nchar).gt.' ') exit
      end do 

1000  format (400a)

      end

      subroutine makecp (inames,mnames,first)
c----------------------------------------------------------------------
c makecp reads the thermodynamic to compute the composition vector of 
c made phases, called by vertex. programs without composition checking
c use smakcp.

c output to console if first = .true.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k,ict, id, incomp(k0), jct, mmeos(k16*k17)

      logical inph(k16*k17), inmk(k16), eof, good, first

      double precision mcp(k16*k17,k0)
      character name*8, mnames(k16*k17)*8

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0 

      do i = 1, nmak

         inmk(i) = .true.

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do 

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j) 

         end do 

      end do 
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do 
c                                 now get the composition vectors for 
c                                 the mnames phases:
      do 

         call getphi (name,.false.,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then 

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do 

               inph(i) = .true.
               mmeos(i) = ieos

               exit
            end if 
         end do 

      end do 
c                                 find valid makes:
      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames

               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then

                  inmk(i) = .false.

                  if (first) then 
                     if (iam.ne.3.and.iam.lt.6.or.iam.eq.15) 
     *                  call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
                  end if 

                  exit

               else if (mnames(k).eq.mknam(i,j)) then

                  mkind(i,j) = k
                  meos(i) = mmeos(k)

               end if 

            end do

            if (.not.inmk(i)) exit

         end do 

      end do 
c                                 compute the composition for each 
c                                 made entitity and check if it's 
c                                 valid
      do i = 1, nmak

         if (inmk(i)) then

            name = mknam(i,mknum(i)+1)

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do 
 
            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do 
            end do 
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k) 
            end do 

            call chkphi (1,name,good)

            if (good) then 
               
               mksat(i) = .false.

            else 
c                                 no thermo componnents, sat comps?
               call chkphi (0,name,good)

               mksat(i) = .true.

               if (.not.good) then

                  inmk(i) = .false.

c                 call warn (52,tot,icmpn,mknam(i,mknum(i)+1))

                  cycle

               end if 

            end if 

         end if 
      end do 
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do 
   

      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mksat(ict) = mksat(i)

            mknum(ict) = mknum(i)

            meos(ict) = meos(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
               mkind(ict,j) = mkind(i,j)
            end do 

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do 

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do 

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do 

         end if

      end do  

      jct = 0 

      do j = 1, icmpn

         if (incomp(j).ne.0) then 
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if

      end do 

      nmak = ict

      if (nmak.gt.0.and.first) 
     *               write (*,1010) (cmpnt(incomp(j)),j=1,jct)
c                                remake list of phases required for 
c                                makes:
      inames = 0 

      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do 

            if (k.le.inames) cycle 

            inames = inames + 1
            mnames(inames) = mknam(i,j) 
   
         end do
c                                write list of valid makes:
         if (first) write (*,1000) mknam(i,mknum(i)+1),
     *                            (mcomp(i,incomp(j)),j=1,jct)

      end do 

      if (nmak.gt.0.and.first) write (*,'(/)')

1000  format (a,1x,15(f5.2,1x))
1010  format (/,'Summary of valid make definitions:',//,10x,15(a,1x)) 

      end

      subroutine sattst (ifer,good)
c----------------------------------------------------------------------
c sorts phases into the appropriate saturated phase list called by 
c input2. returns good if data is valid            
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer j,ifer,idc

      logical good

      character name*8
      common/ csta6 /name

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ic
      common/ cst42 /ic(k0)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec

      integer ifct,idfl
      common/ cst208 /ifct,idfl
c-----------------------------------------------------------------------

      good = .false.

      if (ifct.gt.0) then
c                               check for fluid species data
         do j = 1, ispec

            if (name.ne.cmpnt(idspe(j))) cycle 
            ifer = ifer + 1
            good = .true.
            call loadit (j,.false.,.true.)
            return

         end do 

      end if
 
      if (isat.gt.0) then 
c                               check for saturated composants:
c                               reject the phase if it contains
c                               a thermodynamic component:
         do j = 1, icp
            if (comp(ic(j)).ne.0d0) return
         end do
c                               now load the phase if it has
c                               the saturated component idc:
         do j = isat, 1, -1
            idc = icp + j
            if (comp(ic(idc)).ne.0d0) then
               isct(j) = isct(j) + 1 
               if (isct(j).gt.h6) call error (17,1d0,h6,'SATTST')
               iphct = iphct + 1
               if (iphct.gt.k1) call error (180,1d0,k1,
     *                            'SATTST increase parameter k1')
               ids(j,isct(j)) = iphct
               call loadit (iphct,.false.,.true.)
               good = .true.
               return
            end if
         end do 

      end if 

      end 

      double precision function gmake (id)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id, jd

      double precision g, gcpd
 
      external gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer mknum, mkind, meos
      double precision mkcoef, mdqf
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make
      common / cst335 /make(k10)
c-----------------------------------------------------------------------

      jd = make(id)

      g = 0d0
c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         g = g + mkcoef(jd,i) * gcpd (mkind(jd,i),.false.)

      end do 
c                                add the dqf correction
      gmake = g + mdqf(jd,1) + t*mdqf(jd,2) + p*mdqf(jd,3)

      end 

      subroutine smakcp (inames,mnames)
c----------------------------------------------------------------------
c smakcp reads the thermodynamic data file to compute the composition 
c vector of made phases without composition checking. vertex uses
c makecp
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer inames, i, j, k, ict, id, incomp(k0), jct 

      logical inph(k16*k17), inmk(k16), eof

      double precision mcp(k16*k17,k0)

      character name*8, mnames(k16*k17)*8

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer iam
      common/ cst4 /iam

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos
c----------------------------------------------------------------------
c                                 make a list of all definition names:
      inames = 0 

      do i = 1, nmak

         inmk(i) = .true.

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do 

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j) 

         end do 

      end do 
c                                 array to find valid makes:
      do i = 1, inames
         inph(i) = .false.
      end do 
c                                 now get the composition vectors for 
c                                 the mnames phases:
      do 

         call getphi (name,.false.,eof)

         if (eof) exit

         do i = 1, inames
            if (name.eq.mnames(i)) then 

               do j = 1, icmpn
                  mcp(i,j) = comp(j)
               end do 

               inph(i) = .true.
 
               exit
            end if 
         end do 

      end do 
c                                 find valid makes:
      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames

               if (mnames(k).eq.mknam(i,j).and.(.not.inph(k))) then 
                  inmk(i) = .false.
                  if (iam.ne.3.and.iam.lt.6.or.iam.eq.15) 
     *               call warn (51,tot,icmpn,mknam(i,mknum(i)+1))
                  exit
               else if (mnames(k).eq.mknam(i,j)) then 
                  mkind(i,j) = k
               end if 

            end do

            if (.not.inmk(i)) exit 

         end do 

      end do 
c                                 compute the composition for each 
c                                 made entitity and check if it's 
c                                 valid
      do i = 1, nmak

         if (inmk(i)) then

            do j = 1, icmpn
               mcomp(i,j) = 0d0
            end do 
 
            do j = 1, mknum(i)
               id = mkind(i,j)
               do k = 1, icmpn
                  mcomp(i,k) = mcomp(i,k) + mkcoef(i,j)*mcp(id,k)
               end do 
            end do 
c                                 test the composition vector
c                                 is it a normal phase (i.e.,
c                                 non-zero thermodynamic components)
            do k = 1, icmpn
               comp(k) = mcomp(i,k) 
            end do 

         end if 
      end do 
c                                 clean up arrays:
      ict = 0

      do i = 1, icmpn
         incomp(i) = 0
      end do 
   

      do i = 1, nmak

         if (inmk(i)) then

            ict = ict + 1

            mknum(ict) = mknum(i)

            do j = 1, mknum(ict)+1
               mknam(ict,j) = mknam(i,j)
            end do 

            do j = 1, mknum(ict)
               mkcoef(ict,j) = mkcoef(i,j)
            end do 

            do j = 1, k17
               mdqf(ict,j) = mdqf(i,j)
            end do 

            do j = 1, icmpn
               mcomp(ict,j) = mcomp(i,j)
c                                get list of used components
               if (mcomp(ict,j).ne.0d0.and.incomp(j).eq.0) incomp(j) = j

            end do 

         end if

      end do  

      jct = 0 
      do j = 1, icmpn
         if (incomp(j).ne.0) then 
            jct = jct + 1
            incomp(jct) = incomp(j)
         end if
      end do 

      nmak = ict

      if (nmak.gt.0) write (*,1010) (cmpnt(incomp(j)),j=1,jct)
c                                remake list of phases required for 
c                                makes:
      inames = 0 

      do i = 1, nmak

         do j = 1, mknum(i)

            do k = 1, inames
               if (mnames(k).eq.mknam(i,j)) exit
            end do 

            if (k.le.inames) cycle

            inames = inames + 1
            mnames(inames) = mknam(i,j) 
   
         end do
c                                write list of valid makes:
         write (*,1000) mknam(i,mknum(i)+1),(mcomp(i,incomp(j)),j=1,jct)

      end do 

      write (*,1020)

1000  format (a,1x,15(f5.2,1x))
1010  format (/,'Summary of valid make definitions:',//,10x,15(a,1x)) 
1020  format (/)

      end

c       the x-files library for adaptive minimization

c local variables:
 
c insp(jstot/jstot+1) - pointer to original index of first kstot independent
c         endmembers, followed by jstot-kstot dependent endmembers, followed
c         by the ordered species, if present (jstot+1).
c istot - total number of endmembers excluding ordered species used
c         to formulate the solution model in the input file
c jmsol(jstot,msp) - chemical species present on m'th site of jstot'th
c          endmember, original indexing.
c jstot - total number of endmembers (excluding ordered species) used
c         in the computation (i.e., excluding missing endmembers), but
c         including dependent endmembers.
c kdsol(istot/istot+nord) - endmember index, 0 if missing, -2 if dependent,
c          -1 or 0 if ordered species, reordered for missing species.
c kstot - total number of endmembers (excluding ordered species) 
c         used in the computation (i.e., excluding missing endmembers)

c global variables

c indx(i,j,k,l) - for solution i, pointer to the l'th original endmember
c                 index with species k on site j. 
c mstot(i) - istot globally
c jgsol(i,j,k) - k species indices of endmember j in solution i (jmsol globally) 
c knsp(i=1:mstot,ids) - points to the original (?) index of endmember i in ids

      subroutine reform (sname,im,first)
c---------------------------------------------------------------------
c reform - counts the number of species that can be respresented for a 
c solution given the present endmembers.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character*10 sname

      logical first

      integer kill,ikill,jkill,kill1,i,j,kosp(mst,msp),kill2,
     *        k,im,idsp,ksp(mst)

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer ostot
      common/ junk /ostot

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf
c----------------------------------------------------------------------

      if (jsmod.eq.20) then

         call reaqus (sname)
         kstot = jstot
         istot = jstot
         return

      end if 

      kill = 1

      if (first.and.isite.gt.1) call warn (50,wg(1,1),isite,sname)

      do

         do i = 1, isite
            ksp(i) = 1
            do j = 1, isp(i)
               kosp(i,j) = 0
            end do 
            do j = 1, isite
c                                 the number of endmembers that
c                                 should be present to represent
c                                 site i if all endmembers are present. 
               if (i.ne.j) ksp(i) = ksp(i)*isp(j)
            end do 
         end do
c                                 count the number of species
c                                 missing on each site
         do i = 1, istot
            if (kdsol(i).eq.0) then
               do j = 1, isite 
                  kosp(j,jmsol(i,j)) = kosp(j,jmsol(i,j)) + 1
               end do 
            end if 
         end do 
c                                 find the species that is missing
c                                 the most from the model
         kill = 99 
         kill1 = 0 

         do i = 1, isite
            if (isp(i).gt.1) then 
               do j = 1, isp(i)
c                                 idsp is the the number of species
c                                 possible - the number of missing
c                                 endmembers
                  idsp = ksp(i) -kosp(i,j)
                 
                  if (idsp.lt.kill) then 
                     ikill = i
                     jkill = j 
                     kill = idsp
                  else if (idsp.eq.kill.and.kosp(i,j).gt.0) then
c                                 kill the species that will kill the
c                                 most dependent endmembers
                     kill2 = 0 

                     do k = 1, istot
c                                 count the number that will be killed
                        if (jmsol(k,ikill).eq.jkill.and.kdsol(k).eq.-2) 
     *                     kill2 = kill2 + 1
                     end do

                     if (kill2.gt.kill1) then 
c                                 this is more than before (kill1)
                        kill1 = kill2 
                        ikill = i
                        jkill = j 
                        kill = idsp
                     end if 
                  end if 
               end do 
            end if 
         end do
c                                 kill the species jkill on site ikill
c                                 and reformulate the model (this is 
c                                 inefficient, but we don't care here). 
         call killsp (ikill,jkill)
c                                 check exit conditions
         if (istot.lt.2) then
c                                 failed, rejected too many endmembers

            im = im - 1
            if (first) call warn (25,wg(1,1),jstot,sname)
            jstot = 0

            exit 

         else if (istot.le.jstot.or.kill.eq.99) then  
c                                 changed from eq to le. JADC Nov 22, 2016
c                                 succeeded
            exit 
c                                 reorder the insp array so that the 
c                                 the first kstot endmembers are 
c                                 independent, followed by the 
c                                 dependent endmember, followed by
c                                 the ordered species. 
         end if 

      end do 
c                                 eliminate sites with only one
c                                 species
      if (isite.gt.1) call dedsit

      if (jsmod.eq.9.or.jsmod.eq.10) then 
         write (*,*) 'wonk reform'
         stop
      else
         ostot = istot
      end if 

      end 

      subroutine dedsit
c---------------------------------------------------------------------
c dedsit - eliminates chemical mixing sites with only one species
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i,j,itic,iwas(mst)

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip
c                                 local input variables
      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *      nr(j3)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      double precision yin
      common/ cst50 /yin(ms1,mst)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod
c----------------------------------------------------------------------

      itic = 0 

      do i = 1, isite

         if (isp(i).gt.1) then 
            itic = itic + 1
            iwas(itic) = i
         end if 

      end do 

      if (itic.eq.isite) return

      isite = itic
c                                 a dead site, shift the counters and limits:
      do i = 1, isite

         isp(i) = isp(iwas(i))
c                                 shift subdivision ranges
         do j = 1, isp(i) - 1
            xmn(i,j) = xmn(iwas(i),j)
            xmx(i,j) = xmx(iwas(i),j)
            xnc(i,j) = xnc(iwas(i),j)
            imd(j,i) = imd(j,iwas(i))
            if (imd(j,i).gt.0) yin(j,i) = yin(j,iwas(i))
         end do

      end do 

      do i = 1, istot + norder 
c                                 reset the species pointers (jmsol) 
         do j = 1, isite

            jmsol(i,j) = jmsol(i,iwas(j))

         end do

      end do  

      if (isite.eq.1) recip = .false.

      if (order) then 

         if (isite.eq.1) jsmod = 6

      else if (depend) then 

         jsmod = 7

      else 

         jsmod = 2

      end if 

      end

      subroutine killsp (ikill,jkill)
c---------------------------------------------------------------------
c killsp - eliminates species jkill from site ikill in a solution model
c and reformulates the model accordingly
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      logical skip, bad, dead

      integer jsp,jtic,morder,jst,jend,jnc,kst,ldep,ltic,
     *        i,j,ikill,jkill,kill,kdep,jdqf,ktic,jold,
     *        iwas(m4),i2ni(m4),kwas(m4),
     *        k,l,itic,ijkill(m4),
     *        j2oj(msp),j2nj(msp),i2oi(m4),maxord,mord
c                                 dqf variables
      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip
c                                 local input variables
      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer ostot
      common/ junk /ostot

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer nsub,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub(m10,m11,m0,m12),
     *      nttyp(m10,m11,m0)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision yin
      common/ cst50 /yin(ms1,mst)

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)
c----------------------------------------------------------------------
      if ((jsmod.eq.9.or.jsmod.eq.10).and.ikill.gt.isite) then
         jst = ikill
         jend = ikill
         jnc = 1
         kst = istot + 1
         ldep = ostot
      else 
         jst = 1
         jend = isite
         jnc = 0
         kst = 1
         ldep = istot
      end if 
c                                 was 1,.site
      do i = jst, isite + jnc 

         if (i.ne.ikill) then
c                                 nothing happens  
            do j = 1, isp(i)
              j2oj(j) = j
            end do

         else
c                                 on a site where we kill a species
            jsp = isp(i) - 1
c                                 should also check and store subdivsion
c                                 model here (i.e., some ternary model
c                                 choices are not valid as binary):

c                                 make a pointer from the new species
c                                 index to the old index
            jtic = 0 

            do j = 1, isp(i)

               if (j.ne.jkill) then 
                  jtic = jtic + 1 
c                              pointer from new j to old j
                  j2oj(jtic) = j
c                              pointer from old j to new j
                  j2nj(j) = jtic
               end if

            end do
c                              now reload
            isp(i) = jsp

            if (jsp.gt.1-jnc) then   
c                              now shift subdivision ranges
               do j = 1, jsp + jnc - 1
                  xmn(i,j) = xmn(i,j2oj(j))
                  xmx(i,j) = xmx(i,j2oj(j))
                  xnc(i,j) = xnc(i,j2oj(j))
                  imd(j,i) = imd(j2oj(j),i)

                  if (imd(j,i).gt.0) yin(j,i) = yin(j2oj(j),i)

               end do

            else 

               xmn(i,j) = 1d0
               xmx(i,j) = 1d0
               xnc(i,j) = 1d0

            end if
         end if 
      end do

      kdep = 0

      if (depend) then
c                                create an array which gives the
c                                original locations of the dependent
c                                endmembers, need this to be able to
c                                reorder the y2p array:
         do i = 1, istot
            if (kdsol(i).eq.-2) then 
               kdep = kdep + 1
               iwas(kdep) = i 
            end if
         end do

      end if

      do i = kst, ostot
c                                 kill endmembers with the species
c                                 to be deleted:
         if (jmsol(i,ikill).eq.jkill) kdsol(i) = -3
      end do
c                                 check for dependent endmembers
      call redep (-3)
c                                 at this point all ordered endmembers to 
c                                 be killed are flagged by kdsol = -3.

c                                 now check the ordered species
      morder = 0 

      if (order) then 
c                                 first check if the ordered endmember
c                                 may be stable 
         do k = 1, norder
c                                 check if a missing constituent 
            bad = .false.

            do j = 1, nr(k)
               if (kdsol(iddeps(j,k)).eq.-3) then 
                  bad = .true.
                  exit 
               end if 
            end do

            if (bad) then 
c                                 add species to the kill list
               kdsol(ostot+k) = -3
            
            else 

               morder = morder + 1
               kdsol(ostot+k) = -1
               kwas(morder) = k

            end if 

         end do  

      end if 
c                                figure out which dependent endmembers have
c                                been killed:
      if (depend) then 

         if (kdep.ne.mdep) call error (54,dqf(1,1),mdep,'KILLSP')

         mdep = 0 

         do i = 1, kdep
            if (kdsol(iwas(i)).ne.-3) then 
               mdep = mdep + 1
c                                 iwas is now the original index of the 
c                                 dependent endmember, and mdep is the reset
c                                 value of the dependent endmember counter
               iwas(mdep) = i                
            end if
         end do 

      end if 
 
      itic = 0
      jtic = 0
      ktic = 0 
      ltic = 0 
      kill = 0
      
      do i = 1, ostot + norder 

         if (kdsol(i).ge.-2) then 
c                                 replacement for istot (itic)
            itic = itic + 1
            if (i.le.istot) then 
               ktic = ktic + 1
            else if (i.le.ostot) then 
               ltic = ltic + 1 
            end if 
c                                 pointers from new to old endmember index (i2oi)
            i2oi(itic) = i
c                                 pointers from new to old endmember index (i2ni)
            i2ni(i) = itic                                
c                                 pointer to original species index
            iorig(itic) = iorig(i)
c                                 number of missing endmembers (jtic)
            if (kdsol(i).eq.0) jtic = jtic + 1
c                                reset the kdsol array
            kdsol(itic) = kdsol(i)

         else 
c                                 kill records the killed endmembers
            kill = kill + 1
            ijkill(kill) = i

         end if

      end do
c                                reset total and present counters
      istot = ktic

      jstot = istot - jtic

      do i = kst, itic

         if (i2oi(i).gt.ostot) cycle
c                                 reset the species pointers (jmsol)
         do j = jst, jend
            if (j.eq.ikill) then
               jmsol(i,j) = j2nj(jmsol(i2oi(i),j))
            else
               jmsol(i,j) = jmsol(i2oi(i),j)
            end if 
         end do
      end do 

      ostot = ktic + ltic
c                                --------------------------------------
c                                excess terms:
      itic = 0 
      maxord = 0 

      do i = 1, iterm 
c                                check for forbidden terms (i.e., terms
c                                with a missing endmember
         skip = .false. 
c                                 macroscopic formulation
         do j = 1, kill
c                                 check if subscript points to a killed 
c                                 endmember
            do k = 1, iord
               if (isub(i,k,1).eq.0) then
                  cycle 
               else if (isub(i,k,1).eq.ijkill(j)) then
                  skip = .true.
                  exit 
               end if 
            end do 

            if (skip) exit

         end do 

         if (skip) cycle 
c                               the term is acceptable
         itic = itic + 1

         mord = iord

         do j = 1, iord
            if (isub(i,j,1).eq.0) then
               isub(itic,j,1) = 0 
            else 
               isub(itic,j,1) = i2ni(isub(i,j,1))
            end if 
         end do 

         if (xtyp.eq.0) then 
c                                save the coefficient
            do j = 1, m3
               wg(itic,j) = wg(i,j)
            end do 
c                                find highest order term
            if (mord.gt.maxord) maxord = mord

         else
c                                 redlich kistler
            rkord(itic) = rkord(i)

            do j = 1, rkord(itic)
               do k = 1, m16
                  wk(k,j,itic) = wk(k,j,i)
               end do
            end do 

            maxord = 2

         end if 

      end do     
c                                reset counters, iord is not reset
      iterm = itic
      iord = maxord                         
c                                --------------------------------------
c                                van laar volume functions
      if (laar) then
         do i = kst, ostot + morder 
            do j = 1, m3
               vlaar(j,i) = vlaar(j,i2oi(i))
            end do 
         end do  
      end if 
c                                 --------------------------------------
c                                 dqf corrections, this is sloppy since
c                                 uses istot instead of kstot
      if (idqf.gt.0) then 

         jdqf = 0 
c                                 check if a retained species has a dqf
c                                 correction
         do j = 1, idqf
c                                 the itoi index must be in the inner loop
c                                 in case the values of indq are not sequential
            do i = 1, ostot
               if (indq(j).eq.i2oi(i)) then 
c                                 found a dqf'd endmember
                  jdqf = jdqf + 1
                  indq(jdqf) = i
                  do k = 1, m3
                     dqf(k,jdqf) = dqf(k,j)
                  end do 
                  exit 
               end if
            end do 

            if (jdqf.eq.idqf) exit
 
         end do 

         idqf = jdqf 

      end if 
c                                 --------------------------------------
c                                 configurational entropy model

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, nsite 
c                                 for each species, read function to define 
c                                 the site fraction of the species and eliminate
c                                 killed species

c                                 species counter is incremented in advance
c                                 and must be decremented before saving the 
c                                 final value:
         jtic = 1 

         do j = 1, nspm1(i)

            ktic = 0 
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 macroscopic formulation:
c                                 note: 4th index (nttyp) is only used
c                                 for bragg-williams models.
               dead = .false.
               do l = 1, kill 
                  if (nsub(i,j,k,1).eq.ijkill(l)) then 
                     dead = .true.
                     exit 
                  end if 
               end do 

               if (.not.dead) then
c                                 the term has survived (and therefore 
c                                 also the species):
c                                 don't save nttyp since this is always
c                                 1 for the macroscopic formulation
                  ktic = ktic + 1
c                                 but my dear peanut brained friend, do 
c                                 not forget to move the pointer:
                  nsub(i,jtic,ktic,1) = i2ni(nsub(i,j,k,1)) 
                  acoef(i,jtic,ktic) = acoef(i,j,k)
               end if  
            end do 
c                                 ktic is the number of terms representing
c                                 the jth species, we won't count species 
c                                 with no terms because the endmember configurational
c                                 entropy is assumed to be implicit. 
         if (ktic.gt.0) then 
c                                 increment the species counter
            nterm(i,jtic) = ktic 
            a0(i,jtic) = a0(i,j)
            jtic = jtic + 1
         end if 

      end do

      nspm1(i) = jtic - 1

      end do 
c                                 ---------------------------------------
c                                 ordered species:
      if (order) then 

         norder = morder

         if (morder.eq.0) then 
c                                 there are no ordered species left
            order = .false.

            if (depend) then

               jsmod = 7

            else if (jsmod.eq.27) then 
c                                 special case, green et al 2016 melt model
c                                 converts to normal HP melt model
               jsmod = 24

            else 
c                                 why jsmod = 2?
               jsmod = 2

            end if 

         else 
c                                 shift the ordered species pointers 
c                                 and data to eliminate kill ordered
c                                 species.    
            do j = 1, morder 

               jold = kwas(j)

               do i = 1, 3
                  denth(j,i) = denth(jold,i)
               end do 

               nr(j) = nr(jold)

               do i = 1, nr(j)
                  iddeps(i,j) = i2ni(iddeps(i,jold))
                  depnu(i,j) = depnu(i,jold)
               end do

               itic = 1 

               do i = 1, limn(jold)
c                                 eliminate absent species from 
c                                 stoichiometric p0 limits
                  ktic = 0 

                  if (limt(i,jold).gt.0) then 
   
                     do k = 1, limt(i,jold)

                        skip = .false.

                        do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                           if (limid(k,i,jold).eq.ijkill(l)) then
                              skip = .true.
                              exit 
                           end if
  
                        end do  

                        if (skip) cycle

                        ktic = ktic + 1

                        limid(ktic,itic,j) = i2ni(limid(k,i,jold))
                        limc(ktic,itic,j) = limc(k,i,jold)

                     end do 

                     if (ktic.eq.0) cycle

                     limt(itic,j) = ktic

                  else 
c                                 constant bounds
                     limt(itic,j) = -1
                     k = 1

                  end if   

                  limc(ktic+1,itic,j) = limc(k,i,jold)
                  limc(ktic+2,itic,j) = limc(k+1,i,jold)
 
c                                 now check the p terms, this assumes
c                                 there are no p terms if there are no
c                                 p0 terms, which maynot be true?
                  ktic = 0 

                  do k = 1, jimt(i,jold)

                     skip = .false.

                     do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                        if (jimid(k,i,jold).eq.ijkill(l)) then
                           skip = .true.
                           exit 
                        end if
  
                     end do  

                     if (skip) cycle

                     ktic = ktic + 1

                     jimid(ktic,itic,j) = i2ni(jimid(k,i,jold))
                     jimc(ktic,itic,j) = jimc(k,i,jold)

                  end do 

                  jimt(itic,j) = ktic

                  itic = itic + 1

               end do 

               limn(j) = itic - 1

            end do 

         end if  

      end if 
c                                 --------------------------------------
c                                 dependent endmember properties, the
      if (depend) then 
c                                 dependent endmembers have been reordered
c                                 in redep, but are still expressed in
c                                 terms of the old indices, so reset the 
c                                 indices:
         do i = 1, mdep
            jdep(i) = i2ni(jdep(i))
            do j = 1, ndph(i)
               idep(i,j) = i2ni(idep(i,j))
            end do 
         end do 

      end if

      if (istot.eq.ostot) then 
         if (jsmod.eq.9) then 
            jsmod = 8
         else if (jsmod.eq.10) then 
            jsmod = 7 
         end if
      end if 

      end

      subroutine kill01 (site)
c---------------------------------------------------------------------
c reformulates orphan vertex model (jsmod = 9) for missing orphans.
c site is the simplex with the orphan vertices.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ksp, site

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)
c----------------------------------------------------------------------
c                                 count the number of species
c                                 missing on site
      ksp = 0

      do

         do i = 1, isp(site)
            if (kdsol(istot+i).eq.0) then 
               ksp = ksp + 1
               call killsp (site,i)
               exit 
            end if 
         end do

         if (i.gt.isp(site)) exit

      end do

      end

      subroutine reaqus (tname)
c---------------------------------------------------------------------
c reaqus - reformulates aqueous solution models to eliminate missing 
c endmembers and shift ranges.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character*10 tname

      integer i, jq, jn, js, lm

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------

      jq = 0
      jn = 0
      js = 0 
      lm = 0 

      do i = 1, ns

         if (kdsol(i).ne.0) then
            js = js + 1
            iorig(js) = i 
            kdsol(js) = kdsol(i)
            lm = lm + 1
            xmn(1,lm) = xmn(1,i)
            xmx(1,lm) = xmx(1,i)
            xnc(1,lm) = xnc(1,i)
            imd(lm,1) = imd(i,1)
         end if

      end do

      do i = ns+1, ns+nn

         if (kdsol(i).ne.0) then
            jn = jn + 1
            iorig(js+jn) = i 
            kdsol(js+jn) = kdsol(i)
            lm = lm + 1
            xmn(1,lm) = xmn(1,i)
            xmx(1,lm) = xmx(1,i)
            xnc(1,lm) = xnc(1,i)
            imd(lm,1) = imd(i,1)
         end if

      end do

      do i = ns+nn+1, ns+nn+nq

         if (kdsol(i).ne.0) then
            jq = jq + 1
            iorig(js+jn+jq) = i 
            kdsol(js+jn+jq) = kdsol(i)

            if (i.eq.ns+nn+nq) cycle 
            lm = lm + 1
            xmn(1,lm) = xmn(1,i)
            xmx(1,lm) = xmx(1,i)
            xnc(1,lm) = xnc(1,i)
            imd(lm,1) = imd(i,1)

         end if

      end do

      ns = js
      nn = jn 
      nq = jq

      if (ns.eq.0) then 

         call warn (99,0d0,0,'rejecting '//tname//' because no solvent'
     *             //' species were identified')
         jstot = 0
         return

      else if (nq.eq.1) then

         call warn (99,0d0,0,'eliminating ions from '//tname//' because'
     *             //' only one charged species was identified')
         nq = 0

      end if

      jstot = ns + nn + nq

      end

      subroutine cmodel (im,idsol,tname,first)
c---------------------------------------------------------------------
c cmodel - checks to see if solution models contain valid endmembers.
c modified to allow saturated phase/component endmembers, 10/25/05.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      logical first, ok
 
      character*10 tname, missin(m4)*8

      integer imiss, im, idsol, i, j, h, ineg, ipos
 
      integer isoct
      common/ cst79 /isoct

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer eos
      common/ cst303 /eos(k10)

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character*8 names
      common/ cst8 /names(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character mname*8
      common/ cst18a /mname(m4)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer ikp
      common/ cst61 /ikp(k1)

      integer iam
      common/ cst4 /iam

      integer ostot
      common/ junk /ostot

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
      jstot = 0
      ineg = 0 
      ipos = 0 
      ok = .false.
c                              if called by build (iam = 4) skip the
c                              name check:
      if (iam.ne.4) then
c                              check that solution name matches a 
c                              name requested in input from n1
         do i = 1, isoct

           if (tname.eq.fname(i)) then 
c                              got a match, exit
               idsol = i 
               im = im + 1
               ok = .true.
               exit 

            end if 
c
         end do 
c                                 didn't find a match, read a new name:
         if (.not.ok) return 

      end if  

      do i = 1, ostot

         kdsol(i) = 0
         ok = .false.

         if (jsmod.ge.7.and.jsmod.le.9) then
c                              solution with dependent endmembers, if endmember i
c                              is dependent endmember flag it by setting kdsol(i) = -2
            do j = 1, mdep

               if (jdep(j).eq.i) then
                  kdsol(i) = -2
                  ok = .true.
                  exit 
               end if 

            end do

            if (ok) cycle 

         end if 

         if (jsmod.eq.20.and.i.gt.ns) then
c                                 aqueous solute, test against aqnam 
            do h = 1, aqct 

               if (aqnam(h).eq.mname(i)) then
c                                 got a valid endmember, count and                     
                  jstot = jstot + 1 
                  kstot = jstot
c                                 create arrays of convenience, where j = 1, jstot
                  kdsol(i) = h + aqst                                 
c                                 tests for aq solvent models
                  if (i.gt.nn+ns) then 
c                                 must be a charged solute species
                     if (thermo(6,h+aqst).eq.0) then
                        write (*,*) aqnam(h),
     *                  ' is not a charged species '
     *                  ,'remove it from the list of charged species in'
     *                  ,'solution model:',tname
                        call errpau
                        stop
                     else if (thermo(6,h+aqst).gt.0) then
                        ipos = ipos + 1
                     else 
                        ineg = ineg + 1
                     end if

                  else if (i.gt.nq) then 
c                                  must be neutral species
                     if (thermo(6,h+aqst).ne.0) then
                        write (*,*) names(h),' is a charged species'
     *                            ,'remove it from the list of neutral'
     *                            ,' species in solution model:',tname
                        call errpau
                        stop
                     end if

                  end if 

                  exit 

               end if 

            end do 

         else
c                                 must be a real enitity:
            do h = 1, ipoint

               if (names(h).eq.mname(i)) then
c                                 got a valid endmember, count and                     
                  jstot = jstot + 1 
                  kstot = jstot
c                                 create arrays of convenience, where j = 1, jstot
                  kdsol(i) = h                                 
c                                 set kill flag: 
c                                 don't reset ikp if it has been set for 
c                                 another model.
                  if (iend(i).ne.0.and.ikp(h).eq.0) ikp(h) = -1

                  if (eos(h).eq.15.or.eos(h).eq.16) then 

                     write (*,'(a,/)') names(h)//' is not described '//
     *                               'by a solvent EoS remove it from',
     *                               'the list of solvent species in '//
     *                               'solution model:'//tname

                     call errpau

                     stop

                  end if 

                  exit 

               end if 

            end do
            
         end if 
c                                 found all endmembers:
         if (jstot.eq.ostot) exit

      end do

      if (jsmod.eq.20) then 
c                                 for solvent models, check that charge balance
c                                 is possible
         if ((ipos.eq.0.and.ineg.gt.0).or.
     *       (ipos.gt.0.and.ineg.eq.0)) then 
 
             write (*,'(/,a,/,a)') 
     *                         'models with charged species must incl'//
     *                         'ude both postive and negative species',
     *                         'solution model: '//tname//
     *                         ' will be rejected'
  
             jstot = 0 

         end if

      end if 

      call redep (0)
c                                done if nothing is missing:
      if (jstot.eq.ostot) return
c                                missing endmember warnings:
      if (jstot.lt.2) then

         im = im - 1
         if (first) call warn (25,wg(1,1),jstot,tname)
         jstot = 0 

      else 

         imiss = 0

         do i = 1, ostot
            if (kdsol(i).eq.0) then
               imiss = imiss + 1
               missin(imiss) = mname(i)
            end if 
         end do 

         if (first) write (*,1000) tname,(missin(i), i = 1, imiss)

      end if 

1000  format (/,'**warning ver114** the following endmembers',
     *          ' are missing for ',a,/,4(8(2x,a),/),/)

      end

      subroutine redep (jkill)
c----------------------------------------------------------------------
c redep does reordering of dependent endmembers

c jkill is the value of kdsol that indicates invalid (missing) endmembers,
c this is 0 from cmodel and -3 from reform.
c----------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i,j,l,ndep,k,jkill

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip
 
      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer ostot
      common/ junk /ostot
c----------------------------------------------------------------------
c                                check for dependent endmembers, necessary?
      if (depend) then

         ndep = 0  

         do 100 i = 1, mdep

            do j = 1, ndph(i)

               if (idep(i,j).le.ostot) then
c                                looking for a dependent endmember component: 
                  if (kdsol(idep(i,j)).eq.jkill.and.
     *                kdsol(jdep(i)).ne.-3) then
c                                a component is missing in cmodel, reset kdsol to 0
                     kdsol(jdep(i)) = 0
                     goto 100 
                  end if 

               else
c                                looking for an ordered species
                  do l = 1, norder 
                     do k = 1, nr(l)
                        if (kdsol(iddeps(k,l)).eq.jkill) then
                           kdsol(jdep(i)) = 0 
                           goto 100
                        end if 
                     end do 
                  end do 

               end if 
            end do 
c                                dependent endmember is ok
            ndep = ndep + 1
            jdep(ndep) = jdep(i)
            ndph(ndep) = ndph(i)

            do j = 1, ndph(i)
               idep(ndep,j) = idep(i,j)
               nu(ndep,j) = nu(i,j)
            end do 

            jstot = jstot + 1

100      continue 

         mdep = ndep

         if (mdep.eq.0) depend = .false.

      end if 

      end 

      subroutine rmodel (tname,tn1,tn2,bad)
c---------------------------------------------------------------------
c rmodel - reads solution models from LUN n9.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'
 
      character*10 tname, tn1*6, tn2*22, tag*3, char*1, key*22, val*3, 
     *          nval1*12, nval2*12, nval3*12, strg*40, strg1*40, 
     *          estrg*80

      integer nreact,i,j,k,l,m,jlaar,idim

      logical bad

      double precision coeffs(k7), rnums(20), enth(3)

      integer ijk(mst),inds(k7),ict

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer ostot
      common/ junk /ostot

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer nsub,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub(m10,m11,m0,m12),
     *      nttyp(m10,m11,m0)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision yin
      common/ cst50 /yin(ms1,mst)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
      mdep = 0 
      norder = 0 
      istot = 0
      ostot = 0 
      ist(1) = 0 
c DEBUG DEBUG
      do i = 1, m4
         kdsol(i) = 0
      end do 
c                               set logical variables
      macro = .false.
      order = .false.
      laar = .false.
      depend = .false.
      fluid = .false.
      recip = .false.
      bad = .false.

      do 
          read (n9, '(3(a,1x))', iostat = i) tname
          if (i.ne.0) return
          read (tname,'(a1)') char
          if (char.eq.' '.or.char.eq.'|') cycle
          if (tname.ne.'begin_mode') exit
      end do 
c                             initialize strings
      tn1 = tname
      tn2 = 'unclassified'
c                             next look for optional abbreviation and full_name
      do 
          call redcd1 (n9,i,key,val,nval1,nval2,nval3,strg,strg1)

          read (key,'(i3)', iostat = i) jsmod

          if (i.eq.0) then 
c                             read jsmod, model type flag
             exit

          else if (key.eq.'abbreviation') then

             tn1 = strg1

          else if (key.eq.'full_name') then 

             tn2 = strg1

          else

             write (estrg,'(4(a))') 'unrecognized text: ', key,
     *                              ' reading solution model ', tname
             call error (72, enth(1), i, estrg)

          end if  
                   
      end do

      if (jsmod.eq.20) then
c                                 aqueous model reads to different 
c                                 arrays
         call raqmod (tname)
         istot = nq + nn + ns
         ostot = istot

         return

       end if
c                                 check for disabled model types
      if (jsmod.eq.1.or.jsmod.eq.3.or.jsmod.eq.5) 
     *                        call error (68,enth(1),jsmod,tname)

      if (jsmod.eq.0) fluid = .true.
      if (jsmod.eq.6.or.jsmod.eq.8.or.jsmod.eq.9.or.
     *                                jsmod.eq.27) order = .true.
      if (jsmod.ge.7.and.jsmod.le.10) depend = .true.
      if (jsmod.ge.7.and.jsmod.le.10) recip = .true.
c                                 assign non-default props to 
c                                 special models:
      if (jsmod.ge.30.and.jsmod.le.31) recip = .true.
c                                 read number of independent sites:
      if (recip) then 
         call readda (rnums,1,tname)    
         isite = idint(rnums(1)) 
      else
         isite = 1
      end if 
c                               read number of species on each site:
      call readda (rnums,isite,tname)  
      do i = 1, isite
         isp(i) = idint(rnums(i))
      end do 
c                               total number of endmembers:
      istot = 1
      do i = 1, isite
         istot = istot*isp(i)
      end do

      if (jsmod.eq.9.or.jsmod.eq.10) then 
c                               read the number of orphan site endmembers
         call readda (rnums,1,tname)
         isp(isite+1) = idint(rnums(1))
         ostot = istot + isp(isite+1)

      else 

         ostot = istot 

      end if
c                               counter for character read routines
c                               and starting index
      idim = ostot

      if (ostot.gt.m4) call error (1,rnums(1),idim,
     *                 'm4 (maximum number of endmembers)')

      i = 0
c                               read endmember names:
      call readn (i,idim,tname)
c                               compound formation models
      if (order) then 
c                               get the number of ordered species
         call readda (rnums,1,tname)    
         norder = idint(rnums(1))

         if (ostot+norder.gt.m4) call error (1,rnums(1),istot+norder,
     *                           'm4 (maximum number of endmembers)')

         if (norder.gt.j3) call error (5,rnums(1),norder,tname)
c                               get ordering reaction and name
c                               of ordered species:   
         do i = 1, norder   

            nreact = -1

            call readr (coeffs,enth,inds,idim,nreact,tname)

            do j = 1, 3
               denth(i,j) = enth(j)
            end do 

            nr(i) = nreact
  
            do j = 1, nreact
               depnu(j,i) = coeffs(j+1)
               iddeps(j,i) = inds(j+1)
            end do

         end do  
c                               read the limit equations for the 
c                               amount of the ordered endmembers
         call readlm (idim,tname,bad)

      end if 
c                               read dependent endmembers
      if (depend) then
c                               number of dependent endmembers
         call readda (rnums,1,tname)    
         mdep = idint(rnums(1))
         if (mdep.gt.m15) call error (1,enth(1),mdep,'m15')

         do i = 1, mdep
c                               nreact is returned by readr
            nreact = 0 

            call readr (coeffs,enth,inds,idim,nreact,tname)

            jdep(i) = inds(1)
            ndph(i) = nreact - 1
            if (ndph(i).gt.j4) call error (1,enth(1),ndph(i),'j4')

            do j = 1, ndph(i)
               nu(i,j) = coeffs(j+1)
               idep(i,j) = inds(j+1)
            end do 

         end do

      end if
c                               read endmember flags:
      call readda (rnums,ostot,tname)  

      do i = 1, ostot
         iend(i) = idint(rnums(i))
      end do 
c                               read composition limits, subdivision type:
      m = 0

      do i = 1, isite
c                               number of ranges to be read
         m = m + isp(i) - 1
      end do
c                               get the numbers
      do i = 1, isite
         
         do j = 1, isp(i) - 1

            call readda (rnums,4,tname)

            xmn(i,j) = rnums(1)
            xmx(i,j) = rnums(2)
            xnc(i,j) = rnums(3)
            imd(j,i) = idint(rnums(4))

         end do 

      end do
c                                read the orphan vertex subdivision data 
      if (jsmod.eq.9.or.jsmod.eq.10) then

         i = 3

         do j = 1, isp(i)

            call readda (rnums,4,tname)

            xmn(i,j) = rnums(1)
            xmx(i,j) = rnums(2)
            xnc(i,j) = rnums(3)
            imd(j,i) = idint(rnums(4))

         end do

      end if
c                                create bragg-williams indexes
      do i = 2, isite
         ijk(i) = 1
      end do 
 
      ijk(1) = 0  

      do 20 l = 1, istot

         do m = 1, isite

            if (ijk(m).lt.isp(m)) then 

               ijk(m) = ijk(m) + 1
c                                increment only one index per endmember
               do i = 1, isite
                  jmsol(l,i) = ijk(i)
               end do 

               if (ijk(m).eq.isp(m).and.l.lt.istot-isp(1)+1) then
c                                saturated counter, increment first
c                                unsaturated counter and reset all
c                                lower site counters to first index,
                  do j = 2, isite
                     if (ijk(j).lt.isp(j)) then
                        ijk(j) = ijk(j) + 1
                        do k = 2, j-1
                           ijk(k) = 1
                        end do 
                        ijk(1) = 0
                        goto 20 
                     end if
                  end do  
               end if 

               goto 20

            end if
         end do 
20    continue
c                              orphan vertex models
      do l = istot + 1, ostot
         jmsol(l,isite+1) = l - istot 
      end do
c                              read excess function
      call readx (idim,tname)
c                              expansion for S(configurational)
      call readda (rnums,1,tname)    

      nsite = idint(rnums(1))

      if (nsite.gt.m10) call error (1,a0(1,1),nsite,'m10')
c                                 for each site
      do i = 1, nsite 
c                                 read # of species, and site 
c                                 multiplicty.
         call readda (rnums,2,tname) 

         smult(i) = rnums(2)
c                                 if multiplicity is 0, the model 
c                                 site has variable multiplicity 
c                                 and molar site population expressions
c                                 are read rather than site fractions
c                                 in which case we need as many expressions
c                                 as species. nspm1 is the counter for the
c                                 number of expressions
         if (smult(i).gt.0) then 
            nspm1(i) = idint(rnums(1)) - 1
         else 
            nspm1(i) = idint(rnums(1))
         end if 

         if (nspm1(i).gt.m11) call error (1,a0(1,1),nspm1(i),'m11')
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
c
         do j = 1, nspm1(i)
c                                 read expression for site
c                                 fraction of species j on 
c                                 site i.
            call readz (coeffs,inds,ict,idim,tname,tag)

            a0(i,j) = coeffs(1)
            nterm(i,j) = ict - 1
            if (nterm(i,j).gt.m0) call error (33,a0(1,1),m0,tname)
c                                 for each term:
            do k = 2, ict
c                                 all terms 1 order type, this
c                                 saved for compatability with 
c                                 old versions:
               nttyp(i,j,k-1)   = 1
               acoef(i,j,k-1)   = coeffs(k)
               nsub(i,j,k-1,1) = inds(k)
            end do 
         end do 
      end do 
c                              look for van laar and/or dqf parameters
c                              or the end of model marker
      call readop (idim,jlaar,ostot-mdep,reach,stck,norf,tname)

      if (jlaar.ne.0) then

         laar = .true.
c                                 high order terms not allowed for
c                                 van laar.
         if (iord.gt.2.and.laar) call error (999,coeffs(1),800,'RMODEL')

      end if 
c                                 save original indices, need this for 
c                                 melt models etc that have species specific
c                                 equations of state.
      do i = 1, ostot + norder
         iorig(i) = i 
      end do

      end

      subroutine raqmod (tname)
c---------------------------------------------------------------------
c raqmod - reads aq electrolyte solution model, jsmod = 20.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character tname*10

      integer i, j

      double precision rnums(10)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer nq,nn,ns
      common/ cxt337 /nq,nn,ns
c----------------------------------------------------------------------
c                               read number of solvent species:
      call readda (rnums,1,tname)
      ns = idint(rnums(1))
c                               read solvent species names:
      i = 0
      if (ns.gt.0) call readn (i,ns,tname)
c                               read number of neutral species:
      call readda (rnums,1,tname)
      nn = idint(rnums(1))
c                               read neutral species names:
      i = ns
      if (nn.gt.0) call readn (i,ns+nn,tname)
c                               read number of charged species:
      call readda (rnums,1,tname)
      nq = idint(rnums(1))
c                               read charged species names:
      i = nn + ns
      if (nq.gt.0) call readn (i,i+nq,tname)
c                               read composition limits, subdivision type
c                               for (nn + ns + nq) - 2 species
      if (i+nq.eq.2) i = i + 1
c
      do j = 1, i + nq - 1
c                               dummy for the ns'th species
         if (j.eq.ns) then 

            xmn(1,j) = 0d0
            xmx(1,j) = 1d0
            cycle

         end if 

         call readda (rnums,4,tname)

         xmn(1,j) = rnums(1)
         xmx(1,j) = rnums(2)
         xnc(1,j) = rnums(3)
         imd(j,1) = idint(rnums(4))
c                                 don't allow imod > 2
         if (imd(j,1).ge.3) call error (169,rnums(1),imd(j,1),tname)

      end do
c                              look for van laar and/or dqf parameters
c                              or the end of model marker
      call readop (j,j,j,reach,stck,norf,tname)

      do i = 1, nq+nn+ns
         iorig(i) = i 
      end do

      end 

      subroutine nmodel 
c---------------------------------------------------------------------
c nmodel - creates some counters and conversion arrays and stores
c local solution model parameters in global arrays. 
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'
 
      integer itic,i,j,k

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip

      integer ostot
      common/ junk /ostot

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)
c----------------------------------------------------------------------
c                                 make the insp arrays, 

c                                 insp orders endmembers; first independent
c                                 then dependent, last ordered species.
 
c                                 jnsp points only to independent+ordered
c                                 species (used for reciprocal solutions 
c                                 with ordering).

c                                 iy2p points from an endmember in the y
c                                 array (1..istot+norder) to the endmember
c                                 in the p array (1..kstot+norder)

      if (depend) then

         itic = 0 

         do i = 1, ostot
            if (kdsol(i).gt.0) then 
               itic = itic + 1
               insp(itic) = i 
               jnsp(itic) = i 
               iy2p(i) = itic
            end if 
         end do

         do i = 1, mdep 
            insp(itic+i) = jdep(i)
         end do 

         do i = 1, norder 
            insp(ostot+i) = ostot + i
            jnsp(itic+i) = ostot + i 
            iy2p(ostot+i) = itic + i 
         end do 

      else 

         do i = 1, istot + norder
            insp(i) = i 
            jnsp(i) = i 
            iy2p(i) = i 
         end do 
c                                added for reformulated reciprocal
c                                solutions with no dependent endmembers
c                                June 12, 2012, JADC.
c        if (recip) kstot = istot
c                                JADC Nov 22, 2016:
c                                changed to general case, i.e., if there
c                                are no dependent endmembers then
         kstot = istot 

      end if 

      if (depend) then 
c                                 make y2p array       
         do i = 1, kstot + norder 

            do j = 1, mdep

               y2p(i,j) = 0d0

               do k = 1, ndph(j)
                  if (idep(j,k).eq.jnsp(i)) y2p(i,j) = nu(j,k)
               end do 

            end do
         end do

      end if

      end 

      recursive double precision function gproj (id)
c-----------------------------------------------------------------------
c gproj - computes projected free energy of phase id. uproj must be 
c called prior to any call to gproj. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,j

      double precision gphase, gcpd

      external gphase, gcpd 

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn
c--------------------------------------------------------------------- 
      if (id.le.ipoint) then 

         gproj = gcpd (id,.true.)
c                                 if istct > 0 must be some saturated
c                                 components
         if (istct.gt.1) then 
c                                 this is a screw up solution
c                                 necessary cause uf(1) and uf(2)
c                                 are allocated independent of ifct!
            if (ifct.gt.0) then 
               do j = 1, 2
                  if (iff(j).ne.0) gproj = gproj - cp(iff(j),id)*uf(j)
               end do 
            end if 

            do j = 1, isat
               gproj = gproj - cp(icp+j,id) * us(j)
            end do 

         end if 

      else 

         gproj = gphase(id)

      end if 

      end

      double precision function dgdp (id)
c----------------------------------------------------------------------
c subroutine to compute the G derivative of solution id with respect to 
c the concentration of the 1st ordered species. assumes the excess function 
c is not higher than second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,i1,i2,id

      double precision gex,dgex,dsconf,tphi,dtphi

      double precision enth
      common/ cxt35 /enth(j3)

      double precision r,v,tr,pr,ps
      common/ cst5   /v(l2),tr,pr,r,ps
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
c                                 get derivative of excess function
      dgex = 0d0

      if (lexces(id)) then 

         if (.not.llaar(id)) then 
c                                 regular excess function
            do i = 1, jterm(id)
               i1 = jsub(1,i,id)
               i2 = jsub(2,i,id)
               dgex = dgex + w(i)*( pa(i1)*dydy(i2,1,id) 
     *                            + pa(i2)*dydy(i1,1,id))
            end do
         else 
c                                 h&p van laar
            tphi = 0d0
            dtphi = 0d0

            do i = 1, nstot(id)
               tphi = tphi + alpha(i)* y(i) 
               dtphi = dtphi + alpha(i)*dydy(i,1,id)
            end do 

            gex = 0d0 

            do i = 1, jterm(id)
c                                 assume holland powell form, all terms regular
              i1 = jsub(1,i,id)
              i2 = jsub(2,i,id)

              gex = gex + w(i) * y(i1) * y(i2)   
              dgex = dgex + w(i) * (y(i1)*dydy(i2,1,id) 
     *                            + y(i2)*dydy(i1,1,id))

            end do 
c                                note the excess energy is gex/tphi
            dgex = (dgex - dtphi*gex/tphi)/tphi

         end if 

      end if 
c                                 now get dg/dy(jd)
      dgdp = enth(1) + dgex - v(2)*dsconf(id)

      end

      double precision function dsconf (id)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of a dependent species.

c THIS DOES NOT INCLUDE ENDMEMBER CONFIGURATION ENTROPY DERIVATIVES!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision dlnz,dscon,zt,q,dzdy,z,dd
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      dscon = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dd = 0d0 
         q = qmult(i,id)
c                                 get site fractions
         do j = 1, ksp(i,id)

            z = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               z = z + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do 

            zt = zt + z
c                                 sdzdp is (dz(i,j)/dp(k))
            dzdy = sdzdp(1,j,i,id)

            if (dzdy.eq.0d0) cycle

            if (z.gt.0d0) then 
               dlnz = 1d0 + dlog(z)
            else 
c                                 the term should goto -Inf
               dlnz = 1d0
            end if 

            dscon = dscon - q * dzdy * dlnz

            dd = dd + dzdy

         end do 

         z = 1d0 - zt 

         if (z.gt.0d0) then 
            dlnz = 1d0 + dlog(z)
         else 
            dlnz = 1d0
         end if 

         dscon = dscon - q * sdzdp(1,j,i,id) * dlnz

      end do 

      dsconf = dscon

      end

      subroutine gvlaar (jd,id,dg)
c-----------------------------------------------------------------------
c gvlaar evaluates the contribution to the gibbs energy of a pseudocompound
c arising from van laar excess properties. 
c
c now redundant with function gex except gvlaar uses sxs instead of y
c
c input:
c      jd - solution pointer
c      id - compound pointer
c      dg - ideal free energy
c output:
c      dg - ideal free energy + excess energy     
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id,jd,i1,i2

      double precision phi(m4), tphi, dg
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c-----------------------------------------------------------------------
c                                 first compute "volumes"
      tphi = 0d0

      do i = 1, nstot(jd)


c                                 these phi's are the numerator of
c                                 holland & powells "phi's"
         phi(i) =  alpha(i)* xco(jco(id)+i)
c                                 tphi is the denominator of holland & powell's
c                                 phi's
         tphi = tphi + phi(i) 
      end do 
c                                 rearrange tphi
      tphi = 2d0/tphi
c                                 dg is initialized as gph in the calling 
c                                 program
      do i = 1, jterm(jd)
c                                 assume holland powell form, all terms regular
         i1 = jsub(1,i,jd)
         i2 = jsub(2,i,jd)
         dg = dg + w(i)  
     *           * tphi * phi(i1) * phi(i2) / (alpha(i1) + alpha(i2))   

      end do  

      end 

      double precision function hpmelt (im,y)
c----------------------------------------------------------------------
c evaluates the configurational entropy of Holland & Powell's haplogranite 
c melt model, dlnw is S/R.

c modified to use global arrays, 10/26/05.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, im

      double precision dlnw, yfac, yfo, yfa, y(m4)
c                                 global arrays:
      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 special model endmember indexing
      integer jspec
      common/ cxt8 /jspec(h9,m4)
c----------------------------------------------------------------------  

      yfo = 0d0
      yfa = 0d0

      if (jspec(im,1).eq.1.and.y(1).gt.0d0) then
c                                 the quadratic water term
c                                 and temkin renormalization term
         dlnw = -2d0*y(1)*dlog(y(1)) - (1d0 - y(1))*dlog(1d0 - y(1))

      else 
 
         dlnw = 0d0

      end if 
c                                 the molecular entropy
      do i = jspec(im,4), nstot(im)

         if (y(i).le.0d0) cycle
         dlnw = dlnw - y(i) * dlog(y(i))
 
      end do 
c                                 the fe-mg fudge factor

      if (jspec(im,2).ne.0) yfo = y(jspec(im,2))   
      if (jspec(im,3).ne.0) yfa = y(jspec(im,3))

      yfac = yfo + yfa

      if (yfo.gt.0d0.and.yfa.gt.0d0) then 
         dlnw = dlnw - 4d0 * yfo * dlog(yfo/yfac)
         dlnw = dlnw - 4d0 * yfa * dlog(yfa/yfac)
      end if 

      hpmelt = dlnw*r

      end 

      double precision function gmelt (im)
c----------------------------------------------------------------------
c evaluates the configurational entropy of Ghiorso's pMELTS/MELTS model, 
c dlnw is S.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, im
      double precision dlnw, yh2o
c                                 global arrays:
      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 special model endmember indexing
      integer jspec
      common/ cxt8 /jspec(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)
c----------------------------------------------------------------------
      dlnw = 0d0

      if (jspec(im,1).eq.1) then

         yh2o = y(1)
c                                 the quadratic water term, cf hp model. 
         if (yh2o.ne.0d0) dlnw = -(yh2o * dlog(yh2o) 
     *                             + (1d0-yh2o)*dlog(1d0-yh2o))
      
      end if 
c                                 i wouldn't expect to count water here yet again, but that
c                                 is what eqs 1 & 2 of Ghiorso et al. (2002) imply

c                                 in any event equation 1 in Ghiorso et al. 2002
c                                 is almost certainly wrong, so the entropy here
c                                 is computed as in ghiorso & sack 1995 (melts).
c                                 this correction made Nov 4, 2005. This model seems
c                                 to be consistent with Nicholls 1980 CMP 74:211 

c                                 the basic entropy
      do i = 1, mstot(im)

         if (y(i).le.0d0) cycle
         dlnw = dlnw - y(i) * dlog(y(i))
 
      end do 

      gmelt = r*dlnw

      end  

      double precision function slvmlt ()
c----------------------------------------------------------------------
c evaluates the configurational entropy of high T fo-fa-SiO2 melts,
c to use this model, all species (fo, fa, sio2) must be involed in the
c calculation, i.e., y1 is assumed to be fo, etc. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision dlnw, yol, xmg, yq
c                                 global arrays:
      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),x(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)
c----------------------------------------------------------------------
      dlnw = 0d0
c                                 fraction of olivine species
      yol = y(1) + y(2)

      if (yol.ne.0d0) then 
         xmg = y(1)/yol
c                                 entropy within the olivine species
         if (xmg.ne.0d0.and.xmg.ne.1d0) dlnw = -2d0*yol*
     *                  (xmg * dlog(xmg) + (1d0-xmg)*dlog(1d0-xmg)) 
c                                 mixing of the olvine and sio2 species
         yq = 1d0 - yol
         if (yq.gt.1d-15) dlnw = dlnw - 
     *                  (yol * dlog(yol) + (yq)*dlog(yq))              

      end if 

      slvmlt = r*dlnw

      end 

      subroutine ytox (ids)
c----------------------------------------------------------------------
c subroutine to convert endmember fractions (y) to geometric coordinates 
c (x) for solution model ids.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids, i, j
c                                 convert y -> x array
      integer indx

      common/ cxt5i /indx(h9,mst,msp)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c----------------------------------------------------------------------
      do i = 1, istg(ids)
c                                 this conditional is necessary
c                                 because reform can generate a 
c                                 one species site. 
         if (ispg(ids,i).gt.1) then 
c                                 initialize
            do j = 1, ispg(ids,i)
               x(i,j) = 0d0
            end do 
  
            do j = 1, ispg(ids,i)
               x(i,j) = x(i,j) + y(indx(ids,i,j))
            end do 
         else 
c                                 one species site. 
            x(i,1) = 1d0
         end if 
      end do
     
      end 

      subroutine getolx (ids,id)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the static xco array loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, kcoor

      double precision xt 
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 stored x coordinate
      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

c      integer ksmod, ksite, kmsol, knsp
c      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c----------------------------------------------------------------------
      if (lrecip(ids)) then 

         kcoor = ico(id)

         do i = 1, istg(ids)

            xt = 0d0 

            do j = 1, ndim(i,ids)
               kcoor = kcoor + 1
               x(i,j) = xco(kcoor)
               xt = xt + xco(kcoor)
            end do 

            x(i,j) = 1d0 - xt 

         end do 

      else 
c                                 the use if istg(ids) as the site
c                                 index is a hack for reformulated
c                                 reciprocal solutions

c                                 this conditional is not necessary
c                                 for 1-species models cause it will
c                                 never look at the array for composition
c                                 in any case???
c         if (ispg(ids,1).gt.1.or.ksmod(ids).eq.39) then 
         if (ispg(ids,1).gt.1) then 
            i = 1
         else 
            i = 2
         end if 

         xt = 0d0 

         do j = 1, nstot(ids) - 1
            x(i,j) = xco(jco(id)+j)
            xt = xt + x(i,j) 
         end do 

         x(i,j) = 1d0 - xt 

      end if 

      end 

      double precision function gsol1 (id)
c-----------------------------------------------------------------------
c gsol computes the total (excess+ideal) free energy of solution 
c for a solution identified by index ids and composition y(m4) input
c from cxt7, the composition y is the independent endmember fractions
c for all model types except reciprocal solutions, in which case it is 
c the y's for the full reciprocal model.

c gsol assumes the endmember g's have not been calculated by gall and is
c      only called by WERAMI.
c gsol1 is identical to gsol but can only been called after gall and is 
c      only called by VERTEX and MEEMUM. ingsol must be called prior to
c      gsol1 to initialize p-t dependnent model parameters. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k, id

      double precision gg, dg

      double precision omega, hpmelt, slvmlt, gmelt, gfluid, gzero,
     *                 gex, gfesi, gfesic, gfecr1, gerk, ghybrid

      external omega, hpmelt, slvmlt, gmelt, gfluid, gzero, gex, gfesi, 
     *         gfesic, gfecr1, gerk, ghybrid

      integer jend
      common/ cxt23 /jend(h9,m4)

      double precision g
      common/ cst2 /g(k1)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)
c----------------------------------------------------------------------
      gg = 0d0

      if (ksmod(id).eq.2) then 
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
         call gdqf (id,gg,y) 

         gg = gg - t * omega(id,y) + gex(id,y)
c                                 get mechanical mixture contribution
         do k = 1, mstot(id) 
            gg = gg + y(k) * g(jend(id,2+k))
         end do 

      else if (ksmod(id).ge.30.and.ksmod(id).le.31) then 
c                                 -------------------------------------
c                                 Nastia's version of BCC/FCC Fe-Si-C Lacaze and Sundman
c                                 this model has to be called ahead of the standard models
c                                 because it sets lrecip(id) = true.

c                                 initialize p's
         call y2p0 (id)

         gg =  gfesic (y(1),y(3),y(4),
     *                 g(jend(id,3)),g(jend(id,4)),
     *                 g(jend(id,5)),g(jend(id,6)),ksmod(id))

      else if (lrecip(id).and.lorder(id)) then 
c                                 -------------------------------------
c                                 prismatic solution, initialize p's
         call y2p0 (id)
c                                 get the speciation, excess and entropy effects.
         call specis (gg,id)
c                                 decompose the ordered species into 
c                                 the independent disordered species
c                                 i.e., the p0a array becomes the p's if the 
c                                 abundance of the ordered species is 0.
         do k = 1, lstot(id)
c                                 compute mechanical g from these z's, 
c                                 specip adds a correction for the ordered species.
            gg = gg + g(jend(id,2+k)) * p0a(k)
         end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd. gex not neccessary as computed in specip
         call gdqf (id,gg,p0a)

      else if (lorder(id)) then 
c                                 -------------------------------------
c                                 simplicial ordering solutions.
c                                 get mechanical mixture contribution
         do k = 1, lstot(id)  
            pa(k) = y(k)
            p0a(k) = y(k)
            gg = gg + y(k) * g(jend(id,2+k))
         end do 
c                                 get the speciation energy effect
         call specis (dg,id)

         gg = gg + dg 
c                                 get dqf corrections
         call gdqf (id,gg,p0a)

      else if (lrecip(id)) then 
c                                 -------------------------------------
c                                 macroscopic reciprocal solution
c                                 initialize p's
         call y2p0 (id)

         do k = 1, lstot(id)
            gg = gg + g(jend(id,2+k)) * p0a(k) 
         end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd
         call gdqf (id,gg,p0a)

         gg = gg - t * omega(id,p0a) + gex(id,p0a)


      else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos
         gg = gfluid(y(jspec(id,1)))
            
         do k = 1, 2
            gg = gg + gzero(jnd(k))*y(k)
         end do 

      else if (ksmod(id).eq.20) then 
c                                 electrolytic solution, need to check
c                                 that this thing is getting the right 
c                                 partial molar volumes. 
         call slvnt1 (gg)

         call slvnt2 (gg) 

      else if (ksmod(id).eq.25.or.ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp and ghiorso pmelt models 
         if (t.lt.nopt(20)) then 
c                                 t < t_melt, destabilize the melt
            gg = 1d12*p

         else

            call gdqf (id,gg,y) 

            if (ksmod(id).eq.24) then 

               gg = gg - t * hpmelt(id,y) + gex(id,y)

         else 

               gg = gg - t * gmelt(id) + gex(id,y)

            end if 
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               gg = gg + y(k) * g(jend(id,2+k))
            end do 

         end if 

      else if (ksmod(id).eq.26) then 
c                                 ------------------------------------
c                                 andreas salt model
         call hcneos (gg,y(1),y(2),y(3))

         do k = 1, 3
            gg = gg + y(k) * g(jend(id,2+k))
         end do 

      else if (ksmod(id).eq.28) then 
c                                 -------------------------------------
c                                 high T fo-fa-sio2 model  
         call gdqf (id,gg,y) 

         gg = gg - t * slvmlt() + gex(id,y)
c                                 get mechanical mixture contribution
         do k = 1, mstot(id)  
            gg = gg + y(k) * g(jend(id,2+k)) 
          end do 

      else if (ksmod(id).eq.29) then 
c                                 -------------------------------------
c                                 BCC Fe-Si Lacaze and Sundman
         gg =  gfesi(y(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.32) then 
c                                 -------------------------------------
c                                 BCC Fe-Cr Andersson and Sundman
         gg =  gfecr1(y(1),g(jend(id,3)),g(jend(id,4)))

      else if (ksmod(id).eq.39) then
c                                 -------------------------------------
c                                 generic hybrid EoS
c                                 initialize pointer array
         do k = 1, nstot(id)
c                                 sum pure species g's
            gg = gg + g(jnd(k)) * y(k)

         end do
c                                 compute and add in activities
         gg = gg + ghybrid (y)


      else if (ksmod(id).eq.41) then 
c                                 hybrid MRK ternary COH fluid
         call rkcoh6 (y(2),y(1),gg) 

         do k = 1, 3 
            gg = gg + g(jnd(k)) * y(k)
         end do 

      else if (ksmod(id).eq.40) then 
c                                 MRK silicate vapor
         do k = 1, nstot(id) 
            gg = gg + gzero (jnd(k)) * y(k)
         end do 

         gg = gg + gerk(y)

      else 

         write (*,*) 'what the **** am i doing here?'
         call errpau

      end if 

      gsol1 = gg 

      end

      subroutine ingsol (id)
c-----------------------------------------------------------------------
c ingso1 initializes p-t dependent solution model id parameters for gsol1
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c----------------------------------------------------------------------
c                                 evaluate margules coefficients
      call setw (id)
c                                 evaluate dqf coefficients
      call setdqf (id)
c                                 evaluate enthalpies of ordering
      if (lorder(id)) call oenth (id)

      end

      logical function zbad (y,ids)
c----------------------------------------------------------------------
c subroutine to site fractions computed from equations entered by 
c user for configurational entropy (macroscopic form).
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical badz

      external badz

      double precision y(m4),z,zt,n(m10)

      integer i,j,k,ids

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c----------------------------------------------------------------------
c                                 for each site
      do i = 1, msite(ids)

         zt = 0d0
          
         if (qmult(i,ids).ne.0d0) then 
c                                 get site fractions
            do j = 1, ksp(i,ids)

               z = d0(j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)
                  z = z + dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))
               end do 

               if (badz(z)) goto 90

               zt = zt + z
            
            end do
 
            z = 1d0 - zt

            if (badz(z)) goto 90

         else if (ksp(i,ids).gt.1) then 
c                                 variable multiplicity model
            do j = 1, ksp(i,ids)
c                                 molar site population
               n(j) = d0(j,i,ids)
c                                 for each term:
               do k = 1, lterm(j,i,ids)
                  n(j) = n(j) + dcoef(k,j,i,ids) * y(ksub(k,j,i,ids))
               end do

               zt = zt + n(j)
 
            end do 

            if (zt.gt.0d0) then
c                                 if site exists, check fractions
               do j = 1, ksp(i,ids)

                  if (badz(n(j)/zt)) goto 90
            
               end do

            else if (zt.lt.-zero) then 
c                                 negative site?
               goto 90

            end if 

         end if 

      end do

      zbad = .false.

      return

90    zbad = .true.

      end

      double precision function omega (id,y)
c----------------------------------------------------------------------
c subroutine to evaluate the configurational entropy of a solution
c with composition y, including the correction for endmember 
c configurational negentropies. reciprocal end-member composition version.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision z,zt,dlnw,dlnz,y(m4),n(m10)

      integer i,j,k,id
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c----------------------------------------------------------------------
      dlnw = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         dlnz = zt

         if (qmult(i,id).ne.0d0) then 
c                                 standard model with fixed site multiplicity
c                                 get site fractions
            do j = 1, ksp(i,id)
               z = d0(j,i,id) 
c                                 for each term:
               do k = 1, lterm(j,i,id)
                  z = z + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do 

               zt = zt + z
               if (z.gt.0d0) dlnz = dlnz - z * dlog (z)
            
            end do 

            z = 1d0 - zt
            if (z.gt.0d0) dlnz = dlnz - z * dlog(z)
            dlnw = dlnw + qmult(i,id)*dlnz
 
         else if (ksp(i,id).gt.1) then
c                                 variable site multiplicities   
c                                 get site fractions
            do j = 1, ksp(i,id)
               n(j) = d0(j,i,id) 
c                                 for each term:
               do k = 1, lterm(j,i,id)
c                                 n(j) is molar site population
                  n(j) = n(j) + dcoef(k,j,i,id) * y(ksub(k,j,i,id))
               end do 
c                                 zt is the multiplicity
               zt = zt + n(j)
            
            end do   

            if (zt.gt.0d0) then 
c                                 if site has non-zero multiplicity
               do j = 1, ksp(i,id)
c                                 z is site fraction
                  z = n(j)/zt
                  if (z.gt.0d0) dlnz = dlnz - z * dlog(z)
            
               end do       
         
            end if   

            dlnw = dlnw + r*zt*dlnz

         end if 

      end do 
c                                 endmember corrections
      do i = 1, nstot(id)
         dlnw = dlnw - y(i)*scoef(i,id)
      end do

      omega = dlnw 

      end

      subroutine snorm (id,tname)
c------------------------------------------------------------------------
c compute endmember configurational entropies
c site fractions expressed as a function of end-member proportions.
c see deleted snorm0 for bragg-williams site fraction version.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 tname

      logical zbad

      integer h,id,j

      double precision omega

      external omega, zbad
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c----------------------------------------------------------------------
c                                 get normalization constants
c                                 for each endmember
      do h = 1, nstot(id) 
c                                 zero y-array
         do j = 1, nstot(id) 
            y(j) = 0d0 
         end do

         y(h) = 1d0        
c                                 check for valid site fractions
         if (zbad(y,id)) call error (125,z(1),1,tname)
c                                 evaluate S
         scoef(h,id) = omega(id,y)

      end do 

      end 

      subroutine setdqf (id)
c---------------------------------------------------------------------
c setdqf - evaluates dqf coefficients at p and t
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i, id

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 dqf corrections
      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)
c----------------------------------------------------------------------

      do i = 1, jdqf(id)
c                                 index points to the endmember in the full
c                                 model:
         iq(i) = jndq(i,id)

         dq(i) = dqfg(1,i,id) + t*dqfg(2,i,id) + p*dqfg(3,i,id)
         
      end do 

      end 

      subroutine gdqf (id,g,y)
c----------------------------------------------------------------------
c subroutine to evaluate the excess G of solution id with macroscopic
c composition y. assumes a previous call to setdqf. When setdqf is 
c not called use 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id
c                                 working arrays
      double precision g, y(m4)
c                                 dqf corrections
      integer jndq, jdqf, iq

      double precision dqfg, dq 

      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)
c----------------------------------------------------------------------
      do i = 1, jdqf(id)

         g = g + y(iq(i))*dq(i)

      end do  

      end 

      subroutine setw (id)
c---------------------------------------------------------------------
c setw - evaluates margules coeffiecients and, for laar models, size
c parameters for solution model id at p and t
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer i, k, l, i1, i2, id, j

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)
c----------------------------------------------------------------------
      if (extyp(id).eq.1) then 
c                                 redlich kistler is a special case

         do i = 1, jterm(id)
            do j = 1, rko(i,id)

               if (wkl(3,j,i,id).eq.0d0.or.wkl(4,j,i,id).eq.0d0.or.
     *             wkl(5,j,i,id).eq.0d0) then

                   wl(j,i) = wkl(1,j,i,id) + t*wkl(2,j,i,id)
               else
                   wl(j,i) = wkl(1,j,i,id) + t*wkl(2,j,i,id) + 
     *                   4d0*((-1d0*wkl(5,j,i,id) - 
     *                   dsqrt((2d0*wkl(5,j,i,id)*p + 
     *                   wkl(4,j,i,id))/wkl(4,j,i,id)))*
     *                   wkl(3,j,i,id)*wkl(4,j,i,id)*dexp(-(-1d0 + 
     *                   dsqrt((2d0*wkl(5,j,i,id)*p + 
     *                   wkl(4,j,i,id))/wkl(4,j,i,id)))/
     *                   wkl(5,j,i,id)) + 
     *                   wkl(3,j,i,id)*wkl(4,j,i,id)*
     *                   (wkl(5,j,i,id)+1d0))
               end if

            end do
         end do

         return

      end if 

      do i = 1, jterm(id)
         w(i) = wgl(1,i,id) + t*wgl(2,i,id) + p*wgl(3,i,id)
      end do

      if (llaar(id)) then 

         do i = 1, nstot(id) 

            alpha(i) = vlar(1,i,id) 
     *               + t * vlar(2,i,id) + p * vlar(3,i,id)

         end do

         do i = 1, jterm(id)
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)
            w(i) = 2d0 * w(i)
     *                 * alpha(i1)*alpha(i2) / (alpha(i1) + alpha(i2))
         end do

      end if 

      if (lorder(id)) then 
c                                 set higher order derivatives for 
c                                 speciation
         do k = 1, nord(id) 
            dt(k) = 0d0
            do l = 1, nord(id)
               d2gx(l,k) = 0d0
            end do  
         end do
c                                 for both laar and regular need
c                                 the d(gex)/dy(k)/dy(l)
         do i = 1, jterm(id)
            do k = 1, nord(id) 
               do l = 1, nord(id)
                  d2gx(l,k) = d2gx(l,k) + w(i) * dppp(l,k,i,id)
               end do  
            end do 
         end do


         if (llaar(id)) then 
c                                 for laar also need:
            do i = 1, nstot(id)
               do k = 1, nord(id)  
c                                 dt, derivative of sum(phi)
                  dt(k) = dt(k) + alpha(i)*dydy(i,k,id)
               end do 
            end do 

         end if 

      end if 

      end 

      double precision function gex (ids,y)
c-----------------------------------------------------------------------
c evaluate the excess function for solution model ids. assuming no prior
c call to set coefficients (as in function gexces). 
c input:
c      ids - solution pointer
c      y - composition array

c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,ids

      double precision y(m4), tphi, xpr, lex(m17,m18)

      double precision z, pa, p0a, x, w, yy, wl
      common/ cxt7 /yy(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 local alpha
      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c----------------------------------------------------------------------
      gex = 0d0 

      if (extyp(ids).eq.1) then 
c                                 redlich kistler; expand polynomial
         do i = 1, jterm(ids)
            do j = 1, rko(i,ids)
               lex(j,i) = 0d0
            end do
         end do

         do i = 1, jterm(ids)
            do j = 1, rko(i,ids)
c                                 interchanged subscripts, G Helffrich, 4/8/16.
               lex(j,i) = lex(j,i) + wl(j,i)*
     *                        (y(jsub(1,i,ids))-y(jsub(2,i,ids)))**(j-1)
            end do
         end do

         do i = 1, jterm(ids)
            do j = 1, rko(i,ids)
               gex = gex + lex(j,i)*y(jsub(1,i,ids))*y(jsub(2,i,ids))
            end do
         end do

      else if (lexces(ids)) then 

         if (llaar(ids)) then 
c                                 holland & powell's version of the van laar
c                                 first compute "volumes"
            tphi = 0d0

            do i = 1, nstot(ids)
c                                 tphi is the sum of holland & powell's
c                                 phis
               tphi = tphi + alpha(i) * y(i)

             end do 
c                                 dg is initialized as gph in the calling 
c                                 program
            do i = 1, jterm(ids)
c                                 holland powell form, all terms regular
               gex = gex + w(i) * y(jsub(1,i,ids)) * y(jsub(2,i,ids))  

            end do 

            gex = gex/tphi 

         else 
c                                 macroscopic margules formulation by default
            do i = 1, jterm(ids)

               xpr = 1d0

               do j = 1, jord(ids)
                  if (jsub(j,i,ids).eq.0d0) exit
                  xpr = xpr * y(jsub(j,i,ids))
               end do 

               gex = gex + xpr * w(i) 

            end do  

         end if

      end if 

      end 

      subroutine endcp (jd,id,ids)
c------------------------------------------------------------------------
c compute the composition of endmember id, for solution ids and load it
c into the jd'th position of the x3 array.
c------------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jd, h, i, j, id, ids
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
c                                 figure out which endmember we
c                                 are looking at:
      do h = 1, mstot(ids)
         if (id.eq.jend(ids,2+h)) exit
      end do
c                                 zero x-array
      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            x3(jd,i,j) = 0d0
         end do
c                                 now assign endmember fractions
         x3(jd,i,kmsol(ids,knsp(h,ids),i)) = 1d0

      end do

      end 

      subroutine gmodel (im,tname,wham)
c---------------------------------------------------------------------
c qmodel - stores ALL solution model parameters in global arrays
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      character tname*10, sname*10

      logical add, bad, wham, zbad

      integer im, nloc, i, j, ind, id, jd, k, l,itic,ii,imatch, killct,
     *        killid(20), inc

      double precision dinc,dzt,dx,gcpd

      external gcpd, zbad

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      logical badend, sck, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),badend(m4,h9),sck(h9),nrf(h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
                             
      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer eos
      common/ cst303 /eos(k10)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer nsub,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub(m10,m11,m0,m12),
     *      nttyp(m10,m11,m0)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision t, p, xco2, u1, u2, tr, pr, r, ps
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip
c                                 GLOBAL SOLUTION PARAMETERS:
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
      
      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)
c                                 convert y -> x array
      integer indx
      common/ cxt5i /indx(h9,mst,msp)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m4)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *           xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)
c                                 special model endmember indexing
      integer jspec
      common/ cxt8 /jspec(h9,m4)

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 dqf parameters
      integer idqf,indq
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf

      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c                                 temporary stretching coordinate
c                                 parameters
      double precision yin
      common/ cst50 /yin(ms1,mst)
c                                 parameters for autorefine
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)

      logical refine
      common/ cxt26 /refine
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)

      integer ostot
      common/ junk /ostot

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(5,2)

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      character specie*4
      integer jsp, ins
      common/ cxt33 /jsp,ins(nsp),specie(nsp)

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      character mname*8
      common/ cst18a /mname(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer iam
      common/ cst4 /iam

      integer tnq,tnn,tns
      common/ cxt337 /tnq,tnn,tns

      double precision stch
      common/ cst47 /stch(h9,mst,msp,4)

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c----------------------------------------------------------------------
c                                 auto_refine changes
      if (refine) then
c                                 check for consistent auto-refine data
         read (n10,'(a)') sname
         if (tname.ne.sname) call error (63,r,i,'GMODEL')

      end if 
c                                 initialize autorefine arrays
      stable(im) = .false.
      limit(im) = .false.
      relax(im) = .true.
c                                 initialize compositional distances
      do i = 1, icp
         dcp(i,im) = 0d0
      end do 
c                                 check endmember counters:
      if (im.gt.h9) call error (52,dq(1),idqf,'GMODEL')
c                                 check for inconsistent model reformation
      if (kstot+mdep.gt.jstot) call error (76,dq(1),idqf,tname)
c                                 set up simple counters for
c                                 charge balance models
      if (jsmod.eq.20) then

         isite = 1

         nq = tnq
         ns = tns 
         nn = tnn 

         nqs = nn + nq + ns
         nqs1 = nqs - 1
         nq1 = nq - 1
         sn = nn + ns
         ns1 = sn - 1
         sn1 = ns + 1
         qn = nn + nq

         isp(1) = nqs

         j = 0 

         do i = 1, nqs
            jmsol(i,1) = i
         end do

      end if 
c                                 number of dependent + independent - ordered endmembers
c                                 prismatic space
      mstot(im) = istot
c                                 number of dependent + independent - ordered endmembers
c                                 prismatic space + orphans
      pstot(im) = ostot
c                                 number of dependent + independent - ordered endmembers
c                                 prismatic space - orphans
      qstot(im) = ostot
c                                 number of independent + ordered endmebers
      nstot(im) = kstot + norder
c                                 number of independent disordered endmembers
      lstot(im) = kstot
c                                 chemical mixing sites
      istg(im) = isite
c                                 site check override
      sck(im) = stck
c                                 don't refine endmembers
      nrf(im) = norf
c                                 number of ordered species
      nord(im) = norder 
c                                 number of species and multiplicity and
c                                 site ranges
      inc = 0
      if (jsmod.eq.9.or.jsmod.eq.10) inc = 1

      ostg(im) = isite + inc
      mcoor(im) = 0
      ncoor(im) = 0
      nsum(im) = 0

      do i = 1, ostg(im)

         imlt(im,i) = ist(i)

         if (i.le.isite) then
            ndim(i,im) = isp(i) - 1
            ispg(im,i) = isp(i)
            nsum(im) = nsum(im) + ndim(i,im)
         else 
            ndim(i,im) = isp(i)
            ispg(im,i) = isp(i) + 1 
         end if

         ncoor(im) = ncoor(im) + ispg(im,i)
         mcoor(im) = mcoor(im) + ndim(i,im)

         do j = 1, ndim(i,im)
c                                 subdivision override (iopt(13))
            if (iopt(13).eq.1) then
c                                 make linear
               imd(j,i) = 0

            else if (iopt(13).eq.2) then 
c                                 make all stretch (if not already)
               if (imd(j,i).eq.0) imd(j,i) = 1

            end if 
c                                 for aqueous solutions always use
c                                 model specified increments
            if (jsmod.ne.20.and.nopt(13).gt.0d0) then 
c                                 otherwise:
               if (imd(j,i).eq.0) then 
c                                 cartesian
                  xnc(i,j) = nopt(13)

               else 
c                                 conformal
                  xnc(i,j) = 1d0/nopt(13)

               end if

            else if (jsmod.eq.20) then 

               if (imd(j,i).eq.1) xnc(i,j) = 1d0/xnc(i,j)

            end if
c                                 perturb xmn by a per-mil scale increment to 
c                                 reduce compositional degeneracies. 
            if (icopt.le.3.and.nopt(13).gt.0d0) 
     *         xmn(i,j) = (1d0 + nopt(15)*float(im-5)) * xmn(i,j)
c                                 save solution model values as hard limits for 
            xmno(im,i,j) = xmn(i,j)
            xmxo(im,i,j) = xmx(i,j)
c                                 set stretch parameters according to xmn specified 
c                                 in the solution model:
            if (imd(j,i).ne.0) then

               if (xmn(i,j).eq.0) then

                  stch(im,i,j,1) = bm1
                  stch(im,i,j,2) = bp1
                  stch(im,i,j,3) = bpm
                  stch(im,i,j,4) = lbpm

               else 
c                                 this makes smallest increment comparable, but less 
c                                 than, xmn. they become equal in the limit xmn->0.
                  stch(im,i,j,1) = xmn(i,j)
                  stch(im,i,j,2) = xmn(i,j) + 2
                  stch(im,i,j,3) = stch(im,i,j,2)/stch(im,i,j,1)
                  stch(im,i,j,4) = dlog(stch(im,i,j,3))
                  xmn(i,j) = 0d0

               end if

            end if 
c                                 set initial resolution
c                                 ------------------------------------
c                                 auto_refine segment
            if (refine) then 
c                                 new values from autorefine file
               read (n10,*) xmn(i,j),xmx(i,j)

               if (icopt.lt.4) then 
c                                 set slop to the initial spacing
                  dinc = xnc(i,j)
                
                  xnc(i,j) = xnc(i,j)/nopt(17)

               else 
c                                 adaptive minimization:
c                                 set slop to the final compositional
c                                 resolution of the exploratory stage
                  dinc = rid(4,1)
c                                 adaptive use refine factor I
                  xnc(i,j) = xnc(i,j)/nopt(17) 

               end if 
c                                 widen the range by the exploratory resolution
               xmx(i,j) = xmx(i,j) + dinc
               xmn(i,j) = xmn(i,j) - dinc

               if (xmx(i,j).gt.1d0) xmx(i,j) = 1d0
               if (xmn(i,j).lt.0d0) xmn(i,j) = 0d0
 
            end if
c                                 -------------------------------------
c                                 stretching stuff
            if (imd(j,i).gt.0) then
c                                 check for old subdivision schemes
               if (imd(j,i).gt.1.or.xnc(i,j).lt.1d0) 
     *            call error (62,nopt(13),imd(j,i),tname)

            end if

            imdg(j,i,im) = imd(j,i)
            xmng(im,i,j) = xmn(i,j)
            xmxg(im,i,j) = xmx(i,j)
            xncg(im,i,j) = xnc(i,j)

         end do

      end do
c                                 initialize high/low ranges
      do i = 1, isite
         do j = 1, ndim(i,im)

            xlo(j,i,im) = 1d0
            xhi(j,i,im) = 0d0

         end do
      end do
c                                 set reach factors
      if (.not.refine.and.(iam.eq.1.or.iam.eq.15).and.iopt(20).ne.2.or.
     *                                 iopt(20).eq.0) then
c                                 if vertex and not in the refine stage
c                                 shut off reach increments
         reachg(im) = nopt(21)/2d0

      else if (reach.le.nopt(23)) then

         reachg(im) = nopt(21)*(1d0 + nopt(23))/2d0

      else 

         reachg(im) = nopt(21)*(1d0 + reach)/2d0

      end if 
c                                 -------------------------------------
c                                 classify the model
      ksmod(im) = jsmod
c                                 -------------------------------------
c                                 save the excess terms.                              
      jterm(im) = iterm
      jord(im) = iord
      extyp(im) = xtyp

      do i = 1, iterm 


         if (xtyp.eq.0) then 
c                                 arbitrary expansion
            do j = 1, m3
               wgl(j,i,im) = wg(i,j)
            end do

         else 
c                                 redlich-kistler
            rko(i,im) = rkord(i)

            do k = 1, rkord(i)
               do j = 1, m16
                  wkl(j,k,i,im) = wk(j,k,i)
               end do  
            end do

         end if  

         do j = 1, iord
c                                 isub points to the position in the list
c                                 of endmembers potentially including dependent 
c                                 species. use iy2p to convert to independent
c                                 endmember pointers.
            if (isub(i,j,1).eq.0) then
c                                 term may be of order < iord
               jsub(j,i,im) = 0

            else

               jsub(j,i,im) = iy2p(isub(i,j,1))

               if (kdsol(isub(i,j,1)).eq.-2) then 

                  call error (77,r,i,'dependent endmember '
     *                 //mname(iorig(isub(i,j,1)))//' in solution '
     *                 //tname//'appears in an excess term.')

               end if 

            end if 

         end do

      end do 


      do i = 1, pstot(im)
c                                 initialize invalid speciation flag
         badend(i,im) = .false.
c                                 save global copy of kdsol
         ldsol(i,im) = kdsol(i) 
c                                 insp points to the original position 
c                                 of endmember i in the solution model input:
         knsp(i,im) = insp(i)
c                                 kmsol points to the species on the j'th site
c                                 of the i'th endmember, used for the xtoy
c                                 conversion
         do j = 1, isite
            kmsol(im,i,j) = jmsol(i,j)
         end do 

      end do 
c                                 ----------------------------------------------
c                                 configurational entropy models

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      nloc = 0

      do i = 1, nsite 
c                                 eliminate sites with 1 species
         if (nspm1(i).eq.0) cycle 

         nloc = nloc + 1
c                                 # of species, and site r*multiplicty.
         qmult(nloc,im) = r*smult(i)
         ksp(nloc,im) = nspm1(i) 
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
         do j = 1, nspm1(i)
c                                 # of terms in the 
c                                 site fraction function and a0.
            lterm(j,nloc,im) = nterm(i,j)
            d0(j,nloc,im) = a0(i,j)
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 term coefficient amd species index:
               dcoef(k,j,nloc,im) = acoef(i,j,k)
               ksub(k,j,nloc,im) = iy2p(nsub(i,j,k,1))

               if (kdsol(nsub(i,j,k,1)).eq.-2) then 

                  call error (77,r,k,'dependent endmember '
     *              //mname(iorig(nsub(i,j,k,1)))//' in solution '
     *              //tname//'appears in a site fraction expression.')

               end if 

            end do 
         end do
      end do  
c                                 number of distinct identisites for entropy
      msite(im) = nloc
c                                 --------------------------------------
c                                 van laar volumes, and pointers for "in" endmembers
      do i = 1, nstot(im)
c                                 if the solution is van laar save
c                                 the "volume" function.
         if (laar) then  
            do l = 1, m3
               vlar(l,i,im) = vlaar(l,jnsp(i))
            end do 
         end if 
c                                 initialize scoef's to zero for config 
c                                 entropy calculation (done by snorm).
         scoef(i,im) = 0d0

      end do   
c                                 -------------------------------------
      if (depend) then 
c                                 march, 2017: deleted y2p4z routine that converted
c                                 z(y) expressions to z(p), i.e., by simply
c                                 eliminating dependent endmembers. 

c                                 save y -> p array
         ndep(im) = mdep

         do i = 1, nstot(im)
            do j = 1, mdep
               y2pg(j,i,im) = y2p(i,j)
            end do
         end do
c                                 check for invalid dependent endmembers, these
c                                 are occasionally used as place holders:
         bad = .false.

         if (ksmod(im).eq.9.or.ksmod(im).eq.10) then 

            do j = 1, ndim(3,im)
               x(3,j) = 0d0
            end do 

            x(3,j) = 1d0

         end if 

         do j = 1, mdep

            do i = 1, mstot(im)
               y(i) = 0d0
            end do 

            y(knsp(lstot(im)+j,im)) = 1d0 

            call y2p0 (im)
c                                 check for invalid site fractions, this is only necessary
c                                 for H&P models that assume equipartition (which is not 
c                                 implemented). 
            if (zbad(pa,im)) then

               bad = .true.

               if (iam.lt.3.or.iam.eq.4.or.iam.eq.15) 
     *            call warn (59,y(1),i,
     *            mname(iorig(knsp(lstot(im)+j,im)))
     *            //' in solution model '//tname)

               badend(knsp(lstot(im)+j,im),im) = .true.

            end if

         end do

         if (bad.and.stck) call error (78,y(1),i,tname)

      end if 
c                                 -------------------------------------
c                                 dqf parameters
      jdqf(im) = idqf

      do i = 1, idqf 
c                                 shift pointer from y array to p array
         jndq(i,im) = iy2p(indq(i))
         do j = 1, m3
            dqfg(j,i,im) = dqf(j,i)
         end do 
      end do 
c                                 -------------------------------------
c                                 if nsite ne 0 get "normalization" constants (endmember
c                                 configurational entropies) for entropy model:
      if (nsite.ne.0) call snorm (im,tname)
c                                 -------------------------------------
      if (order) then 
c                                 models with speciation: 
         do j = 1, norder 

            do i = 1, 3 
               deph(i,j,im) = denth(j,i) 
            end do 

            nrct(j,im) = nr(j)

            do i = 1, nr(j)
               ideps(i,j,im) = iy2p(iddeps(i,j))
            end do 
c                                 stoichiometric limits on ordered species
            if (limn(j).gt.0) then

               ln(j,im) = limn(j)
 
               do k = 1, limn(j)
c                                 for each limit, number of p0 terms
                  lt(k,j,im) = limt(k,j)

                  do i = 1, limt(k,j)
c                                 for each p0 term
                     lc(i,k,j,im) = limc(i,k,j)
c                                 convert index from y to p pointer and save
                     lid(i,k,j,im) = iy2p(limid(i,k,j))

                  end do 
c                                 the constant and delta
                  l0c(1,k,j,im) = limc(i,k,j)
                  l0c(2,k,j,im) = limc(i+1,k,j)
c                                 and the number of p terms...
                  jt(k,j,im) = jimt(k,j)

                  do i = 1, jimt(k,j)
c                                 for each p0 term
                     jc(i,k,j,im) = jimc(i,k,j)
c                                 convert index from y to p pointer and save
                     jid(i,k,j,im) = iy2p(jimid(i,k,j))

                  end do 

               end do 

            end if 

         end do 
c                                 classify multiple species models according
c                                 to whether the disordered reactants are 
c                                 partially or completely correlated, assume
c                                 anti-correlation is not a possible case.
         icase(im) = 0

         if (norder.gt.1) then 

            imatch = 0

            do j = 1, nr(1) 

               id = ideps(j,1,im)
     
               do i = 1, nr(2)          
                  if (id.eq.ideps(i,2,im)) then
                     imatch = imatch + 1
                     exit
                  end if 
               end do 

            end do 
 
            if (imatch.eq.1) then 
c                                 if match = 1 one species didn't match
c                                 assume partial correlation
               icase(im) = 2
            else if (imatch.ge.2) then 
               icase(im) = 1
            end if 

         end if  
c                                first create derivatives of endmember
c                                fractions with respect to the ordered
c                                species:
         do j = 1, norder 

            do i = 1, nstot(im)
               dydy(i,j,im) = 0d0
            end do
c                                derivative of the ordered species with
c                                respect to itself:
            dydy(kstot+j,j,im) = 1d0
c                                each ordered species decomposes to
c                                two disordered species iddeps(1-2,j)
c                                depnu is the stoichiometric coefficient
c                                of the disordered species in the ordered
c                                species.

c                                derivatives of the consituent species 
c                                with respect to the ordered species
            do i = 1, nr(j)
               dydy(ideps(i,j,im),j,im) = dydy(ideps(i,j,im),j,im) 
     *                                  - depnu(i,j)
            end do

         end do 
c                                evaluate the second derivative of each
c                                pi*pj term in the excess function with
c                                respect to kth species
         do i = 1, iterm 
            do j = 1, norder
               do k = 1, norder

                  dppp(k,j,i,im) =  dydy(jsub(1,i,im),k,im)
     *                             *dydy(jsub(2,i,im),j,im)
     *                           + 
     *                              dydy(jsub(2,i,im),k,im)
     *                             *dydy(jsub(1,i,im),j,im)  
               end do                     
            end do 
         end do 

         if (depend) then 
c                                make an array to get the the fractions
c                                of the disordered species in the ordered
c                                species. this is essentially identical to 
c                                the dydy array formed above.                                 
            do i = 1, kstot

               do k = 1, norder

                  dvnu(i,k,im) = 0d0

                  do j = 1, nr(k)
                     if (i.ne.ideps(j,k,im)) cycle
                     dvnu(i,k,im) = dvnu(i,k,im) + depnu(j,k)
                  end do

               end do  
            end do  

         end if
c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
         do i = 1, msite(im)
c                                 for each species, read
c                                 function to define the
c                                 site fraction of the species:
            do k = 1, norder
               sdzdp(k,ksp(i,im)+1,i,im) = 0d0 
            end do 

            do j = 1, ksp(i,im)
c                                 # of terms in the 
c                                 site fraction function and a0.
               do l = 1, norder
                  sdzdp(l,j,i,im) = 0d0
               end do 
c                                 for each term:
               do k = 1, lterm(j,i,im)
c                                 endmember indexes  
                  ind = ksub(k,j,i,im) 
c                                 get derivatives, of species fractions
c                                 with respect to ordered species
                  do l = 1, norder 
                     itic = 0 
                     do ii = 1, nr(l)
                        if (ind.eq.ideps(ii,l,im)) then
                           sdzdp(l,j,i,im) = sdzdp(l,j,i,im) 
     *                     + dydy(ideps(ii,l,im),l,im)*dcoef(k,j,i,im)  
                           itic = itic + 1
c                                 high order terms not allowed
                           if (itic.gt.1) 
     *                        call error (999,r,801,'GMODEL')
                        end if                       
                     end do
c                                 the derivative of a term with the 
c                                 ordered species.
                     if (ind.eq.kstot+l) 
     *                  sdzdp(l,j,i,im) = sdzdp(l,j,i,im) 
     *                                  + dcoef(k,j,i,im) 

                  end do   
               end do 
            end do 
         end do
c                                 multiply each dzdp by qmult (R*site
c                                 multiplicity) to reduce operation 
c                                 count in evaluating derivatives. 
         do k = 1, norder 
            do i = 1, nsite 

               dzt = 0d0 

               do j = 1, ksp(i,im) 
                  if (dabs(sdzdp(k,j,i,im)).lt.zero) 
     *                     sdzdp(k,j,i,im) = 0d0
                  dzt = dzt + sdzdp(k,j,i,im)
               end do 

               if (dabs(dzt).lt.zero) dzt = 0d0
               sdzdp(k,j,i,im) = -dzt

            end do 
         end do 

      end if 
c                                 ----------------------------------------------
c                                 models with special endmember indexing:  
      if (jsmod.eq.24.or.jsmod.eq.25.or.jsmod.eq.27) then 
c                                 hp & ghiroso models:
         do i = 1, 3
            jspec(im,i) = 0 
         end do 
c                                 set start index assuming no water:
         jspec(im,4) = 1  

         if (iorig(1).eq.1) then 
c                                 h2o is present:
            jspec(im,1) = 1
c                                 set start index to avoid h2o:
            jspec(im,4) = 2
            if (iorig(2).eq.2) then
               jspec(im,2) = 2
               if (iorig(3).eq.3) jspec(im,3) = 3
            else if (iorig(2).eq.3) then 
               jspec(im,3) = 2
            end if 
         else if (iorig(1).eq.2) then
c                                 h2o absent, fo (in hp) is first endmember
            jspec(im,2) = 1
            if (iorig(2).eq.3) jspec(im,3) = 2
         else if (iorig(1).eq.3) then 
c                                 h2o and fo absent, fa (in hp) is first endmember
            jspec(im,3) = 1
         end if
c                                 set the ifp flag for t_melt option
         do i = 1, mstot(im)
c                                 insp points to the original position 
c                                 of endmember i in the solution model input:
            ifp(kdsol(knsp(i,im))) = -1

         end do 

      else if (jsmod.eq.0) then
c                                 fluid eos, make pointer to co2
         do i = 1, 2
            id = kdsol(insp(i))
            if (cp(2,id).ne.0d0) then
               jspec(im,1) = i
               exit
            end if 
         end do 

      else if (jsmod.eq.39.or.jsmod.eq.20) then
c                                 generic hybrid fluid EoS, convert endmember EoS
c                                 flag to fluid species indices
         if (jsmod.eq.20) then 
            j = ns
         else 
            j = istot
         end if 

         do i = 1, j
            k = eos(kdsol(insp(i)))
            if (k.gt.200) then
               k = k - 200
               call warn (99,0d0,0,'using special components '//
     *               'with solution model '//tname//' may lead '//
     *               'to incorrect results')
            else if (k.gt.100) then
               k = k - 100
            else
               call error (72,0d0,0,'invalid endmember EoS '//
     *               'specification in solution model '//tname)
            end if 

            jspec(im,i) = k

         end do

      else 
c                                 save original endmember indexes for hard-wired 
c                                 solution models
         do i = 1, ostot
            jspec(im,i) = iorig(i)
         end do 

      end if
c                                  set fluid flags for non-special case melts
      if (lname(im).eq.'liquid'.or.lname(im).eq.'fluid') then 
c                                  set ifp for t_melt and melt_is_fluid
         do i = 1, lstot(im)
c                                 of endmember i in the solution model input:
            if (lname(im).eq.'liquid') then 
               ifp(kdsol(knsp(i,im))) = -1
            else 
               ifp(kdsol(knsp(i,im))) = 1
            end if 

         end do

      end if

      if (jsmod.ne.20) then 
c                                 -------------------------------------
c                                 create a y -> x array, this array is 
c                                 to be used to convert endmember fractions (y's)
c                                 back to geometric coordinates (x's).
         do i = 1, isite
            do j = 1, isp(i)
c                                 now find all endmembers with
c                                 species j on site i, this method
c                                 is inefficient but idependent of
c                                 endmember order.   

c DEBUG DEBUG just patched in ostot.
               do k = 1, ostot
                  if (jmsol(k,i).eq.j) indx(im,i,j) = k
               end do 
            end do 
         end do 

      end if 
      
      if (ostot+norder.gt.m4) call error (39,0d0,m4,'INPUT9')    

      smod(im) = .true.
      pmod(im) = .true.
      killct = 0
c                                 set fluid flag from the full name of the 
c                                 phase:
      if (lname(im).eq.'fluid') then 
         fp(im) = .true.
      else 
         fp(im) = .false.
      end if 
c                                 classify liquid model as fluid/not fluid
c                                 according to the melt_is_fluid option, this
c                                 is only relevant for WERAMI
      if (lname(im).eq.'liquid'.and.iam.eq.3.and.lopt(6))
     *                                                 fp(im) = .true.

      do i = 1, kstot
c                                 pointer to endmember
         id = kdsol(insp(i))
c                                 figure out the compositional distance between
c                                 the endmembers, this is used to scale the solvus
c                                 tolerance
         do j = i+1, kstot

            jd = kdsol(insp(j))
c                                 

            do k = 1, icp
                
               dx = dabs(cp(k,id) - cp(k,jd))
               
               if (dx.lt.nopt(5)) then 
                  cycle
               else if (dcp(k,im).lt.dx) then 
                  dcp(k,im) = dx
               end if 
               
            end do 

         end do 
c                                 set ifp for models with an 
c                                 endmember identified as "fluid"
         if (ifp(id).gt.0) fp(im) = .true.  
c
         jend(im,2+i) = id
c                                 set shear/bulk moduli flags
         if (iemod(id).eq.0) smod(im) = .false.
         if (iemod(id).lt.2) pmod(im) = .false.
c                                 look for endmembers to be killed
         if (iend(insp(i)).ne.2) cycle

         add = .true. 

         do j = 1, killct

            if (killid(j).eq.id) then
               add = .false.
               exit
            end if 
         end do 

         if (add) then 

            killct = killct + 1
            if (killct.gt.10) call error (999,wg(1,1),killct,tname)
            killid(killct) = id

         end if 

      end do 
c                                 this looks like bad news, for laar/recip
c                                 or laar/order, but appears to be overridden
c                                 by use of logical classification variables,
c                                 in which case, why is it here????
      if (laar) then 

         if (iterm.eq.0) laar = .false.

      end if 
c                                 set type flags, presently no provision for 
c                                 bw summation
      llaar(im) = .false.
      lexces(im) = .false.
      lorder(im) = .false.
      lrecip(im) = .false.
      
      if (iterm.gt.0) then 
         lexces(im) = .true.
         if (laar) then 
            llaar(im) = .true.
            extyp(im) = 2
         end if    
      end if 

      if (order) lorder(im) = .true.
c                                 the ksmod(im) test is made because
c                                 reform may dump the dependent endmembers
c                                 setting depend = .false., while retaining
c                                 a dummy site with no mixing. reform should
c                                 be redone to truly reformulate multiple
c                                 models to single site models. 
c                                 a non-reciprocal model (ksmod=5) with 
c                                 dependent endmembers is also classified
c                                 as lrecip.
      if (recip.or.depend) lrecip(im) = .true. 

      if (.not.lopt(3)) then 
c                                 hard limits are off, set limits to 0/1
         do i = 1, isite 
            do j = 1, isp(i) - 1

               xmxo(im,i,j) = 1d0
               xmno(im,i,j) = 0d0

            end do 
         end do
      end if   

      if (ksmod(im).eq.0. or.ksmod(im).eq.20.or.ksmod(im).eq.39.or.
     *    ksmod(im).eq.40.or.ksmod(im).eq.41) then
c                                 routines that invoke fluid EoS, as currently 
c                                 configured these will set ins/isp arrays only 
c                                 once. therefore some parameters and indices 
c                                 can be saved in simple arrays
         call setsol (im,wham) 

      end if 
c                                 set independent species names and counters for output
c                                 special cases first:
      if (ksmod(im).eq.0.or.ksmod(im).eq.40.or.ksmod(im).eq.41) then

         spct(im) = jsp

         do i = 1, jsp
            spnams(i,im) = specie(ins(i))
         end do 

      else 

         spct(im) = nstot(im) 
c                                 independent disordered species
         do i = 1, lstot(im)
            spnams(i,im) = mname(iorig(knsp(i,im)))
         end do 
c                                 ordered species
         do i = lstot(im)+1, nstot(im)
            spnams(i,im) = mname(iorig(ndep(im)+i))
         end do 

         if (ksmod(im).eq.29.or.ksmod(im).eq.32) then 
c                                 BCC Fe-Si Lacaze and Sundman (29) 
c                                 BCC Fe-Cr Andersson and Sundman (32)
            spct(im) = 4
            spnams(3,im) = 'osp'
            spnams(4,im) = 'asp'

         end if 

      end if 

      end 

      double precision function strtch (y)
c----------------------------------------------------------------------
c get the x-value from the unstretched coordinate y
c----------------------------------------------------------------------
      implicit none

      double precision y, t

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      t = bpm**(1d0-y)

      strtch = (bp1 - bm1*t)/(1d0 + t)

      end 

      double precision function unstch (x)
c----------------------------------------------------------------------
c get the y-value from the stretched coordinate x
c----------------------------------------------------------------------
      implicit none

      double precision x

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm

      unstch = 1d0 - dlog((bp1-x)/(bm1+x))/lbpm

      end

      subroutine factor (a,n,ipvt,ier)
c-----------------------------------------------------------------------
c factor is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision a(k8,k8),d(k8),rmax,tmax,temp,ratio

      integer ipvt(k8),i,j,k,ier,ip1,n,istr,nm1

      double precision wmach(9)
      common /ax02za/wmach
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n

         ipvt(i) = i
         rmax = 0d0

         do j = 1, n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.wmach(3)) goto 9000

         d(i) = rmax

      end do 
c                            begin decomposition: 
      nm1 = n-1
c
      do i = 1, nm1
         ip1 = i+1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/d(i)
         istr = i
         do j = ip1, n
            tmax = dabs(a(j,i))/d(j)
            if (tmax.gt.rmax) then
               rmax = tmax
               istr = j
            end if 
         end do

         if (dabs(rmax).lt.wmach(3)) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = d(istr)
            d(istr) = d(i)
            d(i) = temp

            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 
         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
         a(j,i) = a(j,i)/a(i,i)
         ratio = a(j,i)
            do k = ip1,n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do
      end do 

      if (dabs(a(n,n)).lt.wmach(3)) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine y2p0 (id)
c-----------------------------------------------------------------------
c y2p0 converts the y array of disordered dependent and independent 
c species abundance to the p0 array of the independent (ordered and 
c disordered) species. the p0 array gives the minimum possible 
c concentrations of the ordered species (the stable abundances are 
c determined by solving the speciation problem). 

c for non-reciprocal solutions the y and p0 arrays are identical.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,k,l
c DEBUG DEBUG
      logical zbad
      external zbad

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)
c-----------------------------------------------------------------------
c                                 for orphan vertex models could check
c                                 that the prismatic fraction is > 0 
c                                 first, this would save setting/resetting
c                                 p0a's as currently done here/
c                                 --------------------------------------
c                                 convert y's to p's
c                                 initialize ordered species
      do k = 1, nord(id)
         p0a(lstot(id)+k) = 0d0
      end do

      do k = 1, nstot(id)
c                                 initialize the independent species
c                                 other then the ordered species
         if (k.le.lstot(id)) p0a(k) = y(knsp(k,id))
c                                 convert the dependent species to
c                                 idependent species
         do l = 1, ndep(id)
            p0a(k) = p0a(k) + y2pg(l,k,id) * y(knsp(lstot(id)+l,id))
         end do 

         pa(k) = p0a(k)

      end do

      if (ksmod(id).eq.9.or.ksmod(id).eq.10) then
c                                 orphan vertex model
         do k = 1, nstot(id)
c                                 renormalize the prismatic endmembers
            p0a(k) = p0a(k)*x(ostg(id),ispg(id,ostg(id)))
            pa(k) = p0a(k)
         end do 

         if (zbad(pa,id)) then 
         write (*,*) 'z2 wonka wonaa'
         end if

         do k = 1, ndim(ostg(id),id)
c                                 add the orphan fractions
            l = lstot(id) - ndim(ostg(id),id) + k 
            p0a(l) = p0a(l) + x(ostg(id),k)
            pa(l) = p0a(l)
         end do

         if (zbad(pa,id)) then 
         write (*,*) 'z3 wonka wonaa'
         end if 

      end if 

      end 

      subroutine p0dord (id)
c-----------------------------------------------------------------------
c decomposes p0a values that specify the abundances of independent 
c ordered endmembers to their stoichiometric equivalent disordered
c species. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,k,l,ind 

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c-----------------------------------------------------------------------

      do k = 1, nord(id)
         do l = 1, nrct(k,id)
            ind = ideps(l,k,id)
            p0a(ind) = p0a(ind) + dvnu(ind,k,id) * p0a(lstot(id)+k)
         end do 
      end do 

      end 

      subroutine specis (g,id)
c----------------------------------------------------------------------
c subroutine to speciation of a solution with disordered composition p0a. 
c the speciated composition is returned in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id

      logical error

      double precision g, gdord, omega, gex

      external omega, gex

      double precision enth
      common/ cxt35 /enth(j3)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
      g = 0d0

      if (ksmod(id).eq.27) then 
c                                 green et al 2016 melt model,
c                                 special case because of non-equimolar
c                                 speciation reaction.
         if (t.lt.nopt(20)) then
            g = 1d6 
         else  
            call gpmelt (g,id)
         end if 

         return

      else if (.not.lrecip(id)) then
c                                 non-reciprocal, initialize p0/pa
         do i = lstot(id)+1, nstot(id)
            p0a(i) = 0d0
            pa(i) = p0a(i)
         end do

      end if
c                                 initialize limit expressions
      call p0limt (id)
c                                 as most models are single species and
c                                 there is so much overhead in computing
c                                 multiple speciation, use a special routine
c                                 for single species models:
      if (nord(id).gt.1) then
 
         call speci2 (g,id,error)

      else 

         call speci1 (g,id,1,error)

      end if 

      if (error.or.iopt(17).ne.0) then 
c                                 if speciation returns error, or order_check is on,
c                                 i.e., iopt(17).ne.0, compute disordered g.
         gdord =  gex(id,p0a) - t*omega(id,p0a)
 
         if (lrecip(id)) then 
            do i = 1, nord(id)
               gdord = gdord + p0a(lstot(id)+i)*enth(i)
            end do 
         end if 

         if (gdord.lt.g) g = gdord

      end if 
c                                 convert the ordered endmember fractions to 
c                                 disordered fractions (stored in the p0a array).
      if (lrecip(id)) call p0dord (id)

      end

      subroutine oenth (id)
c----------------------------------------------------------------------
c subroutine to the enthalpy of ordering for speciation models
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,id

      double precision enth
      common/ cxt35 /enth(j3)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      do i = 1, nord(id) 
         enth(i) = deph(1,i,id) + t * deph(2,i,id) + p * deph(3,i,id) 
      end do 

      end 

      subroutine gderiv (id,g,dp,error)
c----------------------------------------------------------------------
c subroutine to compute the g of a solution (id) and it's 1st and 2nd 
c derivatives with respect to the oncentrations of nord(id) ordered 
c species. the formulation assumes atomic site fractions are linear 
c functions of the ordered species concentrations (p's) and that the 
c excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer i,k,l,i1,i2,id,norder,ipvt(j3)

      double precision g,dp(j3),t,s,ds(j3),d2s(j3,j3),dg(j3),d2g(j3,j3)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision p,tk,xc,u1,u2,tr,pr,r,ps
      common/ cst5 /p,tk,xc,u1,u2,tr,pr,r,ps

      double precision enth
      common/ cxt35 /enth(j3)

      logical pin
      common/ cyt2 /pin(j3)
c----------------------------------------------------------------------
c                                 initialize, d2gx has been set in setw
      g = 0d0

      norder = nord(id)

      do k = 1, norder
         dg(k) = 0d0
         do l = k, norder
            d2g(l,k) = d2gx(l,k)
         end do 
      end do 

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 assuming regular terms
            i1 = jsub(1,i,id)
            i2 = jsub(2,i,id)

            g = g + w(i) * pa(i1) * pa(i2)

            do k = 1, norder

               if (.not.pin(k)) cycle

               dg(k) = dg(k) + w(i) * (pa(i1)*dydy(i2,k,id) 
     *                               + pa(i2)*dydy(i1,k,id))
            end do 

         end do  
c                                 get derivative of excess function
         if (llaar(id)) then 
c                                 so far this is unessecary because
c                                 t does not vary in h&p ordering models.
            t = 0d0 
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do 
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions 
c                                 because of the "tphi" term in the 
c                                 van laar.

            do k = 1, norder

               if (.not.pin(k)) cycle
c                                 convert dg and d2g to the full derivative
               dg(k) = (dg(k) - g*dt(k)/t)/t
               do l = k, norder
                  d2g(l,k) = (d2g(l,k) - 2d0*dt(k)*dg(k))/t
               end do 
            end do 
c                                 and the full excess energy 
            g = g/t

         end if 

      end if 
c                                 get the delta configurational entropy and derivatives
      call sderiv (id,s,ds,d2s)

      do k = 1, norder

         if (.not.pin(k)) cycle

         g = g + enth(k) * pa(lstot(id)+k)
c                                 dg is the negative of the differential of g 
c                                 with respect to the kth species.
         dg(k) = -(enth(k) + dg(k) - tk*ds(k))
         do l = k, norder 
            d2g(l,k) = d2g(l,k) - tk*d2s(l,k)
         end do 
      end do
c                                 determininats, to check for a saddle point
c      if (norder.eq.2) then 
c         detg = d2g(1,1)*d2g(2,2)-d2g(2,1)**2
c      else 
c         detg = d2g(1,1)*(d2g(2,2)*d2g(3,3)-d2g(3,2)**2) 
c     *        - d2g(2,2)*d2g(1,3)**2
c     *        + 2d0*d2g(2,1)*d2g(3,2)*d2g(3,1)-d2g(2,1)**2*d2g(3,3)
c      end if 

      g = g - tk*s 
c                                 copy dg and d2g into dp and d2s
      do k = 1, norder
         if (pin(k)) then 
            dp(k) = dg(k)
            d2s(k,k) = d2g(k,k)
            do l = k+1, norder
               if (pin(l)) then  
                  d2s(l,k) = d2g(l,k)
                  d2s(k,l) = d2g(l,k) 
               end if 
            end do 
         end if 
      end do 

      do k = 1, norder
         if (.not.pin(k)) then
            dp(k) = 1d0 
            d2s(k,k) = 1d0 
            do l = 1, norder
               if (l.eq.k) cycle 
               d2s(l,k) = 0d0
               d2s(k,l) = 0d0 
            end do 
         end if 
      end do     
c                                 get the newton-raphson increments:
c                                 this is a general factorization routine, should
c                                 exploit that d2g is symmetric.
      call factr2 (d2s,j3,norder,ipvt,error)
c                                 solve for the increments by back-substitution, 
c                                 this routine is also not efficient and should 
c                                 be re written. 
      if (.not.error) call subst2 (d2s,ipvt,j3,norder,dp,error) 
c                                 substitute replaces the values of dg with the 
c                                 newton-raphson increments for the ordered species
c                                 compositions. 
      end

      subroutine sderiv (id,s,dsy,dsyy)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of a dependent species.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      logical inf

      double precision zt,dzdy,s,dsy(j3),dsyy(j3,j3),q,zl,
     *                 z(m11,m10),s0,ztemp,zlnz,
     *                 dsinf(j3),d2sinf(j3,j3)
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical badend, sck, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),badend(m4,h9),sck(h9),nrf(h9)
c DEBUG
      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      s = 0d0
c                                 for each site
      do i = 1, msite(id)

         zt = 0d0
         s0 = zt 
c                                 get site fractions
         do j = 1, ksp(i,id)

            ztemp = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               ztemp = ztemp + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do  

            zt = zt + ztemp
            if (ztemp.gt.0d0) s0 = s0 - ztemp * dlog (ztemp)
            z(j,i) = ztemp

         end do 

         ztemp = 1d0 - zt
         if (ztemp.gt.0d0) s0 = s0 - ztemp * dlog (ztemp)
         z(j,i) = ztemp
         s = s + qmult(i,id) * s0

      end do 
c                                 initialize derivatives:
      inf = .false.

      do k = 1, nord(id)
         dsy(k) = 0d0 
         dsinf(k) = 0d0 
         do l = k, nord(id) 
            dsyy(l,k) = 0d0 
            d2sinf(l,k) = 0d0 
         end do 
      end do  
c                                 evaluate derivatives:
      do i = 1, msite(id)

         q = qmult(i,id)
 
         do j = 1, ksp(i,id) + 1   

            zl = z(j,i) 

            if (zl.gt.0d0) then 
               zlnz = 1d0 + dlog(zl)
            else
               zlnz = 1d0 
            end if 

            do k = 1, nord(id)        
c                                 skip species not in the model
               if (.not.pin(k)) cycle
c                                 sdzdp is (dz(i,j)/dp(k))
               dzdy = sdzdp(k,j,i,id)

               if (dzdy.eq.0d0) cycle 

               if (zl.gt.0d0) then 
c                                 the first derivative is
                  dsy(k) = dsy(k) - q * dzdy * zlnz
c                                 and the jacobians are
                  do l = k, nord(id)

                     if (.not.pin(l)) cycle 
                     dsyy(l,k) = dsyy(l,k) 
     *                         - q * dzdy * sdzdp(l,j,i,id) / zl
                  end do 

               else 

                   if (zl.lt.-nopt(5)) then 
                      write (*,*) 'wacka boom',zl,j,i
                      write (*,*) (p0a(l),l=1,nstot(id))
                      write (*,*) (pa(l),l=1,nstot(id))
                      write (*,*) p, t
                      write (*,*) 'please report this error'
                      write (*,*) 
                      inf = .true.

                   end if 
c                                 a species with a non-zero
c                                 derivative is zero, the 
c                                 first will be sign(dzdy)*infinity
                  dsinf(k) = dsinf(k) + dsign(q,dzdy)

                  do l = k, nord(id)
c                                 the 2nd will be -sign of 
c                                 cross term * infinity
                     if (.not.pin(l)) cycle 
                     d2sinf(l,k) = dsinf(k) -  
     *                             dsign(q,dzdy * sdzdp(l,j,i,id))
                  end do 
               
               end if 
            end do  

         end do 

      end do 

      if (inf) then

         do k = 1, nord(id)

            if (.not.pin(k)) cycle 
            if (dabs(dsinf(k)).gt.zero) dsy(k) = 1d4*dsinf(k)

            do l = k, nord(id)
               if (.not.pin(l)) cycle 
               if (dabs(d2sinf(l,k)).gt.zero) 
     *                                  dsyy(l,k) = 1d5*d2sinf(l,k)
            end do  
 
         end do 

      end if
c                                 endmember corrections
      do i = 1, nstot(id)

         s = s - pa(i) * scoef(i,id)

         do k = 1, nord(id) 
           dsy(k) = dsy(k) - dydy(i,k,id) * scoef(i,id)
         end do 

      end do

      end

      subroutine factr2 (a,m,n,ipvt,error)
c-----------------------------------------------------------------------
c factr2 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the actual dimension of matrix a
c           m- the phsical dimension of matrix a
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c       error- false if a is of rank = n, and true if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      logical error 

      include 'perplex_parameters.h'

      integer m,ipvt(m),i,j,k,ip1,n,istr
    
      double precision a(m,m),d(m),rmax,tmax,temp,ratio

      double precision wmach(9)
      common /ax02za/wmach
c-----------------------------------------------------------------------
      error = .false.
c                            initialize ipvt,d
      do i = 1, n

         ipvt(i) = i
         rmax = 0d0

         do j = 1, n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.wmach(3)) then 
            error = .true.
            return 
         end if 

         d(i) = rmax

      end do 
c                            begin decomposition: 
      do i = 1, n-1
c                            determine pivot row (istr).
         ip1 = i+1
         rmax = dabs(a(i,i))/d(i)
         istr = i

         do j = ip1, n

            tmax = dabs(a(j,i))/d(j)

            if (tmax.gt.rmax) then
               rmax = tmax
               istr = j
            end if 

         end do

         if (dabs(rmax).lt.wmach(3)) then 
            error = .true.
            return
         end if 
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 

            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = d(istr)
            d(istr) = d(i)
            d(i) = temp

            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 

         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1, n

            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)

            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 

         end do

      end do 

      if (dabs(a(n,n)).lt.wmach(3)) error = .true.
 
      end

      subroutine subst2 (a,ipvt,m,n,b,error)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the actual dimension of the matrix a.
c           m- the physical dimension of a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error 
 
      integer m,ipvt(m),ip,i,j,n,ii

      double precision a(m,m),b(m),x(m),sum
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j) + sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous. should check
c                                 what's with factor. 
         error = .true.
         return
      end if 

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do 

         if (a(i,i).eq.0d0) then
c                                 as above.
            error = .true.
            return
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 

      b(n) = x(n)

      end

      subroutine speci0 (g,h,w,n,fac,c0,f)
c----------------------------------------------------------------------
c subroutine to solve speciation of a 0-d solution with 1 ordering parameter
c by halving. assumes an ordered species in which A is on 1 site and B is on 
c n sites, and a disordered state in which A and B are distributed over all 
c n+1 sites. 

c    h   - is the enthalpy of complete disordering
c    w   - is the interaction energy
c    fac - is an empirical correction to the entropy, supposedly accounting for SRO.
c    g   - is the change in G for the stable speciation relative to the ordered state. 
c    y   - is the fraction of the ordered species

c                                                  JADC, Aug 29, 2010. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision g,h,w,sign,dy,odg,ndg,n,f,y,dgdy,rt,c0,c1,c2,fac

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 check ordered state
      y = 1d0 - nopt(5)
      rt = r*t*fac

      odg = dgdy(h,w,n,f,y,rt)
c                                 if dgdy < 0 must be fully ordered
c                                 (not really, non-zero w could make
c                                 a zero at intermediate y).
      if (odg.lt.0d0) then 
c                                 fully ordered (y=1) 
         y = 1d0 

      else 
c                                 initialize at halfway point
         dy = -0.5d0 
c                                 iteration loop:
         do 

            y = y + dy 
            if (y.le.0d0) y = nopt(5)

            ndg = dgdy(h,w,n,f,y,rt)

            sign = odg*ndg

            if (sign.lt.0d0) then 
c                                 crossed the zero, flip the search
               odg = ndg
               dy = -dy/2d0

            else if (dabs(dy).lt.nopt(5)) then
c                                 refined to tolerance
               exit

            else if (y.le.nopt(5)) then 
c                                 fully disordered, y=0, c1 = c2
               y = 0d0
               exit 

            end if 

         end do 

      end if 

      c1 = (n+y)/c0

      if (c1.lt.1d0-nopt(5).and.c1.gt.nopt(5)) then 
         g = rt*n*(c1*dlog(c1)+(1d0-c1)*dlog(1d0-c1))
      else  
         g = 0d0 
      end if 

      c2 = (1d0-y)*n/c0

      if (c2.lt.1d0-nopt(5).and.c2.gt.nopt(5)) 
     *   g = g + rt*(c2*dlog(c2) + (1d0-c2)*dlog(1d0-c2))

      g = g + (1d0-y)*( w*y + h)
      
      end 

      double precision function dgdy (h,w,n,f,y,rt)
c----------------------------------------------------------------------
c function to compute dg/dy for subroutine speci0
c----------------------------------------------------------------------
      implicit none

      double precision h,w,n,f,y,rt

      dgdy = (1d0-2d0*y)*w - h
     *       - rt*f*dlog(n*(1d0-y)**2/(n+y)/(1d0+n*y))

      end 

      subroutine speci1 (g,id,k,error)
c----------------------------------------------------------------------
c subroutine to speciation of a solution with a single ordering parameter
c composition is returned 
c in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, jd, k, itic, nr, ind(m4)

      logical error, done

      double precision g,pmax,pmin,dp,dpmax,omega,gex,dy(m4)

      external gex, omega

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)
c----------------------------------------------------------------------
c                                 number of reactants to form ordered species k
      nr = nrct(k,id)

      do i = 1, nrct(k,id)
c                                 dependent disordered species
         ind(i) = ideps(i,k,id)
c                                 stoichiometric coefficients
         dy(i) = dydy(ind(i),k,id)
      end do 

      jd = lstot(id) + k

      error = .false.
c                                 starting point
      if (lrecip(id)) then 
c                                 reciprocal
         call plimit (pmin,pmax,k,id) 
         dpmax = pmax - pmin

      else 
c                                 find the maximum proportion of the 
c                                 ordered species cannot be > the amount 
c                                 of reactant initially present.
         dpmax = 1d0
         do i = 1, nr
            if (-p0a(ind(i))/dy(i).lt.dpmax) dpmax = -p0a(ind(i))/dy(i)
         end do 

      end if 
c                                 to avoid singularity set the initial 
c                                 composition to the max - nopt(5), at this
c                                 condition the first derivative < 0, 
c                                 and the second derivative > 0 (otherwise
c                                 the root must lie at p > pmax - nopt(5).               
      pin(k) = .true.
      dp = dpmax - nopt(5)
      pmax = p0a(jd) + dp
      pmin = p0a(jd) + nopt(5)

      if (pmax-pmin.lt.nopt(5)) return 
c                                 get starting end for the search
c                                 first try the maximum
      call pincs (dp,dy,ind,jd,nr)

      call gderi1 (k,id,dp)

      if (dp.ge.0d0) then 
c                                 at the maximum concentration
c                                 and the increment is positive,
c                                 the solution is fully ordered
c                                 or a local minimum, try the 
c                                 the disordered case:
         call pincs (pmin,dy,ind,jd,nr)

         call gderi1 (k,id,dp)

         if (dp.le.0d0) then 
c                                 neither min nor max starting point
c                                 is possible. setting error to 
c                                 true will cause specis to compare 
c                                 the min/max order cases, specis 
c                                 computes the min case g, therefore
c                                 the case is set to max order here:
            error = .true.

         end if

      end if 

      if (.not.error) then 
c                                 increment and check p
         call pcheck (pa(jd),pmin,pmax,dp,done)
c                                 set speciation
         call pincs (pa(jd)-p0a(jd),dy,ind,jd,nr)
c                                 iteration counter
         itic = 0
c                                 newton raphson iteration
         do 

            call gderi1 (k,id,dp)

            call pcheck (pa(jd),pmin,pmax,dp,done)
c                                 done means the search hit a limit
c                                 or dp < tolerance.
            if (done) then

               goodc(1) = goodc(1) + 1d0
               jcount(1) = jcount(1) + itic 
c                                 use the last increment
               call pincs (pa(jd)-p0a(jd),dy,ind,jd,nr)

               exit

            else 

               itic = itic + 1
c                                 apply the increment
               call pincs (pa(jd)-p0a(jd),dy,ind,jd,nr)

               if (itic.le.iopt(21)) cycle
c                                 failed to converge. exit
               error = .true.
               badc(1) = badc(1) + 1
               jcount(1) = jcount(1) + itic 

               exit

            end if

         end do

      end if
c                                 didn't converge or couldn't 
c                                 find a starting point, set 
c                                 ordered speciation, specis will
c                                 compare this the disordered case.
      if (error) call pincs (dpmax,dy,ind,jd,nr)

      g = pa(jd)*enth(k) - t*omega(id,pa) + gex(id,pa)

      end

      subroutine speci2 (g,id,error)
c----------------------------------------------------------------------
c subroutine to multiple speciation of a solution with disordered composition 
c p0a. the speciated composition is returned in array pa. 
c    id identifies the solution.
c    g is the change in G for the stable speciation relative to a mechanical
c      mixture of the endmembers.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error

      integer i,k,id,lord,itic

      double precision g,dp(j3),tdp,gold,xtdp

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)
c----------------------------------------------------------------------
c                                 get initial p values
      call pinc0 (id,lord)
c                                 lord is the number of possible species
      if (lord.eq.1) then 

         do i = 1, nord(id)
            if (pin(i)) then 
               call speci1 (g,id,i,error)
               exit 
            end if 
         end do

      else if (lord.gt.1) then 

         itic = 0  
         gold = 0d0
         xtdp = 0d0

         do 

            call gderiv (id,g,dp,error)

            if (error) then
               badc(1) = badc(1) + 1
               exit
            end if 

            tdp = 0d0 

            do k = 1, nord(id)
               
               if (.not.pin(k)) cycle

               call pinc (dp(k),k,id)

               tdp = tdp + dabs(dp(k))

            end do 
c                                 nov 23, 2016 added exit if diverging 
c                                 g > gold, itic > 2
            if (tdp.lt.nopt(5).or.
     *          dabs((gold-g)/g).lt.nopt(5).or.tdp.eq.xtdp.or.
     *          itic.gt.2.and.gold.le.g) then

               goodc(1) = goodc(1) + 1d0
               jcount(1) = jcount(1) + itic
               exit

            end if 
c                                 before nov 23, 2016, xtdp was only set 
c                                 if g > gold (i.e., diverging).
            xtdp = tdp

            gold = g 

            itic = itic + 1

            if (itic.gt.iopt(21)) then
c                                 not converging, under the assumption that 
c                                 this happens at low T use pinc0 to set an ordered
c                                 composition and exit
               badc(1) = badc(1) + 1
               jcount(1) = jcount(1) + itic 
               error = .false.
               call pinc0 (id,lord)
               exit 

            end if 

         end do

      else 
c                                 no speciation possible, but still need
c                                 to calculate g (setting error will do this).
c                                 set g to a large number to assure specis
c                                 selects the disordered configuration.
         g = 1d9
         error = .true.

      end if

      end


      subroutine pincs (dp,dy,ind,jd,nr)
c----------------------------------------------------------------------
c subroutine to increment the proportions of endmembers involved in a
c predefined ordering reaction (called only by speci1, see pinc and speci2
c for the general case).

c this routine replicates dpinc, but the p's are computed from p0 and
c uses simple arrays. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jd, nr, ind(m4)

      double precision dp, dy(m4)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c----------------------------------------------------------------------
      pa(jd) = p0a(jd) + dp

      do i = 1, nr
         pa(ind(i)) = p0a(ind(i)) + dy(i)*dp
      end do

      end 

      subroutine pinc (dp,k,id)
c----------------------------------------------------------------------
c subroutine to increment the k'th species of solution id, if the increment
c violates a stoichiometric limit, it's set to half it's maximum value.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k,id,jd

      double precision dp,pmx,pmn

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 given dp check if it violates
c                                 stoichiometric constraints
      jd = lstot(id) + k 

      call plimit (pmn,pmx,k,id)       

      if (pa(jd)+dp.gt.pmx-nopt(5)) then 
         dp = pmx - pa(jd) - nopt(5)
      else if (pa(jd)+dp.lt.pmn+nopt(5)) then 
         dp = pmn - pa(jd) + nopt(5)
      end if  
c                                 adjust the composition by the increment
      call dpinc (dp,k,id,jd)

      end 

      subroutine dpinc (dp,k,id,jd)
c----------------------------------------------------------------------
c subroutine to increment the k'th species of solution id, if the increment
c violates a stoichiometric limit, it's set to half it's maximum value.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,k,id,jd

      double precision dp

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c----------------------------------------------------------------------
c                                 adjust the composition by the increment
      do i = 1, nrct(k,id)

         pa(ideps(i,k,id)) = pa(ideps(i,k,id)) 
     *                     + dydy(ideps(i,k,id),k,id)*dp

      end do

      pa(jd) = pa(jd) + dp

      end 

      subroutine pinc0 (id,lord)
c----------------------------------------------------------------------
c subroutine set initial species concentrations to half their 
c stoichiometric limit.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id,jd,lord,iout,ibad(m4)

      double precision dp,pmn,pmx,dpp(j3),dinc,tinc
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      lord = 0

      if (icase(id).eq.1) then 
c                                 case 1: fully correlated
         dinc = 0.9d0/dfloat(nord(id))
         tinc = dinc

         do k = 1, nord(id)

            call plimit (pmn,pmx,k,id) 

            if (pmn.ge.pmx.or.pmx-pmn.lt.nopt(5)) then 
               pin(k) = .false.
               cycle 
            else 
               pin(k) = .true.
               lord = lord + 1
            end if 

            jd = lstot(id) + k 

            dp = pmn + (pmx - pmn) * tinc - pa(jd)
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

            tinc = tinc + dinc

         end do

      else if (icase(id).eq.2.or.icase(id).eq.0) then 

         if (icase(id).eq.2) then 
c                                 case 2: positive partial correlation
c                                         iterate 4 times for increments
            iout = 5
         else 
c                                 case 0: no correlation/iteration
            iout = 1
         end if

         do i = 1, iout

            do k = 1, nord(id)

               call plimit (pmn,pmx,k,id) 

               if (i.eq.1) then 

                  if (pmn.ge.pmx.or.pmx-pmn.lt.nopt(5)) then 
                     pin(k) = .false.
                     cycle 
                  else 
                     pin(k) = .true.
                     lord = lord + 1
                 end if 

               end if 
c                                 adjust the composition by the first increment
               jd = lstot(id) + k 
               dp = pmx - pa(jd) 
               pa(jd) = pa(jd) + dp
               dpp(k) = pa(jd) - p0a(jd)

            end do 
c                                 no species possible
            if (lord.eq.0) return

         end do 
c                                 back off from maximum for final assignements
         do k = 1, nord(id)

            if (.not.pin(k)) cycle

            jd = lstot(id) + k 
            pa(jd) = p0a(jd)

            dp = dpp(k)*0.9d0
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

         end do 

      else if (nord(id).eq.1) then
c                                 only one order parameter, as currently programmed
c                                 this will never be called. 
         call plimit (pmn,pmx,1,id) 

         if (pmn.ge.pmx) then 
            pin(1) = .false.
         else 

            pin(1) = .true.
            lord = 1
            jd = lstot(id) + 1

            dp = pmn + (pmx - pmn) * 0.9d0 - pa(jd)
c                                 adjust the composition by the first increment
            call dpinc (dp,k,id,jd)

         end if 

      else 
c                                 unanticipated case?
         call error (999,p0a(1),i,
     *              'unanticipated correlation between ordered species')

      end if
c                                 check for degenerate compositions
      if (lord.gt.0) then 

         iout = 0 

         do i = 1, lstot(id)
            if (p0a(i).eq.0d0) then
               iout = iout + 1
               ibad(iout) = i 
            end if
         end do 
c                                 the indices of the present components are igood(1..in)
         if (iout.gt.0) then 
            do k = 1, nord(id)
               if (pin(k)) then    
c                                 check that the ordered species are in the subcomposition
                 do j = 1, nrct(k,id)                                    
                     do i = 1, iout
                        if (ideps(j,k,id).eq.ibad(i)) then 
                           lord = 0 
                           return 
                        end if 
                     end do 
                  end do 
               end if
            end do 
         end if 
      end if  

      end 

      subroutine pcheck (x,xmin,xmax,dx,quit)
c-----------------------------------------------------------------------
c subroutine to increment x for a 1-d root search between xmin and xmax.
c modified to redefine limits according to gradient, Dec 20, 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit 

      double precision x, xmin, xmax, dx, xt

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      quit = .false.

      xt = x + dx

      if (xt.eq.xmin.or.xt.eq.xmax) then 
c                                 hit the limit, don't set x to 
c                                 the limit to save revaluating x
c                                 dependent variables.

        write (*,*) 'this should not happen!!'

        x = xt
        quit = .true.

        return

      end if 

      if (dx.lt.0d0) then 
c                                 narrow limit to avoid oscillating
         if (x.lt.xmax) xmax = x
c                                 increment would make x < xmin
c                                 revise the increment
         if (xt.lt.xmin) dx = (xmin - x)/2d0 

      else if (dx.gt.0d0) then 
c                                 narrow limit to avoid oscillating
         if (x.gt.xmin) xmin = x
c                                 increment would make x > xmax
c                                 revise the increment 
         if (xt.gt.xmax) dx = (xmax - x)/2d0

      end if

      x = x + dx
c                                 check if dx has dropped below 
c                                 threshold for convergence
      if (dabs(dx).lt.nopt(5)) quit = .true.

      end 

      subroutine gderi1 (k,id,dg)
c----------------------------------------------------------------------
c subroutine computes the newton raphson increment dg from the 1st and 2nd 
c derivatives of the g of solution (id) with respect to the concentration 
c of the kth ordered species. 

c on input dg is the last increment. 

c the formulation assumes:

c  1) the speciation reaction is equimolar (see gpmelt for non-equimolar 
c     case.

c  2) atomic site fractions are linear functions of the ordered species 
c     concentrations (p's).

c  3) the excess function is second order.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,k,i1,i2,id

      logical inf 

      double precision g,dg,d2g,t,ds,d2s,dp
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision alpha,dt
      common/ cyt0  /alpha(m4),dt(j3)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      double precision enth
      common/ cxt35 /enth(j3)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize, d2gx has been set in setw
      g = 0d0
      dp = dg 
      dg = g
      d2g = d2gx(k,k)

      if (lexces(id)) then

         do i = 1, jterm(id)
c                                 assuming regular terms
           i1 = jsub(1,i,id)
           i2 = jsub(2,i,id)

           g = g + w(i) * pa(i1) * pa(i2)
           dg = dg + w(i) * (pa(i1)*dydy(i2,k,id) 
     *                     + pa(i2)*dydy(i1,k,id))

         end do  
c                                 get derivative of excess function
         if (llaar(id)) then 
c                                 for h&p van laar, this is unnecessary because
c                                 t is constant.
            t = 0d0 
c                                 h&p van laar
            do i = 1, nstot(id)
               t = t + alpha(i)* pa(i)
            end do 
c                                 coming out of this loop g, dg, and
c                                 d2g  are not the complete functions 
c                                 because of the "tphi" term in the 
c                                 van laar.

c                                 convert dg and d2g to the full derivative
            dg = (dg - g*dt(k)/t)/t
            d2g = (d2g - 2d0*dt(k)*dg)/t

         end if 

      end if 
c                                 get the configurational entropy derivatives
      call sderi1 (k,id,ds,d2s,inf)

      dg  = dg + enth(k)  - v(2)*ds
      d2g = d2g - v(2)*d2s
c                                 dg becomes the newton raphson increment
      if (inf) then 
         dg = dsign(dp/2d0,-dg/d2g)
      else
         dg = -dg/d2g
      end if  

      end

      subroutine sderi1 (l,id,ds,d2s,inf)
c----------------------------------------------------------------------
c subroutine to the derivative of the configurational entropy of a 
c solution with respect to the proportion of the lth ordered species.

c ASSUMES ORDERED SPECIES HAVE NO CONFIGURATIONAL ENTROPY!
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,id

      logical inf

      double precision zt,dzdy,dzy,dzyy,zl,ds,d2s,zlnz,dsinf
c                                 working arrays
      double precision zz, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),zz(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 configurational entropy variables:
      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      double precision dppp,d2gx,sdzdp
      common/ cxt28 /dppp(j3,j3,m1,h9),d2gx(j3,j3),sdzdp(j3,m11,m10,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)
c----------------------------------------------------------------------

      inf = .false.
      ds = 0d0 
      d2s = 0d0

      do i = 1, msite(id)

         dzy = 0d0  
         dzyy = 0d0
 
         zt = 0d0 
         dsinf = 0d0 
  
         do j = 1, ksp(i,id)    

            zl = d0(j,i,id)
c                                 for each term:
            do k = 1, lterm(j,i,id)
               zl = zl + dcoef(k,j,i,id) * pa(ksub(k,j,i,id))
            end do 

            if (zl.gt.0d0) then 
               zt = zt + zl 
               zlnz = 1d0 + dlog(zl)
            else
               zlnz = 1d0 
            end if         
c                                 sdzdp is (dz(i,j)/dp(l))
            dzdy = sdzdp(l,j,i,id)

            if (dzdy.eq.0d0) cycle 

            if (zl.gt.0d0) then 
c                                 the first derivative is
               dzy = dzy - dzdy * zlnz
c                                 and the jacobians are

               dzyy = dzyy  - dzdy**2 / zl

            else 
c                                 a species with a non-zero
c                                 derivative is zero, the s
c                                 derivative may be +/-infinite
               dsinf = dsinf + dsign(1d0,dzdy)
                             
            end if 

         end do 
c                                 add the contibution from the ksp(i,id)+1th
c                                 species:
         zl = 1d0 - zt

         if (zl.gt.0d0) then 
            zlnz = 1d0 + dlog(zl)
         else
            zlnz = 1d0 
         end if 

         dzdy = sdzdp(l,j,i,id)

         if (dzdy.ne.0d0) then 
            if (zl.gt.0d0) then 
c                                 the first derivative is
               dzy = dzy - dzdy * zlnz
c                                 and the second is 
               dzyy = dzyy  - dzdy**2 / zl

            else 
c                                 
c                                 a species with a non-zero
c                                 derivative is zero, the s
c                                 derivative may be +/-infinite
               dsinf = dsinf + dsign(1d0,dzdy) 
               
            end if 

         end if 

         if (dabs(dsinf).lt.zero) then 
            ds = ds + qmult(i,id)*dzy
            d2s = d2s + qmult(i,id)*dzyy
         else 
            inf = .true.
            ds = ds + qmult(i,id)*dsinf*1d4
            d2s = d2s - qmult(i,id)*dabs(dsinf)*1d5
         end if 

      end do
c                                 for models with disordered 
c                                 endmembers, correct first derivative for the
c                                 change in endmember configurational entropy
      do i = 1, nstot(id)
         ds = ds - dydy(i,l,id)*scoef(i,id)
      end do 

      end

      subroutine p0limt (id)
c----------------------------------------------------------------------
c subroutine to compute the sums of the p0 terms in ordered species 
c limit expressions for solution id.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      do k = 1, nord(id) 
c                                 for ordered species k
         do i = 1, ln(k,id)
c                                 for limit i
            tsum(i,k) = l0c(1,i,k,id)  

            do j = 1, lt(i,k,id)
c                                 for term j
               tsum(i,k) = tsum(i,k) + lc(j,i,k,id)*p0a(lid(j,i,k,id))

            end do
 
         end do 

      end do 

      end  

      subroutine plimit (pmn,pmx,k,id)
c----------------------------------------------------------------------
c subroutine to compute minimum and maximum concentration of ordered
c species k in solution id from site fraction constraints, assumes the
c p0 terms have been accumulated in tsum (routine p0limt)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision pmn,pmx,mini
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)
c----------------------------------------------------------------------
      pmx = 1d99
      pmn = -1d99

      do i = 1, ln(k,id)

         mini =  tsum(i,k) 

         do j = 1, jt(i,k,id)

            mini = mini + jc(j,i,k,id)*pa(jid(j,i,k,id))

         end do 

         if (mini.gt.pmn) pmn = mini
         if (l0c(2,i,k,id)+mini.lt.pmx) pmx = mini + l0c(2,i,k,id)

      end do 

      end  

      subroutine readlm (idim,tname,bad)
c---------------------------------------------------------------------
c readlm - reads stoichiometric limits on ordered species concentrations
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer j,k,l,jd,len,inds(k7),ier,ict, idim

      double precision coeffs(k7)

      logical bad

      character begin*5, tag*3, tname*10

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(lchar)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer ostot
      common/ junk /ostot

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                               initialize limit counter
      do j = 1, norder
         limn(j) = 0
      end do 

      if (jsmod.eq.27) return

      call readcd (n9,len,ier,.true.)

      write (begin,'(5a)') (chars(j),j=1, 5)

      if (begin.eq.'begin') then 

         do 
c                                 read the limit equations for the 
c                                 amounts of the ordered endmembers
            call readz (coeffs,inds,ict,idim,tname,tag)

            if (tag.eq.'end') then 
               exit 
c                                 commented to allow fixed bounds
c           else if (ict.eq.1) then
c               bad = .true.
c               exit
            end if 
c                                 convert the endmember index to the 
c                                 ordered species index
            jd = inds(1) - ostot
            limn(jd) = limn(jd) + 1
            k = limn(jd) 
            if (k.gt.j5) call error (999,coeffs(1),k,tname)
                   
            limt(k,jd) = 0
            jimt(k,jd) = 0 
c DEBUG DEBUG just patched in ostot
            do l = 2, ict
c                                 four cases
               if (    (inds(l).le.ostot) 
     *             .or.(inds(l) - ostot.eq.jd.and.
     *                  inds(l+1).ne.inds(l))
     *             .or.(inds(l).eq.inds(l-1).and.
     *                  inds(l).gt.ostot.and.
     *                  inds(l) - ostot.ne.jd)) then
c                                 1) disordered species go in p0 array
c                                 2) the species is limit species, the 
c                                    species must be in the p0 array.
c                                 3) the species is an ordered species, 
c                                    but not the limit species,
c                                    and it's the second occurence
                  limt(k,jd) = limt(k,jd) + 1
                  j = limt(k,jd)
                  if (j.gt.j6) call error (33,coeffs(1),j,tname)
                  limid(j,k,jd) = inds(l)
                  limc(j,k,jd) = coeffs(l)

               else if (inds(l).gt.ostot.and.
     *                  inds(l) - ostot.ne.jd) then
c                                 4) is an ordered species p-term
c                                    if the coefficient is zero it's just
c                                    a place holder and can be skipped.
                  if (coeffs(l).eq.0d0) cycle
                  jimt(k,jd) = jimt(k,jd) + 1
                  j = jimt(k,jd)
                  if (j.gt.j3) call error (33,coeffs(1),j,tname)
                  jimid(j,k,jd) = inds(l)
                  jimc(j,k,jd) = coeffs(l)

               else 

                  bad = .true.

               end if 

            end do 
c                                 the constant and delta (max-min) are:
            j = limt(k,jd) + 1
c                                 set number of terms to negative if constant bounds
            if (limt(k,jd).eq.0) then
               limt(k,jd) = -1 
               jimt(k,jd) = -1
            end if 

            limc(j  ,k,jd) = coeffs(1)
            limc(j+1,k,jd) = coeffs(ict+1)

         end do 

      else if (jsmod.eq.6.and.norder.eq.1) then 

         backspace (n9)

      else 

         bad = .true.

      end if 

      if (bad) then 
         if (iam.lt.3.or.iam.eq.15) then 
            write (*,1000) tname,(chars(j),j=1,len)
            write (*,1010)
         end if  
         backspace (n9)
      end if 

1000  format ('**warning ver203** READLM missing or invalid format for '
     *       ,'stoichiometric limit of ordered species',/,'currently ',
     *        'reading (and rejecting) solution model: ',a,
     *      /,'last record was:',/,400a)
1010  format (/,'This error may be due to an out-of-date '
     *         ,'solution model file.',/
     *         ,'The current version is: '
     *         ,'www.perplex.ethz.ch/perplex/datafiles/solution_model'
     *         ,'.dat',/)

      end 

      subroutine input9 (first,output)
c-----------------------------------------------------------------------
c given a list of solution phase names (fname(h9)) input9 searches a
c data file (on unit n9) for the relevant data and subdivides the
c solutions into pseudo-compounds.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer icoct,h,i,j,im,icky,id,icpct,idsol,ixct

      logical output, first, bad, chksol, wham
 
      character*10 tname, sname(h9), new*3, tn1*6, tn2*22

      double precision zt

      external chksol

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      character prject*100,tfname*100
      common/ cst228 /prject,tfname

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character*8 names
      common/ cst8 /names(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      character mname*8
      common/ cst18a /mname(m4)

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer isoct
      common/ cst79 /isoct

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer ikp
      common/ cst61 /ikp(k1)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iam
      common/ cst4 /iam

      integer ostot
      common/ junk /ostot

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      logical badend, sck, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),badend(m4,h9),sck(h9),nrf(h9)

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)
c-----------------------------------------------------------------------
c                                 initialize counters
      ixct = 0 
c                                 gloabl coordinate counter for xcoor (cxt10)
      icoct = 0  
c                                 initialize model counter
      im = 0
c                                 a flag to check if more than one solution model 
c                                 references an internal molecular EoS.
      wham = .false.
c                                 pointer to aq solution model
      idaq = 0
c                                 no request for solutions
      if (io9.eq.1.or.isoct.eq.0) then 

         isoct = 0 
c                                 identify the fluid for aqrxdo
         call aqidst

         return 

      end if 
c                                 open pseudocompund list file
      if (output.and.lopt(10)) then

         call mertxt (tfname,prject,'_pseudocompound_list.txt',0)
         open (n8,file=tfname)

      end if 
c                                 format test line
      read (n9,'(a)') new
c                                 check version compatability
      if (.not.chksol(new)) call error (3,zt,im,new)

      do while (im.lt.isoct) 
c                                 -------------------------------------
c                                 read the solution name
         call rmodel (tname,tn1,tn2,bad)

         if (bad) cycle 
c                                 ostot is zero, if eof: 
         if (ostot.eq.0.and.isoct-im.gt.0) then 
c                                 then at least one solution phase referenced
c                                 in the input is not present in the
c                                 solution phase data file, write warning:
            if (iam.lt.3.or.iam.eq.15) 
     *                               call warn (43,zt,isoct-im,'INPUT9')
            exit

         end if 
c                                 -------------------------------------
c                                 check the solution model:
         call cmodel (im,idsol,tname,first)

         if (jstot.eq.1.and.jsmod.eq.39.and.lopt(32)) then 
c                                  lagged aqueous speciaton with a pure water solvent.
         else 
c                                  normal solution. 
            if (jstot.lt.2) cycle 
c                                 -------------------------------------
c                                 reformulate the model so that it has 
c                                 no missing endmembers:
            if (jsmod.eq.9.or.jsmod.eq.10) call kill01 (isite+1)

            if (jstot.lt.ostot) call reform (tname,im,first)

            if (ostot.lt.2) cycle

         end if 
c                                 -------------------------------------
c                                 make various book keeping arrays (y2p,
c                                 jmsol, dydz, .....)
         call nmodel
c                                 check that the name has not already been found, i.e., 
c                                 that the name is duplicated in the solution model file
         if (first) then
            do i = 1, im - 1
                  if (tname.eq.sname(i)) call error (75,0d0,i,tname)
            end do
         end if
c                                 save solution name
         sname(im) = tname
c                                 abbreviation
         aname(im) = tn1
c                                 long name
         lname(im) = tn2
c                                 save found solutions in global solution 
c                                 model arrays
         call gmodel (im,tname,wham)
c                                 generate pseudocompound compositions.
c                                 subdiv returns the total
c                                 number of pseudocompounds (ipcps) and 
c                                 array y, of which element y(h,i,j) is
c                                 the site fraction of the jth species on
c                                 the ith site of the hth pseudocompound.
         if (iam.lt.3.or.iam.eq.15) then  
c                                 vertex/meemum need static pseudocompounds
            call subdiv (im,.false.)       
c                                 subdiv generates ntot compositions,
c                                 generate the compound data for each solution:
c                                 save the identities of the endmembers

c                                 global pseudo-cpd counter for sxs
            icpct = 0

            do i = 1, kstot

               id = kdsol(knsp(i,im))
               if (iend(knsp(i,im)).eq.0) ikp(id) = im

            end do

c                                 normal solutions
            do h = 1, ntot
c                                 load the composition into the
c                                 the simple prismatic arrays x,y:
               call prs2xy (h,im,.false.,bad)

               if (bad) cycle
c                                 generate the pseudocompound:
               call soload (im,icoct,icpct,tname,icky,im)

            end do

            if (icpct.gt.0) then 
c                                 write pseudocompound count
               write (*,1100) icpct, tname
c                                 write reach_increment
               if (int(reachg(im)*2d0/nopt(21)-1d0).gt.0) 
     *            write (*,1030) int(reachg(im)*2d0/nopt(21)-1d0), tname
c                                 indicate site_check_override and refine endmembers
               if (.not.sck(im)) write (*,1080) tname
               if (nrf(im)) write (*,1090) tname

               if (output.and.lopt(10)) then

                  if (jsmod.ne.6.and.jsmod.ne.8) then
                     write (n8,1060) tname, 
     *                         (names(jend(im,2+i)), i =1, lstot(im))
                  else if (jsmod.eq.6) then 
                     write (n8,1060) tname, 
     *                         (mname(iorig(knsp(i,im))), 
     *                         i = 1, lstot(im)),'* see footnote 1'
                  else if (jsmod.eq.8) then 
                     write (n8,1060) tname, 
     *                         (mname(iorig(knsp(i,im))), 
     *                         i = 1, lstot(im)),'* see footnotes 1 & 2'
                  end if

                  do i = iphct-icpct+1, iphct
                     write (n8,1070) names(i),
     *                            (xco(jco(i)+j), j = 1, lstot(im))
                  end do

                  if (jsmod.eq.6) then 
                     write (n8,1120)
                  else if (jsmod.eq.8) then
                     write (n8,1120)
                     write (n8,1130)
                  end if 
               end if 

            end if 

            jend(im,2) = icpct

         end if 
c                               read next solution
      end do 

      if (isoct.gt.0) then 

         if (iam.lt.3.or.iam.eq.15) write (*,1110) iphct - ipoint
c                               flush for paralyzer's piped output
         flush (6)
c                               scan for "killed endmembers"
         do i = 1, ipoint
c                               reset ikp
            if (ikp(i).lt.0) ikp(i) = 0
         end do 

         if (io3.eq.0.and.output.and.(iam.lt.3.or.iam.eq.15)
     *                                     .and.isoct.ne.im) then

            write (n3,1020)
            write (n3,1010) (fname(i), i = 1, isoct)
            if (im.gt.0) then 
               write (n3,1000)
               write (n3,1010) (sname(i), i = 1, im)
            else if (output) then 
               write (n3,1040) 
            end if 

         end if 

         do i = 1, im
            fname(i) = sname(i)
         end do 

      end if

      isoct = im
c                              identify the fluid for aqrxdo
      call aqidst
c                              close pseudocompound list
      if (output.and.lopt(10)) close (n8)
c                              close solution model file
      close (n9)

      first = .false.

1000  format (/,'the following solution models will be considered:',/)
1010  format (7(2x,a10))
1020  format (/,'Of the requested solution models:',/)
1030  format (9x,'a reach_increment of ',i2,' is specified for ',a)
1040  format (/,'no models will be considered.',/)
1060  format (/,'Solution: ',a,/,12x,'Endmember fractions:',
     *        /,12x,20(a,1x))
1070  format (a,2x,20(1x,f6.3,2x))
1080  format (9x,'site_check_override is on for ',a)
1090  format (9x,'refine_endmembers is on for ',a)
1100  format (i8,' pseudocompounds generated for: ',a)
1110  format (/,'Total number of pseudocompounds:',i8)
1120  format (/,'1 - Although the bulk composition of pseudocompounds'
     *        ,' for this solution is fixed,',/,' the proportions of'
     *        ,' its endmembers may vary due to respeciation.',/)
1130  format (/,'2 - Proportions output here may sum to <1 ', 
     *          'because the ordered species',/,'may have non-zero ',
     *          'initial proportions.',/)
      end 

      subroutine errk24 (resub)
c---------------------------------------------------------------------
c if an entry will exceed dimension k24 and write apporpriate 
c diagnostic on error.
c---------------------------------------------------------------------
      implicit none

      logical resub

      logical refine
      common/ cxt26 /refine
c---------------------------------------------------------------------
c                                 error diagnostic
         if (resub) then
c                                 adaptive minimization array
            call error (41,0d0,2,'SUBDIV')
         else if (refine) then 
            call error (41,0d0,1,'SUBDIV')
         else 
            call error (41,0d0,0,'SUBDIV')
         end if

      end

      subroutine subdiv (ids,resub)
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical resub

      integer last,i,j,np1,h,index,ids, k, l, n

      double precision sum

      logical refine
      common/ cxt26 /refine

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c---------------------------------------------------------------------
      if (ksmod(ids).eq.20) then
c                                subdivision with charge balance
         call cartaq (ids,resub)
c                                assign to y()?
         return

      else if (ksmod(ids).eq.9.or.ksmod(ids).eq.10) then 

         call oddprs (ids,resub)

         return

      end if 
c                                 do the first site:

      call cartes (1,ids)

      last = ndim(1,ids)
      if (last.eq.0) last = 1

      do h = 1, npairs
c                                 load the k simplicial coordinates
c                                 into as the first k of j prismatic coordinates
         j = (h-1) * mcoor(ids)
         k = (h-1) * ndim(1,ids)

         do i = 1, last
            prism(j+i) = simp(k+i)
         end do 

      end do 

      np1 = npairs
      ntot = np1 

      if (istg(ids).eq.1) return
c                                 do the second site:
      call cartes (2,ids)
 
      last = ndim(2,ids)
      if (last.eq.0) last = 1
c                                 there will be a total of
c                                 (npairs-1)*np1 compositions,
c                                 copy the first site 2 distribution
c                                 into the first np1 compositions:
      do h = 1, np1
c                                 load the k simplicial coordinates
c                                 into as the ndim(1) + k of j prismatic coordinates
         j = (h-1) * mcoor(ids) + ndim(1,ids)

         do i = 1, last
c                                 this could be an invalid compostion for
c                                 a 3 site model.
            prism(j+i) = simp(i)

         end do
      end do

      do h = 2, npairs
c                                 for each site 2 composition,
c                                 duplicate the range of site 1
c                                 compositions.
         n = (h-1) * ndim(2,ids)

         do i = 1, np1

            ntot = ntot + 1
            l = (i-1) * mcoor(ids) 
            k = (ntot-1) * mcoor(ids)

            if (k+mcoor(ids).gt.k24) call errk24 (resub)

            do j = 1, ndim(1,ids)
               prism(k+j) = prism(l+j)
            end do 
c                                 put in the new site 2 compositions:
            k = k + ndim(1,ids)

            do j = 1, ndim(2,ids)
               prism(k+j) = simp(n+j)
            end do

         end do
      end do

      np1 = ntot
c                                 do the third site:
      if (ksmod(ids).eq.9.or.ksmod(ids).eq.10) then
c                                 prism + orphan vertices (the prism is treated
c                                 as a vertex of the simplex including the orphans)
c                                 do the orphan site:

c                                 zero the orphan concentrations in the
c                                 ntot prismatic coordinates generated so far:
         do i = 1, ntot

            k = (i-1)*mcoor(ids) + nsum(ids) + 1
            l = i*mcoor(ids)

            do j = k, l
               prism(j) = 0d0
            end do

         end do
c                                 get the npair simplicial coordinates for the 
c                                 orphan site
         call cartes (3,ids)
c                                 for each simplicial coordinate (skipping the 
c                                 origin) make ntot prismatic compositions, for
c                                 a grand total of ntot*(npair-1) compositions
         do i = 2, npairs

            sum = 1d0
c                                 starting point in simp
            n = (i-1)*ndim(3,ids)

            do j = 1, ndim(3,ids)
               sum = sum - simp(n+j)
            end do 

            do h = 1, np1

               ntot = ntot + 1
c                                 the starting position of the pure prism is
               l = (h-1) * mcoor(ids) 
c                                 the starting position in the final array is
               k = (ntot -1) * mcoor(ids) 

               if (k+mcoor(ids).gt.k24) call errk24 (resub)
c                                 load the diluted prism coordinates
               do j =  1, nsum(ids)
                  prism(k+j) = prism(l+j)*sum
               end do 
c                                 and the simplicial coordinates
               do j = 1, ndim(3,ids)
                  prism(k+nsum(ids)+j) = simp(n+j)
               end do 

            end do

         end do 

         return 

      end if

      if (istg(ids).eq.2) return

      call cartes (3,ids)
c                                 the use of "index" is necessary to avoid
c                                 compiler error if parameter mst < 3
      index = 3 
      last = ndim(index,ids)
      if (last.eq.0) last = 1
c                                 copy the first site 3 distribution
c                                 into the first np1*np2 compositions:
      do h = 1, ntot

         j = (h-1) * mcoor(ids) + nsum(ids)

         do i = 1, last
            prism(j+i) = simp(i)
         end do 

      end do 
c                                 for each site 3 composition,
c                                 duplicate the range of site 1 and 
c                                 site 2 compositions.
      do h = 2, npairs

         l = (h-1) * ndim(index,ids)

         do i = 1, np1

            ntot = ntot + 1

            l = (i-1) * mcoor(ids) 
            k = ntot * mcoor(ids)

            do j = 1, ndim(1,ids)
               prism(k+j) = prism(l+j)
            end do 

            k = k + ndim(1,ids)
            l = l + ndim(1,ids)

            do j = 1, ndim(2,ids)
               prism(k+j) = prism(l+j)
            end do 

            k = k + ndim(2,ids)

            do j = 1, ndim(index,ids)
               prism(k+j) = simp(l+i)
            end do
 
         end do 
      end do 
 
      end

      subroutine satsrt 
c---------------------------------------------------------------------
c routine to sort pseudocompounds consisting entirely of saturated
c components.
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      integer j,idc
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2
 
      do j = isat, 1, -1
         idc = icp + j
         if (cp(idc,iphct).ne.0d0) then
            isct(j) = isct(j) + 1
            if (isct(j).gt.h6) call error (17,cp(1,1),h6,'SATSRT')
            if (iphct.gt.k1) call error (180,cp(1,1),k1,
     *                                  'SATSRT increase parameter k1')
            ids(j,isct(j)) = iphct
            exit
         end if
      end do 

      end

      subroutine soload (isoct,icoct,icpct,tname,icky,im)
c--------------------------------------------------------------------------
c soload - loads/requires solution properties: 

c   jend(h9,m4)  - h9 is the maximum number of solutions
c                   k12 is the maximum number of endmembers pers
c                   solution plus two.
c   jend(i,1)     - OBSOLETE! is the number of endmembers in solution i.
c   jend(i,2)     - is the number of pseudocompounds of solution i.
c   jend(i,3-3+j) - are the indices of the j endmembers in solution i.
c   sxs(k13)      - contains the mole fractions of the endmembers
c                   in the pseudocompounds.
c   ixp(i)        - a pointer that locates the first mole fraction (- 1)
c                   of the ith pseudocompound, the remaining mole fractions
c                   follow sequentially (as in jend(i,3-3+j)).
c   ikp(i)        - the index of the solution corresponding to pseudocompound i.
c   exces(j,i)    - the excess function of pseudocompound i, accounts for 
c                   excess properties and configurational entropy as a function
c                   of pressure and temperature:

c                       gexces(i) = exces(1) + exces(2)*T + exces(3)*P 
c--------------------------------------------------------------------------

      implicit none
  
      include 'perplex_parameters.h'

      character tname*10, znm(2,2)*2, pnm(3)*2 
 
      double precision zpr,hpmelt,slvmlt,gmelt,smix,esum,ctotal,omega,x

      logical zbad

      integer id,im,h,i,j,l,icpct,isoct,icky,index,icoct,icoct0,i228

      external zbad

      double precision ctot
      common/ cst3 /ctot(k1)

      logical refine
      common/ cxt26 /refine

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      character*8 names
      common/ cst8 /names(k1)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer iddeps,norder,nr 
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      double precision qmult, d0, dcoef, scoef      
      common/ cxt1r /qmult(m10,h9),d0(m11,m10,h9),dcoef(m0,m11,m10,h9),
     *               scoef(m4,h9)

      integer jndq, jdqf, iq
      double precision dqfg, dq 
      common/ cxt9 /dqfg(m3,m4,h9),dq(m4),jndq(m4,h9),jdqf(h9),iq(m4)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      double precision wgl, wkl, vlar
      common/ cxt2r /wgl(m3,m1,h9),wkl(m16,m17,m18,h9),vlar(m3,m4,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ikp
      common/ cst61 /ikp(k1)

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip

      double precision pa, p0a, zp, w, y, z, wl
      common/ cxt7 /y(m4),zp(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      double precision y2pg
      common/ cxt4  /y2pg(m15,m4,h9)

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer msite, ksp, lterm, ksub
      common/ cxt1i /msite(h9),ksp(m10,h9),lterm(m11,m10,h9),
     *               ksub(m0,m11,m10,h9)

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      logical badend, sck, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),badend(m4,h9),sck(h9),nrf(h9)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      double precision exces
      common/ cst304 /exces(m3,k1)

      save i228
      data i228/0/
c----------------------------------------------------------------------
c                                 reject special case:
c                                 ternary coh fluids above the CH4-CO join
      if (ksmod(im).eq.41.and.y(1).ge.1d0/3d0+y(2)) return
c                                 move site fractions into p0a = pa arrays indexed 
c                                 only by independent disordered endmembers, this is 
c                                 done even for models without disorder so the pa
c                                 array can be used for all solutions. 
      call y2p0 (im)

c DEBUG DEBUG
      x = 0
      do i = 1, nstot(im)
         x = x + pa(i)
      end do 
      if (dabs(x-1d0).gt.zero) then 
         write (*,*) 'soload 2'
      end if 

      iphct = iphct + 1
      icpct = icpct + 1 

      if (iphct.gt.k1) then
         if (refine) then 
            call error (41,pa(1),1,'SOLOAD')
         else 
            call error (41,pa(1),0,'SOLOAD')
         end if 
      end if

      ikp(iphct) = isoct

      icky = 0 

      do i = 1, m3
         exces(i,iphct) = 0d0
      end do
c                                encode a name
      if (istg(im).gt.1) then 
c                                make character nums for standard cases
c                                this is only to avoid run-time errors
c                                during debugging. 
         do i = 1, istg(im)
            do j = 1, 2 
               h = idint(1d2*z(j,1))
               if (h.eq.100) then 
                  znm(i,j) = '**'
               else 
                  write (znm(i,j),'(i2)') h 
               end if 
            end do
         end do 
      end if 

      do j = 1, 3

         h = idint(1d2*pa(j))

         if (h.eq.100) then 
            pnm(j) = '**'
         else 
            write (pnm(j),'(i2)') h 
         end if 

c use mname array to flag retained absent endmembers

      end do
      
      if (istg(im).eq.2.and.mstot(im).eq.4) then
c                                special case 1, bin-bin reciprocal solution
         write (names(iphct),1020) tname, znm(1,1),znm(2,1)

      else if (istg(im).eq.2.and.mstot(im).eq.6.and.ispg(im,1).eq.3) 
     *        then
c                                special case 2, tern-bin reciprocal solution
         write (names(iphct),1060) tname, znm(1,1),znm(1,2),znm(2,1)

      else if (istg(im).eq.2.and.mstot(im).eq.6.and.ispg(im,1).eq.2) 
     *        then
c                                special case 3, bin-tern reciprocal solution
         write (names(iphct),1060) tname, znm(1,1),znm(2,1),znm(2,2)

      else if (istg(im).eq.2.and.mstot(im).eq.9) then
c                                special case 4, tern-tern reciprocal solution
         write (names(iphct),1060) znm(1,1),znm(1,2),znm(2,1),znm(2,2)

      else if (mstot(im).eq.2) then 
c                                binary solutions
         if (pa(1).gt.0.9999d0) then 
            write (names(iphct),'(a3,a)') names(jend(im,3)),'_100*'
         else if (pa(1).ge.0.98d0) then 
            write (names(iphct),'(a2,a,f5.2)') 
     *             names(jend(im,3)),'_',1d2*pa(1)
         else if (pa(1).lt.1d-6) then 
            write (names(iphct),'(a3,a)') names(jend(im,3)),'_0*'
         else if (pa(1).lt.0.02d0) then 
            write (names(iphct),'(a2,a,f5.4)') 
     *             names(jend(im,3)),'_',1d2*pa(1)
         else 
            write (names(iphct),1070) names(jend(im,3)),1d2*pa(1)
         end if 

      else if (mstot(im).eq.3) then 
c                                ternary solutions
         write (names(iphct),1060) (names(jend(im,2+j)),
     *                             pnm(j), j = 1, 2)
      else if (mstot(im).eq.4) then 
c                                quaternary solutions
         write (names(iphct),1060) tname, (pnm(j), j = 1, 3)

      else
c                                all the rest:
         icky = 1

         if (iphct.lt.1000000) then 
            write (names(iphct),1080) tname, iphct
         else if (iphct.lt.10000000) then
            write (names(iphct),1100) tname, iphct
         else
            write (names(iphct),'(i8)') iphct
         end if 

      end if 
c                                get blanks out of name:
      if (mstot(im).lt.4) then 
         call unblnk (names(iphct)) 
      else 
         call reblnk (names(iphct)) 
      end if   
c                                 initialize constants:
      smix = 0d0
      esum = 0d0

      do i = 1, icomp
         cp(i,iphct) = 0d0
      end do 
c                                 load constants:
      ctotal = 0d0
      icoct0 = icoct
c                                 load xcoors if reciprocal
      if (lrecip(im)) then 

         ico(iphct) = icoct 

         do i = 1, istg(im)
            do j = 1, ndim(i,im)
               icoct = icoct + 1
               xco(icoct) = z(i,j)
            end do
         end do

      end if

      jco(iphct) = icoct

      do h = 1, lstot(im)
c                               do not count the mole
c                               fractions of absent endmembers
         id = jend(im,2+h)

         if (h.lt.lstot(im).or.order.and.depend) then 

            icoct = icoct + 1
            if (icoct.gt.k18) call error (40,y(1),k18,'SOLOAD')

            xco(icoct) = pa(h)

         end if 

         if (pa(h).ne.0d0) then
c                              composition vector
            do l = 1, icomp

               if (ksmod(im).eq.20.and.h.gt.ns) then
                  zpr = pa(h) * aqcp(l,id-aqst)
               else
                  zpr = pa(h) * cp(l,id)
               end if

               cp(l,iphct) = cp(l,iphct) + zpr
               if (l.le.icp) ctotal = ctotal + zpr

            end do 
c                              accumulate endmember configurational entropy
            esum = esum + pa(h) * scoef(h,im)

         end if

      end do

      if (order.and.depend) then 

         do i = 1, norder 

            h = lstot(im) + i

            if (i.lt.norder) then 
                icoct = icoct + 1
                if (icoct.gt.k18) call error (40,y(1),k18,'SOLOAD')
                xco(icoct) = pa(h)
            end if
c                              split these fraction into the fractions of the
c                              consituent disordered species:
            do j = 1, nr(i)

               x = depnu(j,i)*pa(h)
               id = jend(im,2+ideps(j,i,im)) 
c                              composition vector
               do l = 1, icomp
                  cp(l,iphct) = cp(l,iphct) + x * cp(l,id)
                  if (l.le.icp) ctotal = ctotal + x * cp(l,id)
               end do 

            end do 

         end do  

      end if 

      do l = 1, icomp
         if (cp(l,iphct).gt.-nopt(5).and.cp(l,iphct).lt.0d0) then 
            cp(l,iphct) = 0d0 
         else if (cp(l,iphct).lt.0d0.and.im.ne.i228) then
            i228 = im
            call warn (228,cp(l,iphct),l,tname)
         end if 
      end do 
c                                 check if the phase consists
c                                 entirely of saturated components:
      if (ctotal.eq.0d0) then
c                                 only allowed for unconstrained minimization
         if (icopt.ge.5.or.
     *       jmct.gt.0.and.jmuct.lt.jmct) then 

            call warn (55,cp(l,iphct),l,tname)

            iphct = iphct - 1
            icoct = icoct0

            return
            
         end if 

         call satsrt 
c                                 to prevent nan's in the compositional coordinates:
         ctot(iphct) = 1d0

      else

         ctot(iphct) = ctotal

      end if 
c                                 compute ideal configurational negentropy:
      if (order) then
c                                 for cpd formation models, configurational entropy
c                                 is evaluated from speciation.
         smix = 0d0

      else if (jsmod.eq.24) then 
c                                 hp melt model, use internal routine to get entropy
         smix = -hpmelt(im,pa)

      else if (jsmod.eq.25) then 
c                                 ghiorso melt model, use internal routine to get entropy
         smix = -gmelt(im)

      else if (jsmod.eq.28) then 

         smix = -slvmlt()

      else if (msite(im).ne.0) then 

         smix = -omega(im,pa)

      end if 
c                                 save it:
      exces(2,iphct) = smix 
c                                 load excess terms, if not Laar or ordered:
      if (extyp(im).eq.0.and.(.not.order)) then 

         do i = 1, jterm(im)

            zpr = 1d0

            do j = 1, jord(im)
               if (jsub(j,i,im).ne.0) zpr = zpr * pa(jsub(j,i,im))
            end do  

            do j = 1, m3
               exces(j,iphct) = exces(j,iphct) + zpr * wgl(j,i,im)
            end do 

         end do  

      else if (extyp(im).eq.1) then
c                                 redlich kister; expand polynomial
c                                 G Helffrich, 4/16
         do i = 1, jterm(im)
            do j = 1, rko(i,im)
               zpr = y(jsub(1,i,im))*y(jsub(2,i,im))
     *             * (y(jsub(1,i,im))-y(jsub(2,i,im)))**(j-1)
               do l = 1, m3
                  exces(l,iphct) = exces(l,iphct) +
     *               zpr * wkl(l,j,i,im)
               end do
            end do
         end do

      end if
c                              dqf corrections are also be saved in the
c                              exces array this implies that speciation
c                              does not effect the amount of the dqf'd
c                              endmembers.

c                              p0dord converts the p0 of any ordered species
c                              to it's disordered equivalents, as necessary
c                              for the dqf.
      if (depend.and.order) call p0dord (im)

      do i = 1, jdqf(im)
c                              index points to the endmember in the full
c                              model:
         index = jndq(i,im)
c                              user has made a dqf to an ordered species
c                              or a dependent endmember
         if (kdsol(knsp(index,im)).lt.0) 
     *                        call error (227,exces(1,1),index,tname)

         if (depend) then
            do j = 1, m3
               exces(j,iphct) = exces(j,iphct) + p0a(index)*dqfg(j,i,im)
            end do 
         else 
            do j = 1, m3
               exces(j,iphct) = exces(j,iphct) + y(index)*dqfg(j,i,im)
            end do 
         end if 

      end do

1020  format (a2,a2,'_',a2)
1060  format (a2,a2,a2,a2)
1070  format (a3,'_',f4.1)
1080  format (a2,i6)
1100  format (a1,i7)

      end

      double precision function gkomab (id,jd,vdp)
c---------------------------------------------------------------------
c evaluate g for iron according to the EoS of Komabayashi & Fei (JGR,2010)
c id points to a phase of iron
c jd points to its parameters in thermo
c vdp is the vdp integral for all phases except HCP
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id,jd

      double precision  g,vdp

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      g = thermo(4,jd) + t*(thermo(5,jd) + thermo(6,jd)*dlog(t)
     *  + t*(thermo(7,jd) + t*thermo(8,jd))) + thermo(9,jd)/t

      if (id.eq.600) then
c                                 BCC iron
         if (t.gt.1811d0) then 
            g = -25383.581d0 + t*(299.31255d0 - 46d0*dlog(t)) 
     *                       + 2.29603d31*t**(-9)
         end if 

      else if (id.eq.601) then 
c                                 FCC iron
           g = g - 2476.28 * dsqrt(t)

      else if (id.eq.602) then 
c                                 HCP iron 
           g = g - 2476.28 * dsqrt(t)
c                                 vdp from daewaele EoS
      else if (id.eq.603) then
c                                 liquid iron, destabilize at T < 1811
      end if 

      gkomab = g + vdp

      end 

      double precision function hserfe (t)
c-----------------------------------------------------------------------
c hserfe returns the hser(fe) function of Lacaze & Sundman 1990.
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1811d0) then 
         hserfe  = 1224.83d0 + (124.134d0 -23.514d0*dlog(t)
     *                  + (-.439752d-2-.5892691d-7*t)*t)*t + 77358.5d0/t
      else 
         hserfe = -25384.451d0 + (299.31255d0 - 46d0*dlog(t))*t 
     *          + 2.2960305e31/t**9
      end if  

      end 

      double precision function hsersi (t)
c-----------------------------------------------------------------------
c hserfe returns the hser(si) function of Lacaze & Sundman 1990.
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1687d0) then 
         hsersi = -8162.61d0 + ((137.227d0-22.8318d0*dlog(t)) + 
     *                     (-.191129d-2-.355178d-8*t)*t)*t+176667d0/t
      else 
         hsersi = -9457.64d0 + t*(167.272d0 - 27.196d0*dlog(t)) 
     *        - .420369e31/t**9
      end if 

      end 


      double precision function fefcc (t)
c-----------------------------------------------------------------------
c fefcc returns the gfefcc function after Andersson and Sundman, 1987
c for calculation of Fe-Cr phase diagram
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.1811d0) then 
         fefcc  = -0.23757d3 + 0.132416d3 * t - 0.246643d2 * t * dlog(t)
     *            - 0.375752d-2 * t ** 2 - 0.589269d-7 * t ** 3 + 
     *              0.773585d5 / t

      else 
         fefcc = -0.27098266d5 + 0.30025256d3 * t - 0.46d2 * t * dlog(t)
     *           + 0.278854d32 / t ** 9

      end if

      end 

      double precision function crbcc (t)
c-----------------------------------------------------------------------
c crbcc returns the gcrbcc function after Andersson and Sundman, 1987
c for calculation of Fe-Cr phase diagram
c-----------------------------------------------------------------------
      implicit none

      double precision t
c----------------------------------------------------------------------

      if (t.lt.2180d0) then 
         crbcc = -0.885193d4 + 0.15748d3 * t - 0.26908d2 * t * dlog(t) +
     *            0.189435d-2 * t ** 2 - 0.147721d-5 * t ** 3 + 
     *            0.139250d6 / t

      else 
         crbcc = -0.34864d5 + 0.34418d3 * t - 0.50d2 * t * dlog(t) - 
     *            0.288526d33 / t ** 9

      end if 


      end 

      double precision function hserc (t)
c-----------------------------------------------------------------------
c hserc returns the reference Gibbs energy of C
c-----------------------------------------------------------------------
      implicit none
      double precision t


      if (t.ge.1d-2.and.t.lt.103d0) then
         hserc = -0.104914084D4 - 0.9009204D-1 * t - 0.275D-4 * t ** 3

      else if (t.ge.103d0.and.t.le.350d0) then
         hserc = -0.98825091D3 - 0.739898691D1 * t + 0.176583D1 * t *
     #           dlog(t) - 0.1706952D-1 * t ** 2

      else
         hserc = -0.17368441D5 + 0.17073D3 * t - 0.243D2 * t * dlog(t) 
     #           - 0.4723D-3 * t ** 2 + 0.2562600D7 / t - 
     #             0.2643D9 / t ** 2 + 0.12D11 / t ** 3

      end if
      end

      double precision function glacaz (id)
c---------------------------------------------------------------------
c evaluate various CALPHAD g functions
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      double precision hserfe, hsersi, crbcc, hserc, fefcc

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c---------------------------------------------------------------------
      if (id.eq.610) then
c                                 Fe-bcc
         glacaz = hserfe(t)

      else if (id.eq.611) then
c                                 Si-bcc 
         glacaz = 0.47D5 - 0.225D2*t + hsersi(t)

      else if (id.eq.612) then
c                                 Fe-fcc
         glacaz = fefcc(t)

      else if (id.eq.613) then
c                                 Si-fcc
         glacaz = 0.51d5 - 0.218D2*t + hsersi(t) 

      else if (id.eq.614) then
c                                 Fe-liq
         if (t.lt.1811d0) then
            glacaz = 0.1204017d5 - 0.655843d1 * t - 0.36751551d-20 * 
     *               t ** 7 + hserfe(t)
         else
            glacaz = -0.108397d5 + 0.291302d3 * t - 0.46d2 * t * 
     *                dlog(t)
         end if

      else if (id.eq.615) then
c                                 Si-liq
         if (t.lt.1687d0) then
            glacaz = 0.506964D5 - 0.300994D2*t + 0.209307d-20*t**7 
     *             + hsersi(t) 
         else
            glacaz = 0.49828D5 - 0.295591D2*t + 0.420369d31/t**9 + 
     *               hsersi(t)

         end if 

      else if (id.eq.616) then
c                                 Fe2Si 
         glacaz = -0.237522D5 - 0.354D1*t + 0.67D0*hserfe(t) 
     *          + 0.33D0 * hsersi(t)

      else if (id.eq.617) then
c                                 Fe5si3
         glacaz = -0.30143D5 + 0.27D0*t + 0.625D0*hserfe(t) 
     *          + 0.375D0 * hsersi(t)

      else if (id.eq.618) then
c                                 FeSi
         glacaz = -0.363806D5 + 0.222D1*t + hserfe(t)/2D0 
     *           + hsersi(t)/2D0

      else if (id.eq.619) then
c                                 FeSi2 
         glacaz = -0.27383D5 + 0.348D1*t + 0.33D0*hserfe(t) 
     *           + 0.67D0*hsersi(t)

      else if (id.eq.620) then
c                                 Fe3Si7
         glacaz = -0.19649D5 - 0.92D0*t + 0.3D0 * hserfe(t) 
     *           + 0.7D0 * hsersi(t)

      else if (id.eq.621) then
c                                 Si-diamond
         glacaz = hsersi(t)

      else if (id.eq.622) then
c                                 FeC-BCC
         glacaz = hserfe(t) + 0.269943d6 + 0.587857d3 * t
     *            - 0.729d2 * t * dlog(t) - 0.14169d-2 * t ** 2 
     *            + 0.76878d7 / t - 0.7929d9 /t ** 2 + 
     *            0.360d11 / t ** 3

      else if (id.eq.623) then
c                                 SiC-BCC
         glacaz = 0.47D5 - 0.225D2*t + hsersi(t) + 0.269944677d6 + 
     *            0.436523d3 * t - 0.729d2 * t * dlog(t) - 0.14169d-2 *
     *            t ** 2 + 0.76878d7 / t - 0.7929d9 / t ** 2 + 
     *            0.36d11 / t ** 3 

      else if (id.eq.624) then
c                                 FeC-FCC
         if (t.lt.1811d0) then
               
            glacaz = 0.58376159d5 + 0.163135d3 * t - 0.2545D2 * 
     *               t * dlog(t) + 0.1677d-3 * t ** 2 + 0.256260d7
     *               / t - 0.2643d9 / t ** 2 + 0.12D11 / t ** 3 + 
     *               hserfe(t)

         else

            glacaz = 0.32740293d5 + 0.45510556d3 * t - 0.703d2 * t
     *               * dlog(t) - 0.4723d-3 * t ** 2 + 
     *               0.25626d7 / t - 0.2643D9 / t ** 2 + 0.12D11 /
     *               t ** 3 + 0.278854D32 / t ** 9

         end if

      else if (id.eq.625) then
c                                 SiC-FCC
         glacaz = hsersi(t) - 0.37879d5 + 0.20943d3 * t - 0.243d2 * t
     *          * dlog(t) - 0.4723d-3 * t ** 2 + 0.25626d7 / t - 
     *          0.2643d9 / t ** 2 + 0.120d11 / t ** 3 

         else if (id.eq.626) then
c                                 C-LIQ
         glacaz = 117369d0 - 24.63*t + hserc(t)

      else if (id.eq.627) then
c                                 C-GPH
         glacaz = -0.17368441d5 + 0.17037d3 * t - 0.243d2 * t * dlog(t)
     *            - 0.4723d-3 * t ** 2 + 0.25626d7 / t - 0.2643d9 / 
     *            t ** 2 + 0.12d11 / t ** 3

      else if (id.eq.628) then
c                                 SiC
         if (t.lt.700d0) then

            glacaz = -0.85572264D5 + 0.1732005D3 * t - 0.25856D2 * t * 
     *                dlog(t) - 0.2107D-1 * t ** 2 + 0.32153D-5 * t ** 3
     *                + 0.438415D6 / t

         else if (t.gt.700d0.and.t.lt.2100) then

            glacaz = -0.95145902d5 + 0.300346d3 * t - 0.45093d2 * t * 
     *                dlog(t) - 0.367d-2 * t ** 2 + 0.22d-6 * t ** 3 + 
     *                0.1341065d7 / t

         else

            glacaz = -0.105007971d6 + 0.360309d3 * t - 0.53073d2 * t * 
     *                dlog(t) - 0.74525d-3 * t ** 2 + 0.173167d-7 * 
     *                t ** 3 + 0.3693345d7 / t

         end if

      else if (id.eq.629) then
c                                 Cementite
         glacaz = -0.10745d5 + 0.70604d3 * t - 0.1206d3 * t * dlog(t)


      else if (id.eq.630) then
c                                 Fe8Si2C
         glacaz = -0.210043d5 + 0.506d0 * t + 0.91d-1 * (-0.17368441d5
     *            + 0.17037d3 * t - 0.243d2 * t * dlog(t) - 0.4723d-3 
     *            * t ** 2 + 0.25626d7 / t - 0.2643d9 / t ** 2 + 
     *            0.12d11 / t ** 3) + 0.727d0 * hserfe(t) + 
     *            0.182d0 * hsersi(t)

      else if (id.eq.631) then
c                                 C diam
         glacaz = -0.16359441D5 + 0.17561D3 * t - 0.2431D2 * t * 
     *             dlog(t) - 0.4723D-3 * t ** 2 + 0.2698D7 / t - 
     *             0.261D9 / t ** 2 + 0.111D11 / t ** 3 
      
      else if (id.eq.632) then
c                                 Cr-BCC
           glacaz = crbcc(t)

      else if (id.eq.633) then
c                                 Cr-FCC
           glacaz = crbcc(t) + 7284d0 + 0.163d0*t

      else if (id.eq.634) then
c                                 Cr_LIQ
         if (t.lt.2180d0) then

            glacaz = crbcc(t) + 0.2433593D5 - 0.1142D2 * t + 
     *               0.237615D-20 * t ** 7

         else

             glacaz = -0.16459d5 + 0.335618d3 * t - 
     *                0.50d2 * t * dlog(t)

         end if

      else if (id.eq.635) then
c                                 FeCr-sig after Hertzman, 1980; used by Nastia
c          glacaz = (183802.0638d0 + 26d0*hserfe(t) +
c    *               4d0*crbcc(t))/30d0
c                                 FeCr-sig after Andersson & Sundman 1987; used by George
         glacaz = (8d0*fefcc(t) + 4d0*crbcc(t) + 18d0*hserfe(t) +
     *            117 300d0 - 95.96d0*t) / 30d0

      else if (id.eq.636) then
c                                 CrFe-sig after Hertzman, 1980; used by Nastia
c        glacaz =(-22624.05686d0 + 20d0*crbcc(t)
c    *            + 10d0*hserfe(t))/30d0
c                                 CrFe-sig after Andersson & Sundman 1987; used by George
         glacaz = (8d0*fefcc(t) + 22d0*crbcc(t)
     *            + 92 300d0 - 95.96d0*t) / 30d0

      else if (id.eq.637) then
c                                 Fe7C3 after Djurkovic et al., 2011
         glacaz =-0.2345062954D5 + 0.1761006488D4 * t - 
     *            0.2975999679D3 * t * dlog(t) - 0.3148668241D-3 * 
     *            t ** 2 + 0.1708400854D7 / t - 0.1762000881D9 / 
     *            t ** 2 + 0.8000004D10 / t ** 3

      end if

      end 

      double precision function gmag (x)
c-----------------------------------------------------------------------
c gmag returns the magnetic contribution to G for BCC Fe in FeSi alloy.
c after Lacaze & Sundman 1990.
c     x - bulk mole fraction of Fe 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision b,t0,tc,f,x 

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      if (x.eq.0d0) then 
         gmag = 0d0
         return
      end if 

      tc = ((-0.1008D4 * x + 0.1512D4) * x + 0.539D3) * x

      t0 = t/tc

      if (t0.lt.1d0) then 

         f = 1d0 - 0.905299383D0 / t0 - (0.153008346D0 
     *       + (0.680037095D-2 + 0.153008346D-2 * t0**6) * t0**6)
     *       * t0**3

      else 

         f = -(0.641731208D-1 + (0.203724193D-2 + 0.42782080051D-3 / 
     *        t0**10) / t0**10) / t0**5

      end if 

      b = 2.22d0 * x

      gmag = r*t*f*dlog(b+1)

      end 

      double precision function gmags (tc,b,pee)
c-----------------------------------------------------------------------
c gmags returns the magnetic contribution to G parameterized as in
c Sundman 1991 J. Ph. Equil. v12, 127-140.
c     tc - transition temperature
c     b - Bohr magneton value
c     pee - structural magnetic parameter
c
c Not described or referenced in the original ref, but in an obscure reference,
c viz. Hertzman and Sundman (1982), CALPHAD v6. 67-80, negative tc means the
c material is antiferromagnetic, and the magnetic strength will also be a
c negative number by convention.  As explained in Xiong et al. (2012)
c CALPHAD v39, 11-20, to get actual values for evaluation, one examines the
c structural parameter pee.  For bcc structure p = 0.4, while for fcc & hcp
c structures p = 0.28.  The convention is for bcc (p=0.4) antiferromagnetic
c structures, divide tc and b by -1; for fcc & hcp (p=0.28) antiferromagnetic
c structures, divide tc and b by -3.

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision a0,a1,t5,t15,t25,f1,f0,f3,f9,f15

      parameter (a0=518d0/1125d0, a1=11692d0/15975d0)
      parameter (t5=1d0/10d0, t15=1d0/315d0, t25=1d0/1500d0)
      parameter (f1=79d0/140d0, f0=474d0/497d0, f3=1d0/6d0,
     *           f9=1d0/135d0, f15=1d0/600d0)

      double precision a,b,t0,tc,f,pee,bc

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      if (tc.lt.0d0) then
         if (pee.lt.0.4d0) then
            bc = -b/3
            t0 = -3*t/tc
         else
            bc = -b
            t0 = -t/tc
         endif
      else
         t0 = t/tc
         bc = b
      endif

      a = a0 + a1*(1d0/pee - 1d0)

      if (t0.lt.1d0) then 

         f = t - (f1*tc/pee + t * f0 * (1d0/pee - 1d0) *
     *       (f3 + t0**6 * (f9 + t0**6 * f15)) * t0**3) / a

      else 

         f = -t*(t5 + (t15 + t25 / t0**10) / t0**10) / t0**5 / a

      end if 

      gmags = r*f*dlog(bc+1)

      end 

      subroutine mtrans (gval,vdp,ndu,id)
c----------------------------------------------------------------------
c mtrans sorts through and evaluates transition functions
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      double precision gval, dg, vdp, ndu, tc, b, pee, gmags

      external gmags

      integer eos
      common/ cst303 /eos(k10)

      double precision therdi,therlm
      common/ cst203 /therdi(m8,m9),therlm(m7,m6,k9)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c---------------------------------------------------------------------- 
         if (ltyp(id).eq.1) then
c                                 ubc-type transitions
            call lamubc (p,t,dg,lmda(id),lct(id))
            gval = gval + dg
 
         else if (ltyp(id).eq.2) then
c                                 standard transitions
            call lamhel (p,t,gval,vdp,lmda(id),lct(id))

         else if (ltyp(id).eq.3) then
c                                 supcrt q/coe lambda transition
            call lamqtz (p,t,gval,ndu,lmda(id),id)
 
         else if (ltyp(id).eq.4) then

            if (eos(id).ne.8.and.eos(id).ne.9) then 
c                                 putnis landau model as implemented in hp98 
               call lamla0 (dg,vdp,lmda(id))
               
            else
             
               call lamla1 (dg,vdp,lmda(id))
               
            end if 

            gval = gval + dg 

         else if (ltyp(id).eq.5) then
c                                 holland and powell bragg-williams model
            call lambw (dg,lmda(id))
            gval = gval + dg

         else if (ltyp(id).eq.7) then
c                                 George's Hillert & Jarl magnetic transition model
            if (lct(id).gt.1) write(0,*)'**>1 type = 7 trans.!?'
            tc = therlm(1,1,lmda(id))
            b = therlm(2,1,lmda(id))
            pee = therlm(3,1,lmda(id))
            gval = gval + gmags (tc,b,pee)

         else 

            write (*,*) 'no such transition model'
            call errpau
 
         end if

      end 

      double precision function gfesi (y,g1,g2)
c-----------------------------------------------------------------------
c gfesi returns the Gibbs free energy for BCC FeSi alloy after 
c Lacaze & Sundman 1990. See FeSiBCC.mws.

c    y   - the bulk Fe mole fraction
c    g01 - free energy of Bcc Fe, without Gmag
c    g02 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical done

      integer itic 

      double precision g1, g2, y, x, w0, w1, w2, rt, dg, xmin, 
     *                 d2g, gord, xmax, dx, gfesi0, g0, g12, gmag

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save w1, w2, gord
      data w1, w2, gord/-11544d0, 3890d0, -10475.64d0/
c----------------------------------------------------------------------

c      g1 = g1p - gmag(1d0)

      if (y.le.nopt(5).or.y.ge.1d0-nopt(5)) then 
c                                 endmember compositions, no order possible
         gfesi = y*g1 + (1d0-y)*g2 + gmag(y)
         return
      end if 

c!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!

      w0  = -27809d0 + 11.62d0 * t
      gord = -10475.64d0*2 + ( (g1 + g2)/2d0 + w0)
      rt  = r*t
      g12 = 2d0*(gord - w0) - g1 - g2
c                                 max concentration of ordered species
      if (y.gt.0.5d0) then 
         xmax = 1d0
c                                 the true xmin (commented) allows for
c                                 anti-ordering, but because the model 
c                                 is symmetric, i up xmin to y
c        xmin = 2d0*(y-.5d0)
      else
         xmax = 2d0*y
c        xmin = 0d0
      end if 

      xmax = xmax - nopt(5)
      xmin = y + nopt(5)
      x = xmax      
c                                 get 1st and 2nd derivatives
      call dgfesi (dg,d2g,y,x,g12,rt)

      done = .false.
c                                 find starting point for newton-raphson
c                                 search
      if (dg.gt.0d0.and.d2g.gt.0d0) then 
c                                 the max order concentration is a
c                                 good starting point
         dx = -dg/d2g

      else if (dg.lt.0d0) then
c                                 the max order is a minimum
         x = y
         done = .true.

      else 
c                                 try the max disordered concentration
         x = xmin

         call dgfesi (dg,d2g,y,x,g12,rt)         

         if (dg.lt.0d0.and.d2g.gt.0d0) then 
c                                 ok
            dx = -dg/d2g
            
         else                
c                                 full disordered
            done = .true.            

         end if

      end if 
c                                 iteration loop
      if (.not.done) then                   
c                                 increment and check bounds 
         call pcheck (x,xmin,xmax,dx,done)
c                                 iteration counter 
         itic = 0

         do

            call dgfesi (dg,d2g,y,x,g12,rt)

            dx = -dg/d2g 

            call pcheck (x,xmin,xmax,dx,done)

            if (done) then 

               exit

            else 

               itic = itic + 1
               if (itic.gt.iopt(21)) exit

            end if

         end do 

      end if 

      gfesi = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)

      if (iopt(17).ne.0) then 
c                                 order check, compare to the 
c                                 max order g
         g0 = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
         if (gfesi.gt.g0) gfesi = g0 
c                                 min order g
         g0 = gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
         if (gfesi.gt.g0) gfesi = g0 

      end if 
c                                 add magnetic component
      gfesi = gfesi + gmag(y) 

      end 

      double precision function gfesi0 (y,x,gord,g2,g12,w0,w1,w2,rt)
c-----------------------------------------------------------------------
c gfesi0 computes the G for BCC FeSi alloy once the speciation has been
c computed if function gfesi. See FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c    g1 - free energy of Bcc Fe, without Gmag
c    g2 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none 

      double precision g2, g12, y, x, w0, w1, w2, xy, yx, x1, rt, 
     *                 gord
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      gfesi0  = ( dlog(x/x1*xy/yx)*x/2d0 
     *          + dlog(yx/xy)*y 
     *          + dlog(xy*x1)/2d0)*rt 
     *      - g12*yx*x 
     *      - 64d0*w2*y**4
     *      + 16d0*(8d0*w2 - w1)*y**3 
     *      + 4d0 *(6d0*w1 - 20d0*w2 - w0)*y**2 
     *      + 2d0 *(8d0*w2 + gord + w0 - 4d0*w1 - g2)*y + g2

      end 

      double precision function gfesi1 (y,x,w0,w1,w2,rt)
c-----------------------------------------------------------------------
c gfesi0 computes the G - Gmech for BCC FeSi alloy once the speciation has been
c computed if function gfesi. See FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c    g1 - free energy of Bcc Fe, without Gmag
c    g2 - free energy of Bcc Si
c-----------------------------------------------------------------------
      implicit none 

      double precision y, x, w0, w1, w2, xy, yx, x1, rt, gcon, gex
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      gfesi1  = ( dlog(x/x1*xy/yx)*x/2d0 
     *          + dlog(yx/xy)*y 
     *          + dlog(xy*x1)/2d0)*rt 
     *          + (((-64d0*w2*y + 128d0*w2 - 16d0*w1)*y 
     *              + 24d0*w1 - 80d0*w2 - 4d0*w0)*Y + 4d0*x*w0 + 2d0*w0 
     *              + 16d0*w2 - 8d0*w1)*y - 2d0*x**2*w0


      gcon = ( dlog(x/x1*xy/yx)*x/2d0 
     *          + dlog(yx/xy)*y 
     *          + dlog(xy*x1)/2d0)*rt 


      gex  = 
     *          + (((-64d0*w2*y + 128d0*w2 - 16d0*w1)*y 
     *              + 24d0*w1 - 80d0*w2 - 4d0*w0)*Y + 4d0*x*w0 + 2d0*w0 
     *              + 16d0*w2 - 8d0*w1)*y - 2d0*x**2*w0


      end 

      subroutine dgfesi (dg,d2g,y,x,g12,rt)
c-----------------------------------------------------------------------
c dgfesi first and second derivatives of gfesi with respect to the ordered
c species concentration (x). After Lacaze & Sundman 1990, see FeSiBCC.mws.

c    y  - the bulk Fe mole fraction
c-----------------------------------------------------------------------
      implicit none

      double precision y, x, xy, yx, x1, rt, dg, d2g, g12
c-----------------------------------------------------------------------
      yx  = 2d0*y - x
      x1  = 1d0 - x
      xy  = 1d0 - 2d0*y + x

      dg  = -2d0*(y - x)*g12 + dlog(x*xy/x1/yx)*rt/2d0

      d2g =  2d0*g12 + (xy/x1/yx + x/x1/yx + x*xy/x1**2/yx 
     *               + x*xy/x1/yx**2)/x/xy*x1*yx*rt/2d0

      end 

      double precision function gfesic (y1,y3,y4,g1,g2,g3,g4,id)
c---------------------------------------------------------------------
c fesic4 returns the free energy change for Fe-Si-C alloy after 
c Lacaze & Sundman 1990. 

c    id = 30 -> BCC
c    id = 31 -> FCC

c    y1..y4 - mole fractions of Fe, Si, FeC and SiC, respectively
c    g1..g4 - free energies of Fe, Si, FeC and SiC, respectively
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision g1, g2, g3, g4, y1, y3, y4, gmech,
     *                 logu, logx, gconf, gex, x, y, u, v, gmag 

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c     x is the site fraction of Fe on the first site, y is the site 
c     fraction of Si on the first site, u is the site fraction of C 
c     on the second site, v is the site fraction of vacancies on the
c     second site

      x = y1 + y3
      u = y3 + y4
      y = 1d0 - x
      v = 1d0 - u

      gmech = x * v * g1 + y * v * g2 + x * u * g3 + y * u * g4

      if (x.gt.0d0.and.x.lt.1d0) then 

         logx = (x * dlog(x) + y * dlog(y))

      else

         logx = 0d0

      end if 

      if (u.gt.0.and.u.lt.1d0) then 

         logu = (u * dlog(u) + v * dlog(v))

      else

         logu = 0d0

      end if 

      if (id.eq.30) then 
c                                 BCC
         gconf = r*t*(logx + 3d0*logu)

         gex = x*y*v*(-0.153138560d6 + 0.4648d2 * t - 0.92352d5 * x +
     *         0.92352d5*y + 0.62240d5*(x-y)**2) + 0.78866d5 *x*y*u
     *         - 0.190d3*x*u*v*t + gmag(1d0)

      else if (id.eq.31) then 
c                                 FCC
         gconf = r * t * (logx + logu)

         gex = x*y*v* (-0.1252477d6 + 0.41116d2 * t - 0.1427076d6 * x
     *         + 0.1427076D6 * y + 0.899073D5 * (x - y) ** 2) + 
     *         x * y * u * (0.1432199d6 + 0.3931d2 * t - 0.2163205d6 * x
     *         + 0.2163205D6 * y) - 0.34671d5 * x * u * v

      end if 

      gfesic = gmech + gconf + gex

      end

      double precision function gfecr1 (y,g1,g2)
c-----------------------------------------------------------------------
c     gfecr1 returns the free energy change for BCC Fe-Cr alloy after 
c     Andersson and Sundman, 1987, updated by Xiong et al. 2011
c     (solution model id = 32) 
c     The only reason this model needs to be built in is due to the
c     continuous change in magnetic Tc and B through the FM-AFM transition.
c     Otherwise it could be a normal solution model

c    y      - mole fractions of Fe
c    g1, g2 - free energies of Fe-bcc and Cr-bcc

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision g1, g2, y, gmech, 
     *                 gconf, gex, gmag2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      gmech = y * g1 + (1-y) * g2 

      if (y.lt.1d0.and.y.gt.0d0) then

          gconf = r * t * (y*dlog(y) + (1-y)*dlog(1-y)) 
      else
          gconf = 0d0
      
      end if

c     gex below is from Xiong et al. (2011) and uses an alternative magnetic
c     model than is implemented in gmag2/gmags.
      gex = (1 - y) * y * (0.2421206D5 - 0.15507D2 * t + 
     *      (1 - 2 * y) * (0.166469D4 + 0.286D0 * t) + 
     *      ((1 - 2 * y) ** 2) * (-0.1325088D5 + 0.8252D1 * t))
c     gex below is from Andersson & Sundman (1987) whose magnetic model is
c     consistent with gmag2/gmags.  This doesn't work due to unknown errors
c     in the Andersson & Sundman (1987) data listed in the article.
c     gex = y*(1d0-y)*(20 500d0 - 9.68d0*t)

      gfecr1 = gmech + gconf + gex + gmag2(y)

      end

      double precision function gmag2 (x)
c-----------------------------------------------------------------------
c gmag2 returns the magnetic contribution to G for BCC Fe in FeCr alloy
c after Andersson and Sundman (1987).
c     x - bulk mole fraction of Fe

c                                      G. Helffrich, ELSI, 8 Apr. 2016.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision b,tc,x,gmags,xfe,xcr

      if (x.eq.0d0) then 

         gmag2 = 0d0

      else

c        tc = 0.13545D4 * x - 0.3115D3 + (1d0 - x) * x * 
c    *        (0.2200D4 - 0.11D4 * x)
c        b = 0.2228D1 * x - 0.8D-2 - 0.85D0 * (1d0 - x) * x

         xfe = x
         xcr = 1d0-x
         tc = 1043d0 * xfe + (-311.5d0) * xcr +
     *        xfe * xcr * (1650d0 + 550d0*(xcr-xfe))
         b = 2.22d0 * xfe + (-0.008d0) * xcr + xfe*xcr*(-0.008d0)
         gmag2 = gmags (tc,b,0.40d0)

      end if

      end 

      double precision function gmet (id)
c----------------------------------------------------------------
c function reads SGTE data format for reference Gibbs free energy
c and evaluates thermal and pressure EoS
c
c polynomial for Gibbs free energy function
c g = a + b*T + c*T*lnT + d/T + e/T**2 + f/T**3 + g/T**9 + 
c         h*T**2 + i*T**3 + j*T**4 + k*T**7
c
c EOS after Brosh et al., 2007, 2008:
c v0 volume at pr,tr; nn number of atoms; gam0 Grüneisen parameter;
c tet0 Enstein temperature; b1,dd1,b0,dd0 fitting coefficients;
c Bo bulk modulus; Bpo pressure derivative of bulk modulus
c                                 -------------------------------
c in terms of coefficient tags nastia's polynomial is 

c g = c1 + c2*T + c3*T*lnT + c4/T + c5/T**2 + c6/T**3 + c7/T**9 + 
c         c8*T**2 + c9*T**3 + c10*T**4 + c11*T**7

c to accomodate saxena & eriksson (2015, Fe-S) with minimal effort
c for eleanor i use the G0 and S0 tags for the sqrt(T) and ln(T)
c SGTE terms. creating a set of tags specifically for SGTE would 
c be preferable, but require rewriting nastia's data (which might, 
c or might not, be worthless), so the SGTE polynomial is now
c
c g = c1 + c2*T + c3*T*lnT + c4/T + c5/T**2 + c6/T**3 + c7/T**9  
c        + c8*T**2 + c9*T**3 + c10*T**4 + c11*T**7 + G0*sqrt(T)
c        + S0*ln(T)
c
c G0 and S0 are loaded into thermo(31...) and thermo(32...) via
c the ic2p pointer array.

c i didn't follow the cp/s/h stuff through the contorted 'cst1'
c 'cst2' mess to see if it's ever relevant, if it is then the
c code must be modified to include the appropriate G0/S0 terms.
c if that stuff really is necessary it's probably better to clean
c this up rather than to apply more plasters. 

c                                           JADC, 12/3/2017
c-----------------------------------------------------------------
      implicit none
      include 'perplex_parameters.h'

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer id, nn

      double precision a,b,c,d,e,f,g,h,i,j,k
      double precision v0,gam0,tet0,tet02,gam02,b1,dd1,b0,dd0,Bo,Bpo
      double precision gsgte,gsgte0,sr,hr,cpr
      double precision gqh,gqhr,sqhr,hqhr,cpqhr
      double precision intp,difc,gbrosh
      double precision colcom,harter
      double precision tc,tc1,t0,beta,pp,ff,gmagn,gmagn0,vv,cst1,cst2
c----------------------------------------------------------------------
c                             allocate polynomial coeffiecients for
c                             reference Gibbs free energy function
      a = thermo(1,id)
      b = thermo(2,id)
      c = thermo(3,id)
      d = thermo(4,id)
      e = thermo(5,id)
      f = thermo(6,id)
      g = thermo(7,id)
      h = thermo(8,id)
      i = thermo(9,id)
      j = thermo(10,id)
      k = thermo(11,id)
c                             allocate coefficients for EoS
c                             and magnetic term
      gam0 = thermo(12,id)
      nn = thermo(13,id)
      tet0 = thermo(14,id)
      b0 = thermo(15,id)
      dd0 = thermo(16,id)
      b1 = thermo(17,id)
      dd1 = thermo(18,id)
      Bo = thermo(19,id)
      Bpo = thermo(20,id)
      v0 = thermo(22,id)
      tc = thermo(23,id)
      beta = thermo(24,id)
      pp = thermo(25,id)
      vv = thermo(26,id)
c                           thermodynamic properties calculated at tr=0 K
c                           cst1: hsgte(at reference T) - hqh(at ref T) -
c                                 - tr/2*(Cpsgte(at ref T) - Cpqh(at ref T)) = hsgte(at reference T)
c                           cst2: (sqh(at ref T) - ssgte(at ref T)) + 
c                                 + (Cpsgte(at ref T)-Cpqh(at ref T)) = - ssgte(at ref T)
      cst1 = thermo(27,id)
      cst2 = thermo(28,id)
c                           additional Grueneisen parameter and Einstein temperature for graphite
      gam02 = thermo(29,id)
      tet02 = thermo(30,id)

c                          read SGTE data 
      gsgte = a + b*t + c*t*dlog(t) + d/t + e/t**2 + f/t**3 
     *                + g/t**9 + h*t**2 + i*t**3 + j*t**4 + k*t**7
     *                + thermo(31,id)*dsqrt(t) + thermo(32,id)*dlog(t)
c                          check for transitions:
      if (ltyp(id).ne.0) call calpht (t,gsgte,lmda(id),lct(id))
        
c                          quasi-harmonic term
      if (nn.eq.0) then
c                          quasi-harmonic term for graphite
          gqh = 1d0*r*t*dlog(1d0 - dexp(-tet0/t)) + 
     *          2d0*r*t*dlog(1d0 - dexp(-tet02/t))

      else

          gqh = 3d0*nn*r*t*dlog(1d0 - dexp(-tet0/t))
      
      end if

c                             interpolation function
      intp = 1d0/(1D0 + b1)*(b1 + dsqrt(1D0 + 
     *              2D0*b1*(1D0 + dd1)*p/Bo))*dexp(0.10D1/b1 - 
     *              1d0/b1*dsqrt(1D0 + 2D0*b1*(1D0 + dd1)*p/Bo))

c                    evaluate S, Cp, and H at reference T and P
c                    at tr and pr defined in the thermodynamic data file (tr=298.15 K)
      if (cst1.eq.0d0.or.cst2.eq.0d0) then

          gsgte0 = a + b*tr + c*tr*dlog(tr) + d/tr + e/tr**2 + f/tr**3 +
     *                 g/tr**9 + h*tr**2 + i*tr**3 + j*tr**4 + k*tr**7
c                                            sr = -dg/dt at reference t
          sr = -b - c*dlog(tr) - c + d/tr**2 + 2d0*e/tr**3 + 
     *              3d0*f/tr**4 + 0.9D1*g/tr**10 - 2d0*h*tr - 
     *              3d0*i*tr**2 - 0.4D1*j*tr**3 - 0.7D1*k*tr**6
c                                            hr = gr+tr*sr at reference t
          hr = gsgte0 + tr*sr
c                                            cpr = t*ds/dt at reference t
          cpr = -c - 2d0*d/tr**2 - 6d0*e/tr**3 - 12d0*f/tr**4 
     *          - 9d1*g/tr**10  
     *          - 2d0*h*tr - 6d0*i*tr**2 - 12d0*j*tr**3 - 42d0*k*tr**6

          gqhr = 3d0*nn*r*tr*dlog(1d0 - dexp(-tet0/tr))

          sqhr = 3d0*nn*r*tet0/tr/(dexp(tet0/tr) - 1d0) - 
     *                     3d0*nn*r*dlog(1d0 - dexp(-tet0/tr))

          hqhr = 3d0*nn*r*tet0/(dexp(tet0/tr) - 1d0)

          cpqhr = 3d0*nn*r*tet0**2/tr**2*dexp(-tet0/tr)/
     *                    (1d0 - dexp(-tet0/tr))**2

c                         Cp(SGTE) - Cp(QH) with low temperature correction
c                         for the case tr=298.15
          if (t.lt.tr) then

             difc = t**2/(2D0*tr)*(cpr - cpqhr)
        
          else
          
             difc = -(gsgte - hr + t*sr) + (gqh - hqhr + t*sqhr) 
     *                + (t-tr/2D0)*(cpr-cpqhr)

          end if


      else

         difc = -gsgte + gqh + cst1 + t*cst2 

      end if


      gbrosh = colcom(Bo,v0,Bpo,p) + 
     *         harter(nn,r,t,p,tet0,tet02,Bo,b0,dd0,gam0,gam02) - 
     *         gqh + difc*(1D0 - intp)


c                     magnetic contribution using Inden-Hillert-Jarl model

      if (tc.eq.0D0.or.pp.eq.0D0) then
c                     no magnetic contribution
         gmagn = 0D0

      else 
c                     pressure dependence of Tc (as modelled for cementite)
         if (vv.eq.0D0) then 
            tc1 = tc
         else
            tc1 = tc*dexp(vv*p)
         end if

          t0 = t/tc1

          if (pp.eq.0.28D0) then
c                           fcc,hcp metals and cementite
c                           for fcc-Fe beta is a negative number,
c                           therefore, gmagn can't be calculated
c                           dlog(<0)

         if (t0.lt.1D0) then
                  ff = 1d0 - 0.8603387544D0/t0 - 
     #                0.1744912404D0*t0**3 - 
     #                0.7755166236D-2*t0**9 - 0.1744912404D-2*t0**15
              else
                  ff = -0.4269022681D-1/t0**5 - 
     #                  0.1355245296D-2/t0**15 - 
     #                  0.2846015121D-3/t0**25
              end if

          gmagn0 = r*t*dlog(beta+1D0)*ff

          else if (pp.eq.0.4D0) then
c                           bcc metals
             if (t0.lt.1D0) then
                 ff = 1d0 - 0.9052993829D0/t0 - 
     #                0.1530083464D0*t0**3 - 
     *                0.6800370949D-2*t0**9 - 
     *                0.1530083464D-2*t0**15
             else
                 ff = -0.6417312080D-1/t0**5 - 
     *                 0.2037241930D-2 /t0**15 - 
     *                 0.4278208053D-3/t0**25
             end if
          gmagn0 = r*t*dlog(beta+1D0)*ff
          end if

          gmagn=gmagn0

      end if

      gmet = gsgte + gbrosh + gmagn

      end function gmet


      double precision function colcom (Bo,v0,Bpo,p)
c------------------------------------------------------------------------
c integral for cold compression path used in EOS Brosh et al., 2007, 2008
c------------------------------------------------------------------------
      implicit none

      double precision          Bo,V0,Bpo,p
      double precision          a,I4XC,G42,G41,G4L,G4M1,G4PC,n

         n = 4d0
         a = (n-1D0)/(3.*Bpo-1d0)
         I4XC = 1d0 - a + a*(1d0 + n/a*p/Bo/3d0)**(1D0/n)

         G42 = 0.15D1*Bpo**3 - 0.6D1*Bpo**2 + 
     *                         0.8D1*Bpo - 0.3555555555D1

         G41 = -9D0*Bpo**3 + 27D0*Bpo**2 - 
     *                       24D0*Bpo + 0.5333333333D1

         G4L = 9D0*Bpo**3 - 18D0*Bpo**2 + 
     *                       9D0*Bpo - 0.1333333333D1

         G4M1 = 3D0*Bpo**3 - 3D0*Bpo**2 + 
     *                           Bpo - 0.111111111D0

         G4PC = G42*I4XC**(-2) + G41*I4XC**(-1) - G4L*dlog(I4XC) + 
     *          G4M1*I4XC - G42 - G41 - G4M1

        colcom = Bo*v0*G4PC

      end function colcom


      double precision function harter (nn,r,t,p,tet0,tet02,Bo,b0,dd0,
     *                                  gam0,gam02)
c---------------------------------------------------------------------
c integrated "thermal" volume used in EOS Brosh et al., 2007, 2008
c---------------------------------------------------------------------
      implicit none

      double precision          r,t,p,tet0,tet02,Bo,b0,dd0,gam0,gam02
      double precision          b,tet,tet1,tet2,I2XT,G2XT
      integer                   m, nn

         m = 2
         b = (m-1D0)/(3d0*b0-1d0)

         I2XT = 1d0 - b + b*(1d0 + m/b*(1D0 + dd0)*
     *          P/Bo/3d0)**(1D0/m)

         G2XT = (4.5D0*b0 - 3D0)*I2XT**(-2) + (-9D0*b0+3D0)*I2XT**(-1)
     *                                                      + 4.5D0*b0
         
         if (nn.eq.0) then

             tet1 = tet0*dexp(gam0/(1D0+dd0)*G2XT)
             tet2 = tet02*dexp(gam02/(1D0+dd0)*G2XT)

             harter = 1D0*r*t*dlog(1D0-dexp(-tet1/t)) + 
     *               2D0*r*t*dlog(1D0-dexp(-tet2/t))
         
         else

             tet = tet0*dexp(gam0/(1D0+dd0)*G2XT)

             harter = 3D0*nn*r*t*dlog(1D0-dexp(-tet/t))
         
         end if

      end function harter

      double precision function gterm2(id)
c----------------------------------------------------------------
c function evaluates 2nd part of modified EoS for liquid C after Brosh
c
c EOS includes 2 terms: G1 (gsgte + gbrosh) and G2; G2 converges to
c a constant at high p; the term is added to the eos of liquid C
c to model graphite-like behavior of liquid at lower p and diamond-
c like at high p
c
c parameters for the second part of eos: v02, b20 (B2), b2(B'2),
c tet2 (theta2), al2 (alpha2); Murnaghan EoS was used for the
c additional term
c in thermodynamic data file: c1 - v02; c2 - b20; c3 - b2; c4 - theta2;
c c5 - alpha2
c-----------------------------------------------------------------
      implicit none
      
      include 'perplex_parameters.h'

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer id

      double precision tet2, b2,t1, t2, b20
c----------------------------------------------------------------------
c                             allocate coefficients for 2nd term in EoS

      b2 = thermo(3,id)
      b20 = thermo(2,id)
      tet2 = thermo(4,id)
      
      t1  = dexp(- b2*thermo(5,id) * (t - dlog(1d0+t/tet2)*tet2) )
      t2  = 1d0-1d0/b2

      gterm2 = b20*thermo(1,id)/(b2-1d0)*((t1+b2*p/b20)**t2 - t1**t2)

      end function gterm2

      subroutine gpmelt (g,id)
c----------------------------------------------------------------------
c subroutine to speciation of the green et al (JMG, 2016) melt model. this
c model is a special case (ksmod(id)=27) because of the peculiar configurational
c entropy expression and because the model has a single ordering parameter, which 
c green et al take as the fraction of the ordered species (an). this formulation is 
c unfortunate because p(an) is not orthogonal to the disordered speciation
c (p0, because the moles of the species is not constant with changing speciation). 
c here the model is recast as g(p0,q) where q is the number of moles of an that can be 
c formed give p0. 

c    id identifies the solution.
c    g  is the change in G for the stable speciation relative to a mechanical
c       mixture of the endmembers.
c    pc is the mass normalization factor, sum(p0*ctot)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, i1, i2, id, jd, itic

      logical error, done 

      double precision g, qmax, qmin, q, dq, rqmax, d0

      double precision hpmelt, gex
      external hpmelt, gex

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dvnu,deph,dydy
      common/ cxt3r /dvnu(m4,j3,h9),deph(3,j3,h9),dydy(m4,j3,h9)

      logical pin
      common/ cyt2 /pin(j3)

      double precision enth
      common/ cxt35 /enth(j3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)
c----------------------------------------------------------------------
      i1 = ideps(1,1,id)
      i2 = ideps(2,1,id)
      jd = lstot(id) + 1 
      error = .false.
c                                 starting point:
c                                 find the maximum amount of the 
c                                 ordered species that can be formed
c                                 from the fully disordered species
c                                 fractions 
      rqmax = dmin1(p0a(i1),p0a(i2))
c                                 to avoid singularity set the initial 
c                                 composition to the max - nopt(5), at this
c                                 condition the first derivative < 0, 
c                                 and the second derivative > 0 (otherwise
c                                 the root must lie at p > pmax - nopt(5).
      if (rqmax.gt.nopt(5)) then

         pin(1) = .true.
         qmax = rqmax - nopt(5)
         qmin = nopt(5)
c                                 the p's are computed in gpderi
         call gpderi (i1,i2,jd,id,qmax,dq,g)

         if (dq.lt.0d0) then 
c                                 at the maximum concentration, the 
c                                 first derivative is positive, if 
c                                 the second is also > 0 then we're 
c                                 business
            q = qmax

         else
c                                 try the min

            call gpderi (i1,i2,jd,id,qmin,dq,g)

            if (dq.gt.0d0) then 
c                                 ok
               q = qmin

            else
c                                 no search from either limit possible
c                                 set error .true. to compare g at the 
c                                 limits.
               error = .true.
               goto 90

            end if 
         end if 
c                                 increment and check p
         call pcheck (q,qmin,qmax,dq,done)   
c                                 iteration counter to escape
c                                 infinite loops
         itic = 0
c                                 newton raphson iteration
         do 

            call gpderi (i1,i2,jd,id,q,dq,g)

            call pcheck (q,qmin,qmax,dq,done)
c                                 done is just a flag to quit
            if (done) then

               goodc(1) = goodc(1) + 1d0
               jcount(1) = jcount(1) + itic 
c                                 in principle the p's could be incremented
c                                 here and g evaluated for the last update.
               return

            end if

            itic = itic + 1

            if (itic.gt.iopt(21)) then
c                                 fails to converge.
               error = .true.
               badc(1) = badc(1) + 1
               jcount(1) = jcount(1) + itic
               exit

            end if 

         end do

      else
c                                 speciation is not stoichiometrically possible
         g = -t*hpmelt(id,p0a) + gex(id,p0a)
         return

      end if

90    if (error) then 
c                                 didn't converge or couldn't find a good 
c                                 starting point compare the fully ordered 
c                                 and disordered g's and choose the lowest
c                                 full disorder:
         dq = -t*hpmelt(id,p0a) + gex(id,p0a)
c                                 full order:
         do i = 1, jd

            if (i.eq.i1.or.i.eq.i2) then
               d0 = (p0a(i) - rqmax)
            else if (i.lt.jd) then 
               d0 = p0a(i)
            else
               d0 = rqmax
            end if 

            pa(i)  = d0/(1d0 - rqmax)

         end do 

         g = (pa(jd)*enth(1) - t*hpmelt(id,pa) + gex(id,pa)) * 
     *       (1d0 - rqmax)
c                                 compare
         if (g.gt.dq) then 
            g = dq
            do i = 1, jd 
               pa(i)  = p0a(i)
            end do 
         end if 

      end if

      end

      subroutine gpderi (kd,ld,jd,id,q,dg,g)
c----------------------------------------------------------------------
c subroutine to compute the newton-raphson increment (dg) in the ordering
c parameter from the 1st and 2nd derivatives of the g of the 
c green et al (JMG, 2016) melt model. kd and kd are the indices of the
c true dependent endmembers (woL and qL), jd is the index of the ordered
c species (anL). id is the solution. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, kd, ld, i1, i2, id, jd

      double precision g, dg, d2g, s, ds, d2s, pfac,
     *                 q, ph2o, dph2o, d2ph2o, pfo, pfa, pnorm, pnorm2,
     *                 dp(m4), d2p(m4), lpa(m4), lpfac, dng,
     *                 dsfo, dsfa, gnorm, dgnorm
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 excess energy variables
      integer jterm, jord, extyp, rko, jsub
      common/ cxt2i /jterm(h9),jord(h9),extyp(h9),rko(m18,h9),
     *               jsub(m2,m1,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision enth
      common/ cxt35 /enth(j3)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
c                                 initialize
      g   = 0d0 
      dg  = 0d0
      d2g = 0d0
 
      pfo = 0d0
      pfa = 0d0
c                                 the difficulty in this model is the
c                                 non-equimolar speciation reaction, this 
c                                 causes the number of moles of the components
c                                 in a mole of solution to change as a function
c                                 of the order parameter even if composition is
c                                 held constant. 

c                                 to keep the number of moles of the components
c                                 in the solution constant the gibbs energy
c                                 is multiplied by gnorm = 1 + sum(nu(i)), where
c                                 the nu(i) are the stoichiometric coefficients of
c                                 the endmembers in the ordering reaction (it being
c                                 assumed that nu(jd) = 1 and p0(jd) = 0). this gives 
c                                 the solutions g when it has the same amounts of the 
c                                 components as in the disordered limit (p = p0). the
c                                 amounts of the species (p) for a partially or completely
c                                 disordered state are p(i) = (p0(i) + nu(i))*q/gnorm.
c                                 q is the molar amount of the ordered species formed
c                                 by the ordering reaction from the amounts of the 
c                                 reactant species in the disordered limit. 

c                                 for the green et al melt model sum(nu(i)) for the
c                                 reaction wo + als = an is -1, therefore 
c                                 gnorm = (1 - q) and pnorm = 1/(gnorm)
      gnorm  = 1d0 - q
      dgnorm = -1d0
      pnorm = 1d0/gnorm
c                                 diff(pnorm,q,q)
      pnorm2 = 2d0*pnorm

      do i = 1, jd
c                                 calculate dp(i)/dq, d2p(i)/dq.
         if (i.eq.kd.or.i.eq.ld) then
c                                 reactive species
            pa(i)  = (p0a(i) - q) * pnorm
            dp(i)  = (pa(i) - 1d0) * pnorm

         else if (i.lt.jd) then
c                                 inert species
            pa(i)  = p0a(i) * pnorm
            dp(i)  = pa(i)  * pnorm

         else 
c                                 the ordered species
            pa(i)  = q*pnorm
            dp(i)  = (pa(i) + 1d0) * pnorm

         end if

         d2p(i) = dp(i)  * pnorm2

         if (pa(i).gt.0d0) lpa(i) = dlog(pa(i)) + 1d0 

      end do

      do i = 1, jterm(id)
c                                 assuming regular terms
        i1 = jsub(1,i,id)
        i2 = jsub(2,i,id)

        g = g + w(i) * pa(i1) * pa(i2)
        dg = dg + w(i) * (pa(i1)*dp(i2) + pa(i2)*dp(i1))
        d2g = d2g + w(i) * (      d2p(i1)* pa(i2) 
     *                      + 2d0*dp(i2) * dp(i1)
     *                      +     d2p(i2)* pa(i1) )

      end do  
c                                 get the configurational entropy derivatives
      if (jspec(id,1).eq.1.and.pa(1).gt.0d0) then
c                                 quadratic water term:
         ph2o = pa(1)
         dph2o = dp(1)
         d2ph2o = d2p(1)
c                                 pfac renormalizes the non-volatile
c                                 fractions (i.e., the multiplicity of
c                                 the non-volatike site is pfac):
         pfac = 1d0 - ph2o
         lpfac = dlog(pfac)  
c                                 entropy expression rewritten
c                                 to collect the pfac term
         s   =  -2d0 * ph2o * dlog(ph2o) - pfac*lpfac
         ds  =  (-2d0 * lpa(1) + (lpfac + 1d0)) * dph2o
         d2s =  (-2d0*lpa(1) + lpfac + 1d0) * d2p(1)
     *        - (2d0/ph2o + 1d0/pfac)*dph2o*dph2o

      else 

         s = 0d0 
         ds = 0d0 
         d2s = 0d0 

      end if 
 
c                                 the configurational entropy
      do i = jspec(id,4), jd

         if (pa(i).le.0d0) cycle
         s = s - pa(i) * dlog(pa(i))
         ds = ds - lpa(i)*dp(i)
         d2s = d2s - (lpa(i)*d2p(i) + dp(i)*dp(i)/pa(i))
 
      end do 
c                                 the fe-mg fudge factor
      if (jspec(id,2).ne.0) pfo = pa(jspec(id,2))   
      if (jspec(id,3).ne.0) pfa = pa(jspec(id,3))

      pfac = pfo + pfa

      if (pfo.ne.0d0.and.pfa.ne.0d0) then 

         dsfo = 4d0 * dlog(pfo/pfac)
         dsfa = 4d0 * dlog(pfa/pfac)

         s  = s  - pfo*dsfo
         s  = s  - pfa*dsfa
         ds  = ds  - dsfo * dp(jspec(id,2))
         ds  = ds  - dsfa * dp(jspec(id,3))
         d2s = d2s - dsfo * d2p(jspec(id,2))
         d2s = d2s - dsfa * d2p(jspec(id,3))

      end if 

      g   = g   + enth(1)*pa(jd)  - r*v(2)*s
      dg  = dg  + enth(1)*dp(jd)  - r*v(2)*ds
c                                 the normalized g derivative
      dng  = g * dgnorm + gnorm * dg
c                                 and second derivative
      d2g = gnorm * (d2g + enth(1)*d2p(jd) - r*v(2)*d2s) 
     *       + 2d0 * dg * dgnorm
c                                 dg becomes the newton-raphson increment:
      dg = -dng/d2g
c                                 and g the normalized g:

      g   = g * gnorm

      end

      subroutine chopit (ycum,fac,jst,jsp,lsite,ids,jump,extra)
c---------------------------------------------------------------------
c subroutine to do cartesian or transform subdivision of species
c jst+1 through jsp on site k of solution ids. ycum is the smallest
c fraction possible (i.e., if the minimum bound for some species 
c is > 0). the fractions are loaded into simp.

c extra - save space for an extra coordinate in simp (ksmod 20 or 9)
c fac   - factor to modify default resolution (1, except ksmod 9)
c jump  - offset for storing coordinates in simp (ksmod 9)
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer mres
 
      parameter (mres=12000)

      logical extra

      integer mode, ind(ms1), iy(ms1), jsp, lsite, indx, iexit, 
     *        ieyit, i, j, k, ids, ico, jst, jump

      double precision y(ms1,mres), ycum, ymax, dy, ync, 
     *                 x, unstch, strtch, delt, dx, fac

      external unstch, strtch

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
c----------------------------------------------------------------------
      if (.not.extra) then
c                                 chopit always generates jsp coordinates
         ico = jsp
      else 
c                                 but in the case of charge balance save 
c                                 space for the dependent coordinate.
         ico = jsp + 1
      end if

      do i = 1, jsp

         k = jst + i

         if (ksmod(ids).eq.20) then 
c                                 electrolyte model skip limit on ns'th
c                                 species
            if (k.eq.ns) k = k + 1

         end if 
c                                 generate coordinates for i'th component
         iy(i) = 1
         y(i,1) = xmn(lsite,k)

         ync = xnc(lsite,k)/fac

         if (ync.eq.0d0) cycle

         mode = imdg(k,lsite,ids)
c                                 avoid impossible compositions 'cause a min > 0
         if (i.gt.1) then 

            ycum = ycum + xmn(lsite,k-1)
c                                 1-ycum is the smallest fraction possible
            if (ycum.gt.1d0) then 
c                                 inconsistent limits
               write (*,'(/,a,/)') '#########BOOM WACKA BOOM###########'
               write (*,*) ycum,ids,ksmod(ids),lsite,k,i,mode
               write (*,*) (xmx(1,j),j=1,jsp)
               write (*,*) (xmn(1,j),j=1,jsp)
               write (*,*) (xnc(1,j),j=1,jsp)
               write (*,*) (iy(j),j=1,jsp)
               call warn (999,ycum,jsp,'cartes')

               cycle

            else
c                                 the smallest fraction possible is lt
c                                 than xmax
               ymax = xmx(lsite,k)

            end if 
         else 
            ymax = xmx(lsite,k)
         end if 
c                                 two means of extracting y-range, cartesian
c                                 imod = 0 and transformation imod = 1
         if (mode.eq.0) then 
c                                 cartesian
            delt = nopt(5)

            do

               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))

               y(i,iy(i)) = y(i,iy(i)-1) + ync

               if (dabs(y(i,iy(i))-ymax).lt.delt.or.
     *             y(i,iy(i)).gt.ymax) then
                  y(i,iy(i)) = ymax
                  exit 
               end if

            end do

         else 
c                                 see 6.8.0 for old multiple interval conformal transformation
c                                 y is the non-linear cartesian coordinate
c                                 x is the linear conformal coordinate.
            call setstc (ids,lsite,k)

            delt = xmno(ids,lsite,k)
            if (delt.gt.nopt(5)) delt = nopt(5)

            x = unstch (xmn(lsite,k))
            dx = 1d0/ync

            do

               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))

               x = x + dx 
               y(i,iy(i)) = strtch (x)

               if (dabs(y(i,iy(i))-ymax).lt.delt.or.
     *             y(i,iy(i)).gt.ymax) then

                  y(i,iy(i)) = ymax
                  exit

               end if

            end do

         end if

      end do
c                                 the first coordinate
      npairs = 1

      do i = 1, jsp
         ind(i) = 1
         simp(jump+i) = y(i,1)
      end do
c                                 now make the array index run over all
c                                 values increasing the last index fastest
      iexit = 0 
      ieyit = 0 
      dy = 0d0

      do while (iexit.eq.0)
c                                 figure out which index to increment
         do i = jsp, 1, -1

            if (ind(i).lt.iy(i).and.ieyit.eq.0) then
c                                 this is the one to increment
               ind(i) = ind(i) + 1
               indx = i 
               exit 

            else if (i.gt.1) then 
c                                 saturated the index
               ind(i) = 1
               ieyit = 0 
               
            else
c                                 saturated first index, done.
                return 

            end if
 
         end do 
c                                 ok now we have the indices, check
c                                 the composition
         ycum = 0d0
         dy = 0d0 

         do i = 1, jsp 
            ycum = ycum + y(i,ind(i))  
         end do 

         if (ycum.gt.1d0) then
c                                 no matter what this is the last point
c                                 to be generated for ind(indx), set ieyit
c                                 to change indx
            ieyit = 1
c                                 but here is where it gets messy:
            if (indx.eq.1) then
c                                 we're at the first point, and already
c                                 over the top
               iexit = 1
               cycle

            else if ( y(indx,ind(indx)) - y(indx,ind(indx)-1)
     *               - ycum + 1d0    .gt. delt ) then
c                                 the excess (ycum-1) is less then the 
c                                 amount the variable was previously incremented
c                                 so it's possible to back off the composition
c                                 to zero excess. this will always be true for
c                                 cartesian transformations, it might not be true
c                                 for conformal stretching (i.e., increments could
c                                 grow in the direction of the nodal index).
               dy =  1d0 - ycum

            else
c                                 must have just hit on the last increment or
c                                 conformal. 
               cycle 

            end if 

         end if 

         npairs = npairs + 1
         j = jump + (npairs-1)*ico

         if (j+jsp.gt.k13) call error (180,ycum,k13,
     *                      'CARTES increase parameter k13')

         do i = 1, jsp
            simp(j+i) = y(i,ind(i))
         end do 

         simp(j+indx) = simp(j+indx) + dy

      end do

      end 

      subroutine cartaq (ids,resub)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on a single site
c solution with charge balance. called by subdiv. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical resub

      integer i, j, k, l, m, n, ids, qpairs, np0

      double precision ycum, sum, q, ratio

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer jend
      common/ cxt23 /jend(h9,m4)
c----------------------------------------------------------------------
c                                 could save only nqs-2 coordinates, but 
c                                 will try nqs-1 first.
      ycum = 0d0

      if (ns1.eq.0) then
c                                 only solvent, test for no solvent has
c                                 already been made in reform
         do j = 1, nqs1
            prism(j) = 0d0
         end do

         npairs = 1

      else
c                                 subdivision of neutral ns+nn-1 species
         call chopit (ycum,1d0,0,ns1,1,ids,0,.true.)

         do i = 1, npairs

            k = (i-1)*sn
            l = (i-1)*nqs1

            sum = 0d0

            do j = 1, ns1
               prism(l+j) = simp(k+j)
               sum = sum + simp(k+j)
            end do

            if (nq.gt.0) prism(l+nqs1) = sum

         end do

      end if

      ntot = npairs

      if (nq.ne.0) then

         np0 = npairs
c                                 do the nq-1 species independently
         ycum = 0d0

         call chopit (ycum,1d0,sn,nq1,1,ids,0,.true.)
c                                 at this point simp contains all 
c                                 possible compositions of the nq-1 species,
c                                 use charge balance to get the nqth species
         qpairs = 1 

         do i = 1, npairs

            q = 0d0
            sum = 0d0

            k = (i-1)*nq
            l = (qpairs - 1)*nq
            m = 2 + sn

            do j = 1, nq1
               q = q + thermo(6,jend(ids,m+j))*simp(k+j)
               sum = sum + simp(k+j)
               simp(l+j) = simp(k+j)
            end do
c                                 charge ratio
            ratio = q/thermo(6,jend(ids,m+j))
c                                 the net charge has the same sign as the nqth
c                                 species or its amount violates closure, reject:
            if (ratio.gt.0d0.or.sum-ratio.ge.1e0) cycle
c                                 the amount of the species determined by charge balance
            simp(l+nq) = -ratio

            qpairs = qpairs + 1

         end do

         qpairs = qpairs - 1
c                                 for every charged species composition
c                                 load all neutral compositions that don't
c                                 violate closure:
         do i = 1, qpairs
c                                 get the sum of the charged species
            k = (i-1)*nq

            sum = 0d0

            do j = 1, nq 
               sum = sum + simp(k+j)
            end do
c                                 now assemble full compositions:
            do j = 1, np0

               l = (j-1)*nqs1
c                                 test for closure
               if (prism(l+nqs1)+sum.ge.1d0) cycle
c                                 acceptable composition
               m = ntot * nqs1 
               if (m+nqs1.gt.k24) call errk24 (resub) 

               ntot = ntot + 1
c                                 load neutral part
               do n = 1, ns1
                  prism(m+n) = prism(l+n)
               end do
c                                 load charged part
               do n = 1, nq
                  prism(m+ns1+n) = simp(k+n)
               end do

            end do

         end do
c                                  zero the charged species coordinates
c                                  for the first np0 neutral compositions
         do i = 1, np0

            l = (i-1)*nqs1

            do j = sn, nqs1
               prism(l+j) = 0d0
            end do

         end do

      end if 

      end

      subroutine setsol (ids, wham)
c-----------------------------------------------------------------------
c load species indices, charges, etc for aqueous model (ksmod = 20) into
c simple arrays
c-----------------------------------------------------------------------
      implicit none

      logical wham
 
      include 'perplex_parameters.h'

      integer i, ids

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
      if (wham) then 
c                                 an internal molecular eos has already 
c                                 been invoked 
         call error (77,rt,i,'only one solution model may invoke an '//
     *                       'internal molecular fluid EoS.')

      else 

         wham = .true.

      end if 
c                                 load endmember indices into a simple array
      do i = 1, mstot(ids)
         jnd(i) = jend(ids,2+i)
      end do 
c                                 model uses an internal speciation routine
c                                 set the ins array and isp pointer
      if (ksmod(ids).eq.0) then
c                                 fluid eos specified via ifug 
         call setins (ifug)

      else if (ksmod(ids).eq.20) then 
c                                 electrolyte
         do i = 1, nqs
            q2(i) = thermo(6,jnd(i))**2
         end do 

         isp = ns

         do i = 1, isp
            ins(i) = jspec(ids,i)
         end do

         na1 = 1

      else if (ksmod(ids).eq.39) then
c                                 hybrid molecular
         isp = mstot(ids)
c                                 set ns in case of aqrxdo or aqlagd
         ns = isp 
         ns1 = ns - 1
         sn1 = ns + 1
         nsa = ns + aqct
         na1 = nsa + 1
         na2 = nsa + 2
         na3 = nsa + 3
         nat = nsa + 8

         do i = 1, isp
            ins(i) = jspec(ids,i)
         end do

      else if (ksmod(ids).eq.40) then 
c                                 MRK silicate vapor (40), EoS code 26
         call setins (26)

      else if (ksmod(ids).eq.41) then 
c                                 MRK COH fluid (41), EoS code 27
         call setins (27)

      end if 

      end 

      subroutine slvnt1 (gsolv)
c-----------------------------------------------------------------------
c computes solvent p-t-composition dependent properties: dielectric cst (eps), 
c molar mass, debye-hueckel (adh), and gHKF function.
c
c assumes pure species volumes and fugacity coefficients have been calculated.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, k

      double precision gsolv, vsolv, cdh, ysum, ysolv(nsp), hyvol

      double precision gfunc, ghybrid, gcpd

      external gfunc, ghybrid, gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

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

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision gh,dvhy,gmrk0
      common/ csthyb /gh(nsp),dvhy(nsp),gmrk0(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

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
c                                 compute and add in solvent activities by calling
c                                 mrkmix ghybrid gets pmvs for hybrid volume calculation
      gsolv = gsolv + ysum*( ghybrid (ysolv) + rt*dlog(ysum) )
c                                 compute pmvs required for dielectric constant 
      hyvol = 0d0

      do k = 1, ns

         i = ins(k)
c                                 hybrid pmv 
         vhyb(i) = dvhy(i) + v(i)
c                                 hybrid total volume
         hyvol = hyvol + yf(i) * vhyb(i)

      end do 

      do k = 1, ns

         i = ins(k)
c                                 volume fractions
         vf(i) = yf(i) * vhyb(i) / hyvol

      end do

      vsolv = ysum * hyvol
c                                 get dielectric cst based on hybrid volumetric properties
      call geteps (epsln)
c                                 set reference dielectric cst, quick fix for mistake
c                                 of having made it composition dependent,
      epsln0 = 78.47d0
c                                 Debye-Hueckel factor, A[cgs] = -q^3*sqrt(NA)/(4*Pi*k^(3/2))
c                                 *(msolg/(10*vsoljbar))^(1/2)/(epsilon*T)^(3/2) for ln(gamma) = +A*....
c                                 A = cdh*(msolkg/(vsoljbar))^(1/2)/(epsilon*T)^(3/2)
      adh = cdh * dsqrt(1d1*msol/vsolv/(epsln*t)**3)
c                                 shock et al 1992 g function (cgs solvent density),
c                                 used by hkf
      gf = gfunc (msol*1d3/vsolv) 

      end

      subroutine slvnt2 (gsolv)
c-----------------------------------------------------------------------
c computes: debye-hueckel (adh) and solute contribution of the fluid 
c gibbs energy. assumes the molar speciation is stored in y. called 
c only for solution model type 20. 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k

      double precision gsolv, mo(m4), lng0, is

      double precision gcpd, aqact

      external gcpd, aqact

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

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
      is = 0d0 
c                                 molalities and ionic strength
      do k = sn1, nqs
c                                 ln molality of solutes
         mo(k) = y(k)/msol
         is = is + q2(k) * mo(k)

      end do
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
c                                 Davies extension.
      lng0 = dlog(aqact(is/2d0))
c                                 add in the solute gibbs energies
      do k = sn1, nqs

         if (y(k).le.0d0) cycle

         gsolv = gsolv + y(k) * 
     *           (gcpd(jnd(k),.true.) + rt*(dlog(mo(k)) + lng0*q2(k)))

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

      double precision gso(nsp), gcpd

      external gcpd

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision fwt
      common/ cst338 /fwt(k10)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
c                                 solvent properties
      call slvnt3 (gso,.true.,.true.,id)
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

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision y,g,v
      common/ cstcoh /y(nsp),g(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)
c----------------------------------------------------------------------
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
         rho = 1d0/vhyb(j)
c                                 Eq 5 of H&L 2005 for  polarization/rho
c                                 for polar species need to add po(j,3)
         eps = po(j,1) + po(j,2)*trt + (po(j,4) + po(j,5)*trt)*rho
     *                + (po(j,6) + po(j,7)*trt)*rho**po(j,8)

         if (j.eq.6) then 
c                                 fake for h2s = 1/3 h2o
            eps = epsh2o (vhyb(j)/1d1) / 3d0

         else if (po(j,3).eq.0d0) then 
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
      epsln  = epsln + vf(ins(i)) * epsh2o (vhyb(ins(i))/1d1)**(1d0/3d0)

      epsln  = epsln **3

      end 

      subroutine slvntg (gso,mu)
c-----------------------------------------------------------------------
c given chemical potentials compute solvent species partial molar gibbs
c energies for aqrxdo/aqrxlg
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j

      double precision gso(nsp), mu(k8)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------

         do i = 1, ns 

            gso(i) = 0d0

            do j = 1, icp

               if (isnan(mu(j))) cycle

               gso(i) = gso(i) + mu(j)*cp(j,jnd(i))

            end do 

         end do

      end 

      subroutine aqrxdo (jd,lu)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation

c   jd - is the pointer to the solvent in the local assemblage 
c   lu - is the LUN for output, also signals type of output, console
c        vs tab (<0).
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, l, ind(l9), lu, jd, badct

      logical bad, output

      character text*200

      double precision mo(l9), blk(k5), dn, smot, err, is, gamm0, totm,
     *                 g0(l9), lnkw, gso(nsp), ph0, tmass, tsmas, tsmol,
     *                 smol(k5), errkw, smo, gcpd

      external gcpd

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

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

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

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

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)

      save badct
      data badct/0/
c----------------------------------------------------------------------
      if (.not.mus) then 

         call warn (99,0d0,0,
     *       'no chemical potentials, cannot back-calculate solute '//
     *       'speciation')

         return

      end if 

      output = .true.

      if (jdaq.ne.0.and.lopt(32)) then
c                                 a multi species solvent is present: 
         do i = 1, ns
c                                 doing lagged speciation, set the 
c                                 solvent mole fractions to the true 
c                                 values
            ysp(i,jd) = caq(jd,i)

         end do

      end if
c                                 set feos = .true. because can't be 
c                                 sure that the last solvent calculation was
c                                 at the present p-t condition.
      call slvnt3 (gso,.true.,.true.,jd)
c                                 slvnt3 doing lots needless calculations
      call slvntg (gso,mu)
c                                 iterate on speciation
      call aqsolv (g0,gso,mo,mu,is,gamm0,lnkw,bad)
c                                 back calculated bulk composition
      if (bad) then 

         if (badct.lt.11) then 

            badct = badct + 1 

            call warn (99,0d0,0,'AQRXDO did not converge on solute '//
     *                          'speciation')

            if (badct.eq.10) call warn (49,0d0,99,'AQRXDO')

         end if 

         do i = 1, iprop
            prop(i) = nopt(7)
         end do 
 
         return

      end if 
c                                 compute charge balance error
      err = 0d0

      do i = 1, ichg 
         err = err + q(jchg(i)) * mo(jchg(i))
      end do
c                                 neutral pH
      ph0 = -lnkw/2d0/2.302585d0

      do i = 1, icp
         blk(i) = 0d0
      end do
c                                 total molality
      smot = 0d0
c                                 compute mole fractions, total moles first
      do i = 1, ns 
c                                 moles/kg-solvent 
         do j = 1, icp 
            blk(j) = blk(j) + y(i)*cp(j,jnd(i))/msol
         end do

        smot = smot + y(i)/msol

      end do

      smo = 0d0 

      do i = 1, aqct
c                                  total solute molality
         smo = smo + mo(i)
         ind(i) = i 

         do j = 1, icp
            blk(j) = blk(j) + mo(i)*aqcp(j,i)
         end do 

      end do

      smot = smot + smo
c                                 bulk fluid composition 
      tmass = 0d0
      totm = 0d0 

      do i = 1, icp
         totm = totm + blk(i)
         tmass = tmass + atwt(i)*blk(i)
      end do 
c                                error in log10(K_w)
      errkw = (-lnkw + dlog(mo(ihy)*mo(ioh)*gamm0**2
     *             ))/2.302585d0

      if (output.and.lu.lt.0) then
c                                 WERAMI props on a grid
         do i = 1, icp
c                                 bulk composition
            if (lopt(23)) then

               prop(i) = blk(i)/totm*1d2

            else 

               prop(i) = blk(i)*atwt(i)/tmass*1d2

            end if 

         end do

         k = icp

         do i = 1, ns 
c                                 solvent speciation
            k = k + 1

            if (lopt(26)) then
c                                 mole fraction
               prop(k) = y(i)/msol/smot
            else 
c                                 molality
               prop(k) = y(i)/msol
            end if  

         end do 

         do i = 1, aqct 
c                                 solute speciation
            k = k + 1

            if (lopt(27)) then
c                                 molality
               prop(k) = mo(i)
            else 
c                                 mole fraction
               prop(k) = mo(i)/smot 
            end if  

         end do 
c                                  other properties
         prop(k+2) = -dlog10(mo(ihy)*gamm0)
         prop(k+1) = prop(k+2) - ph0
         prop(k+3) = dabs(errkw)
         prop(k+4) = epsln
         prop(k+5) = is
         prop(k+6) = smo

      else if (output) then 
c                                 WERAMI/MEEMUM console output 
         call rankem (mo,ind,aqct,iopt(32))

         write (lu,1000)

         write (text,1050) -dlog10(mo(ihy)*gamm0),dabs(errkw),
     *                     -dlog10(mo(ihy)*gamm0)-ph0,is,gamm0
         call deblnk (text)
         write (lu,'(400a)') (chars(j), j = 1, length)
         write (text,1070) epsln,msol*1d3,smo
         call deblnk (text)
         write (lu,'(400a)') (chars(j), j = 1, length)

         if (jdaq.eq.20.or.lopt(32).and.jdaq.eq.39) then

            write (lu,1100)

         else

            write (lu,1040)

         end if 

         do i = 1, iopt(32)

            k = ind(i)

            if (jdaq.eq.20) then 
c                                 compare to forward speciation:
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
     *                        mo(k),mo(k)/smot,ysp(l,jd),
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))
               else

                  write (lu,1090) aqnam(k),int(thermo(6,k+aqst)),
     *                        mo(k),mo(k)/smot,
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))
               end if 

            else if (lopt(32)) then
c                                 compare to lagged speciation:
                  write (lu,1080) aqnam(k),int(thermo(6,k+aqst)),
     *                        mo(k),mo(k)/smot,caq(jd,k+ns)/caq(jd,na2),
     *                        int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                        int(g0(k))

            else 
c                                 only back-calculated result:
               write (lu,1010) aqnam(k),int(thermo(6,k+aqst)),
     *                         mo(k),mo(k)/smot,
     *                         int(g0(k)+rt*(dlog(mo(k)*gamm0**q2(k)))),
     *                         int(g0(k))

            end if 
         end do 

         if (lopt(32)) then 
c                               lagged
            write (lu,1140)

            do i = 1, ns 

               write (lu,1150) names(jnd(i)), y(i)/msol, y(i)/msol/smot,
     *                         ysp(i,jd), vf(ins(i)), vhyb(ins(i)), 
     *                         int(gso(i)), int(gcpd(jnd(i),.true.))

            end do

         else 

            if (jdaq.eq.20) then 
c                               direct 
               write (lu,1020)
            else 
c                               back-calculated, indicate normalized
               write (lu,1160)

            end if 

            do i = 1, ns 

               write (lu,1030) names(jnd(i)), ysp(i,jd)/msol, ysp(i,jd),
     *                         vf(ins(i)), vhyb(ins(i)), 
     *                         int(gso(i)), int(gcpd(jnd(i),.true.))

            end do 

         end if 
c                                 footnotes:
         if (lopt(32).or.jdaq.eq.20) then 

            write (lu,1170)

         else 

            write (lu,1060)

         end if 
c                                 bulk chemistry:
         if (jdaq.eq.20.or.lopt(32)) then 

            tsmol = 0d0
            tsmas = 0d0
 
            do i = 1, icp
               smol(i) = 0d0
            end do 

            if (lopt(32)) then
c                                 lagged model
               do i = 1, nsa

                  do j = 1, icp

                     if (i.lt.sn1) then 
                        dn = y(i)/caq(jd,na3) * cp(j,jnd(i))
                     else 
                        dn = caq(jd,i) * aqcp(j,i-ns)
                     end if 

                     smol(j) = smol(j) + dn

                     tsmol = tsmol + dn
                     tsmas = tsmas + dn*atwt(j)

                  end do

               end do 

            else 
c                                 forward model
               do i = 1, nqs
               
                  k = jnd(i)

                  do j = 1, icp

                     if (i.lt.sn1) then 
                        dn = ysp(i,jd) * cp(j,k)
                     else 
                        dn = ysp(i,jd) * aqcp(j,k-aqst)
                     end if 

                     smol(j) = smol(j) + dn

                     tsmol = tsmol + dn
                     tsmas = tsmas + dn*atwt(j)

                  end do

               end do 

            end if 

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
     *        9x,'molality',3x,'mol_fraction',2x,'vol_fraction#',
     *        1x,'v,cm3/mol*',2x,'g,J/mol*',3x,'g0,J/mol***')
1030  format (a8,2x,f7.4,5x,f7.5,8x,f7.5,5x,f7.3,3x,i8,4x,i8)
1040  format (/,'Solute endmember properties:',//,10x,'charge',3x,
     *       'molality',5x,'mol_fraction',6x,'g,J/mol*',5x,'g0,J/mol**')
1050  format ('pH = ',f6.3,'+/-',f5.3,
     *        '; Delta_pH = ',f6.3,'; ionic_strength = ',
     *        g10.4,'; gamma/q^2 = ', g10.4)
1060  format (/,'*partial molar, **molal ref. state, ***molar ref. ',
     *        'state.',/)
1070  format ('dielectric constant = ', g10.4,
     *        '; solvent molar mass, g/mol = ',f8.4,
     *        '; solute molality = ',g10.4)
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
1140  format (/,'Solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'opt_mol_frac',
     *        2x,'vol_fraction# v,cm3/mol*',2x,'g,J/mol*',
     *        3x,'g0,J/mol***')
1150  format (a8,2x,f7.4,5x,2(f7.5,8x),f7.5,5x,f7.3,3x,i8,4x,i8)
1160  format (/,'Normalized solvent endmember properties:',//,
     *        9x,'molality',3x,'mol_fraction',2x,'vol_fraction',
     *        1x,'v,cm3/mol*',2x,'g,J/mol*',3x,'g0,J/mol***')
1170  format (/,'#normalized, *partial molar, **molal ref. state, ',
     *        '***molar ref. state.',/)
      end 

      subroutine rankem (a,ind,right,n)
c-----------------------------------------------------------------------
c rank the n largest values of array a(left:right) in array ind. assumes ind has 
c been initialized (left:right). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, imax, ind(*), left, right, n

      double precision a(*), amax
c----------------------------------------------------------------------
      left = 1

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

         f = 0d0
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

         if (x.le.0d0.or.it.gt.iopt(21)) then 
            bad = .true.
            exit
         else if (dx/x.lt.nopt(5)) then 
            bad = .false.
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

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

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

      subroutine aqidst
c-----------------------------------------------------------------------
c identify the aqueous phase for aqrxdo/aqlagd
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k

      logical lagged

      double precision tot

      character name*100

      character prject*100, tfname*100
      common/ cst228 /prject,tfname

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer isoct
      common/ cst79 /isoct

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer eos
      common/ cst303 /eos(k10)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
c                                 set option flags if necessary
      if (lopt(25).or.lopt(32)) then 

         if (jpot.ne.0) then 

            call warn (99,0d0,0,'aq_output or aq_lagged_speciation is T'
     *                        //', but dependent_potentials is off, '//
     *                         'dependent_potentials set = on (AQIDST)')

            jpot = 0

         end if 
c                                reset iopt(32) [# aq species output]
         if (iopt(32).gt.aqct) iopt(32) = aqct

      else 

         aqct = 0
         iopt(32) = 0
         return 

      end if 

      jdaq = 0
      lagged = .false.
c                                 look among solutions:
      do i = 1, isoct

         if (ksmod(i).eq.20.or.ksmod(i).eq.39) then

            idaq = i
            jdaq = ksmod(i)

            if (lopt(32)) then

               do j = 1, ns
c                                 set quack flag so the pure endmembers
c                                 won't be speciated by aqlagd
                  quack(jnd(j)) = .true.

               end do
c                                 identify non-solvent components
               isolc = 0 

               do j = 1, icp

                  tot = 0d0 

                  do k = 1, ns 
                     tot = tot + cp(j,jnd(k)) 
                  end do 

                  if (tot.gt.0d0) cycle

                  isolc = isolc + 1
                  solc(isolc) = j

               end do

               lagged =.true.

            end if

          end if

      end do

      if (jdaq.eq.0) then
c                                 no solution model found:
c                                 turn off lagged speciation just to be sure
         lopt(32) = .false. 

        if (.not.lopt(25)) aqct = 0

c                                 else look for H2O
         do i = 1, ipoint

            if (eos(i).eq.101) then
               idaq = -i
               jnd(1) = i
c                                 set solvent/species pointers on the
c                                 off chance they will be used
               ns = 1
               ins(1) = 1
               isp = 1
               return
            end if
         end do

      end if
c                                open a bad point file for lagged and
c                                back-calculated speciation calculations
      if (lagged.and.iam.le.2) then

         if (iam.eq.1) then 
            call mertxt (name,prject,'.pts',0)
         else 
            call mertxt (name,prject,'_MEEMUM.pts',0)
         end if 
         open (n13,file=name)
c                                meemum/vertex
       else if (.not.lagged.and.iam.eq.3.and.lopt(25)) then 
c                                werami back-calc
          call mertxt (name,prject,'_WERAMI.pts',0)
          open (n13,file=name)

      end if

      end 

      subroutine slvnt0 (gsolv,vsolv)
c-----------------------------------------------------------------------
c sets solvent p-t dependent properties for pure water: dielectric cst (eps), 
c molar mass, debye-hueckel (adh), and gHKF function.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision cdh, gsolv, vsolv

      double precision gfunc, gcpd, epsh2o, lnfpur

      external gfunc, gcpd, epsh2o, lnfpur

      double precision vol
      common/ cst26 /vol

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision vmrk0, vhyb, vf
      common/ cxt38 /vmrk0(nsp),vhyb(nsp),vf(nsp)

      integer iam
      common/ cst4 /iam

      save cdh
      data cdh/-5661800.47810d0/
c-----------------------------------------------------------------------
      if (iam.ne.5) then 
c                                 not frendly
         gsolv  = gcpd(jnd(1),.true.)

      else 
c                                 frendly
         gsolv = lnfpur(101)

      end if 
c                                 calling gcpd will get the molar volume (cm3/mol) 
c                                 in cst26 as vol
      vsolv = vol 
      msol   = 18.01528d-3
      epsln0 = 78.47d0

      epsln  = epsh2o (vol/1d1)
c                                 Debye-Hueckel factor, A[cgs] = -q^3*sqrt(NA)/(4*Pi*k^(3/2))
c                                 *(NH2O/(10*vh2o))^(1/2)/(epsilon*T)^(3/2) for ln(gamma) = +A*....
c                                 for reasons of stupidity this is set up for v in j/bar
      adh = cdh/dsqrt(vol/1d1*(epsln*t)**3)
c                                 shock et al 1992 g function (cgs solvent density),
c                                 used by hkf
      gf = gfunc (msol*1d3/vol)

      yf(1) = 1d0
      vf(1) = 1d0 

      end

      subroutine outtit
c-----------------------------------------------------------------------
c outtit writes title information and a brief description of the
c chemical system for each calculation requested.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i, j

      character*162 title
      common/ csta8 /title(4)
 
      integer iasmbl
      common/ cst27  /iasmbl(j9)

      double precision ctot
      common/ cst3   /ctot(k1)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character cname*5
      common/ csta4  /cname(k5) 

      character*8 names
      common/ cst8 /names(k1)

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common/ cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
c----------------------------------------------------------------------
      write (n3,1000)
c                          title:
      write (n3,1190) title(1)
c                          data base
      write (n3,1210) dname
c                          fluid
      if (ifct.gt.0.or.gflu) call rfluid (2,ifug) 
c                          independent potentials:
      write (n3,1070) (vname(jv(j)), j = 1, ipot)
c                          saturated phase components:
      if (ifct.gt.0) then 
         j = icp + isat
         write (n3,1200) (cname(j+i), i = 1, ifct)
      end if 
c                          saturated components:
      if (isat.gt.0) then 
         j = icp + isat
         write (n3,1180) (cname(i), i = icp1, j)
      end if 
c                          unconstrained components
      write (n3,1080) (cname(i), i = 1, icp)
c                          phases
      if (icp.gt.3) then 
         write (n3,1150) (cname(i), i = 1, icp)
         do i = istct, iphct
            write (n3,'(3x,a,12(1x,f5.3,1x))') 
     *            names(i), (cp(j,i)/ctot(i), j = 1, icp)
         end do 
      else if (icp.eq.3) then 
         write (n3,1090) (cname(i), i = 2, 3)
         write (n3,'(3(1x,a,2x,f5.3,2x,f5.3,5x))') 
     *         (names(i), cp(2,i)/ctot(i), cp(3,i)/ctot(i),
     *                                             i = istct, iphct)
      else if (icp.eq.2) then 
         write (n3,1040) cname(2)
         write (n3,'(4(2x,a,1x,f5.3))') 
     *         (names(i), cp(2,i)/ctot(i), i = istct, iphct)
      else if (icp.eq.1) then 
         write (n3,1130)
         write (n3,'(7(1x,a,1x))') (names(i), i = istct, iphct)
      end if
c                          saturation composant phases
      if (isat.ne.0) write (n3,'(/,a,/)')
     *          'Phases on saturation and buffering surfaces:'

      do i = 1, isat

         write (n3,'(/,3a,/)') ' for component ', cname(i+icp),':'
         write (n3,'(7(1x,a,1x))') (names(ids(i,j)), j = 1, isct(i))

      end do 
c                          excluded phases
      if (ixct.ne.0) then 
         write (n3,1140)
         write (n3,'(7(1x,a,1x))') (exname(i), i = 1, ixct)
      end if 
 
      write (n3,1000)
 
1000  format (/,80('-'),/)
1040  format (/,'Phases and (projected) mol fraction ',a,':',/)
1070  format (/,'Independently constrained potentials:',//,3x,8(a,1x))
1080  format (/,'Components with unconstrained potent'
     *       ,'ials:',//,3x,10(a5,3x))
1090  format (/,'Phases and (projected) composition with respect to '
     *         ,a5,' and ',a5,':',/)
1130  format (/,'Phases:',/)
1140  format (/,'Excluded phases:',/)
1150  format (/,'Phases and (projected) compositions:',//,
     *        11x,12(1x,a5,1x),/)
1180  format (/,'Saturated or buffered components:',//,3x,7(a,3x))
1190  format (/,'Problem title: ',a,/)
1200  format (/,'Saturated phase components:',//,3x,5(a,3x))
1210  format ('Thermodynamic data base from: ',a)
      end

      recursive double precision function gphase (id)
c-----------------------------------------------------------------------
c gphase computes the gibbs free energy of static compound id.
c gphase does not assume that the properties of a pseudocompound endmembers
c are known (unlike gall) it is thus less efficient than gall.

c used only (i think) for mixed variable and schreinemaker's projections.
c alloy solution models are commented out. these must be reinstated for the 
c aforementioned diagram types. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k,id,ids

      double precision gzero, dg, gerk, gph, gproj, gcpd, gfesi,
     *                 gfecr1, gfesic

      external gzero, gerk, gproj, gcpd, gfesi, gfecr1, gfesic

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer ikp
      common/ cst61 /ikp(k1)

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 new global arrays, 10/25/05:
c                                 bookkeeping variables
      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
      ids = ikp(id)
 
      if (id.le.ipoint) then
c                                 a pure compound
         gph = gcpd (id,.true.)

      else if (lorder(ids)) then
c                                 get composition
         call setxyp (ids,id)
c                                 initialize gph and any dqf corrections
         call gexces (id,gph)
c                                 compute margules coefficients
         call setw (ids)
c                                 evaluate enthalpies of ordering
         call oenth (ids)
c                                 get the speciation energy effect     
         call specis (dg,ids)
c                                 get gmech
         do k = 1, lstot(ids)
            gph = gph + gproj(jend(ids,2+k)) * p0a(k)
         end do 

         gph = gph + dg

      else if (ksmod(ids).eq.0) then
c                              get the excess and/or ideal mixing effect
c                              and/or dqf corrections:
         call fexces (id,gph)

         call setxyp (ids,id)
c                              excess props don't include vdp:
         do k = 1, nstot(ids) 
            gph = gph + gzero (jend(ids,2+k)) * y(k)
         end do 

      else if (ksmod(ids).eq.40) then
c                                 si-o mrk fluid
         gph = 0d0 

         call setxyp (ids,id)

         do k = 1, nstot(ids)
            gph = gph + gzero(jnd(k)) * y(k)
         end do 

         gph = gph + gerk (p0a)

      else if (ksmod(ids).ge.29.and.ksmod(ids).le.32) then 
c                                 nastia's models:
          if (ksmod(ids).eq.29) then 
c                                 BCC Fe-Si Lacaze and Sundman
             gph = gfesi(xco(jco(id)+1), gproj (jend(ids,3)), 
     *                                   gproj (jend(ids,4)) )

          else if (ksmod(ids).eq.32) then 
c                                 BCC Fe-Cr Andersson and Sundman
             gph = gfecr1(xco(jco(id)+1), gproj (jend(ids,3)), 
     *                                    gproj (jend(ids,4)) )

          else  
c                                 Nastia's version of BCC/FCC Fe-Si-C Lacaze and Sundman
c                                 this model has to be called ahead of the standard models
c                                 because it sets lrecip(id) = true.
            gph =  gfesic (xco(jco(id)+1),xco(jco(id)+3),xco(jco(id)+4),
     *                     gproj (jend(ids,3)), gproj (jend(ids,4)),
     *                     gproj (jend(ids,5)), gproj (jend(ids,6)),
     *                     ksmod(ids))

          end if 

      else
c                                 normal models (configurational 
c                                 entropy fixed, excess function
c                                 linear in p-t) and special models 
c                                 with normal gmech term
         call setxyp (ids,id)

         if (ksmod(ids).eq.41) then
c                                 ternary coh fluid deltag
            call rkcoh6 (y(2),y(1),gph)

         else if (ksmod(ids).eq.26) then

            call hcneos (gph,y(1),y(2),y(3))

         else 
c                                 get the excess and/or ideal mixing effect
c                                 and/or dqf corrections:
            call gexces (id,gph)

         end if 
c                                 add gmech
         do k = 1, nstot(ids) 
            gph = gph +  gproj (jend(ids,2+k)) * y(k)
         end do 
c                                 for van laar get fancier excess function
         if (llaar(ids)) then

            call setw (ids) 

            call gvlaar (ikp(id),id,gph)

         end if 

      end if 

      gphase = gph

      end

      subroutine factr1 (n,ier)
c-----------------------------------------------------------------------
c factr1 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,n,ip1,istr,ier

      double precision temp,ratio,tmax,rmax

      integer ipvt
      double precision a,d,x(k5)
      common/ cst301 /a(k5,k5),d(k5),ipvt(k5)
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1,n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         x(i) = rmax
      end do 
c                            begin decomposition:
      do i = 1, n - 1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/x(i)
         istr = i
         ip1 = i + 1

         do j = ip1, n
            tmax = dabs(a(j,i))/x(j)
            if (tmax.le.rmax) cycle
            rmax = tmax
            istr = j
         end do 

         if (dabs(rmax).lt.1d-5) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = x(istr)
            x(istr) = x(i)
            x(i) = temp
            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 
         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)
            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do 
 
      end do 
     
      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine uproj
c----------------------------------------------------------------------
c subroutine uproj computes the potentials of saturated phase components
c and saturated components.
  
c the energies of saturated components are projected through
c saturated volatile components.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i,j,k,l,ict,ll,i1,id

      double precision uss(h6), fo2, gph, u, gphase

      external gphase

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      fo2 = 0d0
c                                 compute the chemical potentials of
c                                 saturated phase components.
      if (ifct.gt.0) call ufluid (fo2)

      do i = 1, isat
c                                 determine stable saturated composants
c                                 and the corresponding chemical potentials
         ict = isct(i)

         ll = icp+i

         do j = 1, ict

            k = ids(i,j)
            gph = gphase (k)
            
            if (ifct.gt.0) then 
               do l = 1, 2
c                                 legendre transform for saturated phase
c                                 component potentials
                  if (iff(l).ne.0) gph = gph - cp(iff(l),k)*uf(l)
               end do 
            end if 

            uss(j) = gph 

            if (i.gt.1) then 
c                                 if multiple component saturation constraints
c                                 apply saturation hierarchy legendre transform:
               i1 = i-1
               do l = 1, i1
                  uss(j) = uss(j)-cp(icp+l,k)*us(l)
               end do
            end if 

            g(k) = uss(j)
            uss(j) = uss(j)/cp(ll,k)
         end do 
c                                 if O2, check if fo2 has been 
c                                 determined by a fluid phase routine,
c                                 if so, add the transform:
         if (io2.eq.i) then 
            do j = 1, ict 
               uss(j) = uss(j) + r*t*fo2
            end do 
         end if 
c                           now find stable "composant":

         u = uss(1)

         id = 1

         if (ict.ne.1) then 
            do j = 2, ict
               if (uss(j).gt.u) cycle  
               id = j
               u = uss(j)
            end do 
         end if 
c                               save the id of the stable composant.
         idss(i) = ids(i,id)
c                               and its chemical potential.
         us(i) = u
c                               in case a phase in the component
c                               saturation space is an endmember of
c                               a solution transform the endmember G's:
         do j = 1, ict
            k = ids(i,j)
            g(k) = g(k) - cp(icp+i,k)*u
         end do 

      end do 

      end

      subroutine gall 
c-----------------------------------------------------------------------
c subroutine gall computes molar free energies of all static compounds.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, id

      double precision gval, dg, g0(m4)

      double precision gex, gfesi, gfesic, gerk, gproj, ghybrid, gzero,
     *                 gfecr1, gcpd

      external gerk, gzero, gex, gfesi, gfesic, gproj, ghybrid, 
     *         gcpd

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision g
      common/ cst2 /g(k1)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer isoct
      common/ cst79 /isoct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jspec
      common/ cxt8 /jspec(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)
c-----------------------------------------------------------------------
c                                 compute the chemical potential
c                                 of the projected components.
      call uproj
c                                 first do the endmembers:
      do id = kphct+1, ipoint

         g(id) = gproj (id) 

      end do 
c                                 now do solutions:
      do i = 1, isoct
c                                 check if normal solution:
         if (.not.llaar(i).and.(ksmod(i).eq.7.or.ksmod(i).eq.10.or.
     *       ksmod(i).eq.2.or.ksmod(i).eq.24.or.ksmod(i).eq.25.or.
     *       ksmod(i).eq.28)) then 
c                                 it's normal margules or ideal:
            do j = 1, jend(i,2)
c                                 initialize with excess energy, dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               call setxyp (i,id)

               do k = 1, lstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * y(k)
               end do 

               id = id + 1

            end do 

         else if (ksmod(i).eq.0) then
c                                 it's a fluid compound, the way 
c                                 things are now it must have two
c                                 components.
            do j = 1, lstot(i)
               g0(j) = gzero(jend(i,2+j))
            end do 

            do j = 1, jend(i,2)

               call fexces (id,gval)

               call setxyp (i,id)

               g(id) = g0(1) * y(1) + g0(2) * y(2) + gval

               id = id + 1

            end do 

         else if (lorder(i)) then
c                                 compute margules coefficients
            call setw (i)
c                                 compute enthalpy of ordering
            call oenth (i)
c                                 now for each compound:
            do j = 1, jend(i,2)
c                                 for speciation models gexces
c                                 evaluates only endmember sconf 
c                                 and internal dqf's
               call gexces (id,g(id))

               call setxyp (i,id)

               call specis (dg,i)
c                                 add in g from real endmembers, this
c                                 must include the g for the disordered equivalent
c                                 of the ordered species
               do k = 1, lstot(i)

                  g(id) = g(id) + g(jend(i,2+k)) * p0a(k)

               end do 

               g(id) = g(id) + dg

               id = id + 1

            end do 

         else if (llaar(i)) then 
c                                 compute margules coefficients
            call setw (i)
c                                 because the hp van laar may have p-t
c                                 dependent volumes, the full expression
c                                 must be evaluated here:                   
            do j = 1, jend(i,2)
c                                 initialize with dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               call setxyp (i,id)

               do k = 1, lstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * y(k)
               end do 
c                                 add the real excess energy
               g(id) = g(id) + gex(i,y)

               id = id + 1

            end do 

         else if (ksmod(i).eq.20) then 
c                                 electrolytic solution, assumes:
c                                 1) molal electrolyte standard state
c                                 for solutes.
c                                 2) water is the last species
c                                 solvent species Gibbs energies:
c                                 solvent Gibbs energies
            rt = r*t

            do k = 1, ns
               aqg(k) = g(jnd(k))
            end do 
c                                 compute compound properties
            do j = 1, jend(i,2)
c                                 get the composition
               call setxyp (i,id)
c                                 solvent properties
               call slvnt1 (g(id))
c                                 add in solute properties
               call slvnt2 (g(id))

               id = id + 1

            end do

         else if (ksmod(i).eq.26) then
c                                 H2O-CO2-Salt:
            do j = 1, jend(i,2)

               call setxyp (i,id)
        
               call hcneos (g(id),y(1),y(2),y(3))

               do k = 1, nstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * y(k)
               end do 

               id = id + 1

            end do

         else if (ksmod(i).eq.39) then
c                                 generic hybrid EoS
            do j = 1, jend(i,2)

               g(id) = 0d0
c                                 load composition array and pointers 
               call setxyp (i,id)

               do k = 1, nstot(i)
c                                 sum pure species g's
                  g(id) = g(id) + g(jnd(k)) * y(k)

               end do
c                                 compute and add in activities
               g(id) = g(id) + ghybrid (y)

               id = id + 1

            end do 

         else if (ksmod(i).ge.29.and.ksmod(i).le.32) then 
c                                 nastia's models:
            do j = 1, jend(i,2)
            
               if (ksmod(i).eq.29) then 
c                                 BCC Fe-Si Lacaze and Sundman
                  g(id) = gfesi(xco(jco(id)+1),g(jend(i,3)),
     *                                         g(jend(i,4)))
               else if (ksmod(i).eq.32) then 
c                                 BCC Fe-Cr Andersson and Sundman
                  g(id) = gfecr1 (xco(jco(id)+1),g(jend(i,3)),
     *                                           g(jend(i,4)))
               else 

                  g(id) = gfesic (xco(jco(id)+1),
     *                            xco(jco(id)+3),xco(jco(id)+4),
     *                            g(jend(i,3)),g(jend(i,4)),
     *                            g(jend(i,5)),g(jend(i,6)),ksmod(i))
               end if 
  
               id = id + 1

            end do

         else if (ksmod(i).eq.41) then 

            do j = 1, jend(i,2)
c                                 hybrid MRK ternary COH fluid
               call setxyp (i,id)

               call rkcoh6 (y(2),y(1),g(id)) 

               do k = 1, 3 
                  g(id) = g(id) + g(jnd(k)) * y(k)
               end do 

               id = id + 1

            end do 

         else if (ksmod(i).eq.40) then 

            do j = 1, jend(i,2)
c                                 MRK silicate vapor
               g(id) = 0d0

               call setxyp (i,id)

               do k = 1, lstot(i) 
                  g(id) = g(id) + gzero(jnd(k)) * y(k)
               end do 

               g(id) = g(id) + gerk(y)

               id = id + 1

            end do  

         end if 

      end do 

      end

      subroutine ufluid (fo2)
c----------------------------------------------------------------------
c subroutine ufluid computes the potential of the components
c of a saturated fluid phase. if the mole fraction of a component is les
c less than 1.d-38 the chemical potential is set to -9.9d09.
c ufluid may call one of three molecular fluid equations of state, or
c alternatively users may supply their own routines, however,
c the routines currently in use return the log of a components fugacity
c which is then added to the reference state potential computed by the
c function gphase.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
 
      double precision xf(2), fo2, fs2, gcpd, gzero 

      external gcpd, gzero

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug
c-----------------------------------------------------------------------
c                           compute the chemical potentials of
c                           fluid components in fluid saturated
c                           systems.
      call cfluid (fo2,fs2)

      if (idfl.ne.0) then
   
         uf(idfl) = gcpd (idfl,.false.) + r * t * f(idfl)

      else

         xf(1) = 1d0 - xco2
         xf(2) = xco2
 
         do i = 1, 2

            if (iff(i).ne.0) then 

               if (xf(i).lt.1d-38) then 

                  uf(i) = -1d10

               else 

                  uf(i) = gzero (i) + r * t * f(i)

               end if

            end if 

         end do

      end if 
 
      end


      subroutine outlim 
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer i, j, k, l, ibad1, ibad2, ibad3, igood

      logical bad1, bad2, good, reach

      double precision num

      character char8*8
c                                 -------------------------------------
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c                                 global variables:
c                                 working arrays
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m4)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      logical refine
      common/ cxt26 /refine

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
      
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct
c----------------------------------------------------------------------
      ibad1 = 0 
      ibad2 = 0 
      ibad3 = 0 
      igood = 0 
      rewind (n10)
      if (lopt(11)) rewind (n11)

      if (isoct.eq.0) goto 99

      bad1 = .false.
      bad2 = .false.
      good = .false.

      do i = 1, isoct

         if (.not.stable(i)) then
            bad1 = .true.
            ibad1 = ibad1 + 1
         else
            good = .true.
            igood = igood + 1
         end if
 
         if (limit(i)) then 
            bad2 = .true.
            if (relax(i)) then 
               ibad2 = ibad2 + 1
            else 
               ibad3 = ibad3 + 1
            end if 
         end if 

      end do 

      if (.not.refine) write (n10,*) ibad1,0,igood
c                                 write solutions present that are 
c                                 not stable
      if (bad1) then 

         write (*,1000)
         if (lopt(11)) write (n11,1000)
     
         do i = 1, isoct
            if (.not.stable(i)) then 
               write (*,'(5x,a)') fname(i) 
               if (.not.refine) write (n10,'(a)') fname(i)
               if (lopt(11)) write (n11,'(5x,a)') fname(i) 
            end if 
         end do
      end if 

      if (.not.good) goto 99
c                                 write solutions that are on an internal
c                                 limit
      if (bad2.and.icopt.gt.3) then 
c                                 adaptive minimization
         if (ibad2.gt.0) then 
c                                 solutions whose limits could be relaxed
            write (*,1080) 
            if (lopt(11)) write (n11,1080) 
            do i = 1, isoct
               if (limit(i).and.relax(i)) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

         if (ibad3.gt.0) then 
c                                 solutions whose limits could NOT be relaxed
            write (*,1090) 
            if (lopt(11)) write (n11,1090) 
            do i = 1, isoct
               if (limit(i).and.(.not.relax(i))) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

      else if (bad2) then 
c                                 non-adaptive minimization,
c                                 solutions on internal limits 
         write (*,1010) 
         if (lopt(11)) write (n11,1010) 
         do i = 1, isoct
            if (limit(i)) then
               write (*,'(5x,a)') fname(i) 
               if (lopt(11)) write (n11,'(5x,a)') fname(i)
            end if  
         end do

      end if 

      reach = .false.

      do i = 1, isoct

         if (.not.stable(i)) cycle

         if (int(reachg(i)*2d0/nopt(21)-1d0).gt.0) reach = .true.

         if (.not.refine) then

            write (n10,'(a)') fname(i)

            do j = 1, ostg(i)
               do k = 1, ndim(j,i) 
                  write (n10,*) xlo(k,j,i),xhi(k,j,i)
               end do
            end do 

         end if 
c                                 special case (1 component solution). 
         if (ksmod(i).eq.39.and.ndim(1,i).eq.0) cycle

         if (istg(i).eq.1) then 
c                                 single site solution
            write (*,1020) fname(i)
            if (lopt(11)) write (n11,1020) fname(i)

            if (ksmod(i).eq.20) then 
c                                 charge balance model:
               do j = 1, ndim(1,i)

                  if (j.lt.ns) then
                     char8 = names(jnd(j))
                  else if (j.eq.ns) then
                     cycle 
                  else 
                     char8 = aqnam(jnd(j)  - aqst)
                  end if 

                  write (*,1030) char8, xlo(j,1,i), xhi(j,1,i)

                  if (lopt(11)) write (n11,1030) char8, xlo(j,1,i),
     *                                                  xhi(j,1,i)

               end do     
        
            else 

               do j = 1, ndim(1,i)

                  write (*,1030) names(jend(i,2+j)),
     *                              xlo(j,1,i),xhi(j,1,i)


                  if (lopt(11)) write (n11,1030) 
     *                          names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
     
               end do

            end if

         else 

            if (ostg(i).eq.istg(i)) then
c                                 prismatic 
               write (*,1040) 'prismatic model: '//fname(i)
               if (lopt(11)) write (n11,1040) 
     *                        'prismatic model: '//fname(i)

            else 
c                                 prismatic + orphan vertices
               write (*,1040) 'prism + orphan model: '//fname(i)
               if (lopt(11)) write (n11,1040)
     *                        'prism + orphan model: '//fname(i)

            end if 

            do j = 1, istg(i)
               
               write (*,1050) j
               if (lopt(11)) write (n11,1050) j

               if (ispg(i,j).eq.1) then 

                  write (*,1060)
                  if (lopt(11)) write (n11,1060) 

               else

                  do k = 1, ispg(i,j) - 1

                     write (*,1070) k,xlo(k,j,i),xhi(k,j,i)
                     if (lopt(11)) write (n11,1070) 
     *                              k,xlo(k,j,i),xhi(k,j,i)

                  end do 

               end if 

            end do

            if (ostg(i).gt.istg(i)) then
c                                 prismatic + orphan vertices
               l = ostg(i)

               write (*,1150)
               if (lopt(11)) write (n11,1150)

               do j = 1, ndim(l,i)

                  k = 2 + lstot(i) - ndim(l,i) + j

                  write (*,1030) names(jend(i,k)),xlo(j,l,i),xhi(j,l,i)
                  if (lopt(11)) write (n11,1030) 
     *                           names(jend(i,k)),xlo(j,l,i),xhi(j,l,i)

               end do

            end if

         end if

      end do

      if (reach) then 

         write (*,1100)

         do i = 1, isoct
            if (int(reachg(i)*2d0/nopt(21)-1d0).eq.0) cycle
            write (*,1110) fname(i), int(reachg(i)*2d0/nopt(21)-1d0)
         end do 

      end if 

99    if (goodc(1)+badc(1).gt.0d0) then
         num = badc(1)/(badc(1)+goodc(1))*1d2
         write (*,1120) num, badc(1) + goodc(1)
         if (num.gt.1d-1) call warn (53,num,i,'OUTLIM')
         write (*,1140) dfloat(jcount(1))/(badc(1)+goodc(1))
         goodc(1) = 0d0
         badc(1) = 0d0 
      end if 

      jcount(3) = 0
      jcount(4) = 0
      jcount(1) = 0

      close (n10)
      if (lopt(11)) close (n11)

1000  format (/,'**warning ver992** The following solutions were input,'
     *         ,' but are not stable:',/)
1010  format (/,'**warning ver991** The following solutions have ',
     *          'compositions at',
     *          ' an internal limit (i.e., 0<x<1)',/,'(see ranges ',
     *          'below to determine which limits should be relaxed or',
     *        /,'if executing in auto_refine mode inrease auto_refine',
     *          '_slop in perplex_option.dat):',/)
1020  format (/,'Endmember fractions for model: ',a,//,5x,
     *          'Endmember     Minimum         Maximum')
1030  format (5x,a8,4x,g12.5,4x,g12.5)
1040  format (/,'Compositions for ',a)
1050  format (/,'  Simplex ',i1,/,5x,
     *               'Composition   Minimum         Maximum')
1060  format (8x,'Dummy site generated by model reformulation',/)
1070  format (8x,i1,8x,g12.5,4x,g12.5)
1080  format (/,'**warning ver993** The compositions of the following',
     *        ' solutions reached internal',/,
     *        'limits that were automatically relaxed:',/)
1090  format (/,'**warning ver991** The compositions of the following',
     *        'solutions reached internal',/,'limits that could not be',
     *        ' automatically relaxed. To avoid this problem change',/,
     *        'the subdivision mode for ',
     *        'the solutions or increase auto_refine_slop.',/)
1100  format (/,'The following solution models have non-zero reach_',
     *          'increment:',//,t30,'reach_increment')
1110  format (4x,a,t35,i2)
1120  format (/,'The failure rate during speciation (order-disorder) ',
     *        'calculations is ',f7.3,'%',/,'out of a total of ',f12.0,
     *        ' calculations.',/)
1130  format (/,'Metastability filtering eliminated ',f7.3,'% of the ',
     *        'trial compositions prior',/,'to optimization.',/)
1140  format (/,'Average number of iterations per speciation ',
     *          'calculation:',f5.1,/)
1150  format (/,2x,'Orphans',/,5x,
     *          'Endmember     Minimum         Maximum')

      end 

      subroutine subst (a,ipvt,n,b,ier)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision a(k8,k8),b(k8),x(k8),sum

      integer ipvt(k8),ip,i,j,n,ii,ier
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)
      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j)+sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous,
c                                 but reopt requires it. should check
c                                 what's with factor. 
         ier = 1
         goto 99
      end if 

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do 

         if (a(i,i).eq.0d0) then
c                                 as above.
            ier = 1
            goto 99
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 
      b(n) = x(n)
 
99    end

      logical function solvus (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two static pseudocompounds of
c solution ids.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp 

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      solvus = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle 

         if (dabs(cp(i,id1)-cp(i,id2))/dcp(i,ids).gt.soltol) then 
            solvus = .true.
            exit
         end if 

      end do 

      end 

      subroutine aqlagd (id,bad,recalc)
c-----------------------------------------------------------------------
c given chemical potentials solve for rock dominated aqueous speciation
c configured to be called from resub with output to the (molar normalized)
c arrays g2/cp2, if recalc is true then aqlagd is being used to recover 
c speciation and arrays g2/cp2 are not loaded.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, id, badct

      logical bad, recalc, lmus, feos

      double precision mo(l9), blk(k5), gamm0, totm, g0(l9), lmu(k8),
     *                 tmu(k8),is, gso(nsp), lnkw, gtot, smo, err,
     *                 slvmo(nsp), solmol, negox, posox

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

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

      double precision yf,g,v
      common/ cstcoh /yf(nsp),g(nsp),v(nsp)

      double precision gh,dvhy,g0mrk
      common/ csthyb /gh(nsp),dvhy(nsp),g0mrk(nsp) 

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

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision sel, cox
      logical hscon, hsc, oxchg
      common/ cxt45 /sel(k0),cox(k0),hscon,oxchg,hsc(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision p0, dz
      common/ cxt46 /p0, dz

      save badct, lmu, lmus
      data badct/0/
c----------------------------------------------------------------------
      if ((.not.mus .and..not.recalc).or.
     *    (.not.lmus.and.recalc)) then 

         lmus = .false.
         bad = .true.
         return

      else 

         if (recalc) then 
c                                 use lagged chemical potentials
            do i = 1, icp 
               tmu(i) = lmu(i)
            end do 
c                                 set flag so slvnt3 evaluates pure 
c                                 fluid eos (not clear if this is necessary).
            feos = .true.

         else 

            lmus = .true.

            do i = 1, icp

               lmu(i) = mu(i)
               tmu(i) = mu(i)

               if (cblk(i).eq.0d0.and..not.lopt(36)) then 
c                                 check that the solvent does not contain
c                                 the absent component
                  do j = 1, ns 

                     if (y(j).gt.0d0.and.cp(i,jnd(j)).gt.0d0) then

                        bad = .true.
                        return

                     end if 

                  end do

               end if

               feos = .false.

            end do

         end if 

         call slvnt3 (gso,.false.,feos,id)

         if (epsln.lt.nopt(34)) then
c                                 eos is predicting vapor phase 
c                                 solvent densities
            bad = .true.
            return
         end if 

         bad = .false.

      end if
c                                 iterate on speciation
      call aqsolv (g0,gso,mo,tmu,is,gamm0,lnkw,bad)

      if (bad) then 

         badct = badct + 1

         if (badct.lt.11) 
     *        call warn (99,0d0,0,'AQLAGD did not converge on solute '//
     *                       'speciation')

         if (badct.eq.10) call warn (49,0d0,99,'AQLAGD')
 
         return

      end if
c                                 back calculated bulk composition
      do j = 1, icp
         blk(j) = 0d0
      end do  

      smo = 0d0 
      gtot = 0d0
      err = 0d0 
c                                 everything on a molal basis
c                                 first the solutes
      do i = 1, aqct

         if (mo(i).eq.0d0) cycle
c                                 charge balance error 
         err = err + q(i) * mo(i)
c                                 total g
         gtot = gtot + mo(i) * (g0(i) + rt*dlog(mo(i)*gamm0**q2(i)))
c                                 total molality
         smo = smo + mo(i)
c                                 accumulate component moles
         do j = 1, icp
            blk(j) = blk(j) + mo(i)*aqcp(j,i)
         end do 

      end do

      solmol = smo
c                                 for the solvent mole fractions
c                                 need to accumulate total
c                                 molality first
      do i = 1, ns 
c                                 solvent molality:
         slvmo(i) = yf(ins(i))/msol
c                                 total molality
         smo = smo + slvmo(i)
c                                 moles/kg-solvent 
         do j = 1, icp 
            blk(j) = blk(j) + slvmo(i)*cp(j,jnd(i))
         end do

      end do

      do i = 1, ns
c                                 solvent bulk mole fraction:
         caq(id,i) = slvmo(i)/smo
         if (caq(id,i).eq.0d0) cycle
          gtot = gtot + slvmo(i) * (gso(i) + rt*dlog(caq(id,i)))

      end do 
c                                 bulk fluid composition 
      totm = 0d0
      err = dabs(err)*1d1

      if (lopt(36).and.oxchg) then 
c                                check on charge imbalance
         posox = 0d0
         negox = 0d0 

         do j = 1, icp
            if (cox(j).gt.0) then
               posox = posox + cox(j)*blk(j)
            else
               i = j
               negox = negox + cox(j)*blk(j)
            end if
         end do

         blk(i) = blk(i) - (posox+negox)/cox(i)

c         write (*,*) 'err, posox+negox0',err,posox+negox

         posox = 0d0
         negox = 0d0 

         do j = 1, icp
            if (cox(j).gt.0) then
               posox = posox + cox(j)*blk(j)
            else
               i = j
               negox = negox + cox(j)*blk(j)
            end if
         end do

c         write (*,*) 'err, posox+negox1',err,posox+negox

      end if 

      do j = 1, icp
c                                zero bulk compositions below chg balance error
         if (blk(j).lt.err) blk(j) = 0d0 
c                                totm is the total number of moles of the
c                                components in a solution of smo moles of 
c                                species
         totm = totm + blk(j)
      end do

      if (recalc) then
c                                 stuff needed for output:
         do i = 1, aqct 
            caq(id,ns+i) = mo(i)
c                                 output bad result point to n13
            if (mo(i).ge.nopt(35).and.q(i).eq.0d0) then

               write (n13,'(i3,1x,4(g14.6,1x),a,9(1x,f9.5))')
     *                                i, p0, dz, t, p, aqnam(i),
     *                                epsln,(yf(ins(j)),j=1, ns)
               write (*,'(i3,1x,4(g14.6,1x),a,9(1x,f9.5))') 
     *                                i, p0, dz, t, p, aqnam(i),
     *                                epsln,(yf(ins(j)),j=1, ns)

            end if

         end do 
c                                 ionic strength
         caq(id,na1) = is
c                                 total molality
         caq(id,na2) = smo
c                                 solvent mass
         caq(id,na3) = msol
c                                 error in log10(Kw)
         caq(id,na3+1) = (-lnkw + dlog(mo(ihy)*mo(ioh)*gamm0**2
     *                    ))/2.302585d0
c                                  pH
         caq(id,na3+2) = -dlog10(mo(ihy)*gamm0)
c                                  Delta_pH
         caq(id,na3+3) = caq(id,na3+2) + lnkw/4.605170d0
c                                  solute molality
         caq(id,na3+4) = solmol
c                                  dielectric cst
         caq(id,nat) = epsln

      else 
c                                 stuff need for optimization:
c                                 load into molar normalized arrays 
c                                 used by resub, g per mole of components
         g2(jphct) = gtot/totm

         do j = 1, icp
c                                 bulk composition per mole of components
            cp2(j,jphct) = blk(j)/totm
         end do
c                                c2tot is the number of moles of the 
c                                components in a solution with 1 mole of
c                                species, this is needed for consistent
c                                ouput (i.e., a mol of the phase is per
c                                mol of species rather than per mole of 
c                                components). at the cost of k21 real vars
         c2tot(jphct) = totm/smo

      end if 

      end

      subroutine slvnt3 (gso,whysp,feos,id)
c-----------------------------------------------------------------------
c for back and lagged speciation calculations get solvent properties.
c   if whysp -> take the fluid fractions from ysp(:,id) and renormalize.
c   if feos -> call fluid eos (i.e., pure g's haven't been calculated).
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical whysp, feos

      integer i, k, id

      double precision gso(nsp), dum, gcpd, ysum

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

      double precision yf,gmrk,v
      common/ cstcoh /yf(nsp),gmrk(nsp),v(nsp)

      double precision g
      common/ cst2 /g(k1)

      double precision gh,dvhy,gmrk0
      common/ csthyb /gh(nsp),dvhy(nsp),gmrk0(nsp)

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer spct
      double precision ysp
      character spnams*8
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
      rt  = r*t

      if (ns.gt.1) then
c                                 a multi species solvent is present: 
         if (whysp) then 

            ysum = 0d0

            do i = 1, ns 
               y(i) = ysp(i,id)
               ysum = ysum + y(i)
            end do
c                                 renormalize
            do i = 1, ns 
               y(i) = y(i)/ysum
            end do

         end if

         if (feos) then

            do k = 1, ns
c                                 call to gcpd will put the pure species molar
c                                 volumes arrays vhyb0 and vmrk0 (cxt38) and 
c                                 the mrk pure fugacity coeff in gp.
               aqg(k) = gcpd(jnd(k),.false.)

            end do

         else 

            do k = 1, ns
c                                 use previously computed g
               aqg(k) = g(jnd(k))

            end do

         end if
c                                  dum is just a dummy.
         call slvnt1 (dum)

         do k = 1, ns
c                                 add in the solvent species activity coefficients
c                                 under the assumption that these are independent of the 
c                                 solute speciation the partial g of a solvent species will be 
c                                 g(i) = gs0(i) + RT ln x(i).
            i = ins(k)

            gso(k) = aqg(k) + rt * dlog(gmrk(i)/gmrk0(i))

         end do

      else
c                                  solvent is pure water
         y(1) = 1d0
         ysp(1,id) = 1d0

         call slvnt0 (gso(1),dum)

      end if

      end 

      subroutine aqsolv (g0,gso,mo,mu,is,gamm0,lnkw,bad)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, it, jt, iexp

      logical bad, kill

      double precision c(l9), mo(*), mu(*), dg, xis, dn, xdix, 
     *                 d(l9), is, gamm0, g0(*), lnkw, dix, 
     *                 gso(*), xdn, qb(l9)

      double precision gcpd, solve, aqact

      external gcpd, solve, aqact

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)
c----------------------------------------------------------------------
      if (epsln.lt.nopt(34)) then
c                                  vapor, same as checking lnkw
         bad = .true.
         return
      end if
c                                  set default charge balance ion (aq_ion_H+, lopt(44)
      if (lopt(44)) then 
         ion = ihy
      else 
         ion = ioh
      end if 
c                                  if default choice fails switch to back-up choice
      do k = 1, 2
c                                 set up coefficients for mo(ion) equation
         g0(ion) = gcpd(aqst+ion,.false.) 
c                                 compute solute properties 
         do i = 1, aqct 
c                                 dg is the solvent oxide potentials - g
            g0(i) = gcpd(aqst + i, .false.)
            qr(i) = q(i)/q(ion)
            qb(i) = (q(ion)-q(i))*q(i)
            dg = -g0(i) + qr(i)*g0(ion)

            kill = .false.

            do j = 1, jbulk 

               dn = aqcp(j,i) - qr(i)*aqcp(j,ion)

               if (dn.eq.0d0) cycle

               if (isnan(mu(j))) then 
c                                 kill any species that depends on
c                                 an undetermined chemical potential
c                                 unless NOT lopt(36), this allows
c                                 oxide components without redox because
c                                 g(H+) is not a function of mu(O2) [but
c                                 is not fail-safe].
                  if (aqcp(j,i).ne.0d0.and..not.lopt(36)) then 
                     kill = .true.
                     exit
                  else
                     cycle
                  end if

               else if (cblk(j).eq.0d0.and..not.lopt(36)) then
c                                 this check is necessary because lp may
c                                 give a zero-amount solution for the chemical
c                                 potential of an absent component. the test
c                                 cannot be made with oxide components.
                  if (aqcp(j,i).ne.0d0) then 
                     kill = .true.
                     exit
                  end if

               end if

               dg = dg + dn * mu(j)

            end do 
c                                 normalize by RT
            if (kill) then 
               dg = 0d0 
            else 
               dg = dexp(dg/rt)
            end if 

            if (q(i).ne.0d0) then 
c                                  this is now c(i)*a(ion)^(q(i)) = mo(i)*gamma(i)*q(i)
c                                  the rhs q(i) is because the eq to be solved will be
c                                  sum (q(i)*m(i)) = 0. 
               d(i) = q(i)*dg
               c(i) = d(i)

            else 
c                                  neutral species assumed to be ideal, molality is
               if (dg.gt.nopt(35)) then 
                  dg = nopt(35)
               end if 

               mo(i) = dg

            end if 

         end do 
c                                  initialize iteration loop
         lnkw = (gso(ns) - g0(ioh))/rt

         if (c(ioh).eq.0d0) then
c                                 no hydrogen or no oxygen
            bad = .true.
            return

         end if

         mo(ion) = dexp(lnkw/2d0)
         gamm0 = 1d0 
         is = 0d0 

         xdn = 1d0
         iexp = 1
         it = 0
         jt = 0
         xdix = 1d99
         bad = .false.
c                                  iteration loop for ionic strength
         do 
c                                  solve charge balance for ion
            mo(ion) = solve(c,qr,mo(ion),jchg,ichg,bad)

            if (bad) then
               xis = is
               exit
            end if
c                                  back calculate charged species molalities
c                                  and ionic strength
            xis = is
            is = 0d0

            do i = 1, ichg 
               j = jchg(i)
               mo(j) = c(j) / q(j) * mo(ion)**qr(j) 
               is = is + q2(j) * mo(j)
            end do

            is = is / 2d0 

            dn = is - xis

            if (dabs(dn).gt.1d0/2d0**iexp) then 
               dn = dn/dabs(dn)/2d0**iexp 
               if (dn*xdn.lt.0d0) iexp = iexp + 1
               is = xis + dn
            end if 
c                                 DH law activity coefficient factor (gamma = gamm0^(q^2))
            gamm0 = aqact(is)
c                                 check for convergence
            dix = dabs(xis-is)/is

            if (dix.lt.nopt(5)) then
c                                 converged
               return

            else if (it.gt.iopt(21)) then

               if (xdix.gt.dix.and.jt.lt.10) then
c                                 try again?
                  it = 0
                  jt = jt + 1
                  xdix = dix

               else 
c                                 diverging
                  bad = .true.
                  exit

               end if

            end if

            xdn = dn 

            it = it + 1
c                                 update coefficients
            do i = 1, ichg 
               j = jchg(i)
               c(j) = d(j)*gamm0**qb(j)
            end do

         end do
c                                 switch to the backup ion
         if (lopt(44)) then
            ion = ioh
         else 
            ion = ihy
         end if

      end do

      end


      double precision function aqact (is)
c-----------------------------------------------------------------------
c compute Debye-Hueckel-type activity coefficient factor
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision is

      integer ihy, ioh
      double precision gf, epsln, epsln0, adh, msol
      common/ cxt37 /gf, epsln, epsln0, adh, msol, ihy, ioh
c----------------------------------------------------------------------
      aqact = dexp(adh*dsqrt(is)/(1d0 + dsqrt(is)) + 0.2d0*is)

      end 

      subroutine tabhed (n,vmn,dv,nv,nvar,n5name,n6name)
c----------------------------------------------------------------------
c  write header for nvar-dimensional table output
c     vmn - minimum values of the indendpendent variables
c      dv - the increments
c      nv - number of nodes
c     lop - werami option
c    nvar - number of independent variables
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, nv(2), nvar, ivar, n

      double precision vmn(2), dv(2)

      character*100 n6name, n5name, vname(l3)*14

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      character vnm*8
      common/ cxt18a /vnm(l3)

      integer iam
      common/ cst4 /iam

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio
c------------------------------------------------------------------------
c                                 generate a file name and
c                                 open the file on n5
      if (iam.eq.1) then 
c                                 frac2d file names: on input n5name
c                                 contains the string to be appended 
c                                 to the project name
         call fopenv (n,n5name)
      else 
         call fopenn (n,nvar,n5name,n6name)
      end if 
c                                 set flag for cumulative properties
c                                 generated by allprp/allmod
      first = .true.
c                                 ctr file tab format header
      write (n,'(a)') '|6.6.6'
c                                 a title
      write (n,'(a)') n5name
c                                 number of indepedent variables
      write (n,*) nvar
c                                 independent variable names, 
c                                 value, increment & nodes
      do i = 1, nvar
         write (n,'(a)') vnm(i)
         write (n,*) vmn(i)
         write (n,*) dv(i)
         write (n,*) nv(i)
      end do 
c                                 number of pseudo-dependent variables,
      ivar = 2
      if (icopt.eq.7.and.fileio) then
         ivar = 3
      else if (icopt.eq.9.and.iam.eq.1) then 
         ivar = 1
      else 
         ivar = 2
      end if 
c                                 convert a8 names to a14
      do i = 1, ivar
         vname(i) = vnm(i)
         call unblnk (vname(i))
      end do 
c                                 output the dependent variable counter and list
      if (kcx(1).eq.999) then 
c                                  phemgp file
         write (n,*) ivar + iprop + 2
         write (n,'(200(a20,1x))') 'Name','Counter',
     *                             (vname(i), i = 1, ivar),
     *                             (dname(i),i = 1, iprop)

      else 
c                                  tab file
         if (lopt(15).or.nvar.eq.1) then

            if (.not.lopt(29)) then
c                                  with pseudo-dependent variables 
               write (n,*) ivar + iprop 
               write (n,'(200(a14,1x))') (vname(i), i = 1,ivar),
     *                                   (dname(i), i = 1,iprop)
            else 
               write (n,*) ivar + iprop + 5

               write (n,'(200(a14,1x))') (vname(i), i = 1,ivar),
     *                                   (dname(i), i = 1,iprop),
     *         'y_SiO2','y_SiO','y_O','y_O2','y_Si'
            end if 

         else 
c                                  terse format
            write (n,*) iprop 
            write (n,'(200(a14,1x))') (dname(i), i = 1,iprop)
          
         end if 

      end if

      end 

      subroutine fopenn (n,dim,n5name,n6name)
c----------------------------------------------------------------------
c decide on file name and type for werami output files, open n5name on 
c LUN n, n6name is opened if necessary by prphed or modhed
c dim - integer 1 = 1d, 2 = 2d, 0 = text echo
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, dim, n

      character*100 n5name, n6name, num*4

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl, first
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),
     *               kop(i11),kcx(i11),k2c(i11),iprop,
     *               first,kfl(i11),tname

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c----------------------------------------------------------------------
c                                 make plot file
      do i = 1, 1000
c                                 loop to find an unused name made of 
c                                 project + "i"
         write (num,'(a1,i3)') '_',i

         call unblnk (num) 

         call mertxt (tfname,prject,num,0)
 
         if ((kop(1).eq.36.or.kop(1).eq.38).and.kcx(1).eq.999) then 
c                                 phemgp format
            call mertxt (n5name,tfname,'.phm',0)

         else

            if (dim.eq.0) then 
               call mertxt (n5name,tfname,'.txt',0)
            else
               call mertxt (n5name,tfname,'.tab',0)
            end if 

            if (kop(1).eq.25) call mertxt (n6name,tfname,'.plt',0) 

         end if 
           
         open (n, file=n5name, status='new', iostat=ier)
c                                 presume error means a file with name
c                                 n5name already exists
         if (ier.eq.0) exit 
 
         if (i.gt.999) call error (999,0d0,i,tfname)

      end do

      if (dim.eq.0) write (*,1000) n5name

1000  format (/,'Console output will be echoed in file: ',a,/)

      end

      subroutine fopenv (n,string)
c----------------------------------------------------------------------
c decide on file name and type for werami output files, open n5name on 
c LUN n, n6name is opened if necessary by prphed or modhed
c dim - integer 1 = 1d, 2 = 2d, 0 = text echo
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier, n

      character string*(*)

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c----------------------------------------------------------------------

      call mertxt (tfname,prject,string,0)

      call mertxt (string,tfname,'.tab',0)

      open (n, file=string, status='replace', iostat=ier)

      if (ier.ne.0) call error (999,0d0,i,
     *             'file '//tfname//' is in use by another application')

      end

      subroutine dumper (iclos,id,hkp,jkp,amt,lambda) 
c----------------------------------------------------------------------
c dump phase data from yclos routines:
c     iclos = 0 => yclos0 etc
c     hkp refinement point
c     jkp solution model
c     amt - amount
c     lambda - lambda
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, iclos, id, hkp, jkp

      double precision amt, lambda

      character name*14

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      call getnam (name,jkp)

      if (iclos.eq.1) then
c                                 static
         write (*,1000) id,hkp,jkp,name,amt,lambda,c(id),
     *                  (a(i,id),i=1,jbulk)
      else
c                                 dynamic 
         write (*,1000) id,hkp,jkp,name,amt,lambda,g2(id),
     *                  (cp2(i,id),i=1,jbulk)
      end if

1000  format (i6,1x,i3,1x,i4,1x,a,20(g14.6,1x))

      end 

      double precision function ginc (dt,dp,id)
c-----------------------------------------------------------------------
      implicit none

      double precision dt,dp,gee,gsol,gfrnd

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 put NaN check to prevent a NaN increment 
c                                 from setting p/t to a NaN, this should only 
c                                 happen if a phase has no volumetric and/or 
c                                 caloric EoS. JADC 9/18/2016.
      if (isnan(dp)) dp = 0d0
      if (isnan(dt)) dt = 0d0

      p = p + dp 
      t = t + dt 

      if (iam.eq.5) then 
c                                 frendly 
         gee = gfrnd(-id)

      else 
c                                 meemum/werami
         gee = gsol(id)

      end if 

      p = p - dp 
      t = t - dt

      ginc = gee 

      end

      double precision function gsol (id)
c-----------------------------------------------------------------------
c gsol computes the total (excess+ideal) free energy of solution 
c for a solution identified by index ids and composition y(m4) input
c from cxt7, the composition y is the independent endmember fractions
c for all model types except reciprocal solutions, in which case it is 
c the y's for the full reciprocal model.

c gsol assumes the endmember g's have not been calculated by gall and is
c      called by WERAMI and MEEMUM.
c gsol1 is identical to gsol but can only been called after gall and is 
c      called by VERTEX and MEEMUM. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, k, id

      double precision dg, g, gso(nsp), gamm0

      double precision omega, hpmelt, gmelt, gfluid, gzero, aqact,
     *                 gex, slvmlt, gfesi, gcpd, gerk, gfecr1, ghybrid

      external gphase, omega, hpmelt, gmelt, gfluid, gzero, gex, slvmlt,
     *         gfesi, gerk, gfecr1, ghybrid, gcpd, aqact

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jnd
      double precision aqg,q2s,rt
      common/ cxt2 /aqg(m4),q2s(m4),rt,jnd(m4)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer jspec
      common/ cxt8 /jspec(h9,m4)

      double precision yf,gh,v
      common/ cstcoh /yf(nsp),gh(nsp),v(nsp)

      integer jd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,jd

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer ion, ichg, jchg
      double precision q, q2, qr
      common/ cstaq /q(l9),q2(l9),qr(l9),jchg(l9),ichg,ion

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c----------------------------------------------------------------------
      if (id.lt.0) then 

         gsol = gcpd (-id,.true.)

      else 

         g = 0d0
c                                 initialize p-t dependent coefficients
         call ingsol (id)

         if (ksmod(id).eq.2) then 
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
            call gdqf (id,g,y) 
c                                 add entropy and excess contributions
            g = g - t * omega (id,y) + gex (id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id) 
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

         else if (lrecip(id).and.lorder(id)) then 
c                                 -------------------------------------
c                                 convert y coordinates to independent p coordinates
            call y2p0 (id)
c                                 get the speciation, excess and entropy effects.
            call specis (g,id)

            do k = 1, lstot(id) 
c                                 compute mechanical g from these z's, 
c                                 specip adds a correction for the ordered species.
               g = g + gcpd (jend(id,2+k),.true.) * p0a(k)
            end do 
c                                 get the dqf, this assumes the independent reactants
c                                 are not dqf'd. gex not neccessary as computed in specip
            call gdqf (id,g,p0a)

         else if (lorder(id)) then 
c                                 -------------------------------------
c                                 non-reciprocal speciation.
            do k = 1, lstot(id)  
               pa(k) = y(k)
               p0a(k) = y(k)
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

            call specis (dg,id)

            g = g + dg 
c                                 get dqf corrections
            call gdqf (id,g,p0a)

         else if (lrecip(id)) then 
c                                 -------------------------------------
c                                 macroscopic reciprocal solution w/o order-disorder

c                                 convert y's to p's (p0a here).
            call y2p0 (id)

            do k = 1, lstot(id)
               g = g + gcpd (jend(id,2+k),.true.) * p0a(k) 
            end do 
c                                 get the dqf
            call gdqf (id,g,p0a)
c                                 and excess contributions
            g = g - t * omega (id,p0a) + gex (id,p0a)

         else if (ksmod(id).eq.20) then 
c                                 electrolytic solution
c                                 -------------------------------------
c                                 compute solvent mass and gibbs energy:
            rt = r*t

            do k = 1, ns
c                                 solvent species gibbs energy and volumes
               if (y(k).le.0d0) cycle
               aqg(k) = gcpd(jnd(k),.true.)

            end do 
c                                 solvent properties
            call slvnt1 (g)
c                                 add in solute properties
            call slvnt2 (g)

         else if (ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp melt model         
            call gdqf (id,g,y) 

            g = g - t * hpmelt (id,y) + gex (id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

         else if (ksmod(id).eq.25) then 
c                                 -------------------------------------
c                                 ghiorso pmelt model  
            call gdqf (id,g,y) 

            g = g - t * gmelt (id) + gex (id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

         else if (ksmod(id).eq.26) then 
c                                 ------------------------------------
c                                 andreas salt model
            call hcneos (g,y(1),y(2),y(3))

            do k = 1, 3
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

         else if (ksmod(id).eq.28) then 
c                                 -------------------------------------
c                                 high T fo-fa-sio2 model  
            call gdqf (id,g,y) 

            g = g - t * slvmlt() + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gcpd (jend(id,2+k),.true.)
            end do 

         else if (ksmod(id).eq.29) then 
c                                 -------------------------------------
c                                 BCC Fe-Si Lacaze and Sundman
            g = gfesi(y(1), gcpd (jend(id,3),.true.), 
     *                      gcpd (jend(id,4),.true.) )

         else if (ksmod(id).eq.32) then 
c                                 -------------------------------------
c                                 BCC Fe-Cr Andersson and Sundman
            g =  gfecr1(y(1), gcpd (jend(id,3),.true.), 
     *                        gcpd (jend(id,4),.true.) )

         else if (ksmod(id).eq.39) then
c                                 -------------------------------------
c                                 generic hybrid EoS
            if (lopt(32).and.caq(jd,na1).ne.0d0) then 
c                                 lagged speciation
               call slvnt3 (gso,.false.,.true.,id)
c                                 DH law activity coefficient factor (ln[g] = lng0*q^2)
               gamm0 = aqact(caq(jd,na1))
c                                 solvent species (caq => mole fraction)
               do k = 1, ns
                  if (caq(jd,k).eq.0d0) cycle
                  g = g + caq(jd,k) * (gso(k) + rt * dlog(caq(jd,k)))
               end do 
c                                 solute species (caq => molality)
               do k = sn1, nsa 

                  i = k - ns

                  if (caq(jd,k).eq.0d0) cycle

                  g = g + caq(jd,k)/caq(jd,na2) * (gcpd(aqst+i,.false.) 
     *                  + rt*dlog(caq(jd,k)*gamm0**q2(i)))

               end do 

            else 

               do k = 1, mstot(id)
c                                 sum pure species g's
                  g = g + gcpd(jnd(k),.true.) * y(k)

               end do
c                                 compute and add in activities
               g = g + ghybrid (y)

            end if

         else if (ksmod(id).eq.41) then 
c                                 hybrid MRK ternary COH fluid
            call rkcoh6 (y(2),y(1),g) 

            do k = 1, nstot(id) 
               g = g + gcpd(jnd(k),.true.) * y(k)
            end do 

         else if (ksmod(id).eq.40) then 
c                                 MRK silicate vapor
            do k = 1, nstot(id) 
               g = g + gzero(jnd(k)) * y(k)
            end do 

            g = g + gerk(y)

         else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos. hardwired to special
c                                 component choices
            do k = 1, 2

               g = g + gzero (jnd(k))*y(k)

            end do 
c                                 don't know whether it's a speciation routine
c                                 so set the fluid species fractions just in case,
c                                 this is only necessay for species output by 
c                                 WERAMI/MEEMUM, these will be reset if it actually
c                                 is a speciation routine.
            yf(2) = y(jspec(id,1))
            yf(1) = 1d0 - yf(2)
c
            g = g + gfluid (yf(2))

         else 

            write (*,*) 'what the **** am i doing here?'
            stop

         end if 

         gsol = g 

      end if 

      end


      double precision function gfrnd (id)
c-----------------------------------------------------------------------
c function to get g's for frendly. differs from gphase in that it checks
c for special components O2, H2O, CO2. sloppy but who cares?
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision gee, fo2, fs2, gcpd
 
      external gcpd

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer eos
      common/ cst303 /eos(k10)

      integer ifct,idfl
      common/ cst208 /ifct,idfl
c-----------------------------------------------------------------------

      gee = gcpd (id,.false.) + r * t * dlog(act(id))

      if (ifct.gt.0.and.eos(id).lt.100) then 
c                                 this is a quick fix that will
c                                 call the fluid routine way more 
c                                 than necessary.
         call cfluid (fo2,fs2)

         if (id.eq.idf(3)) then 

            gee = gee + r*t*fo2

         else if (id.eq.idf(1)) then 
        
            gee = gee + r*t*fh2o
         
         else if (id.eq.idf(2)) then 
        
            gee = gee + r*t*fco2

         end if
 
      end if 

      gfrnd = gee

      end 

      subroutine setstc (ids,i,j) 
c----------------------------------------------------------------------
c set stretching transformation for the site i, species j, of solution ids,
c it seems doubtful this is worth the effort, functions strtch and unstch
c could just use the stch array directly.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids,i,j

      double precision stch
      common/ cst47 /stch(h9,mst,msp,4)

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm
c----------------------------------------------------------------------
      bm1  = stch(ids,i,j,1)
      bp1  = stch(ids,i,j,2)
      bpm  = stch(ids,i,j,3)
      lbpm = stch(ids,i,j,4)

      end

      double precision function stinc (x,dy,ids,i,j) 
c----------------------------------------------------------------------
c stinc increments cartesian x by the conformal increment dy for solution 
c ids, site i, species j. this is used to set limits for the conformal
c y in terms of the cartesian x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids,i,j

      double precision x, dy, strtch, unstch

      external strtch, unstch

      double precision stch
      common/ cst47 /stch(h9,mst,msp,4)

      double precision bp1,bm1,bpm,lbpm
      common/ cst46 /bp1,bm1,bpm,lbpm
c----------------------------------------------------------------------
c                                 set stretching parameters
      call setstc (ids,i,j) 
c                                 unstretch x and decrement/increment it
c                                 by +/- one conformal increament, then 
c                                 restretch.
      stinc = strtch ( unstch (x) + dy )

      end

      subroutine prs2xy (i,ids,dynam,bad)
c----------------------------------------------------------------------
c convert the raw compositional coorinates stored for the ith entry
c of prism to the prismatic compositional array x and convert the x
c array to y.

c if dynam then store the prismatic coordinates in the the zcoor array/
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer i, j, k, l, m, ids, kcoct

      logical bad, dynam

      double precision ysum

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)
c----------------------------------------------------------------------
      bad = .false.

      if (dynam) then

         jcoor(jphct) = jcoct - 1
         kcoct = jcoct + mcoor(ids)
c                                 counter for number of non 0 or 1 compositions
         if (kcoct.gt.k20) call error (59,x(1,1),k20,'resub')

      end if 

      l = (i-1)*mcoor(ids)
      m = 0
       
      if (ksmod(ids).ne.20) then

         do j = 1, ostg(ids)

            ysum = 0d0

            do k = 1, ndim(j,ids)

               m = m + 1

               x(j,k) = prism(l+m)
               ysum = ysum + x(j,k)

               if (dynam) then 

                  zcoor(jcoct) = x(j,k)
 
                  if (x(j,k).lt.xmno(ids,j,k).and.
     *                x(j,k).gt.xmxo(ids,j,k)) then 
c                                 the composition is out of range
                     jphct = jphct - 1
                     jcoct = kcoct - mcoor(ids)
                     bad = .true.
                     return

                  end if

                  jcoct = jcoct + 1

               end if 

            end do

            x(j,k) = 1d0 - ysum

         end do

      else 
c                                 charge balance models: a wierd shuffle to put
c                                 the first nqs - 1 species in zcoor
         ysum = 0d0

         do k = 1, nqs

            if (k.eq.ns) cycle 

            m = m + 1

            x(1,k) = prism(l+m)
            ysum = ysum + x(1,k)

            if (k.eq.nqs) exit

            if (dynam) then 

               zcoor(kcoct-nqs+k) = x(1,k)

               if (x(1,k).lt.xmno(ids,1,k).and.
     *             x(1,k).gt.xmxo(ids,1,k)) then 
c                                 the composition is out of range
                  jphct = jphct - 1
                  jcoct = kcoct - mcoor(ids)

                  bad = .true.
                  return

               end if

            end if

         end do

         x(1,ns) = 1d0 - ysum

         if (dynam) then 
            zcoor(kcoct-qn) = x(1,ns)
            jcoct = jcoct + nqs1
         end if

      end if

      call xtoy (ids,ids,.true.,bad)

      if (bad.and.dynam) then 

         jphct = jphct - 1
         jcoct = kcoct - mcoor(ids)

      end if

      end

      subroutine xtoy (ids,id,usex,bad)
c----------------------------------------------------------------------
c subroutine to convert prismatic solution compositions (x(i,j))
c to geometric endmember fractions (y) for solution model ids.

c usex - use x coordinates, otherwise use x3 coordinates (only done by getloc
c        when called by meemum).
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids, id, k, l, m

      logical bad, zap, zbad, usex
 
      external zbad

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical badend, sck, nrf
      integer ldsol
      common/ cxt36 /ldsol(m4,h9),badend(m4,h9),sck(h9),nrf(h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c----------------------------------------------------------------------

      bad = .false.
      zap = bad
      k = 0

      if (usex) then

         do l = 1, mstot(ids)
            y(l) = 1d0
            do m = 1, istg(ids)
               y(l) = y(l)*x(m,kmsol(ids,l,m))
            end do

            if (badend(l,ids).and.y(l).gt.zero) zap = .true.
            if (y(l).gt.one) k = l

         end do

      else

         do l = 1, mstot(ids)
            y(l) = 1d0
            do m = 1, istg(ids)
               y(l) = y(l)*x3(id,m,kmsol(ids,l,m))
            end do
         end do

      end if

      if (ksmod(ids).eq.9.or.ksmod(ids).eq.10) then
c                                 necessary? probably none of this is. 
         do l = mstot(ids) + 1, pstot(ids)
            y(l) = 0d0
         end do

         if (.not.usex) then
c                                 load the orphan fractions into local x
c                                 for y2p0
            do l = 1, ispg(ids,ostg(ids))
               x(ostg(ids),l) = x3(id,ostg(ids),l)
            end do 

         end if

      end if
c                                 checks for x array compositions:
      if (k.ne.0) then 
c                                 reject pure independent endmember compositions. 
         if (ldsol(k,ids).gt.0.and.nrf(ids)) then 

            bad = .true.

            return

         end if 

         y(k) = 1d0 

         do l = 1, mstot(ids)
            
            if (l.eq.k) cycle
            
            y(l) = 0d0

         end do

      end if
c                                 invalid dependent endmember
      if (lopt(43).and.zap) then 
c                                 convert y's to p's
         call y2p0 (ids)
c                                 check for bad z's
         if (zbad(pa,ids)) bad = .true.

      end if

      end 

      subroutine oddprs (ids,resub)
c---------------------------------------------------------------------
c subdivision for a regular prism + orphan vertices
c---------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical resub

      integer i, j, h, ids, k, l, m, n, np0, np1, np2, nst1, nst2, i1,
     *        nco

      double precision sum, ycum

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c---------------------------------------------------------------------
c                                 get the npair simplicial coordinates for the 
c                                 orphan site:
      ycum = 0d0

      call chopit (ycum,1d0,0,ndim(ostg(ids),ids),ostg(ids),ids,0,
     *             .true.)

      np0 = npairs
c                                 save an extra coordinate for prismatic vertex (sum)
      nco = ispg(ids,ostg(ids))
      nst1 = nco * np0
      ntot = 0
c                                 for each simplicial coordinate (skipping the 
c                                 origin) make prismatic compositions, adjust 
c                                 the subdivision resolution of the prism so 
c                                 that the renormalized resolution remains constant
      do i = 1, np0
c                                 starting point in simp
         n = (i-1)*nco
c                                 sum is the factor by which the resolution
c                                 on the prismatic site is changed, i.e., 
c                                 inc = inc0/sum [sum is the fraction of the
c                                 prism vertex]
         sum = 1d0

         do j = 1, ndim(ostg(ids),ids)
            sum = sum - simp(n+j)
         end do

         simp(n+nco) = sum

         if (sum.ne.0d0) then 
c                                 the prism vertex fraction is > 0
            ycum = 0d0
c                                 first simplex of prism
            call chopit (ycum,sum,0,ndim(1,ids),1,ids,nst1,.false.)

            np1 = npairs
            nst2 = nst1 + np1*ndim(1,ids)

            ycum = 0d0
c                                 second simplex of prism
            call chopit (ycum,sum,0,ndim(2,ids),2,ids,nst2,.false.)

            np2 = npairs
c                                 ready to rock
            do j = 1, np1 
c                                 site 1 coordinate starting point in simp
               m = nst1 + (j-1)*ndim(1,ids)
c                                  for each population on site 1
               do k = 1, np2
c                                 site 2 coordinate starting point in simp
                  l = nst2 + (k-1)*ndim(2,ids)
c                                  count the coordinate
                  ntot = ntot + 1
c                                 starting point of the coordinate
                  h = (ntot-1)*mcoor(ids)

                  if (h+mcoor(ids).gt.k24) call errk24 (resub)
c                                 load the coordinate
                  do i1 = 1, ndim(1,ids)
                     prism(h+i1) = simp(m+i1)
                  end do

                  do i1 = 1, ndim(2,ids)
                     prism(h+ndim(1,ids)+i1) = simp(l+i1)
                  end do

                  do i1 = 1, ndim(3,ids)
                     prism(h+nsum(ids)+i1) = simp(n+i1)
                  end do 

               end do
            end do

         else
c                                 the prism vertex is zero, load only the
c                                 fractions of the orphans, the other 
c                                 fractions will not be used.

c                                 count the coordinate
            ntot = ntot + 1
c                                 starting point of the coordinate
            h = (ntot-1)*mcoor(ids)

            do i1 = 1, ndim(3,ids)
               prism(h+nsum(ids)+i1) = simp(n+i1)
            end do 

         end if

      end do 

      end