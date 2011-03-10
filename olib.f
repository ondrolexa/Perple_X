c--------------olib.f---------------------------------------------------
c output routines called only by werami/meemum
c-----------------------------------------------------------------------
      subroutine calpr0 (lu)
c----------------------------------------------------------------------
c calpr0 - output properties of an assemblage, can be called by either
c meemum or werami. if meemum, prints chemical potentials.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      character cprop*18

      integer i,j,l,lu

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5) 

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      double precision atwt
      common/ cst45 /atwt(k0)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision mu
      common/ cst330 /mu(k8)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer jtest,jpot
      common/ debug /jtest,jpot

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      character pname*14
      common/ cxt21a /pname(k5)

      integer iam
      common/ cst4 /iam
c---------------------------------------------------------------------- 
                                    
      write (lu,1000)  

      if (iam.eq.2) then 
         write (lu,1120) (vname(jv(i)),v(jv(i)), i = 1, ipot)
         write (lu,1120) (vname(jv(i)),v(jv(i)), i = 3, ipot)
      else 
         write (lu,1120) (vnm(i), var(i), i = 1, jvar)
      end if 

      if (iopt(2).eq.0) then 
         cprop = 'molar  proportions'
      else
         cprop = 'weight percentages'
      end if

      write (lu,1020) cprop, (cname(i), i = 1, icomp)

      do i = 1, ntot

         write (lu,1030) pname(i), 
c                                 weight %
     *                   props(17,i)*props(16,i)/psys(17)*1d2,
c                                 vol %
     *                   props(1,i)*props(16,i)/psys(1)*1d2,
c                                 mol %
     *                   props(16,i)/psys(16)*1d2,
c                                 mol
     *                   props(16,i),
c                                 molar or weight composition
     *                   (pcomp(l,i), l = 1, icomp)
      end do 

      write (lu,1160)
c                                 phase/system summary, normal thermo:
      do i = 1, ntot
c                                 N, H, V, Cp, alpha, beta, density
         write (lu,1170) pname(i),props(17,i),props(2,i),props(15,i),
     *                   props(1,i),(props(j,i),j=12,14),props(10,i)
      end do

      write (lu,1170) 'System        ',psys(17),psys(2),psys(15),
     *                psys(1),(psys(j),j=12,14),psys(10)
      if (aflu) write (lu,1170) 'System - fluid',psys1(17),psys1(2),
     *                psys1(15),psys1(1),(psys1(j),j=12,14),psys1(10)

      write (lu,1190)
c                                 phase/system summary, seismic:
      do i = 1, ntot
         write (lu,1200) pname(i), (props(j,i), j = 3, 8) 
      end do

      write (lu,1200) 'System        ',(psys(j), j = 3, 8) 
      if (aflu) write (lu,1200) 'System - fluid',(psys1(j), j = 3, 8) 

      write (lu,1240)
c                                 phase/system summary, seismic derivatives:
      do i = 1, ntot
         write (lu,1250) pname(i),props(18,i),props(20,i),props(19,i),
     *                  props(21,i),props(22,i),props(25,i),
     *                  props(23,i),props(26,i),props(24,i),props(27,i)
      end do

      write (lu,1250) 'System        ',psys(18),psys(20),psys(19),
     *                psys(21),psys(22),psys(25),psys(23),psys(26),
     *                psys(24),psys(27)
      if (aflu) write (lu,1250) 'System - fluid',psys1(18),psys1(20),
     *                psys1(19),psys1(21),psys1(22),psys1(25),psys1(23),
     *                psys1(26),psys1(24),psys1(27)

      if (.not.aflu.or.(aflu.and.psys1(1).eq.0)) then 
c                                 no fluid is present, or the system consists
c                                 entirely of fluid (psys1(1)=0):
         write (lu,1210)

         write (lu,1040)

         do i = 1, icomp

            write (lu,1110) cname(i),fbulk(i), fbulk(i)/gtot*1d2,
     *                      fbulk(i)*atwt(i)/psys(17)*1d2
         end do

         write (lu,1220)

         write (lu,1060) 
c                                 enthalpy, specific enthalpy
     *                   psys(2)/psys(1)*1d5/psys(10),
     *                   psys(2)/psys(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys(15)/psys(1)*1d5/psys(10),
     *                   psys(15)/psys(1)*1d5,
c                                 cp, specific cp 
     *                   psys(12)/psys(1)*1d5/psys(10),
     *                   psys(12)/psys(1)*1d5

      else 
c                                 fluid is present
         write (lu,1210)

         write (lu,1080)

         do i = 1, icomp

            write (lu,1110) cname(i),fbulk(i),fbulk(i)/gtot*1d2,
     *               fbulk(i)*atwt(i)/psys(17)*1d2,fbulk1(i)/gtot1*1d2,
     *               fbulk1(i)*atwt(i)/psys1(17)*1d2          
         end do

         write (lu,1220)
c                                 true bulk properties:
         write (lu,1060) 
c                                 enthalpy, specific enthalpy
     *                   psys(2)/psys(1)*1d5/psys(10),
     *                   psys(2)/psys(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys(15)/psys(1)*1d5/psys(10),
     *                   psys(15)/psys(1)*1d5,
c                                 cp, specific cp 
     *                   psys(12)/psys(1)*1d5/psys(10),
     *                   psys(12)/psys(1)*1d5 

c                                 solid only bulk properties:
         write (lu,1100) 
c                                 enthalpy, specific enthalpy
     *                   psys1(2)/psys1(1)*1d5/psys1(10),
     *                   psys1(2)/psys1(1)*1d5, 
c                                 entropy, specific entropy 
     *                   psys1(15)/psys1(1)*1d5/psys1(10),
     *                   psys1(15)/psys1(1)*1d5,
c                                 cp, specific cp 
     *                   psys1(12)/psys1(1)*1d5/psys(10),
     *                   psys1(12)/psys1(1)*1d5

         write (lu,1230) 

      end if 
c                                 chemical potentials variance
      if (jpot.ne.1) then 
         write (lu,1130) (cname(i), i = 1, hcp)
         write (lu,1140) (mu(i), i = 1, hcp)
         write (lu,1071) 2, jbulk - ntot + 2 
      else 
         write (lu,1070) 2, jbulk - ntot + 2 
      end if 

1000  format (/,40('-'),//,'Stable phases at:')
1020  format (/,'Phase Compositions (',a,'):',
     *        /,19x,'wt %',6x,'vol %',5x,'mol %',5x,'mol  ',
     *          5x,20(1x,a,3x))
1030  format (1x,a,3x,3(f6.2,4x),g9.3,1x,20(f7.3,2x))
1040  format (/,14x,'mol',7x,'mol %',6x,'wt %')
1060  format (/,' Enthalpy (J/kg) = ',g12.6,/,
     *          ' Specific Enthalpy (J/m3) = ',g12.6,/,
     *          ' Entropy (J/K/kg) = ',g12.6,/,
     *          ' Specific Entropy (J/K/m3) = ',g12.6,/,
     *          ' Heat Capacity (J/K/kg) = ',g12.6,/,
     *          ' Specific Heat Capacity (J/K/m3) = ',g12.6,/)
1070  format ('Variance (c-p+',i1,') = ',i2,//,40('-'),/)
1071  format (/,'Variance (c-p+',i1,') = ',i2,//,40('-'),/)
1080  format (/,16x,'Complete Assemblage',15x,'Solid+Melt Only',
     *        /,14x,'mol',7x,' mol %',6x,'wt %',9x,' mol %',6x,'wt %')
1100  format (/,' Solid Enthalpy (J/kg) = ',g12.6,/,
     *          ' Solid Secific Enthalpy (J/m3) (2) = ',g12.6,/,
     *          ' Solid Entropy (J/K/kg) = ',g12.6,/,
     *          ' Solid Specific Entropy (J/K/m3) = ',g12.6,/,
     *          ' Solid Heat Capacity (J/K/kg) (1) = ',g12.6,/,
     *          ' Solid Specific Heat Capacity (J/K/m3) (1) = ',g12.6,/)
1110  format (1x,a8,2x,f8.3,5x,2(f6.2,4x),5x,2(f6.2,4x))
1120  format (29x,a8,' = ',g12.6)
1130  format (/,'Chemical Potentials (J/mol):',/,2x,20(4x,a,5x))
1140  format (2x,20(1x,g13.6))
1160  format (/,'Molar Properties and Density:'
     *        /,20x,'N(g)',8x,'H(J)',6x,'S(J/K)',6x,'V(J/bar)',6x,
     *         'Cp(J/K)'
     *         ,6x,'Alpha(1/K)',2x,'Beta(1/bar)',2x,'Density(kg/m3)')
1170  format (1x,a,1x,f9.2,3x,13(g12.5,1x))
1190  format (/,'Seismic Properties:'
     *        /,17x,'Gruneisen',7x,'Ks(bar)',7x,'Mu(bar)',
     *        4x,'V0(km/s)',5x,'Vp(km/s)',5x,'Vs(km/s)')
1200  format (1x,a,3x,12(g12.5,1x))
1210  format (/,'Bulk Composition:')
1220  format (/,'Other Bulk Properties:')
1230  format (/,'N.B.: Aggregate properties represent the entire stable'
     *         ,' assemblage.',/,'Solid aggregate properties represent '
     *         ,'solid and melt properties,',/,'but do not include '
     *         ,'molecular fluid properties.',/)
1240  format (/,'Isochemical Seismic Derivatives:',
     *        /,16x,'Ks_T(bar/K)',2x,'Ks_P',2x,'Mu_T(bar/K)',
     *           2x,'Mu_P',2x,'Vphi_T(km/s/K)',1x,'Vphi_P(km/s/bar)',1x,
     *          'Vp_T(km/s/K)',2x,'Vp_P(km/s/bar)',2x,
     *          'Vs_T(km/s/K)',2x,'Vs_P(km/s/bar)')
1250  format (1x,a,3x,2(f8.2,1x,f8.4,2x),12(g12.5,3x))

      end 

      subroutine getloc (itri,jtri,ijpt,wt,nodata)
c-----------------------------------------------------------------------
c getloc computes local properties requested by either meemum or werami

c if called by werami and ijpt > 1 properties are computed as a 
c weighted mixture (wt(i)) of the assemblages present at nodes (itri(i),
c jtri(i), i = 1, ijpt). otherwise (ijpt=1) the assemblage is that 
c present at node itri(1), jtri(1).

c    ncpdg  -> ncpd
c    npg    -> np
c    idbulk -> kkp
c    xcoor  -> x3

c also sets flags (could set a solvus flag):

c    aflu  -> fluid present
c    fluid -> the phase is a fluid (indexed)

c see getphp got contents of props/psys/psys1 arrays
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,m,ids,jds,jd,kd,jcoor,kcoor,itri(4),jtri(4),ijpt

      double precision wt(3), cst

      logical sick(i8), nodata, ssick, ppois
c                                 x-coordinates for the assemblage solutions
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 
      integer ifp
      common/ cxt32 /ifp(k1)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 global assemblage data
      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer igrd
      common/ cst311/igrd(l7,l7)

      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision bg
      common/ cxt19 /bg(k5,k2)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision mus
      common/ cst48 /mus(k8,k2)

      double precision mu
      common/ cst330 /mu(k8)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------
c                                 logarithmic_p option
      if (lopt(14)) p = 1d1**p 

      if (iam.ne.2) then 

         jd = igrd(itri(1),jtri(1))
         ias = iap(jd)
c                                 no data test
         if (ias.eq.k3) then 
            nodata = .true.
            goto 99 
         end if 

         np = iavar(1,ias)
         ncpd = iavar(2,ias)
         ntot = iavar(3,ias)

         do i = 1, ntot
            kkp(i) = idasls(i,ias)
         end do 

         jcoor = icoor(jd)
c                                 get the dependent potentials
         if (jpot.ne.1) then
 
            do i = 1, hcp
               mu(i) = 0d0
            end do 

            do i = 1, ijpt

               kd = igrd(itri(i),jtri(i))

               do j = 1, hcp
                  mu(j) = mu(j) + wt(i) * mus(j,kd)
               end do 

            end do 

         end if
c                                 if s/v independent variables, 
c                                 then set p and t 
         if (hcp.gt.icp) then

            if (lopt(14)) then 
               write (*,*) 'ERROR: logarithmic_p must be false for USV'
               stop
            end if 

            p = -mu(pindex)
            t = mu(tindex)

            if (p.lt.0d0.or.t.lt.1d2) then 
               nodata = .true.
               goto 99
            end if 

         end if 

      end if 

      aflu = .false.
      shear = .true.
      volume = .true.
      nodata = .false.
      ssick = .false.
      ppois = .false.
c                                 flag for bulk bad bulk properties

c                                 initialize bulk properites
c                                 total mass
      gtot = 0d0
      gtot1 = 0d0

      do i = 1, i8
         psys(i) = 0d0
         psys1(i) = 0d0
         pgeo(i) = 0d0 
         pgeo1(i) = 0d0
      end do 

      do i = 1, icomp
c                                 total molar amounts
         fbulk(i) = 0d0
         fbulk1(i) = 0d0

      end do

      do i = 1, ntot

         ids = kkp(i)

         if (i.le.np) then 

            if (ksmod(ids).eq.0.or.ksmod(ids).gt.20.and.lopt(6)) then 
               aflu = .true.
               fluid(i) = .true.
            else 
               fluid(i) = .false.
            end if

            if (iam.ne.2) then
 
               do j = 1, istg(ids)
                  do k = 1, ispg(ids,j)
                     jcoor = jcoor + 1
                     cst = bg(i,jd)
c                                 in case zero mode is not on, allow
c                                 composition of zero phase
                     if (ijpt.eq.1.and.cst.eq.0d0) cst = 1d0
                     x3(i,j,k) = wt(1)*cst*xcoor(jcoor)

                  end do 
               end do 

            end if 

         else 

            if (ifp(-ids).eq.1) then 
               aflu = .true.
               fluid(i) = .true.
            else 
               fluid(i) = .false.
            end if

         end if
c                                 molar amounts
         if (iam.eq.2) then 

            props(16,i) = amt(i)
c                                 convert x3 to y for calls to gsol            
            if (ids.gt.0) call x3toy (i,ids)

         else 
c                                 if werami with interpolation, average
c                                 compositions:
            props(16,i) = wt(1)*bg(i,jd)
c                                 now average in other assemblages
            do j = 2, ijpt
            
               kd = igrd(itri(j),jtri(j))
               ias = iap(kd)

               kcoor = icoor(kd)

               do k = 1, i

                  jds = idasls(k,ias)

                  if (k.le.np) then 

                     do l = 1, istg(jds)
                        do m = 1, ispg(jds,l)
                           kcoor = kcoor + 1
                           if (k.lt.i) cycle
c                                 this is done so as to count
c                                 the coordinates as well as make
c                                 the composition.   
                           x3(i,l,m) = x3(i,l,m) 
     *                               + wt(j)*bg(k,kd)*xcoor(kcoor)

                        end do 
                     end do 

                  end if 

                  if (k.lt.i) cycle 
c                                 molar amounts
                  props(16,i) = props(16,i) + wt(j)*bg(k,kd)

               end do  

            end do
c                                 renormalize the composition
            if (i.le.np) then 

               do l = 1, istg(ids)
                  do m = 1, ispg(ids,l)    
                     cst = props(16,i)
                     if (cst.eq.0d0) cst = 1d0    
                     x3(i,l,m) = x3(i,l,m)/cst
                  end do 
               end do 

            end if
             
            if (ids.gt.0) then 
c                                 revover x from x3, 2nd arg has no meaning.
               call getxz (i,i,ids)
c                                 convert x to y for calls to gsol
               call xtoy (ids)

            end if 

         end if 
c                                 get and sum phase properties
         call getphp (ids,i,sick,ssick,ppois)     

      end do 
c                                 compute aggregate properties:
      call gtsysp (sick,ssick)

99    if (lopt(14)) p = dlog10(p)

      end

      subroutine x3toy (id,ids)
c----------------------------------------------------------------------
c subroutine to convert geometric reciprocal solution compositions (x3(id,i,j))
c to geometric endmember fractions (y) for solution model ids. replicate 
c of subroutine xtoy, but for the x3 array (used only by getloc from meemum).
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer ids, l, m, ld, id
c                                 -------------------------------------
c                                 global variables:
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------

      do l = 1, mstot(ids)
c                                 the endmembers may have been
c                                 rearranged from the original order,
c                                 use knsp(l,ids) to assure correct
c                                 indexing
         ld = knsp(l,ids) 

         y(ld) = 1d0

         do m = 1, istg(ids)
            y(ld) = y(ld)*x3(id,m,kmsol(ids,ld,m))
         end do

      end do   

      end 

      subroutine gtsysp (sick,ssick)
c-----------------------------------------------------------------------
c computes aggregate (system) properties from sums accumulated by
c prior calls to getphp
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical sick(i8), ssick

      integer i, iwarn

      double precision chi, chi1, units, root, r43

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      character pname*14
      common/ cxt21a /pname(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      save iwarn
      data iwarn /0/
c----------------------------------------------------------------------
c                                 correct for proportional wt 
c                                 on intensive properties (alpha, beta). 
      do i = 3, 21

         if (i.gt.5.and.i.lt.13.or.i.gt.14.and.i.lt.18) cycle 
         psys(i) = psys(i)/psys(1)
         pgeo(i) = pgeo(i)/psys(1)
         if (psys1(1).ne.0d0) then  
            psys1(i) = psys1(i)/psys1(1)
            pgeo1(i) = pgeo1(i)/psys1(1)
         end if 

      end do 
c                                 weighting scheme for seismic velocity
c                                 chi = 1 -> voigt 0 -> reuss 0.5 -> VRH
      chi = nopt(6)
      chi1 = 1d0 - chi
      units = dsqrt(1d5)/1d3
      r43   = 4d0/3d0
c                                 aggregate properties
      if (psys1(1).ne.0d0) psys1(10) = psys1(17)/psys1(1)*1d2
c                                 density, kg/m3
      psys(10) = psys(17)/psys(1)*1d2

      if (volume.and..not.rxn) then 
c                                 gruneisen
         psys(3) = psys(3) + chi1/pgeo(3)
c                                 bulk modulus
         psys(4) = psys(4) + chi1/pgeo(4) 
c                                 bulk modulus T-deriv
         psys(18) = psys(18) + chi1/pgeo(18) 
c                                 bulk modulus P-deriv
         psys(20) = psys(20) + chi1/pgeo(20) 

         root = psys(4)/psys(10)

         if (root.gt.0d0) then 
c                                 sound velocity
            psys(6) = dsqrt(root) * units
c                                 sound velocity T derivative
            psys(22) = (psys(18) + psys(4) * psys(13)) 
     *                  / dsqrt(root) / psys(10) / 2d0 * units
c                                 sound velocity P derivative
            psys(25) = (psys(20) - psys(4) * psys(14)) 
     *                  / dsqrt(root) / psys(10) / 2d0 * units
         end if 

      end if 

      if (volume.and.shear.and..not.rxn) then 
c                                 aggregate seismic properties, if a 
c                                 fluid is present the reuss mean is 
c                                 is infinite, signaled by pgeo = 0.
         if (pgeo(5).gt.0d0) then 
c                                 shear modulus
            psys(5) = psys(5) + chi1/pgeo(5)
         else 
c                                 fluid present
            psys(5) = psys(5)
         end if 

         if (pgeo(19).ne.0d0) then 
c                                 shear modulus P-derivative
            psys(19) = psys(19) + chi1/pgeo(19)
         else 
c                                 fluid present
            psys(19) = psys(19)
         end if

         if (pgeo(21).ne.0d0) then 
c                                 shear modulus T-derivative
            psys(21) = psys(21) + chi1/pgeo(21)
         else 
c                                 fluid present
            psys(21) = psys(21)
         end if

         root = (psys(4)+r43*psys(5))/psys(10)

         if (root.gt.0d0) then 
c                                 p-wave velocity
            psys(7) = dsqrt(root)*units
c                                 p-wave velocity T derivative
            psys(23) = (psys(18) + r43*(psys(19) + psys(13) * psys(5)) 
     *                 + psys(4) * psys(13)) / 
     *                 dsqrt(root) / psys(10) / 2d0 * units
c                                 p-wave velocity P derivative
            psys(26) = (psys(20) + r43*(psys(21) - psys(14) * psys(5)) 
     *                 - psys(4) * psys(14)) /
     *                 dsqrt(root) / psys(10) / 2d0 * units
         end if 

         root = psys(5)/psys(10)

         if (root.gt.0d0) then 
c                                 s-wave velocity
            psys(8) = dsqrt(root) * units
c                                 T-derivative
            psys(24) = (psys(19) + psys(5) * psys(13)) 
     *                 / dsqrt(root) / psys(10) / 2d0 * units
c                                 P-derivative 
            psys(27) = (psys(21) - psys(5) * psys(14)) 
     *                 / dsqrt(root) / psys(10) / 2d0 * units

         end if 
c                                 vp/vs
         if (psys(8).gt.0d0) then 
            psys(9) = psys(7)/psys(8)
         else
            psys(9) = nopt(7)
         end if 

      else 

         do i = 5, 9
            psys(i) = nopt(7)
         end do 

      end if 
c                                 the psys1(1) condition is for the 
c                                 special case of a system consisting 
c                                 only of fluid. 
      if (aflu.and.(.not.ssick).and..not.rxn.and.psys1(1).gt.0d0) then 
c                                 fluid absent properties:
c                                 gruneisen T
         psys1(3) = psys1(3) + chi1/pgeo1(3)
c                                 adiabatic bulk modulus
         psys1(4) = psys1(4) + chi1/pgeo1(4) 
c
         root = psys1(4)/psys1(10)

         if (root.gt.0d0) then 
c                                 sound velocity
            psys1(6) = dsqrt(root) * units
c                                 sound velocity T derivative
            psys1(22) = (psys1(18) + psys1(4) * psys1(13)) 
     *                  / dsqrt(root) / psys(10) / 2d0 * units
c                                 sound velocity P derivative
            psys1(25) = (psys1(20) - psys1(4) * psys1(14)) 
     *                  / dsqrt(root) / psys1(10) / 2d0 * units
         end if 

         if (shear.and.(pgeo1(5).gt.0d0.or.chi.gt.0d0)) then 

            if (pgeo1(5).gt.0d0) then 
c                                 shear modulus
               psys1(5) = psys1(5) + chi1/pgeo1(5)
            else 
c                                 fluid present, use arithmetic mean
               psys1(5) = psys1(5)/chi
            end if 

            if (pgeo1(19).ne.0d0) then 
c                                 shear modulus P-derivative
               psys1(19) = psys1(19) + chi1/pgeo1(19)
            else 
c                                 fluid present, use arithmetic mean
               psys1(19) = psys1(19)/chi
            end if

            if (pgeo1(21).ne.0d0) then 
c                                 shear modulus T-derivative
               psys1(21) = psys1(21) + chi1/pgeo1(21)
            else 
c                                 fluid present, use arithmetic mean
               psys1(21) = psys1(21)/chi
            end if

            root = (psys1(4)+r43*psys1(5))/psys1(10)

            if (root.gt.0d0) then 
c                                 p-wave velocity
               psys1(7) = dsqrt(root)*units
c                                 p-wave velocity T derivative
               psys1(23) = (psys1(18) + r43*(psys1(19) 
     *                    + psys1(13) * psys1(5)) 
     *                    + psys1(4) * psys1(13)) / 
     *                    dsqrt(root) / psys1(10) / 2d0 * units
c                                 p-wave velocity P derivative
               psys1(26) = (psys1(20) + r43*(psys1(21)
     *                    - psys1(14) * psys1(5)) 
     *                    - psys1(4) * psys1(14)) /
     *                    dsqrt(root) / psys1(10) / 2d0 * units
            end if 

            root = psys1(5)/psys1(10)

            if (root.gt.0d0) then 
c                                 s-wave velocity
               psys1(8) = dsqrt(root) * units
c                                 T-derivative
               psys1(24) = (psys1(19) + psys1(5) * psys1(13)) 
     *                    / dsqrt(root) / psys1(10) / 2d0 * units
c                                 P-derivative 
               psys1(27) = (psys1(21) - psys1(5) * psys1(14)) 
     *                    / dsqrt(root) / psys1(10) / 2d0 * units
            end if 

c                                 vp/vs
            if (psys1(8).gt.0d0) then 
               psys1(9) = psys1(7)/psys1(8)
            else
               psys1(9) = nopt(7)
            end if 

         else 

            do i = 5, 9
               psys1(i) = nopt(7)
            end do 

         end if 
         
      end if 

      if ((.not.volume.or..not.shear).and.(iwarn.lt.11)) then

         iwarn = iwarn + 1

         if (.not.shear.and.volume) then
             write (*,1000) t,p
          else if (.not.volume) then 
            do i = 1, ntot
               if (sick(i)) write (*,1010) t,p,pname(i)
            end do 
         end if 

         if (iwarn.eq.11) call warn (49,r,177,'GETSYP') 
                  
      end if  

1000  format (/,'**warning ver177** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'aggregate seismic properties',/,'cannot be computed ',
     *        'because of a missing/invalid shear modulus.',/)
1010  format (/,'**warning ver177** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'aggregate seismic properties ',/,'cannot be computed ',
     *        'because of missing/invalid properties for: ',a,/)

      end 

      subroutine getphp (id,jd,sick,ssick,ppois)
c-----------------------------------------------------------------------
c gets properties of phase id and saves them in props(1:i8,i); 
c if called by werami/meemum id is a general phase pointer; if 
c called by frendly id is an endmember pointer. 

c the properties are saved prop as follows

c 1  - molar volume
c 2  - molar enthalpy
c 3  - gruneisen thermal parm
c 4  - K_S
c 5  - Mu_S
c 6  - v_phi
c 7  - v_p
c 8  - v_s
c 9  - v_p/v_s
c 10 - rho
c 11 - G
c 12 - cp
c 13 - alpha
c 14 - beta
c 15 - S
c 16 - molar amount
c 17 - molar weight

c 18 - KS_T
c 19 - MuS_T
c 20 - KS_P
c 21 - MuS_P

c 22 - vphi_T
c 23 - vp_T
c 24 - vs_T
c 25 - vphi_P
c 26 - vs_P
c 27 - vp_P

c getphp computes isostatic props of the phase identified by id as 
c computed by centered finite differences from the Gibbs energy
c as stored in props(i8,jd)

c the difference increments are

c dt0, dp0 for 1st order derivatives (entropy,volume and enthalpy)
c dt1, dp1 for 2nd order derivatives (heat capacity, expansivity*, 
c          compressibility*)
c dt2, dp2 for 3rd order derivatives (gptt, gppt, gppp, gttt)
c          used for the T derivative of the bulk modulus. 

c *expansivity (alpha) as returned here is 1/v*dv/dt
c *compressibility (beta) as returned here is -1/v*dv/dp

c corrected to check for negative pressure, in which case 
c forward differences are used. june 22, 2004, JADC.
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      logical ok, sick(i8), ssick, pois, ppois

      integer id,jd,iwarn1,iwarn2,j,itemp

      double precision dt0,dt1,dt2,g0,g1a, g2a, dg, ss,alpha1,alpha2,
     *                 dp0,dp1,dp2,e,alpha,v,ginc,beta,cp,s,rho,gtt,r43,
     *                 g1,g2,g3,g4,g5,g7,gppp,gppt,gptt,gttt,mols,units,
     *                 root,chi,chi1

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      character pname*14
      common/ cxt21a /pname(k5)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer iam
      common/ cst4 /iam

      double precision atwt
      common/ cst45 /atwt(k0)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct

      save dt0,dt1,dt2
      data dt0,dt1,dt2/0.5d0,5d0,5d1/

      save iwarn1, iwarn2
      data iwarn1, iwarn2 /2*0/
c----------------------------------------------------------------------
      sick(jd) = .false.
      pois = .false.
c                                 make name and composition, 
c                                 redundant for frendly
      call getnam (pname(jd),id)
c                                 composition, don't call if meemum
      if (iam.ne.2) call getcmp (jd,jd,id)
c                                 component counter for frendly is different
c                                 than for all other programs
      if (iam.ne.5) then 
         itemp = icomp
      else
         itemp = k0
      end if 
c                                 formula weight
      props(17,jd) = 0d0

      do j = 1, itemp
c                                 formula weight
         props(17,jd) = props(17,jd) + cp3(j,jd) * atwt(j) 
c                                 molar amounts of the components
         mols = props(16,jd)*cp3(j,jd)
c                                 mass of the components
         fbulk(j) = fbulk(j) + mols
         gtot = gtot + mols

         if (.not.fluid(jd)) then 
            fbulk1(j) = fbulk1(j) + mols
            gtot1 = gtot1 + mols
         end if 
c                                 molar phase composition
         pcomp(j,jd) = cp3(j,jd)

      end do
c                                 an entity with no mass signals that 
c                                 frendly is using a make definition
c                                 which corresponds to a balanced reaction
      if (gtot.eq.0d0) rxn = .true.

      if (iopt(2).eq.1) then 
c                                 convert molar phase composition to 
c                                 mass % composition:
         do j = 1, itemp
            pcomp(j,jd) = pcomp(j,jd)*atwt(j)*1d2/props(17,jd)
         end do  

      end if 
c                                 shear modulus
      if (.not.fluid(jd)) then 

         call moduli (id,props(5,jd),props(19,jd),props(21,jd),ok)
         if (.not.ok.and.iopt(16).eq.0) shear = .false.  

      else

         props(5,jd)  = 0d0
         props(19,jd) = 0d0    
         props(21,jd) = 0d0

      end if   
c                                 compute g-derivatives for isostatic 
c                                 thermodynamic properties
      dp2 = 5d-2 * p
      dp1 = dp2/1d1
      dp0 = dp1/1d1
            
      g0 = ginc(0d0,0d0,id)
c                                 g0 used only by frendly
      props(11,jd) = g0 
c                                 straight derivatives:
c                                 first order
      if (p-dp0.le.0d0) then 

         v = (ginc(0d0,dp0,id) - g0)/dp0
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment if invalid v
     *   v = (ginc(0d0,dp1,id) - g0)/dp1
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment more if invalid v
     *   v = (ginc(0d0,dp2,id) - g0)/dp2

      else 

         v = (ginc(0d0,dp0,id) - ginc(0d0,-dp0,id))/dp0/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp1.gt.0d0)  
c                                 expand increment if invalid v
     *   v = (ginc(0d0,dp1,id) - ginc(0d0,-dp1,id))/dp1/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp2.gt.0d0)  
c                                 expand increment more if invalid v
     *   v = (ginc(0d0,dp2,id) - ginc(0d0,-dp2,id))/dp2/2d0

      end if 
c                                 in case the evaluating routine fails
c                                 on both calls to ginc 
      s = (ginc(-dt0,0d0,id) - ginc(dt0,0d0,id))/dt0/2d0
c                                 this crap is necessary because 
c                                 optimization or my bad programming
c                                 corrupts ginc with compaq visual fortran.
      g1a = ginc(-dt0,0d0,id)
      g2a =  ginc(dt0,0d0,id)
      dg = g1a-g2a
      ss = dg/dt0/2d0
      s = ss
c
c 
c     write (*,*) s, ss, dg, ginc(-dt0,0d0,id) - ginc(dt0,0d0,id), dt0

      e = g0 + t * s
c                                 second order
      gtt = (ginc(dt1,0d0,id) + ginc(-dt1,0d0,id) - 2d0*g0)/dt1/dt1
      cp = -t*gtt

      if (cp.lt.0d0.or.dabs(cp).gt.1d9)  
c                                 expand increment if invalid cp
     *   cp = -t*(ginc(dt2,0d0,id) + 
     *            ginc(-dt2,0d0,id) - 2d0*g0)/dt2/dt2

      if (cp.lt.0d0.or.dabs(cp).gt.1d9)  
c                                 shrink increment if invalid cp
     *   cp = -t*(ginc(dt0,0d0,id) + 
     *            ginc(-dt0,0d0,id) - 2d0*g0)/dt0/dt0

c                                 volumetric properties only if v is ok:
      if (v.gt.0d0) then 
   
         if (p-dp1.le.0d0) then 
c                                 use forward difference at small p's
            beta = (ginc(0d0,2d0*dp1,id) + g0 - 2d0*ginc(0d0,dp1,id))
     *             /dp1/dp1
            if (beta.lt.-v.or.beta.ge.0d0)
c                                 expand increment if invalid beta
     *      beta = (ginc(0d0,2d0*dp2,id) + g0 - 2d0*ginc(0d0,dp2,id))
     *             /dp2/dp2                                 
            if (beta.lt.-v.or.beta.ge.0d0)
c                                 shrink increment if invalid beta
     *      beta = (ginc(0d0,2d0*dp0,id) + g0 - 2d0*ginc(0d0,dp0,id))
     *             /dp0/dp0   

            alpha = ( ginc( dt1,dp1,id) - ginc( dt1,0d0,id)
     *               -ginc(-dt1,dp1,id) + ginc(-dt1,0d0,id))/dp1/dt1/2d0
            if (alpha.gt.v.or.alpha.le.0d0)
c                                 expand increment if invalid alpha
     *      alpha = ( ginc( dt2,dp2,id) - ginc( dt2,0d0,id)
     *               -ginc(-dt2,dp2,id) + ginc(-dt2,0d0,id))/dp2/dt2/2d0
            if (alpha.gt.v.or.alpha.le.0d0)
c                                 shrink increment if invalid alpha
     *      alpha = ( ginc( dt0,dp0,id) - ginc( dt0,0d0,id)
     *               -ginc(-dt0,dp0,id) + ginc(-dt0,0d0,id))/dp0/dt0/2d0

         else

            beta = (ginc(0d0,dp1,id) + ginc(0d0,-dp1,id) - 2d0*g0)
     *             /dp1/dp1
            if (beta.lt.-v.and.p-dp2.ge.0d0.or.beta.ge.0d0)
c                                 expand increment if invalid beta
     *      beta = (ginc(0d0,dp2,id) + ginc(0d0,-dp2,id) - 2d0*g0)
     *             /dp2/dp2
     
            if (beta.lt.-v.or.beta.ge.0d0)
c                                 shrink increment if invalid beta
     *         beta = (ginc(0d0,dp0,id) 
     *               + ginc(0d0,-dp0,id) - 2d0*g0)/dp0/dp0

            alpha = ( ginc( dt1,dp1,id) - ginc( dt1,-dp1,id)
     *            -ginc(-dt1,dp1,id) + ginc(-dt1,-dp1,id))/dp1/dt1/4d0

            if (alpha.gt.v.or.alpha.le.0d0) then 
c                                 expand increment if invalid alpha
               alpha1 = ( ginc( dt2,dp2,id) - ginc( dt2,-dp2,id)
     *                  - ginc(-dt2,dp2,id) 
     *                  + ginc(-dt2,-dp2,id))/dp2/dt2/4d0

               if (alpha1.gt.v.or.alpha1.le.0d0) then
c                                 shrink increment if invalid alpha
                  alpha2 = ( ginc( dt0,dp0,id) - ginc( dt0,-dp0,id)
     *                      -ginc(-dt0,dp0,id) + ginc(-dt0,-dp0,id))
     *                     /dp0/dt2/4d0

                  if (alpha2.lt.v.and.alpha2.ge.0d0) then 
                     alpha = alpha2
                  end if 

               else 

                  alpha = alpha1

               end if 
            end if      
     
         end if  
c                                 third order derivatives, only need for
c                                 derivatives of seismic props.
         if (p-2d0*dp2.le.0d0) then 

            g1 = ginc(-dt2,0d0,id)
            g2 = ginc( dt2,2d0*dp2,id)
            g3 = ginc( dt2,0d0,id)
            g4 = ginc(-dt2,2d0*dp2,id)
            g5 = ginc(0d0,dp2,id)
            g7 = g3 - g1

            gppp = ((ginc(0d0,4d0*dp2,id) - g0)/2d0
     *              - ginc(0d0,3d0*dp2,id) + g5)/dp2**3
            gppt = (g2 + g3 + 2d0*(ginc(-dt2, dp2,id)-ginc(dt2,dp2,id))
     *              -g4 - g1)/dp2/dp2/dt2/2d0
            gptt = (g2 + g4  + 2d0*(g0 - ginc(0d0,2d0*dp2,id))
     *              -g3 - g1)/2d0/dp2/dt2/dt2

         else 

            g1 = ginc(-dt2,-dp2,id) - ginc( dt2,dp2,id)
            g3 = ginc( dt2,-dp2,id)
            g4 = ginc(-dt2,dp2,id)
            g5 = ginc(0d0,dp2,id) - ginc(0d0,-dp2,id)
            g7 = ginc(-dt2,0d0,id) - ginc(dt2, 0d0,id)
   
            gppp = ((ginc(0d0,2d0*dp2,id) - ginc(0d0,-2d0*dp2,id))/2d0 
     *           - g5)/dp2**3
            gppt = (g3 - g4 + 2d0*g7 - g1)/dp2/dp2/dt2/2d0
            gptt = (g4 - g3 - 2d0*g5 - g1)/2d0/dp2/dt2/dt2
 
         end if 

         gttt = ((ginc(dt2*2d0,0d0,id) - ginc(-dt2*2d0,0d0,id))/2d0 
     *           + g7)/dt2**3

         g7 = (gtt*beta-alpha**2)**2

         if (g7.ne.0d0) then 
c                                 temperature derivative of the adiabatic bulk modulus:
            props(18,jd) = (((v*gppt-alpha*beta)*gtt
     *                      -(2d0*v*gptt-alpha**2)*alpha)*gtt
     *                      +v*gttt*alpha**2)/g7
c                                 pressure derivative of the adiabatic bulk modulus:
            props(20,jd) = (((v*gppp-beta**2)*gtt
     *                   +(alpha*beta-2d0*v*gppt)*alpha)*gtt
     *                   +v*gptt*alpha**2)/g7
         else 

            props(18,jd) = nopt(7)
            props(20,jd) = nopt(7)

         end if 

      end if 
c                                 -------------------------------------
c                                 up to this point beta = d2g/dp2
c                                 and alpha = d2g/dp2 now convert 
c                                 to their normal forms:
      if (v.le.0d0) then 

         sick(jd) = .true.
         v = nopt(7)
         beta = nopt(7)
         alpha = nopt(7)
         rho = nopt(7)
      
      else 

         beta = -beta/v
         alpha = alpha/v
         rho = props(17,jd)/v*1d2

c                                 ideal gas beta = 1/p           
         if (beta.gt.v.or.beta.lt.0d0) then
            beta = nopt(7)
            sick(jd) = .true.
         end if
c                                 aug 28, 2007, removed check on alpha to 
c                                 accomodate -alpha's generated by landau 
c                                 transition models. ideal gas alpha = 1/t
         if (alpha.gt.v) then 
            alpha = nopt(7)
            sick(jd) = .true.
         end if 

      end if  

      if (cp.gt.1d9.or.cp.lt.0d0) then
         cp = nopt(7)
         sick(jd) = .true.
      end if 

      props(2,jd) = e
      props(12,jd) = cp
      props(15,jd) = s
      props(1,jd) = v
      props(13,jd) = alpha
      props(14,jd) = beta 
      props(10,jd) = rho  

      if (.not.sick(jd).and..not.rxn) then
c                                 gruneisen parameter
         props(3,jd) = v/(cp*beta/alpha - t*alpha*v)
c                                 aug 28, 2007, removed check on gruneisen to 
c                                 accomodate -alpha's generated by landau 
c                                 transition models (prop(13,jd))
c        if (props(3,jd).le.0d0) sick(jd) = .true.

c                                 adiabatic bulk modulus
         props(4,jd) = (1d0 + t*alpha*props(3,jd))/beta

         if (props(4,jd).le.0d0) sick(jd) = .true.

         if (.not.fluid(jd).and.iopt(16).gt.0) then 
c                                 use poisson ratio estimates if iopt(16).ne.0
            if ((iopt(16).eq.1.and..not.ok).or.iopt(16).eq.2) then

               if (volume.and..not.sick(jd)) then
 
                  props(5,jd)  = nopt(16)*props(4,jd)
                  props(19,jd) = nopt(16)*props(18,jd) 
                  props(21,jd) = nopt(16)*props(20,jd)

                  if (iopt(16).eq.1) then 
                     ppois = .true.
                     pois = .true.
                  end if 

               else

                  shear = .false.
                  props(5,jd) = nopt(7)
                  props(19,jd) = nopt(7)    
                  props(21,jd) = nopt(7)

               end if 

            end if

         end if  

      end if 

      if (sick(jd)) then

         props(3,jd) = nopt(7)
         props(4,jd) = nopt(7)

         volume = .false.

         if (.not.fluid(jd).and..not.ok) shear = .false.

         if (.not.fluid(jd)) ssick = .true.

      end if 
c                                 seismic properties
      if (.not.sick(jd).and..not.rxn) then 

         units = dsqrt(1d5)/1d3
         r43   = 4d0/3d0
c                                 sound velocity
         root = props(4,jd)/rho
         props(6,jd) = dsqrt(root) * units
c                                 sound velocity T derivative
         props(22,jd) = (props(18,jd) + props(4,jd) * alpha) 
     *                  / dsqrt(root) / rho / 2d0 * units
c                                 sound velocity P derivative
         props(25,jd) = (props(20,jd) - props(4,jd) * beta) 
     *                  / dsqrt(root) / rho / 2d0 * units

c                                 p-wave velocity
         root = (props(4,jd)+r43*props(5,jd))/rho

         if (root.ge.0d0) then 

            props(7,jd) = dsqrt(root)*units
c                                 p-wave velocity T derivative
            props(23,jd) = (props(18,jd) + r43 * 
     *                     (props(19,jd) + alpha * props(5,jd)) 
     *                    + props(4,jd) * alpha) / 
     *                      dsqrt(root) / rho / 2d0 * units
c                                 p-wave velocity P derivative
            props(26,jd) = (props(20,jd) + r43 *
     *                     (props(21,jd) - beta * props(5,jd)) 
     *                    - props(4,jd) * beta) /
     *                      dsqrt(root) / rho / 2d0 * units

         else 

            props(7,jd) = nopt(7)
            props(23,jd) = nopt(7)
            props(26,jd) = nopt(7)
            shear = .false.

         end if 

         if (.not.fluid(jd)) then 
c                                 s-wave velocity
            root = props(5,jd)/rho

            if (root.gt.0d0) then 

               props(8,jd) = dsqrt(root)*units
               props(24,jd)= (props(19,jd) + props(5,jd) * alpha)
     *                     / dsqrt(root) / rho / 2d0 * units
               props(27,jd)= (props(21,jd) - props(5,jd) * beta)
     *                     / dsqrt(root) / rho / 2d0 * units
            else

               props(8,jd) = nopt(7)
               props(24,jd) = nopt(7)
               props(27,jd) = nopt(7)
               shear = .false.

            end if 
c                                 vp/vs
            if (props(8,jd).ne.0d0) then 
               props(9,jd) = props(7,jd)/props(8,jd)
            else
               props(9,jd) = nopt(7)
            end if 

         else 

            props(8,jd) = 0d0
            props(24,jd) = 0d0 
            props(27,jd) = 0d0 

         end if 

      else 

         do j = 3, 9
            props(j,jd) = nopt(7)
         end do 

      end if 
c                                 check and warn if necessary for negative
c                                 expansivity
      if (.not.sick(jd)) then 

         if (alpha.le.0d0.and.iwarn1.lt.11) then

            write (*,1030) t,p,pname(jd)
            iwarn1 = iwarn1 + 1
            if (iwarn1.eq.11) call warn (49,r,179,'GETPHP') 
         end if 

      end if

      if (ppois.and.iwarn2.lt.11) then

         if (pois) then 
            iwarn2 = iwarn2 + 1
            write (*,1040) t,p,pname(jd)
         end if
 
         if (iwarn2.eq.11) call warn (49,r,178,'GETPHP')

      end if 
c                                 accumulate aggregate totals, some
c                                 totals may be incomplete if volume or 
c                                 shear is false for an individual phase
c                                 weighting scheme for seismic velocity
c                                 chi = 1 -> voigt 0 -> reuss 0.5 -> VRH
      chi = nopt(6)
      chi1 = 1d0 - chi

      if (iam.ne.5) then
c                                 weighting factor for molar properties
         mols = props(16,jd)

      else 
c                                 if frendly use reaction coefficients
         mols = vnu(jd)

      end if 
c                                 vol of phase per mole of system
      v = v*mols
c                                 system molar volume
      psys(1)  = psys(1)  + v
c                                 molar enthalpy
      psys(2)  = psys(2)  + e*mols 
c                                 molar gibbs energy
      psys(11) = psys(11) + props(11,jd)*mols 
c                                 molar heat capacity
      psys(12) = psys(12) + cp*mols 
c                                 expansivity
      psys(13) = psys(13) + alpha*v 
c                                 compressibility
      psys(14) = psys(14) + beta*v
c                                 molar entropy
      psys(15) = psys(15) + props(15,jd)*mols 
c                                 moles of assemblage
      psys(16) = psys(16) + mols
c                                 mass of assemblage 
      psys(17) = psys(17) + props(17,jd)*mols
       
      if (volume.and..not.rxn) then 
c                                 gruneisen
         psys(3) = psys(3) + v*props(3,jd)*chi
         pgeo(3) = pgeo(3) + v/props(3,jd)
c                                 Aggregate Bulk Modulus                                
         psys(4) = psys(4) + v*props(4,jd)*chi
         pgeo(4) = pgeo(4) + v/props(4,jd)
c                                 Aggregate Bulk Modulus T-derivative                                
         psys(18) = psys(18) + v*props(18,jd)*chi
         pgeo(18) = pgeo(18) + v/props(18,jd)
c                                 Aggregate Bulk Modulus P-derivative                                
         psys(20) = psys(20) + v*props(20,jd)*chi
         pgeo(20) = pgeo(20) + v/props(20,jd)
c                                 aggregate shear props only if
c                                 shear mod is available for all phases.
         if (shear) then 
c                                 Aggregate Shear Modulus
            psys(5) = psys(5) + v*props(5,jd)*chi
            if (.not.aflu) pgeo(5) = pgeo(5) + v/props(5,jd)
c                                 Aggregate Shear Modulus T-derivative
            psys(19) = psys(19) + v*props(19,jd)*chi
            if (.not.aflu.and.props(19,jd).ne.0d0) 
     *                        pgeo(19) = pgeo(19) + v/props(19,jd)
c                                 Aggregate Shear Modulus P-derivative
            psys(21) = psys(21) + v*props(21,jd)*chi
            if (.not.aflu.and.props(21,jd).ne.0d0) 
     *                        pgeo(21) = pgeo(21) + v/props(21,jd)

         end if 

      end if

      if (aflu.and..not.rxn) then 
c                                 assemblage includes fluid
         if (.not.fluid(jd)) then
c                                 get total without fluid
            psys1(1)  = psys1(1) + v 
c                                 molar enthalpy
            psys1(2)  = psys1(2) + e*mols 
c                                 molar heat capacity
            psys1(12) = psys1(12) + cp*mols 
c                                 molar expansivity
            psys1(13) = psys1(13) + alpha*v 
c                                 molar compressibility
            psys1(14) = psys1(14) + beta*v 
c                                 molar entropy
            psys1(15) = psys1(15) + s*mols 
c                                 total number of moles of phases
            psys1(16) = psys1(16) + mols

            psys1(17) = psys1(17) + props(17,jd)*mols

            if (.not.ssick) then 
c                                 gruneisen
               psys1(3) = psys1(3) + v*props(3,jd)*chi
               pgeo1(3) = pgeo1(3) + v/props(3,jd)
c                                 Aggregate Bulk Modulus                                
               psys1(4) = psys1(4) + v*props(4,jd)*chi
               pgeo1(4) = pgeo1(4) + v/props(4,jd)
c                                 Aggregate Bulk Modulus T-derivative                                
               psys1(18) = psys1(18) + v*props(18,jd)*chi
               pgeo1(18) = pgeo1(18) + v/props(18,jd)
c                                 Aggregate Bulk Modulus P-derivative                                
               psys1(20) = psys1(20) + v*props(20,jd)*chi
               pgeo1(20) = pgeo1(20) + v/props(20,jd)

               if (shear) then 
c                                 Aggregate Shear Modulus
                  psys1(5) = psys1(5) + v*props(5,jd)*chi
                  pgeo1(5) = pgeo1(5) + v/props(5,jd)
c                                 Aggregate Shear Modulus T-derivative
                  psys1(19) = psys1(19) + v*props(19,jd)*chi
                  if (props(19,jd).ne.0d0) 
     *                       pgeo1(19) = pgeo1(19) + v/props(19,jd)
c                                 Aggregate Shear Modulus P-derivative
                  psys1(21) = psys1(21) + v*props(21,jd)*chi
                  if (props(21,jd).ne.0d0) 
     *                       pgeo1(21) = pgeo1(21) + v/props(21,jd)

               end if 

            end if

         end if
 
      end if 

1030  format (/,'**warning ver179** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'the effective expansivity of: ',a,/,'is negative. ',
     *        'Most probably this is because of a Landau ordering ',
     *        'model. The Gruneisen',/,'thermal parameter and seismic',
     *        ' velocities for this phase should be considered ',
     *        'with caution.',/)

1040  format (/,'**warning ver178** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'the shear modulus of: ',a,/,'is missing or invalid ',
     *        'and has been estimated from the default poisson ',
     *        'ratio ',/)

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
c      only called by WERAMI.
c gsol1 is identical to gsol but can only been called after gall and is 
c      only called by VERTEX and MEEMUM. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k,id

      double precision omega, gproj, hpmelt, gmelt, gfluid, gzero, g, 
     *                 dg, gex

      external gex

      integer jend
      common/ cxt23 /jend(h9,k12)

      double precision r,tr,pr,ps,p,t,xco2,u1,u2
      common/ cst5   /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ispec
      common/ cxt8 /ispec(h9,m4)
c----------------------------------------------------------------------
      if (id.lt.0) then 

         call gphase (-id,g)
         gsol = g

      else 

         g = 0d0
c                                 evaluate dqf coefficients
         call setdqf (id)

         call setw (id) 

         if (ksmod(id).eq.2.or.ksmod(id).eq.3) then 
c                                 -------------------------------------
c                                 macroscopic formulation for normal solutions.
            call gdqf (id,g,y) 
c                                 add entropy and excess contributions
            g = g - t * omega(id,y) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id) 
               g = g + y(k) * gproj (jend(id,2+k))
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
               g = g + gproj(jend(id,2+k)) * p0a(k)
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
               g = g + y(k) * gproj (jend(id,2+k))
            end do 
c                                 get the speciation energy effect
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
               g = g + gproj (jend(id,2+k)) * p0a(k) 
            end do 
c                                 get the dqf
            call gdqf (id,g,p0a)
c                                 and excess contributions
            g = g - t * omega(id,p0a) + gex(id,p0a)

         else if (ksmod(id).eq.23) then 

             write (*,*) 'toop samis model not coded'

         else if (ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp melt model         
            call gdqf (id,g,y) 
            g = g - t * hpmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.25) then 
c                                 -------------------------------------
c                                 ghiorso pmelt model  
            call gdqf (id,g,y) 

            g = g - t * gmelt(id) + gex(id,y)
c                                 get mechanical mixture contribution
            do k = 1, mstot(id)  
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.26) then 
c                                 ------------------------------------
c                                 andreas salt model
            call hcneos (g,y(1),y(2),y(3))

            do k = 1, 3
               g = g + y(k) * gproj (jend(id,2+k))
            end do 

         else if (ksmod(id).eq.27) then 

            do k = 1, mstot(id)
               if (y(k).gt.0d0)   
     *            g = g + (gproj(jend(id,2+k))+r*t*dlog(y(k)))*y(k) 
            end do 

         else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos
            do k = 1, 2
               g = g + gzero(jend(id,2+k))*y(k)
            end do 

            g = g + gfluid(y(ispec(id,1)))

         else 

            write (*,*) 'what the **** am i doing here?'
            stop

         end if 

      gsol = g 

      end if 

      end

      subroutine shearm (mu,mut,mup,id)
c-----------------------------------------------------------------------
c shearm returns a linear model for the adiabatic shear modulus
c relative to the current pressure and temperature.

c three cases:

c make(id) = non-zero, use make definition to compute shear modulus.

c iemod = 1, linear model is input

c iemod = 2, shear modulus is known as a function of (V,T), then
c computed by centered finite differences.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision mu,mut,mup,mu2,dt,dp,g,ginc

      double precision smu
      common/ cst323 /smu

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer make
      common / cst335 /make(k10)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      save dt,dp
      data dt,dp/5d0,50d0/
c-----------------------------------------------------------------------

      if (make(id).ne.0) then 

         call makmod (id,mu,mut,mup)

      else if (iemod(id).eq.1) then 

         mu  = emod(1,id) + (p-pr)*emod(2,id) + (t-tr)*emod(3,id)
         mut = emod(3,id)
         mup = emod(2,id)

      else if (iemod(id).eq.2) then 
c                                 by calling ginc a call to
c                                 stixrudes EoS for the adiabatic
c                                 shear modulus is implicit (cst323)
         g = ginc(0d0,0d0,-id)
         mu = smu
c                                 temperature derivative
         g = ginc(dt,0d0,-id)
         mu2 = smu
         g = ginc(-dt,0d0,-id)         
         mut = (mu2 - smu)/dt/2d0

         if (p-dp.gt.0d0) then 
c                                 centered pressure derivative
            g = ginc(0d0,dp,-id)
            mu2 = smu
            g = ginc(0d0,-dp,-id)         
            mup = (mu2 - smu)/dp/2d0

         else 

            g = ginc(0d0,dp,-id)
            mu2 = smu
            g = ginc(0d0,2d0*dp,-id)         
            mup = (mu2 - smu)/dp/2d0

         end if 

      end if          

      end 

      subroutine makmod (id,mu,mut,mup)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id, jd

      double precision mu, pmu, mut, pmut, mup, pmup

      double precision mkcoef, mdqf

      integer mknum, mkind
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16)

      integer make
      common / cst335 /make(k10)
c-----------------------------------------------------------------------

      jd = make(id)

      mu = 0d0
      pmut = 0d0
      pmup = 0d0 

c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         call shearm (pmu,pmut,pmup,mkind(jd,i))

         mu = mu + mkcoef(jd,i) * pmu
         mut = mut + mkcoef(jd,i) * pmut
         mup = mup + mkcoef(jd,i) * pmup

      end do 

      end

      subroutine moduli (ids,mu,mut,mup,ok) 
c-----------------------------------------------------------------------
c subroutine moduli determines shear moduli (mods) for entity ids, returns
c ok = false if moduli are unavailable.

c jmod is the number of cpds for which it is possible to calculate coeffs.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision mu, pmu, mut, pmut, mup, pmup

      integer i, ids

      logical ok
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer jend
      common/ cxt23 /jend(h9,k12)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer iemod,kmod
      logical smod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),iemod(k10),kmod

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
  
      ok = .true.

      mu = 0d0
      mut = 0d0
      mup = 0d0

      if (ids.le.0) then 

         if (iemod(-ids).ne.0) then

            call shearm (mu,mut,mup,-ids)

         else

            ok = .false.

         end if 

      else 

         if (smod(ids)) then 

            if (lrecip(ids)) then
c                                 get the p0a coordinates (amounts of 
c                                 the independent disordered endmembers)     
               call getpp (ids) 

               do i = 1, lstot(ids)

                  call shearm (pmu,pmut,pmup,jend(ids,2+i))

                  mu = mu + p0a(i) * pmu
                  mut = mut + p0a(i) * pmut
                  mup = mup + p0a(i) * pmup

               end do

            else 

c                                 for solutions with no dependent endmembers
c                                 the y coordinates can be used to compute 
c                                 the composition. for speciation models
c                                 (ksmod = 6) this assumes xtoy has been 
c                                 called before moduli, so that y is the 
c                                 unspeciated composition (this will not be 
c                                 the case if gsol has been called, as might
c                                 happen if shearm calls gsol to evaluate a 
c                                 speciation model using stixrude's EoS).
               do i = 1, mstot(ids)

                  call shearm (pmu,pmut,pmup,jend(ids,2+i))

                  mu  = mu + y(i) * pmu
                  mut = mut + y(i) * pmut
                  mup = mup + y(i) * pmup

               end do
 
            end if 

            if (mu.lt.0d0) then 
               mu = nopt(7)
               ok = .false.
            end if 

            if (mut.gt.0d0) then 
               mut = nopt(7)
               ok = .false.
            end if

            if (mup.lt.0d0) then 
               mup = nopt(7)
               ok = .false.
            end if
       
         else

            ok = .false.

         end if 

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

      double precision gee, fo2, fs2

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn

      double precision fh2o,fco2
      common/ cst11 /fh2o,fco2

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      call gphase (id,gee)

      gee = gee + r * t * dlog(act(id))

      if (ifyn.eq.0) then 
c                                 this is a quick fix that will
c                                 call the fluid routine way more 
c                                 than necessary.
         call cfluid (fo2,fs2)

         if (id.eq.idf(3)) then 

            gee = gee + r*t*fo2

         else if (id.eq.idf(1)) then 
        
            gee = gee + r*t*fh2o
         
         else if (id.eq.idf(2)) then 
        
            gee = gee + r*t*fh2o

         end if
 
      end if 

      gfrnd = gee

      end 

