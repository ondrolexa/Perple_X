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
      common/ cxt81 /gtot,fbulk(k5),gtot1,fbulk1(k5)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k5,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(2),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      double precision atwt
      common/ cst45 /atwt(k0)

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume

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
      common/ cst324 /pcomp(k5,k5)

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
      common/ cxt15 /cp3(k5,k5),amt(k5),kkp(k5),np,ncpd,ntot
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

      logical gflu,aflu,fluid,shear,lflu,volume
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume

      double precision bg
      common/ cxt19 /bg(k5,k2)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k5),gtot1,fbulk1(k5)

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
