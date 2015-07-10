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
 
      character cprop*18, text*240

      integer i,j,l,lu,id

      double precision poiss

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

      integer spct
      double precision ysp
      character*8 spnams
      common/ cxt34 /ysp(m4,k5),spct(h9),spnams(m4,h9)

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(240)

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

         if (iopt(2).eq.0) then 

            write (lu,1030) pname(i), 
c                                 weight %
     *                      props(17,i)*props(16,i)/psys(17)*1d2,
c                                 vol %
     *                      props(1,i)*props(16,i)/psys(1)*1d2,
c                                 mol %
     *                      props(16,i)/psys(16)*1d2,
c                                 mol
     *                      props(16,i),
c                                 molar or weight composition
     *                      (pcomp(l,i), l = 1, icomp)

         else 

            write (lu,1031) pname(i), 
c                                 weight %
     *                      props(17,i)*props(16,i)/psys(17)*1d2,
c                                 vol %
     *                      props(1,i)*props(16,i)/psys(1)*1d2,
c                                 mol %
     *                      props(16,i)/psys(16)*1d2,
c                                 mol
     *                      props(16,i),
c                                 molar or weight composition
     *                      (pcomp(l,i), l = 1, icomp)

         end if 

      end do 

      if (lopt(21).and.np.gt.0) then 
c                                 species_ouput
         write (lu,'(/,a,/)') 'Phase speciation (molar proportions):'

         do i = 1, np 

            id = kkp(i) 

            write (text,'(20(a,a,f7.5,a))')
     *            (spnams(j,id),': ',ysp(j,i),', ', j =1,spct(id))

            call deblnk (text)

            write (lu,'(1x,a,4x,240a)') pname(i), 
     *                                 (chars(j), j = 1, length)

         end do 

      end if 

      write (lu,1160)
c                                 phase/system summary, normal thermo:
      do i = 1, ntot
c                                 N, H, S, V, Cp, alpha, beta, density
         write (lu,1170) pname(i),props(17,i),props(2,i),props(15,i),
     *                   props(1,i),(props(j,i),j=12,14),props(28,i),
     *                   props(10,i)
      end do

      write (lu,1170) 'System        ',psys(17),psys(2),psys(15),
     *                psys(1),(psys(j),j=12,14),psys(28),psys(10)
      if (aflu) write (lu,1170) 'System - fluid',psys1(17),psys1(2),
     *                psys1(15),psys1(1),(psys1(j),j=12,14),psys1(28),
     *                psys1(10)


      if (iopt(14).gt.0) then 
c                                 phase/system summary, seismic:
         write (lu,1190)     

         do i = 1, ntot
 
            write (lu,1200) pname(i), (props(j,i), j = 3, 8),
     *                      poiss(props(7,i),props(8,i))
         end do

         write (lu,1200) 'System        ',(psys(j), j = 3, 8),
     *                                    poiss(psys(7),psys(8))

         if (aflu) write (lu,1200) 'System - fluid',(psys1(j), j = 3, 8)
     *                             ,poiss(psys1(7),psys1(8))

         if (iopt(14).eq.2) then 
c                                 phase/system summary, seismic derivatives:
            write (lu,1240)

            do i = 1, ntot
               write (lu,1250) pname(i),props(18,i),props(20,i),
     *                         props(19,i),props(21,i),props(22,i),
     *                         props(25,i),props(23,i),props(26,i),
     *                         props(24,i),props(27,i)
            end do

            write (lu,1250) 'System        ',psys(18),psys(20),psys(19),
     *                     psys(21),psys(22),psys(25),psys(23),psys(26),
     *                     psys(24),psys(27)

            if (aflu) write (lu,1250) 'System - fluid',psys1(18),
     *                psys1(20),psys1(19),psys1(21),psys1(22),psys1(25),
     *                psys1(23),psys1(26),psys1(24),psys1(27)

          end if 

      end if 

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
         write (lu,1130) (cname(i), i = 1, jbulk)
         write (lu,1140) (mu(i), i = 1, jbulk)
         write (lu,1071) 2, jbulk - ntot + 2 
      else 
         write (lu,1070) 2, jbulk - ntot + 2 
      end if 

1000  format (/,40('-'),//,'Stable phases at:')
1020  format (/,'Phase Compositions (',a,'):',
     *        /,19x,'wt %',6x,'vol %',5x,'mol %',5x,'mol  ',
     *          5x,20(1x,a,3x))
1030  format (1x,a,3x,3(f6.2,4x),g9.3,1x,20(f8.5,1x))
1031  format (1x,a,3x,3(f6.2,4x),g9.3,1x,20(f8.3,1x))
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
     *         'Cp(J/K)',6x,'Alpha(1/K)',2x,'Beta(1/bar)',4x,'Cp/Cv',5x,
     *         'Density(kg/m3)')
1170  format (1x,a,1x,f9.2,3x,13(g12.5,1x),3x,f7.4)
1190  format (/,'Seismic Properties:'
     *        /,17x,'Gruneisen_T',7x,'Ks(bar)',7x,'Mu(bar)',
     *        4x,'V0(km/s)',5x,'Vp(km/s)',5x,'Vs(km/s)',5x,
     *        'Poisson ratio')
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

      logical sick(i8), nodata, ssick, ppois, bulkg, bsick
c                                 x-coordinates for the assemblage solutions
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 
      integer ifp
      logical fp
      common/ cxt32 /ifp(k10), fp(h9)
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

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iam
      common/ cst4 /iam

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)
c----------------------------------------------------------------------
c                                 logarithmic_p option
10    if (lopt(14)) p = 1d1**p 

      nodata = .false. 

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
 
            do i = 1, jbulk
               mu(i) = 0d0
            end do 

            do i = 1, ijpt

               kd = igrd(itri(i),jtri(i))

               do j = 1, jbulk
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
c                                 initialize system props/flags
      call insysp (ssick,ppois,bulkg,bsick)

      do i = 1, ntot

         ids = kkp(i)

         if (i.le.np) then 

            if (fp(ids)) then 
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

            if (ifp(-ids).gt.0.or.ifp(-ids).lt.0.and.lopt(6)) then 
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
         call getphp (ids,i,sick,ssick,ppois,bulkg,bsick)   

c         if ((ssick.or.bsick.or.sick(i)).and.ijpt.gt.1) then 
c            wt(1) = 1d0
c            ijpt = 1
c            write (*,*) 'turned off at :', t, p, dlog10(p)
c            goto 10 
c         end if  

      end do 
c                                 compute aggregate properties:
      call gtsysp (sick,ssick,bulkg,bsick)

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
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
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

      double precision omega, hpmelt, gmelt, gfluid, gzero, g, x1(5), 
     *                 dg, gex, slvmlt, gfesi, gcpd, gerk, gfecr1

      external gphase, omega, hpmelt, gmelt, gfluid, gzero, gex, slvmlt,
     *         gfesi, gerk, gfecr1

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
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer jspec
      common/ cxt8 /jspec(h9,m4)
c----------------------------------------------------------------------
      if (id.lt.0) then 

         gsol = gcpd (-id,.true.)

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
               g = g + gcpd (jend(id,2+k),.true.) * p0a(k) 
            end do 
c                                 get the dqf
            call gdqf (id,g,p0a)
c                                 and excess contributions
            g = g - t * omega (id,p0a) + gex (id,p0a)


         else if (ksmod(id).eq.24) then 
c                                 -------------------------------------
c                                 hp melt model         
            call gdqf (id,g,y) 

            g = g - t * hpmelt (id) + gex (id,y)
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

         else if (ksmod(id).eq.41) then 
c                                 hybrid MRK ternary COH fluid
            call rkcoh6 (y(2),y(1),g) 

            do k = 1, nstot(id) 
               g = g + gcpd(jend(id,2+k),.true.) * y(k)
            end do 

         else if (ksmod(id).eq.40) then 
c                                 MRK silicate vapor
            do k = 1, nstot(id) 
               g = g + gzero(jend(id,2+k)) * y(k)
               x1(k) = y(k)
            end do 

            g = g + gerk(x1)

         else if (ksmod(id).eq.0) then 
c                                 ------------------------------------
c                                 internal fluid eos
            do k = 1, 2
               g = g + gzero (jend(id,2+k))*y(k)
            end do 

            g = g + gfluid (y(jspec(id,1)))

         else 

            write (*,*) 'what the **** am i doing here?'
            stop

         end if 

         gsol = g 

      end if 

      end

      subroutine getspc (id,jd)
c-----------------------------------------------------------------------
c getspc loads independent endmember fractions of a solution id into a 
c single compositional array y(m4). it is used only for output and 
c must be called after function gsol (id) and routine setspc (id).
c jd is the local phase counter from the calling routine. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer k, id, jd
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer spct
      double precision ysp
      character*8 spnams
      common/ cxt34 /ysp(m4,k5),spct(h9),spnams(m4,h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      character specie*4
      integer isp, ins
      common/ cxt33 /isp,ins(nsp),specie(nsp)

      double precision xs,g,v
      common/ cstcoh /xs(nsp),g(nsp),v(nsp)
c----------------------------------------------------------------------

      if (ksmod(id).eq.2.or.ksmod(id).eq.3.or.ksmod(id).eq.24.or.
     *    ksmod(id).eq.25.or.ksmod(id).eq.26.or.ksmod(id).eq.28) then 
c                                 macroscopic formulation for normal solutions (2,3) and
c                                 hp melt model (24)
c                                 ghiorso melt model (25)
c                                 andreas salt model (26)
c                                 high T melt model (28)
         do k = 1, spct(id) 
            ysp(k,jd) = y(k)
         end do 

      else if ((lrecip(id).and.lorder(id)).or.lorder(id)) then 

         do k = 1, spct(id) 
            ysp(k,jd) = pa(k)
         end do 

      else if (lrecip(id)) then 

         do k = 1, spct(id) 
            ysp(k,jd) = p0a(k)
         end do 

      else if (ksmod(id).eq.29.or.ksmod(id).eq.32) then 
c                                 BCC Fe-Si Lacaze and Sundman (29) 
c                                 BCC Fe-Cr Andersson and Sundman (32)
         spct(id) = 4
c                                 need to correct routines to give o/d
         do k = 1, spct(id) 
            ysp(k,jd) = 0d0
         end do 

      else if (ksmod(id).eq.40.or.ksmod(id).eq.41.or.ksmod(id).eq.0) 
     *        then 
c                                 fluid speciation with known routines:
c                                 MRK silicate vapor (40) 
c                                 hybrid MRK ternary COH fluid (41)
         do k = 1, spct(id) 
            ysp(k,jd) = xs(ins(k))
         end do          

      end if 

      end

      recursive subroutine shearm (mu,mut,mup,ks,kst,ksp,id)
c-----------------------------------------------------------------------
c shearm returns a linear model for the adiabatic shear/bulk modulus
c relative to the current pressure and temperature.

c three cases:

c make(id) = non-zero => use make definition to compute shear modulus.

c eos = 5 or 6 => shear modulus is known as a function of (V,T), then
c computed by centered finite differences.

c iemod = 1 or 2 => linear model is input
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision mu,mut,mup,mu2,ks,kst,ksp,dt,dp,g,ginc

      double precision smu
      common/ cst323 /smu

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer make
      common / cst335 /make(k10)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      integer eos
      common/ cst303 /eos(k10)

      save dt,dp
      data dt,dp/5d0,50d0/
c-----------------------------------------------------------------------

      if (make(id).ne.0) then 

         call makmod (id,mu,mut,mup,ks,kst,ksp)

      else if (eos(id).eq.5.or.eos(id).eq.6) then 
c                                 by calling ginc a call to
c                                 stixrudes EoS for the adiabatic
c                                 shear modulus is implicit (cst323)
         g = ginc(0d0,0d0,-id)
         mu = smu
c                                 set ks = 0 in case implicit 
c                                 bulk modulus is T.
         ks = 0d0
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

      else if (iemod(id).ne.0) then 

         mu  = emod(1,id) + (p-pr)*emod(2,id) + (t-tr)*emod(3,id)
         mut = emod(3,id)
         mup = emod(2,id)

         ks  = emod(4,id) + (p-pr)*emod(5,id) + (t-tr)*emod(6,id)
         kst = emod(6,id)
         ksp = emod(5,id)

      end if          

      end 

      subroutine makmod (id,mu,mut,mup,ks,kst,ksp)
c-----------------------------------------------------------------------
c gmake computes and sums the component g's for a make definition.
c the component g's may be calculated redundantly because gmake is
c called by gcpd, which in turn may be called by routines that call
c for a single g (e.g., gphase). 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id, jd

      double precision mu, pmu, mut, pmut, mup, pmup, 
     *                 ks, pks, kst, pkst, ksp, pksp

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

      ks = 0d0
      pkst = 0d0
      pksp = 0d0 
c                                compute the sum of the component g's
      do i = 1, mknum(jd)

         call shearm (pmu,pmut,pmup,pks,pkst,pksp,mkind(jd,i))

         mu = mu + mkcoef(jd,i) * pmu
         mut = mut + mkcoef(jd,i) * pmut
         mup = mup + mkcoef(jd,i) * pmup

         ks = ks + mkcoef(jd,i) * pks
         kst = kst + mkcoef(jd,i) * pkst
         ksp = ksp + mkcoef(jd,i) * pksp

      end do 

      end

      subroutine moduli (ids,mu,mut,mup,ks,kst,ksp,ok) 
c-----------------------------------------------------------------------
c subroutine moduli determines shear moduli (mods) for entity ids, returns
c ok = false if moduli are unavailable.

c jmod is the number of cpds for which it is possible to calculate coeffs.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      double precision mu, pmu, mut, pmut, mup, pmup, ks, ksp, kst,
     *                 pks, pksp, pkst

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
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
  
      ok = .true.

      mu = 0d0
      mut = 0d0
      mup = 0d0

      ks = 0d0
      kst = 0d0
      ksp = 0d0

      if (ids.le.0) then 

         if (iemod(-ids).ne.0) then

            call shearm (mu,mut,mup,ks,kst,ksp,-ids)

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

                  call shearm (pmu,pmut,pmup,
     *                         pks,pkst,pksp,jend(ids,2+i))

                  mu = mu + p0a(i) * pmu
                  mut = mut + p0a(i) * pmut
                  mup = mup + p0a(i) * pmup

                  if (.not.pmod(ids)) cycle 

                  ks = ks + p0a(i) * pks
                  kst = kst + p0a(i) * pkst
                  ksp = ksp + p0a(i) * pksp

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

                  call shearm (pmu,pmut,pmup,
     *                         pks,pkst,pksp,jend(ids,2+i))

                  mu  = mu + y(i) * pmu
                  mut = mut + y(i) * pmut
                  mup = mup + y(i) * pmup

                  if (.not.pmod(ids)) cycle

                  ks  = ks + y(i) * pks
                  kst = kst + y(i) * pkst
                  ksp = ksp + y(i) * pksp

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

      double precision gee, fo2, fs2, gcpd
 
      external gcpd

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn

      double precision fh2o,fco2,funk
      common/ cst11 /fh2o,fco2,funk

      integer idf
      double precision act
      common/ cst205 /act(k7),idf(3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer eos
      common/ cst303 /eos(k10)
c-----------------------------------------------------------------------

      gee = gcpd (id,.false.) + r * t * dlog(act(id))

      if (ifyn.eq.0.and.eos(id).lt.100) then 
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

      double precision function poiss (vp,vs)
 
      implicit none

      include 'perplex_parameters.h'

      double precision vp, vs

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      if (isnan(vp).or.isnan(vs)) then 
         poiss = nopt(7)
      else if (vs.eq.0d0) then 
         poiss = 0.5d0
      else 
         poiss =  0.5d0*((vp/vs)**2-2d0)/((vp/vs)**2-1d0)
      end if 

      end 

      subroutine getphp (id,jd,sick,ssick,ppois,bulkg,bsick)
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

c 28 - heat capacity ratio (cp/cv)

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

      logical ok, sick(i8), ssick, pois, ppois, bulk, bulkg, bsick, 
     *        lshear, fow

      integer id, jd, iwarn1, iwarn2, iwarn3, j, itemp, m

      character*14 wname1, wname2, wname3

      double precision dt0, dt1, dt2, g0, gpt, gpp, dp0, dp1, dp2, e, 
     *                 alpha, v, ginc, dt,
     *                 beta, cp, s, rho, gtt, g1, g2, g3, g4, g5,
     *                 g7, gppp, gppt, gptt, gttt, mols, root, poiss

      external poiss, ginc

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

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)

      double precision units, r13, r23, r43, r59, r1, r2
      common/ cst59 /units, r13, r23, r43, r59, r1, r2

      double precision vrt
      integer irt
      logical sroot
      common/ rkroot /vrt,irt,sroot

      double precision y,g,vsp
      common / cstcoh /y(nsp),g(nsp),vsp(nsp)

      double precision pv, pvv
      integer iroots
      logical switch, rkmin, min
      common/ rkdivs /pv,pvv,iroots,switch,rkmin,min

      double precision vp,vvp
      integer rooti
      common/ srkdiv /vp(3),vvp(3),rooti(3)

      save dt
      data dt /.5d0/

      save iwarn1, iwarn2, iwarn3, wname1, wname2, wname3
      data iwarn1, iwarn2, iwarn3, wname1, wname2, wname3 /3*0,3*' '/
c----------------------------------------------------------------------
      sick(jd) = .false.
      pois     = .false.
      lshear   = .false.
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
c                                 bulk modulus flag, if false use explicit form
      bulk = .true.
c                                 shear modulus
      if (.not.fluid(jd)) then 

         call moduli (id,props(5,jd),props(19,jd),props(21,jd),
     *                   props(4,jd),props(18,jd),props(20,jd),ok)

         if (.not.ok.and.iopt(16).eq.0) shear = .false.  
c                                 explicit bulk modulus is allowed and used
         if (lopt(17).and.props(4,jd).gt.0d0) then
            bulk = .false.
            bulkg = .false.
         end if 

      else

         props(5,jd)  = 0d0
         props(19,jd) = 0d0    
         props(21,jd) = 0d0

      end if   
            
      g0 = ginc(0d0,0d0,id)
c                                 get speciation 
      if (id.gt.0) call getspc (id,jd)
c                                 set flag for multiple root eos's
      sroot = .true.
c                                 save derivative for cp search
c      vp(jd) = pv
c      vvp(jd) = pvv
c      rooti(jd) = iroot
c                                 compute g-derivatives for isostatic 
c                                 thermodynamic properties
      if (p.gt.nopt(26)) then 
         dp0 = nopt(27) * p
      else 
         dp0 = nopt(27) * nopt(26)
      end if 

c     if (fluid(jd)) dp0 = dp0/10d0

      dt0 = dt
c                                 if a reaction, cannot use sign to 
c                                 test behavior, just run through the
c                                 whole list
      if (rxn) then 
c                                 use sign of s and v and second derivatives
c                                 to refine difference increments
         call getdpt (g0,dp0,dp1,dp2,dt0,dt1,dt2,v,gpp,s,gtt,gpt,id,fow)

         e = g0 + t * s

         cp = -t*gtt
c                                  these are the only properties output
c                                  for reactions:
         alpha = gpt/v
         beta = -gpp/v
         props(2,jd) = e
         props(11,jd) = g0 
         props(12,jd) = cp
         props(15,jd) = s
         props(1,jd) = v
         props(13,jd) = alpha
         props(14,jd) = beta 

      else 
c                                 real phase, use sign of s and v and second derivatives
c                                 to refine difference increments
         call getdpt (g0,dp0,dp1,dp2,dt0,dt1,dt2,v,gpp,s,gtt,gpt,id,fow)
c                                 enthalpy
         e = g0 + t * s
c                                 heat capacity
         cp = -t*gtt
c                                 volumetric properties only if v is ok:
         if (v.gt.0d0) then 
c                                 third order derivatives, only need for
c                                 derivatives of seismic props.
            if (fow) then 

               g1 = ginc(-dt2,0d0,id)
               g2 = ginc( dt2,2d0*dp2,id)
               g3 = ginc( dt2,0d0,id)
               g4 = ginc(-dt2,2d0*dp2,id)
               g5 = ginc(0d0,dp2,id)
               g7 = g1 - g3

               gppp = ((ginc(0d0,4d0*dp2,id) - g0)/2d0
     *                 - ginc(0d0,3d0*dp2,id) + g5)/dp2**3
               gppt = (g2 + g3 + 2d0*(ginc(-dt2, dp2,id)
     *                -ginc(dt2,dp2,id)) -g4 - g1)/dp2/dp2/dt2/2d0
               gptt = (g2 + g4  + 2d0*(g0 - ginc(0d0,2d0*dp2,id))
     *                 -g3 - g1)/2d0/dp2/dt2/dt2

            else 

               g1 = ginc(-dt2,-dp2,id) - ginc( dt2,dp2,id)
               g3 = ginc( dt2,-dp2,id)
               g4 = ginc(-dt2,dp2,id)
               g5 = ginc(0d0,dp2,id) - ginc(0d0,-dp2,id)
               g7 = ginc(-dt2,0d0,id) - ginc(dt2, 0d0,id)
   
               gppp = ((ginc(0d0,2d0*dp2,id) 
     *                - ginc(0d0,-2d0*dp2,id))/2d0  - g5)/dp2**3
               gppt = (g3 - g4 + 2d0*g7 - g1)/dp2/dp2/dt2/2d0
               gptt = (g4 - g3 - 2d0*g5 - g1)/2d0/dp2/dt2/dt2

            end if 

            gttt = ((ginc(dt2*2d0,0d0,id) - ginc(-dt2*2d0,0d0,id))/2d0 
     *              + g7)/dt2**3

            g7 = (gtt*gpp-gpt**2)**2

            if (g7.ne.0d0.and.bulk) then 
c                                 temperature derivative of the adiabatic bulk modulus:
               props(18,jd) = (((v*gppt-gpt*gpp)*gtt
     *                         -(2d0*v*gptt-gpt**2)*gpt)*gtt
     *                         +v*gttt*gpt**2)/g7
c                                 pressure derivative of the adiabatic bulk modulus:
               props(20,jd) = (((v*gppp-gpp**2)*gtt
     *                      +(gpt*gpp-2d0*v*gppt)*gpt)*gtt
     *                      +v*gptt*gpt**2)/g7
            else if (bulk) then 

               props(18,jd) = nopt(7)
               props(20,jd) = nopt(7)
               sick(jd) = .true.

            end if 

         else 
c                                 no volume, zero second derivatives
            gpp = 0d0
            gpt = 0d0    
      
         end if 

         if (v.le.0d0) then 

            if (v.lt.0d0.or.(v.eq.0d0.and.iam.ne.5)) then 
c                                 allow zero volume for 1 bar gas reference state 
c                                 in frendly
               sick(jd) = .true.
               v = nopt(7)
               beta = nopt(7)
               alpha = nopt(7)
               rho = nopt(7)
         
            else 

               beta = 0d0
               alpha = 0d0
               rho = 0d0 
      
            end if 

         else 

            beta = -gpp/v
            alpha = gpt/v
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

         if (gpp.eq.0d0.or.gtt.eq.0d0.or.gpt.eq.0d0) sick(jd) = .true.

         props(2,jd) = e
         props(11,jd) = g0 
         props(12,jd) = cp
         props(15,jd) = s
         props(1,jd) = v
         props(13,jd) = alpha
         props(14,jd) = beta 
         props(10,jd) = rho  

      end if 

      if (.not.sick(jd).and..not.rxn.and.v.gt.0d0) then
c                                 heat capacity ratio (cp/cv)
         props(28,jd) = 1d0/(1d0 - gpt**2/gpp/gtt)
c                                 gruneisen parameter
         props(3,jd) = v/t/(gtt*gpp/gpt-gpt)
c                                 aug 28, 2007, removed check on gruneisen to 
c                                 accomodate -alpha's generated by landau 
c                                 transition models (prop(13,jd))
c        if (props(3,jd).le.0d0) sick(jd) = .true.

c                                 adiabatic bulk modulus
         if (bulk) props(4,jd) = -v/(gpp-gpt**2/gtt)

         if (props(4,jd).le.0d0) sick(jd) = .true.

         if (.not.fluid(jd).and.iopt(16).gt.0) then 
c                                 use poisson ratio estimates if iopt(16).ne.0
            if ((iopt(16).eq.1.and..not.ok).or.iopt(16).eq.2) then

               if (.not.sick(jd).or..not.bulk) then
 
                  props(5,jd)  = nopt(16)*props(4,jd)

                  if (.not.sick(jd)) then 
                     props(19,jd) = nopt(16)*props(18,jd) 
                     props(21,jd) = nopt(16)*props(20,jd)
                  else 
                     props(19,jd) = nopt(7)
                     props(21,jd) = nopt(7)
                  end if 

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

      if (sick(jd).and.bulk) then

         props(3,jd) = nopt(7)
         props(4,jd) = nopt(7)

         volume = .false.

         if (.not.fluid(jd).and..not.ok) shear = .false.

         if (.not.fluid(jd)) ssick = .true.

      else if (.not.bulk.and.rho.le.0d0.or.props(4,jd).lt.0d0) then 

         volume = .false.

      else if (.not.bulk.and.sick(jd)) then

         bsick = .true. 

      end if 
c                                 at this point sick(jd) is true if:
c                                 1) v < 0 and not a reaction or an ideal gas in frendly
c                                 2) alpha and beta are undefined
c                                 3) cp or beta is < 0 or unreasonably large

c                                 bulk is true if bulk modulus is computed thermodynamically
c                                 volume is false only if sick and .not.bulk

c                                 seismic properties
      if (volume) then

         if (props(5,jd).gt.0d0.or.fluid(jd)) lshear = .true.
c                                 sound velocity
         root = props(4,jd)/rho

         props(6,jd) = dsqrt(root) * units

         if (.not.sick(jd)) then 
c                                 sound velocity T derivative
            props(22,jd) = (props(18,jd) + props(4,jd) * alpha) 
     *                     / dsqrt(root) / rho / 2d0 * units
c                                 sound velocity P derivative
            props(25,jd) = (props(20,jd) - props(4,jd) * beta) 
     *                  / dsqrt(root) / rho / 2d0 * units
         else 
            props(22,jd) = nopt(7)
            props(25,jd) = nopt(7)
         end if 

         if (lshear) then
c                                 p-wave velocity
            root = (props(4,jd)+r43*props(5,jd))/rho 

            props(7,jd) = dsqrt(root)*units

            if (.not.sick(jd)) then 
c                                 p-wave velocity T derivative
               props(23,jd) = (props(18,jd) + r43 * 
     *                        (props(19,jd) + alpha * props(5,jd)) 
     *                         + props(4,jd) * alpha) / 
     *                         dsqrt(root) / rho / 2d0 * units
c                                 p-wave velocity P derivative
               props(26,jd) = (props(20,jd) + r43 *
     *                        (props(21,jd) - beta * props(5,jd)) 
     *                       - props(4,jd) * beta) /
     *                         dsqrt(root) / rho / 2d0 * units
            else 
               props(23,jd) = nopt(7)
               props(26,jd) = nopt(7)
            end if 
c                                 s-wave velocity
            root = props(5,jd)/rho

            props(8,jd) = dsqrt(root)*units

            if (fluid(jd)) then 

               props(24,jd) = 0d0
               props(27,jd) = 0d0

            else if (.not.sick(jd)) then

               props(24,jd) = (props(19,jd) + props(5,jd) * alpha)
     *                           / dsqrt(root) / rho / 2d0 * units
               props(27,jd) = (props(21,jd) - props(5,jd) * beta)
     *                           / dsqrt(root) / rho / 2d0 * units
            else 
               props(24,jd) = nopt(7)
               props(27,jd) = nopt(7)
            end if 
c                                 vp/vs
            if (isnan(props(8,jd))) then 
               props(9,jd) = nopt(7) 
            else if (props(8,jd).ne.0d0) then 
               props(9,jd) = props(7,jd)/props(8,jd)
            else
               props(9,jd) = nopt(7)
            end if 
c                                 check for negative poisson ratio
            if (poiss(props(7,jd),props(8,jd)).lt.0d0) then

               if (iwarn3.lt.11.and.pname(jd).ne.wname3) then 

                 write (*,1050) t,p,pname(jd)

                 if (lopt(20)) then 
                    write (*,1070)
                 else 
                    write (*,1060)
                 end if 

                 iwarn3 = iwarn3 + 1
                 wname3 = pname(jd)
                 if (iwarn3.eq.11) call warn (49,r,180,'GETPHP') 

               end if

               if (lopt(20)) then 

                  sick(jd) = .true.
                  bsick = .true.
                  ssick = .true.
                  volume = .false.
                  shear = .false.

                  v = nopt(7)
                  beta = nopt(7)
                  alpha = nopt(7)
                  rho = nopt(7)

                  props(13,jd) = nopt(7)
                  props(14,jd) = nopt(7) 
                  props(10,jd) = nopt(7) 

                  props(5,jd) = nopt(7)
                  props(19,jd) = nopt(7)    
                  props(21,jd) = nopt(7)

                  props(1,jd) = nopt(7)
                  props(3,jd) = nopt(7)
                  props(4,jd) = nopt(7)

                  props(6,jd)  = nopt(7)
                  props(22,jd) = nopt(7)
                  props(25,jd) = nopt(7)

                  props(7,jd)  = nopt(7)
                  props(23,jd) = nopt(7)
                  props(26,jd) = nopt(7)

                  props(8,jd) = nopt(7)
                  props(24,jd) = nopt(7)
                  props(27,jd) = nopt(7)

                  props(28,jd) = nopt(7)
               end if 

            end if 

         else

            props(7,jd)  = nopt(7)
            props(23,jd) = nopt(7)
            props(26,jd) = nopt(7)

            props(8,jd) = nopt(7)
            props(24,jd) = nopt(7)
            props(27,jd) = nopt(7)

            shear = .false.

         end if  

      else 

         do j = 3, 9
            props(j,jd) = nopt(7)
         end do 

         do j = 18, 27
            props(j,jd) = nopt(7)
         end do 

      end if 
c                                 get min/max moduli for hashin-strikman
c                                 bounds. also saves the corresponding 
c                                 T and P derivatives.
      if (.not.lopt(16)) then

         do j = 1, 2
c                                 property index 
            m = hs2p(j)
            if (isnan(props(m,jd))) cycle 
c                                 min aggregate prop
            if (props(m,jd).lt.hsb(j,1)) then
               hsb(j,1) = props(m,jd)
               hsb(j+2,1) = props(hs2p(j+2),jd)
               hsb(j+4,1) = props(hs2p(j+4),jd)
            end if 
c                                 max aggregate prop
            if (props(m,jd).gt.hsb(j,2)) then
               hsb(j,2) = props(m,jd)
               hsb(j+2,2) = props(hs2p(j+2),jd)
               hsb(j+4,2) = props(hs2p(j+4),jd)
            end if 

            if (fluid(jd)) cycle 
c                                 min solid prop
            if (props(m,jd).lt.hsb(j,3)) then
               hsb(j,3) = props(m,jd)
               hsb(j+2,3) = props(hs2p(j+2),jd)
               hsb(j+4,3) = props(hs2p(j+4),jd)
            end if 
c                                 max solid prop
            if (props(m,jd).gt.hsb(j,4)) then
               hsb(j,4) = props(m,jd)
               hsb(j+2,4) = props(hs2p(j+2),jd)
               hsb(j+4,4) = props(hs2p(j+4),jd)
            end if 

         end do 

      end if 
c                                 check and warn if necessary for negative
c                                 expansivity
      if (.not.sick(jd).and.v.gt.0d0.and.alpha.le.0d0.and.iwarn1.lt.11
     *    .and.pname(jd).ne.wname1) then

         write (*,1030) t,p,pname(jd)
         iwarn1 = iwarn1 + 1
         wname1 = pname(jd)
         if (iwarn1.eq.11) call warn (49,r,179,'GETPHP') 

      end if

      if (ppois.and.iwarn2.lt.11.and.pname(jd).ne.wname2.and.pois) then 

         iwarn2 = iwarn2 + 1
         wname2 = pname(jd)
         write (*,1040) t,p,pname(jd)
         if (iwarn2.eq.11) call warn (49,r,178,'GETPHP')

      end if 
c                                 accumulate non-seismic totals 
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
c                                 gruneisen T
      psys(3)  = psys(3)  + props(3,jd)*v 
c                                 molar gibbs energy
      psys(11) = psys(11) + g0*mols 
c                                 molar heat capacity
      psys(12) = psys(12) + cp*mols 
c                                 expansivity
      psys(13) = psys(13) + alpha*v 
c                                 compressibility
      psys(14) = psys(14) + beta*v
c                                 molar entropy
      psys(15) = psys(15) + s*mols 
c                                 moles of assemblage
      psys(16) = psys(16) + mols
c                                 mass of assemblage 
      psys(17) = psys(17) + props(17,jd)*mols
c                                 solid only totals:
      if (.not.fluid(jd)) then 

         psys1(1)  = psys1(1)  + v
         psys1(2)  = psys1(2)  + e*mols 
         psys1(3)  = psys1(3)  + props(3,jd)*v 
         psys1(11) = psys1(11) + g0*mols 
         psys1(12) = psys1(12) + cp*mols 
         psys1(13) = psys1(13) + alpha*v 
         psys1(14) = psys1(14) + beta*v
         psys1(15) = psys1(15) + s*mols 
         psys1(16) = psys1(16) + mols
         psys1(17) = psys1(17) + props(17,jd)*mols

      end if 

      sroot = .false. 

1030  format (/,'**warning ver179** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        'the effective expansivity of: ',a,/,'is negative. ',
     *        'Most probably this is because of a Landau ordering ',
     *        'model. The Gruneisen',/,'thermal parameter and seismic',
     *        ' speeds for this phase should be considered ',
     *        'with caution.',/)
1040  format (/,'**warning ver178** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        /,'the shear modulus of: ',a,/,'is missing or invalid ',
     *        'and has been estimated from the default poisson ',
     *        'ratio',/)
1050  format (/,'**warning ver180** at T(K)=',g12.4,' P(bar)=',g12.4,1x,
     *        /,'the poisson ratio of: ',a,' is negative.')
1060  format (/,'the suspect data will be output and used to ',
     *        'compute aggregate properties.',/,'to prevent this set ',
     *        'poisson_check to true in perplex_option.dat',/)
1070  format (/,'the suspect data will be rejected and aggregate ',
     *        'properties will not be output.',/,'to prevent this set',
     *        'poisson_check to false in perplex_option.dat',/)

      end

      subroutine getdpt (g0,dp0,dp1,dp2,dt0,dt1,dt2,v,gpp,s,gtt,gpt,
     *                   id,fow)
c----------------------------------------------------------------------
c getdpt computes finite difference increments for phase id based on finite
c difference estimates of v, gpp (dv/dp), s, gtt (ds/dt) and initial guesses 
c for dp0.  returns:
c    v, gpp (dv/dp), s, gtt (ds/dt)  - finite difference estimates
c    dp0, dp1, dp2 - increments for 1st, 2nd and 3rd order p finite differences
c    dt0, dt1, dt2 - increments for 1st, 2nd and 3rd order p finite differences
c    fow - if true, forward differences are needed on p, otherwise centered.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical fow, okt

      integer i, j, id

      double precision g0, dp0, dp1, dp2, v, gpp, ginc, gpt, gpt1, gpt2,
     *                 s, gtt, dt0, dt1, dt2, xdp, xdt

      external ginc

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
c                                 pressure increments
      okt = .false.

      xdp = dp0
      xdt = dt0

      do i = 1, 2

         do j = 1, 3 

            call getgpp (g0,dp0,dp1,dp2,v,gpp,id,fow)

            if (v.gt.0d0.and.gpp.lt.0d0.and.gpp.gt.-v) then
               okt = .true.
               exit
            end if 

            if (i.eq.1) then 
               dp0 = dp0 * nopt(31)
            else 
               dp0 = dp0 / nopt(31)
            end if 
 
         end do 

         if (okt) exit

         dp0 = dp0/nopt(31)**4

      end do

      if (okt) then 
         dp0 = dabs(1d-5*v/gpp)
      else 
c                                 negative compressibility?
         dp0 = xdp
      end if 
c                                 final values
      call getgpp (g0,dp0,dp1,dp2,v,gpp,id,fow)
c                                 -------------------------------------
c                                 temperature increments
      okt = .false.

      do i = 1, 2

         do j = 1, 3

            call getgtt (g0,dt0,dt1,dt2,s,gtt,id)

            if (s.gt.0d0.and.gtt.lt.0) then
                 okt = .true.
                 exit 
            end if 

            if (i.eq.1) then 
               dt0 = dt0 * nopt(31)
               if (dt0.gt.t) exit 
            else
                dt0 = dt0 / nopt(31)
            end if 

         end do 

         if (okt) exit

         dt0 = dt0/nopt(31)**4

      end do

      if (okt) then 
         dt0 = dabs(1d-5*s/gtt)
      else 
c                                 something has gone horribly wrong! 
         dt0 = xdt
      end if 
c                                 final values
      call getgtt (g0,dt0,dt1,dt2,s,gtt,id)

      if (v.gt.0d0) then 
c                                 get the cross derivative gpt for expansivity
         if (fow) then 

            gpt = ( ginc( dt1,dp1,id) - ginc( dt1,0d0,id)
     *             -ginc(-dt1,dp1,id) + ginc(-dt1,0d0,id))/dt1/dp1/2d0

            if (gpt.gt.v.or.gpt.le.0d0)
c                                 expand increment if invalid alpha
     *         gpt = ( ginc( dt2,dp2,id) - ginc( dt2,0d0,id)
     *                -ginc(-dt2,dp2,id) + ginc(-dt2,0d0,id))
     *                /dp2/dt2/2d0

               if (gpt.gt.v.or.gpt.le.0d0)
c                                 shrink increment if invalid alpha
     *            gpt = ( ginc( dt0,dp0,id) - ginc( dt0,0d0,id)
     *                   -ginc(-dt0,dp0,id) + ginc(-dt0,0d0,id))
     *                  /dp0/dt0/2d0

         else

            gpt = ( ginc( dt1,dp1,id) - ginc( dt1,-dp1,id)
     *             -ginc(-dt1,dp1,id) + ginc(-dt1,-dp1,id))/dt1/dp1/4d0

            if (gpt.gt.v.or.gpt.le.0d0) then 
c                                 expand increment if invalid alpha
               gpt1 = ( ginc( dt2,dp2,id) - ginc( dt2,-dp2,id)
     *                - ginc(-dt2,dp2,id) 
     *                     + ginc(-dt2,-dp2,id))/dp2/dt2/4d0

               if (gpt1.gt.v.or.gpt1.le.0d0) then
c                                 shrink increment if invalid alpha
                  gpt2 = ( ginc( dt0,dp0,id) - ginc( dt0,-dp0,id)
     *                    -ginc(-dt0,dp0,id) + ginc(-dt0,-dp0,id))
     *                     /dp0/dt0/4d0

                  if (gpt2.lt.v.and.gpt2.ge.0d0) then 
                     gpt = gpt2
                  end if 

               else 

                  gpt = gpt1

               end if 

            end if  

         end if  

      end if 

      end 

      subroutine getgpp (g0,dp0,dp1,dp2,v,gpp,id,fow)
c----------------------------------------------------------------------
c for phase id getgpp computes:
c    v, gpp (dv/dp) by finite difference
c    dp1, dp2 - increments for 2nd and 3rd order p finite differences
c    fow - if true, forward differences needed, otherwise centered.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical fow

      integer id 

      double precision g0,dp0,dp1,dp2,v,gpp, ginc

      external ginc

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      dp1 = dp0 * nopt(31)
      dp2 = dp1 * nopt(31)

      if (p-dp2.le.0d0) then
         fow = .true.
      else 
         fow = .false.
      end if  

      if (fow) then 

         v = (ginc(0d0,dp0,id) - g0)/dp0
         gpp = (ginc(0d0,2d0*dp1,id) + g0 
     *                  - 2d0*ginc(0d0,dp1,id))/dp1/dp1

      else             

         v = (ginc(0d0,dp0,id) - ginc(0d0,-dp0,id))/dp0/2d0
         gpp = (ginc(0d0,dp1,id) + ginc(0d0,-dp1,id) - 2d0*g0)/dp1/dp1

      end if 

      end 

      subroutine getgtt (g0,dt0,dt1,dt2,s,gtt,id)
c----------------------------------------------------------------------
c for phase id getgtt computes:
c    s, gtt (ds/dt) by finite difference
c    dt1, dt2 - increments for 2nd and 3rd order t finite differences
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id 

      double precision g0,dt0,dt1,dt2,s,gtt, ginc

      external ginc

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------

      if (nopt(31)**2*dt0.ge.t) dt0 = t/nopt(31)**2*dt0*0.9d0

      dt1 = dt0 * nopt(31)
      dt2 = dt1 * nopt(31)

      s = (ginc(-dt0,0d0,id) - ginc(dt0,0d0,id))/dt0/2d0
      gtt = (ginc(dt1,0d0,id) + ginc(-dt1,0d0,id) - 2d0*g0)/dt1/dt1

      end 

      subroutine gtsysp (sick,ssick,bulkg,bsick)
c-----------------------------------------------------------------------
c computes aggregate (system) properties from phase properties 
c obtained by prior calls to getphp
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical sick(i8), ssick, solid, bad, bulkg, bsick

      integer i, j, iwarn, m

      double precision chi, chi1, root, k, g, v, vs

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      character pname*14
      common/ cxt21a /pname(k5)

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer iam
      common/ cst4 /iam

      double precision units, r13, r23, r43, r59, r1, r2
      common/ cst59 /units, r13, r23, r43, r59, r1, r2

      save iwarn
      data iwarn/0/
c----------------------------------------------------------------------
c                                 check if volume is there, if not assume
c                                 things are really bad
      if (rxn) then 
c                                 frendly return if a reaction
         shear = .false.
         return 

      else if (isnan(psys(1))) then 
     
         bad = .true.
    
      else if (iam.eq.5.and.psys(1).eq.0d0) then 
c                                 frendly but not a reaction
         shear = .false.
         return         
 
      else if (psys(1).lt.0d0) then 

         bad = .true.

      else 

         bad = .false. 

      end if 

      if (bad) then 

         do i = 1, i8
            psys(i) = nopt(7)
            psys1(i) = nopt(7)
         end do

         return

      end if 

c     if (.not.volume) shear = .false.
c                                 not so bad....
      solid = .true.
c                                 weighting used to compute average
c                                 moduli (and derivatives) from bounds:
c                                 chi = 1 -> fast bound 
c                                 chi = 0 -> slow bound 
      chi = nopt(6)
      chi1 = 1d0 - chi
c                                 hashin-shtrikman limiting values, set
c                                 in calphp. ony 4,5,18,19,20,21 are used
c                                 1 - min, 2 - max, 3 - min solid, 4 - max
c                                 solid. 
      if (.not.lopt(16)) then 

         do m = 1, 5, 2
            do i = 1, 4

               k = hsb(m,i)
               g = hsb(m+1,i)
               hsb(m,i)   = r43*g

               if (k.eq.0d0.and.g.eq.0d0) then 
                  hsb(m+1,i) = 0d0
               else 
                  hsb(m+1,i) = g*((9d0*k+8d0*g)/(k+2d0*g))/6d0
               end if 

            end do 
         end do 

      end if 
c                                 compute aggregate props:

c                                 convert volumetrically weighted totals
c                                 to arithmetic means (for aseismic properties)
      do i = 13, 14
c                                 normalize volumetrically weighted alpha/beta
         psys(i) = psys(i)/psys(1)
         if (psys1(1).ne.0d0) psys1(i) = psys1(i)/psys1(1)

      end do 
c                                 heat capacity ratio (cp/cv)
      psys(28) = 1d0/(1d0-t*psys(1)*psys(13)**2/psys(14)/psys(12))
c                                 density, kg/m3
      psys(10) = psys(17)/psys(1)*1d2
c                                 gruneisen T
      psys(3) = psys(3)/psys(1)

      if (psys1(1).gt.0d0) then
c                                 if the system is not entirely fluid,
c                                 solid aggregate props, cp/cv 
         psys1(28) = 1d0 / 
     *               (1d0-t*psys1(1)*psys1(13)**2/psys1(14)/psys1(12)) 

         if (.not.isnan(psys1(3))) then 
c                                 gruneisen T
            psys1(3) = psys1(3)/psys1(1)
         else 
            psys1(3) = nopt(7)
         end if 
c                                 density
         psys1(10) = psys1(17)/psys1(1)*1d2

      end if 
c                                 if a reaction (frendly) return
      if (rxn) return
c                                 if not solid, don't compute solid only
c                                 properties
      if (psys1(1).le.0d0.or..not.aflu) then 
         solid = .false.
         aflu = .false.
      end if 
c                                 accumulate aggregate totals, some
c                                 totals may be incomplete if volume or 
c                                 shear is false for an individual phase
c                                 weighting scheme for seismic velocity
      do i = 1, ntot 
c                                 total volume fraction
         v = props(1,i)*props(16,i)/psys(1)
c                                 solid volume fraction
         if (.not.fluid(i).and.aflu) 
     *                    vs = props(1,i)*props(16,i)/psys1(1)
c                                 for elastic properties use
c                                 VRH if lopt(16), else HS
         if (lopt(16)) then 

            if (volume) then 

               do j = 1, 5, 2

                  m = hs2p(j)

                  if (j.gt.1.and.bsick .or. props(m,i).eq.0d0) cycle 
c                                 Aggregate Bulk Modulus, T-derivative, P-derivative                                
                  psys(m) = psys(m) + v*props(m,i)
                  pgeo(m) = pgeo(m) + v/props(m,i)

               end do 

            end if 
c                                 aggregate shear props only if
c                                 shear mod is available for all phases.
            if (shear) then

               do j = 2, 6, 2

                  m = hs2p(j) 
                  if (j.gt.2.and.bsick .or. props(m,i).eq.0d0) cycle 
c                                 Aggregate shear Modulus, T-derivative, P-derivative                                
                  psys(m) = psys(m) + v*props(m,i)
                  if (.not.aflu) pgeo(m) = pgeo(m) + v/props(m,i)

               end do 

            end if 

            if (aflu.and..not.fluid(i).and.(bulkg.and..not.ssick .or.
     *                                                .not.bulkg)) then
c                                 assemblage includes fluid, solid only
c                                 totals:
               do j = 1, 5, 2

                  m = hs2p(j)
                  if (j.gt.1.and.bsick .or. props(m,i).eq.0d0) cycle 
c                                 Aggregate Bulk Modulus, T-derivative, P-derivative                                
                  psys1(m) = psys1(m) + vs*props(m,i)
                  pgeo1(m) = pgeo1(m) + vs/props(m,i)

               end do 

            end if 

            if (aflu.and..not.fluid(i).and.shear) then
               
               do j = 2, 6, 2

                  m = hs2p(j)
                  if (j.gt.2.and.bsick .or. props(m,i).eq.0d0) cycle 
c                                 Aggregate shear Modulus, T-derivative, P-derivative                                
                  psys1(m) = psys1(m) + vs*props(m,i)
                  pgeo1(m) = pgeo1(m) + vs/props(m,i)

               end do 

            end if 

         else 
c                                 HS sums: 
c                                 psys is used for upper bound
c                                 pgeo used for lower bound
            if (volume) then 
c                                 Aggregate bulk modulus, T-derivative, P-derivative 
               do j = 1, 5, 2

                  m = hs2p(j)
                  if (j.gt.1.and.bsick .or. props(m,i).eq.0d0) cycle 
                  psys(m) = psys(m) + v/(props(m,i)+hsb(j,2))
                  pgeo(m) = pgeo(m) + v/(props(m,i)+hsb(j,1))

               end do 

            end if 

            if (shear) then 
c                                 Aggregate shear modulus, T-derivative, P-derivative 
               do j = 2, 6, 2

                  m = hs2p(j) 

                  if (j.gt.2.and.bsick .or. props(m,i).eq.0d0) cycle 
                  psys(m) = psys(m) + v/(props(m,i)+hsb(j,2))
                  if (.not.aflu) pgeo(m) 
     *                              = pgeo(m) + v/(props(m,i)+hsb(j,1))

               end do 

            end if 

            if (aflu.and..not.fluid(i).and..not.ssick) then
c                                 assemblage includes fluid, solid only totals:
c                                 Aggregate bulk modulus, T-derivative, P-derivative 
               do j = 1, 5, 2

                  m = hs2p(j)

                  if (j.gt.1.and.bsick .or. props(m,i).eq.0d0) cycle 
                  psys1(m) = psys1(m) + vs/(props(m,i)+hsb(j,4))
                  pgeo1(m) = pgeo1(m) + vs/(props(m,i)+hsb(j,3))

               end do 

            end if 

            if (aflu.and..not.fluid(i).and.shear) then
c                                 assemblage includes fluid, solid only totals:
c                                 Aggregate shear modulus, T-derivative, P-derivative 
               do j = 2, 6, 2

                  m = hs2p(j)

                  if (j.gt.2.and.bsick .or. props(m,i).eq.0d0) cycle 
                  psys1(m) = psys1(m) + vs/(props(m,i)+hsb(j,4))
                  pgeo1(m) = pgeo1(m) + vs/(props(m,i)+hsb(j,3))

               end do 

            end if 

         end if 

      end do 
c                                 seismic moduli and derivatives
      do j = 1, 6
c                                 property index
         m = hs2p(j)
c                                 combine as VRH or HS means
         if (lopt(16)) then 
c                                 VRH
            psys(m) = chi*psys(m)
            if (pgeo(m).ne.0d0) psys(m) = psys(m)+chi1/pgeo(m) 
            psys1(m) = chi*psys1(m)
            if (pgeo1(m).ne.0d0) psys1(m) = psys1(m)+chi1/pgeo1(m) 

         else 
c                                 HS 
            if (psys(m).ne.0d0) psys(m) = 1d0/psys(m) 
            psys(m) = chi*(psys(m) - hsb(j,2))

            if (pgeo(m).ne.0d0) pgeo(m) = 1d0/pgeo(m)
            psys(m) = psys(m) + chi1*(pgeo(m) - hsb(j,1))

            if (psys1(m).ne.0d0) psys1(m) = 1d0/psys1(m) 
            psys1(m) = chi*(psys1(m) - hsb(j,4))

            if (pgeo1(m).ne.0d0) pgeo1(m) = 1d0/pgeo1(m)
            psys1(m) = psys1(m) + chi1*(pgeo1(m) - hsb(j,3))

         end if 
 
      end do 
c                                 ----------------------------------
c                                 aggregate velocities
      if (volume) then
c                                 sound velocity
         root = psys(4)/psys(10)

         if (root.gt.0d0) then 

            psys(6) = dsqrt(root) * units

            if (.not.bsick) then 
c                                 sound velocity T derivative
               psys(22) = (psys(18) + psys(4) * psys(13)) 
     *                     / dsqrt(root) / psys(10) / 2d0 * units
c                                 sound velocity P derivative
               psys(25) = (psys(20) - psys(4) * psys(14)) 
     *                     / dsqrt(root) / psys(10) / 2d0 * units
            else 
               psys(22) = nopt(7)
               psys(25) = nopt(7) 
            end if 
         end if 

         if (shear) then 
c                                 s-wave velocity
            root = psys(5)/psys(10)

            if (root.gt.0d0) then 

               psys(8) = dsqrt(root) * units

               if (.not.bsick) then 
c                                 T-derivative
                  psys(24) = (psys(19) + psys(5) * psys(13)) 
     *                          / dsqrt(root) / psys(10) / 2d0 * units
c                                 P-derivative 
                  psys(27) = (psys(21) - psys(5) * psys(14)) 
     *                          / dsqrt(root) / psys(10) / 2d0 * units
               else 
                  psys(24) = nopt(7)
                  psys(27) = nopt(7) 
               end if 
            end if 

            root = (psys(4)+r43*psys(5))/psys(10)

            if (root.gt.0d0) then 
c                                 p-wave velocity
               psys(7) = dsqrt(root)*units

               if (.not.bsick) then 
c                                 p-wave velocity T derivative
                  psys(23) = (psys(18) + r43*(psys(19) 
     *                       + psys(13) * psys(5)) 
     *                       + psys(4) * psys(13)) / 
     *                         dsqrt(root) / psys(10) / 2d0 * units
c                                 p-wave velocity P derivative
                  psys(26) = (psys(20) + r43*(psys(21) 
     *                       - psys(14) * psys(5)) 
     *                       - psys(4) * psys(14)) /
     *                         dsqrt(root) / psys(10) / 2d0 * units
               else 
                  psys(23) = nopt(7)
                  psys(26) = nopt(7) 
               end if 
            end if 
c                                 vp/vs
            if (psys(8).gt.0d0) then 
               psys(9) = psys(7)/psys(8)
            else
               psys(9) = nopt(7)
            end if
         else 

            psys(7)  = nopt(7)
            psys(8)  = nopt(7)
            psys(9)  = nopt(7)
            psys(23) = nopt(7)
            psys(24) = nopt(7)
            psys(26) = nopt(7)
            psys(27) = nopt(7)

         end if 

      else 
c                                 set missing data
         do j = 3, 9
            psys(j) = nopt(7)
         end do 

         do j = 18, 27
            psys(j) = nopt(7)
         end do 

      end if 
c                                 ----------------------------------
c                                 fluid-absent properties:
c                                 the psys1(1) condition is for the 
c                                 special case of a system consisting 
c                                 only of fluid. 
      if (solid.and.(.not.ssick.and.bulkg).or..not.bulkg) then 
c                                 fluid absent properties:
         root = psys1(4)/psys1(10)

         if (root.gt.0d0) then 
c                                 sound velocity
            psys1(6) = dsqrt(root) * units

            if (.not.bsick) then 
c                                 sound velocity T derivative
               psys1(22) = (psys1(18) + psys1(4) * psys1(13)) 
     *                     / dsqrt(root) / psys(10) / 2d0 * units
c                                 sound velocity P derivative
               psys1(25) = (psys1(20) - psys1(4) * psys1(14)) 
     *                     / dsqrt(root) / psys1(10) / 2d0 * units
            else 
               psys1(22) = nopt(7)
               psys1(25) = nopt(7) 
            end if 
         end if

         if (shear.and.solid) then 

            root = psys1(5)/psys1(10)

            if (root.gt.0d0) then 
c                                 s-wave velocity
               psys1(8) = dsqrt(root) * units

               if (.not.bsick) then 
c                                 T-derivative
                  psys1(24) = (psys1(19) + psys1(5) * psys1(13)) 
     *                       / dsqrt(root) / psys1(10) / 2d0 * units
c                                 P-derivative 
                  psys1(27) = (psys1(21) - psys1(5) * psys1(14)) 
     *                       / dsqrt(root) / psys1(10) / 2d0 * units
               else
                  psys1(24) = nopt(7)
                  psys1(27) = nopt(7) 
               end if 
            end if 

            root = (psys1(4)+r43*psys1(5))/psys1(10)

            if (root.gt.0d0) then 
c                                 p-wave velocity
               psys1(7) = dsqrt(root)*units

               if (.not.bsick) then 
c                                 p-wave velocity T derivative
                  psys1(23) = (psys1(18) + r43*(psys1(19) 
     *                       + psys1(13) * psys1(5)) 
     *                       + psys1(4) * psys1(13)) / 
     *                       dsqrt(root) / psys1(10) / 2d0 * units
c                                 p-wave velocity P derivative
                  psys1(26) = (psys1(20) + r43*(psys1(21)
     *                       - psys1(14) * psys1(5)) 
     *                       - psys1(4) * psys1(14)) /
     *                       dsqrt(root) / psys1(10) / 2d0 * units
               else 
                  psys1(23) = nopt(7)
                  psys1(26) = nopt(7) 
               end if 
            end if 
c                                 vp/vs
            if (psys1(8).gt.0d0) then 
               psys1(9) = psys1(7)/psys1(8)
            else
               psys1(9) = nopt(7)
            end if 
         else

         end if 
      else 
c                                 set missing data
         do j = 3, 9
            psys1(j) = nopt(7)
         end do 

         do j = 18, 27
            psys1(j) = nopt(7)
         end do 
         
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

      subroutine insysp (ssick,ppois,bulkg,bsick)
c-----------------------------------------------------------------------
c insysp initializes system properties accumulated in getphp and gtsysp
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      logical ssick,ppois,bulkg,bsick

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer hs2p
      double precision hsb
      common/ cst84 /hsb(i8,4),hs2p(6)
c----------------------------------------------------------------------
c                                 assemblage flags
c                                 ----------------
c                                 aflu   = T if a fluid is present
      aflu   = .false.
c                                 shear  = T if shear modulus can be computed
      shear  = .true.
c                                 volume = T if volume and bulk modulus can be computed
      volume = .true.
c                                 ssick  = T if one or more phases has problems
      ssick  = .false.
c                                 ppois  = T if a shear modulus has been estimated from poisson ratio
      ppois  = .false.
c                                 rxn    = T if a reaction (frendly) => no associated mass
      rxn    = .false.
c                                 bulkg  = F if a explicit bulk modulus function has been used 
      bulkg  = .true.
c                                 bsick  = T if bulkg and volume and shear and ssick
      bsick  = .false.
c                                 initialize sums
      do i = 1, i8
         psys(i) = 0d0
         psys1(i) = 0d0
         pgeo(i) = 0d0 
         pgeo1(i) = 0d0
      end do 
c                                 initialize bulk properites
c                                 total mass
      gtot = 0d0
      gtot1 = 0d0

c                                 HS limiting moduli
      do i = 1, 6
         hsb(i,1) = 1d99
         hsb(i,2) = 0d0         
         hsb(i,3) = 1d99
         hsb(i,4) = 0d0 
      end do     

      do i = 1, icomp
c                                 total molar amounts
         fbulk(i) = 0d0
         fbulk1(i) = 0d0

      end do

      end 