c this is a main program to illustrate how meemum can be called as a subroutine

      program MAIN 

      implicit none

      include 'perplex_parameters.h'
c                                 some local variables:
      logical bad

      integer i, j, l, ier

      character cprop*18, amount*6
c----------------------------------------------------------------------
c                                 these common blocks are necessary to 
c                                 communicate between the MAIN program
c                                 and perplex:

c                                 jbulk -> the number of chemical components
c                                 cblk -> the molar bulk composition
      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 iam -> a variable indicating the Perple_X program
      integer iam
      common/ cst4 /iam
c                                 perplex option values:
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                 ntot - number of phases stable
c                                 np - number of solution phases
c                                 ncpd - number of compounds
c                                 kkp(ntot) - pointer to cpd or solution model
c                                 cp3(1:jbulk,1:ntot) - molar phase compositions
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 phase and system properties loaded by 
c                                 subroutine getloc, the call to getloc
c                                 is unnecessary if only the molar phase
c                                 proportions and compositions are of interest.
c                                 phase (prop(i,1:ntot)) and system (psys(i))
c                                 are for index i (this list may not be exhaustive):
c                                  1  - molar volume
c                                  2  - molar enthalpy
c                                  3  - gruneisen thermal parm
c                                  4  - K_S
c                                  5  - Mu_S
c                                  6  - v_phi
c                                  7  - v_p
c                                  8  - v_s
c                                  9  - v_p/v_s
c                                  10 - rho
c                                  11 - G
c                                  12 - cp
c                                  13 - alpha
c                                  14 - beta
c                                  15 - S
c                                  16 - molar amount
c                                  17 - molar weight
c                                  18 - KS_T
c                                  19 - MuS_T
c                                  20 - KS_P
c                                  21 - MuS_P
c                                  22 - vphi_T
c                                  23 - vp_T
c                                  24 - vs_T
c                                  25 - vphi_P
c                                  26 - vs_P
c                                  27 - vp_P
c                                  28 - heat capacity ratio (cp/cv)
      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)
c                                 atwt -> g-formula wts of the chemical components
      double precision atwt
      common/ cst45 /atwt(k0)
c                                 v -> the standard perplex potential variables
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c                                 ipot -> the number of potential variables in use
c                                 jv -> the indices of the potential variables
      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c                                 vname -> the name of the potential variables
      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
c                                 cname -> the names of the components
      character*5 cname
      common/ csta4 /cname(k5)
c                                 phase compositions
      double precision pcomp
      common/ cst324 /pcomp(k0,k5)
c                                 phase names
      character pname*14
      common/ cxt21a /pname(k5)
c                                 iwt => 1 mass bulk comp, 0 molar bulk comp
      integer iwt
      common/ cst209 /iwt
c----------------------------------------------------------------------
      iam = 2
c                                 initialization, read files etc. 
      call iniprp
c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      if (iwt.eq.1) then 
         amount = 'weight'
      else 
         amount = 'molar '
      end if
c                                 computational loop
      do 
c                                 get physical conditions
         write (*,'(10a)') 'Enter ',(vname(jv(i)), i = 1, ipot)
         read (*,*,iostat=ier) (v(jv(i)), i = 1, ipot)
c                                 get chemical composition
         do 
             write (*,1060) amount
             write (*,'(12(a,1x))') (cname(i),i=1,jbulk)
             read (*,*,iostat=ier) (cblk(i),i=1,jbulk)
             if (ier.eq.0) exit
         end do  
         
         if (iwt.eq.1) then 
c                                 convert mass to molar composition
            do i = 1, jbulk
               cblk(i) = cblk(i)/atwt(i)
            end do 

         end if
c                                 do the optimization
         call meemum (bad)

         if (bad) then
            write (*,*) 'optimization failed'
            cycle
         end if 
c                                 do some output
         if (iopt(2).eq.0) then 
            cprop = 'molar  proportions'
         else
            cprop = 'weight percentages'
         end if

         write (*,1020) cprop, (cname(i), i = 1, jbulk)
c                           phase proportions and compositions
         do i = 1, ntot

            write (*,'(1x,a,3x,3(f6.2,4x),g9.3,1x,20(f8.5,1x))') 
     *                      pname(i), 
c                                 weight %
     *                      props(17,i)*props(16,i)/psys(17)*1d2,
c                                 vol %
     *                      props(1,i)*props(16,i)/psys(1)*1d2,
c                                 mol %
     *                      props(16,i)/psys(16)*1d2,
c                                 mol
     *                      props(16,i),
c                                 wieght or molar composition
     *                      (pcomp(l,i), l = 1, jbulk)
         end do
c                           phase/system elementary properties
         do i = 1, ntot
c                                 N, G, S, V, Cp, alpha, beta, density
            write (*,1170) pname(i),props(17,i),int(props(11,i)),
     *                      props(15,i),
     *                      props(1,i),(props(j,i),j=12,14),props(28,i),
     *                      props(10,i)
         end do

         write (*,1170) 'System        ',psys(17),int(psys(11)),
     *                                   psys(15),psys(1),
     *                               (psys(j),j=12,14),psys(28),psys(10)

      end do 

1020  format (/,'Phase Compositions (',a,'):',
     *        /,19x,'wt %',6x,'vol %',5x,'mol %',5x,'mol  ',
     *          5x,20(1x,a,3x))
1060  format (/,'Enter ',a,' amounts of the components:')
1160  format (/,'Molar Properties and Density:'
     *        /,20x,'N(g)',10x,'G(J)',5x,'S(J/K)',5x,'V(J/bar)',6x,
     *         'Cp(J/K)',7x,'Alpha(1/K)',2x,'Beta(1/bar)',4x,'Cp/Cv',4x,
     *         'Density(kg/m3)')
1170  format (1x,a,1x,f9.2,3x,i12,1x,12(g12.5,1x),3x,f7.4)

      end


      subroutine meemum (bad)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, idead

      logical nodata, bad

      integer itri(4),jtri(4),ijpt

      double precision wt(3)

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 initialization
      rxn = .false.
c                                 normalize the composition vector, this 
c                                 is necessary for reasons of stupidity (lpopt0). 
      ctotal = 0d0

      do i = 1, icp
          ctotal = ctotal + cblk(i)
      end do 

      do i = 1, icp
         b(i) = cblk(i)/ctotal
      end do
c                                 set dependent variables
      call incdp0
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
      call lpopt0 (idead)

      if (idead.eq.0) then
c                                 compute derivative properties
         call getloc (itri,jtri,ijpt,wt,nodata)

         bad = .false.

      else 

         bad = .true.

      end if

      end 