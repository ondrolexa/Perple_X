      program cohsrk
c----------------------------------------------------------------------
c                       ************************
c                       *                      *
c                       *    cohsrk.6.1993     *
c                       *                      *
c                       ************************
c----------------------------------------------------------------------

c COHSRK is a fortran program to call various subroutines to calculate  
c C-O-H-S fluid speciation as a function of X(O), X(CO2), X(S), X(C) 
c a(C), f(S2), f(O2) or O2-S2-buffer assemblage. 

c The variables X(O), X(S), and X(C) are defined as:

c  X(O) = n(O)/{n(O)+n(H)}
c  X(S) = n(S)/{n(S)+n(C)}
c  X(C) = n(C)/{n(S)+n(C)+n(O)+n(H)}

c where n(C), n(S), n(O) and n(H)  are the total number of moles of 
c carbon, sulfur, oxygen and hydrogen (as opposed to the amounts
c of the species) in the fluid.

c The user may choose from the following routines, identified by
c number:

c    0 - Modified Redlich-Kwong (MRK)
c    1 - Kerrick & Jacobs 1981 hard sphere MRK (HSMRK)
c    2 - Hybrid MRK/HSMRK
c    3 - Saxena & Fei 1987 pseudo-virial expansion
c    4 - Bottinga & Richet 1981 (CO2 RK)
c    5 - Holland & Powell 1990 (CORK)
c    6 - Hybrid Haar et al 1979/HSMRK (TRKMRK)
c    7 - Graphite buffered COH MRK fluid
c    8 - Graphite buffered COH hybrid-RK fluid
c    9 - Maximum X(H2O) GCOH fluid Cesare & Connolly 1993
c   10 - X(O) GCOH-fluid hybrid-MRK Connolly & Cesare 1993
c   11 - X(O) GCOH-fluid MRK Connolly and Cesare 1993
c   12 - X(O) GCOHS Connolly & Cesare 1993
c   13 - X(H) H2-H2O-hybrid
c   14 - hogbrd, don't use this if you dont know what it is.
c   15 - X(H) low T H2-H2O-hybrid
c   16 - X(O) H-O HSMRK/MRK hybrid
c   17 - X(O) H-O-S HSMRK/MRK hybrid
c   18 - Delany/HSMRK/MRK hybrid, for P > 10 kb
c   19 - X(O)-X(S) GCOHS Connolly & Cesare 1993
c   20 - X(O)-X(C) GCOHS Connolly & Cesare 1993
c   24 - f(O2)-N/C graphite saturated COHN MRK fluid

c Routines 0-3, 5-6, and 18 are for conventional P-T-X(CO2) 
c calculations, where X(CO2) is the mole fraction of CO2 in 
c a binary H2O-CO2 mixture. The Saxena & Fei routine was
c programmed by someone else, and I suspect it has been 
c entered incorrectly. 

c Routine 4 is for pure CO2 fluids.

c Routines 7 and 8 are for C-O-H fluids as a function of 
c P-T and f(O2) or f(O2)-buffer at specified 
c graphite activity (usually 1). I recommend choice 8.

c Routine 9 returns the fugacities of O2, H2O, and CO2
c for a H:O = 2 graphite saturated C-O-H fluid, i.e.,
c the fugacities computed with routine 10 at X(O) = 1/3
c and a(graphite) = 1.

c Routines 10 and 11 are for C-O-H fluids as a function 
c of X(O) at specified graphite activity. I recommend 10.

c Routines 12 and 17 are for C-O-H-S and H-O-S fluids as
c a function of P-T-X(O) conditions at specified sulfur
c fugacity or buffer, for C-O-H-S fluids graphite activity
c must also be specified. 

c Routines 13 and 15 are for binary H2/H2O mixtures 
c (effectively H-O fluids at X(O) < 1/3), routine 15 uses
c parameters optimized for the solvus region.

c Routine 16 is for H-O fluids as a function of P-T-X(O).

c Routines 2, 8, 12, 13, 14, 16 and 17 are a hybrids of the 
c HSMRK EoS (Kerrick & Jacobs) and the MRK (Redlich & Kwong/
c DeSantis et al/Holloway) EoS as described by Connolly & 
c Cesare, J Met Geol, 11:379-388, and for most purposes I would
c recommend these equations.

c Routine 19 calculates C-O-H-S graphite-UNDERSATURATED fluid 
c properties as a function of X(O) and the atomic S/C ratio 
c (expressed by the variable X(S)) at specified sulfur fugacity. 
c This routine often may not converge. The EoS is identical to
c that used in Routine 12.

c Routine 20 calculates C-O-H-S graphite-UNDERSATURATED fluid 
c properties as a function of X(O) and the atomic fraction of
c carbon in the fluid (X(C)) at specified sulfur fugacity.
c This routine can be made to calculate simple COH fluid 
c speciation as a function of bulk composition by setting 
c sulfur fugacity to such a low value that the concentration of 
c sulfur becomes negligible.
c This routine often may not converge. The EoS is identical to
c that used in Routine 12.

c Routine 24 calculates C-O-H-N graphite saturated fluid properties
c as a function of f(O2) and the fluid N/C ratio using a modified
c Redlich-Kwong EoS. The routine will fail in the limits x(H2O)->0
c and N/C->0.

c COHSRK is primarily intended to provide an example of how the
c various speciation routines in PERPLEX can be called. If you
c can follow the logic of the program then you shall probably be
c able to customize the program to use input, and generate output,
c more in line with your specific needs. I do have another program
c that generates a plot file with various species, concentrations,
c and fugacities plotted against X(O), and I shall be happy to 
c provide it to interested users.

c----------------------------------------------------------------------

c Limitations/warnings:

c The C-O-H-S and H-O-S routines do not consider the species C2H6,
c O2, and SO3 and pure S species, these will not be significant at
c the conditions realized in crustal metamorphic environments. 
c It is relatively easy to incorporate new species in the routines
c so if you suspect the species are important, please contact me
c and I will add the species.

c The HSMRK equation of state for water seems unreliable at P > 
c 15 Kbar, and hence all the hybrid EoS speciation routines
c as well. This has not bothered me because I am mostly concerned  
c with crustal metamorphic conditions; however, should you wish
c to make calculations at high pressure you need only replace the
c calls which get pure fluid fugacities with calls to the high
c pressure routine of your choice (e.g., DeSantis-Holloway MRK,
c Delany & Helgeson, Bottinga & Richet, etc. etc.).

c Functions for equilibrium constants are fit in the range 400-1400K.

c No warnings are issued if P-T-X(O) conditions are within a solvus.

c All X(O) routines fail at X(O) = 0 or 1, to avoid this problem
c the routines reset extreme values of X(O) to 1d-10 or 0.9999999999
c as appropriate.

c All X(O) routines may fail at X(O) = 1/3 if the fluid becomes an
c essentially mono-species H2O fluid, this is most likely to occur
c in simple H-O fluids.

c The C-O-H-S and H-O-S routines permit calculations for a fluid
c in equilibrium with pyrrhotite with a fixed atomic Fe/S ratio,
c this is not the same as N(FeS), the mole fraction of troilite
c relative to S2, as often reported for pyrrhotite analyses.

c Routines 19 and 20 may not converge, or may converge to 
c invalid roots, in the vicinity of graphite saturated and
c supersaturated conditions; and at very carbon under-saturated 
c conditions numerical slop may lead to errors on the order of
c 0.4 log units in calculated properties such as f(O2).

c In speciation routines, if the concentration of a species becomes 
c negligible (as determined by the numerical precision of the
c computer being used), its properties may be undefined or
c may have meaningless fluctuations.
c----------------------------------------------------------------------

c NOTES ON THIS SOURCE AND COMPILING:

c The MAIN program COHSRK and a BLOCK DATA are included in this file,
c all the subroutines called by COHSRK are in the file flib.f.
c The files may be compiled separately and the objects linked together,
c or you may concatenate the two files. 

c This program was put together from routines used in PERPLEX (see
c below) and for this reason it is not as compact or as flexible
c as might be desired. The size of the flib.f may cause problems
c for DOS compilers, these can be overcome by splitting the source
c into two or three blocks that can then be combined during linking.
c Alternatively, routines that are not considered may be useful
c can be eliminated if the calls to these routines from subprogram
c cfluid are also eliminated.

c Specifically I would recommend eliminating:

c      subroutine brmrk
c      subroutine simps 
c      subroutine qromb 
c      subroutine polint 
c      subroutine trapzd 
c      function brvol 
c      function vdpdv 
c      subroutine hosrk5 
c      subroutine cohfit 
c      subroutine haar
c      subroutine psat2
c      subroutine aideal 
c      subroutine trkmrk
c      subroutine saxfei
c      subroutine hprk
c      subroutine cohgra 
c      subroutine hh2ork 
c      subroutine lohork
c      subroutine lomrk 

c this requires elimination of the calls to: brmrk, cohfit,
c trkmrk, saxfei, hprk, cohgra, hh2ork, and lohork from cfluid.

c subroutines warn and error can also be eliminated if the 
c calls to these routines, from subroutines rfluid, brmrk,
c cohgra, cohsgr and cohhyb, are replaced by statements that write 
c an appropriate warning or error message. 

c----------------------------------------------------------------------

c These routines are incorporated in the PERPLEX programs for calculating
c phase equilibria and diagrams, if you would like a copy of these 
c programs, or if you have any problems with this program, please
c contact me at:

c                     James Connolly
c                     IMP-ETHZ
c                     CH-8092 Zuerich

c by e-mail at:

c                     jamie@erdw.ethz.ch

c or by telephone/fax at:

c                     0041-1-632-7804/0041-1-632-1088

c-----------------------------------------------------------------------
c I/O: most input and output is done through the common blocks below,
c the significance of the variables as named in the main program is:

c ifug   - number indexing the requested EoS.
c p      - pressure, bars.
c t      - temperature, kelvins.
c xo     - X(O) for multispecies routines, and X(CO2) or X(H2) for
c          binary routines.
c vol    - molar volume for all multispecies routines, and binary
c          routines 0, 1, 13, 15.
c fh2o   - natural log (f(H2O)) for all routines.
c fco2   - natural log of the species other than H2O (i.e., CO2 or H2)
c          in all binary routines.
c fo2    - natural log (f(O2)).
c fs2    - 1/2 natural log (f(S2)).

c the following variables are only used for multispecies calculations:

c ins(i) - pointers indicating the species are to be calculated.
c isp    - number of species to be calculated.
c nsp    - dimensioning for the maximum number of species. 
c ibuf   - pointer indicating method of calculating f(O2) for
c          routines 7 and 8, or f(S2) for routines 12 and 17.
c dlnfo2 - for routines 7 and 8 the displacement of the f(O2)
c          relative to a buffer (in log units) or the absolute
c          ln(f(O2)) (if ibuf = 3). For routines 12 and 17,
c          if ibuf = 2 then dlnfo2 is the atomic Fe/S of pyrrhotite,
c          if ibuf = 3 then dlnfo2 is 1/2 the natural log (f(S2)). 
c elag   - natural log of graphite activity.
c g(i)   - fugacity coefficient of the ith species.
c x(i)   - mole fraction of the ith species.

c the indices of eleven species presently defined are:

c         1 = H2O
c         2 = CO2
c         3 = CO
c         4 = CH4 
c         5 = H2
c         6 = H2S
c         7 = O2
c         8 = SO2
c         9 = COS
c        10 = N2
c        11 = NH3

c O2 should be replaced by SO3, and ethane should be added.
c-----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      character specie(nsp)*3, vname*8, y*1

      integer ier, igo, ins(nsp), i, isp, j, k, kmax

      double precision nc, nh, no, ns, nn, tentoe, fo2, fs2, xfh, 
     *                 xfc, ag, tot, totx

      double precision fh2o,fco2
      common / cst11 /fh2o,fco2

      double precision p,t,xo,u
      common/ cst5 /p,t,xo,u(6)

      double precision xs,g,v
      common/ cstcoh /xs(nsp),g(nsp),v(nsp)

      double precision vol
      common/ cst26 /vol

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10 /iff(2),idss(h5),ifug,ifyn,isyn

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      integer iam
      common/ cst4 /iam

      data tentoe, fo2, fs2, specie /2.302585093d0, 0.d0, 0d0,
     *      'H2O','CO2','CO ','CH4','H2 ','H2S','O2 ',
     *      'SO2','COS','N2 ','NH3'/
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 11
c                                 version info
      call vrsion

50    igo = 0
      fs2 = -9999d0*tentoe/2d0
      elag = 0d0
c                                  initialize species fractions/volumes
      do i = 1, nsp
         xs(i) = 0d0
      end do 

      vname = ' X(CO2) '
c                                  get the users choice of EoS:   
      call rfluid (1, ifug)
c                                  for multispecies fluids set
c                                  up species indices:
      if (ifug.gt.6.and.ifug.lt.13.or.
     *    ifug.eq.19.or.ifug.eq.20.or.ifug.eq.24) then
         vname = '  X(O)  '
         if (ifug.eq.7.or.ifug.eq.8.or.ifug.eq.24) vname = 'log(fo2)'
         isp = 5

         do i = 1, 6
            ins(i) = i
         end do
 
         if (ifug.eq.19.or.ifug.eq.20) then
            isp = 8
            ins(7) = 8
            ins(8) = 9
         else if (ifug.ge.12.and.ifug.lt.19) then
            isp = 9
            ins(7) = 7
            ins(8) = 8
            ins(9) = 9
         else if (ifug.eq.24) then
            isp = 7
            ins(6) = 10
            ins(7) = 11
         end if

      else if (ifug.eq.16) then
         vname = '  X(O)  '
         isp = 3
         ins(1) = 1
         ins(2) = 5
         ins(3) = 7
      else if (ifug.eq.17) then
         vname = '  X(O)  '
         isp = 4
         ins(1) = 1
         ins(2) = 5
         ins(3) = 6
         ins(4) = 8
      else if (ifug.eq.13.or.ifug.eq.15) then 
         vname = ' X(H2)  '
      end if 

      write (*,1120)

10    if ((igo.eq.0.and.(ifug.eq.7.or.ifug.eq.8)).or.
     *    (ibuf.ne.3.and.(ifug.eq.7.or.ifug.eq.8)).or.
     *    ifug.eq.4.or.ifug.eq.9) then 
c                                  get P-T conditions:
         write (*,1010) 
         xo = 1d0
         igo = 1
         read (*,*,iostat=ier) p, t
         call rerror (ier, *10)
      
      else if (ifug.eq.19) then

         write (*,1340)
         read (*,*,iostat=ier) p, t, xo, elag
         call rerror (ier, *10)

      else if (ifug.eq.20) then

         write (*,1360)
         read (*,*,iostat=ier) p, t, xo, elag
         call rerror (ier, *10)

      else if (ifug.eq.24) then 
         
         if (ibuf.ne.3) then 
            write (*,1390)
            read (*,*,iostat=ier) p, t, gz
            call rerror (ier, *10)
         else 
            write (*,1410)
            read (*,*,iostat=ier) p, t, dlnfo2, gz
            call rerror (ier, *10)
            dlnfo2 = tentoe*dlnfo2
         end if 

         if (igo.eq.0) then 
            igo = 1
            write (*,1430) 
            read (*,*,iostat=ier) elag
            call rerror (ier, *10)
         end if 

      else 
c                                  or get P-T-X/f conditions:
         write (*,1000) vname
         read (*,*,iostat=ier) p, t, xo
         call rerror (ier, *10)
         if (ifug.eq.7.or.ifug.eq.8.or.ifug.eq.24) dlnfo2 = tentoe * xo

      end if 
c                                  quit if p = 0
      if (p.eq.0d0) goto 99
c                                  if sulfur dependent, get
c                                  fs2 if user hasn't opted for
c                                  a buffer:
      if ( (ifug.eq.12.and.ibuf.eq.3.and.igo.ne.0).or.
     *     (ifug.eq.17.and.ibuf.eq.3.and.igo.ne.0).or.
     *     (ifug.eq.19.and.ibuf.eq.3.and.igo.ne.0).or.
     *     (ifug.eq.20.and.ibuf.eq.3.and.igo.ne.0) ) then 

30       write (*,1110)
         read (*,*,iostat=ier) dlnfo2
         call rerror (ier,*30)
         dlnfo2 = tentoe * dlnfo2
      end if 
c                                  call fluid routine:
      call cfluid (fo2, fs2)
      write (*,*) ' '
      call rfluid (3, ifug)
      write (*,1280) p,t
c                                  output results:
      igo = 1
      xfh = fh2o
      xfc = fco2
      fh2o = dexp(fh2o) 
      fco2 = dexp(fco2) 
      fo2 = fo2 / tentoe

      if (ifug.lt.4.or.ifug.eq.5.or.ifug.eq.6.or.ifug.eq.18.or.
     *    ifug.eq.14.or.ifug.gt.20.and.ifug.lt.24) then

         write (*,1130) fh2o, fco2
c                                  increment pressure for
c                                  finite difference estimate of
c                                  volume:
         p = p + 1d0
         call cfluid (fo2, fs2)
         write (*,1300) 83.14d0*t*((1d0-xo)*(fh2o-xfh)+xo*(fco2-xfc))
      else if (ifug.eq.4) then
         write (*,1140) fco2
         write (*,1300) vol
      else if (ifug.eq.9) then
         write (*,1150) fh2o, fco2, fo2
      else if (ifug.eq.13.or.ifug.eq.15) then
         write (*,1160) fh2o,fco2,fo2
         write (*,1300) vol
      else 
         if (ifug.eq.16.or.ifug.eq.17) then
            ag = 0d0
         else if (ifug.eq.19.or.ifug.eq.20) then
            ag = dexp (gz)
         else
            ag = dexp (elag)
         end if 
c                                  routine cfluid returns ln(fs2)/2
         write (*,1170) fo2, 2d0*fs2/tentoe, ag
c                                  output speciation:
         write (*,1230)
         do 20 j = 1, isp, 4
            kmax = j + 3
            if (kmax.gt.isp) kmax = isp
            write (*,1180) (specie(ins(k)), k = j, kmax)
            write (*,1190) (xs(ins(k)), k = j, kmax)
20          write (*,1200) (g(ins(k))*p*xs(ins(k)), k = j, kmax)
c                                  total species fractions:
         totx = 0d0
         do 21 k = 1, isp
21          totx = totx + xs(ins(k))

         write (*,1370) totx
         if (totx.gt.1.001d0.or.totx.lt.0.999d0) write (*,1380)

c                                  output bulk properties and V:
         write (*,1240)

         ns = xs(6) + xs(8) + xs(9) 
         no = xs(1) + xs(2)*2d0 + xs(3) + xs(7)*2d0 
     *               + xs(8)*2d0 + xs(9) 
         nc = xs(2) + xs(3) + xs(4) + xs(9)
         nh = (xs(1) + xs(5) + xs(6))*2d0 + xs(4)*4d0 + xs(11)*3d0 
         nn = 2d0*xs(10) + xs(11)

         tot = ns + no + nh + nc + nn     

         write (*,1250) nc/tot, nh/tot, no/tot, ns/tot, nn/tot
         write (*,1270) no/(no+nh)
         if (ns.ne.0d0) write (*,1310) ns/(ns+nc)
         write (*,1350) nc/(nc+no+nh+ns)
         if (ifug.eq.24) write (*,1400) nn/nc
         if (nh.ne.0d0.and.ns.ne.0d0) write (*,1290) ns/nh
         if (ifug.eq.19.or.ifug.eq.20.and.ag.gt.1d0) write (*,1330)

         write (*,1420) vol

      end if 

      goto 10 

1000  format (/,'Enter p(bar), T(K), ',a,': ')
1010  format (/,'Enter p(bar), T(K): ')
1110  format (/,'Enter the log10[f(S2)]:',/)
1120  format (/,'Enter a zero for pressure to quit.',/)
1130  format (/,10x,'f(H2O) = ',g12.5,/,10x,'f(CO2) = ',g12.5,/)
1140  format (/,10x,'f(CO2) = ',g12.5,/)
1150  format (/,10x,'f(H2O)     = ',g12.5,/
     *         ,10x,'f(CO2)     = ',g12.5,/
     *         ,10x,'log[f(O2)] = ',g12.5,/)
1160  format (/,10x,'f(H2O)     = ',g12.5,
     *        /,10x,'f(H2 )     = ',g12.5,
     *        /,10x,'log[f(O2)] = ',g12.5,/)
1170  format (/,10x,'log[f(O2)] = ',g12.5,
     *        /,10x,'log[f(S2)] = ',g12.5,
     *        /,10x,'a(gph/dia) = ',g12.5,/)
1180  format (10x,4(5x,a3,5x))
1190  format (5x,'x',4x,4(g12.5,1x))
1200  format (5x,'f',4x,4(g12.5,1x),/)
1210  format (/,'Change EoS, buffer, or graphite activity (y/n)?',/)
1220  format (a1)
1230  format (22x,'Speciation/Fugacities',/)
1240  format (/,22x,'Atomic Proportions',//,
     *        10x,'C',12x,'H',12x,'O',12x,'S',12x,'N')
1250  format (4x,5(g12.5,1x),/)
1270  format (5x,'Back-calculated X(O) = ',g16.9)
1280  format (/,10x,'p(bar)     = ',g12.5,/,10x,'T(K)       = ',g12.5)
1290  format (5x,'S/H = ',g12.5,/)
1300  format (10x,'V(cm3/mol) = ',g12.5,/)
1310  format (5x,'Back-calculated X(S) = ',g16.9)
1330  format (/,5x,'COMPOSITION IS SUPERSATURATED WITH ',
     *           'RESPECT TO CARBON!!',/)
1340  format (/,'Enter p(bar), T(K), X(O), X(S): ')
1350  format (5x,'Back-calculated X(C) = ',g16.9)
1360  format (/,'Enter p(bar), T(K), X(C):')
1370  format (/,5x,'Sum of species fractions: ',f14.9,/)
1380  format (/,5x,'INVALID SPECIIATION!!',/)
1390  format (/,'Enter p(bar), T(K), molar N/C:')
1400  format (5x,'Back-calculated N/C  = ',g16.9)
1410  format (/,'Enter p(bar), T(K), log10[f(O2)], molar N/C:')
1420  format (/,5x,'Molar Volume (cm3/mol) = ',g12.5,/)
1430  format (/,'Enter log10[a(gph/dia)]:')

99    write (*,1210) 
      read (*,1220) y
      if (y.eq.'y'.or.y.eq.'Y') goto 50

      end 
