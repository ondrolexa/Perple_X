c      include 'nlib.f'
c      include 'olib.f'
c      include 'clib.f'
c      include 'resub.f'
c      include 'rlib.f'
c      include 'tlib.f'
c      include 'flib.f'

      program meemm       
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,ier,idead

      logical nodata, bulk

      character amount*6, yes*1

      integer itri(4),jtri(4),ijpt

      double precision wt(3), num

      integer iwt
      common/ cst209 /iwt

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k5),ctotal,jdv(k19),npt,fulrnk

      double precision atwt
      common/ cst45 /atwt(k0)

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

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)

      integer iam
      common/ cst4 /iam

      double precision tlv, dt, rho, tst, tlv1, tlv2, rho1, rho2, cp1, 
     *                 cp2, x1, dp, rhoc, x2, pv1, pv2, pvv1, pvv2,
     *                 spec1(5,2),n(2),x(2),molwt(2),specwt(5),nat,
     *                 prps(7,2),xb(2),no,nsi,tot,p,t,lnk1,lnk2,lnk3,
     *                 lnk4,lnk5,ravg,savg

      data specwt/60.084, 44.085, 15.999, 31.998, 28.086/

      integer imax,imin,ir1,ir2,tic,j,lprops(7),lun
      logical go, quit
c                                 N, H, S, Cp, Cp/Cv, rho, vphi
      data lprops/17,2,15,12,28,10,7/
      integer idspec
      double precision spec
      common/ tspec /spec(nsp,k5),idspec

      common/ rcrt /rhoc

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      double precision vp,vvp
      integer rooti
      common/ srkdiv /vp(3),vvp(3),rooti(3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 2
c                                 version info
      call vrsion

      lun = 112

      open (lun,file='lv.tab')
      open (111,file='scratch')

      write (lun,'(a)') '|6.6.6'
      write (lun,*) 'output from LVmeemum'
      write (lun,*) 1
      write (lun,'(a)') 'log(p)'
      write (lun,*) 0.
      write (lun,*) 1.
      write (lun,*) 1
      write (lun,*) 47
      write (lun,'(9a)') 'log(p) p(bar) T(K) dtL(K) dtG(K) xL xG ',
     *    'nL(mol) nG(mol) NL(g-mol) NG(g-mol) NA(g-at) DH(J/kg) ',
     *    'SL SG CpL CpG Cp/CvL Cp/CvG rhoL rhoG v_phiL v_phiG ',
     *    'yL_SiO2 yL_SiO yL_O yL_O2 yL_Si ', 
     *    'yG_SiO2 yG_SiO yG_O yG_O2 yG_Si ',
     *    'pvL pvvL irL pvG pvvG irG ln(K4) ln(C4) ln(K5) ln(C5) ravg ',
     *    'savg ssG rsG'

c      * v(1),10d0**v(1),tlv,tlv-tlv1,tlv2-tlv,
c      * x, n, molwt, prps(1,1), 
c      * prps(2,2) - prps(2,1), 
c                                 S, Cp, Cp/Cv, rho, vphi
c      * ((prps(j,i),i=1,2),j=3,7),
c      * ((spec1(j,i),i=1,2),j=1,5)
c                                 initialization, read files etc. 
      rxn = .false.
      call iniprp

      write (*,1000) 
c     read (*,'(a)') yes
      yes = 'no'

      if (yes.eq.'y'.or.yes.eq.'Y') then 
c                                 bulk is true, user enters composition and p-t conditions
         bulk = .true.

      else 
c                                 else user enters only p-t and composition read from input file.
         bulk = .false.

      end if 

c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      amount = 'molar '

      if (iwt.eq.1) amount = 'weight'
c                                 computational loop
      do 
c                                 read potential variable values    
c                                 v(1) is P(bar), v(2) is T(K) the pointer jv used 
c                                 for general problems but can be eliminated for calculations 
c                                 simply as a f(P,T)       
         write (*,1070) (vname(jv(i)), i = 1, ipot)
         read (*,*,iostat=ier) (v(jv(i)), i = 1, ipot)
         if (ier.ne.0) cycle
         if (v(jv(1)).eq.0d0) exit 
          
         if (bulk) then 
c                                 load the composition into b, the component names are  
c                                 in cname, if iwt = 1 the composition is in mass fractions
c                                 otherwise in molar units. 
            do 
               write (*,1060) amount
               write (*,'(12(a,1x))') (cname(i),i=1,jbulk-1)
               read (*,*,iostat=ier) (cblk(i),i=1,jbulk-1)
               cblk(jbulk) = 1
               write (*,'(g12.6)') cblk(2)/(cblk(2)+cblk(1))
               if (ier.eq.0) exit
            end do  
         
            if (iwt.eq.1) then 
c                                 convert mass to molar 
               do i = 1, jbulk
                  cblk(i) = cblk(i)/atwt(i)
               end do 

            end if
c                                 normalize the composition vector, this 
c                                 is necessary for reasons of stupidity (lpopt0). 
            ctotal = 0d0

            do i = 1, icp
               ctotal = ctotal + cblk(i)
            end do 

            do i = 1, icp
               b(i) = cblk(i)/ctotal
            end do

         end if 
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
         tlv1 = v(2)
         quit = .false.
         rhoc = 1035d0
         dp = nopt(30)

         do

         dt = 10d0
         go = .true.
         v(2) = tlv1 - 5d0
         b(1) = 2d0/3d0
         b(2) = 1d0 - b(1)
         x(2) = 0d0

         do 

            call lpopt0 (idead)

            if (idead.eq.0) then 
    
               call getloc (itri,jtri,ijpt,wt,nodata)

               if (np.gt.1.or..not.go.or.
     *             np.eq.1.and.props(10,1).gt.rhoc) then 

                  tlv = v(2)

                  if (np.gt.1) go = .true.

                  if (go) then

                     if (dt.lt.0d0) dt = -dt/2d0

                     rho = 0d0
                     imax = 0

                     do i = 1, np
                        if (props(10,i).gt.rho.and.
     *                      pcomp(2,i).gt.1d0/3d0) then
                           rho = props(10,i)
                           imax = i
                        end if 
                     end do 

                     do i = 1, icomp
                        cblk(i) = pcomp(i,imax)
                        b(i) = cblk(i)
                     end do 

                      if (b(1).gt.0.99) then
                         ir1 = 0
                      end if 

                     if (np.gt.1) then 
                        if (imax.eq.1) then 
                            x(2) = pcomp(2,2) 
                        else
                            x(2) = pcomp(2,1)
                        end if 
                     end if      

                  end if                         
                   
                  if (dabs(dt).lt.nopt(29)) then

                     if (np.gt.1) then 
                        if (imax.eq.1) then 
                            x(2) = pcomp(2,2) 
                        else
                            x(2) = pcomp(2,1)
                        end if 
                     end if 

                           call getloc (itri,jtri,ijpt,wt,nodata)

                     if (rho.lt.rhoc.or.rho.gt.2d3.and.v(1).gt.2.5d0) 
     *                  then
c                                 either the O-rich phase or the faux-gas, back off 
c                                 to find the plausible root
                        dt = -dabs(dt)
                        tic = 0 

                        do

                           v(2) = v(2) + dt 
                           call lpopt0 (idead)
                           if (idead.ne.0) cycle
c                           if (idead.ne.0.or.np.eq.1) cycle
                           call getloc (itri,jtri,ijpt,wt,nodata)

            if (go) then 
c            call getloc (itri,jtri,ijpt,wt,nodata)
c            call calpr0 (6)
            end if 
 
                           imax = 0 
c comp 1 is O.
                           do i = 1, np
                              if (pcomp(2,i).lt.1d0/3d0) cycle
                              imax = i
                           end do 


                     if (imax.eq.0) then
                        imax = iabs(imax) 
                     end if 

                           if (tic.gt.10.and.dabs(dt).lt.1) then
                              tic = 0
                              dt = 10d0*dt
                              if (dabs(dt).gt.1d0) dt = -1d0
                           end if 

                           tic = tic + 1


                           if (imax.eq.0.or.isnan(props(10,imax))) then
                              write (*,*) 'terrible'
                              stop
                           end if 

                           if (props(10,imax).gt.rhoc.and.
     *                         props(10,imax).lt.2d3.and.v(1).ge.2.5d0
     *                         .or.
     *                      props(10,imax).gt.rhoc.and.
     *                      props(10,imax).lt.2.5d3.and.v(1).lt.2.5d0) 
     * then

            if (go) then 
c            call getloc (itri,jtri,ijpt,wt,nodata)
c            call calpr0 (6)
            end if 

                              exit 
                           end if

                        end do  

                     end if 
c                                 save the liq props:
                     tlv1 = v(2)
                     rho1 = props(10,imax)
                     cp1 = props(12,imax)
                     pv1  = vp(imax)
                     pvv1 = vvp(imax)
                     ir1  = rooti(imax)

                     do j = 1, 5
                        spec1(j,1) = spec(j,imax)
                     end do 

                     do j = 1, 7
                        prps(j,1) = props(lprops(j),imax)
                     end do 

c                                 now increase t to get into the gas composition
                      do i = 1, icomp
                            cblk(i) = pcomp(i,imax) 
                            b(i) = cblk(i)
                      end do  

                      if (b(1).gt.0.99) then
                         ir1 = 0
                      end if 

                      x(1) = pcomp(2,imax) 
                      if (np.gt.1) then 
                         if (imax.eq.1) then 
                            x(2) = pcomp(2,2) 
                         else
                            x(2) = pcomp(2,1)
                         end if 
                      end if 

                      dt = dabs(dt)
                      v(2) = tlv
                      tic = 0 

                      do 

                         v(2) = v(2) + dt 
                         call lpopt0 (idead)
                         if (idead.ne.0.or.np.gt.1) cycle

                         call getloc (itri,jtri,ijpt,wt,nodata)

                         tlv2 = v(2)
                         rho2 = props(10,1)

                         if (isnan(rho2)) then 
                            write (*,*) 'oink'
                            cycle
                         else if (rho2.lt.1d0) then
c                           cycle
                         end if 

                         if (tic.gt.10.and.dt.lt.1) then
                            tic = 0
                            dt = 10d0*dt
                            if (dt.gt.1) dt = 1d0
                         end if 

                         tic = tic + 1

                         if (rho2.gt.0.99*rho1) cycle 
                         cp2 = props(12,1)
                         pv2  = vp(1)
                         pvv2 = vvp(1)
                         ir2  = rooti(1)
                         do j = 1, 5
                           spec1(j,2) = spec(j,1)
                         end do 
                         do j = 1, 7
                           prps(j,2) = props(lprops(j),1)
                         end do 


                         exit
 
                      end do




                      exit

                  end if 

                  if (v(2).gt.62d3) then 

                     quit = .true. 
                     exit 

                  end if 
            
               else if (np.eq.1.and.go.and.dt.gt.0d0) then

                  dt = -dt/2d0

               end if 

             end if 

             v(2) = v(2) + dt

         end do 

         if (quit.and.v(1).gt.3.874d0) exit 
c                                 compute true molar weight
         do i = 1, 2
c                                 atomic weight 
            nat = x(1)*atwt(1) + (1d0-x(1))*atwt(2)
            molwt(i) = 0d0
            
            do j = 1, 5
               molwt(i) = molwt(i) + spec1(j,i) * specwt(j)
            end do

            n(i) = molwt(i)/nat

         end do 





          if (v(2).lt.6200d0)  
     *    write (111,'(120(g14.8,1x))') v(1),tlv,tlv-tlv1,tlv2-tlv,
     *                                x,
     *                                rho1,rho2,cp1,cp2,
     *                                pv1,pvv1,ir1,pv2,pvv2,ir2

          write (*,'(12(g12.6,1x))') v(1),tlv,tlv-tlv1,tlv2-tlv,x,
     *                                rho1,rho2,cp1,cp2
          write (*,'(12(g12.6,1x))') pv1,pvv1,ir1,pv2,pvv2,ir2

          do j = 1, 2 

             write (*,'(12(g12.6,1x))') (prps(i,j),i=1,7)
c                         back calculate compositions
             nsi = spec1(1,j) + spec1(2,j) + spec1(5,j)
             no  = 2d0*spec1(1,j) + spec1(2,j) + 2d0*spec1(4,j)
     *             + spec1(3,j)

             xb(j) = nsi/(nsi+no)
             tot = 0d0
             do i = 1,5
                tot = tot + spec1(i,j)
             end do 
           
             do i = 2,4
                prps(i,j) = prps(i,j)/prps(1,j)*1d3
             end do 

             write (*,'(12(g12.6,1x))') xb(j),x(j),tot
             write (*,'(12(g12.6,1x))') n(j),molwt(j),(spec1(i,j),i=1,5)

          end do 

      ravg = (rho1 + rho2)/2d0
      savg = (prps(3,1) + prps(3,2))/2d0

      t = v(2)
      p = 1d1**(v(1))
      lnk1 = (-0.9214495D6/t + 0.6234471D5)/t - 0.1631235D2
      lnk2 = (-1.133204d+06/t - 5.491882d+04)/t + 1.710990d+01
      lnk3 = (1.906315d6/t - 1.005993d5)/t + 1.664069d1
      lnk4 = lnk1 + lnk2 + lnk3
c      write (*,'(a,12(g12.6,1x))') 'k4 ',lnk4,lnk4+dlog(p)
      lnk5 = lnk1/2d0 + lnk3
c      write (*,'(a,12(g12.6,1x))') 'k5 ',lnk5,lnk5*dlog(p)/2d0    
c      write (113,'(12(g12.6,1x))') 
c     *  lnk4,lnk4+dlog(p),lnk5,lnk5*dlog(p)/2d0  
       write (*,'(a,12(g12.6,1x))') 'rho-s crit ',ravg,savg
c                                 get entropy/rho of stoichiomentric composition
         b(1) = 2d0/3d0
         b(2) = 1d0 - b(1)



            call lpopt0 (idead)

            if (idead.eq.0) call getloc (itri,jtri,ijpt,wt,nodata)

            




          write (lun,'(120(g14.8,1x))') 
     *  v(1),10d0**v(1),tlv,tlv-tlv1,tlv2-tlv,
     *  x, n, molwt, prps(1,1), 
     *  prps(2,2) - prps(2,1), 
     *  ((prps(j,i),i=1,2),j=3,7),
     *  ((spec1(j,i),j=1,5),i=1,2),
     *   pv1,pvv1,ir1,pv2,pvv2,ir2,
     *   lnk4,lnk4+dlog(p),lnk5,lnk5*dlog(p)/2d0,ravg,savg,
     *   props(15,1)/props(17,1)*1e3,props(10,1)


         if (v(1).ge.2.795d0.and.dp.ge.0.09) then
            dp = 0.01
         else if (v(1).ge.3.77d0.and.dp.ge.0.009) then
            dp = 0.00001
c         else if (v(1).ge.3.86d0.and.dp.ge.0.0009) then 
c            dp = 0.0001
         end if                       
 
         v(1) = v(1) + dp 

         end do 

         if (idead.gt.0) then
 
            write (*,*) 'minimization failed'

         else 
c                                 compute derivative properties
            call getloc (itri,jtri,ijpt,wt,nodata)
c                                 print summary to LUN 6
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)

         end if 

         if (goodc(1)+badc(1).gt.0d0) then
            num = badc(1)/(badc(1)+goodc(1))*1d2
            if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')
            goodc(1) = 0d0
            badc(1) = 0d0 
         end if 

      end do

1000  format (/,'Interactively enter bulk compositions (y/n)?',/,
     *          'If you answer no, MEEMUM uses the bulk composition',
     *         ' specified in the input file.',/)
1060  format (/,'Enter ',a,' amounts of the components:')
1070  format (/,'Enter (zeroes to quit) ',7(a,1x))

      end 
    
      subroutine iniprp
c----------------------------------------------------------------------
c iniprp - read data files and initialization for meemum
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, output, err 

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod
c----------------------------------------------------------------------- 
      first = .true.
      output = .false.
      err = .false.
c                                 elastic modulii flag
      kmod = 0 
c                                 -------------------------------------------
c                                 open statements for units n1-n5 and n9
c                                 are in subroutine input1
      call input1 (first,output,err)
c                                 for meemum turn auto_refine OFF
      iopt(6) = 0 
c                                 read thermodynamic data on unit n2:
      call input2 (first)
c                                 allow reading of auto-refine data 
      call setau1 (output)
c                                 read data for solution phases on n9:
      call input9 (first,output)
c                                 call initlp to initialize arrays 
c                                 for optimization.
      call initlp     

      end
