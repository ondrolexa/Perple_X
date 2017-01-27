c routines common to all programs? could be in tlib.f?

      subroutine setau1 (output)
c----------------------------------------------------------------------
c setau1 sets autorefine dependent parameters. called by vertex, werami,
c pssect and meemum.

c output is set to false if autorefine mode is not auto (i.e., iopt(6) = 2) 
c or it is auto and in the second cycle.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical output
 
      character*8 y*1

      character*10 badnam(h9)

      integer ibad2,ibad1,igood,i,j,ierr

      character*100 n10nam,n11nam,n12nam

      character*100 prject,tfname
      common/ cst228 /prject,tfname
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      logical refine
      common/ cxt26 /refine

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(4,2)

      integer iam
      common/ cst4 /iam
c-----------------------------------------------------------------------
      refine = .false.
c                                 only use autorefine if solutions
c                                 are present and it is requested.
      if (isoct.ne.0) then 

         call mertxt (n10nam,prject,'.arf',0)
         open (n10, file = n10nam, iostat = ierr, status = 'old')

         call mertxt (n12nam,prject,'.tof',0)

         if (iam.eq.1.or.iam.eq.2) then
c                                 VERTEX or MEEMUM:
            if (iam.eq.1) then 

               open (n8, file = n12nam, status = 'unknown')
c                                 user friendly text version 
               if (lopt(11)) then 
                  call mertxt (n11nam,prject,'_auto_refine.txt',0)
                  open (n11, file = n11nam, status = 'unknown')
               end if 

            end if 

            ibad1 = 0 

            if (ierr.ne.0.and.iam.eq.1) then 
c                                 no auto_refine data
               write (*,1020) n10nam
               open (n10, file = n10nam, status = 'unknown')

            else if (ierr.eq.0.and.iam.eq.1) then 
                   
               read (n10,*,iostat=ierr) ibad1, ibad2, igood
               if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)

               if (iopt(6).ne.2.or.output) write (*,1030) n10nam

               if (iopt(6).eq.1) then 
c                                 manual mode, allow reinitialization
c                                 or suppression.
                  write (*,1060) 
                  read (*,'(a)') y

                  if (y.eq.'y'.or.y.eq.'Y') then

                     iopt(6) = 0

                  else 

                     refine = .true.  

                  end if

                  output = .true.
 
               else if (output) then  
c                                 second cycle of automated mode
                  refine = .true.

               end if  

               write (n8,*) refine
          
            else if (ierr.eq.0.and.iam.eq.2) then 
c                                 MEEMUM, ask the user if he wants
c                                 to use the data 
               write (*,'(/,a,a,/,a)') 'Auto-refine data exists from a',
     *                  ' previous calculation with VERTEX.',
     *                   'Do you want MEEMUM to use this data (y/n)?'
               read (*,'(a)') y

               if (y.ne.'y'.and.y.ne.'Y') then

                  iopt(6) = 0

               else 

                  refine = .true.  
                  read (n10,*,iostat=ierr) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
                  iopt(6) = 1
                  write (*,1030) n10nam

               end if

            end if 
c                                 set cycle dependent parameters
            if (refine) then 

               i = 2 

            else 

               i = 1

            end if
c                                 solvus tolerance 
            if (lopt(9)) nopt(8) = 1.5d0*rid(3,i)
c                                 number of iterations
            iopt(10) = grid(6,i)
c                                 bound relaxation rate
            if (i.eq.1) then 
c                                 the initial resolution of the exploratory stage
               nopt(10) = rid(3,1) * nopt(29)
           
            else 
c                                 the final resolution of the exploratory stage   
               nopt(10) = rid(4,1) * nopt(29)

            end if 

         else 
c                                 werami/pssect if refine, get the 
c                                 solution models to be rejected
            open (n8, file = n12nam, iostat=ierr, status = 'old')
        
            if (ierr.eq.0) then 
c                                 write a flag to indicate if auto-refine
c                                 has been used, this is necessary so that other
c                                 perplex programs know whether to reject the
c                                 badnam phases:
               read (n8,*,iostat=ierr) refine
c                                 read phases to be rejected if in auto-refine
               if (refine) then 
                  read (n10,*,iostat=ierr) ibad1, ibad2, igood
                  if (ibad1.gt.0) read (n10,'(a)') (badnam(i),i=1,ibad1)
               end if 

            end if 

         end if 

      end if 

      close (n8)
c                                 just to be sure
      if (iopt(6).eq.0) refine = .false.

      if (refine) then 
c                                 reject solution models that were 
c                                 not found to be stable and set parameters 
c                                 that depend on refinement
         ibad2 = 0 

         do 50 i = 1, isoct

            do j = 1, ibad1
               if (fname(i).eq.badnam(j)) then
                  if (iam.eq.1) write (*,1070) fname(i)
                  goto 50
               end if 
            end do 

            ibad2 = ibad2 + 1
            fname(ibad2) = fname(i)

50       continue 

         isoct = ibad2 

         write (*,'(/)')

      end if

      if (iopt(6).eq.2.and..not.refine) then 
         output = .false.
      else
         output = .true.
      end if 

1020  format (/,'Writing data for auto-refinement to file: ',a,/)
1030  format (/,'Reading data for auto-refinement from file: ',a,/)
1060  format ('Suppress or reinitialize auto-refinement (y/n)?')
1070  format ('Eliminating solution model: ',a,' in auto-refinement.')

      end 

      subroutine setau2 (output)
c----------------------------------------------------------------------
c setau2 sets/resets autorefine parameters after the solution models have
c been read. setau1 must be called first.

c output is set to true if autorefine mode is auto (i.e., iopt(6) = 2) 
c but no solutions are present (isoct = 0). 
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical output

      integer i,index
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer grid
      double precision rid 
      common/ cst327 /grid(6,2),rid(4,2)

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      logical oned
      common/ cst82 /oned

      logical refine
      common/ cxt26 /refine
c-----------------------------------------------------------------------
      if (isoct.eq.0) then 
     
         index = 2
         output = .true.

      else if (.not.output) then

         index = 1

      else 

          if (refine) then

             index = 2

          else 

             index = 1

          end if 

      end if 
c                                 set auto-refine dependent parameters
      if (icopt.eq.5) then 
c                                 gridded minimization
         if (oned) then 
            jlow = grid(4,index)
            loopx = 1
         else 
            jlow = grid(2,index)
            loopx = grid(1,index) 
         end if

         jlev = grid(3,index) 
          
      else if (icopt.gt.5) then 
c                                 1d/2d phase fractionation
         jlow = grid(4,index)

      else if (icopt.eq.1) then 
c                                 schreinemakers diagrams

c                                 max variance of curves to be traced
          isudo = grid(5,index)
c                                 default variable tracing increment
          do i = 1, 2
             dv(iv(i)) = (vmax(iv(i)) - vmin(iv(i)))*rid(1,index)
          end do 

      else if (icopt.eq.3) then 
c                                 mixed variable diagrams 

c                                 no variance restriction
          isudo = 99
c                                 default search increment
          dv(iv(1)) = (vmax(iv(1)) - vmin(iv(1)))*rid(1,index)

      end if 

      end 

      subroutine input1 (first,output,err)
c-----------------------------------------------------------------------
c input1 reads data from a file on unit n1, this data controls the
c computational options and is modified frequently.

c iam - indicates calling program 1 - vertex
c                                 2 - meemum
c                                 3 - werami
c                                13 - unsplt, global call
c                                14 - unsplt, local call
c                                 any other values no output
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      logical output, eof, first, err

      character*100 blank*1,string(3)*8,rname*5,name*8,strg*80,n2name,
     *              n9name,y*1,sname*10,prt*3,plt*3

      integer idum,nstrg,i,j,ierr,icmpn,jcont,kct

      double precision dip

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio

      character*100 cfname
      common/ cst227 /cfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      character*162 title
      common/ csta8 /title(4)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      character*8 xname,vname
      common/ csta2 /xname(k5),vname(l2)

      character*5 cname
      common/ csta4 /cname(k5) 

      integer icp2
      common/ cst81 /icp2

      character*5 zname
      common/ cst209a /zname

      character tcname*5,xcmpnt*5
      common/ csta9 /tcname(k0),xcmpnt(k0)

      double precision buf
      common/ cst112 /buf(5)

      integer iwt
      common/ cst209 /iwt

      integer ivfl
      common/ cst102 /ivfl

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ibuf,hu,hv,hw,hx 
      double precision dlnfo2,elag,gz,gy,gx
      common/ cst100 /dlnfo2,elag,gz,gy,gx,ibuf,hu,hv,hw,hx

      double precision ctrans
      integer ictr,itrans
      common/ cst207 /ctrans(k0,k5),ictr(k5),itrans

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer isoct
      common/ cst79 /isoct

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ivarrx,ivarip,isudo,ivar
      common/ cst62 /ivarrx(k2),ivarip(k2),isudo,ivar

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      logical oned
      common/ cst82 /oned

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      character*8 eoscmp
      common/ cst98 /eoscmp(2)

      integer iam
      common/ cst4 /iam

      save blank
      data blank/' '/
c-----------------------------------------------------------------------
c                             output = .false. then in 1st cycle of
c                             autorefine.
      if (.not.output) then 
c                                 read computational option file 
         call fopen1 
      
      else 
c                                 create the file name           
         call mertxt (tfname,prject,'.dat',0)
         open (n1, file = tfname, iostat = ierr, status = 'old')
         if (ierr.ne.0) call error (120,r,n1,tfname)

      end if 
c                                 begin reading input:

c                                 read name of thermodynamic data file
      read (n1,'(a)') n2name
      call enblnk (n2name)
c                                 read print and graphic file names
      read (n1,'(a)') prt

      read (n1,'(a)') plt

      read (n1,'(a)') n9name
      call enblnk (n9name)
c
      do i = 1, 4
         title(i) = ' '
      end do 
c                                 read title for the calculation:
      read (n1,'(a)') title(1)
c                                 read computational option or option file name
c                                 use error condition to determine which:
      read (n1,'(a)') tfname
c                                 get first non-blank string 
      call getstg (tfname)

      read (tfname,'(i2)',iostat=ierr) icopt 

      if (ierr.eq.0) then 
c                                 if no error, old version
         tfname = 'perplex_option.dat'

      else
c                                 new version, read icopt
         read (n1,*,err=998) icopt

      end if 
c                                 if meemum, override whatever computational option
c                                 is set in the input file. 
      if (iam.eq.2) icopt = 5
c                                 if fractionation path from data 
c                                 file, get name:
      fileio = .false.

      if (icopt.eq.10.or.icopt.eq.11) then 

         fileio = .true.

         read (n1,'(a)') cfname
         call enblnk (cfname)

         if (icopt.eq.10) then 
            icopt = 7
         else 
            icopt = 9
         end if 

      end if 
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum 
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum

      read (n1,*,err=998) itrans
      read (n1,*,err=998) icmpn
c                                 read new component definitions:
      do i = 1, itrans
         read (n1,'(a,1x,i2)') tcname(i), ictr(i)
         read (n1,*) (ctrans(j,i), j = 1, icmpn)
      end do

      read (n1,*,err=998) iwt
c                                 dummy variable place holders
      read (n1,*,err=998) idum
      read (n1,*,err=998) idum  
      read (n1,*,err=998) idum
c                                 read code for choice of fluid equation
c                                 of state from terminal. 
      read (n1,*,err=998) ifug
      if (ifug.ge.7.and.ifug.le.12.and.ifug.ne.9.and.ifug.ne.14.or.
     *    ifug.eq.19.or.ifug.eq.16.or.ifug.eq.17.or.ifug.eq.24.or.
     *    ifug.eq.20.or.ifug.eq.25) 
     *                  read (n1,*,err=998) ibuf,hu,dlnfo2,elag
      if (ibuf.eq.5) read (n1,*,err=998) buf     

      if (hu.eq.1) then 
c                                 hardwired fluid EoS endmember names
         eoscmp(1) = 'H2      '
         eoscmp(2) = 'O2      '

      else 

         eoscmp(1) = 'H2O     '
         eoscmp(2) = 'CO2     '

      end if 
c                                 no dependent variable
      iind = 0 
c                                 dummy variable
      read (n1,*,err=998) loopx
c                                 here loopx is just a 1d/2d flag for 
c                                 gridded minimization, for backwards 
c                                 compatibility set the to 2d if > 2 or < 1.
      if (loopx.eq.1) then 
         oned = .true.
      else
         oned = .false.
      end if 

      read (n1,*,err=998) idep
      read (n1,*,err=998) c0,c1,c2,c3,c4

      if (idep.eq.1) then 
         iind = 2
      else if (idep.eq.2) then 
         iind = 1
      end if 
c                                 decode thermodynamic components
c                                 read to the beginning of the component list
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 count (icp) and save names (cname)
      icp = 0
      jbulk = 0 
  
      do 

         read (n1,'(a,a)') rname,strg

         if (rname.eq.'end t') then 
c                                 finished, could check for no components
            if (icp.eq.0) then
               write (*,*) 'No thermodynamic components'
               goto 998
            else if (icopt.eq.5.and.jbulk.lt.icp) then 
               write (*,*) 'All thermodynamic components must be ',
     *                     'constrained.'
               goto 998
            end if 
        
            exit 

         else if (rname.eq.blank) then 
 
            cycle 

         else if (rname.eq.'V'.or.rname.eq.'S') then

            usv = .true.

         else

            icp = icp + 1
            cname(icp) = rname
c                                 encode a graphics names for the
c                                 compositional variables, this is kind of
c                                 pointless, but it looks good.
            write (xname(icp),'(a,a,a)') 'x(',rname,')'
c                                 unblank the name
            call unblnk (xname(icp))
            if (icp.gt.k5) call error (197,r,icp,'INPUT1')

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) icont
c                                 for meemum override icont
         if (icont.ne.0.or.iam.eq.2) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, icont)    
         end if 

      end do           

      icp1 = icp + 1
      icp2 = icp + 2

      if (usv) then

         hcp = icp2
         tindex = icp1
         pindex = icp2
         cname(tindex) = 'T(K) '
         cname(pindex) = '-P(b)'

      else

         hcp = icp

      end if 
c                                 decode saturated components    
c                                 isat is the saturated component counter
      isat = 0
      io2  = 0 

      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 

      do 

         read (n1,'(a,a)') rname,strg
         if (rname.eq.blank) cycle 

         if (rname.eq.'end s') then 

            icomp = icp + isat
            exit 

         else if (rname.eq.blank) then 

            cycle 

         end if 
c                                 check for compositional constraints
         read (strg,*,err=998) jcont

         if (jcont.ne.0) then 
            jbulk = jbulk + 1
            read (strg,*,err=998) j, (dblk(i,jbulk), i = 1, icont)    
c                                 override variance flag choice, why here?
            isudo = 0    
         end if  

         isat = isat + 1
         if (isat.gt.h5) call error (15,r,i,'BUILD')
         cname(icp+isat) = rname
         if (rname.eq.'O2') io2 = isat

      end do 
c                                 decode saturated phase components
      do 
         read (n1,'(a)',end=998) rname
         if (rname.eq.'begin') exit
      end do 
c                                 ifct is the saturated phase component counter
      ifct = 0

      do 

         read (n1,'(a)') rname

         if (rname.eq.'end s') then 
            icomp = icomp + ifct
            exit 
         else if (rname.eq.blank) then 
            cycle 
         end if 
      
         ifct = ifct + 1
         if (ifct.gt.2) call error (44,r,i,'BUILD')
c                                 save the component if only one
c                                 for use in input2.
         if (ifct.eq.1) zname = rname
         cname(icomp+ifct) = rname
      end do 
c                                  decode mobile components
c                                  jmct - mobile component counter
      jmct = 0 
      ifact = 0
      jmuct = 0 

      do 

         call rdstrg (n1,nstrg,string,eof)

         if (eof) then 

            goto 998

         else if (string(1).eq.'begin') then

            cycle 

         else if (string(1).eq.'end') then

            icomp = icomp + jmct
            exit 

         else 

            read (string(1),'(a5)') rname
            jmct = jmct + 1
            if (jmct.gt.2) call error (45,r,i,'BUILD')
            cname(icomp+jmct) = rname

            if (nstrg.eq.1) then 
c                                 old format, create variable name
               write (vname(3+jmct),'(a,a)') 'mu_',rname
               imaf(jmct) = 1
               jmuct = jmuct + 1

            else 
c                                 new format
               read (string(2),'(a1)') y
               vname(3+jmct) = string(2)
               afname(jmct) = string(3)

               if (y.eq.'m') then 
c                                 chemical potential
                  imaf(jmct) = 1
                  jmuct = jmuct + 1

               else if (y.eq.'f') then 

                  imaf(jmct) = 2

               else if (y.eq.'a') then 

                  imaf(jmct) = 3

               end if 

               if (imaf(jmct).gt.1) ifact = ifact + 1 

            end if 
               
         end if 

      end do 
c                             the ifct flag can probably be set later if fluid
c                             is in the thermodynamic composition space.   
      jfct = icp + isat 
c                             jprct+1..icomp -> (jmct.ne.0) mobile components 
      jprct = icomp - jmct 
c                             excluded phases
      ixct = 0
c                             decode excluded phases
      do 
         read (n1,'(a)',end=998) name
         if (name.eq.'begin ex') exit
      end do

      do 

        read (n1,'(a)') name

         if (name.eq.'end excl') then
            exit
         else if (name.eq.blank) then
            cycle 
         end if 

         ixct = ixct + 1
         if (ixct.gt.h8) call error (13,r,i,'BUILD')
         exname(ixct) = name

      end do  
c                             solution phases:
      do 
         read (n1,'(a)',end=998) sname
         if (sname.eq.'begin solu') exit
      end do
c                             isoct - solution phase counter,
c                             io9 is a flag = 0 no solution file
      isoct = 0

      do 

         read (n1,'(a)') sname
 
         if (sname.eq.'end soluti') then 
            if (io9.eq.1) isoct = 0 
            exit 
         else if (name.eq.blank) then 
            cycle  
         end if 

         isoct = isoct + 1
         if (isoct.gt.h9) call error (25,r,i,'BUILD')
         fname(isoct) = sname

      end do  
c                             read the maximum pressure, temper-
c                             ature, xco2, u1, and u2; the minimum
c                             pressure temperature, xco2, u1, and u2;
c                             and the default pressure, temperature,
c                             xco2, and chemical
c                             potential increments use kelvins, bars and
c                             joules as units (if no mobile components
c                             enter two zeroes for each read).
      read (n1,*,err=998) vmax
      read (n1,*,err=998) vmin
      read (n1,*,err=998) dv
c                             read the default indices of the
c                             dependent, independent, and secondary
c                             independent intensive variables, p = 1,
c                             t = 2, and xco2 = 3, respectively.
      read (n1,*,err=998) (iv(i), i = 1, 5)
c                             check variable ranges are consistent,
c                             variable iv(1):
      if (icopt.ne.0.and.icopt.ne.4.and.iam.ne.2) then

         if (iv(1).eq.3.and.ifct.eq.0) call error (110,r,i,'I')

         if (iv(1).eq.3.and.ifct.eq.1) then 

            if (icopt.ne.7.and.iv(2).ne.3) call error (111,r,i,'I')

         end if 

         if (vmin(iv(1)).ge.vmax(iv(1)).and.icopt.lt.5) then 

            call error (112,r,i,'less than or equal')

         else if (vmin(iv(1)).eq.vmax(iv(1)).and.
     *            icopt.eq.5.and.icont.lt.3) then

            call error (112,r,i,'equal')

         end if 

         if (vname(iv(1)).eq.blank) call error (116,dip,i,'I')

      end if
c                             variable iv(2):
      if (iam.ne.2.and.(icopt.eq.1.or.
     *                  (icopt.eq.5.and.icont.eq.1.and..not.oned))) then

         if (iv(2).eq.3.and.ifct.eq.0) call error (110,r,i,'INPUT1')

         if (iv(2).eq.3.and.ifct.eq.1) call error (111,r,i,'INPUT1')

         if (icopt.eq.1) then 

            if (vmin(iv(2)).ge.vmax(iv(2))) call error (112,r,i,
     *                                            'less than or equal')

         else 

            if (vmin(iv(2)).eq.vmax(iv(2))) call error (112,r,i,'equal')

         end if 

         if (vname(iv(2)).eq.blank) call error (116,r,i,'INPUT1')

      end if
c                             if a chemical potential is specified as an
c                             independent variable (iv(1-3)), check if
c                             the variable is defined:
      kct = 0

      do i = 1, 3
         if (iv(i).gt.3) kct = kct + 1
      end do 
c                             identify the variable used to determine
c                             which phases lie on the left hand side
c                             of a reaction equation.
      if (icopt.eq.3) then
         ivfl = iv(1)
      else if (iv(1).eq.2.or.iv(2).eq.2) then
c                             choose T
         ivfl = 2
      else if (iv(1).eq.1.or.iv(2).eq.1) then
c                             no T, so choose P
         ivfl = 1
      else
c                             no P or T, choose independent V
         if (iv(2).ne.3) then
            ivfl = iv(2)
         else
            ivfl = iv(1)
         end if
      end if
c                             ok, now find out which variables are
c                             dummies and store the indexes of the
c                             non-dummy variables in jv.
      ipot = 0

      do i = 1, 5
c                             variables v(1) (p) and v(2) (t) are
c                             only dummies if idep is set.
         if ((iv(i).ne.idep.or.icopt.eq.7.or.icopt.eq.9).and.
     *       (iv(i).eq.1.or.iv(i).eq.2)) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variable v(3) is a dummy if ifct = 0:
         else if ((iv(i).eq.3).and.ifct.gt.0) then
            ipot = ipot+1
            jv(ipot) = iv(i)
c                             variables v(4) and v(4) are dummies if
c                             imyn = 1:
         else if (jmct.ne.0) then
            if (iv(i).eq.4) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            else if (iv(i).eq.5.and.jmct.eq.2) then
               ipot = ipot+1
               jv(ipot) = iv(i)
            end if
         end if

      end do 
c                                 if dependent variable add to jv list, could
c                                 increment ipot, but maybe it's better not to.
      if (idep.ne.0) jv(ipot+1) = idep
c                                 set convergence criteria for routine univeq
      if (icopt.le.3) call concrt

      if (icopt.ne.0) close (n1)
c                                 open files requested in input
      call fopen (n2name,prt,plt,n9name,jbulk,icp,err)
c                                 err only set for unsplt (iam.eq.14)
      if (err) return
c                                 read auxilliary input for 2d fractionation
      if (icopt.eq.9) call rdain
c                                 get runtime parameters
      if (first.or.(.not.first).and.(.not.output)) 
     *   call redop1 (first,tfname)

      goto 999
c                                 archaic error trap
998   call error (27,r,i,n2name)

999   end

      subroutine input2 (first)
c----------------------------------------------------------------------
c input2 reads the thermodynamic data file for most perplex programs, 
c a (the?) notable exception being frendly that calls the parallel 
c routine jnput2.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      character*5 mnames(k16*k17)*8

      double precision twt(k5),cst
 
      integer i,j,im, ict, k, ifer,inames, jphct, imak(k16)
 
      logical eof, good, first

      integer iff,idss,ifug
      common / cst10 /iff(2),idss(h5),ifug

      double precision ctot
      common/ cst3  /ctot(k1)

      integer iwt
      common/ cst209 /iwt

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
     
      character*5 zname
      common/ cst209a /zname

      character*5 cname
      common/ csta4 /cname(k5)

      character*8 names
      common/ cst8 /names(k1)

      character*8 name
      common/ csta6 /name

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer ic
      common/ cst42 /ic(k0)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer cl
      character cmpnt*5, dname*80
      common/ csta5 /cl(k0),cmpnt(k0),dname

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision atwt
      common/ cst45 /atwt(k0)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ixct,ifact
      common/ cst37 /ixct,ifact 

      integer aqct
      common/ cst336 /aqct

      character*8 exname,afname
      common/ cst36 /exname(h8),afname(2)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer jfct,jmct,jprct,jmuct
      common/ cst307 /jfct,jmct,jprct,jmuct

      integer imaf,idaf
      common/ cst33 /imaf(i6),idaf(i6)

      double precision mcomp
      character mknam*8
      integer nmak
      logical mksat
      common / cst333 /mcomp(k16,k0),nmak,mksat(k16),mknam(k16,k17)

      double precision mkcoef, mdqf
      integer mknum, mkind, meos
      common / cst334 /mkcoef(k16,k17),mdqf(k16,k17),mkind(k16,k17),
     *                 mknum(k16),meos(k16)

      integer make
      common / cst335 /make(k10)

      integer eos
      common/ cst303 /eos(k10)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision vnumu
      common/ cst44 /vnumu(i6,k10)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer idspe,ispec
      common/ cst19 /idspe(2),ispec
c-----------------------------------------------------------------------
c                               initialization for each data set
c                               for k10 endmembers
      do i = 1, k10
         make(i) = 0 
         names(i) = ' '
      end do
c                               for k1 phases:
      do i = 1, k1
         ikp(i) = 0
      end do 
c                               other counters and flags:
      do i = 1, h5
         isct(i) = 0
      end do 
c                               counters for bounds
      iphct = 0
      lamin = 0 
      idsin = 0 
      idfl = 0
c                               read data base header, do component
c                               transformations, read make definitions.
      call topn2 (0)
c                               general input data for main program

c                               reorder thermodynamic components
c                               if the saturated phase components are 
c                               present
      if (lopt(7)) then

         do k = 1, ispec 
                             
            do i = 1, icp

               if (cname(i).eq.cmpnt(idspe(k))) then 

                  if (i.eq.k) exit 

                  cname(i) = cname(k)

                  do j = 1, 3
                     cst = dblk(j,i)
                     dblk(j,i) = dblk(j,k) 
                     dblk(j,k) = cst
                  end do 

                  cname(k) = cmpnt(idspe(k))

                  exit            

               end if 

            end do 

         end do 

      end if  
c                              load the old cbulk array
      if (ifct.gt.0) iphct = 2
c                               identify nonzero components.
c                               initialize icout(i) = 0
      do i = 1, icmpn
         icout(i) = 0
      end do

      do i = 1, icomp

         im = 0

         do j = 1, icmpn

            if (cname(i).eq.cmpnt(j)) then 

               twt(i) = atwt(j)
               ic(i) = j
               icout(j) = 1

               do k = 1, ispec
                  if (j.eq.idspe(k)) then 
                     iff(k) = i
                     idfl = idfl + 1
                  end if 
               end do 
 
               im = 1

            end if 

         end do 
c                               write error message if a component
c                               was not found:
         if (im.eq.0) then 
            write (*,1230) cname(i), (cmpnt(k), k = 1, icmpn)
            write (*,1240)
            stop
         end if 
 
      end do 
c                                 this segment is to check if
c                                 a possible saturated phase component
c                                 has been made a mobile component,
c                                 if there is also a saturated phase
c                                 component idfl is the identity of the
c                                 mobile component otherwise idfl = 0.
      if (ifct.eq.1.and.idfl.eq.2) then

         do i = 1, ispec
            if (zname.ne.cmpnt(idspe(i))) cycle 
            idfl = i
            exit 
         end do 

      else 
         idfl = 0
      end if
c                                 load atwts in updated order
      do i = 1, icomp
         atwt(i) = twt(i)
      end do 
c                                 convert weight to molar amounts
      if (jbulk.ne.0) then 

         if (iwt.eq.1) then 
            do i = 1, jbulk
               do j = 1, 3
                  dblk(j,i) = dblk(j,i)/atwt(i)
               end do 
            end do 
         end if 

         do i = 1, jbulk
            cblk(i) = dblk(1,i)
         end do   

      end if 
c                                 get composition vectors for entities
c                                 defined by a make definition:
      call makecp (inames,mnames,first)
c                                 loop to read reference phase data for
c                                 activity/fugacity variables
      ict = 0 

      if (ifact.gt.0) then
c                                 rewind and read 'til end of header
         call eohead (n2)

         good = .false.

         do

            call getphi (name,.false.,eof)

            if (eof) then 

               write (*,1000) (afname(i),i=1,jmct)
               write (*,1010)
               stop

            end if 
c                                 now look for a match with the 
c                                 reference phase names
            do i = 1, jmct

               if (name.eq.afname(i)) then 
c                                 got a match, count
                  iphct = iphct + 1

                  ict = ict + 1

                  idaf(i) = iphct
c                                 store thermodynamic parameters:
                  call loadit (iphct,.false.,.true.)
c                                 zero the component
                  vnumu(i,iphct) = 0d0

                  if (imaf(i).eq.2) then 
c                                 if some cretin chooses fugacity, prevent
c                                 gphase from calling the EoS.   
                     eos(iphct) = ieos 

                  else if (lopt(7)) then 
c                                 check for special component names
c                                 this is necessary because loadit 
c                                 will not set isfp if ifct > 0.
                     do k = 1, ispec
                        if (name.ne.cmpnt(idspe(k))) cycle
                        eos(iphct) = 100 + k 
                        exit 
                     end do 
 
                  end if 
c                                 blank the name, this has two purposes,
c                                 it prevents problems if an entry is 
c                                 replicated in the data file, and flags
c                                 tagged entries 
                  afname(i) = ' '

                  if (ict.eq.jmct) good = .true.

                  exit 

               end if 

            end do 

            if (good) exit 

         end do 

      end if 
c                                 begin first read loop for data on
c                                 saturated components.
      if (isat.eq.0.and.ifct.eq.0) goto 40
c                                 read 'til end of header
      call eohead (n2)
c                                 loop to read real saturated
c                                 entities:
      ifer = 0

      do 

         call getphi (name,.false.,eof)

         if (eof) exit
 
         call chkphi (0,name,good)

         if (good) call sattst (ifer,good)

      end do 
c                                 loop to load made saturated entities
      do i = 1, nmak

         if (.not.mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check:
         call chkphi (2,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)
c                                 set eos flag
         ieos = meos(i)

         call sattst (ifer,good)

         if (good) then 
            make(iphct) = i
c                                 pointer used for iemod.
            imak(i) = iphct
         end if 

      end do 
c                                 check that there is data for
c                                 every fluid component.
      if (ifct.gt.0.and.ifer.ne.ifct) call error (36,r,i,'INPUT2')
c                                 check that there is one phase
c                                 for each saturation constraint
40    do i = 1, isat
         if (isct(i).lt.1) call error (15,r,i,cname(icp+i))
      end do
c                                 save endmembers that consist entirely 
c                                 of saturated phase or mobile components:
      kphct = iphct 

      if (ifct+jmct.gt.0) then 

         call eohead (n2)

         do 

            call getphi (name,.false.,eof)

            if (eof) exit

            call chkphi (4,name,good)

            if (.not.good) cycle 
c                                 reject phases already in the list
            do i = 1, kphct
               if (names(i).eq.name) then
                  good = .false.
                  exit
               end if 
            end do 

            if (.not.good) cycle             
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end do

      end if 
c                                 read data for the remaining
c                                 phases of appropriate composition.
      istct = iphct + 1
c                                 read till end of header
      call eohead (n2)
c                                 loop to load normal thermodynamic data:
      do  
    
         call getphi (name,.false.,eof)

         if (eof) exit 
c                                 check if valid phase:
         call chkphi (1,name,good)

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)
         end if 
      end do 

c                                 loop to load made entities
      do i = 1, nmak

         if (mksat(i)) cycle
c                                 load make data 
         do j = 1, icmpn
            comp(j) = mcomp(i,j)
         end do 

         name = mknam(i,mknum(i)+1)
c                                 redundant check, but makes ctot.
         call chkphi (3,name,good)
c                               
         if (.not.good) call error (57,comp(1),iphct,name)

         iphct = iphct + 1
         ctot(iphct) = tot
c                                 set ieos flag to that of the first
c                                 real entity in the make definition
         ieos = meos(i)

         call loadit (iphct,.true.,.true.)

         make(iphct) = i
c                                 pointer used for iemod.
         imak(i) = iphct

      end do 
c                                 load thermodynamic data for solute species, at this
c                                 point iphct points to the last real entity, save
c                                 this value and restore it later.
      jphct = iphct
c                                 read header
      call eohead (n2)
c                                 loop to load solute data:
      do  
    
         call getphi (name,.false.,eof)

         if (eof) exit
c                                 skip non-solute standard state data
         if (ieos.ne.15.and.ieos.ne.16) cycle
c                                 check if valid species:
         call chkphi (1,name,good)

         if (good) then 
c                                 acceptable data, count the phase:
            iphct = iphct + 1
c                                 for normalized composition, probably
c                                 con't need this, but could be used to
c                                 save molar wt or something like that:
            ctot(iphct) = tot
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.true.)

         end if 

      end do
c                                 later aqct is used to determine
c                                 the indices of aqueous species used
c                                 in the calculation. 
      aqct = iphct
c                                 get/save data for makes, this
c                                 data is saved in the arrays thermo
c                                 and cp by loadit, but are not counted,
c                                 i.e., the counters ipoint and iphct
c                                 are reset. soload will then load the
c                                 cp array over the values loaded here,
c                                 but thermo should not be affected. gmake
c                                 then gets the data using the array 
c                                 mkind. the names array will also be 
c                                 overwritten.

c                                 read header
      call eohead (n2)

      do 

         call getphi (name,.true.,eof)

         if (eof) exit

         do i = 1, inames

            if (name.ne.mnames(i)) cycle
c                                 matched a name
            iphct = iphct + 1
c                                 store thermodynamic parameters:
            call loadit (iphct,.false.,.false.)

         end do

      end do 

      do i = 1, nmak
c                                remake pointer array for makes 
         do j = 1, mknum(i)
            do k = jphct + 1, iphct
               if (names(k).ne.mknam(i,j)) cycle
               mkind(i,j) = k
            end do
         end do 
      end do  
c                                reset ipoint counter, but do not 
c                                reset iphct, because the compositions
c                                of the make phases are necessary for
c                                chemical potential variables.

c                                really? then why was it reset here? this
c                                IS going to wipe out the compositions of
c                                aqueous species, so iphct is now set to
c                                aqct (formerly it was jphct). Jan 6, 2017.
      iphct = aqct
      ipoint = jphct

      do i = 1, nmak
c                                make an iemod flag for made
c                                endmembers:   
         do j = 1, mknum(i)
            if (iemod(mkind(i,j)).eq.0) exit 
         end do 

         if (j.le.mknum(i)) cycle

         iemod(imak(i)) = iemod(mkind(i,1))

      end do 

1000  format ('**error ver007** at least one of the reference ',
     *        'phases:',/,5(a,1x))
1010  format ('needed to define an independent fugacity/activity ',
     *    'variable is missing from the',/,'thermodynamic data file',/)
1230  format ('**error ver013** ',a,' is an incorrect component'
     *       ,' name, valid names are:',/,12(1x,a))
1240  format ('check for upper/lower case matches or extra blanks',/)

      close (n2)

      end

      subroutine setvr0 (i,j)
c--------------------------------------------------------------------
c setvr1 computes nodal variables for node ij, three cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont
c----------------------------------------------------------------------

      if (icont.eq.1) then 

         v(iv1) = vmin(iv1) + (i-1)*dv(iv1)
         v(iv2) = vmin(iv2) + (j-1)*dv(iv2)
         call incdp0

      else if (icont.eq.2) then 

         v(iv1) = vmin(iv1) + (j-1)*dv(iv1)
         call incdep (iv1)

         cx(1) =  (i-1)*dvr(1)
         call setblk 

      else 

         cx(1) = (i-1) * dvr(1)
         cx(2) = (j-1) * dvr(2)
         call setblk

      end if 

      end

      subroutine setblk
c-----------------------------------------------------------------------
c for gridded minimization setblk computes the bulk composition
c and initializes the arrays for lpopt.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x0

      integer i,j

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      x0 = 1d0

      if (lopt(1)) then 
c                                 closed composition
         do j = 1, icont-1
            x0 = x0 - cx(j)
         end do 

      end if 

      do j = 1, jbulk
         cblk(j) = x0*dblk(1,j)
      end do 
         
      do j = 1, jbulk
         do i = 2, icont 
            cblk(j) = cblk(j) + cx(i-1)*dblk(i,j)
         end do 
      end do
c                                 modify cblk here to change the 
c                                 composition before minimization.
      ctotal = 0d0 
c                                 get total moles to compute mole fractions             
      do i = 1, hcp
         ctotal = ctotal + cblk(i)
      end do

      do i = 1, hcp 
         b(i) = cblk(i)/ctotal
      end do

      end 

      subroutine setvar 
c--------------------------------------------------------------------
c setvar initializes the variables for gridded minimization, three
c cases:

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j

      double precision rloopy,rloopx

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio
c----------------------------------------------------------------------

      rloopy = dfloat(loopy-1)
      rloopx = dfloat(loopx-1)
c                                 for 1d calculations
      if (rloopx.eq.0) rloopx = rloopy

      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do

      call incdp0

      if (icopt.eq.7.and.fileio) then 
c                                using nodal coordinate system
         dvr(1) = 1

      else if (icont.eq.1) then 
c                                v(iv1) on x, v(iv2) on y
         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopx
         dvr(1) = dv(iv1)

         dv(iv2) = (vmax(iv2) - vmin(iv2))/rloopy
         dvr(2) = dv(iv2)

      else if (icont.eq.2) then 
c                               composition is on x, v(iv1) on y
         dvr(1) = 1d0/rloopx

         dv(iv1) = (vmax(iv1) - vmin(iv1))/rloopy
         dvr(2) = dv(iv1)

      else 
c                                compositions on both axes
         dvr(1) = 1d0/rloopx
         dvr(2) = 1d0/rloopy 
         cx(1) = 0d0
         cx(2) = 0d0

      end if 
c                                set the bulk composition:
      do j = 1, jbulk
         cblk(j) = dblk(1,j)
      end do 

      end 

      subroutine inipot 
c--------------------------------------------------------------------
c setvar initializes the independent potential variables to their 
c minimum values
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
c----------------------------------------------------------------------
c                                 initialize potentials
      do i = 1, ipot
         v(jv(i)) = vmin(jv(i))
      end do 
c                                 set dependent potential, if it exists
      call incdp0

      end 

      subroutine getcmp (jd,id,ids)
c-----------------------------------------------------------------------
c getcmp gets the composition of pseudocompund id, where:
c  if ids < 0, -ids points to the composition of a true compound in array cp
c  if ids > 0, id points to the composition of a solution defined in terms
c              on endmember fractions defined and saved by routine resub
c              in array zcoor.
c the composition is saved in arrays cp3 and x3, entry jd
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,id,jd,ids
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 single site solution coordinates:
      integer jend
      common/ cxt23 /jend(h9,m4)
c                                 refined compositions and solution 
c                                 pointer
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp0
      common/ cst71 /cp0(k0,k5)

      integer iam
      common/ cst4 /iam

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------

      kkp(jd) = ids
      cptot(jd) = 0d0

      if (ids.lt.0) then 
c                                 simple compounds
         if (iam.ne.5) then
c                                 all programs except frendly 
            do i = 1, icomp
               cp3(i,jd) = cp(i,-ids)
            end do 
c                                 check if it's a solution endmember
            if (ikp(-ids).ne.0) call endcp (jd,-ids,ikp(-ids))
   
         else 
c                                 frendly 
            do i = 1, k0
               cp3(i,jd) = cp0(i,-ids)
            end do 

         end if 

      else
c                                 solutions, initialize
         do i = 1, icomp
            cp3(i,jd) = 0d0
         end do
c                                 get the x(i,j) coordinates for the
c                                 composition from the zcoor array,
c                                 this routine also saves a copy of the
c                                 xcoordinates in x3(jd,i,j)
         call getxz (jd,id,ids)
c                                 convert the x(i,j) coordinates to the
c                                 geometric y coordinates
         call xtoy (ids)

         if (lrecip(ids)) then
c                                 get the p' coordinates (amounts of 
c                                 the independent endmembers)     
            call getpp (ids) 

            do i = 1, lstot(ids)
               do j = 1, icomp 
                  cp3(j,jd) = cp3(j,jd) + p0a(i) * cp(j,jend(ids,2+i))
               end do 
            end do          

         else 
c                                 solutions with no dependent endmembers:
c                                 y coordinates used to compute the composition
            do i = 1, mstot(ids)
               do j = 1, icomp
                  cp3(j,jd) = cp3(j,jd) + y(i) * cp(j,jend(ids,2+i))
               end do
            end do

         end if 

      end if 

      do i = 1, icp
         cptot(jd) = cptot(jd) + cp3(i,jd)
      end do 

      end 

      subroutine inblnk (text,char)
c----------------------------------------------------------------------
c inblnk - scan text to last '/' or '\' and insert char after.
 
c     text - character string 
c----------------------------------------------------------------------
      implicit none

      integer i, nchar
 
      character text*(*), bitsy(400)*1, char*1 
c----------------------------------------------------------------------
      nchar = len(text) 
      read (text,1000) (bitsy(i), i = 1, nchar)
c                                 scan for blanks:

      do i = nchar,1,-1
c                                 this line may cause problems
c                                 on some operating systems that 
c                                 recognize the backslash as an escape
c                                 character.
         if (bitsy(i).eq.'/') goto 10
         bitsy(i+1) = bitsy(i)
      end do 

      i = 0

10    bitsy(i+1) = char

      write (text,1000) (bitsy(i), i = 1, nchar)
 
1000  format (400a1)
      end

      subroutine matchj (unnown,itis)
c----------------------------------------------------------------------
 
c matchj - subroutine to determine if the string unnown is a valid
c          solution or compound name.
 
c   itis = -id if compound
c   itis = ikp if solution 
c   itis = 0 if invalid
c----------------------------------------------------------------------
      implicit none

      integer i, itis
 
      character*10 unnown
 
      include 'perplex_parameters.h'
 
      integer isoct
      common/ cst79 /isoct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character names*8
      common/ cst8 /names(k1)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c---------------------------------------------------------------------- 
 
      itis = 0

      do i = 1, isoct
         if (unnown.eq.fname(i)) then
             itis = i
             goto 99
         end if
      end do

      do i = 1, iphct
         if (unnown.eq.names(i)) then
            itis = -i
            goto 99
         end if
      end do 

99    end

      subroutine maktit 
c-----------------------------------------------------------------------
c create a title for graphics output, the title consists of the 
c calculation title + saturation hierarchy (provided one is 
c specified) and is the first two elements of title (csta8).
c if icopt = 1 or 3, also adds a blurb about reaction convention.

c title is max 3 lines, but four lines are written to be consistent
c with old plot file formats written by frendly, pt2curv etc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      character*162 title
      common/ csta8 /title(4)

      character*8 vname,xname     
      common/ csta2  /xname(k5),vname(l2)

      integer ivfl
      common/ cst102 /ivfl

      character*5 cname
      common/ csta4 /cname(k5)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c-----------------------------------------------------------------------
      do i = 2, 4
         title(i) = ' '
      end do                              
c                               saturated and buffered component names:
      if (isat.gt.0) then 
         write (title(2),1070) (cname(i+icp), i= 1, isat)
      else 
         write (title(2),1000) ' '
      end if 
c                                 reaction convention
      if (icopt.eq.1.or.icopt.eq.3) write (title(3),1080) vname(ivfl)

      do i = 1, 3
         call deblnk (title(i))
      end do 

1000  format (a)
1070  format ('Component saturation hierarchy: ',7(a,1x))
1080  format ('Reaction equations are written with the high ',
     *         a,'assemblage to the right of the = sign')

      end

      subroutine rdain
c-----------------------------------------------------------------------
c a subprogram to read auxilliary input file for 2d fractionation 
c calculations, called by VERTEX and WERAMI
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxbox,maxlay
      parameter (maxbox=1760,maxlay=6) 

      integer i,j,k,ier

      double precision zlayer

      character*100 name

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer ipot,jv,iv
      common/ cst24 /ipot,jv(l2),iv(l2)

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc 

      integer gloopy,ilay,irep
      double precision a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk
      common/ cst66 /a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk(maxlay,k5),gloopy,ilay,irep(maxlay)

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio

      character*100 cfname
      common/ cst227 /cfname
c-----------------------------------------------------------------------
c                                 look for input data from a file 
c                                 of type aux
      call mertxt (name,prject,'.aux',0)
      open (n8,file=name,status='old',iostat=ier)

      if (ier.ne.0) call error (51,zbox,gloopy,name) 
c                                 set the number of independent variables
c                                 to 1 (the independent path variable must
c                                 be variable jv(1), and the dependent path
c                                 variable must be jv(2), the path variables
c                                 can only be pressure and temperature
      ipot = 1
c                                 in old versions this was the number of steps 
c                                 along the path, in 07 loopy is now set via
c                                 jlow/1dpath (above). gloopy is not used as a 
c                                 flag if top layer composition is to be refreshed
c                                 after each step (gloopy=999).
      read (n8,*) gloopy
c                                 thickness of a box in column
      read (n8,*) zbox 
c                                 gradient in variable jv(1) with z, jv(1)
c                                 is the independent variable, for subduction
c                                 this is logically pressure, i.e.,
c                                 dp(bar)/dz(m)
      read (n8,*) dv1dz 
c                                 now we need a path function for the dependent
c                                 variable, here we take a function defined in
c                                 terms of the absolute depth of the top of the
c                                 column (z0) and the relative depth (dz) within
c                                 the column
c                                 
c                                 v2 = a(z0)*dz^2 + b(z0)*dz + c(z0)

c                                 e.g., T(K) =  a(z0)*dz^2 + b(z0)*dz + c(z0)

c                                 where a(z0) = a0 + a1*z0 + a2*z0^2 + a3*z0^3 + ...
c                                 b(z0) = b0 + b1*z0 + b2*z0^2 + b3*z0^3 + ...
c                                 c(z0) = c0 + c1*z0 + c2*z0^2 + c3*z0^3 + ...
      read (n8,*) a0, a1, a2, a3
      read (n8,*) b0, b1, b2, b3
      read (n8,*) c0, c1, c2, c3
c                                 get the initial global composition array
c                                 consisting of ibox compositions defined 
c                                 in terms of icp components. this read
c                                 statement assumes that H2O an CO2 (if 
c                                 thermodynamic components) are the 1st and
c                                 2nd components (if present). 
      ilay = 0
      ncol = 0
      vmax(2) = 0d0
      vmin(2) = 0d0
c                                 number of nodes with appended composition
c                                 end of data indicated by zero 
      do 

         read (n8,*) zlayer

         if (zlayer.eq.0) exit 

         ilay = ilay + 1

         if (ilay.eq.maxlay) then 
            write (*,*) 'increase maxlay in routine FRAC2D'
            stop
         end if 

         read (n8,*) (iblk(ilay,i),i=1,icp)

         irep(ilay) = idint(zlayer/zbox)

         ncol = ncol + irep(ilay)

         if (ncol.gt.maxbox) then
            write (*,*) 'increase maxbox in routine DUMMY1'
            stop
         end if 
c                                 set the y coodinate to depth below top
         vmin(2) = vmin(2) - irep(ilay)*zbox

      end do 

      close (n8)
c                                 two cases, file input or analytical
      if (fileio) then 
c                                 file input of nodal p-t coordinates
         open (n8,file=cfname,status='old',iostat=ier)
c                                 read header info
         read (n8,*) i, nrow

         if (ncol*nrow.gt.k2) then 
            write (*,'(/,a,i6,a,i6,a)') 
     *     '**error ** too many coordinates, nodes*columns>k2',i*nrow,
     *     'increase k2 (',k2,')'

            stop 

         else if (i.ne.ncol) then 
            write (*,'(/,a,i4,a,a,/,a,i4,a)') 
     *     '** error ** the number of nodes in a column (',i,
     *     ') specified in: ',cfname,'must equal the',
     *     'number of lithological nodes (',ncol,
     *     ')specified in the aux file.'

           stop 

         end if 

         do i = 1, nrow

            k = (i-1) * ncol

            do j = 1, ncol
               read (n8,*) vn(k+j,1),vn(k+j,2)
            end do 

         end do

         close (n8)

      end if 

      end 

      subroutine fr2dpt (p0,dz)
c----------------------------------------------------------------------
c subroutine to set p-t variables from i-j coordinates in 2d-fractionation
c calculations, called by VERTEX and WERAMI
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer maxlay,i,j

      parameter (maxlay=6) 

      double precision p0, z0, dz

      integer gloopy,ilay,irep
      double precision a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk
      common/ cst66 /a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,dv1dz,
     *               zbox,iblk(maxlay,k5),gloopy,ilay,irep(maxlay)

      integer irct,ird
      double precision vn
      common/ cst31 /vn(k2,k7),irct,ird

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio 

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)  

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
c----------------------------------------------------------------------
      if (fileio) then 
c                                convert p0-dz coordinate to nodal 
c                                values
         i = (p0 - vmin(1))/dv(1) + 1
         j = ncol - int(dz/zbox)

         v(1) = vn((i-1)*ncol + j, 1)
         v(2) = vn((i-1)*ncol + j, 2)

      else 
c                                 convert to depth at top of column
         z0 = p0/dv1dz
c                                 set the independent variable
         v(1) = p0 + dz * dv1dz   
c                                 set the dependent variable
         v(2) = (a0 + a1*z0 + a2*z0**2 + a3*z0**3)*dz**2 
     *        + (b0 + b1*z0 + b2*z0**2 + b3*z0**3)*dz 
     *        +  c0 + c1*z0 + c2*z0**2 + c3*z0**3
      end if 
                       
      end 
   
      subroutine getpp (id)
c-----------------------------------------------------------------------
c getpp computes the amounts of the indepdendent edmembers of a reciprocal
c solution in terms of the disordered endmembers (i.e., the p coordinates
c corrected for the amounts of the ordered species if present [ksmod=8]).
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
c                                  first convert the istot disordered
c                                  endmember coordinates to the 
c                                  kstot + nord p0 coordinates
      call y2p0 (id) 
c                                  decompose ordered species
      if (nord(id).gt.0) call p0dord (id)

      end

      subroutine fopen (n2name,prt,plt,n9name,jbulk,icp,err)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical first, err

      integer ierr,jbulk,icp
 
      character*100 blank*1,n2name,prt*3,plt*3,name,n9name

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iam
      common/ cst4 /iam

      save first,blank

      data first,blank/.true.,' '/
c----------------------------------------------------------------------
c                                 open thermodynamic data file
      call fopen2 (0,n2name) 

      if (iam.gt.2.and.iam.ne.13) then 
c                                 perplex programs other than
c                                 meemum and vertex only open
c                                 exisiting files:
c                                 open the plot file
         call mertxt (name,prject,'.plt',0)
         open (n4, file = name, iostat = ierr, status = 'old')
         if (ierr.ne.0) then
            if (iam.eq.14) then
               err = .true.
               return
            else 
               call error (122,0d0,n4,name)
            end if
         end if 
c                                 open solution model file
         if (n9name.ne.blank) then
            io9 = 0 
            open (n9,file = n9name,iostat = ierr,status = 'old')
            if (ierr.ne.0) call error (120,0d0,n9,n9name)
         else
            io9 = 1
         end if
c                                 open assemblage file
         call mertxt (name,prject,'.blk',0)
         open (n5, file = name, iostat = ierr, status = 'old')
         if (ierr.ne.0) call error (122,0d0,n4,name)

         return

      end if 

      if (first) then 
         call mertxt (name,prject,'.dat',0)
         write (*,1160) name
         write (*,1170) n2name
      end if 

      if (n9name.ne.blank) then

         io9 = 0 
c                                 open solution model file
         open (n9,file = n9name,iostat = ierr,status = 'old')
         if (ierr.ne.0) call error (120,0d0,n9,n9name)

         if (first) write (*,1210) n9name

      else

         io9 = 1
         if (first) write (*,1210) 'not requested'

      end if
c                                 open print/plot files if requested
      if (prt.ne.blank.and.prt.ne.'no_'.and.iam.ne.13) then 
         io3 = 0 
         call mertxt (name,prject,'.prn',0)
         open (n3, file = name)
      else
         io3 = 1
         name = 'none requested'
      end if

      if (first.and.iam.ne.13) write (*,1180) name

      if (plt.ne.blank.and.plt.ne.'no_') then
         io4 = 0
         call mertxt (name,prject,'.plt',0)
         open (n4, file = name)
      else
         io4 = 1
         name = 'none requested'
      end if

      if (first) write (*,1190) name

      if (jbulk.ge.icp.and.io4.ne.1) then
c                                 create special plot output file
         call mertxt (name,prject,'.blk',0)
         open (n5, file = name)
         if (first) write (*,1220) name

      else if (jbulk.ge.icp.and.io4.eq.1) then 

         if (first) write (*,1220) 'none requested'

      end if

      first = .false.

1160  format (/,'Reading problem definition from file: ',a)
1170  format ('Reading thermodynamic data from file: ',a)
1180  format ('Writing print output to file: ',a)
1190  format ('Writing plot output to file: ',a)
1210  format ('Reading solution models from file: ',a)
1220  format ('Writing bulk composition plot output to file: ',a)

      end 

      subroutine outgrd (loopx,loopy,jinc)
c----------------------------------------------------------------------
c output grid data to the plot file
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer loopx,loopy,jinc,i,j,jst,kst,kd,ltic,iend

      character string*240

      integer igrd
      common/ cst311 /igrd(l7,l7)

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(240)
c----------------------------------------------------------------------

      write (n4,*) loopx, loopy, jinc
c                                 fill in grid
      do i = 1, loopx

         if (i.ne.1.and.igrd(i,1).eq.0) igrd(i,1) = igrd(i-1,1)

         kst = 1

20       jst = kst
         if (i.ne.1.and.igrd(i,jst).eq.0) igrd(i,jst) = igrd(i-1,jst)
         kd = igrd(i,jst)
         ltic = -1

         do j = jst, loopy

            if (i.ne.1.and.igrd(i,j).eq.0) igrd(i,j) = igrd(i-1,j)

            if (igrd(i,j).eq.0.or.igrd(i,j).eq.kd) then
               ltic = ltic + 1
               if (j.eq.loopy) write (n4,*) ltic,kd
            else 
               write (n4,*) ltic,kd
               kst = j
               goto 20 
            end if 
         end do 
      end do         
c                                 write assemblage list
      write (n4,*) iasct

      do i = 1, iasct
         write (n4,*) iavar(1,i),iavar(2,i),iavar(3,i)
         write (n4,*) (idasls(j,i), j = 1, iavar(3,i))
      end do 
c                                 write assemblages to print file
      if (io3.eq.0) then 

         write (n3,'(/,1x,a,a,/)') 'Stable assemblages identified ',
     *                         'by assemblage index:'
         do i = 1, iasct
            call psbtxt (i,string,iend)
            write (n3,'(i4,a,240a)') i,' - ',(chars(j), j = 1, length)
         end do       

      end if 

      end 

      subroutine psbtxt (id,string,iend)
c----------------------------------------------------------------------
c subprogram to write a text labels for bulk composition output 
c id identifies the assemblage

      implicit none

      include 'perplex_parameters.h'

      character string*(*), pname*14

      integer i, j, ist, iend, id, np, ntot, ids

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer length,iblank,icom
      character chars*1
      common/ cst51 /length,iblank,icom,chars(240)
c----------------------------------------------------------------------
      iend = 0

      string = ' '

      ist = 1
      np = iavar(1,id)
      ntot = iavar(3,id)

      do i = 1, 240
         chars(i) = ' '
      end do
c                                 first solution names:
      do i = 1, ntot
             
         ids = idasls(i,id)

         call getnam (pname,ids) 

         ist = iend + 1
         iend = ist + 14
         read (pname,'(400a1)') (chars(j),j=ist,iend)

         call ftext (ist,iend)

      end do 

      write (string,'(400a1)') (chars(j),j=1,iend) 

      length = iend

      end 
