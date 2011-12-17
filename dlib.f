c routines common to psect and reader

      subroutine bplinp
c-----------------------------------------------------------------------
c read the b-plot file that contains the information on the assemblages
c stable at each grid node
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer nco, jxco, kxco, i, j, ids, ier
c                                 -------------------------------------
c                                 global variables
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 global assemblage data
      integer icog,jcog
      common/ cxt17 /icog(k2),jcog(k2)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)

      double precision bg
      common/ cxt19 /bg(k5,k2)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jtest,jpot
      common/ debug /jtest,jpot

      double precision mus
      common/ cst48 /mus(k8,k2)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 assemblage counter
      ibulk = 0
c                                 pointer to solution compositional coordinates
      jxco = 0 
      kxco = 0

      do 

         ibulk = ibulk + 1

         if (ibulk.gt.k2) call error (183,0d0,k2,'BLINP1')

         read (n5,*,end=99) icog(ibulk),jcog(ibulk),iap(ibulk)

         ias = iap(ibulk)
c                                phase molar amounts
         read (n5,*) (bg(i,ibulk),i=1,iavar(3,ias))

         icoor(ibulk) = jxco

         do i = 1, iavar(1,ias)

            ids = idasls(i,ias)

            nco = 0 

            do j = 1, istg(ids)
               nco = nco + ispg(ids,j)
            end do       

            jxco = jxco + 1
            kxco = jxco + nco - 1

            if (kxco.gt.k18) call error (61,0d0,k18,'BPLINP')

            read (n5,*) (xcoor(j), j = jxco, kxco)
         
            jxco = kxco

         end do 

         jxco = kxco  
c                                 read mu's if available
         if (jpot.ne.1) read (n5,*,iostat=ier) 
     *                       (mus(i,ibulk), i = 1, hcp)

         if (ier.ne.0) then 
c                                 if error on read most probably its
c                                 because of NaN's for the chemical 
c                                 potentials
            do i = 1, hcp
               mus(i,ibulk) = nopt(7)
            end do 
 
         end if 

      end do                 

99    end  

      subroutine getvar  
c--------------------------------------------------------------------
c getvar makes a list of variables to be used for i/o:

c if icopt = 10 -> using nodal coordinates else, 

c if icopt =  9 -> using 2d frac coordinates else:

c one-dimensional diagram (oned = .true.) then 

c the vertical (real) axis is variable iv(2), the horizontal axis
c is dummy.

c two-dimensional diagram (oned = .false.) then 

c icont = 1 -> independent variables are the 1st and 2nd potentials
c icont = 2 -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3 -> independent variables are compositional variables

c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      double precision vmax,vmin,dv
      common/ cst9  /vmax(l2),vmin(l2),dv(l2)

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)   

      character vname*8, xname*8
      common / csta2 /xname(k5),vname(l2)  

      logical oned
      common/ cst82 /oned

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc
c----------------------------------------------------------------------

      if (icopt.lt.9) then 

         jvar = ipot

         if (idep.ne.0) jvar = ipot + 1

         if (icont.eq.1) then 

            do i = 1, jvar
               vnm(i) = vname(jv(i))
               vmx(i) = vmax(jv(i))
               vmn(i) = vmin(jv(i))
               var(i) = vmin(jv(i))
            end do   

         else 

            if (icont.eq.2) then 

               jvar = jvar + 1

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               do i = 2, jvar
                  vnm(i) = vname(jv(i-1))
                  vmx(i) = vmax(jv(i-1))
                  vmn(i) = vmin(jv(i-1))
                  var(i) = vmin(jv(i-1))
               end do   
 
            else 

               jvar = jvar + 2

               vnm(1) = ' X(C1)  '
               vmx(1) = 1d0
               vmn(1) = 0d0

               vnm(2) = ' X(C2)  '
               vmx(2) = 1d0
               vmn(2) = 0d0

               do i = 3, jvar
                  vnm(i) = vname(jv(i-2))
                  vmx(i) = vmax(jv(i-2))
                  vmn(i) = vmin(jv(i-2))
                  var(i) = vmin(jv(i-2))
               end do   

            end if 

         end if 

         if (oned) then 
c                                 make a fake y-axis for 1-d plots
            vmx(2) = 1d0
            vmn(2) = 0d0 

         end if

      else if (icopt.eq.9) then 
c                                 2d fractionation
         vnm(1) = 'P0(bar)   '
         vnm(2) = 'DZ(m)     '    
c                                 switch loopx and loopy     
         do i = 1, 2
            vmx(i) = vmax(jv(i))
            vmn(i) = vmin(jv(i))
            var(i) = vmin(jv(i))
         end do 

         jvar = 4

         do i = 3, 4
            vnm(i) = vname(jv(i-2))
         end do 

      else 
c                                 icopt = 10:
c                                 using nodal coordinates as the x axis
         vnm(1) = 'node #'
         vmn(1) = 1
         vmx(2) = 1d0
         vmn(2) = 0d0 
         vmx(1) = loopy 
         oned = .true.
         
         jvar = ipot + 1

         do i = 2, jvar
            vnm(i) = vname(jv(i-1))
         end do   

      end if 

      end

      subroutine fopen (n2name,pr,pl,n9name,jbulk,icp)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer ierr,jbulk,icp
 
      character*100 blank*1,n2name,pr*3,pl*3,name,n9name

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      data blank/' '/
c----------------------------------------------------------------------
c                                 open thermodynamic data file
      call fopen2 (0,n2name) 
c                                 open the plot file
      call mertxt (name,prject,'.plt',0)
      open (n4, file = name, iostat = ierr, status = 'old')
      if (ierr.ne.0) call error (122,0d0,n4,name)
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

      end 

      subroutine plinp
c---------------------------------------------------------------------- 
c plinp - subroutine to read assemblage info for gridded min calculations.
c if icopt = 7 and fileio also reads nodal coordinates.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, jst, irep, kd, jend, ier

      logical count

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer igrd
      common/ cst311/igrd(l7,l7)

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer idasls,iavar,iasct,ias
      common/ cst75 /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer jcont
      common/ cst315 /jcont

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iap,ibulk
      common/ cst74 /iap(k2),ibulk

      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)

      double precision vip
      common/ cst28 /vip(l2,k2)

      character*100 cfname
      common/ cst227 /cfname

      logical fileio
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab

      integer idsol,nrep,nph
      common/ cst38/idsol(k5,k3),nrep(k5,k3),nph(k3)
c----------------------------------------------------------------------
      if (jcont.ne.0) then 
c                                 turn interpolation off for
c                                 fractionation calcs or compositional
c                                 variables, this could be optional.
         iopt(4) = 0
         write (*,3000) 

      end if 
c                                 top of plot file
      read (n4,*) loopx, loopy, jinc
c                                 decompress the grid data
      do i = 1, loopx
         jst = 1
         do while (jst.le.loopy)
            read (n4,*) irep, kd
            if (kd.eq.0) write (*,*) 'bad un at i, j',i,j
            jend = jst + irep 
            do j = jst, jend
               igrd(i,j) = kd
            end do 
            jst = jend + 1
         end do 
      end do 
c                                 read assemblages
      read (n4,*) iasct

      istab = 0 

      do i = 1, iasct
         read (n4,*) iavar(1,i),iavar(2,i),iavar(3,i)
         read (n4,*) (idasls(j,i), j = 1, iavar(3,i))
c                                 make a cumulative list of stable phases
c                                 first get the number of occurrences of 
c                                 each phase in the assemblage
         nph(i) = 0
         do j = 1, k5
            idsol(j,i) = 0
            nrep(j,i) = 0
         end do 

         do j = 1, iavar(3,i) 

            count = .true.

            if (j.le.iavar(1,i)) then 

               do k = 1, nph(i)
                  if (idsol(k,i).eq.idasls(j,i)) then 
                     count = .false.
                     nrep(k,i) = nrep(k,i) + 1
                     exit 
                  end if 
               end do
 
            end if 

            if (count) then
               nph(i) = nph(i) + 1
               idsol(nph(i),i) = idasls(j,i)
               nrep(nph(i),i) = 1
            end if

         end do  
c                                 make an array in which each id 
c                                 occurs only once
         
c                                 next compare to the existing list
         do k = 1, nph(i)  

            count = .true.

            do j = 1, istab

               if (idsol(k,i).eq.idstab(j)) then 
                  if (nrep(k,i).gt.nstab(j)) nstab(j) = nrep(k,i)
                  count = .false.
                  exit
               end if 

            end do 

            if (count) then 
               istab = istab + 1
               if (istab.gt.k10) call error (999,0d0,istab,'ISTAB ')
               nstab(istab) = nrep(k,i)
               idstab(istab) = idsol(k,i)
            end if 

         end do 

      end do 
c                                 make the "null" assemblage
      iap(k2) = k3
      iavar(1,k3) = 0
      iavar(2,k3) = 0 
      iavar(3,k3) = 0 

      if (icopt.eq.7.and.fileio) then 
c                                 if coodinates from a file, read
c                                 coordinate file.
         open (n8,file=cfname,status='old',iostat=ier)
         if (ier.ne.0) call error (6,vip(1,1),i,cfname)
         if (loopy.gt.k2) call error (1,vip(1,1),loopy,'k2')
         do j = 1, loopy
            read (n8,*,iostat=ier) (vip(i,j), i = 1, ipot)
            if (ier.ne.0) then 
               write (*,1000) cfname
               stop
            end if 
         end do 
         close (n8)

      end if 

1000  format (/,'**error ver635** Coordinate file ',a,/,
     *       'is inconsistent with plot file, re-run VERTEX.',/)
3000  format (/,'**warning ver636** For this computational mode ',
     *          'interpolation of physical',/,'properties has been ',
     *          'turned OFF.',/)
      end

      subroutine  mkcomp (jcomp,ids)
c----------------------------------------------------------------
c mkcomp makes the jcomp'th user defined compositional variable
c the first k5 compositions are reserved for chsprp, the remaining
c k5 compositions are for solvus testing

c the solution ids is associated with the composition.

c   jcx  - the number of components to define the numerator of
c          the composition.
c   jcx1 - the number of components to define the denominator of
c          the composition.
c   icps - the indices of the components (1..jcx,jcx+1...jcx1).
c   rcps - the cofficients on the compenents as indexed by icps.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*5 y*1, units*15, text*195

      integer jcomp, ier, i, ids

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer icps, jcx, jcx1, kds
      logical stol, savg
      double precision rcps
      common/ comps /rcps(k7,2*k5),icps(k7,2*k5),jcx(2*k5),jcx1(2*k5),
     *               kds(2*k5),stol(h9),savg(h9)
c----------------------------------------------------------------------
c                                choose units for composition
      if (iopt(2).eq.0) then
         units = 'mole proportion'
      else 
         units = 'weight fraction'
      end if 
c                                get the composition to be contoured
10    write (*,1100) units
  
      do 

         write (*,1030) 'numerator',k5+1
         read (*,*,iostat=ier) jcx(jcomp)

         if (ier.ne.0.or.jcx(jcomp).lt.1) then
            write (*,1020)
            cycle 
         end if 

         exit

      end do 
c                                define the numerator
      do 

         write (*,1040) 'numerator'
         write (*,1010) (i,cname(i),i = 1, icomp)
         read (*,*,iostat=ier) (icps(i,jcomp),rcps(i,jcomp), 
     *                                     i = 1, jcx(jcomp))
         do i = 1, jcx(jcomp)
            if (icps(i,jcomp).lt.1.or.icps(i,jcomp).gt.icomp) then
               ier = 1
               exit 
            end if 
         end do 

         if (ier.ne.0) then
            write (*,1020)
            cycle 
         end if 

         exit 

      end do  
c                                define the denominator
  
      do 

         write (*,1030) 'denominator',k5+1-jcx(jcomp)
         write (*,1140)
         read (*,*,iostat=ier) jcx1(jcomp)

         if (ier.ne.0.or.jcx1(jcomp).lt.0) then
            write (*,1020)
            cycle 
         end if 
 
         jcx1(jcomp) = jcx(jcomp) + jcx1(jcomp)
        
         exit 

      end do 

      if (jcx1(jcomp).gt.jcx(jcomp)) then 

         do 

            write (*,1040) 'denominator'
            write (*,1010) (i,cname(i),i = 1, icomp)
            read (*,*,iostat=ier) (icps(i,jcomp),rcps(i,jcomp), 
     *                                 i = jcx(jcomp)+1, jcx1(jcomp))

            do i = jcx(jcomp)+1, jcx1(jcomp)
               if (icps(i,jcomp).lt.1.or.icps(i,jcomp).gt.icomp) then
                  ier = 1
                  exit 
               end if 
            end do 

            if (ier.ne.0) then
               write (*,1020)
               cycle 
            end if 
c                                show the user the composition: 
            write (*,1070)           
            write (text,1130) (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                      i = 1, jcx(jcomp))
            call deblnk (text)
            write (*,1150) text 
            write (*,*) '   divided by '

            write (text,1130) (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                 i = jcx(jcomp)+1, jcx1(jcomp))
            call deblnk (text)
            write (*,1150) text 

            exit 

         end do 

      else 

         write (text,1130) (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                    i = 1, jcx(jcomp))
         call deblnk (text)
         write (*,1080) text 

      end if 
 
      write (*,1090)
      read (*,'(a)') y
      if (y.eq.'y'.or.y.eq.'Y') goto 10

      kds(jcomp) = ids

1010  format (2x,i2,' - ',a5)
1020  format (/,'Invalid input, try again:',/)
1030  format (/,'How many components in the ',a,' of the',
     *          ' composition (<',i2,')?')
1040  format (/,'Enter component indices and weighting factors for the '
     *        ,a,':')
1070  format (/,'The compositional variable is:')
1080  format (/,'The compositional variable is: ',a,/)
1090  format ('Change it (y/n)?')
1100  format (/,'Compositions are defined as a ratio of the form:',/,
     *        4x,' Sum {w(i)*n(i), i = 1, c1} / Sum {w(i)*n(i), i',
     *        ' = c2, c3}',/,15x,
     *        ' n(j) = ',a,' of component j',/,15x,
     *        ' w(j) = weighting factor of component j (usually 1)')
1130  format (15('+',1x,f4.1,1x,a5,1x))
1140  format ('Enter zero to use the numerator as a composition.')
1150  format (/,a,/)  

      end

      subroutine rnam1 (iex,xnam)
c----------------------------------------------------------------------
c read a solution/compound name from console, return
c iex = -id if a compound
c iex = ikp if a solution
c iex = 0 if invalid choice
c----------------------------------------------------------------------
      implicit none

      integer iex

      character*10 xnam
c----------------------------------------------------------------------
      iex = 0

      do 

         write (*,1040) 
         read (*,'(a)') xnam

         call matchj (xnam,iex)

         if (iex.ne.0) exit 
         write (*,1100) xnam

      end do 

1040  format (/,'Enter solution or compound name (left justified): ')
1100  format (/,'No such entity as ',a,', try again: ')

      end

      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the x3 array (post vertex) or zcoor array (in vertex). 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, jd, ids
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 xcoordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            x(i,j) = x3(jd,i,j) 
         end do 
      end do 

      end 
