c routines common to psect and reader and NOT called by vertex/meemum

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

      logical fileio, flsh, anneal
      integer ncol, nrow
      common/ cst226 /ncol,nrow,fileio,flsh,anneal

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      if (icopt.eq.7.and.fileio) then 
c                                 1d-fractionation with file input,
c                                 use nodal coordinates:
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

      else if (icopt.lt.9) then 

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
c                                 2d fractionation, this is not gonna 
c                                 work if fileio, oder?
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

      else if (icopt.eq.12) then 

         vnm(1) = 'n(H2O)     '
         vnm(2) = 'node#      '

         vmn(2) = 1d0
         vmx(2) = dfloat(iopt(36)) + 1d0
         var(2) = 1d0
 
         vmn(1) = 0d0
         vmx(1) = nopt(36)*dfloat(iopt(36))
         var(1) = 0d0

         v(1) = vmin(1)
         v(2) = vmin(2)

         jvar = ipot + 2

         do i = 3, jvar
            vnm(i) = vname(jv(i-2))
            vmx(i) = vmax(jv(i-2))
            vmn(i) = vmin(jv(i-2))
            var(i) = vmin(jv(i-2))
         end do

      end if 

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

      character*5 y*1, units*13, text*195, what*9, sym*1

      integer jcomp, ier, i, ids, count

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      character cname*5
      common/ csta4  /cname(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer icps, jcx, jcx1, kds
      logical stol, savg, spec
      double precision rcps, a0
      common/ comps /rcps(k7,2*k5),a0(k7,2),icps(k7,2*k5),jcx(2*k5),
     *               jcx1(2*k5),kds(2*k5),stol(i11),savg(i11),spec(2*k5)

      integer spct
      double precision ysp
      character*8 spnams
      common/ cxt34 /ysp(l10,k5),spct(h9),spnams(l10,h9)
c----------------------------------------------------------------------
c                                choose components vs species
      write (*,1000) fname(ids)
      read (*,'(a)') y

      if (y.eq.'y'.or.y.eq.'Y') then 
         spec(jcomp) = .true.
         what = ' species'
      else
         spec(jcomp) = .false.
         what = 'component'
      end if 
c                                set units for composition
      if (spec(jcomp)) then
         units = 'mole fraction'
         sym = 'y'
      else if (iopt(2).eq.0) then
         units = 'molar amount '
         sym = 'n'
      else 
         units = ' mass amount '
         sym = 'm'
      end if 
c                                get the composition to be contoured
10    if (lopt(22)) then
c                                with moronic constant 
         write (*,1100) sym, sym, units, what, what
      else 
         write (*,1110) sym, sym, units, what, what
c                                zero the constant
         a0(jcomp,1) = 0d0
         a0(jcomp,2) = a0(jcomp,1)
      end if 
  
      do 

         if (spec(jcomp)) then 
            write (*,1030) what,'numerator',k5+1
         else 
            write (*,1030) what//'s','numerator',k5+1
         end if 

         read (*,*,iostat=ier) jcx(jcomp)

         if (ier.ne.0.or.jcx(jcomp).lt.1) then
            write (*,1020)
            cycle 
         end if 

         exit

      end do 
c                                define the numerator
      do 

         write (*,1040) what,'numerator'

         if (spec(jcomp)) then
            write (*,1010) (i,spnams(i,ids), i = 1, spct(ids))
            count = spct(ids)
         else 
            write (*,1010) (i,cname(i), i = 1, icomp)
            count = icomp
         end if 

         read (*,*,iostat=ier) (icps(i,jcomp),rcps(i,jcomp), 
     *                                     i = 1, jcx(jcomp))
         do i = 1, jcx(jcomp)
            if (icps(i,jcomp).lt.1.or.icps(i,jcomp).gt.count) then
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

      if (lopt(22)) then 
         write (*,1050) 'a1'
         call rdnumb (a0(jcomp,1),0d0,i,0,.true.)
      end if 
c                                define the denominator
      do 

         if (spec(jcomp)) then 
            write (*,1030) what,'denominator',k5+1-jcx(jcomp)
         else 
            write (*,1030) what//'s','denominator',k5+1-jcx(jcomp)
         end if 

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

            write (*,1040) what,'denominator'

            if (spec(jcomp)) then
               write (*,1010) (i,spnams(i,ids), i = 1, spct(ids))
            else 
               write (*,1010) (i,cname(i), i = 1, icomp)
            end if 

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

            if (lopt(22)) then 

               write (*,1050) 'a2'
               call rdnumb (a0(jcomp,2),0d0,i,0,.true.)
c                                show the user the composition: 
               write (*,1070)   

               if (spec(jcomp)) then
                  write (text,1120) a0(jcomp,1),(rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = 1, jcx(jcomp))
               else           
                  write (text,1130) a0(jcomp,1),
     *                         (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                          i = 1, jcx(jcomp))
               end if 

            else 

               write (*,1070)  

               if (spec(jcomp)) then
                  write (text,1120) (rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = 1, jcx(jcomp))
               else          
                  write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)),
     *                                            i = 1, jcx(jcomp))
               end if 

            end if  

            call deblnk (text)
            write (*,1150) text    
            write (*,*) '   divided by '

            if (lopt(22)) then 

               if (spec(jcomp)) then
                  write (text,1130) a0(jcomp,1),(rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               else  
                  write (text,1130) a0(jcomp,2),
     *                        (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               end if 

            else

               if (spec(jcomp)) then
                  write (text,1120) (rcps(i,jcomp),
     *                               spnams(icps(i,jcomp),ids),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               else  
                  write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)),
     *                               i = jcx(jcomp)+1, jcx1(jcomp))
               end if 

            end if 

            call deblnk (text)
            write (*,1150) text 

            exit 

         end do 

      else 

         if (lopt(22)) then 

            if (spec(jcomp)) then
               write (text,1130) a0(jcomp,1),(rcps(i,jcomp),
     *                           spnams(icps(i,jcomp),ids),
     *                           i = 1, jcx(jcomp))
            else  
               write (text,1130) a0(jcomp,1),
     *                           (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                           i = 1, jcx(jcomp))
            end if 

         else

            if (spec(jcomp)) then
               write (text,1120) (rcps(i,jcomp),
     *                           spnams(icps(i,jcomp),ids),
     *                           i = 1, jcx(jcomp))
            else  
               write (text,1120) (rcps(i,jcomp),cname(icps(i,jcomp)), 
     *                                       i = 1, jcx(jcomp))
            end if 

         end if  

         call deblnk (text)
         write (*,1080) text 

      end if 
 
      write (*,1090)
      read (*,'(a)') y
      if (y.eq.'y'.or.y.eq.'Y') goto 10

      kds(jcomp) = ids

1000  format (/,'Define the composition in terms of the species/endmem',
     *          'bers of ',a,'(y/n)?',//,'Answer no to define a ',
     *          'composition in terms of the systems components.',/,
     *          'Units (mass or molar) are controlled by the ',
     *          'composition keyword in',/,'perplex_option.dat.')
1010  format (2x,i2,' - ',a)
1020  format (/,'Invalid input, try again:',/)
1030  format (/,'How many ',a,' in the ',a,' of the',
     *          ' composition (<',i2,')?')
1040  format (/,'Enter ',a,' indices and weighting factors for the '
     *        ,a,':')
1050  format (/,'Enter the optional constant ',a,' [defaults to 0]:')
1070  format (/,'The compositional variable is:')
1080  format (/,'The compositional variable is: ',a,/)
1090  format ('Change it (y/n)?')
1100  format (/,'Compositions are defined as a ratio of the form:',/,
     *        4x,'[a1 + Sum {w(i)*n(i), i = 1, c1}] / [a2 + Sum {w(i)*',
     *        a,'(i), i = c2, c3}]',/,15x,
     *        a,'(j)   = ',a,' of ',a,' j',/,15x,
     *        'w(j)   = weighting factor of ',a,' j (usually 1)',/,
     *    15x,'a1, a2 = optional constants (usually 0)')
1110  format (/,'Compositions are defined as a ratio of the form:',/,
     *        4x,' Sum {w(i)*n(i), i = 1, c1} / Sum {w(i)*',
     *        a,'(i), i = c2, c3}',/,15x,
     *        a,'(j)   = ',a,' of ',a,' j',/,15x,
     *        'w(j)   = weighting factor of ',a,' j (usually 1)')
1120  format (15('+',1x,f4.1,1x,a5,1x))
1130  format (f4.1,1x,15('+',1x,f4.1,1x,a5,1x))
1140  format ('Enter zero to use the numerator as a composition.')
1150  format (/,a,/)  

      end

      subroutine rnam1 (iex,xnam,what)
c----------------------------------------------------------------------
c read a solution name (what = 0) compound name (what = 1) or either
c (what = 2) from console, return
c iex = -id if a compound
c iex = ikp if a solution
c iex = 0 if invalid choice
c----------------------------------------------------------------------
      implicit none

      integer iex, what

      character*10 xnam
c----------------------------------------------------------------------
      iex = 0

      do 

         if (what.eq.0) then 
            write (*,1040) 'solution' 
         else if (what.eq.1) then 
            write (*,1040) 'compound' 
         else 
            write (*,1040) 'solution or compound' 
         end if  

         read (*,'(a)') xnam

         call matchj (xnam,iex)

         if (iex.ne.0) exit 
         write (*,1100) xnam

      end do 

1040  format (/,'Enter ',a,' (left justified): ')
1100  format (/,'No such entity as ',a,', try again: ')

      end

      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the x3 array (post vertex) or zcoor array (in vertex). 

c getxz is duplicated in resub.f, and consequently requires argument ID
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, jd, ids
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 xcoordinates for the final solution, a
c                                 leetle witz.
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)
c----------------------------------------------------------------------

      do i = 1, ostg(ids)
         do j = 1, ispg(ids,i)
            x(i,j) = x3(jd,i,j)
         end do
      end do

      end
