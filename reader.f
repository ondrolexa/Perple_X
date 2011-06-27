c       lun restrictions
c       n8 - echo 
c       n7 - geotherm
c       n5 - out
c       n6 - special out
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, output, fake

      integer imode, ierr, i

      character*100 n5name,n6name

      logical oned
      common/ cst82 /oned

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer iam
      common/ cst4 /iam

      integer icps, jcx, jcx1, kds
      logical stol, savg
      double precision rcps
      common/ comps /rcps(k7,2*k5),icps(k7,2*k5),jcx(2*k5),jcx1(2*k5),
     *               kds(2*k5),stol(h9),savg(h9)
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 3
c                                 version info
      call vrsion
c                                 initialize some flags
      first = .true.
      output = .false.
      fake   = .false.
      rxn = .false.
c                                 this could be eliminated by passing first 
c                                 to chsprp.
      do i = 1, h9
         stol(i) = .false.
      end do 
c                                 initialize ind only used for outmod, but
c                                 only set for 1d calcs in getind
      ind = 1
c                                 read input from unit n1 (terminal/disk).
c                                 input1 also initializes:
c                                 equilibrium counters; units n2 n4 and n6;
c                                 and the limits for numerical results.
      call input1 (first,output)
c                                 set ivar flag, this indicates the number
c                                 of possible independent plotting variables, jvar
c                                 indicates the number of thermodynamic variables
      ivar = 2
      if (icopt.eq.10) ivar = 3
c                                 don't allow users to do anything
c                                 other than gridded min
      if (icopt.lt.5) call error (4,1d0,icopt,'PSVDRAW')
c                                 read thermodynamic data on unit n2:
      call input2 (fake)
c                                 read autorefine lists
      call setau1 (output)
c                                 read data for solution phases on n9:
      call input9 (fake,output)

      call setau2 (output)
c                                 read the plot file for grid info
      call plinp
c                                 organize variables 
      call getvar
c                                 initialize the grid parameters
      call setvar 
c                                 read bulk composition data file:
      call bplinp 

      do 

         write (*,1000) 
         if (.not.oned) write (*,1010)
         write (*,1020)
         if (.not.oned) write (*,1025)
         write (*,1026)

         read (*,*,iostat=ierr) imode
         if (ierr.ne.0) cycle 

         if (first.and.imode.eq.1) then 
c                                 make console output echo to rpl file
            call fopenn (n8,0,n5name,n6name)
            first = .false.

         end if 

         if (imode.eq.1) then          
c                                 mode 1 - print assemblage data at a 
c                                 specific p-t condition
            call mode1 

         else if (imode.eq.2) then 

            if (oned) then 
               write (*,1030) 
               cycle
            end if 

            call mode2 

         else if (imode.eq.3) then 
c                                extract the data
            if (oned) then 
               call mode31 
            else 
               call mode3 
            end if 

         else if (imode.eq.4) then 

            call mode4 

         else if (imode.eq.0) then 

            exit 

         end if 

      end do 
c                                 close "echo" file
      if (imode.eq.1.or.imode.eq.3) close (n8)

1000  format ('Select operational mode:',/,
     *        4x,'1 - properties at specified conditions')
1010  format (4x,'2 - properties on a 2d grid')
1020  format (4x,'3 - properties along a 1d path')        
1025  format (4x,'4 - as in 3, but input from file')
1026  format (4x,'0 - EXIT')
1030  format (/,'Invalid choice for 1d grids',/)

      end 

      subroutine mode2 
c----------------------------------------------------------------------
c sample data on an x-y grid 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node

      integer i, j, nxy(2), dim

      double precision tmin(2), tmax(2), dx(2)

      character*100 n6name, n5name, yes*1

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      character vnm*8
      common/ cxt18a /vnm(l3)
c----------------------------------------------------------------------
      node = .false. 
      dim = 2
c                                 select the property
      call chsprp 

      do i = 1, iprop
         prmx(i) = -1d99
         prmn(i) = 1d99
      end do 
c                                 set up coordinates etc
c                                 allow restricted plot limits
      write(*,1040)
      read (*,'(a)') yes 

      if (yes.eq.'y'.or.yes.eq.'Y') then 

         do i = 1, 2
30          write (*,1060) vnm(i),vmn(i),vmx(i)
            read (*,*,err=30) tmin(i),tmax(i)
         end do 

      else 

         tmin(1) = vmn(1)
         tmin(2) = vmn(2)
         tmax(1) = vmx(1)
         tmax(2) = vmx(2)

      end if 
c                                 number of grid points
      write (*,1080) vnm(1),vnm(2)
      read (*,*) nxy
         
      do i = 1, 2
         tmin(i) = tmin(i) + (tmax(i)-tmin(i))*1d-6
         tmax(i) = tmax(i) - (tmax(i)-tmin(i))*1d-6
         dx(i) = (tmax(i)-tmin(i))/dfloat(nxy(i)-1)
      end do 

      call tabhed (tmin,dx,nxy,dim,n5name,n6name)

      do j = 1, nxy(2)

         var(2) = tmin(2) + dx(2)*dfloat(j-1)

         do i = 1, nxy(1) 

            var(1) = tmin(1) + dx(1)*dfloat(i-1)

            call polprp (dim)
 
         end do
 
      end do 
c                                 wrap up the calculation
      call finprp (dim,n5name,n6name,node) 

1040  format (/,'Change default variable range (y/n)?')
1060  format (/,'Current limits on ',a,' are: ',g14.7,'->',g14.7,/,
     *          'Enter new values:')
1080  format (/,'Enter number of nodes in the ',a,' and ',a,
     *          ' directions:')

      end 

      subroutine mode1 
c----------------------------------------------------------------------
c read x-y coordinates for a 2-d section from the console. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit, nodata

      integer itri(4),jtri(4),ijpt

      double precision wt(3)
c----------------------------------------------------------------------
      do 

         call readxy (quit)

         if (quit) exit

         call triang (itri,jtri,ijpt,wt)

         if (ijpt.eq.0) then 
            nodata = .true.
         else 
            call getloc (itri,jtri,ijpt,wt,nodata)
         end if 

         if (nodata) then 
            write (*,1000) 
         else 
            call calpr0 (6)
            call calpr0 (n8)
         end if

      end do 

1000  format (/,'No data at this condition, presumably because',
     *          ' minimization failed.',/)

      end 

      subroutine amiin1 (j,left)
c----------------------------------------------------------------------
c amiin1 - identifies the grid point associated with coordinate x and
c          indicates if x is left of the nodal coordinate.        
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer j

      logical left 

      double precision res

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar
c----------------------------------------------------------------------
c                                 find node associated with condition
      res = var(1)-vmn(1)
      j = int(res/dvr(1))
      res = res - j*dvr(1)

      if (res.lt.0d0) then
         left = .true.
      else 
         left = .false.
      end if 

      if (res.gt.0.5d0*dvr(1)) then
         j = j + 1
         left = .true.
      end if 

      j = j + 1
       
      end 

      subroutine amiin2 (i,j)
c----------------------------------------------------------------------
c amiin - identifies the grid point associated with coordinate x-y 
c         and the sector of the x-y coordinate relative the nodal point        
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j

      logical left, down

      double precision res

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar
c----------------------------------------------------------------------
c                                 find node associated with condition
      res = (var(1)-vmn(1))/dvr(1) + 1d0
      i = int(res) 

      if (res-dfloat(i).gt.0.5d0) then 
         i = i + 1
         left = .true.
      else 
         left = .false.
      end if 

      res = (var(2)-vmn(2))/dvr(2) + 1d0
      j = int(res)

      if (res-dfloat(j).gt.0.5d0) then 
         j = j + 1
         down = .true.
      else 
         down = .false.
      end if 

      end 

      subroutine setval 
c--------------------------------------------------------------------
c setval sets the values of the thermodynamic variables after a call
c to readxy gets the section coordinates in array var, five cases:

c icopt = 10 -> var(1) is the index of coordinate array
c icopt = 9  -> frac2d calculations x-y indirectly related to p-t
c icont = 1  -> independent variables are the 1st and 2nd potentials
c icont = 2  -> 1st independent variable is a composition variable,  
c              2nd independent variable is the 1st potential (iv1)
c icont = 3  -> independent variables are compositional variables
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind, j

      double precision wt 

      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps

      integer ipot,jv,iv1,iv2,iv3,iv4,iv5
      common / cst24 /ipot,jv(l2),iv1,iv2,iv3,iv4,iv5

      integer icont
      double precision dblk,cx
      common/ cst314 /dblk(3,k5),cx(2),icont

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer iind, idep
      double precision c0,c1,c2,c3,c4,c5
      common/ cst316 /c0,c1,c2,c3,c4,c5,iind,idep
 
      double precision vip
      common/ cst28 /vip(l2,k2)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c----------------------------------------------------------------------
      if (icopt.eq.10) then

          ind = idint(var(1))
          wt = var(1) - dfloat(ind)

          do j = 1, ipot
             v(jv(j)) = vip(j,ind)*(1d0-wt)
             if (wt.gt.0d0) v(jv(j)) = v(jv(j)) + vip(j,ind+1)*wt
          end do 

          var(2) = v(jv(1))
          var(3) = v(jv(2))

      else if (icopt.eq.9) then 
c                                 change sign on dz because of downward
c                                 directed depth coordinate.
         call fr2dpt (var(1),-var(2))
         var(3) = v(1)
         var(4) = v(2)

      else if (icont.eq.1) then 

         v(iv1) = var(1)
         v(iv2) = var(2)
         call incdp0
         if (idep.ne.0) var(jvar) = v(idep) 

      else if (icont.eq.2) then 

         cx(1) =  var(1)
         call setblk 

         v(iv1) = var(2)
         call incdep (iv1)
         if (idep.ne.0) var(jvar) = v(idep) 

      else 

         cx(1) = var(1)
         cx(2) = var(2)
         call setblk

      end if 

      end

      subroutine readxy (quit)
c----------------------------------------------------------------------
c read x-y coordinates for a 2-d section from the users console, these
c are then assigned to the thermodynamic variables by calling setval.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, ier
      logical quit

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3) 

      logical oned
      common/ cst82 /oned
c----------------------------------------------------------------------
      if (.not.oned) then

         do 

            quit = .false.

            write (*,1000) vnm(1),vnm(2)
            read (*,*,iostat=ier) var(1),var(2)
            if (ier.ne.0) cycle 

            if (var(1)+var(2).eq.198d0) quit = .true.

            if (.not.quit) then 
                
               quit = .false.

               do i = 1, 2
                  if (vmn(i).lt.vmx(i)) then
                     if (var(i).lt.vmn(i).or.var(i).gt.vmx(i)) then  
                        write (*,1010) vnm(i),vmn(i),vmx(i)
                        quit = .true.
                     end if 
                  else 
                     if (var(i).lt.vmx(i).or.var(i).gt.vmn(i)) then  
                        write (*,1010) vnm(i),vmn(i),vmx(i)
                        quit = .true.  
                     end if 
                  end if 
               end do 

               if (quit) cycle

            end if 

            exit
 
         end do 

      else 

         do 

            quit = .false.

            write (*,1020) vnm(1)
            read (*,*) var(1)

            if (var(1).eq.999d0) quit = .true.

            if (.not.quit) then 
               if (vmn(1).lt.vmx(1)) then
                  if (var(1).lt.vmn(1).or.var(1).gt.vmx(1)) then  
                     write (*,1010) vnm(1),vmn(1),vmx(1)
                     cycle  
                  end if 
               else 
                  if (var(1).lt.vmx(1).or.var(1).gt.vmn(1)) then  
                     write (*,1010) vnm(1),vmn(1),vmx(1)
                     cycle  
                  end if 
               end if 
            end if
            
            exit
 
         end do 

      end if 

      if (.not.quit) call setval

1000  format (/,'Enter ',a,' and ',a,' (99 and 99 to quit):')
1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1020  format (/,'Enter ',a,' (999 to quit):')

      end 

      subroutine triang (itri,jtri,ijpt,wt)
c----------------------------------------------------------------------
c routine to extract interpolation points for an x-y point (set by readxy). 
c the algorithm seeks points over an interval of jinc(1), jnterp 
c should be read from perplex_option.dat but is currently set here to 1. 
c This whole process could be done 
c much more efficiently (and better) if the grid points were mapped 
c to a triangular mesh. 

c    returns:
c              ijpt       - number of good interpolation points
c              itri, jtri - nodal cordinates of the (3) interpolation 
c                           points
c              wt         - weights of interpolation points
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x, y, wt(3), px(3), py(3), div, dst(4), x0, x1

      integer jloc, iloc, j, i, k, l, jmax, np, jd, j0,
     *        itri(4), jtri(4), ijpt, iam, jam, kinc, jmin, 
     *        imin, imax, ktri(4), ltri(4), ibest

      logical rinsid, in, isok, jn, warned, left

      integer pi(4,4)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer jlow,jlev,loopx,loopy,jinc
      common/ cst312 /jlow,jlev,loopx,loopy,jinc

      integer igrd
      common/ cst311/igrd(l7,l7)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                 global assemblage data
      integer icog,jcog
      common/ cxt17 /icog(k2),jcog(k2)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      logical oned
      common/ cst82 /oned

      save warned, pi
      data warned/.false./
      data pi/1, 2, 3, 4, 1, 2, 4, 3, 1, 3, 4, 2, 2, 3, 4, 1/
c----------------------------------------------------------------------

      x = var(1)
      y = var(2)
c                                 get nodal coordinate   
      call xy2ij (iloc,jloc,left)
c                                 identify the assemblage    
      jd = igrd(iloc,jloc)
c                                 the iap(jd) check works if k2 is inconsistent
c                                 between vertex and werami. the question is then
c                                 whether there's any point in checking jd?           
      if (jd.eq.k2.or.iap(jd).eq.0) then
c     if (jd.eq.k2) then 
c                                 no data at point, set ijpt = 0 and return
         ijpt = 0 
         goto 99 
      end if 

      ias = iap(jd)
      np = iavar(1,ias)
      ijpt = 1
      itri(1) = iloc
      jtri(1) = jloc 
c                                 exit if interpolation is off
      if (iopt(4).eq.0) goto 99   
c                                 check for solvus if nopt(8) < 1, 
c                                 i, e., solvus testing is on
      if (nopt(8).lt.1d0) then  

         do i = 1, np-1
            do j = i+1, np
               if (idasls(i,ias).eq.idasls(j,ias)) then
c                                 a solvus, turn interpolation off, warn and return
                  if (.not.warned) then 
                     warned = .true.
                     write (*,1000) 
                  end if

                  goto 99

               end if 
            end do 
         end do 
      end if 

      ijpt = 0

      iam = iap(igrd(iloc,jloc))
c                                 interpolation for 1d grids                             
      if (oned) then 
c                                 set jloc to the real node
         j0 = jloc
         jloc = jcog(jd)
c                                 this is gonna move j even if no interpolation 
c                                 is possible, but we don't care, right?
         jtri(1) = jloc
c                                 find a real point with the same assemblage to 
c                                 the left
         jmin = 0

         do j = jloc - 1, jloc - iopt(4)*jinc, -1   
    
            if (j.lt.1) exit 

            if (iap(igrd(1,j)).eq.ias) then 
c                                 is the point real?
               if (jcog(igrd(1,j)).ne.j) cycle

               jmin = j
c                                  found the assemblage
               exit

            else 
c                                  ran into a new assemblage
               exit

            end if
         end do             
c                                  find a real point to the right
         jmax = 0 

         do j = jloc + 1, jloc + iopt(4)*jinc

            if (j.gt.loopy) exit

            if (iap(igrd(1,j)).eq.ias) then 
c                                 is the point real?
               if (jcog(igrd(1,j)).ne.j) cycle

               jmax = j
c                                 found the assemblage
               exit

            else 
c                                 ran into a new assemblage
               exit

            end if

         end do   

         ijpt = 2
         itri(ijpt) = 1
c                                 check cases, here only interpolation is 
c                                 allowed. 
         if (j0.eq.jloc) then

            if (left.and.jmin.gt.0) then
               jtri(ijpt) = jmin
            else if ((.not.left).and.jmax.gt.0) then 
               jtri(ijpt) = jmax
            else 
               ijpt = 1
            end if 

         else if (j0.gt.jloc) then 

            if (jmax.ne.0) then 
               jtri(ijpt) = jmax
            else 
               ijpt = 1
            end if 

         else 

            if (jmin.ne.0) then 
               jtri(ijpt) = jmin
            else 
               ijpt = 1
            end if

         end if 

         if (ijpt.eq.1) goto 99
c                                 compute weights of the interpolation points   
         x1 = vmn(1) + dfloat(jtri(2)-1)*dvr(1) 
         x0 = vmn(1) + dfloat(jtri(1)-1)*dvr(1)    
   
         wt(1) = (x1-x)/(x1-x0)
         wt(2) = 1d0 - wt(1)

         goto 99 

      end if 
c                                 make a spiral-like search outward
c                                 from the node to find interpolation
c                                 points
      kinc = 0       
c                                 the multiplier on jinc (the increment
c                                 for the lowest level grid) controls
c                                 the search area. 
      do while (kinc.le.iopt(4)*jinc)

         jmin = jloc - kinc
         jmax = jloc + kinc 

         do j = jmin, jmax
c                                 skip out of bounds points           
            if (j.lt.1.or.j.gt.loopy) cycle

            imin = iloc - kinc
            imax = iloc + kinc

            do i = imin, imax
c                                 skip out of bounds points           
               if (i.lt.1.or.i.gt.loopx) cycle
c                                 skip interior points (this is sloppy)
               if (j.ne.jmin.and.j.ne.jmax.and.
     *             i.ne.imin.and.i.ne.imax) cycle 
c                                 is the point the same assemblage?
               if (iap(igrd(i,j)).ne.iam) cycle
c                                 pointer to reference node
               jam = igrd(i,j)
c                                 if so, is the point real?
               if (icog(jam).ne.i.or.jcog(jam).ne.j) cycle
c                                 if here the point is valid, now
c                                 check if it's geometrically feasible
               if (ijpt.lt.2) then 
c                                 always take the two points
                  ijpt = ijpt + 1
                  itri(ijpt) = i
                  jtri(ijpt) = j

               else if (ijpt.eq.2) then 

                  itri(3) = i
                  jtri(3) = j      
c                                 check if on line of previous
c                                 points
                  if (isok(itri,jtri)) then 
                     ijpt = 3
                     in = rinsid(itri,x,jtri,y,dst(1))
                  end if 

               else if (ijpt.eq.3) then 
c                                 four permutations are possible
c                                 check all 
                  itri(4) = i
                  jtri(4) = j 
                  ibest = 1
     
                  do k = 2, 4
                     do l = 1, 3
                        ktri(l) = itri(pi(l,k))
                        ltri(l) = jtri(pi(l,k))
                     end do 


                     if (isok(ktri,ltri)) then

                        jn = rinsid(ktri,x,ltri,y,dst(k))

                        if (jn.and..not.in) then
c                                 no in prior bounding triangle, so
c                                 k is now the best guess  
                           ibest = k 
                           in = jn                    
                        else if (jn.and.in.or..not.jn.and..not.in)
     *                                                             then 
c                                 both are bounding, compare weights
c                                 (total distance to vertices).
                           if (dst(ibest).gt.dst(k)) then
                              ibest = k 
                           end if                              
                        end if 
                     end if 
                  end do    

                  if (ibest.ne.1) then 
c                                 load the best choice
                     do k = 1, 3
                        itri(k) = itri(pi(k,ibest))
                        jtri(k) = jtri(pi(k,ibest))
                     end do 
                  end if 

               end if 
            end do 
         end do 

         kinc = kinc + 1

      end do 
  
      if (ijpt.eq.3) then 

         do j = 1, 3
            px(j) = vmn(1) + (itri(j)-1)*dvr(1)
            py(j) = vmn(2) + (jtri(j)-1)*dvr(2)
         end do
c                                 compute iterpolation coefficients
         div = px(2)*py(3)-px(1)*py(3)-px(2)*py(1)+
     *         px(3)*py(1)-px(3)*py(2)+px(1)*py(2)
c                                 z[1] coef
         wt(1) = (y*(px(3)-px(2))+px(2)*py(3)-
     *            px(3)*py(2)+x*(py(2)-py(3)))/div
c                                 z[2] coef
         wt(2) = (px(3)*py(1)-px(1)*py(3)+
     *            y*(px(1)-px(3))-x*(py(1)-py(3)))/div
c                                 z[3] coef
         wt(3) = (x*(py(1)-py(2))+px(1)*py(2)-px(2)*py(1)
     *           +y*(px(2)-px(1)))/div
c                                 if the triangle is non-bounding
c                                 and extrapolation is off reset 
c                                 counter.
         if (.not.in.and.iopt(5).eq.0) ijpt = 1

      else 

         ijpt = 1 

      end if 

99    if (ijpt.eq.1) wt(1) = 1d0

1000  format (/,'**warning ver637** Immiscibility occurs in one or ',
     *          'more phases ',/,'interpolation will be turned off ',
     *          'at all affected nodes.',/,'To overide this feature ',
     *          'at the risk of computing inconsistent properties',/,
     *          'set solvus_tolerance = 1 and rerun VERTEX',/)

      end 

      logical function rinsid (itri,xp,jtri,yp,dst)
c----------------------------------------------------------------------
c function to determine if triangle i1-j1 ... i3-j3 bounds point xp - yp
c     returns true if internal or on an edge
c     returns false if outside or on a vertex.
c algorithm after jimscott@blackpawn.com
c                                 JADC 1/2005
c----------------------------------------------------------------------

      implicit none

      include 'perplex_parameters.h'

      logical rsmsid

      double precision xp, yp, x(3), y(3), dst, dist

      integer itri(4), jtri(4), j

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar
c----------------------------------------------------------------------
c                                convert nodal coordinates to real
c                                cordinates
      dst = 0d0

      do j = 1, 3
         x(j) = vmn(1) + (itri(j)-1)*dvr(1)
         y(j) = vmn(2) + (jtri(j)-1)*dvr(2)
         dst = dst + dist(x(j),y(j),itri(j),jtri(j))
      end do 
c                                1->2 join
      if (rsmsid(x(2)-x(1),y(2)-y(1),x(3)-x(1),
     *           y(3)-y(1),xp-x(1),yp-y(1)).and.
c                                1->3 join
     *    rsmsid(x(3)-x(1),y(3)-y(1),x(2)-x(1),
     *           y(2)-y(1),xp-x(1),yp-y(1)).and.
c                                2->3 join
     *    rsmsid(x(3)-x(2),y(3)-y(2),x(1)-x(2),
     *           y(1)-y(2),xp-x(2),yp-y(2))) then

          rinsid = .true.

      else 

          rinsid = .false.

      end if 

      end 

      logical function rsmsid (x1,y1,x2,y2,x3,y3)
c----------------------------------------------------------------------
c function to determine if points (x1,y1) (x2,y2) are on the same side
c (or colinear with) the line (0,0)-(x3,y3) (if so true).

      double precision x1,y1,x2,y2,x3,y3
c                                 test the cross product, if 0 a point
c                                 is on the line.
      if ((x1*y3 - y1*x3)*(x1*y2 - y1*x2).ge.0) then
         rsmsid = .true.
      else
         rsmsid = .false.
      end if 

      end 

      double precision function dist (x,y,i,j)
c----------------------------------------------------------------------
c get distance from nodal point i,j to coordinate x-y      

      implicit none
 
      include 'perplex_parameters.h'

      integer i, j

      double precision dely, delx, x, y

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar
c-----------------------------------------------------------------------
c                                 find normalized distance 
      delx = (x - vmn(1))/dvr(1) - (i - 1)
      dely = (y - vmn(2))/dvr(2) - (j - 1)
      dist = sqrt (delx**2 + dely**2)

      end   

      subroutine polprp (dim)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical nodata

      integer itri(4), jtri(4), ijpt, lop, icx, komp, i, dim

      double precision wt(3)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10) 

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer igrd
      common/ cst311/igrd(l7,l7)
c----------------------------------------------------------------------
c                                 set variables to x-y value
      call setval

      do i = 1, iprop

         lop = kop(i)
         icx = kcx(i)
         komp = k2c(i)
         lflu = kfl(i)
c                                 initialize, necessary?
         prop(i) = nopt(7)
c                                 get node(s) to extract value
         call triang (itri,jtri,ijpt,wt)

         if (ijpt.eq.0) then 
c                                 missing data at the node
            call badnum (dim)

            return 

         else 
c                                 compute all properties
            call getloc (itri,jtri,ijpt,wt,nodata)

            if (nodata) then 

               call badnum (dim)

               return  

            else 
c                                 get the properties of interest
               if (lop.eq.25) then 
c                                 get all modes
                  call allmod

                  exit 

               else if (lop.eq.36.or.lop.eq.38) then 
c                                 multiple property lists, PHEMGP
                  call allprp (dim)

                  exit 

               else if (lop.ne.24) then 
c                                 general properties
                  call getprp (prop(i),lop,icx,komp,.false.)

               else 
c                                 assemblage index request, lots
c                                 of redundant calc, but should be 
c                                 moved to getptp.
c                                 no need to call triang/getlow
                  call xy2ij (itri(1),jtri(1),nodata)

                  prop(i) = iap(igrd(itri(1),jtri(1)))

               end if 

            end if 

         end if 

      end do 

      if (lop.ne.36.and.lop.ne.38) call outprp (dim)

      end 

      subroutine badnum (dim)  
c----------------------------------------------------------------
c badnum: writes missing data message; assigns and outputs badnumber
c data.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i 

      character vnm*8
      common/ cxt18a /vnm(l3)   

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot
c----------------------------------------------------------------------
      write (*,1000) vnm(1),var(1),vnm(2),var(2),nopt(7)
c                                 for phemgp format:
      ntot = 0
      tname = 'Missing data  '
c                                 in general
      do i = 1, iprop
         prop(i) = nopt(7)
      end do 

      call outprp (dim)

1000  format ('Missing data at ',a,'=',g12.5,', ',a,'=',g12.5,
     *        ' assigned ',g12.5,' to all properties')

      end 

      subroutine outprp (dim)
c----------------------------------------------------------------------
c outprp outputs properties computed by chsprp for dim-dimensional tables: 
c   dim    = 1 tables include independent variables;
c   dim    = 2, kcx(1) ne 999, tables include independent variables if lopt(15) = T;
c   kcx(1) = 999, write phemgp format that includes a name, counter, and 
c            independent variables. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot
c----------------------------------------------------------------------
      if (kcx(1).eq.999) then 
c                                 write phemgp format
         write (n5,'(a14,1x,7x,i2,6x,200(g14.7,1x))') tname,ntot,
     *                                               (var(i),i=1,ivar), 
     *                                               (prop(i),i=1,iprop)

      else if (lopt(15).or.dim.eq.1) then 
c                                 write spreadsheet tab format
         write (n5,'(200(g14.7,1x))') (var(i),i=1,ivar), 
     *                                (prop(i),i=1,iprop)
      else 
c                                 write compact tab format
         write (n5,'(200(g14.7,1x))') (prop(i),i=1,iprop)
      end if 
c                                 check property ranges       
      do i = 1, iprop
c                                
         if (isnan(prop(i))) cycle
c                                 first check eliminates logical 
c                                 comparisons with NaN's that compaq
c                                 fortran doesn't like.
         if (.not.isnan(nopt(7))) then 
            if (prop(i).eq.nopt(7)) cycle 
         end if 
         
         if (prop(i).gt.prmx(i)) prmx(i) = prop(i)
         if (prop(i).lt.prmn(i)) prmn(i) = prop(i)
 
      end do

      end 

      subroutine getprp (prop,lop,icx,komp,aprp)
c----------------------------------------------------------------
c getprp gets properties:

c   aprp - true if called by allprp, false otherwise

c   jd   - is the pointer to the assemblage
c   jflu - if 1, bulk properties include fluid phase; else 0

c   lop  - flag indicating the property chosen
c   icx  - if lop = 6, the component chosen
c   icx  - if the identity of the phase chosen
c          icx = 0 if bulk properties requested.

c   icps - if lop = 8, the indices of the components. 
c   
c requestable properties (indicated by lop)

c 1                 Specific enthalpy (J/m3)',
c 2                 Density (kg/m3)',
c 3                 Specific Heat capacity (J/K/m3)',
c 4                 Expansivity (1/K, for volume)',
c 5                 Compressibility (1/bar, for volume)
c 6                 Weight percent of a component
c 7                 Mode (Vol %) of a compound or solution',
c 8                 Composition of a solution'
c 9                 Grueneisen thermal ratio',
c 10                Adiabatic bulk modulus (bar)',
c 11                Sound velocity (km/s)
c 12                Shear modulus (bar)',
c 13                P-wave velocity (km/s)',
c 14                S-wave velocity (km/s)',
c 15                Vp/Vs
c 16                Specific Entropy (J/K/m3)'
c 17                Entropy (J/K/kg)'
c 18                Enthalpy (J/kg)'
c 19                Heat Capacity (J/K/kg)'
c 20                Specific mass of phase (kg/m3)
c 21                Poisson's Ratio 
c 22                Molar Volume (J/bar) 
c 23                Chemical potentials (J/mol)
c 24                not passed to getprp
c 25                output all modes
c 26                Sound velocity temperature derivative (km/s/K)
c 27                P-wave velocity temperature derivative (km/s/K)
c 28                S-wave velocity temperature derivative (km/s/K)
c 29                Adiabatic bulk modulus temperature derivative (bar/K)
c 30                Shear modulus temperature derivative (bar/K)
c 31                Sound velocity pressure derivative (km/s/bar)
c 32                P-wave velocity pressure derivative (km/s/bar)
c 33                S-wave velocity pressure derivative (km/s/bar)
c 34                Adiabatic bulk modulus pressure derivative (unitless)
c 35                Shear modulus pressure derivative (unitless)
c 36                All properties of a phase or the system
c 37                Amount of a phase per unit system

c system properties from seismo:

c 1 vsys       - sys molar volume 
c 2 esys       - sys molar enthalpy 
c 3 grun       - sys gruenesien ratio
c 4 bulkm      - sys adiabatic bulk modulus
c 5 mu         - sys shear modulus
c 6 vel0       - sys sound velocity 
c 7 velp       - sys p-wave velocity
c 8 vels       - sys s-wave velocity
c 9 vp/vs      
c 10 rsys      - sys density
c 11 csys      - sys specific heat capacity
c 12 cp
c 13 alpha
c 14 beta 
c 15 entropy
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer lop, icx, id, komp

      logical aprp

      double precision prop, r, gtcomp, mode(3)

      double precision atwt
      common/ cst45 /atwt(k0)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      double precision mu
      common/ cst330 /mu(k8)
c----------------------------------------------------------------------
      if (lop.eq.6) then 
c                                 wt % of component 
         if (aflu.and.lflu.or.(.not.aflu)) then 
c                                 include fluid
            prop = fbulk(icx)*atwt(icx)/psys(17)*1d2
         else 
c                                 exclude fluid
            prop = fbulk1(icx)*atwt(icx)/psys1(17)*1d2
         end if 
      
      else if (lop.eq.23) then 
c                                 chemical potential
         prop = mu(icx)

      else 

         if (icx.eq.0) then 
c                                 a system property is  requested:
c                                 if psys1 = 0, the system is just fluid.
            if (aflu.and.lflu.or.(.not.aflu).or.psys1(1).eq.0d0) then 
c                                 if lflu the property is to include 
c                                 fluid (if present), i.e, use psys/ptot array:
               if (lop.eq.1) then 
c 1                               Specific enthalpy (J/m3)
c                                 psys(2)  is J/mole psys(1) is J/bar/mole (volume)
                  prop = psys(2)/psys(1)*1d5
               else if (lop.eq.2) then 
c 2                               Density (kg/m3)
                  prop = psys(10)
               else if (lop.eq.3) then 
c 3                               Specific Heat capacity (J/K/m3)
                 prop = psys(12)/psys(1)*1d5
               else if (lop.eq.4) then
c 4                               Expansivity (1/K, for volume)
                  prop =  psys(13)
               else if (lop.eq.5) then 
c 5                               Compressibility (1/bar, for volume)
                  prop =  psys(14)
               else if (lop.eq.7) then 
c                                 Mode %
                  prop = 1d2
               else if (lop.ge.9.and.lop.le.15) then 
c                                 gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = psys(lop-6) 
               else if (lop.eq.16) then 
c                                 specific s (j/k/m3)
                  prop = psys(15)/psys(1)*1d5
               else if (lop.eq.17) then 
c                                 S (J/K/kg)
                  prop = psys(15)/psys(1)*1d5/psys(10)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = psys(2)/psys(1)*1d5/psys(10)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = psys(12)/psys(1)*1d5/psys(10)
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (psys(8).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (psys(7)/psys(8))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = psys(1)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = psys(lop-11)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = psys(lop-14)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = psys(lop-4)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = psys(lop-6)
               end if 
            else 
c                                 fluid absent system properties:
               if (lop.eq.1) then 
c 1                               Specific enthalpy (J/m3)
                  prop = psys1(2)/psys1(1)*1d5
               else if (lop.eq.2) then 
c 2                               Density (kg/m3)
                  prop = psys1(10)
               else if (lop.eq.3) then 
c 3                               Specific Heat capacity (J/K/m3)
                  prop = psys1(12)/psys1(1)*1d5
               else if (lop.eq.4) then
c 4                               Expansivity (1/K, for volume)
                  prop =  psys1(13)
               else if (lop.eq.5) then 
c 5                               Compressibility (1/bar, for volume)
                  prop =  psys1(14)
               else if (lop.eq.7) then 
c                                 Mode %
                  prop = 1d2
               else if (lop.ge.9.and.lop.le.15) then 
c                                 m  gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = psys1(lop-6) 
               else if (lop.eq.16) then
c                                 m  specific s (j/k/m3)
                  prop = psys1(15)/psys1(1)*1d5
               else if (lop.eq.17) then
c                                 S (J/K/kg)
                  prop = psys1(15)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = psys1(2)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = psys1(12)/psys1(1)*1d5/psys1(10)
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (psys1(8).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (psys1(7)/psys1(8))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = psys1(1)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = psys1(lop-11)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = psys1(lop-14)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = psys1(lop-4)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = psys1(lop-6)
               end if 
            end if 

         else 

            if (aprp) then
c                                 call from allprp
               id = icx

            else 
c                                 normal call by phase
c                                 find the phase index
               call soltst (id,icx)

            end if 
                                 
            if (id.eq.0.and.lop.ne.7) then 
c                                 the phase is absent, set 
c                                 the property to bad_number
               prop = nopt(7)

            else if (id.ne.0) then 

               if (lop.eq.1) then
c                                 specific enthalpy 
                   prop = props(2,id)/props(1,id)*1d5
               else if (lop.eq.2) then 
c                                 density (kg/m3)
                   prop = props(10,id)
               else if (lop.eq.3) then 
c                                 specific cp 
                   prop = props(12,id)/props(1,id)*1d5
               else if (lop.eq.4) then 
c                                 expansivity
                   prop = props(13,id)
               else if (lop.eq.5) then 
c                                 compressibility
                   prop = props(14,id)
               else if (lop.eq.7) then                           
c                                 mode (%)
                   call gtmode (mode,id)
                   prop = mode(iopt(3)+1)
               else if (lop.eq.8) then 
c                                 composition (external function)
                  prop = gtcomp(id,icx,komp)
               else if (lop.ge.9.and.lop.le.15) then 
c                                 gruneisen T, K, mu, Vphi, vp, vs, vp/vs
                  prop = props(lop-6,id) 
               else if (lop.eq.16) then
c                                 specific s (j/k/m3)
                  prop = props(15,id)/props(1,id)*1d5
               else if (lop.eq.17) then 
c                                 S (J/K/kg)
                  prop = props(15,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.18) then 
c                                 H (J/kg)
                  prop = props(2,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.19) then
c                                 Cp (J/K/kg)
                  prop = props(12,id)/props(1,id)*1d5/props(10,id)
               else if (lop.eq.20) then                           
c                                 specific weight of a phase is the mass per
c                                 m3 of solid+melt

c                                 the number of moles of system/m3 is 1d5/psys(1)
c                                 twt is g/mol phase                                        
                  if (aflu.and.lflu.or.(.not.aflu)) then
c                                 total mode:
                     prop = props(16,id)*props(17,id)*1d2/psys(1)
                  else 
c                                 solid only mode:
                     prop = props(16,id)*props(17,id)*1d2/psys1(1)
                  end if 
               else if (lop.eq.21) then 
c                                 Poisson's ratio
                  if (props(8,id).eq.0d0) then 
                     prop = 0.5d0
                  else 
                     r = (props(7,id)/props(8,id))**2
                     prop = 0.5d0*(r-2d0)/(r-1d0)
                  end if 
               else if (lop.eq.22) then
c                                 molar volume
                  prop = props(1,id)
               else if (lop.eq.29.or.lop.eq.30) then
c                                 (Ks or Mu)_T
                  prop = props(lop-11,id)
               else if (lop.eq.34.or.lop.eq.35) then
c                                 (Ks or Mu)_P
                  prop = props(lop-14,id)
               else if (lop.ge.26.and.lop.le.28) then
c                                 (Vphi, Vp, Vs)_T
                  prop = props(lop-4,id)
               else if (lop.ge.31.and.lop.le.33) then
c                                 (Vphi, Vp, Vs)_P
                  prop = props(lop-6,id)

               else if (lop.eq.37) then                           
c                                 absolute amount of phase
c                                 per unit system
                  if (iopt(3).eq.0) then 
c                                 volume J/bar multiply by 1d-5 to get (m3)
                     prop = props(1,id)*props(16,id)
                  else if (iopt(3).eq.1) then   
c                                 mass (kg) 
                     prop = props(16,id)*props(17,id)/1d3
                  else if (iopt(3).eq.2) then 
c                                 mol 
                     prop = props(16,id)
                  end if 

               end if 

            end if

         end if 

      end if                
 
      end 

      logical function isok (itri,jtri)
c----------------------------------------------------------------------
c check if vertices itri-jtri define a triangle

      implicit none 

      integer itri(4), jtri(4)
      double precision m, b, di
c----------------------------------------------------------------------
      if (itri(1).eq.itri(2).and.itri(1).eq.itri(3).or.
     *    jtri(1).eq.jtri(2).and.jtri(1).eq.jtri(3)) then 

         isok = .false.

      else if (itri(1).ne.itri(2)) then 
c                          if the points are not on the same
c                          row, then they may all be on a 
c                          diagonal
         di = itri(1) - itri(2)
         m = (-jtri(2)+jtri(1))/di
         b = -(jtri(1)*itri(2)-itri(1)*jtri(2))/di + 1d-3
         if (jtri(3).eq.int(m*itri(3)+b)) then
            isok = .false.
         else
            isok = .true.
         end if 

      else

         isok = .true.         

      end if 

      end 

      double precision function gtcomp (id,icx,jcomp)
c-------------------------------------------------------------------
c function comp returns icomp'th composition of phase id
c in an assemblage whose properties have been defined in routine
c seismo, the composition is defined in routine mkcomp. 

c returns -1d99 if composition jcomp is not for phase id.
c ------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision totden, comp

      integer jcomp, j, id, icx
 
      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      integer icps, jcx, jcx1, kds
      logical stol, savg
      double precision rcps
      common/ comps /rcps(k7,2*k5),icps(k7,2*k5),jcx(2*k5),jcx1(2*k5),
     *               kds(2*k5),stol(h9),savg(h9)
c----------------------------------------------------------------------
      if (icx.eq.kds(jcomp)) then 

         comp = 0d0 
         totden = 0d0
c                                 now compute the composition:
c                                 numerator:
         do j = 1, jcx(jcomp)
            comp = comp + rcps(j,jcomp)*pcomp(icps(j,jcomp),id)
         end do
c                                 denominator:
         do j = jcx(jcomp)+1, jcx1(jcomp)
            totden = totden + rcps(j,jcomp)*pcomp(icps(j,jcomp),id)
         end do     
c                                 numerator/denominator:        
         if (totden.ne.0d0) comp = comp / totden

      else 

         comp = -1d99
 
      end if 

      gtcomp = comp

      end

      subroutine soltst (index,icx)  
c-------------------------------------------------------------------
c soltst checks for solvi, if solution icx has a solvus in 
c assemblage jd, soltsts asks the user to define a compositional
c criterion to decide which phase is relevant for property 
c contours. the index of the chosen phase is returned as index
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h' 

      character cprop*6

      integer jdsol(k5), choice, index, kdsol(k5), isol, i, j, icx, 
     *        jsol, ier, phase

      double precision cmin(k5) ,cmax(k5), tcomp, gtcomp

      character fname*10
      common/ csta7 /fname(h9)

      character*5 cname
      common/ csta4 /cname(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer icps, jcx, jcx1, kds
      logical stol, savg
      double precision rcps
      common/ comps /rcps(k7,2*k5),icps(k7,2*k5),jcx(2*k5),jcx1(2*k5),
     *               kds(2*k5),stol(h9),savg(h9)

      save cprop, cmin, cmax
c----------------------------------------------------------------------
      index = 0 
      isol = 0 
c                                 how many times does the phase occur?
      do i = 1, iavar(3,ias)
         if (icx.eq.idasls(i,ias)) then
             isol = isol + 1
             jdsol(isol) = i
         end if 
      end do
 
      if (isol.ne.0) index = jdsol(1)
c                                 the phase doesn't occur or occurs once
      if (isol.lt.2) goto 99 

      phase = idasls(index,ias)
c                                 if here, the phase must be a solution

10    if (.not.stol(phase)) then
c                                 immisicible phases are present (isol>1)
c                                 but there is no criterion (stol = .false.)
         stol(phase) = .true.

         if (iopt(2).eq.0) then 
            cprop = 'molar '
         else
            cprop = 'weight'
         end if

         write (*,1000) isol,fname(phase),cprop
         write (*,1040) (cname(i), i = 1, icomp)
         do i = 1, isol 
            write (*,1050) (pcomp(j,jdsol(i)), j = 1, icomp)
         end do 

         write (*,1030) 
c                                 get choice
         call rdnumb (tcomp,0d0,choice,1,.false.)

         if (choice.eq.2) then
c                                 average the compositions, turn 
c                                 solvus testing off
            nopt(8) = 1d0

         else
c                                 savg is a flag used only by 
c                                 werami to indicate whether 
c                                 solutions should be averaged 
c                                 within an exisiting set of solvus
c                                 criteria, initialized here.
  
            savg(phase) = .false.

            write (*,1010) isol,isol-1
            
            do i = 1, isol-1

               call mkcomp (i+k5,phase)
c                                 get the range for the compositional
c                                 variable:
5020           write (*,1020) i
               read (*,*,iostat=ier) cmin(i), cmax(i)
               call rerror (ier,*5020)

            end do 

         end if 

      end if 

      if (nopt(8).eq.1d0) then 
c                                 average immiscible compositions
c                                 get mole fractions
         call avgcmp (isol,jdsol)

      else 
c                                 identify the immiscible phase of interest     
c                                 test which phase (if any) satisfy
c                                 the compositional criteria:
         jsol = 0

         do 20 i = 1, isol
c                                 for each phase, test the isol-1 
c                                 conditions:
            do j = 1, isol-1
c                                 comp is a function that returns
c                                 the j+1th composition 
               tcomp = gtcomp (i,idasls(jdsol(i),ias),k5+j)
c                                 the composition is not relevant
               if (tcomp.eq.-1d99) cycle
c                                 the composition is out of bounds
               if (tcomp.lt.cmin(j).or.tcomp.gt.cmax(j)) goto 20

            end do
c                                 the solution past all tests
            jsol = jsol + 1
            kdsol(jsol) = jdsol(i)

20       continue 

         if (jsol.gt.1.and..not.savg(phase)) then 
c                                 two or more phases satisfy the
c                                 existing criteria
            write (*,1060) jsol,fname(phase)
            write (*,1040) (cname(i), i = 1, icomp)
            do i = 1, jsol 
               write (*,1050) (pcomp(j,kdsol(i)), j = 1, icomp)
            end do 

            write (*,1070)
c                                 get choice
            call rdnumb (tcomp,0d0,choice,1,.false.)

            if (choice.eq.2) then 
c                                 2 - average within existing criterion 
               savg(phase) = .true.

            else if (choice.eq.3) then 
c                                 3 - ignore and hope for the best
               call avgcmp (jsol,kdsol)

            else 
c                                 not 2 or 3, redefine the criterion
               stol(phase) = .false.
               goto 10             

            end if 

         end if 
c                                 user has elected to average within
c                                 existing criterion
         if (jsol.gt.1.and.savg(phase)) call avgcmp (jsol,kdsol) 

         if (jsol.eq.0) then 
            index = 0 
         else 
            index = kdsol(1)
         end if 

      end if 

1000  format (/,i1,' immiscible phases of ',a,/,'coexist with the ',
     *        'following ',a,' compositions:',/)
1010  format (/'The following prompts define the compositional ',
     *        'variable(s) (C[i]) to be used',/,
     *        'to identify the phase of interest.',//,
     *        'As there are ',i1,' coexisting phases',
     *        ' you will be prompted',/,'for ',i1,
     *        ' compositional variable(s).',/)
1020  format (/,'Enter the range (minimum, maximum) of C[',i1,'] that ',
     *        'defines the phase of interest:',/)
1030  format (/,'Choose an option:',/,3x,
     *        '1 - specify compositional criteria to identify the',
     *        ' phase of interest [default].',/,3x,
     *        '2 - average the compositions of immiscible phases',/)
1040  format (/,4x,20(a,4x))
1050  format (3x,20(f7.3,2x))
1060  format (/,i1,' coexisting phases of ',a,' satisfy your ',
     *        ' compositional criteria',/,'with the compositions:')
1070  format (/,'Choose an option:',/,3x,
     *        '1 - redefine the compositional criteria [default].',/,3x,
     *        '2 - average the compositions of all phases that',
     *        ' meet the existing criterion.',/,3x,
     *        '3 - ignore this instance.',/)
99    end 

      subroutine mode3 
c----------------------------------------------------------------------
c sample data on an x-y path defined by interactive user input 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node, ok

      integer i, j, icurve, ivi, ivd, iord, ipts, jpts, ier, k(2), dim

      double precision coef(0:10), dxy(2), xyp(2,2), s, d

      character*100 n5name, n6name, yes*1, text*200

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      character vnm*8
      common/ cxt18a /vnm(l3)  

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem
c----------------------------------------------------------------------
      node = .false.
      dim = 1
c                                 set up path information
      icurve = 0 
c                                 ask if non-linear path
10    write (*,1200) 
      read (*,'(a)') yes

      if (yes.eq.'y'.or.yes.eq.'Y') then 

         icurve = 1
c                                 select independent variable:
5        write (*,1160) (i,vnm(i),i= 1, 2)
         read (*,*,err=5) ivi

         if (ivi.eq.1) then 
            ivd = 2
         else
            ivi = 2
            ivd = 1
         end if

         write (*,1210) vnm(ivd),vnm(ivi)
         read (*,*) iord

         do i = 0, iord
            write (*,1220) i
            read (*,*) coef(i)
         end do

         write (text,1350) vnm(ivd),(coef(i),vnm(ivi),i,i=0,iord)
         call unblnk (text)
         write (*,1340) text
c                                 ask if ok.
         write (*,1320)
         read (*,'(a)') yes 

         if (yes.eq.'y'.or.yes.eq.'Y') goto 10
c                                 it's ok.
         dxy(ivi) = vmx(ivi)-vmn(ivi)
         dxy(ivd) = vmx(ivd)-vmn(ivd)
         xyp(ivi,1) = vmn(ivi)
         xyp(ivd,1) = vmn(ivd)

      else 
c                                 linear path
         ivi = 1
         ivd = 2
         i = 1 

         do 

            write (*,1140) i, vnm(1), vnm(2)
            read (*,*,iostat=ier) xyp(1,i), xyp(2,i)
            if (ier.ne.0) cycle

            ok = .true.

            do j = 1, 2

               if (vmn(j).lt.vmx(j)) then 

                  if (xyp(j,i).gt.vmn(j).and.xyp(j,i).lt.vmx(j)) cycle 
 
                  write (*,1010) vnm(j),vmn(j),vmx(j)
                  ok = .false. 

               else

                  if (xyp(j,i).gt.vmx(j).and.xyp(j,i).lt.vmn(j)) cycle  
                  write (*,1010) vnm(j),vmn(j),vmx(j)
                  ok = .false.

               end if
 
            end do 

            if (.not.ok) cycle

            i = i + 1

            if (i.eq.3) then 

               do j = 1, 2 
                  dxy(j) = xyp(j,2) - xyp(j,1)
               end do 

               if (dxy(1).eq.0d0.and.dxy(2).eq.0d0) then

                  write (*,*) 
     *               'initial and final coordinates cannot be identical'
                  i = 1
                  cycle

               else if (dxy(1).eq.0d0) then 
c                                 parallel to the x-axis
                     ivi = 2
                     ivd = 1
                     s = 0d0

               else if (dxy(2).eq.0d0) then 
c                                 parallel to y axis
                  s = 0 

               end if 

               exit 

            end if 

         end do 

      end if 
c                                 set up counters, pointers:
      do 

         write (*,1150) 
         read (*,*,iostat=ier) ipts
         if (ipts.lt.2) ipts = 2
         if (ier.eq.0) exit 

      end do 

      d = dxy(ivi)/dfloat(ipts - 1)
c                                 ind could be set by getind              
      ind = ivi

      jpts = 0 
c                                 select properties
      call chsprp 
c                                 write file header
      call tabhed (dxy,dxy,k,dim,n5name,n6name)
c                                 compute properties
      do i = 1, ipts

         var(ivi) = xyp(ivi,1) + dfloat(i-1)*d

         if (icurve.eq.0) then 

            var(ivd) = xyp(ivd,1) + s*(var(ivi)-xyp(ivi,1))
            jpts = 1

         else
 
            var(ivd) = 0d0

            do j = 0, iord
               var(ivd) = var(ivd) + coef(j)*var(ivi)**j
            end do

            if (var(ivd).le.vmx(ivd).and.var(ivd).ge.vmn(ivd)) then 
               jpts = jpts + 1
            else 
               exit  
            end if
 
         end if 

         call polprp (dim)

      end do 

      call finprp (dim,n5name,n6name,node) 

      if (jpts.eq.0) write (*,1330) 

1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1140  format (/,'Enter endpoint ',i1,' (',a,'-',a,') coordinates:')
1150  format (/,'How many points along the profile?')
1160  format (/,'Select independent variable: ',2(/,1x,i1,' - ',a))
1200  format (/,'Construct a non-linear profile (y/n)?')
1210  format (/,'Profile must be described by the function',/,a,
     *        ' = Sum ( c(i) * ',a,' ^i, i = 0..n)',/,'Enter n (<10)')
1220  format (/,'Enter c(',i2,')')
1320  format (/,'Change the profile (Y/N)?')
1330  format (/,'Your polynomial does not yield conditions within',
     *          'computational coordinate frame.',/)
1340  format (/,'Your polynomial is:',/,a)
1350  format (a,'=',5('+(',g12.6,')','*',a,'^',i1))

      end 

      subroutine mode31
c----------------------------------------------------------------------
c sample data in a 1d section 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node

      integer i, ipts, ier, k(2), dim 

      double precision dxy, xyp(2), d

      character*100 n5name, n6name

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      character vnm*8
      common/ cxt18a /vnm(l3)  

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem
c----------------------------------------------------------------------
      node = .false.
      dim = 1
c                                 get the independent output variable
      call getind
c                                 path endpoints
30    do i = 1, 2

20       if (i.eq.1) then 
            write (*,1130) vnm(1)
         else
            write (*,1140) vnm(1)
         end if 

         read (*,*,err=20) xyp(i)

         if (vmn(1).lt.vmx(1)) then 
            if (xyp(i).lt.vmn(1).or.xyp(i).gt.vmx(1)) then  
               write (*,1010) vnm(1),vmn(1),vmx(1)
               goto 20
            end if 
         else
            if (xyp(i).lt.vmx(1).or.xyp(i).gt.vmn(1)) then  
               write (*,1010) vnm(1),vmn(1),vmx(1)
               goto 20
            end if 
         end if 

      end do 

      dxy = xyp(2) - xyp(1)

      if (dxy.eq.0d0) then
         write (*,*) 'initial and final coordinates cannot be identical'
         goto 30
      end if   
c                                 set up counters, pointers:
      do 
         write (*,1150) 
         read (*,*,iostat=ier) ipts
         if (ipts.lt.2) ipts = 2
         if (ier.eq.0) exit 
      end do

      d = dxy/dfloat(ipts - 1)
c                                 select properties:
      call chsprp
c                                 name and open plot file, write header 
      call tabhed (xyp,xyp,k,1,n5name,n6name)

      do i = 1, ipts

         var(1) = xyp(1) + dfloat(i-1)*d

         call polprp (dim)

      end do 

      call finprp (dim,n5name,n6name,node) 

1010  format (/,'The plot file range for ',a,' is ',g12.4,' - ',g12.4,
     *        /,'Try again:',/)
1130  format (/,'Enter the ',a,' coordinate at the beginning of',
     *          ' the profile:')
1140  format (/,'Enter the ',a,' coordinate at the end of',
     *          ' the profile:')
1150  format (/,'How many points along the profile?')

      end 

      subroutine mode4
c----------------------------------------------------------------------
c sample data on an x-y path defined by input from a file
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical node

      integer i, icoors, jmode, ixy, inc, ierr, k(2), dim

      double precision r(2), tmin1, tmax1,a1, b1, c1, d1, x0, x1,  
     *                 dt1, x, y, xx(5*l5), yy(5*l5), pmin, pmax

      character*100 n5name, n6name, dname

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar 

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      logical oned
      common/ cst82 /oned
c----------------------------------------------------------------------
      node = .false.
      dim = 1

      do 

         if (oned) then 

            write (*,1280) 
            jmode = 2 

         else 
c                                 select input type
            write (*,1290) 
            read (*,*) jmode   
            write (*,1300) 

         end if 

         read (*,'(a)') dname
         open (n7,file=dname,iostat=ierr)
         if (ierr.eq.0) exit 
         write (*,*) 'No such data file as: ',dname

      end do 

      if (jmode.eq.1) then 
c                                 points from a polynomial
c                                 the association of "x" and "y"
c                                 with the actual variables of the
c                                 diagram is determined by the flag
c                                 ixy 
         read (n7,*) pmin,pmax,ixy
         read (n7,*,iostat=ierr) tmin1,tmax1,dt1,a1,b1,c1,d1
         if (ierr.ne.0) return

         if (dt1.lt.0d0) then 
            x0 = tmax1
            x1 = tmin1
         else 
            x0 = tmin1
            x1 = tmax1
         end if 

         x = x0
c                                 select properties
         call chsprp 
c                                 write plot file header
         call tabhed (r,r,k,dim,n5name,n6name)

         do 

            y = a1 + b1*x + c1*x**2 + d1*x**3

            if (y.le.pmax.and.y.ge.pmin) then 
               if ((dt1.gt.0d0.and.x.le.x1).or.
     *             (dt1.lt.0d0.and.x.ge.x1)) then
c                                 condition is in bounds
                  if (ixy.eq.0) then 
                     var(1) = x
                     var(2) = y 
                  else
                     var(1) = y
                     var(2) = x 
                  end if 

                  call polprp (dim) 

               end if 

            end if 

            x = x + dt1

            if (x.gt.x1) exit 

         end do  

         call finprp (dim,n5name,n6name,node) 

      else 
c                                   points from a data file:
         icoors = 1

         do 

            if (oned) then
   
               read (n7,*,iostat=ierr) xx(icoors)
               if (ierr.ne.0) exit 

               yy(icoors) = vmn(2)

            else

               read (n7,*,iostat=ierr) xx(icoors),yy(icoors)
               if (ierr.ne.0) exit 

            end if 

            if (xx(icoors).lt.vmn(1).or.xx(icoors).gt.vmx(1)) cycle 
            if (yy(icoors).lt.vmn(2).or.yy(icoors).gt.vmx(2)) cycle 

            icoors = icoors + 1

            if (icoors.gt.5*l5) call error (69,xx(1),i,'MODE4') 

         end do 
               
         icoors = icoors - 1

         close (n7)

         if (icoors.eq.0) then 
            write (*,*) 'file contains no points within data bounds'
            return
         end if 

         write (*,1310) icoors
         read (*,*) inc
c                                 select properties            
         call chsprp
c                                 write plot file header
         call tabhed (r,r,k,dim,n5name,n6name)

         do i = 1, icoors, inc

            var(1) = xx(i)
            var(2) = yy(i)

            call polprp (dim)

         end do  

         call finprp (dim,n5name,n6name,node) 

      end if 

1280  format (/,'Enter the name of the file containing the path',
     *          ' ordinates:',/)
1290  format (/,'Path will be described by:',/,
     *          '   1 - a file containing a polynomial function',/,
     *          '   2 - a file containing a list of x-y points',/,
     *          'Enter 1 or 2:'/) 
1300  format (/,'Enter the file name:',/)
1310  format (/,'File contains ',i5,' points',/,
     *          'every nth plot will be plotted, enter n:',/)

      end

      subroutine avgcmp (isol,jdsol)  
c-------------------------------------------------------------------
c makes the average composition and properties of from isol compositions
c of the phase indexed by the array jdsol
c-------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h' 

      integer jdsol(k5), index, isol, i, j

      double precision x(k5), ntot

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------
c                                 average immiscible compositions
c                                 get mole fractions
      ntot = 0d0 
      index = jdsol(1)

      do i = 1, isol
         ntot = ntot + props(16,jdsol(i))
      end do 

      do i = 1, isol
         x(i) = props(16,jdsol(i))/ntot
      end do 
c                                 set composition
      do j = 1, icomp

         pcomp(j,index) = x(1) * pcomp(j,index)

         do i = 2, isol
            pcomp(j,index) = pcomp(j,index) + x(i)*pcomp(j,jdsol(i))
         end do 
      end do 
c                                 set physical properties assuming molar
c                                 weighting (this is probably wrong for 
c                                 some). 
      do j = 1, 17

         if (j.eq.16) cycle

         props(j,index) = x(1) * props(j,index)

         do i = 2, isol 
            props(j,index) = props(j,index) + x(i)*props(j,jdsol(i))
         end do 

      end do 

      end

      subroutine xy2ij (iloc,jloc,left)
c----------------------------------------------------------------------
c given the current x-y coordinates (loaded in var), return the nodal
c coordinates
c----------------------------------------------------------------------
      implicit none

      integer jloc, iloc

      logical left

      logical oned
      common/ cst82 /oned


      if (oned) then 
c                                 get the nodal coordinate and find if the
c                                 real coordinate lies to the left or right 
c                                 of the nodal coordinate.                             
         call amiin1 (jloc,left)
         iloc = 1

      else
c                                 could get fancy and record up/left here
         call amiin2 (iloc,jloc) 

      end if 

      end 

      subroutine allmod 
c----------------------------------------------------------------
c computes modes of all stable phases, i.e., over the entire range
c of physical conditions. 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, id, jk

      double precision mode(3)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab

      integer idsol,nrep,nph
      common/ cst38/idsol(k5,k3),nrep(k5,k3),nph(k3)
c----------------------------------------------------------------------
      do i = 1, iprop
         prop(i) = 0d0
      end do 

      id = 0 

      do i = 1, nph(ias)

         jk = 0 

         do j = 1, istab

            if (idstab(j).eq.idsol(i,ias)) then 

               do k = 1, nrep(i,ias)

                  id = id + 1

                  call gtmode (mode,id)
                  prop(jk+k) = mode(iopt(3)+1)

               end do 

            end if
c                                mode column pointer
            jk = jk + nstab(j)
 
         end do

      end do

      if (lopt(2)) then
c                                 convert to cumulative modes if
c                                 requested
         do j = 2, iprop
            prop(j) = prop(j) + prop(j-1)
         end do

      end if
  
      end 

      subroutine gtmode (mode,id) 
c----------------------------------------------------------------
c function to extract vol/wt/mol mode (%) of phase id 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      double precision mode(3)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
c----------------------------------------------------------------
      if (aflu.and.lflu.or.(.not.aflu).or.psys1(1).eq.0d0) then
c                     total mode:
c                     volume fraction
         mode(1) = props(1,id)*props(16,id)/psys(1)*1d2
c                     weight fraction 
         mode(2) = props(16,id)*props(17,id)/psys(17)*1d2
c                     mol fraction
         mode(3) = props(16,id)/psys(16)*1d2

      else 
c                     solid only mode:
         if (fluid(id)) then 
            mode(1) = 0d0
            mode(2) = 0d0
            mode(3) = 0d0
         else 
c                     volume fraction
            mode(1) = props(1,id)*props(16,id)/psys1(1)*1d2
c                     wt fraction
            mode(2) = props(16,id)*props(17,id)/psys1(17)*1d2
c                     mol fraction
            mode(3) = props(16,id)/psys1(16)*1d2
         end if 

      end if 

      end 

      subroutine outmod (dim,n6name,node)
c----------------------------------------------------------------
c reformat output from "all_modes" requests according to the 
c options set in perplex_option.dat and dimension dim. 
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*8 ynm, n6name*100, rec*1

      logical node

      integer i, j, k, ipt, dim, ier

      double precision ymx,ymn,xl(i11),yl(i11),x(3),dx,dy(i11),xmx,xmn

      character vnm*8
      common/ cxt18a /vnm(l3) 

      double precision vip
      common/ cst28 /vip(l2,k2)

      character*162 title
      common/ csta8 /title(4)

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      integer inv
      character dname*14, titl1*162
      common/ cst76 /inv(i11),dname(i11),titl1

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
      if (dim.eq.1) then 

         rewind (n5)
c                                 skip header lines
         do i = 1, 9
            read (n5,'(a)') rec
         end do 

         open (n6,file=n6name)
c                                 read file to get limits, this 
c                                 must also be done for x if 
c                                 ind ne 1.
         ymn = 1d30
         ymx = -1d30
         xmn = ymn
         xmx = ymx
         ynm = 'y axis'
         ipt = 0 

         do i = 1, iprop
            dy(i) = -1d30
            yl(i) = -1d30
         end do 
c                                 read the data to get the range  
         do           
                  
            read (n5,*,iostat=ier) (x(i),i=1,ivar), (prop(i),i=1,iprop)

            ipt = ipt + 1
            if (ipt.gt.k2) call error (999,ymn,ipt,'OUTMOD')

            if (node) then 
               vip(4,ipt) = ipt
            else 
               vip(4,ipt) = x(ind)
            end if 

            if (ier.ne.0) exit 

            do i = 1, iprop

               if (x(ind).gt.xmx) xmx = x(ind)
               if (x(ind).lt.xmn) xmn = x(ind)
               if (prop(i).gt.ymx) ymx = prop(i)
               if (prop(i).lt.ymn) ymn = prop(i)
c                                 get a label position
               if (lopt(2)) then 
c                                 cumulative mode
                  if (i.eq.1) then 

                     if (prop(i).gt.2d0*yl(i)) then
                        yl(i) = prop(i)/2d0
                        xl(i) = x(ind)
                     end if 

                  else 

                     if (prop(i)-prop(i-1).gt.dy(i)) then
                        dy(i) = prop(i)-prop(i-1)
                        yl(i) = (prop(i)+prop(i-1))/2d0 
                        xl(i) = x(ind)
                     end if
 
                  end if 

               else

                  if (prop(i).gt.yl(i)) then 
                     yl(i) = prop(i)
                     xl(i) = x(ind)
                  end if 

               end if 

            end do 

         end do

         dx = (xmx-xmn)/1d1
c                                 header of psvdraw file
         write (n6,1000) title,xmx,xmn,ymx,ymn,vnm(ind),ynm
c                                 make each column into a curve use of
c                                 the vip array might cause problems 
c                                 if someone uses chemical potential
c                                 variables and calculates multiple
c                                 property sets, as vip(4,5) will be 
c                                 reset and are not reread. 
         do i = 1, iprop

            rewind (n5)
c                                 skip header lines
            do j = 1, 9
               read (n5,'(a)') rec
            end do 

            do j = 1, ipt       
                   
               read (n5,*,iostat=ier) (x(k),k=1,ivar),
     *                                (prop(k),k=1,iprop)
               vip(5,j) = prop(i)

            end do
c                                 output the curve
            write (n6,1010) ipt*2,i,dname(i)
            write (n6,*) (vip(4,j),vip(5,j),j=1,ipt)

         end do
c                                 write trailers for psipts, pscurv
         write (n6,1010) 1,1,'trailer'
         write (n6,*) '0 0'
c                                 write label coordinates and text
         do i = 1, iprop
            if (lopt(2)) then 
               write (n6,*) xl(i)-dx,yl(i)
            else
               write (n6,*) xl(i),yl(i)
            end if 
            write (n6,'(a)') dname(i)
         end do 

      end if

      close (n6)

1000  format ('1',/,'0 0 0',/,'0 0 0 0 0 0',4(/,a162),/,'2 1 2 0 0',/,
     *        '0 0 0 0. 0. 0. 0. 0.',/,4(g14.7,1x),/,a,/,a)
1010  format (i5,' 1 ',i3,' 1 1 1 1 1 1 ',/,'0. ',a)

      end 

      subroutine allprp (dim)
c----------------------------------------------------------------
c output spreadsheet format properties (lop = 36 or 38).
c   icx = 0   - output just system properties
c   icx = 999 - system + phase properties
c   else      - phase icx properties
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, dim  

      double precision mode(3)

      double precision gtot,fbulk,gtot1,fbulk1
      common/ cxt81 /gtot,fbulk(k0),gtot1,fbulk1(k0)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision props,psys,psys1,pgeo,pgeo1
      common/ cxt22 /props(i8,k5),psys(i8),psys1(i8),pgeo(i8),pgeo1(i8)

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer jvar
      double precision var,dvr,vmn,vmx
      common/ cxt18 /var(l3),dvr(l3),vmn(l3),vmx(l3),jvar

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision pcomp
      common/ cst324 /pcomp(k0,k5)

      double precision mu
      common/ cst330 /mu(k8)

      character pname*14
      common/ cxt21a /pname(k5)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab
c----------------------------------------------------------------------
      if (kop(1).eq.38) then 
c                                 custom property choices
         if (kcx(1).eq.999.or.kcx(1).eq.0) then 
c                                 get system props
            do i = 1, iprop
               call getprp (prop(i),kop(i+1),0,0,.true.)
            end do

            tname = 'system        '

            call outprp (dim)

         end if 

         if (kcx(1).ne.0) then 
c                                 properties of phases
            do j = 1, ntot
             
               if (kcx(1).eq.999) then
                  id = j 
               else
                  id = kcx(1)
               end if 

               do i = 1, iprop
                  call getprp (prop(i),kop(i+1),id,0,.true.)
               end do

               tname = pname(id)

               call outprp (dim)

               if (kcx(1).ne.999) exit 

            end do       

         end if  

      else
c                                 "all property" option, lop.eq.36
         if (kcx(1).eq.999.or.kcx(1).eq.0) then 
c                                 only system properties requested:
            if (aflu.and.lflu.or.(.not.aflu)) then 
c                                 normal properties
               do i = 1, i8
                  prop(i) = psys(i)
               end do 
c                                 modes
               do i = i8+1, i8+3
                  prop(i) = 1d2
               end do 
c                                 bulk composition 
               do i = i8+4, i8+3+icomp
                  prop(i) = fbulk(i-i8-3)
               end do 
c                                 chemical potentials
               do i = i8+icomp+4, i8+3+icomp+ichem
                  prop(i) = mu(i-i8-3-icomp)
               end do

               tname = 'system        '

               call outprp (dim)

            else
c                                 exclude fluid:
c                                 normal properties
               do i = 1, i8
                  prop(i) = psys1(i)
               end do 
c                                 modes
               do i = i8+1, i8+3
                  prop(i) = 1d2
               end do 
c                                 bulk composition 
               do i = i8+4, i8+3+icomp
                  prop(i) = fbulk1(i-i8-3)
               end do 
c                                 chemical potentials
               do i = i8+icomp+4, iprop
                  prop(i) = mu(i-i8-3-icomp)
               end do 

               tname = 'system        '

               call outprp (dim)

            end if 

            if (kcx(1).ne.0) then
c                                 properties of all phases
               do j = 1, ntot

                  if (kcx(1).eq.999) then 
                     id = j 
                  else
c                                 all properties of a specific phase
c                                 find the phase index
                     call soltst (id,kcx(1))

                  end if 
c                                 normal properties
                  do i = 1, i8
                     prop(i) = props(i,id)
                  end do 
c                                 compute modes
                  call gtmode (mode,id)

                  do i = 1, 3
                     prop(i8+i) = mode(i)
                  end do 
c                                 bulk composition 
                  do i = i8+4, i8+3+icomp
                     prop(i) = pcomp(i-i8-3,id)
                  end do 
c                                 chemical potentials
                  do i = i8+icomp+4, iprop
                     prop(i) = mu(i-i8-3-icomp)
                  end do 

                  tname = pname(id)

                  call outprp (dim) 

                  if (kcx(1).ne.999) exit 

               end do

            end if 

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
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

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

      subroutine grxn (g)
c--------------------------------------------------------------------
c a dummy routine to allow rk to be linked with rlib.f
c--------------------------------------------------------------------
      implicit none
      double precision g
      g = g
      end

      subroutine tabhed (vmn,dv,nv,nvar,n5name,n6name)
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

      integer i, nv(2), nvar, ivar

      double precision vmn(2), dv(2)

      character*100 n6name, n5name, vname(l3)*14

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      character vnm*8
      common/ cxt18a /vnm(l3)  
c------------------------------------------------------------------------
c                                 generate a file name and
c                                 open the file on n5
      call fopenn (n5,nvar,n5name,n6name)
c                                 ctr file tab format header
      write (n5,'(a)') '|6.6.6'
c                                 a title
      write (n5,'(a)') n5name
c                                 number of indepedent variables
      write (n5,*) nvar
c                                 independent variable names, 
c                                 value, increment & nodes
      do i = 1, nvar
         write (n5,'(a)') vnm(i)
         write (n5,*) vmn(i)
         write (n5,*) dv(i)
         write (n5,*) nv(i)
      end do 
c                                 number of pseudo-dependent variables,
      ivar = 2
      if (icopt.eq.10) ivar = 3
c                                 convert a8 names to a14
      do i = 1, ivar
         vname(i) = vnm(i)
         call unblnk (vname(i))
      end do 
c                                 output the dependent variable counter and list
      if (kcx(1).eq.999) then 
c                                  phemgp file
         write (n5,*) ivar + iprop + 2
         write (n5,'(200(a14,1x))') 'Name','Counter',
     *                              (vname(i), i = 1, ivar),
     *                              (dname(i),i = 1, iprop)

      else 
c                                  tab file
         if (lopt(15).or.nvar.eq.1) then
c                                  with pseudo-dependent variables 
            write (n5,*) ivar + iprop 
            write (n5,'(200(a14,1x))') (vname(i), i = 1,ivar),
     *                                 (dname(i), i = 1,iprop)

         else 
c                                  terse format
            write (n5,*) iprop 
            write (n5,'(200(a14,1x))') (dname(i), i = 1,iprop)
          
         end if 

      end if

      end 

      subroutine getind
c----------------------------------------------------------------
c get the plotting variable index (ind) for 1-d property plots
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ier

      character vnm*8
      common/ cxt18a /vnm(l3) 

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem
c----------------------------------------------------------------------
c                                 choose plotting variable
      write (*,1000) vnm(1)

      do 
         write (*,1030) (ier,vnm(ier),ier=1,ivar)
         read (*,'(bn,i80)',iostat=ier) ind
         if (ier.ne.0) cycle
         if (ind.ne.2.and.ind.ne.3) ind = 1
         exit 
      end do  

1000  format (/,'The independent variable for this calculation is: ',a)
1030  format (/,'Choose the independent variable for data plots:',/,
     *       4x,i1,' - ',a,' [default]',6(/,4x,i1,' - ',a))
 
      end 

      subroutine chsprp 
c----------------------------------------------------------------
c chsprp asks the user to choose properties to be extracted, it
c then creates a list of names for the properties (dname).

c   lop/kop  - flag indicating the property chosen
c   icx/kcx  - if lop = 6, the component chosen
c   icx/kcx  - if lop > 6, the identity of the solution chosen,
c              icx = -1 if a solution is not chosen.

c   lflu/kfl - .true. include fluids for bulk props.

c 1                 Specific enthalpy (J/m3)',
c 2                 Density (kg/m3)',
c 3                'Specific Heat capacity (J/K/m3)',
c 4                'Expansivity (1/K, for volume)',
c 5                'Compressibility (1/bar, for volume)',
c 6                'Weight percent of a component',
c 7                'Mode (Vol %) of a compound or solution',
c 8                'Composition of a solution'
c 9                 Grueneisen thermal ratio',
c 10               'Adiabatic bulk modulus (bar)',
c 11               'Adiabatic shear modulus (bar)'
c 12                Sound velocity (km/s)
c 13               'P-wave velocity (km/s)',
c 14                S-wave velocity (km/s)',
c 15                Vp/Vs
c 16               'Specific Entropy (J/K/m3)'
c 17               'Entropy (J/K/kg)'
c 18               'Enthalpy (J/kg)'
c 19               'Heat Capacity (J/K/kg)'
c 20               'Specific mass (kg/m3) of a phase'
c 21               'Poisson's Ratio'
c 22               'Molar Volume (J/bar)'
c 23                Dependendent potentials (J/mol)
c 24                Assemblage index
c 25                Modes of all phases (wt or vol%)
c 26                Sound velocity temperature derivative (km/s/K)
c 27                P-wave velocity temperature derivative (km/s/K)
c 28                S-wave velocity temperature derivative (km/s/K)
c 29                Adiabatic bulk modulus temperature derivative (bar/K)
c 30                Shear modulus temperature derivative (bar/K)
c 31                Sound velocity pressure derivative (km/s/bar)
c 32                P-wave velocity pressure derivative (km/s/bar)
c 33                S-wave velocity pressure derivative (km/s/bar)
c 34                Adiabatic bulk modulus pressure derivative (unitless)
c 35                Shear modulus pressure derivative (unitless)
c 36                All properties of a phase or the system
c 37                Absolute amounts
c 38                Multiple property grid for system and phases
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, icx, kprop, ier, lop, komp, mprop

      parameter (kprop=38)

      character propty(kprop)*60, y*1, pname*10

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer ivar,ind,ichem
      common/ cst83 /ivar,ind,ichem

      character cname*5
      common/ csta4  /cname(k5)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer idstab,nstab,istab
      common/ cst34 /idstab(i11),nstab(i11),istab

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      integer icps, jcx, jcx1, kds
      logical stol, savg
      double precision rcps
      common/ comps /rcps(k7,2*k5),icps(k7,2*k5),jcx(2*k5),jcx1(2*k5),
     *               kds(2*k5),stol(h9),savg(h9)

      save propty

      data propty/'Specific Enthalpy (J/m3)',
     *            'Density (kg/m3)',
     *            'Specific heat capacity (J/K/m3)',
     *            'Expansivity (1/K, for volume)',
     *            'Compressibility (1/bar, for volume)',
     *            'Weight (%) of a component',
     *            'Mode (Vol, Mol, or Wt proportion) of a phase',
     *            'Composition (Mol or Wt) of a solution',
     *            'Grueneisen thermal ratio',
     *            'Adiabatic bulk modulus (bar)',
     *            'Adiabatic shear modulus (bar)',
     *            'Sound velocity (km/s)',
     *            'P-wave velocity (Vp, km/s)',
     *            'S-wave velocity (Vs, km/s)',
     *            'Vp/Vs',
     *            'Specific entropy (J/K/m3)',
     *            'Entropy (J/K/kg)',
     *            'Enthalpy (J/kg)',
     *            'Heat Capacity (J/K/kg)',
     *            'Specific mass of a phase (kg/m3-system)',
     *            'Poisson ratio','Molar Volume (J/bar)',
     *            'Dependent potentials (J/mol, bar, K)',
     *            'Assemblage Index',
     *            'Modes of all phases',
     *            'Sound velocity T derivative (km/s/K)',
     *            'P-wave velocity T derivative (km/s/K)',
     *            'S-wave velocity T derivative (km/s/K)',
     *            'Adiabatic bulk modulus T derivative (bar/K)',
     *            'Shear modulus T derivative (bar/K)',
     *            'Sound velocity P derivative (km/s/bar)',
     *            'P-wave velocity P derivative (km/s/bar)',
     *            'S-wave velocity P derivative (km/s/bar)',
     *            'Adiabatic bulk modulus P derivative (unitless)',
     *            'Shear modulus P derivative (unitless)',
     *            'All phase &/or system properties (PHEMGP format)',
     *            'Absolute amount (Vol, Mol, or Wt) of a phase',
     *            'Multiple property output (PHEMGP format)'/
c----------------------------------------------------------------------
      do i = 1, istab

         if (stol(i)) then 
c                                 doing a second run, with an 
c                                 existing solvus criterion, ask
c                                 whether to change.
            write (*,1030)
            read (*,'(a)') y
            if (y.eq.'y'.or.y.eq.'Y') stol(i) = .false.
    
         end if
 
      end do 
c                                 property counter
      iprop = 0
c                                 phase composition counter
      komp = 0
c                                 choose property
      do 

         icx = 0
         lflu = .false.

         write (*,1050)

         do i = 1, kprop

            if (iprop.gt.0.and.i.eq.25.or.i.eq.36) cycle 
            write (*,1060) i,propty(i)

         end do 

         do 

            call rdnumb (prop(1),0d0,lop,999,.false.)
            if (lop.ne.999) exit
            write (*,'(a)') 'Select a property or enter 0 to finish...'

         end do 
  
         if (lop.lt.0.or.lop.gt.kprop) then 

            write (*,1020)
            cycle 

         else if (lop.eq.0) then 

            exit 

         end if

         if (lop.eq.7.or.lop.eq.20.or.lop.eq.37) then 
c                                 modes:
c                                 get phase name
             call rnam1 (icx,pname)
c                                 ask if fluid should be included:
             if (gflu) then 

                write (*,1120) 
                read (*,'(a)') y
                if (y.eq.'y'.or.y.eq.'Y') lflu = .true.

             end if 
c                                 write blurb about units
             if (lop.eq.7) then 
                write (*,1080)
             else if (lop.eq.37) then 
                write (*,1090)
             end if 
             
         else if (lop.eq.25) then 
c                                 all modes
             write (*,1070)
             read (*,'(a)') y

             if (y.eq.'y'.or.y.eq.'Y') then 
                lopt(2) = .true.
             else
                lopt(2) = .false.
             end if 
c                                 ask if fluid should be included:
             if (gflu) then 

                write (*,1120) 
                read (*,'(a)') y
                if (y.eq.'y'.or.y.eq.'Y') lflu = .true.

             end if 
c                                 double loop necessary because solution 
c                                 i may occur as j coexisting phases
             do i = 1, istab
                do j = 1, nstab(i)

                   iprop = iprop + 1
                   call getnam (dname(iprop),idstab(i))

                end do 
             end do 
 
         else if (lop.eq.6.or.lop.eq.23) then
c                                 warn if no potentials
            if (jpot.eq.1.and.lop.eq.23) then

               call warn (31,nopt(1),iopt(1),'CHSPRP')
               cycle

            end if 
c                                 get component to be contoured
5010        write (*,1000)

            if (lop.eq.23) then 
c                                 for usv calculations make names
c                                 for the extra potentials (p,t)                            
               if (hcp.gt.icp) then 
                  ichem = hcp 
                  cname(icp+1) = 'T(K)   '
                  cname(icp+2) = '-P(bar)'
               else 
                  ichem = icp 
               end if 

               write (*,1010) (i, cname(i), i = 1, ichem)

            else 

               write (*,1010) (i, cname(i), i = 1, icomp)

            end if 

            read (*,*,iostat=ier) icx
            call rerror (ier,*5010)   
c                                 ask if fluids included
            if (gflu.and.lop.eq.6) then 

               write (*,1120) 
               read (*,'(a)') y
               if (y.eq.'y'.or.y.eq.'Y') lflu = .true. 

            end if 

         else if (lop.eq.8) then
c                                 get solution identity
            do 

               call rnam1 (icx,pname)
               if (icx.gt.0) exit  
               write (*,1140)

            end do 
c                                 get user defined composition:
            komp = komp + 1

            if (komp.gt.k5) then 
               write (*,1160) k5
               cycle 
            end if 

            call mkcomp (komp,icx)

         else if (lop.eq.36.or.lop.eq.38) then   
c                                 multi-prop options, get case i: 
c                                 1 - system, 2 - phase
c                                 3 - system + phases
            write (*,1130)
            call rdnumb (nopt(1),0d0,i,1,.false.)
c                                 convert kop to the internal 
c                                 icx flag => icx = 0 system prop,
c                                 icx = 999 all props, else phase index
            if (i.eq.3) then 

               icx = 999 

            else if (i.eq.2) then 
c                                 get phase index
               call rnam1 (icx,pname)

            end if 

            if (gflu) then 
               write (*,1120) 
               read (*,'(a)') y
               if (y.eq.'y'.or.y.eq.'Y') lflu = .true. 
            end if         

            if (lop.eq.36) then 

               iprop = i8 + 3 + icomp + ichem
               if (iprop.gt.i11) call error (1,0d0,iprop,'I11')

            else 
c                                 lop.eq.38
               kop(1) = 38
c                                 custom list
               do

                  write (*,1150)
                  read (*,*,iostat=ier) i

                  if (ier.ne.0.or.i.gt.kprop-1.or.i.lt.0) then
                     write (*,1020)
                     cycle
                  else if (i.eq.8.or.i.eq.6.or.i.eq.23.or.i.eq.24.or.
     *                     i.eq.25.or.i.eq.36.or.i.eq.38) then 
                     write (*,1100) 
                     cycle
                  else if (i.eq.0) then 
                     exit
                  end if 
c                                 save property choice
                  iprop = iprop + 1
                  if (iprop.gt.i11) call error (1,0d0,iprop,'I11')
                  kop(iprop+1) = i                     

               end do

            end if 

         else if (lop.ne.6.and.lop.ne.8) then
               
            if (icx.eq.0) then   
c                                 ask if bulk or phase property
               write (*,1110) 
               read (*,'(a)') y

               if (y.ne.'y'.and.y.ne.'Y') then 
c                                 it's a bulk property, ask if fluid
c                                 should be included:
                  if (gflu) then 
                     write (*,1120) 
                     read (*,'(a)') y
                     if (y.eq.'y'.or.y.eq.'Y') lflu = .true. 
                  end if 

               else if (lop.ne.24) then 
c                                 get phase name
                  do 
                     call rnam1 (icx,pname)
                     if (icx.lt.1.and.lop.eq.8) then
                        write (*,1140)
                        cycle 
                     end if 
                     exit 
                  end do  

               end if

            end if 

         end if 
c                                 make dependent variable names:
         if (lop.eq.25.or.lop.eq.36.or.lop.eq.38) then
c                                 multi prop options, only allowed as 
c                                 single choices.
            kop(1) = lop
            kcx(1) = icx
            kfl(1) = lflu
c                                 assign property names
c                                 lop = 25 -> all mode names are assigned above
            if (lop.eq.36) then 
c                                 "all" prop option
               mprop = i8 + 3
c                                 basic props
               do i = 1, mprop
                  call gtname (lop,icx,i,komp,pname)
               end do 
c                                 for all prop option make the
c                                 bulk composition and chemical 
c                                 potential variable names
               do i = mprop + 1, mprop + icomp
c                                 bulk compositions, lop = 6
                  call gtname (6,i-mprop,i,komp,pname)
               end do 

               do i = mprop + icomp + 1, iprop
c                                 chemical potentials, lop = 23
                  call gtname (23,i-mprop-icomp,i,komp,pname)
               end do 

            else if (lop.eq.38) then 
c                                 "custom" prop option, kop
c                                 pointer is offset because 
c                                 kop(1) = 38
               if (icx.eq.999) pname = 'phase     '

               do i = 1, iprop
                  call gtname (kop(i+1),icx,i,komp,pname)
               end do 

            end if 

            exit 

         else 
c                                 save the local choice options in the 
c                                 global arrays
            iprop = iprop + 1
            kop(iprop) = lop
            kcx(iprop) = icx
            kfl(iprop) = lflu
            k2c(iprop) = komp
c                                 make the name of the property and save
c                                 it in array dname
            call gtname (lop,icx,iprop,komp,pname)

         end if 

      end do 

1000  format (/,'Enter a component:')
1010  format (2x,i2,' - ',a5)
1020  format (/,'Invalid input, try again...',/)
1030  format (/,'Retain the compositional criteria you defined ',
     *          'earlier (y/n)?',/,'Answer yes only if you intend ',
     *          'to extract properties for the same phase.',/)
1050  format (/,'Select properties [enter 0 to finish]:')
1060  format (3x,i2,' - ',a60)
1070  format (/,'Output cumulative modes (y/n)?',/
     *         ,'(see www.perplex.ethz.ch/perplex_options.html'
     *         ,'#cumulative_modes)')
1080  format (//,'Fractions are Wt, Vol, or Mol depending on the '
     *         ,'perplex_option.dat proportions keyword.',//)
1090  format (//,'Amounts are kg, m3, or Mol per unit quantity system '
     *         ,'as specified by the',/
     *         ,'perplex_option.dat proportions keyword.',/)
1100  format (/,'Property not allowed for this option, try again...',/)
1110  format (/,'Calculate individual phase properties (y/n)?')
1120  format (/,'Include fluid in computation of aggregate ', 
     *          '(or modal) properties (y/n)?')
1130  format (/,'In this mode you may tabulate:',
     *     /,4x,'1 - properties of the system',
     *     /,4x,'2 - properties of a phase',   
     *     /,4x,'3 - properties of the system and its phases',/
     *         ,'Output for option 1 & 2 can be plotted with '
     *         ,'PSPLOT, PYWERAMI or MatLab.',/,'Output for '
     *         ,'option 3 can only be plotted with PHEMGP.',//
     *         ,'Select an option [default = 1]:')
1140  format (/,'Hey cowboy, that warnt no solution, try again.',/)
1150  format (/,'Specify a property to be computed from the ',
     *          'list above [0 to end]')
1160  format (/,'**warning ver011** only ',i2,' user defined '
     *      'compositions permitted.',/,'do multiple runs with WERAMI',
     *      'or redimension common block comps.',/)
  
      end

      subroutine gtname (lop,icx,jprop,komp,pname)
c----------------------------------------------------------------
c makes the name of property iprop and saves it in dname(iprop)
c called only by chsprp, variable names as in chsprp.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer icx, jprop, lop, komp, l2p(38)

      character prname(44)*14, pname*10

      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname

      integer inv
      character dname*14, title*162
      common/ cst76 /inv(i11),dname(i11),title

      character cname*5
      common/ csta4  /cname(k5)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      save prname,l2p

      data prname/'V,J/bar/mol   ','H,J/mol       ','Gruneisen_T   ',
c                                 4-6
     *      'Ks,bar        ','Gs,bar        ','v0,km/s       ',
c                                 7-9
     *      'vp,km/s       ','vs,km/s       ','vp/vs         ',
c                                 10-12
     *      'rho,kg/m3     ','G,J/mol       ','cp,J/K/mol    ',
c                                 13-15
     *      'alpha,1/K     ','beta,1/bar    ','S,J/K/mol     ',
c                                 16-18
     *      'n,mol         ','N,g           ','Ks_T,bar/K    ',
c                                 19-21
     *      'Gs_T,bar/K    ','Ks_P          ','Gs_P          ',
c                                 22-24
     *      'v0_T          ','vp_T          ','vs_T          ',
c                                 25-27
     *      'v0_P          ','vp_P          ','vs_P          ',
c                                 28-30
     *      'wt,%          ','vol,%         ','mol,%         ',
     *      'h,J/m3        ','cp,J/K/m3     ','blk_comp      ',
     *      'mode          ','composition   ','s,J/K/m3      ',
     *      's,J/K/kg      ','h,J/kg        ','cp,J/K/kg     ',
c                                 40-43
     *      'specific_mass ','poisson_ratio ','chemical_pot  ',
     *      'assemblage_i  ','extent        '/
c                                 l2p points from lop to prname, 
c                                 1-10
      data l2p/34,10,32,13,14,33,37,38, 3, 4,
c                                 11-20            
     *         5 , 6, 7, 8, 9,36,37,38,39,40,
c                                 21-30 
     *         41, 1,42,43, 0,22,23,24,18,19,
c                                 31-38
     *         25,26,27,20,21, 0, 44, 0/
c----------------------------------------------------------------------
c                                 make property name
      if (lop.eq.6) then
c                                 wt% component icx
         write (dname(jprop),'(a,a)') cname(icx),',wt%     '
         call unblnk (dname(jprop))

      else if (lop.eq.7) then
c                                 mode of a phase
         if (iopt(3).eq.0) then 
c                                 vol%
            write (dname(jprop),'(a,a)') pname,',vo%'

         else if (iopt(3).eq.1) then 
c                                 wt%
            write (dname(jprop),'(a,a)') pname,',wt%'

         else  
c                                 mol%
            write (dname(jprop),'(a,a)') pname,',mo%'

         end if 

         call unblnk(dname(jprop))  

      else if (lop.eq.8) then 
c                                phase composition
          write (dname(jprop),'(a,i1,a,a)') 'C',komp,pname
          write (*,1000) dname(jprop)

      else if (lop.eq.23) then 
c                                chemical potential of a component
         write (dname(jprop),'(a,a,a)') 'mu_',cname(icx),',J/mol'
         call unblnk (dname(jprop))

      else if (lop.eq.36) then 
c                                allprp option, lop points directly
c                                to prmame
         dname(jprop) = prname(jprop)

      else if (lop.eq.37) then
c                                extent of a phase
         if (iopt(3).eq.0) then 
c                                volume
            write (dname(jprop),'(a,a)') pname,',m3 '

         else if (iopt(3).eq.1) then 
c                                 mass
            write (dname(jprop),'(a,a)') pname,',kg '

         else  
c                                 mol
            write (dname(jprop),'(a,a)') pname,',mol'

         end if 

         call unblnk(dname(jprop))  

      else if (lop.eq.38) then 
c                                 custom prop (phemgp)
         dname(jprop) = prname(l2p(jprop))

      else 

         dname(jprop) = prname(l2p(lop))
          
      end if    

1000  format (/,'This composition will be designated: ',a,/)

      end  

      subroutine finprp (dim,n5name,n6name,node)
c----------------------------------------------------------------
c wrap up property output, writes blurb on processing.
c----------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer dim, i 

      logical node

      character*100 n5name,n6name

      integer inv
      character dname*14, titl1*162
      common/ cst76 /inv(i11),dname(i11),titl1

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      character*14 tname
      integer kop,kcx,k2c,iprop
      logical kfl
      double precision prop,prmx,prmn
      common/ cst77 /prop(i11),prmx(i11),prmn(i11),kop(i11),kcx(i11),
     *               k2c(i11),iprop,kfl(i11),tname  
c----------------------------------------------------------------------
c                                 write data ranges
      write (*,1090) nopt(7)
      write (*,'(5x,200(g14.7,1x))') (dname(i),i=1,iprop)
      write (*,'(a3,2x,200(g14.7,1x))') 'min',(prmn(i),i=1,iprop)
      write (*,'(a3,2x,200(g14.7,1x))') 'max',(prmx(i),i=1,iprop)      

      if (kop(1).eq.25) then 
c                                 create plt format file for "all
c                                 modes" output option
         call outmod (dim,n6name,node)

         if (dim.eq.1) then 
            write (*,1000) n6name, n5name
            write (*,1030) 
            write (*,1010) dim,'tab'
            write (*,1020) 
         else 
            write (*,1040) dim,'tab',n5name
            write (*,1010) dim,'tab'
            write (*,1060)
         end if

      else if (kcx(1).eq.999) then 

         write (*,1040) dim,'phm',n5name
         write (*,1010) dim,'phm'
         if (dim.eq.1) then 
            write (*,1070)
         else 
            write (*,1080)
         end if

      else 

         write (*,1040) dim,'tab',n5name
         write (*,1010) dim,'tab'
        
         if (dim.eq.1) then 
            write (*,1020)
         else 
            write (*,1060) 
         end if 

      end if 

      close (n5)

1000  format (/,'Output has been written to two files:',//,
     *       5x,'plt format is in file: ',a,/,
     *       5x,'1d tab format is in file: ',a)
1010  format (/,i1,'d ',a,' format files can be processed with:',/)
1020  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PERPLE_X_PLOT - a Matlab plotting script',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on tab format refer to:',/,5x,
     *       'www.perplex.ethz.ch/faq/perple_x_tab_file_format.txt',/)
1030  format (/,'plt format files can be plotted with:',//,
     *       5x,'PSVDRAW')
1040  format (/,'Output has been written to the ',i1,
     *          'd ',a,' format file: ',a)
1060  format (5x,'PSTABLE - a Perple_X plotting program',
     *      /,5x,'PERPLE_X_PLOT - a MATLAB plotting script',
     *      /,5x,'PYWERAMI - petrol.natur.cuni.cz/~ondro/pywerami:home',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on tab format refer to:',
     *      /,5x,'perplex.ethz.ch/faq/perple_x_tab_file_format.txt',/)
1070  format (5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on phm format refer to:',
     *      /,5x,'perplex.ethz.ch/faq/perple_x_phm_file_format.txt',/)
1080  format (5x,'PHEMGP - perplex.ethz.ch/phemgp',
     *      /,5x,'spread-sheet programs, e.g., EXCEL',//,
     *       'for details on phm format refer to:',
     *      /,5x,'perplex.ethz.ch/faq/perple_x_phm_file_format.txt',/)
1090  format (/,'Data ranges excluding values equal to bad_number ',
     *       '(',g10.3,') specified in perplex_option.dat:',/)
      end  