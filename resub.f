      subroutine resub (jd,id,ids,iref,iter)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompound compositions around the
c pseudocompound id of solution ids in iteration iter. ifst is the 
c pointer to the first pseudocompound of the solution ids and ilst
c points to the last. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      logical bad

      double precision xxnc, ysum

      integer i, j, k, ids, id, jd, iter, kcoct, iref
c                                 -------------------------------------
c                                 functions
      double precision gsol1, ydinc
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive g and compositions
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 temporary subdivision limits:
      double precision wg,xmn,xmx,xnc
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,
     *        jstot,kstot
      common/ cst108 /wg(m1,m3),xmn(mst,msp),xmx(mst,msp),xnc(mst,msp),
     *      iend(m4),isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),
     *      isp(mst),isite,iterm,iord,istot,jstot,kstot
c                                 coordinates output by subdiv
      double precision xy,yy
      integer ntot,npairs
      common/ cst86 /xy(mdim,k1),yy(ms1,mst,k1),ntot,npairs
c                                 max number of refinements for solution h9
      integer ncoor
      common/ cxt24 /ncoor(h9)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 option values
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
      
      if (iter.eq.1) then
c                                first iteration id array points to 
c                                original compound arrays:
         call getolx (ids,id)

      else
c                                on subsequent iterations get the y's
c                                stored in the ycoor array by routine 
c                                saver, these are reindexed copies of the
c                                coordinates originally saved in zcoor
c                                below. 
         call getxy0 (ids,id)

      end if
c                                load the subdivision limits into
c                                temporary limit arrays:
      isite = istg(ids)
      
      do i = 1, isite

         isp(i) = ispg(ids,i)

         do j = 1, isp(i) - 1

            imd(j,i) = imdg(j,i,ids)

            xxnc = nopt(14)*xncg(ids,i,j)/iopt(11)**(iter-1)

            if (imd(j,i).eq.0) then 
c                                 cartesian
               xmn(i,j) = x(i,j) - xxnc
               xmx(i,j) = x(i,j) + xxnc

            else
c                                 conformal
               xmn(i,j) = ydinc (x(i,j),-xxnc,imd(j,i),j,i,ids)
               xmx(i,j) = ydinc (x(i,j),xxnc,imd(j,i),j,i,ids)

            end if 

            xnc(i,j) = 2d0*nopt(14)*xncg(ids,i,j)/iopt(11)**iter
            if (xmn(i,j).lt.xmng(ids,i,j)) xmn(i,j) = xmng(ids,i,j)
            if (xmx(i,j).gt.xmxg(ids,i,j)) xmx(i,j) = xmxg(ids,i,j)

         end do 
      end do 
                            
      call subdv1 ('characters',ids) 

      do 10 i = 1, ntot 

         jphct = jphct + 1
         if (jphct.gt.k21) call error (58,x(1,1),k21,'resub')
c                                 convert to compositional corrdinates 
c                                 required by routine gsol, y coordinates
c                                 are placed in first array of cxt7,
c                                 store a copy of x coordinates in 
c                                 1-d array zcoor
         jkp(jphct) = ids
         hkp(jphct) = jd
         jcoor(jphct) = jcoct - 1
         kcoct = jcoct + ncoor(ids)
c                                 counter for number of non 0 or 1 compositions

         if (kcoct.gt.k20) call error (59,x(1,1),k20,'resub')

         do j = 1, isite
            ysum = 0d0
            do k = 1, isp(j) - 1
               x(j,k) = yy(k,j,i)
               ysum = ysum + x(j,k)
               zcoor(jcoct) = x(j,k)

               if (x(j,k).lt.xmno(ids,j,k).and.
     *             x(j,k).gt.xmxo(ids,j,k)) then 
c                                 the composition is out of range
                  jphct = jphct - 1
                  jcoct = kcoct - ncoor(ids)
                  goto 10 
               end if 
               jcoct = jcoct + 1
            end do 
            x(j,isp(j)) = 1d0 - ysum
            zcoor(jcoct) = x(j,isp(j))
            jcoct = jcoct + 1
         end do 

         call xtoy (ids)
c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
         g2(jphct) = gsol1(ids)
c                                 this is to check for invalid site fractions
c                                 arising from h&p's equipartition models.
c                                 this is redundant with the z calculations 
c                                 made in gsol1, but to avoid modifying the code
c                                 the test is repeated here, the test will be 
c                                 obsolete once equipartition is phased out.
         if (lopt(5).and.lrecip(ids)) then 
c                                 check for invalid site fractions, this is only necessary
c                                 for H&P models that assume equipartition (which is not 
c                                 implemented). 
            call zchk (pa,ids,bad)

            if (bad) then
               jphct = jphct - 1
               jcoct = jcoct - ncoor(ids)
               cycle
            end if 
         end if 
c                                 use the coordinates to compute the composition 
c                                 of the solution
         call csol (ids)

         iref = iref + 1
     
10    continue  

      end 

      subroutine csol (id)
c-----------------------------------------------------------------------
c csol computes chemical composition of solution id from the macroscopic
c endmember fraction array y or p0a (cxt7), these arrays are prepared by a prior
c call to function gsol. the composition is loaded into the array cp2 at
c position jphct.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,id

      double precision ctot2
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      integer jend
      common/ cxt23 /jend(h9,k12)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------

      ctot2 = 0d0

      do i = 1, icp
         cp2(i,jphct) = 0d0
      end do  

      if (lrecip(id).or.lorder(id)) then 
c                                 solutions with dependent endmembers, p0a 
c                                 contains the p's. for ksmod=8 these are a 
c                                 reformulation of the p's to eliminate the ordered 
c                                 endmembers. p0a is constructed in function gsol.
         do i = 1, lstot(id) 
            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + p0a(i) * cp(j,jend(id,2+i))
            end do 
            ctot2 = ctot2 + p0a(i)*ctot(jend(id,2+i))
         end do 

      else 
c                                 general case (y coordinates)
         do i = 1, mstot(id)

            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + y(i) * cp(j,jend(id,2+i))
            end do

            ctot2 = ctot2 + y(i)*ctot(jend(id,2+i)) 

         end do 
         
      end if 
c                                  normalize the composition and free energy
      g2(jphct) = g2(jphct)/ctot2

      do j = 1, icp 
         cp2(j,jphct) = cp2(j,jphct)/ctot2
      end do  

      end

      subroutine reopt (idead,jdv,npt)
c-----------------------------------------------------------------------
c reopt - given the results of an initial optimization for lpopt, reopt
c iteratively refines the solution by generating pseudocompounds in the
c neighborhood of the initial optimization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw, lw, iter, iref, i, j, npt, id, idead, ids, jstart, 
     *        inc

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)  

      double precision  ax(k5), x(k21), clamda(k21+k5), w(lw)

      integer is(k21+k5), iw(liw), jdv(k19)
c                                 -------------------------------------
c                                 global variables
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision xa,b,xc
      common/ cst313 /xa(k5,k1),b(k5),xc(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      double precision g
      common/ cst2  /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      iter = 1
      jphct = 0
      iref = 0 
      jcoct = 1
      inc = istct - 1
c                                 --------------------------------------
c                                 first iteration
      do i = 1, npt

         id = jdv(i) + inc

         if (id.le.ipoint) then
c                                 the point is a true compound
            jphct = jphct + 1
c                                 jkp indicates which phase a point is associated with
            jkp(jphct) = -id
c                                 hkp indicates which refinement point
            hkp(jphct) = i

            g2(jphct) = g(id)/ctot(id)

            do j = 1, icp
               cp2(j,jphct) = cp(j,id)/ctot(id)
            end do 

         else 
c                                 the point is a pseudocompound 
c                                 to be refined
            call resub (i,id,ikp(id),iref,iter)

         end if
c                                 reset jdv to point to the new
c                                 indexing, this is only for the
c                                 case of no iteration (i.e., resub
c                                 increments jphct by one.
         jdv(i) = jphct

      end do 

      if (iref.eq.0) goto 90
 
      do 
c                                 iter is incremented before the operations,
c                                 i.e., on the nth iteration, iter is n+1
         iter = iter + 1
c                                 cold start
         jstart = 0 
c                                 set idead = 0 to prevent lpnag from
c                                 overwriting warm start parameters
         idead = 0 
c                                 do the optimization
         call lpnag (jphct,icp,cp2,k5,b,g2,is,x,ax,
     *               clamda,iw,liw,w,lw,idead,l6,jstart)
c                                 warn if severe error
         if (idead.gt.0) then 
            call lpwarn (idead,'REOPT ')
            goto 99
         end if 
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,jdv,npt,iter,idead)

         if (idead.gt.0) goto 99
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the xcoor array.
         call saver (jdv,npt)

         if (iter.gt.iopt(10)) exit 

         jphct = 0 
         iref = 0 
         jcoct = 1
c                                 generate new pseudocompounds
         do i = 1, npt

            ids = lkp(i)

            if (ids.lt.0) then 
c                                 the point is a true compound
               jphct = jphct + 1
               jkp(jphct) = ids
               hkp(jphct) = mkp(i)
               ids = -ids
               g2(jphct) = g(ids)/ctot(ids)

               do j = 1, icp
                  cp2(j,jphct) = cp(j,ids)/ctot(ids)
               end do 

            else
c                                 the point is a pseudocompound 
               call resub (mkp(i),i,ids,iref,iter)

            end if
c                                 reset jdv in case we're going to exit
            jdv(i) = i 

         end do 

         if (iref.eq.0) then
c                                 has this ever happened? 
            do i = 1, npt
               jdv(i) = i 
            end do  
            exit

         end if 

      end do 

90    if (npt.ne.icp) then 
         idead = 1
      else 
c                                 get chemical potentials for output.
         call getmus (cp2,g2,jdv,idead)       
      end if 

99    end 

      subroutine sortin (ind,k,n)
c-----------------------------------------------------------------------
c sort the first k values of ind
c-----------------------------------------------------------------------
      implicit none

      integer i, j, k, n, ind(n), imin

      do j = 1, k-1

         imin = ind(j)

         do i = j+1, k

            if (ind(i).lt.imin) then 
               imin = ind(i)
               ind(i) = ind(j)
               ind(j) = imin
            end if
 
         end do 

      end do 

      end 

      subroutine lpwarn (idead,char)
c----------------------------------------------------------------------
c write warning messages from lpnag as called by routine 'char',
c set flags ier and idead, the optimization is a total fail if
c idead set to 1.
c----------------------------------------------------------------------
      implicit none

      integer idead

      character*6 char     

      double precision c
c----------------------------------------------------------------------
c                                             look for errors                                            
      if (idead.eq.2.or.idead.gt.4) then 
c                                             unbounded solution, or
c                                             other programming error.
         call warn (91,c,idead,char) 

      else if (idead.eq.3) then 
c                                             no feasible solution
         call warn (42,c,idead,char)

      else if (idead.eq.4) then 
c                                             iteration count exceeded,
c                                             probable cause no feasible
c                                             solution.
         call warn (90,c,idead,char) 

      end if

      end 

      subroutine saver (jdv,npt)
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound x(i,j) compositions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      integer jdv(k19), npt, i, j, k, kcoct, id, ids, itic
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c----------------------------------------------------------------------
      kcoct = 0

      do i = 1, npt

         id = jdv(i)
         ids = jkp(id)
         lkp(i) = ids
c                                 cycle on a compound
         if (ids.lt.0) cycle
c                                 it's a solution:
         lcoor(i) = kcoct
         itic = 0

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
               itic = itic + 1
               if (kcoct+itic.gt.k22) 
     *             call error (60,ycoor(1),k22,'saver')
               ycoor(lcoor(i)+itic) = zcoor(jcoor(id)+itic)
            end do 
         end do 

         kcoct = kcoct + itic

      end do 

      end 

      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the zcoor array loaded in resub.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jd, id, ids, icoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                  xcoordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------
      icoor = jcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            icoor = icoor + 1
            x(i,j) = zcoor(icoor)
            x3(jd,i,j) = x(i,j)
         end do 
      end do 

      end 

      subroutine getxy0 (ids,id)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the ycoor array loaded in saver
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c----------------------------------------------------------------------
      jcoor = lcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            jcoor = jcoor + 1
            x(i,j) = ycoor(jcoor)
         end do 
      end do 

      end 

      subroutine rebulk (jdv,idead)
c----------------------------------------------------------------------
c rebulk computes the amounts of the stable compounds and eliminates 
c those with zero modes (if requested).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer jdv(k19),idead,i,j,k,npt

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
c                                 load the transpose of the
c                                 concentration matrix of the pseudo-
c                                 invariant assemblage.
      npt = icp 

      do i = 1, jbulk
         if (i.le.icp) then
c                                 getcmp assigns kkp, so jkp is uneccesary
c                                 from here on
            call getcmp (i,jdv(i),jkp(jdv(i)))
            do j = 1, jbulk
               a(j,i) = cp3(j,i)
            end do
         else
            npt = npt + 1
            k = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(k).eq.0) then 
               kkp(i) = -k
            else 
               kkp(i) = ikp(k)
            end if 

            do j = 1, jbulk
               a(j,i) = cp(j,k)
            end do
         end if
      end do  
c                                 factor the matrix
      call factr1 (jbulk,idead)
      if (idead.eq.1) goto 99 
c                                 test for bounded compositions
      do i = 1, jbulk
c                                 load composition vector, the alpha
c                                 vector is returned in the same array
         b(i) = cblk(i)
      end do
c                                 solve for the alpha vector
      call subst1 (jbulk)
c                                 check for and eliminate zero mode
c                                 phases
      ntot = jbulk

      if (nopt(9).gt.0d0) call zmode (idead)

99    end

      logical function solvs1 (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of solution
c ids, called only for final solution vales by avrger.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,h8),soltol

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c-----------------------------------------------------------------------
      solvs1 = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle 

         if (dabs(cp3(i,id1) - cp3(i,id2))/dcp(i,ids).gt.soltol) then 
            solvs1 = .true.
            exit 
         end if 
      end do 

      end 

      logical function solvs2 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of solution
c ids, intermediate solution vales. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id1,id2,i

      integer icomp,iphct,icp,istct
      common/ cst6 /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,h8),soltol

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c-----------------------------------------------------------------------
      solvs2 = .false.

      do i = 1, icp

         if (dabs(cp2(i,id1) - cp2(i,id2)).gt.soltol) then
            solvs2 = .true.
            exit
         end if  

      end do 

      end 

      subroutine avrger 
c----------------------------------------------------------------------
c avrger combines discretization points into a single solution
c composition. on output

c     np  - is the number of solutions, 
c     ncpd - is the number of true compounds
c     ntot - np+ncpd

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvs1, check
c                                 -------------------------------------
c                                 local variables
      integer idsol(k5),kdsol(k5,k5),ids,isite,xidsol,xkdsol,irep,
     *        i,j,jdsol(k5,k5),jd,k,l,nkp(k5),xjdsol(k5)

      double precision bsol(k5,k5),cpnew(k5,k5),xx,ct,xb(k5), 
     *                 bnew(k5),xnew(k21,mst,msp)
c                                 -------------------------------------
c                                 global variables:
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,h8),soltol

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv
c-----------------------------------------------------------------------
c                                first check if solution endmembers are
c                                among the stable compounds:
      do i = 1, ntot
         if (kkp(i).lt.0) then 
            if (ikp(-kkp(i)).ne.0) then 
c                                we have an endmember
               nkp(i) = ikp(-kkp(i))
            else
               nkp(i) = kkp(i)
            end if 
         else 
            nkp(i) = kkp(i)
         end if 
      end do 
c                                check if any solutions
      do i = 1, ntot
         if (nkp(i).gt.0) goto 10
      end do 

      if (usv) then 
         do i = 1, ntot
            do j = i+1, ntot
               if (nkp(i).eq.nkp(j)) goto 10 
            end do 
         end do
      end if 

      np = 0
      ncpd = ntot
c                                the assemblage consists entirely 
c                                of true compounds, nothing to do, but load
c                                the cp3 array
      do i = 1, ncpd
         l = -kkp(i)
         do j = 1, icomp
            cp3(j,i) = cp(j,l)
         end do 
      end do 

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0
c                                set solvus tolerance, avrger
      soltol = 1.5d0*nopt(8)

      do 30 i = 1, ntot
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = b(i)
            kdsol(ntot,ncpd) = nkp(i)       
            jdsol(ntot,ncpd) = i   
         else 
            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (kdsol(j,1).eq.nkp(i)) then 
c                                 if match check for a solvus
                  if (.not.solvs1(i,jdsol(j,idsol(j)),nkp(i))) then
c                                 the pseudocompound matches a solution
c                                 found earlier.
                     idsol(j) = idsol(j) + 1
                     bsol(j,idsol(j)) = b(i)
                     jdsol(j,idsol(j)) = i  
                     goto 30 
                  end if 
               end if 
            end do
c                                 the pseudocompound is a new solution 
c                                 phase.
            np = np + 1
            idsol(np) = 1
            kdsol(np,1) = nkp(i)
            jdsol(np,1) = i 
            bsol(np,1) = b(i)

         end if    
30    continue  
c                                 check if a solution occurs more than once
c                                 but the occurences are not sequential (this
c                                 can only occur if an endmember is immiscible 
c                                 with a general composition
      if (np.gt.2) then
 
         do i = 1, np

            check = .false.
            irep = 0

            do j = i+1, np
               if (kdsol(j,1).ne.kdsol(i,1)) then
                  check = .true.
               else 
                  irep = irep + 1
               end if 
            end do 

            if (check.and.irep.gt.0) then 

               l = i + 1

               if (kdsol(l,1).ne.kdsol(i,1)) then 
c                                 not in sequence, find the next occurence
                  do j = i+2, np 
                     if (kdsol(i,1).eq.kdsol(j,1)) exit
                  end do 
c                                 swap phase at i+1 with the one at j
                  xidsol = idsol(l)
                  xkdsol = kdsol(l,1)
                  do k = 1, xidsol
                     xb(k) = bsol(l,k)
                     xjdsol(k) = jdsol(l,k)
                  end do 

                  idsol(l) = idsol(j)
                  kdsol(l,1) = kdsol(j,1)
                  do k = 1, idsol(j)
                     bsol(l,k) = bsol(j,k)
                     jdsol(l,k) = jdsol(j,k)
                  end do 

                  idsol(j) = xidsol
                  kdsol(j,1) = xkdsol
                  do k = 1, xidsol
                     bsol(j,k) = xb(k)
                     jdsol(j,k) = xjdsol(k)
                  end do 

               end if 
            end if 
         end do 
      end if 
c                                 if a solution is represented by
c                                 more than one pseudocompound get
c                                 the everage composition
      do i = 1, np 
c                                 initialize
         bnew(i) = 0d0

         do j = 1, icomp
            cpnew(j,i) = 0d0
         end do 

         ids = kdsol(i,1)
         isite = istg(ids)

         do j = 1, isite
            do k = 1, ispg(ids,j)
               xnew(i,j,k) = 0d0
            end do 
         end do 

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + b(jdsol(i,j))
         end do 

         do j = 1, idsol(i)

            jd = jdsol(i,j)
c                                conditional in case zero mode
c                                is off:
            if (bnew(i).gt.0d0) then 

               xx =  b(jd)/bnew(i)
c                                save the new compositions
               do k = 1, icomp
                  cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = xnew(i,k,l) + xx*x3(jd,k,l)
                  end do 
               end do 
            
            else 
c                               
               do k = 1, icomp
                  cpnew(k,i) = cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = x3(jd,k,l)
                  end do 
               end do 

            end if 

         end do 

      end do
c                                now reform the arrays kdv and b
      do i = 1, np

         b(i) = bnew(i)
         kkp(i) = kdsol(i,1)
         ids = kkp(i)

         ct = 0d0 

         do j = 1, icomp
            ct = ct + cpnew(j,i)
            cp3(j,i) = cpnew(j,i)
         end do

         ctot3(i) = ct

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
c                                 set x's for sollim
               x(j,k) = xnew(i,j,k)
c                                 set x's for global storage
               x3(i,j,k) = x(j,k) 

            end do 
         end do 
c                                 check composition against solution model ranges
c                                 if auto_refine is on:
         call sollim (ids)

      end do

      if (.not.usv) then 

         do i = 1, ncpd
            k = np + i
            l = kdsol(ntot,i)
            b(k) = bsol(ntot,i)
            kkp(k) = l
c                                for the sake of completeness load
c                                compound composition into cp3 array
            do j = 1, icomp
               cp3(j,k) = cp(j,-l)
            end do 
         end do
 
      else 

         irep = 1
         idsol(1) = kdsol(ntot,1)

         do i = 2, ncpd

            check = .false.

            do j = 1, irep

               if (kdsol(ntot,i).eq.idsol(j)) then 
                  check = .true.
                  exit 
               end if 
               
               if (check) exit

            end do 

            if (.not.check) then 
               irep = irep + 1
               idsol(irep) = kdsol(ntot,i)
            end if 
             
         end do 
         
         do i = 1, irep

            k = np + i
            l = idsol(i)

            do j = 1, ncpd
               if (kdsol(ntot,j).eq.l) b(k) = b(k) + bsol(ntot,j)
            end do 

            b(k) = bsol(ntot,i)
            kkp(k) = l
c                                for the sake of completeness load
c                                compound composition into cp3 array
            do j = 1, icomp
               cp3(j,k) = cp(j,-l)
            end do 
         end do 

         ncpd = irep

      end if 

      ntot = np + ncpd

99    end 

      subroutine sollim (ids)
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids,i,j

      double precision z, pa, p0a, x, w, y
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)

      character fname*10
      common/ csta7 /fname(h9)
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 set stable flag
      stable(ids) = .true.
c                                 check x-ranges
      do i = 1, istg(ids)
         do j = 1, ispg(ids,i) - 1
c                                 low limit:
            if (x(i,j).lt.xlo(j,i,ids)) then

               xlo(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).gt.xmno(ids,i,j).and.
     *             x(i,j).le.xmng(ids,i,j)) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmng(ids,i,j) = xmng(ids,i,j) - nopt(10)
                     if (xmng(ids,i,j).lt.0d0) xmng(ids,i,j) = 0d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then 
c                                 assymmetric stretching towards xmin
                     yint(1,j,i,ids) = yint(1,j,i,ids) - nopt(10)
                     if (yint(1,j,i,ids).lt.0d0) yint(1,j,i,ids) = 0d0
                     xmng(ids,i,j) =  yint(1,j,i,ids)

                  else 
c                                 symmetric modes, don't reset, but
c                                 set xmn to prevent future warnings
                     xmng(ids,i,j) = 0d0 
                     relax(ids) = .false.

                  end if 

                  limit(ids) = .true.

               end if 
            end if 
c                                 high limit:
            if (x(i,j).gt.xhi(j,i,ids)) then
               xhi(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).lt.xmxo(ids,i,j).and.
     *             x(i,j).ge.xmxg(ids,i,j)) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmxg(ids,i,j) = xmxg(ids,i,j) + nopt(10)
                     if (xmxg(ids,i,j).gt.1d0) xmxg(ids,i,j) = 1d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then 
c                                 assymmetric stretching
                     yint(2,j,i,ids) = yint(2,j,i,ids) + nopt(10)
                     if (yint(2,j,i,ids).gt.1d0) yint(2,j,i,ids) = 1d0
                     xmxg(ids,i,j) = yint(2,j,i,ids)

                  else 
c                                 symmetric modes, don't reset
c                                 set xmx to prevent future warnings
                     xmxg(ids,i,j) = 1d0 
                     relax(ids) = .false.

                  end if 

                  limit(ids) = .true.

               end if 
            end if 
         end do 
      end do  

      end 

      subroutine outlim 
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer i,j,k,ibad1,ibad2,ibad3,igood

      logical bad1,bad2,good
c                                 -------------------------------------
c                                 global variables:
c                                 working arrays
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character*10 fname
      common/ csta7 /fname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      logical refine
      common/ cxt26 /refine

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
      ibad1 = 0 
      ibad2 = 0 
      ibad3 = 0 
      igood = 0 
      rewind (n10)
      if (lopt(11)) rewind (n11)

      if (isoct.eq.0) goto 99

      bad1 = .false.
      bad2 = .false.
      good = .false.

      do i = 1, isoct

         if (.not.stable(i)) then
            bad1 = .true.
            ibad1 = ibad1 + 1
         else
            good = .true.
            igood = igood + 1
         end if
 
         if (limit(i)) then 
            bad2 = .true.
            if (relax(i)) then 
               ibad2 = ibad2 + 1
            else 
               ibad3 = ibad3 + 1
            end if 
         end if 

      end do 

      if (.not.refine) write (n10,*) ibad1,0,igood
c                                 write solutions present that are 
c                                 not stable
      if (bad1) then 

         write (*,1000)
         if (lopt(11)) write (n11,1000)
     
         do i = 1, isoct
            if (.not.stable(i)) then 
               write (*,'(5x,a)') fname(i) 
               if (.not.refine) write (n10,'(a)') fname(i)
               if (lopt(11)) write (n11,'(5x,a)') fname(i) 
            end if 
         end do
      end if 

      if (.not.good) goto 99
c                                 write solutions that are on an internal
c                                 limit
      if (bad2.and.icopt.gt.3) then 
c                                 adaptive minimization
         if (ibad2.gt.0) then 
c                                 solutions whose limits could be relaxed
            write (*,1080) 
            if (lopt(11)) write (n11,1010) 
            do i = 1, isoct
               if (limit(i).and.relax(i)) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

         if (ibad3.gt.0) then 
c                                 solutions whose limits could NOT be relaxed
            write (*,1090) 
            if (lopt(11)) write (n11,1010) 
            do i = 1, isoct
               if (limit(i).and.(.not.relax(i))) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

      else if (bad2) then 
c                                 non-adaptive minimization,
c                                 solutions on internal limits 
         write (*,1010) 
         if (lopt(11)) write (n11,1010) 
         do i = 1, isoct
            if (limit(i)) then
               write (*,'(5x,a)') fname(i) 
               if (lopt(11)) write (n11,'(5x,a)') fname(i)
            end if  
         end do

      end if 

      do i = 1, isoct
         if (.not.stable(i)) cycle

         if (.not.refine) then

            write (n10,'(a)') fname(i)

            do j = 1, istg(i)
               do k = 1, ispg(i,j)-1
                  write (n10,*) xlo(k,j,i),xhi(k,j,i)
               end do
            end do 

         end if 

         if (istg(i).eq.1) then 
c                                 single site solution
            write (*,1020) fname(i)
            if (lopt(11)) write (n11,1020) fname(i)
            do j = 1, ispg(i,1) - 1
               write (*,1030) names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
               if (lopt(11)) write (n11,1030) 
     *                        names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
            end do 
         else
c                                 reciprocal solution
            write (*,1040) fname(i)
            if (lopt(11)) write (n11,1040) fname(i)

            do j = 1, istg(i)
               
               write (*,1050) j
               if (lopt(11)) write (n11,1050) j

               if (ispg(i,j).eq.1) then 
                  write (*,1060)
                  if (lopt(11)) write (n11,1060) 
               else
                  do k = 1, ispg(i,j) - 1
                     write (*,1070) k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                     if (lopt(11)) write (n11,1070) 
     *                              k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                  end do 
               end if 
            end do
         end if 
      end do 

99    close (n10)
      if (lopt(11)) close (n11)

1000  format (/,'WARNING: The following solutions were input, but are',
     *          ' not stable:',/)
1010  format (/,'WARNING: The following solutions have compositions on',
     *          ' an internal limit (i.e., 0<x<1)',/,'(see ranges ',
     *          'below to determine which limits should be relaxed or',
     *        /,'if executing in auto_refine mode inrease auto_refine',
     *          '_slop in perplex_option.dat):',/)
1020  format (/,'Endmember compositional ranges for model: ',a,//,5x,
     *        'Endmember   Minimum   Maximum')
1030  format (5x,a8,4x,f7.5,3x,f7.5)
1040  format (/,'Site fraction ranges for multisite model: ',a)
1050  format (/,'  Site ',i1,/,5x,'Species   Minimum   Maximum   ')
c     *          'Endmember with this species')
1060  format (8x,'Dummy site generated by model reformulation',/)
1070  format (8x,i1,6x,f7.5,3x,f7.5,3x,12(a8,1x))
1080  format (/,'WARNING: The compositions of the following solutions ',
     *        'reached internal limits',/,'that were automatically ',
     *        'relaxed.',/)
1090  format (/,'WARNING: The compositions of the following solutions ',
     *        'reached internal limits',/,'that could not be ',
     *        'automatically relaxed. To avoid this problem change',/,
     *        'the subdivision mode for ',
     *        'the solutions or increase auto_refine_slop.',/)

      end 

      subroutine sorter (kdbulk,ico,jco,output)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,m,kdbulk,ico,jco,ids,ioct,inct

      logical output 

      double precision cpt(k5,k5),xt(k5,mst,msp),bt(k5),ct3(k5)
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
c----------------------------------------------------------------------
c                                 look for a match with known assemblages
      do 110 i = 1, iasct

         if (np.ne.iavar(1,i).or.ncpd.ne.iavar(2,i)) cycle 

         do 120 j = 1, ntot
            do k = 1, ntot 
               if (idasls(k,i).eq.kkp(j)) then 
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                   inct = 0 
                   ioct = 0 
                   do l = 1, ntot
                      if (kkp(l).eq.kkp(j)) inct = inct + 1
                      if (idasls(l,i).eq.kkp(j)) ioct = ioct + 1
                   end do 

                   if (ioct.ne.inct) goto 110 

                   goto 120

               end if 
            end do
c                                 no match with compound j:
c                                 do next assemblage
            goto 110 

120      continue  

         if (ibulk.gt.k2) call error (183,0d0,k2,'SORTER')
         ibulk = ibulk + 1
         iap(ibulk) = i
         kdbulk = ibulk
c                                 reorder the result arrays of the
c                                 current occurence to match initial 
c                                 occurence:
         do j = 1, ntot

            do k = 1, ntot

               if (kkp(k).eq.idasls(j,i)) then
c                                 load temporary array
                  bt(j) = b(k)

                  if (kkp(k).gt.0) then 

                     do l = 1, icomp
                        cpt(l,j) = cp3(l,k)
                     end do

                     ct3(j) = ctot3(k)

                     do l = 1, istg(kkp(k))
                        do m = 1, ispg(kkp(k),l)
                           xt(j,l,m) = x3(k,l,m) 
                        end do 
                     end do 
                  end if 
c                                 this eliminates immiscible phases
                  kkp(k) = 0

                  exit 
 
               end if 

            end do 

         end do
c                                 reload final arrays from temporary
         do j = 1, ntot

            b(j) = bt(j)
            ids = idasls(j,i)
            kkp(j) = ids

            if (ids.gt.0) then 

               do k = 1, icomp
                  cp3(k,j) = cpt(k,j)
               end do

               ctot3(j) = ct3(j)

               do k = 1, istg(ids)
                  do l = 1, ispg(ids,k)
                     x3(j,k,l) = xt(j,k,l) 
                  end do 
               end do 
            end if 
         end do 

         goto 98 

110   continue 
c                                 the assemblage is new:
      iasct = iasct + 1
      if (iasct.gt.k3) call error (184,0d0,k3,'BLKMAT')

      do i = 1, ntot
         idasls(i,iasct) = kkp(i)
      end do

      ibulk = ibulk + 1
      if (ibulk.gt.k2) call error (183,0d0,k2,'BLKMAT')
      kdbulk = ibulk 
      iap(ibulk) = iasct 

      iavar(1,iasct) = np
      iavar(2,iasct) = ncpd
      iavar(3,iasct) = np + ncpd
c                                
98    if (output) call outbl1 (ico,jco)
     
      end 

      subroutine outbl1 (ico,jco)
c----------------------------------------------------------------------
c output data for compositions and phases of assemblage ibulk
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ico,jco,i,j,k,ids
c                                 -------------------------------------
c                                 global variables
      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp, np, ncpd, ntot
      double precision cp3, ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 molar amounts (b)
      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)
c                                 i/o
      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      double precision mu
      common/ cst330 /mu(k8)

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 
c----------------------------------------------------------------------
      if (io4.eq.1) return
c                                graphics output  
      write (n5,1000) ico,jco,iap(ibulk)
c                                phase molar amounts
      write (n5,1010) (b(i),i=1,np+ncpd)
c                                solution phase compositions
      do i = 1, np
         ids = kkp(i)
         write (n5,1010) ((x3(i,j,k),k=1,ispg(ids,j)),j=1,istg(ids))
      end do 
c                                dependent potentials
      if (jpot.ne.1) write (n5,1010) (mu(i),i=1,hcp)

1000  format (10(i8,1x))
1010  format (6(g16.8,1x))

      end 

      subroutine yclos1 (clamda,x,is,jphct,jdv,npt,idead)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration 1.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, jdv(k19), npt, i, j, k, is(k1+k5), imin,
     *        idsol(k5), kdv(h8+1), nsol, mpt, iam, id, idead, is1,

     *        left, right, inc, jdsol(k5,k5), kdsol(k5), max, idm(h8+1)

      external ffirst

      logical solvus

      double precision clamda(k1+k5), cmin, dlamda(h8+1), x(k1), 
     *                 slam(h8+1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isoct
      common/ cst79 /isoct

      double precision dcp,soltol
      common/ cst57 /dcp(k5,h8),soltol

      integer ikp
      common/ cst61 /ikp(k1)
c----------------------------------------------------------------------
      npt = 0 
      nsol = 0 
      cmin = 1d99
      imin = 0 
      inc = istct - 1
      is1 = isoct + 1 
c                                 solvus tolerance, yclos1, was 5
      soltol = 5d0*nopt(8)

      do i = 1, jphct

         if (is(i).eq.0.or.is(i).eq.2) then 
c                                 make a list of found phases:
            id = i + inc

            if (ikp(id).ne.0) then 
               do j = 1, nsol
                  if (ikp(id).eq.idsol(j)) then 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = id
                     goto 10
                  end if 
               end do
c                                 new phase, add to list
               nsol = nsol + 1
               idsol(nsol) = ikp(id)
               jdsol(1,nsol) = id
               kdsol(nsol) = 1

            end if 
c                                 new point, add to list
10          npt = npt + 1
            jdv(npt) = i

         else if (clamda(i).lt.cmin) then 
c                                 find the nearest phase           
            imin = i
            cmin = clamda(i)

         end if 

      end do

      if (iopt(10).gt.1) then
c                                 not done iterating, add the metastable
c                                 phases
         do i = 1, is1
            slam(i) = 1d99
            idm(i) = 0 
         end do 
c                                 perp 6.6.3, make a list of metastable
c                                 phases, this list includes all compounds
c                                 and the least metastable composition of
c                                 each solution.      
         do 20 i = 1, jphct

            if (is(i).ne.1) cycle 

            id = i + inc 
            iam = ikp(id)

            if (iam.ne.0) then   

               if (clamda(i).lt.slam(iam)) then
c                                the composition is more stable
c                                than the previous composition 
c                                of the solution. check if it's 
c                                one of the stable solutions                               
                  do j = 1, nsol
                     if (iam.eq.idsol(j)) then
c                                it's already stable, only accept
c                                it if its further than the solvus
c                                tolerance from any of the stable
c                                compositions.
                        do k = 1, kdsol(j)
                           if (.not.solvus(jdsol(k,j),id,iam)) goto 20
                        end do 
                     end if
                  end do
c                                the composition is stable
                  slam(iam) = clamda(i)
                  idm(iam) = i

               end if

            else 
c                                a compound, save only one
               if (clamda(i).lt.slam(is1)) then 
                  slam(is1) = clamda(i)
                  idm(is1) = i
               end if 
 
            end if 

20       continue 
c                                 load the metastable points into
c                                 kdv
         mpt = 0 

         do i = 1, is1

            if (idm(i).eq.0) cycle
            mpt = mpt + 1
            kdv(mpt) = idm(i)
            dlamda(mpt) = clamda(idm(i))

         end do 

         if (mpt.eq.0) then
c                                 no points other than active
            goto 999

         else if (mpt.le.3) then 
c                                 less metastable refinement points than
c                                 iopt(12)
            max = mpt

         else 
c                                 sort the metastable points to
c                                 find the most stable iopt(12) points
            left = 1
            right = mpt
            max = 3

            call ffirst (dlamda,kdv,left,right,max,h8+1,ffirst)

         end if 
 
         do i = 1, max
            npt = npt + 1
            jdv(npt) = kdv(i)
         end do
c                                 sort the phases, why? don't know, but it's 
c                                 necessary
         call sortin (jdv,npt,k19) 

      else 
c                                 done iterating, try to correct bad results
         if (npt.gt.icp) then 
c                                 look through for zero phases
            do i = 1, npt
               if (x(jdv(i)).eq.0d0) idead = idead + 1
            end do

            if (idead.ge.npt-icp) then 

               idead = npt - icp 
               mpt = 0
c                                 eliminate zero phase
               do i = 1, npt 

                  if (x(jdv(i)).eq.0d0) then 
                     idead = idead - 1
                     cycle
                  end if 

                  mpt = mpt + 1
                  jdv(mpt) = jdv(i)
                  if (idead.eq.0) exit

               end do 

               npt = mpt 

            else 

               idead = 1

            end if 

         else if (icp-npt.eq.1.and.imin.ne.0) then 
c                                 recover by adding the closest point 
            npt = npt + 1
            jdv(npt) = imin

         else if (npt.lt.icp) then 

            idead = 1 

         end if 

      end if 

999   end 

      subroutine ffirst (a, ind, left, right, k, n, dumsub)
c-----------------------------------------------------------------------
c find the k smallest values of array between indices left and right
c from http://en.wikipedia.org/wiki/Selection_algorithm
c-----------------------------------------------------------------------
      implicit none

      integer left, right, k, n, pivot, opivot, partit, ind(n)

      external dumsub

      double precision a(n)

      if (right.gt.left) then 

         opivot = left + (right-left)/2
         pivot = partit (a, ind, left, right, opivot, n)

         if (pivot.gt.k) then 
            call dumsub (a,ind,left,pivot-1,k,n,dumsub)
         else if (pivot.lt.k) then 
            call dumsub (a,ind,pivot+1,right,k-pivot,n,dumsub)
         end if 

      end if 

      end 

      integer function partit (a, ind, left, right, opivot, n)

      implicit none

      integer left, right, n, pivot, opivot, iold, ind(n), i

      double precision a(n), value, oldval

      value = a(opivot)
c                                 swap a(opivot) with a(right)
      iold = ind(opivot)
      a(opivot) = a(right)
      ind(opivot) = ind(right)
      a(right) = value
      ind(right) = iold

      pivot  = left

      do i = left, right-1

         if (a(i).le.value) then
 
            iold = ind(pivot)
            oldval = a(pivot)
            a(pivot) = a(i)
            ind(pivot) = ind(i)
            a(i) = oldval
            ind(i) = iold
            
            pivot = pivot + 1
         end if 
      end do 
c                                 swap a(right) with a(pivot)
      iold = ind(pivot)
      oldval = a(pivot)
      a(pivot) = a(right)
      ind(pivot) = ind(right)
      a(right) = oldval 
      ind(right) = iold
      
      partit = pivot

      end 

      subroutine subst (a,ipvt,n,b,ier)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision a(k8,k8),b(k8),x(k8),sum

      integer ipvt(k8),ip,i,j,n,ii,ier
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)
      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j)+sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous,
c                                 but reopt requires it. should check
c                                 what's with factor. 
         ier = 1
         goto 99
      end if 

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do 

         if (a(i,i).eq.0d0) then
c                                 as above.
            ier = 1
            goto 99
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 
      b(n) = x(n)
 
99    end

      subroutine subst1 (n)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x(k5), sum

      integer n, i, j, im1, ip1, nm1, ii, ip

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

c                            solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n
         sum = 0d0
         im1 = i-1
         do j = 1, im1
            sum = a(i,j)*x(j)+sum
         end do 
         ip = ipvt(i)
         x(i) = b(ip)-sum
      end do 
c                            solve ux = y for x:
      x(n) = x(n)/a(n,n)
      nm1 = n-1

      do ii = 1, nm1
         i = n-ii
         ip1 = i+1
         sum = 0d0
         do j = ip1, n
            sum = a(i,j)*x(j)+sum
         end do
         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)
      end do 

      b(n) = x(n)
 
      end

      subroutine factr1 (n,ier)
c-----------------------------------------------------------------------
c factr1 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,n,ip1,istr,ier

      double precision temp,ratio,tmax,rmax

      integer ipvt
      double precision a,d,x(k5)
      common/ cst301 /a(k5,k5),d(k5),ipvt(k5)
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1,n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         x(i) = rmax
      end do 
c                            begin decomposition:
      do i = 1, n - 1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/x(i)
         istr = i
         ip1 = i + 1

         do j = ip1, n
            tmax = dabs(a(j,i))/x(j)
            if (tmax.le.rmax) cycle
            rmax = tmax
            istr = j
         end do 

         if (dabs(rmax).lt.1d-5) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = x(istr)
            x(istr) = x(i)
            x(i) = temp
            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 
         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)
            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do 
 
      end do 
     
      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine ufluid (fo2)
c----------------------------------------------------------------------
c subroutine ufluid computes the potential of the components
c of a saturated fluid phase. if the mole fraction of a component is les
c less than 1.d-38 the chemical potential is set to -9.9d09.
c ufluid may call one of three molecular fluid equations of state, or
c alternatively users may supply their own routines, however,
c the routines currently in use return the log of a components fugacity
c which is then added to the reference state potential computed by the
c function gphase.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
 
      double precision xf(2),fo2,fs2,gph

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(2)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn
c-----------------------------------------------------------------------
c                           compute the chemical potentials of
c                           fluid components in fluid saturated
c                           systems.
      call cfluid (fo2,fs2)

      if (idfl.ne.0) then
         call gphase (idfl,gph)
         uf(idfl) = gph + r * t * f(idfl)
      else
         xf(1) = 1d0 - xco2
         xf(2) = xco2
 
         do i = 1, 2
            if (iff(i).ne.0) then 
               if (xf(i).lt.1d-38) then 
                  uf(i) = -1d10
               else 
                  call gphase (i,gph)
                  uf(i) = gph + r * t * f(i)
               end if
            end if 
         end do

      end if 
 
      end

      subroutine uproj
c----------------------------------------------------------------------
c subroutine uproj computes the potentials of saturated phase components
c and saturated components.
  
c the energies of saturated components are projected through
c saturated volatile components.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i,j,k,l,ict,ll,i1,id

      double precision uss(h6),fo2,gph,u

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      fo2 = 0d0
c                                 compute the chemical potentials of
c                                 saturated phase components.
      if (ifyn.ne.1) call ufluid (fo2)

      do i = 1, isat
c                                 determine stable saturated composants
c                                 and the corresponding chemical potentials
         ict = isct(i)

         ll = icp+i

         do j = 1, ict

            k = ids(i,j)
            call gphase (k,gph)
            
            if (ifct.gt.0) then 
               do l = 1, 2
c                                 legendre transform for saturared phase
c                                 component potentials
                  if (iff(l).ne.0) gph = gph - cp(iff(l),k)*uf(l)
               end do 
            end if 

            uss(j) = gph 

            if (i.gt.1) then 
c                                 if multiple component saturation constraints
c                                 apply saturation hierarchy legendre transform:
               i1 = i-1
               do l = 1, i1
                  uss(j) = uss(j)-cp(icp+l,k)*us(l)
               end do
            end if 

            g(k) = uss(j)
            uss(j) = uss(j)/cp(ll,k)
         end do 
c                                 if O2, check if fo2 has been 
c                                 determined by a fluid phase routine,
c                                 if so, add the transform:
         if (io2.eq.i) then 
            do j = 1, ict 
               uss(j) = uss(j) + r*t*fo2
            end do 
         end if 
c                           now find stable "composant":

         u = uss(1)

         id = 1

         if (ict.ne.1) then 
            do j = 2, ict
               if (uss(j).gt.u) cycle  
               id = j
               u = uss(j)
            end do 
         end if 
c                               save the id of the stable composant.
         idss(i) = ids(i,id)
c                               and its chemical potential.
         us(i) = u
c                               in case a phase in the component
c                               saturation space is an endmember of
c                               a solution transform the endmember G's:
         do j = 1, ict
            k = ids(i,j)
            g(k) = g(k) - cp(icp+i,k)*u
         end do 

      end do 

      end

      subroutine initlp 
c--------------------------------------------------------------------
c initialize arrays and constants for lp minimizarion
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,inc,kphct,id,jk

      logical bad

      double precision ctotal,u,s,vol,tot 

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 stuff used in lpnag 
      double precision wmach(9)
      common /ax02za/wmach

      integer ldt,ldq
      common /be04nb/ldt,ldq

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c-----------------------------------------------------------------------
      inc = istct - 1

      tloop = 400
      ploop = 400
      dv(1) = (vmax(1)-vmin(1))/(ploop-1)
      dv(2) = (vmax(2)-vmin(2))/(tloop-1)
c                                 load arrays for lp solution
      jphct = iphct - inc

      if (.not.usv) then
c                                 pressure and temperature are allowed 
c                                 EoS variables
         ctotal = 0d0

         do i = 1, icp
            ctotal = ctotal + cblk(i)
         end do 
c                                 composition constraint
         do i = 1, icp
            b(i) = cblk(i)/ctotal
         end do 

         do i = 1, jphct
            id = i + inc
            iam(i) = id
            do j = 1, icp
               a(j,i) = cp(j,id)/ctot(id)
            end do
         end do

         ldt = icp + 1
         ldq = icp + 1

      else
c                                 using usv formulation
         kphct = 0 

         do i = 1, jphct

            id = i + inc
            jk = 0
c                                 for each static chemical cpd
            do j = 1, tloop
c                                 for each temperature
               do k = 1, ploop
c                                 for each pressure
c                 
                  jk = jk + 1  

                  call pt4sv (jk)  

                  call getusv (u,s,vol,id,bad)
  
          
                  kphct = kphct + 1
                  
                  tot = ctot(id) + s + vol
                  
                  do l = 1, jbulk
                    a(l,kphct) = cp(l,id)/tot
                  end do

                  a(icp+1,kphct) = s/tot
                  a(icp+2,kphct) = vol/tot
                  c(kphct) = u/tot
                  iam(kphct) = id
                  jam(kphct) = jk

               end do 

            end do 

         end do 

         jphct = kphct

         ldt = jbulk + 1
         ldq = jbulk + 1

      end if 
c                                 cold start istart = 0
      istart = 0
c                                 stuff for lpnag
      wmach(3) = 1.11022302462516d-16
      wmach(4) = dsqrt(wmach(3))
      wmach(5) = 2.22507385850721d-308
      wmach(7) = 1d0/wmach(5)
      wmach(8) = dsqrt(wmach(7))
      wmach(9) = max(1d0/wmach(4),1d2)
      wmach(2) = wmach(3)**0.8d0
      wmach(1) = wmach(3)**0.9d0

      end 

      subroutine gall 
c-----------------------------------------------------------------------
c subroutine gall computes molar free energies of all phases.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,id

      double precision dg1,gval,dg,gzero,g0(k5),gex

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer isoct
      common/ cst79 /isoct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ifp
      common/ cxt32 /ifp(k1)

      integer ixp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1)

      integer jend
      common/ cxt23 /jend(h9,k12)
c                                 working arrays
      double precision z, pa, p0a, x, w, y
      common/ cxt7 /x(m4),y(m4),pa(m4),p0a(m4),z(mst,msp),w(m1)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c-----------------------------------------------------------------------
c                                 compute the chemical potential
c                                 of the projected components.
      call uproj
c                                 first do the endmembers:
      do id = istct, ipoint

         call gcpd(id,gval)
c                                 this is a screw up solution
c                                 necessary cause uf(1) and uf(2)
c                                 are allocated independent of ifct!
         if (ifct.gt.0) then 
            do j = 1, 2
               if (iff(j).ne.0) gval = gval - cp(iff(j),id)*uf(j)
            end do 
         end if 

         do j = 1, isat
            gval = gval - us(j) * cp(icp+j,id)
         end do

         g(id) = gval 

      end do 
c                                 now do solutions:
      do i = 1, isoct
c                                 check if normal solution:
         if (ifp(id).le.0) then 
c                                 it's normal margules or ideal:
            do j = 1, jend(i,2)
c                                 initialize with excess energy, dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               do k = 1, lstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 

               id = id + 1

            end do 

         else if (ifp(id).eq.1) then
c                                 it's a fluid compound, the way 
c                                 things are now it must have two
c                                 components.
c                                 first get the reference state 
c                                 free energy:
            
            do j = 1, mstot(i)
               g0(j) = gzero(jend(i,2+j))
            end do 

            do j = 1, jend(i,2)

               call fexces (id,gval)

               g(id) = g0(1) * sxs(ixp(id)+1) + g0(2) * sxs(ixp(id)+2) 
     *                                        + gval
               id = id + 1

            end do 

         else if (lrecip(i).and.lorder(i)) then
c                                 reciprocal solution with ordering 

c                                 compute margules coefficients
            call setw (i)
c                                 now for each compound:
            do j = 1, jend(i,2)
c                                 assign x's
               do k = 1, nstot(i) 
                  p0a(k) = sxs(ixp(id)+k)
                  pa(k) = p0a(k)
               end do 
c                                 get the speciation energy effect
               call specis (dg,i)
c                                 and endmember dqf
               call gexces (id,dg1)

               g(id) = dg + dg1 
c                                 add in g from real endmembers, this
c                                 must include the g for the disordered equivalent
c                                 of the ordered species
               do k = 1, lstot(i)

                  g(id) = g(id) + g(jend(i,2+k)) * p0a(k)

               end do 

               id = id + 1

            end do 

         else if (lorder(i)) then
c                                 compute margules coefficients
            call setw (i)
c                                 now for each compound:
            do j = 1, jend(i,2)
c                                 here excess gets the endmember entropies
c                                 and dqf corrections
               call gexces (id,g(id))

               do k = 1, lstot(i)
                  p0a(k) = sxs(ixp(id)+k)
                  pa(k) = p0a(k)
                  g(id) = g(id) + g(jend(i,2+k)) * p0a(k)
               end do 
c                                 get the speciation energy effect
               call specis (dg,i)

               g(id) = g(id) + dg  

               id = id + 1

            end do 

         else if (llaar(i)) then 
c                                 compute margules coefficients
            call setw (i)
c                                 because the hp van laar may have p-t
c                                 dependent volumes, the full expression
c                                 must be evaluated here:                   
            do j = 1, jend(i,2)
c                                 initialize with dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               do k = 1, lstot(i) 
                  x(k) =  sxs(ixp(id)+k)
                  g(id) = g(id) + g(jend(i,2+k)) * x(k)
               end do 
c                                 add the real excess energy
               g(id) = g(id) + gex(i,x)

               id = id + 1

            end do 

         else if (ifp(id).eq.23) then
c                                 toop melt:
            do j = 1, nstot(i)
               g0(j) = gzero(jend(i,2+j))
            end do 

            do j = 1, jend(i,2)
               call toop (id,g(id))
               do k = 1, nstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 
               id = id + 1
            end do 
         else if (ifp(id).eq.26) then
c                                 H2O-CO2-Salt:
            do j = 1, jend(i,2)
        
               call hcneos (g(id),sxs(ixp(id)+1),
     *                      sxs(ixp(id)+2),sxs(ixp(id)+3))

               do k = 1, nstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 

               id = id + 1

            end do
         end if 
      end do 

      end

      logical function solvus (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two static pseudocompounds of
c solution ids.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp 

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,h8),soltol

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      solvus = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle 

         if (dabs(cp(i,id1)-cp(i,id2))/dcp(i,ids).gt.soltol) then 
            solvus = .true.
            exit
         end if 

      end do 

      end 

      subroutine fopen (n2name,prt,plt,n9name,jbulk,icp,icopt,j)
c-----------------------------------------------------------------------
c open files for subroutine input1.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      logical first

      integer ierr,icopt,jbulk,icp,j
 
      character*100 blank*1,n2name,prt*3,plt*3,name,n9name

      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      character*100 prject,tfname
      common/ cst228 /prject,tfname

      save first,blank

      data first,blank/.true.,' '/
c----------------------------------------------------------------------
c                                 open thermodynamic data file
      call fopen2 (0,n2name)
      if (first) write (*,1170) n2name
c                                 open print/plot files if requested
      if (prt.ne.blank.and.prt.ne.'no_') then 
         io3 = 0 
         call mertxt (name,prject,'.prn',0)
         open (n3, file = name)
         if (first) write (*,1180) name
      else
         io3 = 1
         if (first) write (*,1180) 'none requested'
      end if

      if (plt.ne.blank.and.plt.ne.'no_') then
         io4 = 0
         call mertxt (name,prject,'.plt',0)
         open (n4, file = name)
         if (first) write (*,1190) name
      else
         io4 = 1
         if (first) write (*,1190) 'none requested'
      end if

      if (n9name.ne.blank) then
         io9 = 0 
c                                 open solution model file
         open (n9,file = n9name,iostat = ierr,status = 'old')
         if (ierr.ne.0) call error (120,0d0,n9,n9name)
         if (first) then 
            call mertxt (name,prject,'_pseudocompound_glossary.txt',0)
            write (*,1200) name
            write (*,1210) n9name
         end if 
      else
         io9 = 1
         if (first) write (*,1210) 'none requested'
       end if

      if (jbulk.ge.icp.and.io4.ne.1) then
c                                 create special plot output file
         call mertxt (name,prject,'.blk',0)
         open (n5, file = name)
         if (first) write (*,1220) name

      else if (jbulk.ge.icp.and.io4.eq.1) then 
         if (first) write (*,1220) 'none requested'
      end if
c                              if jtest = 3, write a list of reactions for
c                              stefano
      if ((icopt.eq.1.or.icopt.eq.3).and.j.eq.3) then 
         call mertxt (name,prject,'_reaction_list.txt',0)
         open (n6,file=name)
         if (first) write (*,1250) name
      end if 

      first = .false.

1170  format (/,'Reading thermodynamic data from file: ',a)
1180  format ('Writing print output to file: ',a)
1190  format ('Writing plot output to file: ',a)
1200  format ('Writing pseudocompound glossary to file: ',a)
1210  format ('Reading solution models from file: ',a)
1220  format ('Writing bulk composition plot output to file: ',a)
1250  format ('Writing complete reaction list to file: ',a)
      end 

      subroutine grxn (gval) 
c-----------------------------------------------------------------------
c grxn computes the free energy of univariant equilibria
c defined by the data in commonn block cst21 which is initialized
c in the subprogram balanc.  grxn is partially redundant with
c the function gphase but because of the frequency that these
c these routines are used a significant increase in efficiency is
c gained by maintaining separate functions.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j

      double precision gval,gproj

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct
c-----------------------------------------------------------------------
c                                 compute potentials of saturated phases
c                                 and components, note that in this
c                                 version of vertex the stoichiometry of
c                                 such components may vary.

c                                 no saturated phase components and no
c                                 saturated components:
      if (iffr.eq.1.and.isyn.eq.1) goto 10
c                                 note that this call to uproj makes a
c                                 subsequent call in gall redundant if
c                                 sfol1 is used to trace a univariant
c                                 curve.
      call uproj
c                                 compute free energy change of the rxn
10    gval = 0d0

      do j = 1, ivct
         gval = gval + vnu(j) * gproj(idr(j))
      end do 

      end

      subroutine lpopt0 (idead)
c-----------------------------------------------------------------------
c lpopt0 - calls lp minimization after a call to initlp. lpopt0
c does the minimization, writes error messages if necessary.

c this is an utterly stupid formulation of the lp problem because i modified
c the lp code to impose the implicit constraint that phase amounts were between
c 0 and 1, but did not impose the constraint the sum of the amounts must be
c be 1 (which would be unwise in any case). this requires that compositions and
c g's must be normalized to the total number of moles (and therefore lots of 
c extra bookkeeping). 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw,lw,k,idead,npt,inc

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),clamda(k1+k5),w(lw),oldt

      integer is(k1+k5),iw(liw),jdv(k19)
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision g
      common/ cst2 /g(k1)

      integer jphct,istart
      common/ cst111 /jphct,istart

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      save ax, x, clamda, w, is, iw
c-----------------------------------------------------------------------
      inc = istct - 1

      if (.not.usv) then 

         oldt = t

         if (t.lt.nopt(12)) t = nopt(12)

         call gall

         do k = 1, jphct
            c(k) = g(k+inc)/ctot(k+inc)
         end do

      end if 
c                                 idead = -1 tells lpnag to save parameters
c                                 for subsequent warm starts
      idead = -1
c                                 optimize by nag
      call lpnag (jphct,hcp,a,k5,b,c,is,x,ax,
     *            clamda,iw,liw,w,lw,idead,l6,istart)

      if (idead.gt.0) then
c                                 look for severe errors                                            
         call lpwarn (idead,'LPOPT ')
c                                 on severe error do a cold start.
c                                 necessary?
         istart = 0
         return

      end if 

      if (hcp.eq.1.or.iopt(10).eq.0.or.isoct.eq.0.or.usv) then 
c                                 no refinement, final processing
c                                 with yclos0
         call yclos0 (x,is,jphct,jdv,idead) 
c                                 get potentials, this is just for output
c                                 if not usv.
         call getmus (a,c,jdv,idead)

      else
c                                 find discretization points
c                                 for refinement
         call yclos1 (clamda,x,is,jphct,jdv,npt,idead)

         if (idead.gt.0) return 
c                                 reoptimize with refinement
         call reopt (idead,jdv,npt)
c                                 coming out of reopt the stable points
c                                 are indexed by jdv
         if (idead.gt.0) return 
c                                 get amounts of pseudocompounds
         call rebulk (jdv,idead)

      end if 
c                                 test for solvi and average
      if (idead.eq.0) call avrger       

      t = oldt

      end 

      subroutine yclos0 (x,is,jphct,jdv,idead)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jdv(k19),inc,i,j,idead,jphct,mpt,npt,is(k1+k5),id

      double precision u,x(k1)

      logical bad

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision ctot
      common/ cst3  /ctot(k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer ikp
      common/ cst61 /ikp(k1)

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c----------------------------------------------------------------------

      npt = 0 
      inc = istct - 1

      do i = 1, jphct

         if (is(i).ne.0.and.is(i).ne.2) cycle  
c                                 new point, add to list, two acceptable
c                                 cases is = 0 inactive 
c                                       is = 2 active at upper bound (this
c                                            is really only possible for a 
c                                            one component system (though
c                                            conceivably might occur in 
c                                            fractionation calculations).
            npt = npt + 1

            jdv(npt) = i 
 
      end do

      if (hcp.eq.1.and.npt.gt.1) then 

         npt = 0 

         do i = 1, jphct

            if (is(i).ne.2) cycle 

            npt = npt + 1
            jdv(npt) = i

         end do

      end if 

      if (npt.gt.hcp) then 
c                                 look through for zero phases
         do i = 1, npt
            if (x(jdv(i)).eq.0d0) idead = idead + 1
         end do

         if (idead.ge.npt-hcp) then 

            idead = npt - hcp 
            mpt = 0
c                                 eliminate zero phase
            do i = 1, npt 

               if (x(jdv(i)).eq.0d0.and.idead.gt.0) then 
                  idead = idead - 1
c                                 reset the active flag
c                 is(jdv(i)) = 1
                  cycle
               end if 

               mpt = mpt + 1
               jdv(mpt) = jdv(i)

            end do 

            npt = mpt 

         else 

            idead = 1

         end if 

      else if (npt.lt.hcp) then 

         idead = 1 

      end if 

      if (idead.eq.1) then 
        call warn (88,0d0,idead,'LPOPT') 
        goto 99
      end if 
c                                 the assemblage is ok, get the modes and 
c                                 save the data, this means creating the kkp
c                                 array, the x3 array (for output), the cp3
c                                 array (only for meemum)
      do i = 1, jbulk

         if (i.le.hcp) then 
c                                 normal phase
            id = iam(jdv(i))
         else
c                                 saturated phase
            id = idss(i-hcp)
         end if 
c                                 set identifier flag
         if (i.le.ipoint) then
            kkp(i) = -id
         else  
            kkp(i) = ikp(id)
         end if 
c                                 load the composition matrix
c                                 a and b are used by getamt
         do j = 1, jbulk

            b(j) = cblk(j)
            if (j.gt.icp.and.usv) cycle
            a(j,i) = cp(j,id)
            cp3(j,i) = cp(j,id)

         end do

         if (usv) then 

            call pt4sv (jam(jdv(i)))

            call getusv (u,a(icp1,i),a(jbulk,i),id,bad)

            if (bad) then 
               idead = 1
               goto 99
            end if 

         end if 
c                                 total component amounts
         ctot3(i) = ctot(id)
c                                 set the x3 array
         if (ikp(id).ne.0) call setx3 (i,id,ikp(id))

      end do 

      ntot = jbulk
c                                 factor the matrix
      call factr1 (jbulk,idead)

      if (idead.eq.1) goto 99 
c                                 get mode, solve for the alpha vector
      call subst1 (jbulk)
c                                 test for zero modes
      if (nopt(9).gt.0d0) call zmode (idead)

99    end 

      subroutine zmode (idead)
c----------------------------------------------------------------------
c rebulk computes the amounts of the stable compounds and eliminates 
c those with zero modes (if requested).
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,ipvt,jzero(k5), idead,tictoc

      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      save tictoc
      data tictoc/0/
c----------------------------------------------------------------------

      ntot = 0

      do i = 1, jbulk
         if (b(i).lt.-nopt(9)) then
c                                 negative mode
            if (icp.lt.jbulk.and.tictoc.lt.5) then 
c                                 probably because saturated component
               call warn (1,b(1),i,'ZMODE')
               tictoc = tictoc + 1

            else 

               write (*,*) ' a negative mode, serious problem?'
               write (*,*) (kkp(j),j=1,jbulk)
               write (*,*) (b(j),j=1,jbulk)
c                                 this probably is serious because it will
c                                 cause the back-calculated bulk to be incorrect.
               idead = 1

            end if 

         else if (b(i).gt.nopt(9)) then 

            ntot = ntot + 1
            jzero(ntot) = i 

         end if 

      end do 

      if (jbulk.ne.ntot) then

         do i = 1, ntot

            j = jzero(i)
c                                 jkp and jdv are probably unnecessary.
            kkp(i) = kkp(j)
            b(i) = b(j)

            do k = 1, icomp
               cp3(k,i) = cp3(k,j)
            end do 

            ctot3(i) = ctot3(j)

            if (kkp(i).gt.0) then 

               do k = 1, istg(kkp(i))
                  do l = 1, ispg(kkp(i),k)
                     x3(i,k,l) = x3(j,k,l)
                  end do 
               end do 
            end if 
         end do 
      end if

      end

      subroutine setx3 (ind,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the xcoor (reciprocal) or sxs (single site) arrays loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor, ind
c                                 x coordinate description
      integer istg, ispg, imlt, imdg

      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 stored x coordinate
      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)
c                                 single site solution coordinates:
      integer ixp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1)
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ipoint,imyn
      common/ cst60 /ipoint,imyn

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      if (id.gt.ipoint) then 
c                                 a normal solution
         if (istg(ids).eq.1) then 
c                                 one site solution
            do j = 1, nstot(ids)
               x3(ind,1,j) = sxs(ixp(id)+j) 
            end do 

         else if (ispg(ids,1).gt.1) then 
c                                 multi-site solution
            jcoor = icoor(id)

            do i = 1, istg(ids)
               do j = 1, ispg(ids,i)
                  jcoor = jcoor + 1
                  x3(ind,i,j) = xcoor(jcoor)
               end do 
            end do 

         else 
c                                 a dummy site
            x3(ind,1,1) = 1d0
          
            do j = 1, nstot(ids)
               x3(ind,2,j) = sxs(ixp(id)+j) 
            end do 

         end if 

      else 
c                                 an endmember 
         call endcp (ind,id,ids)

      end if

      end 

      subroutine getmus (a,c,jdv,idead)
c-----------------------------------------------------------------------
c getmus computes chemical potentials from a pseudocompound assemblage
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, idead,  ier, ipvt(k8), jerk, jdv(k19)

      double precision  a(k5,*), c(*), comp(k8,k8)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      double precision mu
      common/ cst330 /mu(k8)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jtest,jpot
      common/ debug /jtest,jpot

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      save jerk 
      data jerk/0/ 
c-----------------------------------------------------------------------
      if (.not.usv.and.jpot.eq.1) return

      if (idead.eq.0) then 
c                                 compute chemical potentials if requested:
         do i = 1, hcp
            do j = 1, hcp 
               comp(i,j) = a(j,jdv(i))
            end do 
            mu(i) = c(jdv(i))
         end do 
c                                 compute chemical potentials if requested:
c                                 reload last solution (no test for uniqueness)
         call factor (comp,hcp,ipvt,ier)

         if (ier.eq.1.and.jerk.lt.1) then 

            jerk = jerk + 1
            write (*,*) 'Deader than a doornail, but dont worry' 
     
         else 
 
            call subst (comp,ipvt,hcp,mu,ier)

            if (ier.eq.1.and.jerk.lt.2) then 
               jerk = jerk + 1
               write (*,*) 'Deader than 2 doornails, but dont worry' 
            end if 

         end if 

      end if 

      if (idead.ne.0.or.ier.eq.1) then 
c                                 optimization failed.
         do i = 1, hcp
            mu(i) = nopt(7)
         end do

      end if 

      end 

      double precision function ginc0 (dt,dp,id)
c-----------------------------------------------------------------------
c id indicates a static pseudocompound, use ginc for dynamic compound
c-----------------------------------------------------------------------
      implicit none

      double precision dt,dp,gee

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      p = p + dp 
      t = t + dt 

      call gphase (id,gee)

      p = p - dp 
      t = t - dt

      ginc0 = gee 

      end 

      subroutine getusv (u,s,v,id,bad)
c-----------------------------------------------------------------------
c getusv computes u,s,v by centered finite differences from the Gibbs energy

c the difference increments are

c dt0, dp0 for 1st order derivatives (entropy,volume and enthalpy)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      logical bad

      double precision dt0,dp0,dp1,dp2,g0,v,ginc0,s,u

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save dt0
      data dt0/0.5d0/
c----------------------------------------------------------------------

      dp0 = 0.5d-3 * p 
      dp1 = 0.5d-2 * p
      dp2 = 0.5d-1 * p
            
      g0 = ginc0(0d0,0d0,id)
c                                 straight derivatives:
c                                 first order
      if (p-dp0.le.0d0) then 

         v = (ginc0(0d0,dp0,id) - g0)/dp0
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment if invalid v
     *   v = (ginc0(0d0,dp1,id) - g0)/dp1
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment more if invalid v
     *   v = (ginc0(0d0,dp2,id) - g0)/dp2

      else 

         v = (ginc0(0d0,dp0,id) - ginc0(0d0,-dp0,id))/dp0/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp1.gt.0d0)  
c                                 expand increment if invalid v
     *   v = (ginc0(0d0,dp1,id) - ginc0(0d0,-dp1,id))/dp1/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp2.gt.0d0)  
c                                 expand increment more if invalid v
     *   v = (ginc0(0d0,dp2,id) - ginc0(0d0,-dp2,id))/dp2/2d0

      end if 
c                                 in case the evaluating routine fails
c                                 on both calls to ginc 
      if (v.eq.0d0) v = 1d0

      s = (ginc0(-dt0,0d0,id) - ginc0(dt0,0d0,id))/dt0/2d0

      u = g0 + t*s - p*v
 
      if (s.lt.0d0.or.v.lt.0d0) then 
         bad = .true.
      else
         bad = .false.
      end if 

      end 

      subroutine pt4sv (ind)
c--------------------------------------------------------------------
c subroutine to recover the p-t condition of discreitization in USV
c calculations.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind,tind,pind

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c----------------------------------------------------------------------

      tind = (ind-1)/tloop
      pind = ind - tind*tloop - 1
      v(2) = vmin(2) + tind*dv(2)
      v(1) = vmin(1) + pind*dv(1)

      end 

      subroutine yclos2 (clamda,x,is,jdv,npt,iter,idead)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jdv(k19), npt, i, is(k21+k5), idead, id, gmin,
     *        jmin(k19), mpt, jzero, iter

      logical stable

      double precision clamda(k21+k5), clam(k19), x(k21), glam

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision amt,cptot
      common/ cst78 /amt(k5),cptot(k5)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)
c----------------------------------------------------------------------
c                                 opt is the number of points refined
c                                 from the previous cycle
      mpt = npt

      do i = 1, npt
         jmin(i) = 0 
         clam(i) = 1d99
      end do 

      gmin = 0
      glam = 1d99 
      npt = 0

      do i = 1, jphct
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 check the stability of all points 
         if (is(i).eq.0.or.is(i).eq.2) then 
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            stable = .true.

         else if (clamda(i).lt.clam(id)) then 
c                                 find the nearest phase           
            jmin(id) = i
            clam(id) = clamda(i)
            if (clam(id).gt.glam) cycle
            glam = clam(id)
            gmin = i

         end if 

      end do 

      if (iter.le.iopt(10)) then
c                                 if not done iterating, add the metastable
c                                 phases
         do i = 1, mpt
            if (jmin(i).eq.0) cycle
            npt = npt + 1
            if (npt.gt.k19) then 
c                                 this should never happen
               write (*,*) 'WARNING! this should never happen, too many'
               write (*,*) 'refinement points in YCLOS2, if you see'
               write (*,*) 'this msg please tell me, if you are not me'
               write (*,*) npt,mpt,k19
               npt = npt - 1
               exit 
            end if 
            jdv(npt) = jmin(i)
         end do
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin (jdv,npt,k19)
c                                 make a pointer to the original refinement 
c                                 point
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do 

      else if (npt.ne.icp) then 
c                                 done iterating, try to correct bad results
         if (npt.gt.icp) then 
c                                 look through for zero amount stable phases
            jzero = 0
   
            do i = 1, npt
               if (x(jdv(i)).eq.0d0) jzero = jzero + 1
            end do

            if (jzero.ge.npt-icp) then 

               jzero = npt - icp 
               mpt = npt
               npt = 0
c                                 eliminate zero phase
               do i = 1, mpt 

                  if (x(jdv(i)).eq.0d0.and.jzero.gt.0) then 
                     jzero = jzero - 1
                     cycle
                  end if 

                  npt = npt + 1
                  jdv(npt) = jdv(i)

               end do 

            else 
c                                 no hope, too many non-zero amount phases     
               idead = 1

            end if
 
         else if (icp-npt.eq.1.and.gmin.ne.0) then 
c                                 missing just one point, use gmin
          
            npt = icp
            jdv(icp) = gmin

         else 

            idead = 1

         end if 

         if (idead.eq.0) then
c                                 sort the phases 
            call sortin (jdv,npt,k19)
c                                 save the amounts
            do i = 1, npt
               amt(i) = x(jdv(i))
            end do 

         end if 

      end if 

      end