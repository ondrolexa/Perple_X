c routines only called by vertex/meemum, could be combined with nlib.f

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

      integer i,liw,lw,k,idead,inc,lphct

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),clamda(k1+k5),w(lw),oldt,oldp

      integer is(k1+k5),iw(liw)

      logical quit, abort
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

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

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer tphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct,jpt

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer idegen, idg(k5), jcp, jin(k5)
      common/ cst315 /idegen, idg, jcp, jin

      logical abort1
      common/ cstabo /abort1

      save ax, x, clamda, w, is, iw
c-----------------------------------------------------------------------
      idegen = 0
      jcp = 0

      if (.not.usv) then 
c                                 degeneracy test
         do k = 1, icp 
            if (b(k).eq.0d0) then 
               idegen = idegen + 1
               idg(idegen) = k
            else 
               jcp = jcp + 1
               jin(jcp) = k
            end if 
         end do

         inc = istct - 1

         oldt = t
         oldp = p
c                                logarithmic_p option
         if (lopt(14)) p = 1d1**p
c                                t_stop option
         if (t.lt.nopt(12)) t = nopt(12)

         call gall

         do k = 1, jphct
            c(k) = g(k+inc)/ctot(k+inc)
         end do
c                                load the adaptive refinement cpd g's
         do k = 1, jpt 
            g2(k) = c(k)
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

      else if (hcp.eq.1.or.iopt(10).eq.0.or.isoct.eq.0.or.usv) then 
c                                 no refinement, find the answer
         call yclos0 (x,is,jphct) 
c                                 final processing, .true. indicates static
         call rebulk (.true.,abort)

      else
c                                 save lphct to recover static solution if
c                                 no refinement 
         lphct = jphct 
c                                 find discretization points
c                                 for refinement

         call yclos1 (clamda,x,is,jphct,quit)
c                                 returns quit if nothing to refine
         if (quit) then 
c                                 final processing, .true. indicates static
            call rebulk (.true.,abort)

         else
c                                 initialize refinement point pointers
            do i = 1, ipoint
               hkp(i) = 0 
            end do 
c                                 reoptimize with refinement
            call reopt (idead)
c                                 final processing, .false. indicates dynamic
            if (idead.eq.0) then 

               call rebulk (.false.,abort)

               if (abort) then
c                                 bad solution (lagged speciation) identified
c                                 in avrger
                  call lpwarn (102,'LPOPT0')
                  if (iopt(22).lt.2) idead = 102

               end if 

               if (abort1) then 

                  idead = 104

               end if

            else if (idead.eq.-1) then
c                                 hail mary
               jphct = lphct
               idead = 0

               call yclos0 (x,is,jphct) 
               call rebulk (.true.,abort)

            end if 

         end if 

      end if 
     
      if (.not.usv) then 
         t = oldt
         p = oldp
      end if 

      end 

      subroutine reopt (idead)
c-----------------------------------------------------------------------
c reopt - given the results of an initial optimization for lpopt, reopt
c iteratively refines the solution by generating pseudocompounds in the
c neighborhood of the initial optimization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw, lw, iter, idead, jstart, opt, kter, kitmax 

      logical quit, kterat

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)

      double precision  ax(k5), x(k21), clamda(k21+k5), w(lw)

      integer is(k21+k5), iw(liw)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision xa,b,xc
      common/ cst313 /xa(k5,k1),b(k5),xc(k1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      quit = .false.
      opt = npt
      kterat = .false.
      kitmax = 0
      kter = 0 

c                                 --------------------------------------
c                                 first iteration
      call resub (1,kterat)


      if (jphct.eq.jpt) then
c                                 if nothing to refine, set idead 
c                                 to recover previous solution,
c                                 DEBUG DEBUG set to error 102 because
c                                 likely failed aqlagd
         idead = 102
         return

      end if

      if (kterat) kitmax = iopt(33)
      iter = 2

      do
c                                 iter is incremented before the operations,
c                                 i.e., on the nth iteration, iter is n+1
         if (kter.gt.kitmax) then 
            iter = iter + 1
            kter = 0
         end if
c                                 set quit flag
         if (iter.gt.iopt(10).and.kter.eq.kitmax) quit = .true.
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

            call lpwarn (idead,'REOPT')
            exit

         end if

         kter = kter + 1
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,iter,opt,idead,quit)

         if (idead.gt.0) then 

            call lpwarn (idead,'REOPT')
            exit 

         end if 
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the xcoor array.
         call saver

         if (quit) exit
c                                 generate new pseudocompounds
         call resub (iter,kterat)

      end do

      end

      subroutine resub (iter,kterat)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompounds around each refinenent 
c point during adaptive optimization. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      logical bad, kterat

      integer i, ids, lds, id, kd, iter, igood

      double precision res0

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer ikp
      common/ cst61 /ikp(k1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ntot,npairs
      common/ cst86 /ntot,npairs

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn
c----------------------------------------------------------------------
c                                 iteration dependent resolution
      res0 = nopt(24)/nopt(21)**iter
c                                 set dynamic array counters:
      jphct = jpt
      jcoct = 1
c                                 reset refinement point flags
      do i = 1, jpt
         hkp(i) = 0
      end do 
c                                 loop on previous stable phases
c                                 refine as necessay:
      lds = 0 

      do kd = 1, npt

         if (iter.eq.1) then 
c                                 static array pointer is
            id = jdv(kd) + istct - 1
c                                 solution model pointer is
            ids = ikp(id)
c                                 refine if a solution
            if (ids.eq.0) cycle 
c                                 special (pointless) iterations?
            if (lopt(32).and.iopt(33).gt.0) then
               if (ksmod(ikp(id)).eq.39) kterat = .true.
            end if
c                                 get the refinement point composition
            if (id.gt.ipoint) then 
               call getolx (ids,id)
            else
               if (.not.lopt(39)) cycle
               call endmmx (kd,id,ids)
            end if

         else
c                                 use pointer array lkp this uses 
c                                 negative values to index static
c                                 compounds and positive values to
c                                 point to solution models
            id = lkp(kd)

            if (id.lt.0) then

               ids = ikp(-id)
               if (ids.eq.0.or..not.lopt(39)) cycle
c                                 endmember refinement point:
c                                 get refine point composition
               call endmmx (kd,-id,ids)

            else

               ids = id
c                                 solution refinement point:
               call getxy0 (ids,kd)

            end if 

         end if 
c                                  get the subdivision limits:
         call sublim (ids,res0)
c                                  do the subdivision
         call subdiv (ids,.true.)
c                                  set solution model parameters for
c                                  gsol1, don't call if the previous
c                                  refinement point was the same solution.
         if (ids.ne.lds) call ingsol (ids)

         lds = ids

         igood = 0 

         do i = 1, ntot
c                                   increment jphct, 
c                                   load jkp[ids], hkp[i], local
c                                   and global composition arrays
            call loadgx (kd,i,ids,bad)

            if (bad) cycle

            igood = igood + 1

         end do
c                                    special call to make H2O for
c                                    lagged speciation, this is necessary
c                                    because non-linear stretching can prevent
c                                    fluid composition from reaching pure water.
         if (igood.ne.0.and.ksmod(ids).eq.39) then 
c                                    load prism with water coordinate
            do i = 1, ndim(1,ids)
               prism(i) = 0d0
            end do 

            call loadgx (kd,1,ids,bad)

         end if

      end do

      end

      subroutine endmmx (ld,jd,ids)
c----------------------------------------------------------------------
c generate compositional coordinates (x(i,j) array) for endmembers 
c during outrefine. if iter = 1, id is the static array index of the
c endmember, else iter is the dynamic array index. ids is the solution
c model index. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, ids, jd, kd, ld

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jend
      common/ cxt23 /jend(h9,m4)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer tphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct,jpt

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)
c----------------------------------------------------------------------
c                                set refinement point index
      hkp(jd) = ld 
c                                locate the endmember in the solution
      do i = 1, lstot(ids)
         if (jend(ids,2+i).eq.jd) then
            kd = i
            exit
         end if 
      end do 
c                                initialize all x's
      do i = 1, ostg(ids)
         do j = 1, ispg(ids,i)
            x(i,j) = 0d0
         end do
      end do
c                                 assign endmember x's
      if (kd.le.mstot(ids)) then 
c                                 normal simplex/prism
         do i = 1, istg(ids)
            x(i,kmsol(ids,kd,i)) = 1d0
         end do 

      else 
c                                 simplex with a prismatic vertex
c DEBUG DEBUG
         x(ostg(ids),kd-mstot(ids)) = 1d0 
c       write (*,*) 'groink ',x(ostg(ids),kd-mstot(ids))
c       pause

      end if

      end

      subroutine csol (id,bad)
c-----------------------------------------------------------------------
c csol computes chemical composition of solution id from the macroscopic
c endmember fraction array y or p0a (cxt7), these arrays are prepared by a prior
c call to function gsol. the composition is loaded into the array cp2 at
c position jphct.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, id

      logical bad, degen

      double precision ctot2

      external degen

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 adaptive coordinates
      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer jend
      common/ cxt23 /jend(h9,m4)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer jnd
      double precision aqg,qq,rt
      common/ cxt2 /aqg(m4),qq(m4),rt,jnd(m4)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
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

      else if (ksmod(id).eq.20) then 

         do i = sn1, nqs

            k = jnd(i) - aqst

            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + y(i) * aqcp(j,k)
            end do 

            ctot2 = ctot2 + y(i)*aqtot(k)

         end do 

         do i = 1, ns 

            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + y(i) * cp(j,jnd(i))
            end do 

            ctot2 = ctot2 + y(i)*ctot(jnd(i))

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
c                                  a phase with a null composition may appear
c                                  as an endmember of a solution in a calculation
c                                  with mobile components:

c                                  sept 22 2017: previously null compositions were
c                                  given unstable properties, bad flag added this 
c                                  date along with degeneracy check. 
      bad = .false.

      if (ctot2.ne.0d0) then
c                                  normalize the composition and free energy
         g2(jphct) = g2(jphct)/ctot2
         c2tot(jphct) = ctot2

         do j = 1, icp
            cp2(j,jphct) = cp2(j,jphct)/ctot2
         end do
c                                  degeneracy test
         if (degen(jphct,2)) then 
c            bad = .true.
            return
         end if 

      else 
c                 DANGER DEBUG 
         write (*,*) 'ctot2?', ctot2
         bad = .true.

      end if

      end 

      subroutine sortin 
c-----------------------------------------------------------------------
c sort the first npt values of jdv
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, imin

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
c----------------------------------------------------------------------
      do j = 1, npt-1

         imin = jdv(j)

         do i = j+1, npt

            if (jdv(i).lt.imin) then 
               imin = jdv(i)
               jdv(i) = jdv(j)
               jdv(j) = imin
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

      integer idead, iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, 
     *        iwarn03

      character char*(*)

      double precision c

      save iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, iwarn03

      data iwarn91, iwarn42, iwarn90, iwarn01, iwarn02, iwarn03/6*0/
c----------------------------------------------------------------------
c                                             look for errors
      if (idead.eq.2.or.idead.gt.4.and.idead.lt.100
     *                            .and.iwarn91.lt.6) then 
c                                             unbounded solution, or
c                                             other programming error.
         call warn (91,c,idead,char) 
         iwarn91 = iwarn91 + 1
         if (iwarn91.eq.5) call warn (49,c,91,'LPWARN')

      else if (idead.eq.3.and.iwarn42.lt.6) then 
c                                             no feasible solution
         call warn (42,c,idead,char)
         iwarn42 = iwarn42 + 1
         if (iwarn42.eq.6) call warn (49,c,42,'LPWARN')

      else if (idead.eq.4.and.iwarn90.lt.6) then 
c                                             iteration count exceeded,
c                                             probable cause no feasible
c                                             solution.
         call warn (90,c,idead,char) 
         iwarn90 = iwarn90 + 1
         if (iwarn90.eq.5) call warn (49,c,90,'LPWARN')

      else if (idead.eq.101.and.iwarn01.lt.10) then

          iwarn01 = iwarn01 + 1
          call warn (100,c,101,'under-saturated solute-component.'
     *              //' To output result set aq_bad_result to 102')
          if (iwarn01.eq.10) call warn (49,c,101,'LPWARN')

      else if (idead.eq.102.and.iwarn02.lt.10) then

          iwarn02 = iwarn02 + 1
          call warn (100,c,102,'pure and impure solvent phases '//
     *              'coexist within solvus_tolerance. '//
     *              'To output result set aq_bad_result to 101')
          if (iwarn02.eq.10) call warn (49,c,102,'LPWARN')

      else if (idead.eq.103.and.iwarn03.lt.10) then

          iwarn03 = iwarn03 + 1
          call warn (100,c,103,'pure and impure solvent phases '//
     *              'coexist. To output result set aq_bad_result.')
          if (iwarn03.eq.10) call warn (49,c,103,'LPWARN')

      end if

      end 

      subroutine saver 
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound x(i,j) compositions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      integer i, j, k, kcoct, id, ids, itic
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

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer jpoint, jstct
      common/ cxt60 /jpoint,jstct

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)
c----------------------------------------------------------------------
      kcoct = 0

      do i = 1, npt

         id = jdv(i)

         if (id.lt.jpoint) then 
            lkp(i) = -(id + jstct)
            cycle
         end if 

         ids = jkp(id)
         lkp(i) = ids
c                                 cycle on a compound
         if (ids.lt.0) cycle
c                                 it's a solution:
         lcoor(i) = kcoct
         itic = 0

         do j = 1, ostg(ids)

            do k = 1, ndim(j,ids)

               itic = itic + 1

               if (kcoct+itic.gt.k22) call error (60,ctotal,k22,'saver')

               ycoor(lcoor(i)+itic) = zcoor(jcoor(id)+itic)

            end do

         end do

         kcoct = kcoct + itic

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

      double precision xt 
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)
c----------------------------------------------------------------------
      jcoor = lcoor(id)

      do i = 1, ostg(ids)

         xt = 0d0 

         do j = 1, ndim(i,ids)
            jcoor = jcoor + 1
            x(i,j) = ycoor(jcoor)
            xt = xt + ycoor(jcoor)
         end do 

         x(i,j) = 1d0 - xt

      end do

      end 

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
      common/ cst57 /dcp(k5,k19),soltol
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
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
c ids, intermediate solution values. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id1,id2,i

      integer icomp,iphct,icp,istct
      common/ cst6 /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt
c-----------------------------------------------------------------------
      solvs2 = .false.

      do i = 1, icp

         if (dabs(cp2(i,id1) - cp2(i,id2)).gt.soltol) then
            solvs2 = .true.
            exit
         end if

      end do 

      end 

      logical function solvs4 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of a lagged
c aqueous solution based 

c ids, called only for final solution vales by avrger.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c-----------------------------------------------------------------------
      solvs4 = .false.

      do i = 1, ns

         if (dabs(x3(id1,1,i) - x3(id2,1,i)).gt.soltol) then 
            solvs4 = .true.
            exit 
         end if 
         
      end do 

      end 

      subroutine avrger (abort)
c----------------------------------------------------------------------
c avrger combines discretization points into a single solution
c composition. on output

c     np  - is the number of solutions, 
c     ncpd - is the number of true compounds
c     ntot - np+ncpd

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are)
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical check, bad, quit, notaq, abort

      integer idsol(k5),kdsol(k5,k5),ids,xidsol,xkdsol,irep,
     *        i,j,jdsol(k5,k5),jd,k,l,nkp(k5),xjdsol(k5),kk

      double precision bsol(k5,k5),cpnew(k5,k5),xx,xb(k5),msol,
     *                 bnew(k5),xnew(k5,mst,msp),ncaq(k5,l10),ximp

      logical solvs1, solvs4
      external solvs1, solvs4
c                                 -------------------------------------
c                                 global variables:
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                  x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      double precision fwt
      common/ cst338 /fwt(k10)

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      logical abort1
      common/ cstabo /abort1
c-----------------------------------------------------------------------
      abort = .false.
c                                first check if solution endmembers are
c                                among the stable compounds:
      do i = 1, ntot
c                                initialize kdsol, this was not done 
c                                before nov 17, 2017; god only knows
c                                how it worked...
         kdsol(1,i) = 0
c                                locate solution endmembers:
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

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0
c                                solvus_tolerance
      soltol = nopt(8)

      do i = 1, ntot
          
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = amt(i)
            kdsol(ntot,ncpd) = nkp(i)
            jdsol(ntot,ncpd) = i

         else 

            if (lopt(32).and.ksmod(nkp(i)).eq.39) then

               notaq = .false.
c                                 get lagged speciation
c                                 loaded into caq(i,1:ns+aqct)
               do k = 1, ns
                  y(k) = x3(i,1,k)
               end do 

               if (abort1) then 
                  quit = .true.
                  abort = .true.
                  cycle 
               end if 

               if (quack(jdv(i))) then 
c                                 pure solvent phase
                  msol = 0d0

                  do k = 1, ns

                     caq(i,k) = y(k)
c                                 solvent molar weight
                     msol = msol + y(k) * fwt(jnd(k))

                  end do 

                  do k = sn1, nat
                     caq(i,k) = 0d0
                  end do
c                                 total molality
                  caq(i,na2) = 1d0/msol
                  caq(i,na3) = msol

               else
c                                 impure solvent, get speciation
                  call aqlagd (i,bad,.true.)

               end if

            else

               notaq = .true.

            end if

            quit = .false.

            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (kdsol(j,1).eq.nkp(i)) then

                  kk = jdsol(j,idsol(j))
c                                 if match check for a solvus
                  if (notaq) then 
                     if (solvs1(i,kk,nkp(i))) cycle
                  else
c                                  special solvus test based on solvent 
c                                  speciation for lagged aq model.
                     if (solvs4(i,kk)) cycle

                     if (iopt(22).lt.2) then 
c                                  check pure and impure solvent coexistence

                        if (caq(i,na1).eq.0d0.and.caq(kk,na1).ne.0d0.or.
     *                      caq(i,na1).ne.0d0.and.caq(kk,na1).eq.0d0) 
     *                                                              then 
c                                  pure solvent and impure solvent coexist
                            abort = .true.
                            return

                        end if

                     end if 

                  end if 
c                                 the pseudocompound matches a solution
c                                 found earlier.
                  idsol(j) = idsol(j) + 1
                  bsol(j,idsol(j)) = amt(i)
                  jdsol(j,idsol(j)) = i
                  quit = .true.
                  exit

               end if

            end do

            if (quit) cycle 
c                                 the pseudocompound is a new solution 
c                                 phase.
            np = np + 1
            idsol(np) = 1
            kdsol(np,1) = nkp(i)
            jdsol(np,1) = i
            bsol(np,1) = amt(i)

         end if

      end do
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
         ids = kdsol(i,1)

         do j = 1, icomp
            cpnew(j,i) = 0d0
         end do 

         do j = 1, ostg(ids)
            do k = 1, ispg(ids,j)
               xnew(i,j,k) = 0d0
            end do 
         end do

         if (lopt(32).and.ksmod(ids).eq.39) then 
c                               lagged speciation
            do k = 1, nat
               ncaq(i,k) = 0d0
            end do 

         end if 

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + amt(jdsol(i,j))
         end do
c                                in case pure and impure solvent is going to be averaged
c                                count fraction of impure solvent
         ximp = 0d0

         do j = 1, idsol(i)

            jd = jdsol(i,j)
c                                conditional for zero-mode stable phases
            if (bnew(i).gt.0d0) then 
               xx =  amt(jd)/bnew(i)
            else 
c                                this is gonna be a disaster if idsol(i) > 1
               xx = 1d0
            end if 
c                                save the new compositions
            do k = 1, icomp
               cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
            end do

            do k = 1, ostg(ids)
               do l = 1, ispg(ids,k)
                  xnew(i,k,l) = xnew(i,k,l) + xx*x3(jd,k,l)
               end do
            end do

            if (lopt(32).and.ksmod(ids).eq.39) then

               if (caq(jd,na1).ne.0d0) ximp = ximp + xx
c                                lagged speciation (1:nsa), ionic strength (na1), total
c                                molality (na2), solvent mass (na3), err_log_kw (na4)
c                                pH, Delta_pH, solute molality, epsilon (nat)
               do k = 1, nat
                  ncaq(i,k) = ncaq(i,k) + xx*caq(jd,k)
               end do

            end if

         end do

         if (lopt(32).and.ksmod(ids).eq.39.and.ximp.gt.0d0) then
c                                 renomalize err_log_kw, pH, Delta_pH, epsilon
            ncaq(i,na3+1) = ncaq(i,na3+1)/ximp
            ncaq(i,na3+2) = ncaq(i,na3+2)/ximp
            ncaq(i,na3+3) = ncaq(i,na3+3)/ximp
            ncaq(i,nat) = ncaq(i,nat)/ximp

         end if

      end do
c                                 now reform the arrays kdv and b
      do i = 1, np

         amt(i) = bnew(i)
         kkp(i) = kdsol(i,1)
         ids = kkp(i)

         do j = 1, icomp
            cp3(j,i) = cpnew(j,i)
         end do

         do j = 1, ostg(ids)
            do k = 1, ispg(ids,j)
c                                 set x's for sollim
               x(j,k) = xnew(i,j,k)
c                                 set x's for global storage
               x3(i,j,k) = x(j,k) 

            end do
         end do

         if (lopt(32).and.ksmod(ids).eq.39) then 
c                               lagged speciation, ionic strength, tot molality
c                               and solvent mass.
            do k = 1, nat
               caq(i,k) = ncaq(i,k)
            end do 

         end if 
c                                 check composition against solution model ranges
c                                 if auto_refine is on:
         call sollim (ids)

      end do

      if (.not.usv) then 

         do i = 1, ncpd
            k = np + i
            l = kdsol(ntot,i)
            amt(k) = bsol(ntot,i)
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
               if (kdsol(ntot,j).eq.l) amt(k) = amt(k) + bsol(ntot,j)
            end do 

            amt(k) = bsol(ntot,i)
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

      integer ids, i, j

      double precision stinc

      external stinc 

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c----------------------------------------------------------------------
c                                 set stable flag
      stable(ids) = .true.
c                                 check x-ranges
      do i = 1, ostg(ids)

         do j = 1, ndim(i,ids)

            if (ksmod(ids).eq.20) then 
               if (j.eq.ns) cycle 
            end if 
c                                 low limit:
            if (x(i,j).lt.xlo(j,i,ids)) then

               xlo(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).gt.xmno(ids,i,j).and.
     *            (x(i,j).le.xmng(ids,i,j).and..not.lopt(3))) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmng(ids,i,j) = xmng(ids,i,j) - xncg(ids,i,j)

                  else 
c                                 assymmetric stretching towards xmin
                     xmng(ids,i,j) = 
     *                   stinc (x(i,j),-1d0/xncg(ids,i,j),ids,i,j)

                  end if

                  if (xmng(ids,i,j).lt.0d0) xmng(ids,i,j) = 0d0

                  limit(ids) = .true.

               end if

            end if 
c                                 high limit:
            if (x(i,j).gt.xhi(j,i,ids)) then

               xhi(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).lt.xmxo(ids,i,j).and.
     *            (x(i,j).ge.xmxg(ids,i,j).and..not.lopt(3))) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmxg(ids,i,j) = xmxg(ids,i,j) + xncg(ids,i,j)

                  else 
c                                 assymmetric stretching 
                     xmxg(ids,i,j) = 
     *                   stinc (x(i,j),1d0/xncg(ids,i,j),ids,i,j)

                  end if 

                  if (xmxg(ids,i,j).gt.1d0) xmxg(ids,i,j) = 1d0
                  limit(ids) = .true.

               end if

            end if

         end do

      end do

      end 

      subroutine sorter (kdbulk,ico,jco,output)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,m,kdbulk,ico,jco,ids,ioct,inct

      logical output, reord, match, nomtch, ok 

      double precision cpt(k5,k5),xt(k5,mst,msp),bt(k5),caqt(k5,l10)
c                                 x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 look for a match with known assemblages
      match = .false.
 
      do i = 1, iasct

         if (np.ne.iavar(1,i).or.ncpd.ne.iavar(2,i)) cycle 

         nomtch = .false.

         do j = 1, ntot

            ok = .false.

            do k = 1, ntot 

               if (idasls(k,i).eq.kkp(j)) then 

                  ok = .true.
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                  inct = 0 
                  ioct = 0 

                  do l = 1, np
                     if (kkp(l).eq.kkp(j)) inct = inct + 1
                     if (idasls(l,i).eq.kkp(j)) ioct = ioct + 1
                  end do 

                  if (ioct.ne.inct) then
                     nomtch = .true.  
                     exit 
                  end if 

               end if 

            end do

            if (.not.ok) nomtch = .true. 
            if (nomtch) exit 

         end do 

         if (nomtch) cycle

         match = .true.
c                                 check if reordering is necessary
         reord = .false.

         do j = 1, ntot
            if (kkp(j).eq.idasls(j,i)) cycle
            reord = .true.
            exit
         end do

         if (reord) then 
c                                 reorder the result arrays of the
c                                 current occurence to match initial 
c                                 occurence:
            do j = 1, ntot

               do k = 1, ntot

                  ids = kkp(k)

                  if (ids.eq.idasls(j,i)) then
c                                 load temporary array
                     bt(j) = amt(k)

                     if (ids.gt.0) then 

                        do l = 1, icomp
                           cpt(l,j) = cp3(l,k)
                        end do

                        do l = 1, istg(ids)
                           do m = 1, ispg(ids,l)
                              xt(j,l,m) = x3(k,l,m) 
                           end do 
                        end do

                        if (lopt(32).and.ksmod(ids).eq.39) then 

                           do l = 1, nat
                              caqt(j,l) = caq(k,l)
                           end do 

                        end if

                     end if 
c                                 this eliminates immiscible phases
                     kkp(k) = 0

                     exit 
 
                  end if 

               end do 

            end do
c                                 reload final arrays from temporary
            do j = 1, ntot

               amt(j) = bt(j)
               ids = idasls(j,i)
               kkp(j) = ids

               if (ids.gt.0) then 

                  do k = 1, icomp
                     cp3(k,j) = cpt(k,j)
                  end do

                  do k = 1, istg(ids)
                     do l = 1, ispg(ids,k)
                        x3(j,k,l) = xt(j,k,l) 
                     end do 
                  end do

                 if (lopt(32).and.ksmod(ids).eq.39) then 

                     do l = 1, nat
                        caq(j,l) = caqt(j,l)
                     end do 

                  end if

               end if 
            end do 

         end if 

         if (ibulk.gt.k2) call error (183,0d0,k2,'SORTER')
         ibulk = ibulk + 1
         iap(ibulk) = i
         kdbulk = ibulk

         exit  

      end do 

      if (.not.match) then 
c                                 the assemblage is new:
         iasct = iasct + 1
         if (iasct.gt.k3) call error (184,0d0,k3,'SORTER')

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

      end if 
                                
      if (output) call outbl1 (ico,jco)
     
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
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 x-coordinates for the final solution
      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 i/o
      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c----------------------------------------------------------------------
      if (io4.eq.1) return
c                                graphics output  
      write (n5,'(3(i8,1x))') ico,jco,iap(ibulk)
c                                phase molar amounts
      write (n5,1010) (amt(i),i=1,np+ncpd)
c                                solution phase compositions
      do i = 1, np
         ids = kkp(i)
         write (n5,1010) ((x3(i,j,k),k=1,ispg(ids,j)),j=1,ostg(ids))
c                                lagged speciation
         if (ksmod(ids).eq.39.and.lopt(32)) write (n5,1010) 
     *                                            (caq(i,j),j=1,nat)

      end do 
c                                dependent potentials
      if (jpot.ne.1) write (n5,1010) (mu(i),i=1,jbulk)

1010  format (20(g16.8,1x))

      end 

      subroutine yclos1 (clamda,x,is,jphct,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement. this routine is only called as preparation
c for iterative refinement.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, i, j, k, is(*), idsol(k5), kdv(h9), nsol, ids,
     *        mpt, iam, id, inc, jdsol(k5,k5), kdsol(k5), max

      external ffirst, degen 

      logical degen, solvus, quit, news, solvnt(1)

      double precision clamda(*), x(*), slam(h9)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision wmach(9)
      common /ax02za/wmach

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isoct
      common/ cst79 /isoct

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer ikp
      common/ cst61 /ikp(k1)

      integer tphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct,jpt

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      npt = 0
      nsol = 0
      inc = istct - 1
      quit = .true.
c                                 solvus_tolerance_II, 0.25
      soltol = nopt(25)

      do i = 1, jphct

         if (is(i).ne.1) then

c            if (x(i).lt.zero) cycle
c                                 degeneracy test
c            if (degen(i,1)) cycle 
c                                 make a list of found phases:
            id = i + inc
c                                 currently endmember compositions are not 
c                                 refined (this is probably a mistake, but 
c                                 seems to work fine), so use id > ipoint
c                                 to identify a refineable point
            ids = ikp(id) 

            if (ids.ne.0) then 

               if (id.gt.ipoint.or.
     *             ksmod(ids).eq.39.and.lopt(32)) then 
c                                 a pseudocompound
                  quit = .false.
                  news = .true.

                  do j = 1, nsol

                     if (ids.ne.idsol(j)) cycle 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = id
                     news = .false.
                     exit 

                  end do

                  if (news) then 
c                                 new phase, add to list
                     nsol = nsol + 1
                     idsol(nsol) = ids
                     jdsol(1,nsol) = id
                     kdsol(nsol) = 1

                  end if

               end if 

            end if 
c                                 new point, add to list
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)

            if (lopt(34)) then

               if (npt.eq.1) write (*,'(/,a,i2)') 'iteration ',0
               if (ikp(id).ne.0) then 
                  call dumper (1,i,0,ikp(i),x(i),clamda(i))
               else 
                  call dumper (1,i,0,-id,x(i),clamda(i))
               end if 
            end if

         end if 

      end do 
c                                 get mus for lagged speciation
      call getmus (1,0,solvnt,.false.,.false.)

      do i = 1, isoct
         slam(i) = 1d99
         kdv(i) = 0
      end do 
c                                 perp 6.6.3, make a list of the least 
c                                 metastable composition of each solution.
      do 20 i = jpt+1, jphct

c DEBUG why was this here? added ~6.7.6, removed april 21, 2017
c i think clamda(i).lt.0 allows degenerate compositions (and probably 
c therefore the 6.7.6 version may be better, on the bright side with
c it removed the solution composition gets refined (if endmember comps are
c being allowed, see ldsol code); restored again april 2017. removed sep 2018.

         if (is(i).ne.1) cycle

         if (degen(i,1)) cycle

         id = i + inc 
         iam = ikp(id)

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
c                                the composition is acceptable
            slam(iam) = clamda(i)
            kdv(iam) = i

         end if

20    continue
c                                 load the metastable solutions into kdv
      mpt = 0 

      do i = 1, isoct

         if (kdv(i).eq.0) cycle
         mpt = mpt + 1
         kdv(mpt) = kdv(i)
         slam(mpt) = slam(i)
         quit = .false.

      end do 

      if (mpt.le.iopt(31)) then 
c                                 less metastable refinement points than
c                                 iopt(31)
         max = mpt

      else 
c                                 sort the metastable points to
c                                 find the most stable iopt(31) points
         max = iopt(31)

         call ffirst (slam,kdv,1,mpt,max,h9,ffirst)

      end if 
 
      do i = 1, max

         jdv(npt+i) = kdv(i)

         if (lopt(34)) then

            id = kdv(i)

            if (ikp(id).ne.0) then 
               call dumper (1,id,0,ikp(id),x(id),clamda(id))
            else 
               call dumper (1,id,0,-id,x(id),clamda(id))
            end if 

         end if

      end do

      if (quit) then 
c                                 zero mode filter and 
c                                 save amounts for final processing
         mpt = npt
         npt = 0

         do i = 1, mpt
            if (x(jdv(i)).lt.nopt(9)) cycle 
            npt = npt + 1
            jdv(npt) = jdv(i)
            amt(npt) = x(jdv(i))
         end do 

      else 

         npt = npt + max
c                                 sort the phases, why? don't know, but it's 
c                                 necessary
         call sortin

      end if 

      end 

      recursive subroutine ffirst (a, ind, left, right, k, n, dumsub)
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

      subroutine initlp 
c--------------------------------------------------------------------
c initialize arrays and constants for lp minimization
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,inc,id,jk

      logical bad

      double precision u,s,vol,tot 

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

      integer ldt,ldq
      common /be04nb/ldt,ldq

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn

      integer jpoint, jstct
      common/ cxt60 /jpoint,jstct

      integer tphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),tphct,jpt

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c-----------------------------------------------------------------------
      inc = istct - 1

      tloop = 40
      ploop = 40
      dv(1) = (vmax(1)-vmin(1))/(ploop-1)
      dv(2) = (vmax(2)-vmin(2))/(tloop-1)
c                                 load arrays for lp solution
      jphct = iphct - inc
      jpt = ipoint - inc

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
c                                 static composition array
         do i = 1, jphct
            id = i + inc
            iam(i) = id
            do j = 1, icp
               a(j,i) = cp(j,id)/ctot(id)
            end do
         end do
c                                 load all compounds into the 
c                                 the dynamic composition array
         do id = istct, ipoint

            i = id - inc 
c                                 jkp indicates which phase a point is associated with
            jkp(i) = -id
            hkp(i) = 0

            do j = 1, icp
               cp2(j,i) = a(j,i)
            end do

         end do
c                                 locate last point in dynamic arrays and set increment
c                                 to static array
         jpoint = ipoint - inc
         jstct = inc 

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

      end 

      subroutine yclos0 (x,is,jphct)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,jphct,is(*)

      logical solvnt(1)

      double precision x(*)
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      npt = 0

      do i = 1, jphct

         if (is(i).eq.1.or.x(i).lt.nopt(9)) cycle  
c                                 acceptable cases 0 active, between bounds
c                                                  2 active, upper bound 
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)
 
      end do

      call getmus (1,0,solvnt,.false.,.true.)

      end

      subroutine setx3 (ind,id,ids)
c----------------------------------------------------------------------
c subroutine to recover prismatic solution compositions (x(i,j))
c from the xco array loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, kcoor, ind

      double precision xt

      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer kd, na1, na2, na3, nat
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,nat,kd

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      if (id.gt.ipoint) then 
c                                 a normal solution
         if (istg(ids).eq.1) then 
c                                 one site solution
            xt = 0d0 

            do j = 1, nstot(ids) - 1
               x3(ind,1,j) = xco(jco(id)+j) 
               xt = xt + x3(ind,1,j)
            end do 

            x3(ind,1,j) = 1d0 - xt 

         else if (ispg(ids,1).gt.1) then 
c                                 multi-site solution
            kcoor = ico(id)

            do i = 1, ostg(ids)

               xt = 0d0 

               do j = 1, ndim(i,ids)
                  kcoor = kcoor + 1
                  x3(ind,i,j) = xco(kcoor)
                  xt = xt + xco(kcoor)
               end do 

               x3(ind,i,j) = 1d0 - xt 

            end do 

         else 
c                                 a dummy site
            x3(ind,1,1) = 1d0
            xt = 0d0 

            do j = 1, nstot(ids)
               x3(ind,2,j) = xco(jco(id)+j) 
               xt = xt + x3(ind,2,j)
            end do 

            x3(ind,2,j) = 1d0 - xt 

         end if 

      else 
c                                 an endmember 
         call endcp (ind,id,ids)

      end if

      end 

      double precision function ginc0 (dt,dp,id)
c-----------------------------------------------------------------------
c id indicates a static pseudocompound, use ginc for dynamic compound
c-----------------------------------------------------------------------
      implicit none

      double precision dt, dp, gphase

      external gphase

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      p = p + dp 
      t = t + dt 

      ginc0 = gphase (id)

      p = p - dp 
      t = t - dt

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

      subroutine yclos2 (clamda,x,is,iter,opt,idead,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. quit is true for final
c iteration.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      external ffirst

      integer i, id, is(*), jmin(k19), opt, kpt, mpt, iter, tic, 
     *        idead, j, k

      double precision clamda(*), clam(k19), x(*)

      logical stable(k19), solvnt(k19), quit, abort, test, good, bad

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision wmach(9)
      common /ax02za/wmach

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer idaq, jdaq
      logical laq
      common/ cxt3 /idaq,jdaq,laq

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      save tic
      data tic/0/
c----------------------------------------------------------------------
c                                 npt is the number of points refined
c                                 from the previous cycle, opt is the
c                                 number of points in the original 
c                                 solution.
      do i = 1, k19
         jmin(i) = 0 
         clam(i) = 1d99
         stable(i) = .false.
      end do

      abort = .false.
      test = .false.

      npt = 0

      do i = 1, jphct
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)

         if (is(i).ne.1) then
c                                 these tests destabilize chemical potential
c                                 back-calculation.
c                                 degeneracy test
c           if (degen(i,2)) cycle
c                                 mode test
c           if (x(i).lt.zero) cycle
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            amt(npt) = x(i)
            if (id.gt.0) stable(id) = .true.

            if (lopt(32)) then
c                                 for lagged aq speciation
c                                 classify as solvent/solid
               if (jkp(i).lt.0) then
c                                 setting abort to true signals 
c                                 test in getmus:
                  if (quack(-jkp(i))) then 
                     abort = .true.
                  end if 

                  if (ikp(-jkp(i)).eq.idaq) then
                     solvnt(npt) = .true.
                     test = .true.
                  else 
                     solvnt(npt) = .false.
                  end if

               else if (jkp(i).eq.idaq) then

                  solvnt(npt) = .true.
                  test = .true.

               else

                  solvnt(npt) = .false.

               end if

            end if

            if (lopt(34)) then
c                                 dump iteration details
               if (npt.eq.1) write (*,'(/,a,i2)') 'iteration ',iter-1
               call dumper (2,i,hkp(i),jkp(i),x(i),clamda(i))

            end if 

         else if (id.gt.0) then 
c                                 a metastable solution cpd
            if (clamda(i).lt.clam(id)) then
c DEBUG DEBUG DEBUG               this is a cheap way of eliminating
c                                 compositionally redundant refinement 
c                                 points, the problem is that the critical 
c                                 value of clamda is problem/machine dependent. 
c              if (clamda(i).lt.1d-7) cycle 
c                                 keep the least metastable point
               jmin(id) = i
               clam(id) = clamda(i)

            end if

         end if

      end do
c                                 get mu's for lagged speciation
      if (abort.and.iopt(22).eq.0) then 

          quit = .true.
          idead = 103

      else 

         if (test) abort = .true.

         call getmus (iter,iter-1,solvnt,abort,.false.)

         if (abort) then 
c                                 undersaturated solute component
            quit = .true.
c                                 don't set idead if iopt =1, this
c                                 allows output of the result.
            if (iopt(22).ne.1.and.iopt(22).ne.99) then 
               idead = 101
            else 
               call lpwarn (101,'YCLOS2')
            end if 

         end if

      end if 

      if (.not.quit) then
c                                 if not done iterating:
         kpt = 0
c                                 make a list of the solutions
         do i = 1, opt

            if (jmin(i).eq.0) then
               cycle
            else if (stable(hkp(jmin(i)))) then
c                                 contrary to what you might expect, this
c                                 definitely improves quality, presumably because
c                                 the metastable point will always be the closest
c                                 composition to the stable point and the resulting
c                                 overlap fogs the lp.

c                                 sep 18, nah it stops the list from being clogged
c                                 up with one phase
               cycle

            else
c                                 delete compositionally degenerate refinement points
               bad = .false. 

               do j = 1, npt + kpt 

                  if (jkp(jdv(j)).ne.jkp(jmin(i))) cycle

                  good = .false.
c                                 metastable point matches a refinement point, 
c                                 check composition
                  do k = 1, icp

                     if (dabs(cp2(k,jdv(j))-cp2(k,jmin(i))).gt.nopt(5))
     *                                                              then
                        good = .true.
                        exit

                     end if

                  end do

                  if (good) cycle

                  bad = .true.

                  exit

               end do

               if (bad) cycle

            end if 

            kpt = kpt + 1
            jmin(kpt) = jmin(i)
            clam(kpt) = clam(i)
c                                 temporarily use jdv(i>npt) for degeneracy test
            jdv(npt+kpt) = jmin(i)

         end do

         if (kpt.gt.iopt(31)) then 
c                                 sort if too many
            kpt = iopt(31)

            call ffirst (clam,jmin,1,kpt,iopt(31),k19,ffirst)

         end if 

         do i = 1, kpt

           npt = npt + 1
           jdv(npt) = jmin(i)

           if (lopt(34)) call dumper (2,jdv(npt),hkp(jdv(npt)),
     *                       jkp(jdv(npt)),x(jdv(npt)),clamda(jdv(npt)))

         end do
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin
c                                 make a pointer to the original refinement 
c                                 point
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do

      else
  
         mpt = npt
         npt = 0
c                                 check zero modes the amounts
         do i = 1, mpt

            if (x(jdv(i)).ge.nopt(9)) then

               npt = npt + 1
               amt(npt) = x(jdv(i))
               jdv(npt) = jdv(i)

            else if (x(jdv(i)).lt.-nopt(9).and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

      end if 

      end

      subroutine rebulk (static,abort)
c----------------------------------------------------------------------
c upon successful completion of an optimization with either static or
c dynamic pseudocompounds rebulk:
c     1) loads the generic arrays cp3, cptot, ctot3, kkp and x3
c     2) computes the amounts of saturated component phases.
c     4) checks for solvi and homogenizes miscible phases.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, j, k, id, tictoc

      logical static, bad, abort

      double precision c(k5),u

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c                                 hcp is different from icp only if usv
      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      double precision ctot
      common/ cst3  /ctot(k1)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      save tictoc
      data tictoc/0/
c----------------------------------------------------------------------
      do i = 1, npt

         if (static) then 

            id = iam(jdv(i))
c                                 set identifier flag
            if (id.le.ipoint) then
               kkp(i) = -id
            else  
               kkp(i) = ikp(id)
            end if 

            cptot(i) = ctot(id)

            do j = 1, icomp
               if (j.gt.icp.and.usv) exit
               cp3(j,i) = cp(j,id)
            end do
c                                 set the x3 array
            if (ikp(id).ne.0) call setx3 (i,id,ikp(id))

            if (usv) then 
c                                 usv calculations are currently
c                                 not set up for solutions and use
c                                 only static optimization
               call pt4sv (jam(jdv(i)))

               call getusv (u,cp3(icp1,i),cp3(jbulk,i),id,bad)

            end if 

         else 
c                                 getcmp assigns cp3, cptot, x3, and kkp
            call getcmp (i,jdv(i),jkp(jdv(i)))

         end if
c                                 convert normalized amounts to molar 
c                                 amounts
         amt(i) = ctotal*amt(i)/cptot(i)

      end do 

      if (jbulk.gt.icp.and.(.not.usv)) then  
c                                 get the amounts of the saturated phases:
         do i = icp+1, jbulk
c                                 k is the saturated component pointer
            k = i - icp
c                                 initialize bulk
            c(k) = cblk(i)
c                                 save chemical potentials from gproj
            mu(i) = us(k)
c                                 subtract the amount of the component in the 
c                                 phases in the thermodynamic c-space
            do j = 1, npt 
               c(k) = c(k)- amt(j)*cp3(i,j)
            end do 

         end do 

         do i = jbulk, icp+1, -1
c                                  cycle through the saturated phases
            npt = npt + 1
            id = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(id).eq.0) then 
               kkp(npt) = -id
            else 
               write (*,1000) 
               stop 
            end if 
c                                  amount of the staurated phase
            amt(npt) = c(i-icp)/cp(i,id)
c                                  warn on undersaturation
            if (amt(npt).lt.nopt(9)) then 
               if (amt(npt).lt.-nopt(9).and.tictoc.lt.5) 
     *                             call warn (2,amt(npt),i,'REBULK')
               npt = npt - 1
               tictoc = tictoc + 1
               exit
            end if 
c                                  remove the saturated phase from 
c                                  the bulk composition.
            do j = icp+1, i - 1
               c(j) = c(j) - amt(npt)*cp(j,id)
            end do 
c                                  load the saturated phase composition 
            do j = 1, icomp
               cp3(j,i) = cp(j,id)
            end do           

         end do

      end if 

      ntot = npt
c                                 test for solvi and average
c                                 homogeneous phases.
      call avrger (abort)
      
1000  format (/,'**error ver901** solutions not allowed in saturated ',
     *   'component composition space',/,'in adaptive optimization ',
     *   'calculations, this limitation can be removed upon request.')

      end

      subroutine getgc (lc,lg,iter) 
c----------------------------------------------------------------------
c iter is a flag indicating where the compositions are and is sort of 
c      related to the iteration count during optimization.
c jter is the true iteration count.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, iter

      double precision lc(k8,k8), lg(k8)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer hcp, idv
      common/ cst52  /hcp,idv(k7) 

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------

      do i = 1, npt

         id = jdv(i) 

         if (iter.gt.1) then

            do j = 1, hcp
               lc(i,j) = cp2(j,id)
            end do

            lg(i) = g2(id)

         else

            do j = 1, hcp
               lc(i,j) = a(j,id)
            end do

            lg(i) = c(id)

         end if

      end do

      end 

      subroutine getmus (iter,jter,solvnt,abort,quit) 
c----------------------------------------------------------------------
c iter is a flag indicating where the compositions are and is sort of 
c      related to the iteration count during optimization.
c jter is the true iteration count.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvnt(*), quit, bad, abort, cslut(k19), cslvt(k19)

      integer i, j, ier, ipvt(k8), iter, jter, imu(k8), kcp, lcp, 
     *        inp(k8)

      double precision comp(k8,k8), g, lc(k8,k8), lg(k8)

      character cname*5
      common/ csta4  /cname(k5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision p0, dz
      common/ cxt46 /p0, dz

      integer hcp, idv
      common/ cst52  /hcp,idv(k7) 

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer jtest,jpot
      common/ debug /jtest,jpot

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)

      integer idegen, idg(k5), jcp, jin(k5)
      common/ cst315 /idegen, idg, jcp, jin

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------

      if (.not.lopt(32)) then

         if ((iter.lt.iopt(10).or.quit).and.jpot.ne.0) return

      end if

      ier = 1
c                                 load c and g into a local array to 
c                                 avoid a myriad of conditionals
      call getgc (lc,lg,iter)
c                                 for lagged speciation 
      if (abort) then

         abort = .false.
c                                 isolc is the number of non-solvent components, 
c                                 the indices of which are in solc(isolc)
         do j = 1, isolc
            cslut(j) = .false.
            cslvt(j) = .false.
         end do 

         do i = 1, npt
            do j = 1, isolc

               if (lc(i,solc(j)).eq.0d0) cycle

               if (solvnt(i)) then 
                  cslvt(j) = .true.
               else 
                  cslut(j) = .true. 
               end if 

            end do
         end do

         do j = 1, isolc

            if (cslvt(j).and..not.cslut(j)) then 
c                                 a component is present only in the solvent
c                                 iteration will become unstable
               abort = .true.

               write (n13,'(i4,1x,4(g14.6,1x),a)') 1000+solc(j), 
     *                                             p0, dz, t, p,
     *                                'disolved_non-solvent_component'
               exit

            end if
         end do 

      end if 
c                                 look for a general solution if npt > icp
      if (npt.ge.hcp) then 

         do i = 1, npt

            do j = 1, hcp
               comp(i,j) = lc(i,j)
            end do 

            mu(i) = lg(i)

         end do

         ier = 0 

         if (npt.gt.hcp) then
c                                try filtering out zero mode phases
            lcp = 0

            do i = 1, npt

               if (amt(i).lt.zero) cycle

               lcp = lcp + 1

               mu(lcp) = mu(i)

               do j = 1, jcp 
                  comp(lcp,j) = comp(i,j)
               end do 

            end do

            if (lcp.ne.hcp) ier = 2

         end if 

         if (ier.eq.0) call factor (comp,hcp,ipvt,ier)

         if (ier.eq.0) call subst (comp,ipvt,hcp,mu,ier)

      end if 
c                                 ier ne 0 => look for degenerate solution
      if (idegen.gt.0.and.ier.ne.0) then

         ier = 0 
         kcp = 0

         do i = 1, npt
c                                 check if the phase contains a degenerate component
            bad = .false.

            do j = 1, idegen
               if (lc(i,idg(j)).lt.zero) cycle
c                                 reject
               bad = .true.
               exit
            end do

            if (bad) cycle
c                                 load the phase
            kcp = kcp + 1

            inp(kcp) = i

            do j = 1, jcp 
               comp(kcp,j) = lc(i,jin(j))
            end do 

            mu(kcp) = lg(i)

         end do

         if (kcp.gt.jcp) then
c                                 over-determined, try eliminating
c                                 phases present in zero amount
            lcp = 0

            do i = 1, kcp

               if (amt(inp(i)).lt.zero) cycle

               lcp = lcp + 1

               mu(lcp) = mu(i)

               do j = 1, jcp 
                  comp(lcp,j) =  comp(i,j)
               end do 

            end do

            if (lcp.ne.jcp) ier = 2

         else if (kcp.lt.jcp) then

            ier = 3

         end if

         if (ier.eq.0) then 

            call factor (comp,jcp,ipvt,ier)

            if (ier.eq.0) call subst (comp,ipvt,jcp,mu,ier)

            if (ier.eq.0) then 
c                                 load the chemical potentials 
c                                 into their correct positions
               do i = jcp, 1, -1
                  mu(jin(i)) = mu(i)
               end do 
c                                 and NaN the missing values
               do i = 1, idegen
                  mu(idg(i)) = nopt(7)
               end do

            end if

         end if 

      end if 

      if (ier.eq.0) then 
c                                 output as requested:
         mus = .true. 

         if (lopt(33)) then 
c                                 output iteration bulk G and mu's
            g = 0d0

            do i = 1, hcp

               if (isnan(mu(i))) then

                  imu(i) = 0d0
                  cycle

               else 

                  g = g + cblk(i)*mu(i)
                  imu(i) = idint(mu(i))

               end if 

            end do

            if (jter.eq.0.or.lopt(34)) 
     *                    write (*,1000) (cname(i), i = 1, hcp)
            write (*,1010) jter, g/ctotal, (imu(i), i = 1, hcp)
            if (lopt(34)) write (*,'(/)') 

         end if

      else

         mus = .false.

         if (lopt(33)) then 
c                                 output failure msg
            if (jter.eq.0.or.lopt(34)) 
     *                    write (*,1000) (cname(i), i = 1, hcp)
            write (*,1020) jter
            if (lopt(34)) write (*,'(/)') 

         end if
c                                 failed
         do i = 1, hcp
            mu(i) = nopt(7)
         end do

      end if

1000  format (/,'Iteration',4x,'G(J/mol)',7x,20(5x,a))
1010  format (3x,i2,4x,f15.4,6x,20(i9,1x))
1020  format (3x,i2,4x,'chemical potential back-calculation failed.')

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

      subroutine sublim (ids,res0) 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer i, j, ids

      double precision res0, xxnc, stinc, sum

      external stinc

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 temporary subdivision limits:
      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      integer pstot,qstot,ostg,odim,nsum
      common/ junk1 /pstot(h9),qstot(h9),ostg(h9),odim(mst,h9),nsum(h9)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c----------------------------------------------------------------------
      if (ksmod(ids).ne.20) then 
c                                normal models     
         do i = 1, ostg(ids)

            sum = 0d0 

            do j = 1, ndim(i,ids)

               sum = sum + x(i,j)

               if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                  xnc(i,j) = xncg(ids,i,j)*res0
                  xxnc = xnc(i,j)*reachg(ids)
                  xmn(i,j) = x(i,j) - xxnc
                  xmx(i,j) = x(i,j) + xxnc

               else
c                                 conformal, xnc is the number 
c                                 of intervals for 0->1 resolution
                  xnc(i,j) = xncg(ids,i,j)/res0

                  xxnc = reachg(ids)/xnc(i,j)

                  xmn(i,j) = stinc (x(i,j),-xxnc,ids,i,j)
                  xmx(i,j) = stinc (x(i,j),xxnc,ids,i,j)

               end if

               if (xmx(3,1).lt.0d0) then 
               write (*,*) 'bad news bears knocking on the door ',sum,i
               end if 
c                                 changed feb 6, 2012 from xmng/xmxg
c                                 to allow hardlimits. JADC
               if (xmn(i,j).lt.xmno(ids,i,j)) xmn(i,j) = xmno(ids,i,j)
               if (xmx(i,j).gt.xmxo(ids,i,j)) xmx(i,j) = xmxo(ids,i,j)

            end do 
         end do

      else 
c                                 charge balance model
         do i = 1, nqs1

            if (i.eq.ns) cycle

            if (imdg(i,1,ids).eq.0) then 
c                                 cartesian
               xnc(1,i) = xncg(ids,1,i)*res0
               xxnc = xnc(1,i)*reachg(ids)
               xmn(1,i) = x(1,i) - xxnc
               xmx(1,i) = x(1,i) + xxnc

            else
c                                 conformal
               xnc(1,i) = xncg(ids,1,i)/res0

               xxnc = reachg(ids)/xnc(1,i)

               xmn(1,i) = stinc (x(1,i),-xxnc,ids,1,i)
               xmx(1,i) = stinc (x(1,i), xxnc,ids,1,i)

            end if 

            if (xmn(1,i).lt.xmno(ids,1,i)) then
               xmn(1,i) = xmno(ids,1,i)
            else if (xmx(1,i).gt.xmxo(ids,1,i)) then 
               xmx(1,i) = xmxo(ids,1,i)
            end if 

         end do

      end if

      end 

      subroutine loadgx (kd,i,ids,bad) 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer i, j, kd, ids, kcoct

      logical bad, recalc

      double precision gsol1

      external gsol1

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      logical quack
      integer solc, isolc
      common/ cxt1 /solc(k5),isolc,quack(k21)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      save recalc
      data recalc/.false./
c----------------------------------------------------------------------
      jphct = jphct + 1
      recalc = .false.

      if (jphct.gt.k21) call error (58,1d0,k21,'loadgx')

      jkp(jphct) = ids
      hkp(jphct) = kd

      call prs2xy (i,ids,.true.,bad)

      if (bad) return
c                                June 8, 2018
      kcoct = jcoct + mcoor(ids)

      if (lopt(32).and.ksmod(ids).eq.39) then

         if (lopt(46)) then
c                                 set as aq_solvent_solvus:
c                                 solute free cpd
            g2(jphct) = gsol1(ids)

            call csol (ids,bad)

            if (bad) then 
               jphct = jphct - 1
               jcoct = kcoct - mcoor(ids)
               return 
            end if 

            quack(jphct) = .true.
c                                 now pad out counters for 
c                                 a solute cpd
            jphct = jphct + 1
            if (jphct.gt.k21) call error (58,1d0,k21,'resub')

            jkp(jphct) = ids
            hkp(jphct) = kd

            jcoor(jphct) = jcoct - 1

            do j = 1, mcoor(ids)
               zcoor(jcoct) = zcoor(jcoct-mcoor(ids))
               jcoct = jcoct + 1
            end do 

            kcoct = kcoct + mcoor(ids)

         end if 
c                                  solute-bearing compound
         call aqlagd (1,bad,recalc)

         if (bad) then

            jphct = jphct - 1
            jcoct = kcoct - mcoor(ids)
            return

         else
 
            quack(jphct) = .false.

         end if

      else 
c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
         g2(jphct) = gsol1(ids)
c                                 use the coordinates to compute the composition 
c                                 of the solution
         call csol (ids,bad)

         if (bad) then 
            jphct = jphct - 1
            jcoct = kcoct - mcoor(ids)
            return
         end if 

      end if 

      end

      logical function degen (index,array)
c----------------------------------------------------------------------
c check compounds for degeneracy
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, index, array

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer idegen, idg(k5), jcp, jin(k5)
      common/ cst315 /idegen, idg, jcp, jin 

      integer jphct, jpt
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct,jpt

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      degen = .false.
c                                 turn test off if aq_oxide_components, 
c                                 in principle this would allow disproportionation
c                                 for non-elemental components
      if (lopt(36)) return

      do i = 1, idegen

         if (array.eq.1) then 
            if (a(idg(i),index).gt.zero) then
               degen = .true.
               exit
            end if
         else if (array.eq.2) then 
            if (cp2(idg(i),index).gt.zero) then
               degen = .true.
               exit
            end if
         end if 

      end do

      end 