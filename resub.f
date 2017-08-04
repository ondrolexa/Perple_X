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

      integer liw,lw,k,idead,inc

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),clamda(k1+k5),w(lw),oldt,oldp

      integer is(k1+k5),iw(liw)

      logical quit
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

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      save ax, x, clamda, w, is, iw
c-----------------------------------------------------------------------
      if (.not.usv) then 

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
         call rebulk (.true.)

      else
c                                 find discretization points
c                                 for refinement
         call yclos1 (clamda,x,is,jphct,quit)
c                                 returns quit if nothing to refine
         if (quit) then 
c                                 final processing, .true. indicates static
            call rebulk (.true.)

         else 
c                                 reoptimize with refinement
            call reopt (idead)
c                                 final processing, .false. indicates dynamic
            if (idead.eq.0) call rebulk (.false.) 

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

      integer liw, lw, iter, iref, i, j, id, idead, ids, jstart, inc, 
     *        opt

      logical first, wad1,wad2

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)  

      double precision  ax(k5), x(k21), clamda(k21+k5), w(lw)

      integer is(k21+k5), iw(liw)
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

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer jnd
      double precision aqg,q2,rt
      common/ cxt2 /aqg(m4),q2(m4),rt,jnd(m4)
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      iter = 1
      jphct = 0
      iref = 0
      wad1 = .false.
      wad2 = .true.
      jcoct = 1
      inc = istct - 1
      opt = npt
c                                 --------------------------------------
c                                 first iteration
      first = .true.

      do i = 1, npt

         id = jdv(i) + inc


c         if (ikp(id).eq.0) then #refine endmember change

         if (id.le.ipoint) then
c                                 the point is a true compound
            jphct = jphct + 1
c                                 jkp indicates which phase a point is associated with
            jkp(jphct) = -id
c                                 hkp indicates which refinement point
            hkp(jphct) = i
c                                 flag for h2o if lagged speciation
            if (lopt(32)) then
               if (id.eq.jnd(ns)) wad1 = .true.
            end if 

            g2(jphct) = g(id)/ctot(id)

            do j = 1, icp
               cp2(j,jphct) = cp(j,id)/ctot(id)
            end do 

         else 
c                                 the point is a pseudocompound, refine it
            call resub (i,id,ikp(id),iref,iter,first,wad1,wad2)

         end if
c                                 reset jdv in case of exit
         jdv(i) = jphct

      end do 

      if (iref.eq.0) return
 
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
            exit

         end if 
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,iter,opt)
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the xcoor array.
         call saver 

         if (iter.gt.iopt(10)) exit 

         jphct = 0
         wad1 = .false.
         wad2 = .true.
         iref = 0
         jcoct = 1

         first = .true.
c                                 generate new pseudocompounds
         do i = 1, npt

            ids = lkp(i)

            if (ids.lt.0) then 
c                                 the point is a true compound
               jphct = jphct + 1
               jkp(jphct) = ids
               hkp(jphct) = mkp(i)
               ids = -ids
               if (lopt(32)) then
                  if (ids.eq.jnd(ns)) wad1 = .true.
               end if 

               g2(jphct) = g(ids)/ctot(ids)

               do j = 1, icp
                  cp2(j,jphct) = cp(j,ids)/ctot(ids)
               end do 

            else
c                                 the point is a pseudocompound 
               call resub (mkp(i),i,ids,iref,iter,first,wad1,wad2)

            end if
c                                 reset jdv in case of exit
            jdv(i) = i 

         end do 

         if (jphct.lt.icp) then 
c                                 something weird, could recover the
c                                 last good solution instead of setting 
c                                 idead.
            call warn (99,0d0,0,'something has gone terribly wrong'//
     *                    ', not enough refinement points in reopt')

            idead = 99

            exit 

         end if 

      end do 

      end

      subroutine resub (jd,id,ids,iref,iter,first,wad1,wad2)
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
      logical bad, first, quack1, quack2, quack3, wad1, wad2, wad3

      double precision xxnc, ysum, res0

      integer i, j, k, l, m, ids, id, jd, iter, kcoct, iref, last 
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
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 temporary subdivision limits:
      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)
c                                 coordinates output by subdiv
      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

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

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ineg
      common/ cst91 /ineg(h9,m15)

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      save last
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
      res0 = nopt(24)/nopt(21)**iter

      if (ksmod(ids).ne.20) then 
c                                normal models     
         do i = 1, istg(ids)

            do j = 1, ndim(i,ids)

               xnc(i,j) = xncg(ids,i,j)*res0
               xxnc = xnc(i,j)*reachg(ids)

               if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                  xmn(i,j) = x(i,j) - xxnc
                  xmx(i,j) = x(i,j) + xxnc

               else
c                                 conformal
                  xmn(i,j) = ydinc (x(i,j),-xxnc,imdg(j,i,ids),j,i,ids)
                  xmx(i,j) = ydinc (x(i,j),xxnc,imdg(j,i,ids),j,i,ids)

               end if 
c                                 changed feb 6, 2012 from xmng/xmxg
c                                 to allow hardlimits. JADC
               if (xmn(i,j).lt.xmno(ids,i,j)) xmn(i,j) = xmno(ids,i,j)
               if (xmx(i,j).gt.xmxo(ids,i,j)) xmx(i,j) = xmxo(ids,i,j)

            end do 
         end do

      else 
c                                 charge balance model
         j = -1

         do i = 1, nqs1

            if (i.eq.ns) cycle

            xnc(1,i) = xncg(ids,1,i)*res0
            xxnc = xnc(1,i)*reachg(ids)

            if (imdg(i,1,ids).eq.0) then 
c                                 cartesian
               xmn(1,i) = x(1,i) - xxnc
               xmx(1,i) = x(1,i) + xxnc

            else
c                                 conformal
               xmn(1,i) = ydinc (x(1,i),-xxnc,imdg(i,1,ids),i,1,ids)
               xmx(1,i) = ydinc (x(1,i),xxnc,imdg(i,1,ids),i,1,ids)

            end if 

            if (xmn(1,i).lt.xmno(ids,1,i)) then
               xmn(1,i) = xmno(ids,1,i)
            else if (xmx(1,i).gt.xmxo(ids,1,i)) then 
               xmx(1,i) = xmxo(ids,1,i)
            end if 

         end do

      end if 
c                                  set solution model parameters for
c                                  gsol1, don't call if the previous
c                                  refinement point was the same solution.
      if (ids.ne.last) call ingsol (ids) 

      call subdiv (ids,.true.)

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
         kcoct = jcoct + mcoor(ids)
c                                 counter for number of non 0 or 1 compositions
         if (kcoct.gt.k20) call error (59,x(1,1),k20,'resub')

         l = (i-1)*mcoor(ids)
         m = 0
       
         if (ksmod(ids).ne.20) then 

            do j = 1, istg(ids)

               ysum = 0d0

               do k = 1, ndim(j,ids)

                  m = m + 1

                  x(j,k) = prism(l+m)
                  ysum = ysum + x(j,k)
                  zcoor(jcoct) = x(j,k)
 
                  if (x(j,k).lt.xmno(ids,j,k).and.
     *                x(j,k).gt.xmxo(ids,j,k)) then 
c                                 the composition is out of range
                     jphct = jphct - 1
                     jcoct = kcoct - mcoor(ids)
                     goto 10

                  end if

                  jcoct = jcoct + 1

               end do

               x(j,ispg(ids,j)) = 1d0 - ysum

            end do 

         else 
c                                 charge balance models: a wierd shuffle to put
c                                 the first nqs - 1 species in zcoor
            ysum = 0d0

            do k = 1, nqs

               if (k.eq.ns) cycle 

               m = m + 1

               x(1,k) = prism(l+m)
               ysum = ysum + x(1,k)

               if (k.eq.nqs) cycle 

               zcoor(kcoct-nqs+k) = x(1,k)
 
               if (x(1,k).lt.xmno(ids,1,k).and.
     *             x(1,k).gt.xmxo(ids,1,k)) then 
c                                 the composition is out of range
                  jphct = jphct - 1
                  jcoct = kcoct - mcoor(ids)
                  goto 10

               end if

            end do

            x(1,ns) = 1d0 - ysum

            zcoor(kcoct-qn) = x(1,ns)
            jcoct = jcoct + nqs1

         end if 

         call xtoy (ids,bad)

         if (bad) then 

            jphct = jphct - 1
            jcoct = kcoct - mcoor(ids)
            cycle

         else if (ksmod(ids).eq.5) then
c                                 this is an el cheapo filter for redundant
c                                 compositions, a better method would be to
c                                 do the subdivision properly.
            do j = 1, ndep(ids)

               if (y(knsp(lstot(ids)+j,ids)).gt.0d0.and.
     *             y(knsp(lstot(ids)+j,ids)).le.y(ineg(ids,j))) then
c                                 reject composition 
                  jphct = jphct - 1
                  jcoct = kcoct - mcoor(ids)
                  bad = .true. 
                  exit 

               end if 

            end do 
      
            if (bad) cycle 

         end if

         if (lopt(32).and.ksmod(ids).eq.39) then

            quack1 = .false.
            quack2 = .true.

            if (quack1) then 
c                                 solute free cpd
            g2(jphct) = gsol1(ids)
            call csol (ids)
c                                 now pad out counters for 
c                                 a solute cpd
            jphct = jphct + 1
            jkp(jphct) = ids
            hkp(jphct) = jd
            jcoor(jphct) = jcoct - 1

            do j = 1, mcoor(ids)
               zcoor(jcoct) = zcoor(jcoct-mcoor(ids))
               jcoct = jcoct + 1
            end do 

            kcoct = kcoct + mcoor(ids)

            end if 

            if (quack2) then 

            call aqlagd (1,bad,.false.)

            if (bad) then

               jphct = jphct - 1
               jcoct = kcoct - mcoor(ids)
               cycle

            end if  

            end if

            if (wad1.and.wad2) then 
c                                 make water, ha ha
               wad2 = .false.

               jphct = jphct + 1
               jkp(jphct) = ids
               hkp(jphct) = jd
               jcoor(jphct) = jcoct - 1

               do j = 1, mcoor(ids)
                  zcoor(jcoct) = 0d0
                  y(j) = 0d0
                  jcoct = jcoct + 1
               end do 

               y(ns) = 1d0

               kcoct = kcoct + mcoor(ids)

               call aqlagd (1,bad,.false.)

               if (bad) then

                  jphct = jphct - 1
                  jcoct = kcoct - mcoor(ids)
                  cycle

               end if

            end if 

         else 
c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
            g2(jphct) = gsol1(ids)
c                                 use the coordinates to compute the composition 
c                                 of the solution
            call csol (ids)

         end if 

         iref = iref + 1

10    continue

      first = .false.

      last = ids 

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

      integer i, j, k, id

      double precision ctot2

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
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
      if (ctot2.ne.0d0) then 
c                                  normalize the composition and free energy
         g2(jphct) = g2(jphct)/ctot2

         do j = 1, icp 
            cp2(j,jphct) = cp2(j,jphct)/ctot2
         end do

      else 

         g2(jphct) = 1d4*p
         cp(1,jphct) = 1d0

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

      integer idead, iwarn91, iwarn42, iwarn90

      character*6 char     

      double precision c

      save iwarn91, iwarn42, iwarn90

      data iwarn91, iwarn42, iwarn90/0,0,0/
c----------------------------------------------------------------------
c                                             look for errors                                            
      if (idead.eq.2.or.idead.gt.4.and.iwarn91.lt.6) then 
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

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt
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

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)
c----------------------------------------------------------------------
      jcoor = lcoor(id)

      do i = 1, istg(ids)

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

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
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

      logical solvs1, check, bad
c                                 -------------------------------------
c                                 local variables
      integer idsol(k5),kdsol(k5,k5),ids,isite,xidsol,xkdsol,irep,
     *        i,j,jdsol(k5,k5),jd,k,l,nkp(k5),xjdsol(k5)

      double precision bsol(k5,k5),cpnew(k5,k5),xx,xb(k5),
     *                 bnew(k5),xnew(k5,mst,msp),ncaq(k5,l10)
c                                 -------------------------------------
c                                 global variables:
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                  x-coordinates for the final solution
      integer kd, na1, na2, na3, na4
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,na4,kd

      integer nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa
      common/ cst337 /nq,nn,ns,ns1,sn1,nqs,nqs1,sn,qn,nq1,nsa

      integer iaq, aqst, aqct
      character aqnam*8
      double precision aqcp, aqtot
      common/ cst336 /aqcp(k0,l9),aqtot(l9),aqnam(l9),iaq(l9),aqst,aqct

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

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0
c                                solvus_tolerance
      soltol = nopt(8)

      do 30 i = 1, ntot
          
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = amt(i)
            kdsol(ntot,ncpd) = nkp(i)
            jdsol(ntot,ncpd) = i

         else 

            if (lopt(32).and.ksmod(nkp(i)).eq.39) then 
c                                 get lagged speciation
c                                 loaded into caq(i,1:ns+aqct)
               do k = 1, ns
                  y(k) = x3(i,1,k)
               end do 

               call aqlagd (i,bad,.true.)

            end if

            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (kdsol(j,1).eq.nkp(i)) then 
c                                 if match check for a solvus
                  if (.not.solvs1(i,jdsol(j,idsol(j)),nkp(i))) then
c                                 the pseudocompound matches a solution
c                                 found earlier.
                     idsol(j) = idsol(j) + 1
                     bsol(j,idsol(j)) = amt(i)
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
            bsol(np,1) = amt(i)

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
         ids = kdsol(i,1)
         isite = istg(ids)

         do j = 1, icomp
            cpnew(j,i) = 0d0
         end do 

         do j = 1, isite
            do k = 1, ispg(ids,j)
               xnew(i,j,k) = 0d0
            end do 
         end do

         if (lopt(32).and.ksmod(ids).eq.39) then 
c                               lagged speciation
            do k = 1, na4
               ncaq(i,k) = 0d0
            end do 

         end if 

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + amt(jdsol(i,j))
         end do 

         do j = 1, idsol(i)

            jd = jdsol(i,j)
c                                conditional in case zero mode
c                                is off:
            if (bnew(i).gt.0d0) then 

               xx =  amt(jd)/bnew(i)

            else 

               xx = 1d0

            end if
c                                save the new compositions
            do k = 1, icomp
               cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
            end do

            do k = 1, isite
               do l = 1, ispg(ids,k)
                  xnew(i,k,l) = xnew(i,k,l) + xx*x3(jd,k,l)
               end do
            end do

            if (lopt(32).and.ksmod(ids).eq.39) then 
c                               lagged speciation (1:nsa), ionic strength (na1), total
c                               molality (na2), solvent mass (na3), err_log_kw (na4)
               do k = 1, na4
                  ncaq(i,k) = ncaq(i,k) + xx*caq(jd,k)
               end do

            end if

         end do

      end do
c                                now reform the arrays kdv and b
      do i = 1, np

         amt(i) = bnew(i)
         kkp(i) = kdsol(i,1)
         ids = kkp(i)

         do j = 1, icomp
            cp3(j,i) = cpnew(j,i)
         end do

         do j = 1, istg(ids)
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
            do k = 1, na4
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
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)

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
      do i = 1, istg(ids)

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
                     if (xmng(ids,i,j).lt.0d0) xmng(ids,i,j) = 0d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then
c                                 assymmetric stretching towards xmin, this
c                                 probably isn't even correct for imd = 1, but
c                                 only a lunatic would use this not going from 
c                                 zero. 
                     yint(1,j,i,ids) = yint(1,j,i,ids) - xncg(ids,i,j)
     *                                   *(xmxg(ids,i,j)-xmng(ids,i,j))
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
     *            (x(i,j).ge.xmxg(ids,i,j).and..not.lopt(3))) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmxg(ids,i,j) = xmxg(ids,i,j) + xncg(ids,i,j)
                     if (xmxg(ids,i,j).gt.1d0) xmxg(ids,i,j) = 1d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then
c                                 assymmetric stretching, this is probably
c                                 only correct for imd = 1. 

c                                 by not changing
c                                 the relative increment, expanding the limits
c                                 lowers the resolution?

c                                 for aq electrolytes i had this, because i was
c                                 worried it was expanding the limit and decreasing 
c                                 resolution. this needs to be checked out. 5/17

c                     yint(2,j,i,ids) = yint(2,j,i,ids) + xncg(ids,i,j)
c     *                                   *(xmxg(ids,i,j)-xmng(ids,i,j))
c                     xncg(ids,i,j) = 1d0/(1d0/xncg(ids,i,j) + 1d0)
c                                  revert to the original form. 5/17
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

      subroutine sorter (kdbulk,ico,jco,output)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,m,kdbulk,ico,jco,ids,ioct,inct

      logical output, reord, match, nomtch, ok 

      double precision cpt(k5,k5),xt(k5,mst,msp),bt(k5)
c                                 x-coordinates for the final solution
      integer kd, na1, na2, na3, na4
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,na4,kd

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

                  if (kkp(k).eq.idasls(j,i)) then
c                                 load temporary array
                  bt(j) = amt(k)

                     if (kkp(k).gt.0) then 

                        do l = 1, icomp
                           cpt(l,j) = cp3(l,k)
                        end do

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
      integer kd, na1, na2, na3, na4
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,na4,kd
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
         write (n5,1010) ((x3(i,j,k),k=1,ispg(ids,j)),j=1,istg(ids))
c                                lagged speciation
         if (ksmod(ids).eq.39.and.lopt(32)) write (n5,1010) 
     *                                            (caq(i,j),j=1,na4)

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

      integer jphct, i, j, k, is(*), idsol(k5), kdv(h9), nsol, 
     *        mpt, iam, id, inc, jdsol(k5,k5), ldv1, ldv2,
     *        kdsol(k5), max, mcpd

      external ffirst

      logical solvus, quit

      double precision clamda(*), x(*),  slam(h9), clam

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

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus
c----------------------------------------------------------------------
      npt = 0 
      nsol = 0
      inc = istct - 1
      quit = .true.
c                                 solvus_tolerance_II, 0.25
      soltol = nopt(25)

      do i = 1, jphct

         if (is(i).ne.1) then 
c                                 make a list of found phases:
            id = i + inc
c                                 currently endmember compositions are not 
c                                 refined (this is probably a mistake, but 
c                                 seems to work fine), so use id > ipoint
c                                 to identify a refineable point

c                                 modified 3/4/2015 to refine endmember compositions. JADC
c           if (ikp(id).ne.0) then  #refine endmember change
            if (ikp(id).ne.0.and.id.gt.ipoint) then

               quit = .false.

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

         end if 

      end do 
c DEBUG DEBUG:
c                                 get mus for lagged speciation
      call getmus (1,0,.false.)

      ldv1 = 0
      ldv2 = 0 
      clam = 1d99

      do i = 1, isoct
         slam(i) = 1d99
         kdv(i) = 0 
      end do 
c                                 perp 6.6.3, make a list of metastable
c                                 phases, this list includes two compounds
c                                 and the least metastable composition of
c                                 each solution.      
      do 20 i = 1, jphct
c DEBUG why was this here? added ~6.7.6, removed april 21, 2017
c i think clamda(i).lt.0 allows degenerate compositions (and probably 
c therefore the 6.7.6 version may be better, on the bright side with
c it removed the solution composition gets refined (if endmember comps are
c being allowed, see ldsol code); restored again april 2017.
         if (is(i).ne.1.or.clamda(i).lt.wmach(3)) cycle 
c        if (is(i).ne.1) cycle 

         id = i + inc 
         iam = ikp(id)
c                                modified to allow endmembers, Mar 4, 2015. JADC
c        if (iam.ne.0.and.id.gt.0*ipoint) then   #refine endmember change
         if (iam.ne.0.and.id.gt.ipoint) then  

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

         else 
c                                a compound, save lowest 2, changed from 1
c                                march 4, 2015, JADC
            if (clamda(i).lt.clam) then 
c                                put the old min into ldv2
               ldv2 = ldv1
c                                and save the current in ldv1
               ldv1 = i

               clam = clamda(i)

            end if 
 
         end if 

20    continue
c                                 load the metastable compounds directly
c                                 into jdv
      if (ldv2.ne.0) then 
         npt = npt + 2
         jdv(npt-1) = ldv2
         jdv(npt) = ldv1
         mcpd = 2
      else if (ldv1.ne.0) then 
         npt = npt + 1
         jdv(npt) = ldv1
         mcpd = 1
      else 
         mcpd = 0 
      end if 
c                                 load the metastable solutions into kdv
      mpt = 0 

      do i = 1, isoct

         if (kdv(i).eq.0) cycle
         mpt = mpt + 1
         kdv(mpt) = kdv(i)
         slam(mpt) = slam(i)

      end do 

      if (mpt.le.iopt(12)) then 
c                                 less metastable refinement points than
c                                 iopt(12)
         max = mpt

      else 
c                                 sort the metastable points to
c                                 find the most stable iopt(12) points
         max = iopt(12)

         call ffirst (slam,kdv,1,mpt,max,h9,ffirst)

      end if 
 
      do i = 1, max

         jdv(npt+i) = kdv(i)
c                                 a metastable solution to be refined
c        if (kdv(i)+inc.gt.ipoint) quit = .false.
         if (ikp(kdv(i)+inc).ne.0) quit = .false.
      end do

      if (quit) then 
c                                 zero mode filter and 
c                                 save amounts for final processing
         mpt = npt - mcpd 
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
c initialize arrays and constants for lp minimizarion
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
c-----------------------------------------------------------------------
      inc = istct - 1

      tloop = 40
      ploop = 40
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

      end 

      subroutine yclos0 (x,is,jphct)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,jphct,is(*)

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

      call getmus (1,0,.true.)

      end 

      subroutine setx3 (ind,id,ids)
c----------------------------------------------------------------------
c subroutine to recover prismatic solution compositions (x(i,j))
c from the xco array loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, kcoor, ind
c                                 x coordinate description
      integer istg, ispg, imlt, imdg

      double precision xt

      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 stored x coordinate
      double precision xco
      integer ico,jco
      common/ cxt10 /xco(k18),ico(k1),jco(k1)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)
c                                  x-coordinates for the final solution
      integer kd, na1, na2, na3, na4
      double precision x3, caq
      common/ cxt16 /x3(k5,mst,msp),caq(k5,l10),na1,na2,na3,na4,kd

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

            do i = 1, istg(ids)

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

      subroutine yclos2 (clamda,x,is,iter,opt)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      external ffirst

      integer i, is(*), id, jmin(k19), kmin(k19), opt, kpt, mpt, iter, 
     *        tic

      double precision clamda(*), clam(k19), x(*)

      logical stable(k19)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision wmach(9)
      common /ax02za/wmach

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

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

      npt = 0 
      kpt = 0

      do i = 1, jphct
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 check the stability of all points 
c DEBUG                           5/9/2017 removed the x(i)>0 conditional
c                                 to get chemical potentials. unfortunately
c                                 i forgot why removing it is a bad idea.
c                                 so we'll try it one more time (july 3, 2017): 
c        if (is(i).ne.1.and.x(i).gt.0d0) then
         if (is(i).ne.1) then  
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            stable(id) = .true.

         else if (clamda(i).lt.clam(id)) then
c DEBUG wish i wrote was this is about, but it may be in yclos0/1
            if (clamda(i).lt.wmach(3)) cycle 

            if (jkp(id).gt.0) then
c                                 it's a solution, keep the  
c                                 least metastable point
               jmin(id) = i
               clam(id) = clamda(i)
  
            else 
c                                 it's a compound, keep all 
c                                 as they hardly cost anything. 
               kpt = kpt + 1
               kmin(kpt) = i

            end if 

         end if 

      end do

c      write (*,*) 'npt opt kpt', npt, opt, kpt

      if (npt.gt.icp) then 

c         write (*,*) 'too many', npt, icp, iter
c         do i = 1, npt
c            write (*,*) jdv(i),hkp(jdv(i)),jkp(idv(i)),
c     *                  x(jdv(i)),clamda(jdv(i)),is(jdv(i))
c         end do 
c                                 july 10, 2017. 
c                                 too many stable points, see if there
c                                 is a zero-mode stable phase. attempt
c                                 to recover pre-july 3, 2017 behavior;
c                                 actually there is only a small chance 
c                                 that npt > icp cause's getmus to fail, 
c                                 and if it did this could be solved by 
c                                 reordering the stable phases until a 
c                                 non-degenerate case is found... comment again
         mpt = npt
         npt = 0 

         do i = 1, mpt
            stable(i) = .false.
         end do

         do i = 1, mpt

            if (x(jdv(i)).gt.0d0) then

               npt = npt + 1
               jdv(npt) = jdv(i)
               stable(hkp(jdv(i))) = .true.

            end if 
c                                  because i expect that this is occurring 
c                                  when the endmember is identical to the 
c                                  solution composition, don't save the phase
c                                  in the list of metastable phases (as was 
c                                  done pre-july 3, 2017.

         end do 

c        if (npt.lt.icp) write (*,*) 'now too few', npt, icp, iter

      end if 
c                                 get mu's for lagged speciation
      call getmus (iter,iter-1,.false.)  

c      write (*,*) 'npt opt kpt', npt, opt, kpt     

      if (iter.le.iopt(10)) then
c                                 if not done iterating, add the metastable
c                                 phases, first the compounds:
         do i = 1, kpt
            npt = npt + 1 
            jdv(npt) = kmin(i)
         end do 
c                                 make a list of the solutions
         kpt = 0 

         do i = 1, opt

            if (jmin(i).eq.0) then
               cycle
            else if (stable(hkp(jmin(i)))) then
c                                 contrary to what you might expect, this
c                                 definitely improves quality, presumably because
c                                 the metastable point will always be the closest
c                                 composition to the stable point and the resulting
c                                 overlap fogs the lp. 
               cycle 
            end if 

            kpt = kpt + 1
            jmin(kpt) = jmin(i)
            clam(kpt) = clam(i)

         end do
c DEBUG DEBUG 
c should be gt 0
         if (kpt.gt.iopt(31)) then 
c                                 sort if too many

            kpt = iopt(31)

            call ffirst (clam,jmin,1,kpt,iopt(31),k19,ffirst)

         end if 

         do i = 1, kpt

           npt = npt + 1
           jdv(npt) = jmin(i)

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
            else if (lopt(13).and.x(jdv(i)).lt.-nopt(9)
     *                             .and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

      end if 

c      write (*,*) 'npt opt ', npt, opt

      end

      subroutine rebulk (static)
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

      logical static, bad

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
      call avrger
      
1000  format (/,'**error ver901** solutions not allowed in saturated ',
     *   'component composition space',/,'in adaptive optimization ',
     *   'calculations, this limitation can be removed upon request.')

      end

      subroutine getmus (iter,jter,quit) 
c----------------------------------------------------------------------
c iter is a flag indicating where the compositions are and is sort of 
c      related to the iteration count during optimization.
c jter is the true iteration count.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical quit 

      integer i, j, id, ier, ipvt(k8), iter, jter, ic(k5), jc(k5), kcp

      double precision comp(k8,k8), g

      character cname*5
      common/ csta4  /cname(k5)

      integer hcp, idv
      common/ cst52  /hcp,idv(k7) 

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

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
c----------------------------------------------------------------------

      mus = .false.

      if (.not.lopt(32)) then
 
         if ((iter.lt.iopt(10).or.quit).and.jpot.ne.0) return

      end if 
         
      if (npt.eq.hcp) then 

         do i = 1, hcp

            id = jdv(i) 

            if (iter.gt.1) then

               do j = 1, hcp 
                  comp(i,j) = cp2(j,id)
               end do 

               mu(i) = g2(id)

            else 

               do j = 1, hcp 
                  comp(i,j) = a(j,id)
               end do 

               mu(i) = c(id)

            end if 


         end do

         call factor (comp,hcp,ipvt,ier)

         if (ier.eq.0) call subst (comp,ipvt,hcp,mu,ier)

      else 
c                                 npt < hcp, look for a zero component
         kcp = 0 

         do i = 1, hcp 

            if (cblk(i).gt.nopt(11)) then 

               kcp = kcp + 1
               ic(kcp) = i
               jc(i) = kcp

            else 

               jc(i) = 0

            end if

         end do

         if (npt.eq.kcp) then
c                                 try to salvage by elimating the zero
c                                 component 
            do i = 1, kcp

               id = jdv(i) 

               if (iter.gt.1) then

                  do j = 1, kcp 
                     comp(i,j) = cp2(ic(j),id)
                  end do 

                  mu(i) = g2(id)

               else 

                  do j = 1, kcp 
                     comp(i,j) = a(ic(j),id)
                  end do 

                  mu(i) = c(id)

               end if 

            end do

            call factor (comp,kcp,ipvt,ier)

            if (ier.eq.0) call subst (comp,ipvt,kcp,mu,ier)

         else
c                                 hopeless, well it could be done, but i won't now. 
            ier = 1

         end if

      end if 

      if (ier.eq.0) then 

         mus = .true. 

         if (npt.lt.hcp) then
c                                 degenerate, shift the mu's back to 
c                                 the original positions
            do i = hcp, 1, -1

               if (jc(i).eq.0) then 
                  mu(i) = nopt(7)
               else 
                  mu(i) = mu(jc(i))
               end if 

            end do

         end if 

         if (lopt(33)) then 
c                                 output iteration bulk G and mu's
            g = 0d0 

            do i = 1, hcp

               if (isnan(mu(i))) cycle

               g = g + cblk(i)*mu(i)
       
            end do 

            if (jter.eq.0) write (*,1000) (cname(i), i = 1, hcp)
            write (*,1010) jter, g/ctotal, (idint(mu(i)), i = 1, hcp)

         end if

      else 

         if (lopt(33)) then 
c                                 output failure msg
            if (jter.eq.0) write (*,1000) (cname(i), i = 1, hcp)
            write (*,1020) jter

         end if
c                                 failed
         do i = 1, hcp
            mu(i) = nopt(7)
         end do

      end if

1000  format (/,'Iteration',4x,'G(J/mol)',7x,20(4x,a))
1010  format (3x,i2,4x,f15.4,6x,20(i8,1x))
1020  format (3x,i2,4x,'chemical potential back-calculation failed.')

      end