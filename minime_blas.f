      subroutine dummy
      end 

      subroutine minfrc (jter)
c-----------------------------------------------------------------------
c minimize the omega function for the independent endmember fractions
c of solution ids subject to site fraction constraints

c     number of independent endmember fractions -> nstot-1 (<m19)
c     number of independent endmember fractions -> nz (<m20)
c     closure is forced in the objective function (gsol2)

c ingsol MUST be called prior to minfrc to initialize solution/p-t
c specific properties!

c endmember gibbs energies must be computed (presumably by gall)
c prior to the call to minfxc!
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical inp, tic, zbad, toc

      integer i, nvar, iter, iwork(m22), jter, itic,
     *        istuff(10), istate(m21), idead, nclin, ntot
c DEBUG691
     *        ,iprint,mode,j,k,ind(m14)

      double precision ggrd(m19), lapz(m20,m19),gsol1, pinc,
     *                 bl(m21), bu(m21), gfinal, ppp(m19), fac,
     *                 clamda(m21),r(m19,m19),work(m23),stuff(2)
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),yt(m4),fdint,ftol,
     *                 zsite(m10,m11), pinc0,sum


      character ctol*20,cdint*20

      external gsol2, gsol1, dummy

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision wmach
      common/ cstmch /wmach(9)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      data fac,pinc0,iprint,inp,toc/1d-4,1d-2,0,.false.,.false./

      save fac,pinc0,iprint,inp,toc
c-----------------------------------------------------------------------
      yt = pa

      tic = .true.

      nclin = nz(rids)
      ntot = nstot(rids)

      if (dnu(rids).eq.0d0) then
         nvar = ntot - 1
      else 
         nvar = ntot
      end if 
c                                 finite difference increments
c                                 will be estimated at this 
c                                 coordinate, so choose a feasible 
c                                 composition
      ppp(1:nvar) = pa(1:nvar)

c        istuff(5) = 1
c                                 save the original point result
c        call gsol2 (mode,nvar,ppp,gfinal,ggrd,idead,istuff,stuff)

      istuff(5) = 0

      istuff(6) = 0

c                                 initialize bounds
      if (boundd(rids)) then 
c                                 the endmember fractions are bounded
         bu(1:nvar) = 1d0
         bl(1:nvar) = 0d0
         if (.not.lorder(rids)) nclin = 0

      else 
c                                 the model has site fractions
         bu(1:nvar) = 1d0
         bl(1:nvar) = -1d0

      end if 
c                                 load the local constraints 
c                                 from the global arrays
      lapz(1:nclin,1:nvar) = apz(rids,1:nclin,1:nvar)

      bl(nvar+1:nvar+nclin) = zl(rids,1:nclin)
      bu(nvar+1:nvar+nclin) = zu(rids,1:nclin)

      if (nvar.eq.ntot) then
c                                 closure for non-equimolar ordering
         nclin = nclin + 1
         bl(nvar+nclin) = 1d0
         bu(nvar+nclin) = 1d0
         lapz(nclin,1:nvar) = 1d0

      else if (nclin.eq.0) then 
c                                 closure for molecular models
         nclin = nclin + 1
         bl(nvar+nclin) = 0d0
         bu(nvar+nclin) = 1d0
         lapz(nclin,1:nvar) = 1d0

      end if


      itic = 0

10    idead = -1
c                                 obj call counter
      istuff(3) = 0
c                                 saved obj value counter
      istuff(4) = 0

      if (inp) then

         write (*,*) 'ftol,fdint'
         read (*,*) ftol, fdint
         write (ctol,'(g14.7)') ftol
         write (cdint,'(g14.7)') fdint
         CALL E04UEF ('optimality tolerance = '//ctol)
         CALL E04UEF ('difference interval = '//cdint)
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      else

         CALL E04UEF ('nolist')
c        CALL E04UEF ('optimality tolerance =  1d-4')
c        CALL E04UEF ('difference interval = 1d-3')
c        write (ctol,'(i4)') iprint
c        CALL E04UEF ('print level = '//ctol)

c                                 in NLPSOL:
c                                 EPSRF is function precision
c                                 CTOL  is feasibility tolerance
c                                 FTOL  is optimality tolerance
c                                 none of these are allowed to go below epsmch, if so they are 
c                                 reset to their defaults in terms of epsmch, this leads to 
c                                 the possibility that decreasing fac increases the result...
c                                 to stop this behavior modify
         CALL E04UEF ('verify level 0')
         write (ctol,'(g14.7)') (wmach(1)*fac)
         CALL E04UEF ('function precision = '//ctol)
         write (ctol,'(g14.7)') (wmach(1)*fac)**(0.8)
         CALL E04UEF ('optimality tolerance = '//ctol)
         write (ctol,'(g14.7)')  (wmach(1)*fac)**(0.5)
         CALL E04UEF ('feasibility tolerance = '//ctol)
c step limit < nopt(5) leads to bad results, coincidence?
         write (ctol,'(g14.7)') 0.5
         CALL E04UEF ('step limit = '//ctol)
c low values -> more accurate search -> more function calls
c                              0.05-.4 seem best
         CALL E04UEF ('linesearch tolerance = 0.225')
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      end if

      if (deriv(rids)) then

         if (itic.le.1) CALL E04UEF ('derivative level = 3')

         if (itic.eq.1) then
            CALL E04UEF ('verify level 1')
c           CALL E04UEF ('print level = 10')
         else if (itic.eq.2) then
            CALL E04UEF ('verify level 0')
            CALL E04UEF ('derivative level = 0')
c           CALL E04UEF ('print level = 10')
         end if 


      else

         CALL E04UEF ('verify level 0')
         CALL E04UEF ('derivative level = 0')

      end if

      call nlpsol (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol2,iter,
     *            istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,m22,work,
     *            m23,istuff,stuff,idead,iprint)

      if (iter.eq.0.and.idead.eq.0.and.itic.le.1.and.deriv(rids)) then

         pa = yt
         itic = itic + 1
         if (itic.eq.2) istuff(6) = 1
         goto 10

      else if (idead.ne.0) then 

         write (*,*) 'woana woaba, wanka?'

      else if (jter.lt.400) then 

         if (iter.eq.0) then 

            write (*,*) 'zapra',itic

         end if 


         if (toc) then 

         yt = pa

         do j = 1, 2

            pinc = pinc0/4**(j-1)

         ind = 0


         do i = 1, 2**ntot

c           pa = yt/pinc

            call binind (ind,ntot)

            sum = 0d0

            pa = yt

            do k = 1, ntot
               pa(k) = pa(k)*(1d0 + ind(k)*pinc)
               sum = sum + pa(k)
            end do

            pa = pa/sum

c           write (*,'(10(1x,i1))') ind(1:ntot)

c           pa(i) = pa(i) + (1d0 - 1d0/pinc)

c           write (*,*) rids
c           write (*,1000) pa(1:ntot)
c           write (*,1000) yt(1:ntot)
c           write (*,1000) (pa(k)-yt(k),k=1,ntot)

1001  format (i5,1x,g12.6,12(1x,f7.4))
1000  format (18x,12(1x,f7.4))


            if (zbad(pa,rids,zsite,fname(rids),.false.,fname(rids))) 
     *         then

c                write (*,*) 'oik0'

               cycle 

            end if

            call makepp (rids)
c                                 if logical arg = T use implicit ordering
            gfinal = gsol1 (rids,.true.)
c                                 get the bulk composition from pp
            call getscp (rcp,rsum,rids,rids,.false.)
c                                 increment the counter
            call savrpc (gfinal,jphct)

         end do

         end do

         else

            yt = pa

            do j = 1, 2

            pinc = 1d0 + pinc0/2**(j-1)

            do i = 1, ntot

               pa = yt/pinc
 
               pa(i) = pa(i) + (1d0 - 1d0/pinc)

            if (zbad(pa,rids,zsite,fname(rids),.false.,fname(rids))) 
     *         then

c                write (*,*) 'oik0'

               cycle 

            end if

            call makepp (rids)
c                                 if logical arg = T use implicit ordering
            gfinal = gsol1 (rids,.true.)
c                                 get the bulk composition from pp
            call getscp (rcp,rsum,rids,rids,.false.)
c                                 increment the counter
            call savrpc (gfinal,jphct)

         end do

         end do

         end if

      end if


      end

      subroutine binind (ind,ntot)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      implicit none

      integer ind(*), i, j, ntot
c                                 find the highest non-zero index
      do j = ntot, 1, -1
         if (ind(j).eq.1) exit
      end do
c                                 find the first zero index
      do i = 1, j + 1

         if (ind(i).eq.0) then

            if (i.gt.j) then
               ind(1:j) = 0
               ind(i) = 1
            else
               ind(i) = 1
            end if

            return

         end if

      end do

      end

      subroutine gsol2 (mode,nvar,ppp,gval,dgdp,istart,istuff,stuff)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for minfrc. can call 
c either gsol1 with order true or false, true seems to give better results
c presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, qwak, toc

      integer i, j, nvar, mode, istuff(*), istart, iwarn

      double precision ppp(*), gval, dgdp(*), stuff(*),
     *                 gsol1, g, sum1

      external gsol1

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      save iwarn, toc

      data iwarn, toc/0,.false./
c-----------------------------------------------------------------------
      sum1 = 0d0

      do i = 1, nvar
         sum1 = sum1 + ppp(i)
         pa(i) = ppp(i)
      end do

      if (nvar.lt.nstot(rids)) pa(nstot(rids)) = 1d0 - sum1

      if (ksmod(rids).eq.39) then

         do i = 1, nstot(rids)

            if (pa(i).gt.1d0.or.pa(i).lt.0d0) then

               if (pa(i).gt.1d0.and.pa(i).lt.1d0+zero) then 
                  pa(i) = 1d0
               else if (pa(i).lt.0d0.and.dabs(pa(i)).lt.zero) then
                  pa(i) = 0d0
               else
                  mode = -1
                  return
               end if

            end if

         end do

      end if

      call makepp (rids)

      if (deriv(rids).and.istuff(6).eq.0) then

         call getder (g,dgdp,rids)
c                                 get the bulk composition from pa
         call getscp (rcp,rsum,rids,rids,.false.)
c                                 convert dgdp to dg'dp
         do i = 1, nvar
            do j = 1, icp
               dgdp(i) = dgdp(i) - dcdp(j,i,rids)*mu(j)
            end do
         end do

         rkwak = .true.

      else if (ksmod(rids).eq.39.and.lopt(32)) then 
c                                 the last argument cancels recalc, in
c                                 which case i is a dummy. smo the total
c                                 species molality it is necessary for 
c                                 renormalization.
         call gaqlgd (g,rcp,rsum,rsmo,i,bad,.false.)

         qwak = .false.

         if (bad) then 
c                                 on failure revert to molecular fluid
            g = gsol1 (rids,.false.)
            call getscp (rcp,rsum,rids,rids,.false.)
            qwak = .true.

            if (iwarn.lt.11) then
               write (*,1000) fname(rids)
               call prtptx
               if (iwarn.eq.10) call warn (49,0d0,205,'MINFRC')
               iwarn = iwarn + 1
            end if

         end if

      else
c                                 if logical arg = T use implicit ordering
         g = gsol1 (rids,.false.)
c                                 get the bulk composition from pp
         call getscp (rcp,rsum,rids,rids,.false.)

         rkwak = .true.

      end if

      gval = g

      do i = 1, icp
         gval = gval - rcp(i)*mu(i)
      end do

      istuff(3) = istuff(3) + 1

      if (istuff(2).ne.0.and.(nvar.lt.nstot(rids).or.
     *    sum1.ge.one.and.sum1.le.1d0+zero).and.rsum.gt.zero) then
c                                 save the composition
         istuff(4) = istuff(4) + 1
c                                 increment the counter
         call savrpc (g,jphct)

c        if (toc) write (*,2000) gval,jphct
2000  format (g14.6,1x,i6)

      end if

1000  format (/,'**warning ver205** lagged speciation failed, ',
     *       'for ',a,'. The molecular',/,'speciation will be ',
     *       'output.',/)

      end

      subroutine savrpc (g,phct)
c-----------------------------------------------------------------------
c save a dynamic composition for the lp solver
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical ok

      integer phct, i, j

      double precision g, diff

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c-----------------------------------------------------------------------
c                                 increment the counter
      phct = phct + 1
c                                 normalize and save the composition
      cp2(1:icomp,phct) = rcp(1:icomp)/rsum
c                                 check if duplicate
      do i = 1, phct - 1

         if (jkp(i).eq.rids) then

            ok = .false.

            do j = 1, icp

               if (dcp(j,rids).eq.0d0) cycle

               if (dabs((cp2(j,phct) - cp2(j,i)) / dcp(j,rids)).gt.1d-4)
     *            then

                  ok = .true.
                  exit

               end if

            end do

            if (.not.ok) then

               phct = phct - 1

               return

            end if

         end if

      end do
c                                 the solution model pointer
      jkp(phct) = rids
c                                 the refinement point pointer
      hkp(phct) = rkds
c                                 save the normalized g
      g2(phct) = g/rsum
c                                 sum scp(1:icp)
      if (ksmod(rids).eq.39.and.lopt(32).and..not.rkwak) then 
         c2tot(phct) = rsum/rsmo
      else
         c2tot(phct) = rsum
      end if

1000  format (i5,1x,g12.6,12(1x,f7.4))
1010  format (18x,12(1x,f7.4))

      quack(phct) = rkwak
c                                 save the endmember fractions
      icoz(phct) = zcoct

      zco(zcoct+1:zcoct+nstot(rids)) = pa(1:nstot(rids))

      zcoct = zcoct + nstot(rids)

      end 

      subroutine gsol3 (mode,nvar,ppp,gval,dgdp,istart,istuff,stuff)
c-----------------------------------------------------------------------
c gsol3 - a shell to call gsol1 from minfxc, ingsol must be called
c         prior to minfxc to initialize solution specific paramters. only
c         called for equimolar explicit o/d models. non-equimolar o/d models
c         (currently melt(G,HGP)) involve non-linear constraints that are not
c         currently implemented in minfxc/nlpsol.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, istart, mode, istuff(*)

      double precision ppp(*), gval, psum, dgdp(*), stuff(*)

      double precision gsol1, omega0

      external gsol1, omega

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
c-----------------------------------------------------------------------
      jds = istuff(1)

      istuff(4) = istuff(4) + 1

      psum = 0d0

      do i = 1, nvar
         psum = psum + ppp(i)
         pa(i) = ppp(i)
      end do

      pa(nstot(jds)) = 1d0 - psum

      if (istuff(3).eq.0d0) then

         if (deriv(jds)) then 

            call getder (gval,dgdp,jds)

         else 
c                                 free energy minimization
            gval = gsol1 (jds,.false.)

         end if

      else
c                                 entropy maximization
         gval = -omega0 (jds,pa)

      end if

      end

      subroutine gsol4 (mode,nvar,ppp,gval,dgdp,istart,istuff,stuff)
c-----------------------------------------------------------------------
c gsol4 - a shell to call gsol1 from minfxc, ingsol must be called
c         prior to minfxc to initialize solution specific paramters. only
c         called for equimolar explicit o/d models. non-equimolar o/d models
c         (currently melt(G,HGP)) involve non-linear constraints that are not
c         currently implemented in minfxc/nlpsol.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical error 

      integer i, ids, nvar, istart, mode, istuff(*)

      double precision ppp(*), gval, dgdp(*), stuff(*), d2s(j3,j3)
c-----------------------------------------------------------------------
      ids = istuff(1)
c                                   on input ppp(1:nord) contains the 
c                                   proportions of the ordered species
c                                   pa(lstot+1:nstot).
c                                   -----------------------------------
c                                   initialize the proportions
      call ppp2pa (ppp,ids)

      if (istuff(3).eq.0) then 
c                                   free energy minimization
c                                   get g and dgdp
        call gderiv (ids,gval,dgdp,.true.,error)

      else 
c                                   negentropy minimization
         call sderiv (ids,gval,dgdp,d2s,.true.)

      end if

      istuff(4) = istuff(4) + 1

      end

      subroutine ppp2pa (ppp,ids)
c-----------------------------------------------------------------------
c set pa from p0a given current proportions of the ordered species in
c ppp
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ids, jd, k

      double precision ppp(*)

      logical pin
      common/ cyt2 /pin(j3)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      pa(1:nstot(ids)) = p0a(1:nstot(ids))
c                                   update pa for the change in the 
c                                   proportions of the ordered species
      do k = 1, nord(ids)

         if (.not.pin(k)) cycle

         jd = lstot(ids) + k

         call dpinc (ppp(k)-p0a(jd),k,ids,jd)

      end do

      end


      subroutine p2yx (id,bad)
c-----------------------------------------------------------------------
c converts the independent endmember fractions to 0-1 bounded barycentric 
c coordinates:

c     number of bounding vertices -> mstot (<m4)
c     number of independent fractions -> nstot (<m14)
c     number of linear constraints -> the number of independent
c        site fractions + closure (<m20)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, site, comp, clos, inv, zbad

      integer liw, lw, mvar, mcon, nvar, i, jter, iprint, iwarn,
     *        iwarn1, iwarn2

      character cit*4, ctol*14

      double precision scp(k5), tol

      parameter (mvar=m4, mcon=m20, liw=2*mvar+3, 
     *           lw=2*(mcon+1)**2+7*mvar+5*mcon)

      integer ncon, id, is(mvar+mcon), iw(liw), idead, istart

      double precision ax(mcon), clamda(mvar+mcon), wrk(lw), c(mvar),
     *                 a(mcon,mvar), bl(mvar+mcon), bu(mvar+mcon), 
     *                 gopt, sum, b(mcon)

      double precision wmach
      common/ cstmch /wmach(9)

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      double precision ayc
      common/ csty2c /ayc(h9,k5,m4)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jend
      common/ cxt23 /jend(h9,m14+2)
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      external zbad

      save iwarn, iwarn1, iwarn2

      data iwarn, iwarn1, iwarn2/3*0/
c-----------------------------------------------------------------------
      bad = .false.
      inv = .false.

      tol = 1d2*zero
c                                 primatic, need to invert to vertex fractions
      if (lstot(id).lt.mstot(id)) inv = .true.
c                                 choose constraints:
      if (lorder(id)) then
c                                 decompose to stoichiometric equivaluents
         call makepp (id)

         if (inv) then
c                                 prism
            site = .true.
            comp = .false.
c                                 explicit closure definitely helps
            clos = .true.

            if (dnu(id).ne.0d0) 
     *         call errdbg ('unanticipated prism/non-eq molar/py2x')

c                                 get the disordered p's
            call minfxc (gopt,id,.true.)

         else
c                                get sum (needed for non-eq molar case):
            sum = 0d0

            do i = 1, lstot(id)
c DEBUG691
               if (pp(i).lt.-1d-2) then 
                  write (*,*) 'wtf, p2yx 2',fname(id),' pp ',
     *                        pp(1:lstot(id))
                  bad = .true.
                  return
               end if 
               if (pp(i).lt.zero) pp(i) = 0d0
               sum = sum + pp(i)
            end do

            x(1,1,1:lstot(id)) = pp(1:lstot(id))/sum

            if (pop1(id).gt.1) 
     *         call errdbg ('houston we have a problem, p2yx 1')

         end if

      else

         if (inv) then
c                                 reciprocal and/or relict
c                                 equipartition
            comp = .true.
            clos = .false.
            site = .false.

         else

            x(1,1,1:lstot(id)) = pa(1:lstot(id))

            if (pop1(id).gt.1) 
     *         call errdbg ('houston we have a problem, p2yx 1')

         end if

      end if

      if (.not.inv) return

      nvar = mstot(id)
      ncon = 0
c                                 dummy objective function coefficients
c                                 (only 1 feasible point?)
      c(1:nvar) = 1d0
      bl(1:nvar) = 0d0
      bu(1:nvar) = 1d0

      if (site) then 
c                                 get the site fraction constraints
         call p2zind (pa,b,ncon,id)
c                                 load the fractions
         bl(nvar+1:nvar+ncon) = b(1:ncon)
         bu(nvar+1:nvar+ncon) = b(1:ncon)
c                                 load the ayz constraint matrix
         a(1:ncon,1:nvar) = ayz(id,1:ncon,1:nvar)

      end if

      if (comp) then 
c                                 load the ayc constraint matrix
         a(ncon+1:ncon+icp,1:nvar) = ayc(id,1:icp,1:nvar)
c                                 get the bulk 
         call getscp (scp,sum,id,1,.true.)
c
         bl(nvar+ncon+1:nvar+ncon+icp) = scp(1:icp)
         bu(nvar+ncon+1:nvar+ncon+icp) = scp(1:icp)
         ncon = ncon + icp

      end if

      if (clos) then 
c                                 add the closure constraint
         ncon = ncon + 1
         a(ncon,1:nvar) = 1d0
         bl(nvar+ncon) = 1d0
         bu(nvar+ncon) = 1d0

      end if
c                                 cold start
      istart = 0
      idead = -1
      iprint = 0

c     if (lopt(28)) call begtim (9)

      write (ctol,'(g14.7)') tol
      write (cit,'(i4)') l6

      call e04mhf ('nolist')
      call e04mhf ('iteration limit = '//cit)
      call e04mhf ('feasibility tolerance = '//ctol)
      call e04mhf ('print level = 0')
      call e04mhf ('cold start')
      call e04mhf ('problem type = fp')

      call lpsol (nvar,ncon,a,mcon,bl,bu,c,is,y,jter,gopt,ax,
     *            clamda,iw,liw,wrk,lw,idead,iprint)

c     if (lopt(28)) call endtim (9,.true.,'p2y inversion')

c                                 reset ldt, ldq, istart for phase eq
      istart = 0
c DEBUG691 to account for the unmodified lpsol ifail setting
      if (idead.le.3) then

         idead = 0

      else
c                                 really bad inversion result
         if (iwarn.lt.11) then

            write (*,1010) fname(id),idead

            call prtptx

            if (iwarn.eq.10) call warn (49,0d0,202,'P2YX')

            iwarn = iwarn + 1

         end if

         badinv(id,1) = badinv(id,1) + 1

         bad = .true.

         return

      end if
c                                 the inversion is generally weak, take any answer
c                                 within 10% of closure or positivity
      sum = 0d0

      do i = 1, mstot(id)

         sum = sum + y(i)

      end do

      if (sum.gt.1.1.or.sum.lt.0.9) then
c                                 closure violation
         if (iwarn1.lt.11) then

            write (*,1000) fname(id),(sum-1d0)*1d2

            call prtptx

            if (iwarn1.eq.10) call warn (49,0d0,201,'P2YX')
            
            iwarn1 = iwarn1 + 1

         end if

         bad = .true.

         badinv(id,1) = badinv(id,1) + 1

         return

      end if

      sum = 0d0

      do i = 1, mstot(id)

         if (y(i).lt.0d0) then
c                                 could do another inversion without
c                                 positivity constraint to see if the
c                                 answer really is outside the prism.
            if (y(i).lt.-0.05) bad = .true.

            if (iwarn2.lt.11.and.y(i).lt.-tol) then

                write (*,1020) i,y(i),fname(id)

                if (bad) then
                   write (*,1040)
                else
                   write (*,1030) i
                end if

                call prtptx

                if (iwarn2.eq.10) call warn (49,0d0,203,'P2YX')

                iwarn2 = iwarn2 + 1

            end if

            if (bad) then

               badinv(id,1) = badinv(id,1) + 1 

               return

            end if

            y(i) = 0d0

         else 

            sum = sum + y(i)

         end if

      end do
c                                 renormalize
      y(1:mstot(id)) = y(1:mstot(id))/sum

      badinv(id,2) = badinv(id,2) + 1
c                                 convert the y's to x's
      call sety2x (id)

1000  format (/,'**warning ver201** p2y inversion for ',a,' violates ',
     *       'closure by ',f5.1,'%, the result',/,'will not be used t',
     *       'o compute compositional ranges, large violations may ind',
     *       'icate that',/,'the compositional polyhedron for the mode',
     *       'l does not span all possible model compositions.',/)
1010  format (/,'**warning ver202** p2y inversion for ',a,' failed, ',
     *       'idead = ',i2,', the result',/,'will not be used t',
     *       'o compute compositional ranges.',/)
1020  format (/,'**warning ver203** negative vertex fraction y(',i2,
     *       ') = ',g8.1,' for ',a,'.',/,'Large negative values may ',
     *       'indicate that the compositional polyhedron for the model',
     *     /,'does not span all possible model compositions.',/)
1030  format ('y(',i2,') will be zeroed for computing compositional ',
     *       'ranges.',/)
1040  format ('The composition will not be used to compute',
     *       ' compositional ranges.',/)

      end

      subroutine xmnfxc (gfinal,ids,maxs)
c-----------------------------------------------------------------------
c optimize solution gibbs energy or configurational entropy at constant 
c composition subject to site fraction constraints.

c     number of compositional constraints -> icomp (<k5)
c     number of independent endmember fractions -> nstot-1 (<m19)
c     number of site fractions -> nz (<m20)
c     closure is forced in the objective function (gsol2) and
c        by using the complete set of site fraction.
c     requires that pp has been loaded in cxt7

c ingsol MUST be called prior to minfxc to initialize solution/p-t
c specific properties!
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical maxs, inp

      integer ids, i, j, k, nvar, iter, iwork(m22), iprint,
     *        istuff(10),istate(m21), idead, nclin, ntot

      double precision sum, ggrd(m19), b(k5),
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(2),
     *                 lapz(m20,m19)
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),xp(m14), ftol, fdint
      character*14 cdint, ctol

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision wmach
      common/ cstmch /wmach(9)

      external gsol3, dummy
c DEBUG691 minfxc
      data iprint,inp/0,.false./

      save iprint,inp
c-----------------------------------------------------------------------
c                                 create constraint matrix z
c                                 --------------------------
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 initialize bounds
      bu(1:nvar) = 1d20
      bl(1:nvar) = -1d20

      ppp(1:nvar) = pp(1:nvar)
      xp(1:ntot) = pp(1:ntot)
      inp = .false.
c                                 load the local constraints 
c                                 from the global arrays
      nclin = nz(ids)
c                                 first the site fraction constraints
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      do i = 1, nclin
         bl(nvar+i) = zl(ids,i)
         bu(nvar+i) = zu(ids,i)
      end do
c                                 get the normalized bulk composition of the solution
c     call getxcp (b,stuff(1),ids)

      nclin = nclin + icomp

      call getscp (b,sum,ids,ids,.false.)
c                                 --------------------------------
c                                 the bulk composition constraints
      do k = 1, icomp

         i = nz(ids) + k

         sum = 0d0
         do j = 1, nvar
            lapz(i,j) = apc(ids,k,j)
            sum = sum + lapz(i,j) * xp(j)
         end do

         bl(nvar+i) = b(k) - apc(ids,k,ntot)
c        if (dabs(bl(nvar+i)).lt.zero) bl(nvar+i) = 0d0
         bu(nvar+i) = bl(nvar+i)

c        write (*,*) 'sum, b ',bl(nvar+i),sum

      end do

      idead = -1

      iprint = 0
c                                 solution model index
      istuff(1) = ids
c                                 istuff(2) is set by NLP and is
c                                 irrelevant.
c                                 istuff(3) = 1 min g, 0 max entropy
      if (maxs) then 
         istuff(3) = 1
      else 
         istuff(3) = 0
      end if
c                                 obj call counter
      istuff(4) = 0

10    if (inp) then

         istuff(4) = 0

         iprint = 10

         ppp(1:nvar) = xp(1:nvar)

         write (*,*) 'ftol,fdint'
         read (*,*) ftol, fdint
         write (ctol,'(g14.7)') ftol
         write (cdint,'(g14.7)') fdint

         CALL E04UEF ('optimality tolerance = '//ctol)
         CALL E04UEF ('difference interval = '//cdint)
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      else

         iprint = 0 
         CALL E04UEF ('nolist')
         write (ctol,'(g14.7)') (wmach(1)*1d3)
         CALL E04UEF ('function precision = '//ctol)
         write (ctol,'(g14.7)') (wmach(1)*1d3)**(0.8)
         CALL E04UEF ('optimality tolerance = '//ctol)
         write (ctol,'(g14.7)') zero
         CALL E04UEF ('feasibility tolerance = '//ctol)
c step limit < nopt(5) leads to bad results, coincidence?
         write (ctol,'(g14.7)') nopt(5)/1d-1
         CALL E04UEF ('step limit = '//ctol)
c low values -> more accurate search -> more function calls
         CALL E04UEF ('linesearch tolerance = 0.1')
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      end if

      if (deriv(ids)) then

c        CALL E04UEF ('difference interval = -1')
         if (ids.eq.4) CALL E04UEF ('verify level 1')
         CALL E04UEF ('derivative level = 3')

      else

         CALL E04UEF ('difference interval = 1d-3')
         CALL E04UEF ('verify level 0')
         CALL E04UEF ('derivative level = 0')

      end if

      call nlpsol (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol3,iter,
     *            istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,m22,work,
     *            m23,istuff,stuff,idead,iprint)

      if (idead.eq.2) then 
         write (*,*) 'minfxc infeasible initial conditions'
c        call errpau
      end if

      sum = 0d0

      do i = 1, nvar
         pa(i) = ppp(i)
         if (dabs(pa(i)).lt.1d2*zero) then 
            pa(i) = 0d0
         else 
            sum = sum + ppp(i)
         end if 
      end do

      pa(i) = 1d0 - sum
c                                 you toucha pa, makepp
      call makepp (ids)

      if (inp) then 
         write (*,*) istuff(4),gfinal,istuff(1)
         goto 10
      end if

      end


      subroutine minfxc (gfinal,ids,maxs)
c-----------------------------------------------------------------------
c optimize solution gibbs energy or configurational entropy at constant 
c composition subject to site fraction constraints.

c     number of independent endmember fractions -> <j3
c     number of constraints -> < j3*j5 * 2
c     requires that pp has been loaded in cxt7

c ingsol MUST be called prior to minfxc to initialize solution/p-t
c specific properties!
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical maxs, inp, tic 

      integer ids, i, j, k, nvar, iter, iwork(m22), iprint,
     *        istuff(10),istate(m21), idead, nclin, ntot, lord

      double precision ggrd(m19), gordp0,g0,
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(2),
     *                 lapz(m20,m19),gord,gsol1
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),xp(m14), ftol, fdint
      character*14 cdint, ctol

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      logical pin
      common/ cyt2 /pin(j3)

      integer ln,lt,lid,jt,jid
      double precision lc, l0c, jc
      common/ cxt29 /lc(j6,j5,j3,h9),l0c(2,j5,j3,h9),lid(j6,j5,j3,h9),
     *               ln(j3,h9),lt(j5,j3,h9),jc(j3,j5,j3,h9),
     *               jid(j3,j5,j3,h9),jt(j5,j3,h9)

      double precision tsum
      common/ cxt31 /tsum(j5,j3)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      double precision wmach
      common/ cstmch /wmach(9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      external gsol4, gsol1, gordp0, dummy
c DEBUG691 minfxc
      data iprint,inp/0,.false./

      save iprint,inp
c-----------------------------------------------------------------------
c                                 initialize limit expressions from p0
      call p0limt (ids)
c                                 set initial p values and count the
c                                 the number of non-frustrated od
c                                 variables.
      call pinc0 (ids,lord)

      if (icase(ids).eq.0) then 
c                                 o/d reactions are independent and
c                                 pin settings from pinc0 are valid
c                                 regardless if whether p0 is fully 
c                                 disorderd
         if (lord.eq.0) then 
            gfinal = gord(ids)
            return
         end if

      else if (maxs) then
c                                 if maxs then p0 is likely partially 
c                                 ordered, the pin settings from pinc0
c                                 can't be relied upon, a routine could 
c                                 be added to check, but given that the 
c                                 maxs inversion is mostly likely to be
c                                 called for a general composition, the
c                                 lazy solution here is to keep everything
c                                 in:
         pin = .true.
         lord = nord(ids)

      else if (icase(ids).eq.1) then 
c                                 p0 for ~maxs will always correspond to
c                                 the disordered limit, in this case 
c                                 pin for uncorrelated o/d reactions can 
c                                 be relied up. currently the only partly
c                                 correlated case (Omph(GHP)) can also be
c                                 relied upon, but this may change. a special
c                                 test for this and fully correlated cases 
c                                 could be made, but since ~maxs calls are
c                                 only for backup from specis the lazy solution
c                                 is adopted here for the fully correlated case.
         pin = .true.
         lord = nord(ids)

      end if

      nvar = nord(ids)
c                                 variable bounds and local (ppp) variable
c                                 initialization
      do k = 1, nord(ids)

         if (pin(k)) then
            bu(k) = 1d0
            bl(k) = -1d0
         else
            bu(k) = pa(lstot(ids)+k)
            bl(k) = pa(lstot(ids)+k)
         end if

         ppp(k) = pa(lstot(ids)+k)

      end do
c                                constraints
      nclin = 0

      do k = 1, nord(ids)
c                                for each constraint
         do i = 1, ln(k,ids)

            nclin = nclin + 1
c                                bounds
            bu(nvar+nclin) = -tsum(i,k)
            bl(nvar+nclin) = -tsum(i,k) - l0c(2,i,k,ids)
c                                coefficients
            lapz(nclin,1:nvar) = 0d0

            do j = 1, jt(i,k,ids)

               lapz(nclin,jid(j,i,k,ids)-lstot(ids)) = jc(j,i,k,ids)

            end do

            if (lapz(nclin,k).ne.0d0) then 
               write (*,*) 'wtf',lapz(nclin,k)
            end if

            lapz(nclin,k) = -1d0

         end do

      end do

      idead = -1

      tic = .true.

      iprint = 0
c                                 solution model index
      istuff(1) = ids
c                                 istuff(2) is set by NLP and is
c                                 irrelevant.
c                                 istuff(3) = 1 min g, 0 max entropy
      if (maxs) then 
         istuff(3) = 1
      else 
         istuff(3) = 0
      end if
c                                 obj call counter
      istuff(4) = 0

10    if (inp) then

         istuff(4) = 0

         iprint = 10

         ppp(1:nvar) = xp(1:nvar)

         write (*,*) 'ftol,fdint'
         read (*,*) ftol, fdint
         write (ctol,'(g14.7)') ftol
         write (cdint,'(g14.7)') fdint

         CALL E04UEF ('optimality tolerance = '//ctol)
         CALL E04UEF ('difference interval = '//cdint)
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      else

         iprint = 0 
         
         CALL E04UEF ('nolist')
         write (ctol,'(g14.7)') (wmach(1)*1d1)
         CALL E04UEF ('function precision = '//ctol)
         write (ctol,'(g14.7)') (wmach(1)*1d1)**(0.8)
         CALL E04UEF ('optimality tolerance = '//ctol)
         write (ctol,'(g14.7)') zero
         CALL E04UEF ('feasibility tolerance = '//ctol)
c step limit < nopt(5) leads to bad results, coincidence?
         ftol = 0.2
         write (ctol,'(g14.7)') ftol
         CALL E04UEF ('step limit = '//ctol)
c low values -> more accurate search -> more function calls
         ftol = 0.2
         write (ctol,'(g14.7)') ftol
         CALL E04UEF ('linesearch tolerance = '//ctol)
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      end if

      if (deriv(ids)) then

c        CALL E04UEF ('difference interval = -1')
         if (.not.tic) CALL E04UEF ('verify level 1')
         CALL E04UEF ('derivative level = 3')

      else

         CALL E04UEF ('difference interval = 1d-3')
         CALL E04UEF ('verify level 0')
         CALL E04UEF ('derivative level = 0')

      end if

      call nlpsol (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol4,iter,
     *            istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,m22,work,
     *            m23,istuff,stuff,idead,iprint)
c                                   set pa to correspond to the final 
c                                   values in ppp.
      call ppp2pa (ppp,ids)

      if (iter.eq.0.and.idead.eq.0.and.tic.and.deriv(ids)) then 

         pa = p0a
         tic = .false.
         goto 10

      else if (.not.tic.and.iter.gt.0) then 

c        write (*,*) 'maialino josefino!',fname(ids),idead,iter

      end if 

      if (idead.eq.2) then 
         write (*,*) 'minfxc infeasible initial conditions'
      else if (idead.eq.7) then
         write (*,*) 'weak solution'
      else if (idead.eq.7) then
         write (*,*) 'bad derivatives'
      else if (idead.ne.0) then 
         write (*,*) 'sommat else bad',idead
      end if

      if (.not.maxs) then

         g0 = gordp0 (ids)

         if (gfinal.gt.g0) then 
            gfinal = g0
            pa = p0a
         end if

      end if 

      end

