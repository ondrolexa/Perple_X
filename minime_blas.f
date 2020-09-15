      subroutine dummy
      end 

      subroutine minfrc (ids,kds)
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

      logical inp, tick

      integer ids, kds, nvar, iter, iwork(m22),
     *        istuff(10), istate(m21), idead, nclin, ntot
c DEBUG691
     *         , i, j, iprint,mode

      double precision ggrd(m19), lapz(m20,m19),
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(2)
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),yt(m4),zp,zt(m10,m11),
     *                 ftol,fdint

      character ctol*20,cdint*20

      external gsol2, dummy

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

      data iprint,inp/0,.false./

      save iprint,inp
c-----------------------------------------------------------------------
      zp = 0d0
      ftol = 1d-1

      do i = 1, nstot(ids)
         if (pa(i).eq.0d0) then 
            pa(i) = ftol
         end if
         zp = zp + pa(i)
      end do

      pa(1:nstot(ids)) = pa(1:nstot(ids)) / zp

      yt = pa

10    nclin = nz(ids)
      ntot = nstot(ids)

      if (dnu(ids).eq.0d0) then
         nvar = ntot - 1
      else 
         nvar = ntot
      end if 
c                                 finite difference increments
c                                 will be estimated at this 
c                                 coordinate, so choose a feasible 
c                                 composition
c DEBUGXXX
c     pa(1:ntot) = 0.25d0
      ppp(1:nvar) = pa(1:nvar)
c                                 initialize bounds
      if (nclin.gt.0) then 
c                                 the model has site fractions
         bu(1:nvar) = 1d20
         bl(1:nvar) = -1d20

      else 
c                                 a molecular model, the endmember fractions 
c                                 are bounded
         bu(1:nvar) = 1d0
         bl(1:nvar) = 0d0

      end if 
c                                 load the local constraints 
c                                 from the global arrays
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      bl(nvar+1:nvar+nclin) = zl(ids,1:nclin)
      bu(nvar+1:nvar+nclin) = zu(ids,1:nclin)

      tick = .false.

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
         tick = .false.

      end if

      idead = -1
c                                 the solution model index
      istuff(1) = ids
c                                 print/save obj value 
c     istuff(2) set by nlpopt
c                                 obj call counter
      istuff(3) = 0
c                                 saved obj value counter
      istuff(4) = 0
c                                 refinement point index
      istuff(5) = kds

      if (inp) then

         ppp = 0

         iprint = 10

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
         if (tick.or.deriv(ids)) iprint = 0

         CALL E04UEF ('nolist')
         CALL E04UEF ('optimality tolerance =  1d-4')
         CALL E04UEF ('difference interval = 1d-3')
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      end if

      if (deriv(ids)) then

c        CALL E04UEF ('verify level 1')
         CALL E04UEF ('derivative level = 3')

      else

         CALL E04UEF ('verify level 0')
         CALL E04UEF ('derivative level = 0')

      end if

c     call nlpopt (nvar,nclin,m20,m19,lapz,bl,bu,gsol2,
c    *             iter,istate,clamda,gfinal,ggrd,r,ppp,iwork,
c    *             m22,work,m23,istuff,stuff,idead,iprint)

c     write (*,*) '*****with derivatives****',fname(ids)

      call nlpsol (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol2,iter,
     *            istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,m22,work,
     *            m23,istuff,stuff,idead,iprint)

c     call gsol2 (mode,nvar,ppp,gfinal,ggrd,idead,istuff,stuff)

c     if (deriv(ids)) then 
c        deriv(ids) = .false.
c        call gsol2 (mode,nvar,ppp,gfinal,ggrd,idead,istuff,stuff)
c        deriv(ids) = .true.
c     end if 


c     if (deriv(ids)) then 
c     write (*,*) '*****without derivatives*****',fname(ids)
c     pa = yt
c        deriv(ids) = .false.
c        CALL E04UEF ('optimality tolerance =  1d-4')
c        CALL E04UEF ('difference interval = 1d-3')
c        CALL E04UEF ('verify level 0')
c        CALL E04UEF ('derivative level = 0')
c     call nlpsol (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol2,iter,
c    *            istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,m22,work,
c    *            m23,istuff,stuff,idead,iprint)
c        deriv(ids) = .true.

c     end if




      if (inp) then 

         deriv(ids) = .false.

         pa = yt

         write (*,*) istuff(3),gfinal,kds

         goto 10

      end if

1000  format (a10,10(g14.7,1x))
1010  format (i5,1x,10(g14.7,1x))
      end

      subroutine gsol2 (mode,nvar,ppp,gval,dgdp,istart,istuff,stuff)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for minfrc. can call 
c either gsol1 with order true or false, true seems to give better results
c presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad, qwak

      integer i, j, jds, nvar, mode, istuff(*), istart, iwarn

      double precision ppp(*), gval, dgdp(*), stuff(*),
     *                 gsol1, g, sum, scp(k5), sum1, smo

      external gsol1

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      save iwarn

      data iwarn/0/
c-----------------------------------------------------------------------
      jds = istuff(1)

      sum1 = 0d0

      do i = 1, nvar
         sum1 = sum1 + ppp(i)
         pa(i) = ppp(i)
      end do

      if (nvar.lt.nstot(jds)) pa(nstot(jds)) = 1d0 - sum1

      if (ksmod(jds).eq.39) then

         do i = 1, nstot(jds)

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

      call makepp (jds)

      if (deriv(jds)) then

         call getder (g,dgdp,jds)
c                                 get the bulk composition from pp
         call getscp (scp,sum,jds,jds,.false.)
c                                 convert dgdp to dg'dp
         do i = 1, nvar
            do j = 1, icp
               dgdp(i) = dgdp(i) - dcdp(j,i,jds)*mu(j)
            end do
         end do

         qwak = .true.

      else if (ksmod(jds).eq.39.and.lopt(32)) then 
c                                 the last argument cancels recalc, in
c                                 which case i is a dummy. smo the total
c                                 species molality it is necessary for 
c                                 renormalization.
         call gaqlgd (g,scp,sum,smo,i,bad,.false.)

         qwak = .false.

         if (bad) then 
c                                 on failure revert to molecular fluid
            g = gsol1 (jds,.false.)
            call getscp (scp,sum,jds,jds,.false.)
            qwak = .true.

            if (iwarn.lt.11) then
               write (*,1000) fname(jds)
               call prtptx
               if (iwarn.eq.10) call warn (49,0d0,205,'MINFRC')
               iwarn = iwarn + 1
            end if

         end if

      else
c                                 if logical arg = T use implicit ordering
         g = gsol1 (jds,.false.)
c                                 get the bulk composition from pp
         call getscp (scp,sum,jds,jds,.false.)

         qwak = .true.

      end if

      gval = g

      do i = 1, icp
         gval = gval - scp(i)*mu(i)
      end do
c                                  normalize for appearances
c     gval = gval/sum

      istuff(3) = istuff(3) + 1

      if (istuff(2).ne.0.and.(nvar.lt.nstot(jds).or.
     *    sum1.ge.one.and.sum1.le.1d0+zero).and.sum.gt.zero) then
c                                 save the composition
         istuff(4) = istuff(4) + 1
c                                 increment the counter
         jphct = jphct + 1
c                                 the solution model pointer
         jkp(jphct) = jds
c                                 the refinement point pointer
         hkp(jphct) = istuff(5)
c                                 save the normalized g
         g2(jphct) = g/sum
c                                 save the normalized bulk
         cp2(1:icomp,jphct) = scp(1:icomp)/sum
c                                 sum scp(1:icp)
         if (ksmod(jds).eq.39.and.lopt(32).and..not.qwak) then 
            c2tot(jphct) = sum/smo
         else
            c2tot(jphct) = sum
         end if

         quack(jphct) = qwak
c                                 save the endmember fractions
         icoz(jphct) = zcoct

         zco(zcoct+1:zcoct+nstot(jds)) = pa(1:nstot(jds))

         zcoct = zcoct + nstot(jds)

      end if

1000  format (/,'**warning ver205** lagged speciation failed, ',
     *       'for ',a,'. The molecular',/,'speciation will be ',
     *       'output.',/)

      end

      subroutine gsol3 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
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

      double precision ppp(*), gval, psum, ggrd(*), stuff(*)

      double precision gsol1, omega0

      external gsol1, omega

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1
      save / cst59 /
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
c                                 free energy minimization
         gval = gsol1 (jds,.false.)

      else
c                                 entropy maximization
         gval = -omega0 (jds,pa)

      end if

c     write (*,1000) 0, (pa(i),i=1,nstot(jds))

c     write (*,1000) 0, (cp2(i,tphct),i=1,icp)
c     write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g12.6,1x,12(f8.5,1x))
1010  format (2(g14.7,2x))

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

      logical bad, site, comp, clos, inv

      integer liw, lw, mvar, mcon, nvar, i, j, jter, iprint, iwarn,
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
               if (pp(i).lt.-1d-2) call errdbg ('wtf, p2yx 2')
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

      if (lopt(28)) call begtim (2)

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

      if (lopt(28)) call endtim (2,.true.,'p2y inversion')

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
1040  format ('The composition will not be will not be used to compute',
     *       ' compositional ranges.',/)

      end

      subroutine minfxc (gfinal,ids,maxs)
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
      save / cst59 /

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
      nclin = nz(ids) + icomp
c                                 first the site fraction constraints
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      do i = 1, nclin
         bl(nvar+i) = zl(ids,i)
         bu(nvar+i) = zu(ids,i)
      end do
c                                 get the normalized bulk composition of the solution
      call getxcp (b,stuff(1),ids)

c     call getscp (b,sum,ids,ids,.false.)
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
         if (dabs(bl(nvar+i)).lt.zero) bl(nvar+i) = 0d0
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

         write (ctol,'(g14.7)') zero
         CALL E04UEF ('nolist')
         CALL E04UEF ('optimality tolerance = '//ctol)
         CALL E04UEF ('feasibility tolerance = '//ctol)
         CALL E04UEF ('difference interval = 1d-3')
         write (ctol,'(i4)') iprint
         CALL E04UEF ('print level = '//ctol)

      end if

      CALL E04UEF ('derivative level = 0')

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
c                                 if you touch pa, makepp
      call makepp (ids)

      if (inp) then 
         write (*,*) istuff(4),gfinal,istuff(1)
         goto 10
      end if

      end
