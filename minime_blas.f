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

      logical inp

      integer ids, kds, nvar, iter, iwork(m22),
     *        istuff(10), istate(m21), idead, nclin, ntot
c DEBUG691
     *         , i, j, iprint

      double precision ggrd(m19), lapz(m20,m19),
     *                 bl(m21), bu(m21), gfinal, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(2)
c DEBUG691                    dummies for NCNLN > 0
     *                 ,c(1),cjac(1,1),yt(m14),zp,zt(m10,m11),
     *                 ftol,fdint

      character ctol*20,cdint*20

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      double precision units, r13, r23, r43, r59, zero, one, r1
      common/ cst59 /units, r13, r23, r43, r59, zero, one, r1

      external gsol2, dummy

      data iprint,inp/0,.false./

      save iprint,inp
c-----------------------------------------------------------------------
      if (.not.mus) then 
         write (*,*) 'no mus'
         call errpau
      end if

10    nclin = nz(ids)
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 finite difference increments
c                                 will be estimated at this 
c                                 coordinate, so choose a feasible 
c                                 composition
      ppp(1:nvar) = pa(1:nvar)
c                                 initialize bounds
      bu(1:nvar) = 1d20
      bl(1:nvar) = -1d20
c                                 load the local constraints 
c                                 from the global arrays
      lapz(1:nclin,1:nvar) = apz(ids,1:nclin,1:nvar)

      bl(nvar+1:nvar+nclin) = zl(ids,1:nclin)
      bu(nvar+1:nvar+nclin) = zu(ids,1:nclin)

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

      else

         CALL E04UEF ('nolist')
         CALL E04UEF ('optimality tolerance =  1d-4')
         CALL E04UEF ('difference interval = 1d-3')

      end if

      CALL E04UEF ('derivative level = 0')

      call nlpopt (nvar,nclin,m20,m19,lapz,bl,bu,gsol2,
     *             iter,istate,clamda,gfinal,ggrd,r,ppp,iwork,
     *             m22,work,m23,istuff,stuff,idead,iprint)

c     call nlpopt (nvar,nclin,0,m20,1,m19,lapz,bl,bu,dummy,gsol2,
c    *             iter,istate,c,cjac,clamda,gfinal,ggrd,r,ppp,iwork,
c    *             m22,work,m23,istuff,stuff,istart)

c     write (*,1000) 'pas',pa(1:3)
c     write (*,1000) 'p0as',p0a(1:3)
1000  format (a10,10(g14.7,1x))

c     do i = 1, nz(ids)
c        zp = 0d0 
c        do j = 1, nvar
c           zp = zp + lapz(i,j) * ppp(j)
c        end do

c        write (*,1010) i, bl(i),zp,bu(i)

c     end do

c     call p2z (pa,zt,ids,.true.)
      if (inp) write (*,*) istuff(3),gfinal,kds
      if (inp) goto 10

      zp = 0d0
      do i = 1, nstot(ids)
         zp = zp + pa(i)
c        if (dabs(pa(i)).lt.zero) pa(i) = 0d0
      end do

      if (zp.lt.0.9999) then 
         write (*,*) 'low sum'
         call p2z (pa,zt,ids,.true.)
      end if


1010  format (i5,1x,10(g14.7,1x))
      end

      subroutine gsol2 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for minfrc. can call 
c either gsol1 that does o/d or gsol4 that does not, gsol1 seems to give
c better results presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, mode, istuff(*), istart

      double precision ppp(*), gval, gsol4, ggrd(*), stuff(*),
     *                 gsol1, g, sum, scp(k5), sum1

      external gsol4, gsol1

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
c-----------------------------------------------------------------------
      jds = istuff(1)

      sum = 0d0

      do i = 1, nvar
         sum = sum + ppp(i)
         pa(i) = ppp(i)
      end do

      pa(nstot(jds)) = 1d0 - sum

      call makepp (jds)

      g = gsol4 (jds,1d0)

c                                 get the bulk composition from pp
      call getscp (scp,sum,jds,jds,.false.)

      gval = g

      do i = 1, icp
         gval = gval - scp(i)*mu(i)
c         if (scp(i).lt.-1d2*zero) then
c            write (*,*) 'gsol2, large -composition',jds,i,scp(i)
c            call p2z (pa,zt,jds,.true.)
c            scp(i) = 0d0
c         end if
      end do
c                                  normalize for appearances
      gval = gval/sum

c     write (*,1000) gval, g

      istuff(3) = istuff(3) + 1

      if (istuff(2).ne.0) then
c                                 save the composition
         istuff(4) = istuff(4) + 1
c        write (*,1000) gval, g/sum, (pa(i),i=1,nstot(jds))
c                                 there is a small possibility 
c                                 that the composition of a phase
c                                 may move entirely into the 
c                                 constrained composition space
         if (sum.gt.zero) then
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
            c2tot(jphct) = sum
c                                 save the endmember fractions
            icoz(jphct) = zcoct

            zco(zcoct+1:zcoct+nstot(jds)) = pa(1:nstot(jds))

            zcoct = zcoct + nstot(jds)

            sum = 0d0
            do i = 1, nstot(jds)
               sum = sum + pa(i)
            end do 
            if (sum.lt.0.9999) then
               write (*,*) 'wug'
            end if

         end if

      end if

1000  format (2(g12.6,1x),12(f8.5,1x))
1010  format (2(g14.7,2x))

      end

      subroutine gsol3 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, istart, mode, istuff(*)

      double precision ppp(*), gval, psum, ggrd(*), stuff(*)

      double precision gsol4, omega

      external gsol4, omega

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      double precision apc, endt, endc
      common/ cstp2c /apc(h9,k5,m14), endt(h9,m14), endc(h9,m14,k5)

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
c                                 for non-equimolar ordering
c                                 get the total moles
      stuff(2) = 0d0

      do i = 1, nstot(jds)
         stuff(2) = stuff(2) + pa(i) * endt(jds,i)
      end do

      if (istuff(3).eq.0d0) then
c                                 free energy minimization
c                                 gsol4 returns the g for stuff(1) moles
         gval = gsol4 (jds,stuff(2)/stuff(1))
c        gval = gsol4 (jds,1d0)
      else
c                                 entropy maximization
         gval = -stuff(2)/stuff(1) * omega(jds,pa)
      end if

c     write (*,1000) 0, (pa(i),i=1,nstot(jds))

c     write (*,1000) 0, (cp2(i,tphct),i=1,icp)
c     write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g12.6,1x,12(f8.5,1x))
1010  format (2(g14.7,2x))

      end

      double precision function gsol4 (id,norm)
c-----------------------------------------------------------------------
c gsol4 computes the total (excess+ideal) free energy of solution
c for a solution identified by index ids for the speciation input
c via cxt7, i.e., in contrast to gsol1 it does not compute the 
c speciation of o/d models. the energy is scaled to stuff1 moles
c by norm.

c ingsol MUST be called prior to gsol4 to initialize solution
c specific parameters! 

c gsol4 assumes the endmember g's have been calculated by gall.

c gsol4 is called only for implicit order-disorder models by minfxc.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer k, id

      double precision g1, g2, norm

      double precision omega, gex, gfesi

      external omega, gex

      integer jend
      common/ cxt23 /jend(h9,m4)

      double precision g
      common/ cst2 /g(k1)

      double precision enth
      common/ cxt35 /enth(j3)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                 for non-equimolar ordering properties
c                                 computed with pa are for stuff2 moles
c                                 and those computed with pp are for 
c                                 stuff1 moles
      g1 = 0d0
c                                 dqf contribution
      call gdqf (id,g1,pp)
c                                 mechanical contribution
      do k = 1, lstot(id)
         g1 = g1 + g(jend(id,2+k)) * pp(k)
      end do
c                                 entropic + excess o/d effect
      g2 = - t * omega(id,pa) + gex(id,pa)
c                                 enthalpic effect o/d effct
      do k = 1, nord(id)
         g2 = g2 + pa(lstot(id)+k)*enth(k)
      end do
c                                 gsol4 is the energy per stuff1 moles. 
      gsol4 = g1 + g2 / norm

      end

      subroutine p2yx (id)
c-----------------------------------------------------------------------
c converts the independent endmember fractions to 0-1 bounded barycentric 
c coordinates:

c     number of bounded coordinates -> mstot (<m4)
c     number of independent fractions -> nstot (<m14)
c     number of linear constraints -> the number of independent
c        site fractions + closure (<m20)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw,lw,mvar,mcon

      parameter (mvar=m4, mcon=m20, liw=2*mvar+3, 
     *           lw=2*(mcon+1)**2+7*mvar+5*mcon)

      integer ncon, id, is(mvar+mcon), iw(liw), idead, istart

      double precision ax(mcon), clamda(mvar+mcon), wrk(lw), c(mvar),
     *                 a(mcon,mvar), b(mcon), gopt

      double precision ayz
      common/ csty2z /ayz(h9,m20,m4)

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer ldq,ldt,ncolt
      common/ be04nb /ldt,ncolt,ldq
c-----------------------------------------------------------------------
c                                 get the disordered p's
      call minfxc (gopt,id,.true.)
c                                 get the site fraction constraints
      call p2zind (pa,b,id)
c                                 dummy objective function coefficients
c                                 (only 1 feasible point?)
      c(1:mstot(id)) = 1d0
c                                 load the ayz constraint matrix
      a(1:nz(id),1:mstot(id)) = ayz(id,1:nz(id),1:mstot(id))
c                                 add the closure constraint
      ncon = nz(id) + 1
      a(ncon,1:mstot(id)) = 1d0
      b(ncon) = 1d0
c                                 cold start
      istart = 0
      ldt = ncon
      ldq = ldt

      if (lopt(28)) call begtim (2)
c                                 optimize by nag
      call lpnag (mstot(id),ncon,a,mcon,b,c,is,y,ax,
     *            clamda,iw,liw,wrk,lw,idead,l6,istart)

      if (lopt(28)) call endtim (2,.true.,'p2y inversion')

      if (idead.gt.0) then
c                                 look for severe errors
         call lpwarn (idead,'LPOPT ')
         call errpau

      end if
c                                 reset ldt, ldq, istart for phase eq
      ldq = icp + 1
      ldt = ldq
      istart = 0
c                                 convert the y's to x's
      call sety2x (id,y)

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
c DEBUG691
     *                 ,xp(m14), ftol, fdint
      character*14 cdint, ctol

      integer nz
      double precision apz, zl, zu
      common/ cstp2z /apz(h9,m20,m19), zl(h9,m20), zu(h9,m20), nz(h9)

      double precision apc, endt, endc
      common/ cstp2c /apc(h9,k5,m14), endt(h9,m14), endc(h9,m14,k5)

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

      external gsol3
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

      ppp(1:nvar) = pa(1:nvar)
      xp(1:ntot) = pa(1:ntot)
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

      else

         iprint = 0 

         CALL E04UEF ('nolist')
         CALL E04UEF ('optimality tolerance =  1d-8')
         CALL E04UEF ('difference interval = 1d-3')

      end if

      CALL E04UEF ('derivative level = 0')

      call nlpopt (nvar,nclin,m20,m19,lapz,bl,bu,gsol3,
     *             iter,istate,clamda,gfinal,ggrd,r,ppp,iwork,
     *             m22,work,m23,istuff,stuff,idead,iprint)

      if (idead.eq.2) then 
         write (*,*) 'minfxc infeasible initial conditions'
         call errpau
      end if

      sum = 0d0

      do i = 1, nvar
         pa(i) = ppp(i)
         sum = sum + ppp(i)
      end do

      pa(i) = 1d0 - sum

      if (inp) then 
         write (*,*) istuff(4),gfinal,istuff(1)
         goto 10
      end if

c     do i = 1, nstot(ids)
c        if (dabs(pa(i)).lt.zero) pa(i) = 0d0
c     end do 

      end
