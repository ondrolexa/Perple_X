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

c DEBUG691 gall
      double precision deph,dydy,dnu
      common/ cxt3r /deph(3,j3,h9),dydy(m4,j3,h9),dnu(h9)

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

      if (dnu(ids).lt.zero) then
         nvar = ntot - 1
      else 
         nvar = ntot
      end if 
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

      if (nvar.eq.ntot) then
         nclin = nclin + 1
         bl(nvar+nclin) = 1d0
         bu(nvar+nclin) = 1d0
         lapz(nclin,1:nvar) = 1d0
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

      if (inp) then 

         zp = 0d0
         do i = 1, nstot(ids)
            zp = zp + ppp(i)
         end do

         write (*,1000) 'pas',pa(1:3)
         write (*,1000) 'p0as',p0a(1:3)

         do i = 1, nclin
            zp = 0d0 
            do j = 1, nvar
               zp = zp + lapz(i,j) * ppp(j)
            end do

            write (*,1010) i, bl(i),zp,bu(i)

         end do

         call p2z (pa,zt,ids,.true.)

         write (*,*) istuff(3),gfinal,kds

         goto 10

      end if

c      zp = 0d0
c      do i = 1, nstot(ids)
c         zp = zp + ppp(i)
c        if (dabs(pa(i)).lt.zero) pa(i) = 0d0
c      end do

c      if (zp.lt.one.or.zp.gt.1d0+zero) then 
c         write (*,*) 'low sum'
c         call p2z (ppp,zt,ids,.true.)
c      end if

1000  format (a10,10(g14.7,1x))
1010  format (i5,1x,10(g14.7,1x))
      end

      subroutine gsol2 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
c function to evaluate gibbs energy of a solution for minfrc. can call 
c either gsol1 with order true or false, true seems to give better results
c presumably because it's using analytical gradients.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, mode, istuff(*), istart

      double precision ppp(*), gval, ggrd(*), stuff(*),
     *                 gsol1, g, sum, scp(k5), sum1

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
c DEBUG691 gall
      double precision deph,dydy,dnu
      common/ cxt3r /deph(3,j3,h9),dydy(m4,j3,h9),dnu(h9)

c-----------------------------------------------------------------------
      jds = istuff(1)

      sum1 = 0d0

      do i = 1, nvar
         sum1 = sum1 + ppp(i)
         pa(i) = ppp(i)
      end do

      if (nvar.lt.nstot(jds)) pa(nstot(jds)) = 1d0 - sum1

      call makepp (jds)
c                                 T use explicit ordering
      g = gsol1 (jds,.false.)
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

c      if (sum1.lt.one.or.sum1.ge.1d0+zero) then 
c         write (*,*) 'whoa nelly, ',sum1, istuff(2)
c      else 
c         write (*,*) 'heigh ho silver',istuff(4)
c      end if 

      if (istuff(2).ne.0.and.sum1.ge.one.and.sum1.le.1d0+zero) then
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

            if (sum.lt.one.or.sum.gt.1d0+zero) then
               write (*,*) 'wug'
            end if

         end if

      end if

1000  format (2(g12.6,1x),12(f8.5,1x))
1010  format (2(g14.7,2x))

      end

      subroutine gsol3 (mode,nvar,ppp,gval,ggrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
c gsol3 - a shell to call gsol1 from minfxc, ingsol must be called
c         prior to minfxc to initialize solution specific paramters. only
c         called for equimolar explicit o/d models. non-equimolar o/d models
c         (currently melt(G,HGP)) involve non-linear constraints that are not
c         currently implemented in minfxc/nlpopt.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jds, nvar, istart, mode, istuff(*)

      double precision ppp(*), gval, psum, ggrd(*), stuff(*)

      double precision gsol1, omega

      external gsol1, omega

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

      if (istuff(3).eq.0d0) then
c                                 free energy minimization
         gval = gsol1 (jds,.false.)

      else
c                                 entropy maximization
         gval = -omega (jds,pa)

      end if

c     write (*,1000) 0, (pa(i),i=1,nstot(jds))

c     write (*,1000) 0, (cp2(i,tphct),i=1,icp)
c     write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g12.6,1x,12(f8.5,1x))
1010  format (2(g14.7,2x))

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
