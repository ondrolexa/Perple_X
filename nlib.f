
c these are nag routines modified to solve the phase eq problem, whereby
c the basic constraints are 0<=x<=1 and the general constraints are ax=b

c many of the routines no longer perform the functions indicated in the 
c nag documentation or source comments. 

      subroutine lpnag (n,nclin,a,lda,bl,cvec,istate,x,ax,
     *                  clamda,iw,leniw,w,lenw,idead,itmax,istart)
c----------------------------------------------------------------------
c     a  is a constant  nclin by n  matrix.
c     n  is the number of variables (dimension of x).
c     nclin  is the number of general linear constraints (rows of  a).

      implicit none

      double precision  obj
      integer           idead, iter, lda, leniw, lenw, n, nclin
      double precision  a(lda,*), ax(*), bl(nclin), 
     *                  clamda(n+nclin), cvec(*), w(lenw), x(n)
      integer           istate(n+nclin), iw(leniw)

      double precision  amin,errmax,feamax,feamin,xnorm,dnrm2 

      integer           ianrmj, it, j, jinf, jmax, istart,
     *                  nact1, nactiv, nartif,
     *                  nfree, ngq, nmoved, nrejtd, nrz, numinf,
     *                  nviol, nz, itmax

      logical           done, found, halted, rowerr, unitq

      character*6       msg

      integer loclc
      common/ ae04mf /loclc(15)

      double precision  wmach(9)
      common/ ax02za /wmach

      integer ldq,ldt
      common/ be04nb /ldt,ldq

      double precision  asize,dtmax, dtmin
      common/ de04nb /asize, dtmax, dtmin

      integer nctotl,lkx,lanorm,ld,lgq,lcq,lrlam,lt,lq,lwtinf,lwrk

      save nctotl,lkx,lanorm,ld,lgq,lcq,lrlam,lt,lq,lwtinf,lwrk
c----------------------------------------------------------------------   
      iter = 0

      if (istart.eq.0.or.idead.eq.-1) then 
c                                 cold start (istart = 0) or reusing 
c                                 matrices from a prior optimization 
c                                 (idead = -1)
         nctotl = n + nclin
         lkx = 4 + n
         lanorm = 1 + nclin + n
         ld = lanorm + nclin + nclin
         lgq = ld + n
         lcq = lgq + n
         lrlam = lcq + n
         lt = lrlam + n
         lq = lt + ldt*ldt
         lwtinf = lq + ldq*ldq
         lwrk = lwtinf + n + nclin

         loclc(4) = lanorm
         loclc(5) = lanorm + nclin
         loclc(7) = lanorm + nclin + nclin
         loclc(8) = lgq
         loclc(9) = lcq
         loclc(10) = lrlam
         loclc(12) = lt
         loclc(13) = lq
         loclc(14) = lwtinf
         loclc(15) = lwrk

      end if 

      idead = 0 
c                                     define the initial feasibility tolerances in clamda.
      call f06fbf(n+nclin,wmach(4),w(1),1)
      call e04mfr('i',n,nclin,nmoved,iter,numinf,istate,clamda,w(1),x)

      if (istart.eq.0) then
c                                     true cold start initialization
         ianrmj = lanorm
         do j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
         end do 

         call f06flf(nclin,w(lanorm),1,asize,amin)
         call f06flf(nctotl,w(1),1,feamax,feamin)
         call dcopy(nctotl,w(1),1,w(lwtinf),1)
         call dscal(nctotl,1d0/feamin,w(lwtinf),1)

c        define the initial working set.

         call lpinit (nctotl,nactiv,nartif,nfree,n,
     *                istate,iw(4),iw(lkx),x,w(lgq))

c        compute the tq factorization of the working set matrix.

         unitq = .true.
         nz = nfree

         if (nactiv.gt.0) then
            it = nactiv + 1
            nact1 = nactiv
            nactiv = 0
            ngq = 0

            call e04nfq(unitq,1,nact1,it,nactiv,nartif,nz,nfree,
     *                  nrejtd,ngq,n,ldq,lda,ldt,istate,iw(4),
     *                  iw(lkx),wmach(9),a,w(lt),w(lgq),w(lq),w(lwrk),
     *                  w(ld),w(lrlam))

            if (nactiv.eq.0) then 

               write (*,*) 'i am gonna crash, set msg and die, but why?'

               idead = 7 

               return

            end if   
          
         end if

      else 
c                                 warm start
         unitq = iw(1) .eq. 1
         nfree = iw(2)
         nactiv = iw(3)
         nz = nfree - nactiv

      end if

c        install the transformed linear term in cq.

      call dcopy(n,cvec,1,w(lcq),1)
      call e04nbw(6,n,nz,nfree,ldq,unitq,iw(lkx),w(lcq),w(lq),w(lwrk))

      jinf = 0
      nrz = 0

c     repeat               (until working set residuals are acceptable)

c     move x onto the constraints in the working set.

   40 call e04mfj(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,ldt,
     *            istate,iw(4),iw(lkx),jmax,errmax,xnorm,a,ax,bl,
     *            w(1),w(lt),x,w(lq),w(ld),w(lwrk))

      if (rowerr) then
         msg = 'infeas'
         obj = errmax
         go to 60
      end if

      call e04mfz(msg,unitq,iter,itmax,
     *            jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,nz,istate,
     *            iw(4),iw(lkx),obj,numinf,xnorm,a,ax,bl,
     *            cvec,clamda,w(1),x,w)

      found = msg .eq. 'feasbl' .or. msg .eq. 'optiml' .or. msg .eq.
     *        'weak  ' .or. msg .eq. 'unbndd' .or. msg .eq. 'infeas'
      halted = msg .eq. 'itnlim'

      if (found) call e04mfr('o',n,nclin,nmoved,iter,numinf,istate,
     *               clamda,w(1),x)

      done = found.and.nviol .eq. 0.and.nmoved .eq. 0

      if (.not. (done .or. halted)) go to 40

      call e04mfk(nfree,n,nctotl,nactiv,iw(4),iw(lkx),
     *            clamda,w(lrlam))

      iw(1) = 0
      if (unitq) iw(1) = 1
      iw(2) = nfree
      iw(3) = nactiv

   60 if (msg.eq.'optiml') then
         idead = 0
      else if (msg.eq.'feasbl') then
         idead = 0
      else if (msg.eq.'weak  ') then
         idead = 0
      else if (msg.eq.'unbndd') then
         idead = 2
      else if (msg.eq.'infeas') then
         idead = 3
      else if (msg.eq.'itnlim') then
         idead = 4
      else if (msg.eq.'errors') then
         idead = 6
      end if

      istart = 1 

      end

      subroutine e04mfk(nfree,n,nctotl,nactiv,kactiv,kx,
     *           clamda,rlamda)
c----------------------------------------------------------------------
c     mark 16 release. nag copyright 1992.
c     e04mfk   creates the expanded lagrange multiplier vector clamda.

      integer           n, nactiv, nctotl, nfree
      double precision  clamda(nctotl), rlamda(n)
      integer           kactiv(n), kx(n)
      double precision  rlam
      integer           j, k, nfixed, nz
c----------------------------------------------------------------------


      nz = nfree - nactiv

c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      call f06fbf (nctotl,0d0,clamda,1)

      nfixed = n - nfree

      do k = 1, nactiv + nfixed
         if (k.le.nactiv) then
            j = kactiv(k) + n
            rlam = rlamda(nactiv-k+1)
         else
            j = kx(nz+k)
            rlam = rlamda(k)
         end if
         clamda(j) = rlam
      end do

      end

      subroutine e04mfz(msg,unitq,iter,itmax,
     *                  jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,
     *                  nz,istate,kactiv,kx,obj,numinf,xnorm,a,
     *                  ax,bl,cvec,featol,featlu,x,w)
c----------------------------------------------------------------------
c     e04mfz  is a subroutine for linear programming.
c     on entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.  the arrays  istate,
c     kactiv  and  kx  will have been set accordingly
c     and the arrays  t  and  q  will contain the tq factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  the tq factorization of the resulting (nactiv by nfree)
c     matrix is  a(free)*q = (0 t),  where q is (nfree by nfree) and t
c     is upper-triangular.

c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.

c     after kdegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.


c     values of istate(j) for the linear constraints.......

c     istate(j)
c     ---------
c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.

c     constraint j may be violated by as much as featol(j).
c----------------------------------------------------------------------
      implicit none 

      character*6       empty
      parameter         (empty='      ')
      double precision  obj, xnorm
      integer           iter, jinf, lda, n, nactiv, nclin, nfree,
     *                  nrz, numinf, nviol, nz
      logical           unitq

      character*6       msg

      double precision  a(lda,*), ax(*), bl(nclin), 
     *                  cvec(*), featlu(n+nclin), featol(n+nclin), w(*),
     *                  x(n)
      integer           istate(n+nclin), kactiv(n), kx(n)

      double precision  alfa, asize,dtmax, dtmin
      integer           jadd, jdel,ldq, ldt, itmax

      double precision  wmach(9)
      integer           loclc
      double precision  alfap, alfhit, bigalf, biggst, condmx, 
     *                  dinky, dnorm, errmax, flmax,
     *                  gfnorm, grznrm, gznorm, objsiz, smllst,
     *                  suminf, tinyst, trubig, trusml, wssize, zerolm
      integer           iadd, ifix, inform, is, it, j, jbigst,
     *                  jmax, jsmlst, jtiny, kbigst, kdel, ksmlst, lad,
     *                  lanorm, lcq, ld, lgq, lq, lrlam, lt,
     *                  lwrk, lwtinf,
     *                  nctotl, ngq, nmoved, notopt, ntfixd
      logical           firstv, hitlow, move, onbnd, overfl,
     *                  unbndd, fail

      double precision  ddot, dnrm2, adivb

      common            /ae04mf/loclc(15)
      common            /ax02za/wmach
      common            /be04nb/ldt,ldq
      common            /de04nb/asize, dtmax, dtmin

      save              firstv, alfa, jdel, jadd
c----------------------------------------------------------------------

c     specify the machine-dependent parameters.

      flmax = wmach(7)

      ngq = 2
      it = 1
      lanorm = loclc(4)
      lad = loclc(5)
      ld = loclc(7)
      lgq = loclc(8)
      lcq = loclc(9)
      lrlam = loclc(10)
      lt = loclc(12)
      lq = loclc(13)
      lwtinf = loclc(14)
      lwrk = loclc(15)

c     we need a temporary array when changing the active set.
c     use the multiplier array.

      if (iter.eq.0) then
         jadd = 0
         jdel = 0
         firstv = .false.
         alfa = 0d0
      end if

      nctotl = n + nclin
      nviol = 0
      condmx = flmax

      call e04mfh(n,nclin,lda,istate,numinf,suminf,bl,a,
     *            featol,w(lgq),x,w(lwtinf))

      if (numinf.gt.0) then
         call e04nbw(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),w(lwrk))
      else 
         call dcopy(n,w(lcq),1,w(lgq),1)
      end if

      if (numinf.eq.0) then
         obj = ddot(n,cvec,1,x)
      else
         obj = suminf
      end if

      msg = empty

c*    ======================start of main loop==========================
c     +    do while (msg .eq. empty)
   20 if (msg.eq.empty) then

         gznorm = 0d0
         if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)

         if (nrz.eq.nz) then
            grznrm = gznorm
         else
            grznrm = 0d0
            if (nrz.gt.0) grznrm = dnrm2(nrz,w(lgq),1)
         end if

         gfnorm = gznorm
         if (nfree.gt.0.and.nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),
     *       1)

         if (numinf.gt.0) then
            dinky = wmach(2)*dabs(suminf)
         else
            objsiz = 1d0 + dabs(obj)
            wssize = 0d0
            if (nactiv.gt.0) wssize = dtmax
            dinky = wmach(2)*max(wssize,objsiz,gfnorm)
         end if

c        if the reduced gradient z'g is small enough,
c        lagrange multipliers will be computed.
         if (grznrm.le.dinky) then

c           the point  x  is a constrained stationary point.
c           compute lagrange multipliers.
c           define what we mean by 'tiny' and non-optimal multipliers.

            notopt = 0
            jdel = 0
            zerolm = -dinky
            smllst = -dinky
            biggst = dinky + 1d0
            tinyst = dinky

            call e04mfm(n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,w(lanorm),w(lgq),w(lrlam),w(lt),
     *                  w(lwtinf))

            if (nrz.lt.nz) call e04mfl(n,nrz,nz,zerolm,notopt,
     *                                 numinf,trusml,smllst,jsmlst,
     *                                 tinyst,jtiny,w(lgq))

            if (abs(jsmlst).gt.0) then

c              delete a constraint.

c              e04mfm  or  e04mfl  found a non-optimal multiplier.

               jdel = jsmlst

               if (jsmlst.gt.0) then

c                 regular constraint.

                  kdel = ksmlst
                  istate(jdel) = 0
               end if
            else
               if (numinf.gt.0.and.jbigst.gt.0) then

c                 no feasible point exists for the constraints but
c                 the sum of the constraint violations can be reduced
c                 by moving off constraints with multipliers greater
c                 than 1.

                  jdel = jbigst
                  kdel = kbigst
                  if (trubig.le.0d0) is = -1
                  if (trubig.gt.0d0) is = -2
                  istate(jdel) = is
                  firstv = .true.
                  numinf = numinf + 1
               end if
            end if

            if (jdel.eq.0) then
               if (numinf.gt.0) then
                  msg = 'infeas'
               else
                  msg = 'optiml'
               end if
               go to 20
            end if

c           constraint  jdel  has been deleted.
c           update the  tq  factorization.

            call e04nfp(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,w(lt),w(lgq),w(lq)
     *                  ,w(ld),w(lrlam),fail)

            if (fail) then
c                                 this error appeared dec 07,
c                                 possibly due to parameter edits?
c                                 JADC 12/07 
               msg = 'infeas'
               goto 20 
            end if 

         else

c           compute a search direction.

            if (iter.ge.itmax) then
               msg = 'itnlim'
               go to 20
            end if

            iter = iter + 1

            call dcopy(nrz,w(lgq),1,w(ld),1)
            call dscal(nrz,-1d0,w(ld),1)

            dnorm = dnrm2(nrz,w(ld),1)

            call e04nbw(1,n,nrz,nfree,ldq,unitq,kx,w(ld),w(lq),w(lwrk))
            call dgemv('n',nclin,n,a,lda,w(ld),0d0,w(lad))

c           find the constraint we bump into along d.
c           update  x  and  ax  if the step alfa is nonzero.

c           alfhit is initialized to bigalf. if it remains that value
c           after the call to  e04mfs, it is regarded as infinite.

            bigalf = adivb(1d20,dnorm,overfl)

            call e04mfs(firstv,n,nclin,istate,bigalf,dnorm,
     *                  hitlow,move,onbnd,unbndd,alfhit,alfap,jadd,
     *                  w(lanorm),w(lad),ax,bl,featol,w(ld),x)

            if (unbndd) then
               msg = 'unbndd'
               go to 20
            end if

            alfa = alfhit
            call daxpy(n,alfa,w(ld),1,x)

            call daxpy(nclin,alfa,w(lad),1,ax)
            xnorm = dnrm2(n,x,1)

c           add a constraint to the working set.
c           update the  tq  factors of the working set.
c           use  d  as temporary work space.

            if (jadd.gt.n) then
               iadd = jadd - n
               istate(jadd) = 3              
            else
               if (hitlow) then
                  istate(jadd) = 1
               else
                  istate(jadd) = 2
               end if

               if (alfa.ge.0d0) then
                  if (hitlow) then
                     x(jadd) = 0d0
                  else
                     x(jadd) = 1d0
                  end if
               end if
               do ifix = 1, nfree
                  if (kx(ifix).eq.jadd) exit
               end do 
            end if

            call e04nfr(unitq,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,ngq,n,lda,ldq,ldt,kx,condmx,a,
     *                  w(lt),w(lgq),w(lq),w(lwrk),w(lrlam),w(ld))

            nz = nz - 1
            nrz = nrz - 1

            if (jadd.le.n) then

c              a simple bound has been added.

               nfree = nfree - 1
            else

c              a general constraint has been added.

               nactiv = nactiv + 1
               kactiv(nactiv) = iadd
            end if

c           increment featol.

            call daxpy(nctotl,0d0,featlu,1,featol)

            if (mod(iter,50).eq.0) then

c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  set inform to force iterative
c              refinement and a switch to phase 1.

               call e04mfq(n,nclin,istate,nviol,jmax,errmax,ax,
     *                     bl,featol,x)

            end if

            if (mod(iter,999999).eq.0) then

c              every  kdegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.

               call e04mfr('e',n,nclin,nmoved,iter,
     *                     numinf,istate,featol,featlu,x)

               nviol = nviol + nmoved
            end if

            if (nviol.gt.0) then
               msg = 'resetx'
               go to 20
            end if

            if (numinf.ne.0) then
               call e04mfh(n,nclin,lda,istate,numinf,suminf,bl,
     *                     a,featol,w(lgq),x,w(lwtinf))

               if (numinf.gt.0) then
                  call e04nbw(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),
     *                        w(lwrk))
               else 
                  call dcopy(n,w(lcq),1,w(lgq),1)
               end if
            end if

            if (numinf.eq.0) then
               obj = ddot(n,cvec,1,x)
            else
               obj = suminf
            end if
         end if
         go to 20
c        +    end while
      end if
c     ======================end of main loop============================

      if (msg.eq.'optiml') then

            if (nrz.lt.nz) then
               msg = 'weak  '
            else
               ntfixd = 0
               do 80 j = 1, n
                  if (istate(j).eq.4) ntfixd = ntfixd + 1
   80          continue
               if (ntfixd.gt.0) msg = 'weak  '
            end if
            if (abs(jtiny).gt.0) msg = 'weak  '

      else if (msg.eq.'unbndd'.and.numinf.gt.0) then
         msg = 'infeas'
      end if

      end

      subroutine e04mfj(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,
     *                  ldt,istate,kactiv,kx,jmax,errmax,xnorm,a,ax,bl,
     *                  featol,t,x,q,p,work)
c     mark 16 release. nag copyright 1992.

c     ******************************************************************
c     e04mfj  computes the point on a working set that is closest in the
c     least-squares sense to the input vector x.

c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable rowerr
c     is set.

c     original version derived from lssetx january-1987.
c     this version of  e04mfj  dated   5-jul-1989.
c     ******************************************************************

      integer           ntry
      parameter         (ntry=5)

      double precision  errmax, xnorm
      integer           jmax, lda, ldq, ldt, n, nactiv, nclin, nfree, nz
      logical           rowerr, unitq

      double precision  a(lda,*), ax(*), bl(nclin), 
     *                  featol(n+nclin), p(n), q(ldq,*), t(ldt,*),
     *                  work(n), x(n)
      integer           istate(n+nclin), kactiv(n), kx(n)
c     .. scalars in common ..
      

c     .. arrays in common ..

c     .. local scalars ..
      double precision  bnd
      integer           i, is, j, k, ktry
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           idamax

c     ------------------------------------------------------------------
c     move  x  onto the simple bounds in the working set.
c     ------------------------------------------------------------------
      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = 0d0
         if (is.ge.2) bnd = 1d0
         if (is.ne.4) x(j) = bnd
   20 continue
c     ------------------------------------------------------------------
c     move  x  onto the general constraints in the working set.
c     ntry  attempts are made to get acceptable row errors.
c     ------------------------------------------------------------------
      ktry = 1
      jmax = 1
      errmax = 0d0

c     repeat

c        set work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.

40       do 60 i = 1, nactiv
            k = kactiv(i)
            work(nactiv-i+1) = bl(k) - ddot(n,a(k,1),lda,x)
   60    continue

         call dtrsv('n',nactiv,t(1,nz+1),ldt,work)
         call f06fbf(n,0d0,p,1)
         call dcopy(nactiv,work,1,p(nz+1),1)
         call e04nbw(2,n,nz,nfree,ldq,unitq,kx,p,q,work)
         call daxpy(n,1d0,p,1,x)


c     ---------------------------------------------------------------
c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.
c     ---------------------------------------------------------------
      xnorm = dnrm2(n,x,1)
      call dgemv('n',nclin,n,a,lda,x,0d0,ax)

c     ---------------------------------------------------------------
c     check the row residuals.
c     ---------------------------------------------------------------

         do k = 1, nactiv
            i = kactiv(k)
            is = istate(n + i)
            if (is.ge.1) work(k) = bl(i) - ax(i)
         end do

         jmax = idamax(nactiv,work)
         errmax = dabs(work(jmax))

      ktry = ktry + 1

      if (.not. (errmax.le.featol(jmax) .or. ktry.gt.ntry)) go to 40

      rowerr = errmax .gt. featol(jmax)

      end

      subroutine e04nbw(mode,n,nz,nfree,nq,unitq,kx,v,zy,wrk)

      integer           mode, n, nfree, nq, nz
      logical           unitq
      double precision  v(n), wrk(n), zy(nq,*)
      integer           kx(n)
      integer           j, j1, j2, k, l, lenv, nfixed


      nfixed = n - nfree
      j1 = 1
      j2 = nfree
      if (mode.eq.1 .or. mode.eq.4) j2 = nz
      if (mode.eq.2 .or. mode.eq.5 .or. mode.eq.7) j1 = nz + 1
      lenv = j2 - j1 + 1

      if (mode.le.3) then

         if (nfree.gt.0) call f06fbf(nfree,0d0,wrk,1)

c        copy  v(fixed)  into the end of  wrk.

         if (mode.ge.2.and.nfixed.gt.0) call dcopy(nfixed,v(nfree+1),
     *       1,wrk(nfree+1),1)

c        set  wrk  =  relevant part of  zy * v.

         if (lenv.gt.0) then
            if (unitq) then
               call dcopy(lenv,v(j1),1,wrk(j1),1)
            else
               call dgemv('n',nfree,j2-j1+1,zy(1,j1),nq,v(j1),1d0,wrk)
            end if
         end if

c        expand  wrk  into  v  as a full n-vector.

         call f06fbf(n,0d0,v,1)
         do 20 k = 1, nfree
            j = kx(k)
            v(j) = wrk(k)
   20    continue

c        copy  wrk(fixed)  into the appropriate parts of  v.

         if (mode.gt.1) then
            do 40 l = 1, nfixed
               j = kx(nfree+l)
               v(j) = wrk(nfree+l)
   40       continue
         end if

      else
c        ===============================================================
c        mode = 4, 5, 6, 7  or  8.
c        ===============================================================
c        put the fixed components of  v  into the end of  wrk.

         if (mode.eq.5 .or. mode.eq.6) then
            do 60 l = 1, nfixed
               j = kx(nfree+l)
               wrk(nfree+l) = v(j)
   60       continue
         end if

c        put the free  components of  v  into the beginning of  wrk.

         if (nfree.gt.0) then
            do 80 k = 1, nfree
               j = kx(k)
               wrk(k) = v(j)
   80       continue

c           set  v  =  relevant part of  zy' * wrk.

            if (lenv.gt.0) then
               if (unitq) then
                  call dcopy(lenv,wrk(j1),1,v(j1),1)
               else
                  call dgemv('t',nfree,j2-j1+1,zy(1,j1),nq,wrk,
     *                       0d0,v(j1))
               end if
            end if
         end if

c        copy the fixed components of  wrk  into the end of  v.

         if (nfixed.gt.0.and.(mode.eq.5 .or. mode.eq.6))
     *       call dcopy(nfixed,wrk(nfree+1),1,v(nfree+1),1)
      end if

      end


      subroutine e04nfq(unitq,k1,k2,it,nactiv,nartif,nz,nfree,
     *                  nrejtd,ngq,n,ldq,lda,ldt,istate,kactiv,kx,
     *                  condmx,a,t,gqm,q,w,c,s)
c     mark 16 release. nag copyright 1992.

c     ******************************************************************
c     e04nfq  includes general constraints  k1  thru  k2  as new rows of
c     the  tq  factorization:
c              a(free) * q(free)  = ( 0 t)
c                        q(free)  = ( z y)

c     a) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt)  of the array  t.

c     original version written by peg,  october-31-1984.
c     this version of  e04nfq  dated  7-jul-1989.
c     ******************************************************************

      implicit none

      double precision  condmx
      integer           it, k1, k2, lda, ldq, ldt, n, nactiv,
     *                  nartif, nfree, ngq, nrejtd, nz
      logical           unitq
      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n),
     *                  t(ldt,*), w(n)
      integer           istate(*), kactiv(n), kx(n)

      double precision  asize, dtmax, dtmin
      double precision  wmach(9)
      double precision  cndmax, cond, delta, dtnew, rnorm, rowmax,
     *                  rtmax, tdtmax, tdtmin
      integer           i, iadd, iartif, ifix, inform, iswap, j, jadd,
     *                  jt, k, l, nzadd
      logical           overfl
      double precision  dnrm2, adivb
      common            /ax02za/wmach
      common            /de04nb/asize, dtmax, dtmin


      rtmax = wmach(8)

      jt = nz + 1

c     estimate the condition number of the constraints already
c     factorized.

      if (nactiv.eq.0) then
         dtmax = 0d0
         dtmin = 1d0

         if (unitq) then

c           first general constraint added.  set  q = i.

            call f06qhf(nfree,nfree,q,ldq)
            unitq = .false.
         end if
      else
         call f06flf(nactiv,t(it,jt),ldt+1,dtmax,dtmin)
      end if

      do k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv.lt.nfree) then

c           transform the incoming row of  a  by  q'.

            call dcopy(n,a(iadd,1),lda,w,1)
            call e04nbw(8,n,nz,nfree,ldq,unitq,kx,w,q,s)

c           check that the incoming row is not dependent upon those
c           already in the working set.

            dtnew = dnrm2(nz,w,1)
            if (nactiv.eq.0) then

c              this is the first general constraint in the working set.

               cond = adivb(asize,dtnew,overfl)
               tdtmax = dtnew
               tdtmin = dtnew
            else

c              there are already some general constraints in the working
c              set. update the estimate of the condition number.

               tdtmax = max(dtnew,dtmax)
               tdtmin = min(dtnew,dtmin)
               cond = adivb(tdtmax,tdtmin,overfl)
            end if

            if (cond.ge.condmx .or. overfl) then

c              this constraint appears to be dependent on those already
c              in the working set.  skip it.

               istate(jadd) = 0
               kactiv(k) = -kactiv(k)
            else
               if (nz.gt.1) then

c                 use a single column transformation to reduce the first
c                 nz-1  elements of  w  to zero.

c                 apply the householder reflection  i  -  w w'.
c                 the reflection is applied to  z  and gqm so that
c                    y  =    z  * w,   z    =  z    -  y w'  and
c                    y  =  gqm' * w,   gqm  =  gqm  -  w y',
c                 where  w = wrk1 (from householder),
c                 and    y = wrk2 (workspace).

c                 note that delta  has to be stored after the reflection
c                 is used.

                  delta = w(nz)

                  call f06frf(nz-1,delta,w,w(nz))

                  if (w(nz).gt.0d0) then

                     call dgemv ('n',nfree,nz,q,ldq,w,0d0,s)
                     call dger (nfree,nz,s,w,q,ldq)

                     if (ngq.gt.0) then
                        call dgemv('t',nz,ngq,gqm,n,w,0d0,s)
                        call dger (nz,ngq,w,s,gqm,n)
                     end if

                  end if

                  w(nz) = delta

               end if

               it = it - 1
               jt = jt - 1
               nactiv = nactiv + 1
               nz = nz - 1
               call dcopy(nactiv,w(jt),1,t(it,jt),ldt)
               dtmax = tdtmax
               dtmin = tdtmin

            end if

         end if

      end do 

      if (nactiv.lt.k2) then

c        some of the constraints were classed as dependent and not
c        included in the factorization.  re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.

         l = k1 - 1
         do k = k1, k2
            i = kactiv(k)
            if (i.ge.0) then
               l = l + 1
               if (l.ne.k) then
                  iswap = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
         end do 

c        if a vertex is required,  add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.

            cndmax = rtmax
            nzadd = nz

            do 80 iartif = 1, nzadd

               if (unitq) then
                  ifix = nfree
                  jadd = kx(ifix)
               else

                  rowmax = 0d0

                  do 60 i = 1, nfree
                     rnorm = dnrm2(nz,q(i,1),ldq)
                     if (rowmax.lt.rnorm) then
                        rowmax = rnorm
                        ifix = i
                     end if
   60             continue

                  jadd = kx(ifix)

                  call e04nfr(unitq,inform,ifix,iadd,jadd,it,
     *                        nactiv,nz,nfree,ngq,n,lda,ldq,ldt,
     *                        kx,cndmax,a,t,gqm,q,w,c,s)
               end if

               nfree = nfree - 1
               nz = nz - 1
               nartif = nartif + 1
               istate(jadd) = 4

   80       continue

         if (it.gt.1) then

c           if some dependent constraints were rejected,  move  t  to
c           the top of the array  t.

            do 120 k = 1, nactiv

               j = nz + k

               do 100 i = 1, k
                  t(i,j) = t(it+i-1,j)
  100          continue

  120       continue
         end if
      end if

      nrejtd = k2 - nactiv

      end

      subroutine lpinit (nctotl,nactiv,nartif,nfree,n,
     *                   istate,kactiv,kx,x,wx)

      integer is, j, jfix, jfree, n, nactiv, nartif, nctotl, nfree,
     *        istate(nctotl), kactiv(n), kx(n)

      double precision  wx(n), x(n)

c     initialize variables 

      do j = 1, n
         x(j) = 0d0
         wx(j) = 0d0
      end do

      nfree = n
      nactiv = 0
      nartif = 0

      do j = 1, n
         istate(j) = 0
      end do 

      do j = n+1, nctotl
         istate(j) = 3
      end do

c     define nfree and kactiv.

      do j = n+1, nctotl
         nactiv = nactiv + 1
         kactiv(nactiv) = j - n
      end do 
 
c        see if any bounds are violated or nearly satisfied.
c        if so,  add these bounds to the working set and set the
c        variables exactly on their bounds.

         j = n
100      if (j.ge.1.and.nactiv.lt.nfree) then

            if (istate(j).eq.0) then
               is = 0
               if (wx(j).le.1d-2) is = 1
               if (1d0-wx(j).le.2d-2) is = 2
               if (is.gt.0) then
                  istate(j) = is
                  if (is.eq.1) wx(j) = 0d0
                  if (is.eq.2) wx(j) = 1d0
                  nfree = nfree - 1
               end if
            end if

            j = j - 1
            go to 100

         end if

      jfree = 1
      jfix = nfree + 1

      do j = 1, n

         if (istate(j).le.0) then

            kx(jfree) = j
            jfree = jfree + 1

         else

            kx(jfix) = j
            jfix = jfix + 1

         end if

      end do 

      end

      subroutine e04mfr(job,n,nclin,nmoved,iter,numinf,istate,
     *                  featol,featlu,x)

      implicit none
      double precision  point6
      parameter         (point6=0.6d+0)
      integer           iter, n, nclin, nmoved, numinf
      character*1       job
      double precision  featlu(n+nclin), featol(n+nclin), x(n)
      integer           istate(n+nclin)
      double precision  tolx0
      integer           itnfix, ndegen
      double precision  wmach(9)
      integer           nfix(2)
      double precision  d, tolz
      integer           is, j, maxfix
      common            /ax02za/wmach/ce04mf/tolx0, ndegen, itnfix, nfix
      save              tolz

      nmoved = 0

      if (job.eq.'i') then

c        job = 'initialize'.
c        initialize at the start of each linear problem.
c        kdegen  is the expand frequency      and
c        featlu  are the user-supplied feasibility tolerances.
c        they are not changed.

         ndegen = 0
         itnfix = 0
         nfix(1) = 0
         nfix(2) = 0
         tolx0 = 0.5d0
         tolz = wmach(3)**point6

         do j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
         end do 

      else

c        job = 'end of cycle' or 'optimal'.
c        initialize local variables maxfix and tolz.

         maxfix = 2

         if (job.eq.'o') then

c           job = 'optimal'.
c           return with nmoved = 0 if the last call was at the same
c           iteration,  or if there have already been maxfix calls with
c           the same state of feasibility.

            if (itnfix.eq.iter) return
            if (numinf.gt.0) then
               j = 1
            else
               j = 2
            end if

            if (nfix(j).ge.maxfix) return
            nfix(j) = nfix(j) + 1
         end if

c        reset featol to its minimum value.

         do j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
         end do 

c        count the number of times a variable is moved a nontrivial
c        distance onto its bound.

         itnfix = iter

         do j = 1, n
            is = istate(j)
            if (is.gt.0.and.is.lt.4) then
               if (is.eq.1) then
                  d = dabs(x(j))
               else
                  d = dabs(x(j)-1d0)
               end if

               if (d.gt.tolz) nmoved = nmoved + 1
            end if
         end do 

      end if

      end

      subroutine e04mfq(n,nclin,istate,nviol,jmax,errmax,ax,bl,
     *                  featol,x)

      implicit none

c     e04mfq  checks the residuals of the constraints that are believed
c     to be feasible.  the number of constraints violated by more than
c     featol is computed, along with the maximum constraint violation.

      double precision  errmax
      integer           jmax, n, nclin, nviol

      double precision  ax(*), bl(nclin), featol(n+nclin), x(n)
      integer           istate(n+nclin)

      double precision  con, feasj, res
      integer           is, j


c     compute the number of constraints (nviol) violated by more than
c     featol and  the maximum constraint violation (errmax).
c     (the residual of a constraint in the working set is treated as if
c     it were an equality constraint fixed at that bound.)

      nviol = 0
      jmax = 0
      errmax = 0d0

      do j = 1, n + nclin
         is = istate(j)

         if (is.ge.0) then
            feasj = featol(j)

            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if
c           check for constraint violations.

            if (j.le.n) then 
               res = -con
               if (res.gt.feasj) then
                  nviol = nviol + 1
                  go to 20
               end if

               res = 1d0 - con
               if (res.lt.(-feasj)) then
                  nviol = nviol + 1
                  res = -res
                  go to 20
               end if
c           this constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.
               res = 0d0

               if (is.eq.1) then
               res = dabs(con)

               else if (is.eq.2) then
               res = dabs(1d0-con)

               else if (is.eq.3) then
               res = dabs(1d0-con)
               end if
            else

               res = bl(j-n) - con

               if (res.gt.feasj) then
                  nviol = nviol + 1
                  go to 20
               end if

c           this constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.

               res = 0d0

               if (is.eq.1) then
               res = dabs(bl(j-n)-con)

               else if (is.eq.2) then
               res = dabs(bl(j-n)-con)

               else if (is.eq.3) then
               res = dabs(bl(j-n)-con)
               end if
            end if 

            if (res.gt.feasj) nviol = nviol + 1

   20       if (res.gt.errmax) then
               jmax = j
               errmax = res
            end if

         end if

      end do 

      end

      subroutine e04mfs(firstv,n,nclin,istate,bigalf,pnorm,
     *                  hitlow,move,onbnd,unbndd,alfa,alfap,jhit,anorm,
     *                  ap,ax,bl,featol,p,x)

      implicit none 

c     e04mfs  finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).

c     in this version of e04mfs, when x is infeasible, the number of
c     infeasibilities will never increase.  if the number stays the
c     same, the sum of infeasibilities will decrease.  if the number
c     decreases by one or more,  the sum of infeasibilities will usually
c     decrease also, but occasionally it will increase after the step
c     alfa  is taken.  (convergence is still assured because the number
c     has decreased.)

c     three possible steps are computed as follows:

c     alfaf = the maximum step that can be taken without violating
c              one of the constraints that are currently satisfied.

c     alfai = reaches a linear constraint that is currently violated.
c              usually this will be the furthest such constraint along
c              p, subject to the angle between the constraint normal and
c              p being reasonably close to the maximum value among
c              infeasible constraints,  but if firstv = .true. it will
c              be the first one along p.  the latter case applies only
c              when the problem has been determined to be infeasible,
c              and the sum of infeasibilities are being minimized.
c              (alfai is not defined when x is feasible.)

c     alfai is needed occasionally when infeasible, to prevent
c     going unnecessarily far when alfaf is quite large.  it will
c     always come into effect when x is about to become feasible.
c     (the sum of infeasibilities will decrease initially as alfa
c     increases from zero, but may start increasing for larger steps.
c     choosing a large alfai allows several elements of  x  to
c     become feasible at the same time.

c     in the end, we take  alfa = alfaf  if x is feasible, or if
c     alfai > alfap (where  alfap  is the perturbed step from pass 1).
c     otherwise,  we take  alfa = alfai.

c     input parameters
c     ----------------
c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            if alfa reaches a bound as large as  it is
c            classed as an unbounded step.
c     featol is the array of current feasibility tolerances used by
c            e04mfh.  typically in the range 0.5*tolx to 0.99*tolx,
c            where tolx is the featol specified by the user.
c     tolinc (in common) is used to determine stepmn (see below),
c            the minimum positive step.
c     istate is set as follows:
c            istate(j) = -2  if a'x .lt. bl - featol
c                      = -1  if a'x .gt. bu + featol
c                      =  0  if a'x is not in the working set
c                      =  1  if a'x is in the working set at bl
c                      =  2  if a'x is in the working set at bu
c                      =  3  if a'x is in the working set (an equality)
c                      =  4  if x(j) is temporarily fixed.
c            values -2 and -1 do not occur once feasible.
c     bl     the lower bounds on the variables.
c     x      the values of       ditto.
c     p      the search direction.


c     output parameters
c     -----------------
c     hitlow  = true  if a lower bound restricted alfa.
c             = false otherwise.
c     move    = true  if  exact ge stepmn  (defined at end of code).
c     onbnd   = true  if  alfa = exact.  this means that the step  alfa
c                     moves x  exactly onto one of its constraints,
c                     namely  bound.
c             = false if the exact step would be too small
c                     (exact .lt. stepmn).
c               (with these definitions,  move = onbnd).
c     unbndd  = true  if alfa = bigalf.  jhit may possibly be zero.
c               the parameters hitlow, move, onbnd, bound and exact
c               should not be used.
c     jhit    = the index (if any) such that constraint jhit reaches
c               a bound.
c     bound   = the bound value bl(jhit) or bu(jhit) corresponding
c               to hitlow.
c     exact   = the step that would take constraint jhit exactly onto
c               bound.
c     alfa    = an allowable, positive step.
c               if unbndd is true,  alfa = stepmx.
c               otherwise,          alfa = max(stepmn, exact).


c     e04mfs is based on minos 5.2 routine m5chzr, which implements the
c     expand procedure to deal with degeneracy. the step alfaf is
c     chosen as in the two-pass approach of paula harris (1973), except
c     that this version insists on returning a positive step, alfa.
c     two features make this possible:

c        1. featol increases slightly each iteration.

c        2. the blocking constraint, when added to the working set,
c           retains the value ax(jhit) + alfa * ap(jhit),
c           even if this is not exactly on the blocking bound.

c     for infeasible variables moving towards their bound, we require
c     the rate of change of the chosen constraint to be at least gamma
c     times as large as the biggest available.  this still gives us
c     freedom in pass 2.
c     gamma = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
c     gamma = 0.001 seems to be safe.

      double precision  gamma
      parameter         (gamma=1d-3)

      double precision  alfa, alfap, bigalf, pnorm
      integer           jhit, n, nclin
      logical           firstv, hitlow, move, onbnd, unbndd

      double precision  anorm(*), ap(*), ax(*), bl(nclin),
     *                  featol(n+nclin), p(n), x(n)
      integer           istate(n+nclin)

      double precision  tolx0
      integer           itnfix, ndegen
      integer           nfix(2)
      double precision  alfai, atp, atpabs, atpmxf, atpmxi, atpscd, atx,
     *                  bound, delta, exact, res,
     *                  stepmn, tolpiv
      integer           i, j, jhitf, jhiti, js
      logical           blockf, blocki

      common            /ce04mf/tolx0, ndegen, itnfix, nfix

      double precision wmach(9)
      common/ax02za/wmach

c     tolpiv is a tolerance to exclude negligible elements of a'p.

      tolpiv = wmach(1)*pnorm

c     first pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

      atpmxi = 0d0
      alfap = bigalf

      do 20 j = 1, n + nclin
         js = istate(j)

         if (js.le.0) then
            delta = featol(j)

            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = dabs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = dabs(atp)
               atpscd = atpabs/(1d0+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

c              res = -one

            else if (atp.le.0d0.and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

c              test for smaller alfap.
c              if the upper bound is violated. test for bigger atp.

               if (j.le.n) then 
                   res = atx + delta                    
                  if (res.lt.alfap*atpabs) alfap = res/atpabs
               else                  
                  res = atx - bl(j-n) + delta 
                  if (res.lt.alfap*atpabs) alfap = res/atpabs
               end if 

               if (js.eq.-1) atpmxi = max(atpmxi,atpscd)

            else if (atp.gt.0d0.and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfap.
c              if the lower bound is violated. test for bigger atp.

               if (j.le.n) then 
                  res = 1d0 - atx + delta
                  if (res.lt.alfap*atp) alfap = res/atp
               else 
                  res = bl(j-n) - atx + delta
                  if (res.lt.alfap*atp) alfap = res/atp
               end if 

               if (js.eq.-2) atpmxi = max(atpmxi,atpscd)
            end if

         end if
   20 continue

c     second pass.
c     for feasible variables, recompute steps without perturbation.
c     amongst constraints that are closer than alfap, choose the one
c     that makes the largest angle with the search direction.
c     for infeasible variables, find the largest step subject to a'p
c     being no smaller than gamma * max(a'p).

      if (firstv) then
         alfai = bigalf
      else
         alfai = 0d0
      end if

      atpmxf = 0d0
      atpmxi = gamma*atpmxi
      jhitf = 0
      jhiti = 0

      do 40 j = 1, n + nclin
         js = istate(j)

         if (js.le.0) then

            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = dabs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = dabs(atp)
               atpscd = atpabs/(1d0+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

c              res = -one

            else if (atp.le.0d0.and. js.ne.-2) then

c              a'x  is decreasing.

c              test for bigger a'p if the lower bound is satisfied.
c              test for smaller alfaf.

               if (atpscd.gt.atpmxf) then

                  if (j.le.n) then 
                     res = atx 

                     if (res.le.alfap*atpabs) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                   else 
                     res = atx - bl(j-n)

                     if (res.le.alfap*atpabs) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                   end if 
               end if

               if (js.eq.-1) then
c                 the upper bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.
                  if (j.le.n) then 
                  if (firstv) then
                     res = atx - 1d0

                     if (res.le.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if

                  else if (atpscd.ge.atpmxi) then
                     res = atx - 1d0

                     if (res.gt.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
                  end if
                  else 
                  if (firstv) then
                     res = atx - bl(j-n)

                     if (res.le.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if

                  else if (atpscd.ge.atpmxi) then
                     res = atx - bl(j-n)

                     if (res.gt.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
                  end if
                  
                  end if 
               end if

            else if (atp.gt.0d0.and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfap.

               if (atpscd.gt.atpmxf) then

                  if (j.le.n) then 
                     res = 1d0 - atx
                     if (res.le.alfap*atp) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                  else 
                     res = bl(j-n) - atx
                     if (res.le.alfap*atp) then
                        atpmxf = atpscd
                        jhitf = j
                     end if

                  end if 
               end if

               if (js.eq.-2) then
c                 the lower bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.
               if (j.le.n) then       
                  if (firstv) then
                     res = -atx

                     if (res.le.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  else if (atpscd.ge.atpmxi) then
                     res = - atx

                     if (res.gt.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  end if       
               else
                  res = bl(j-n) - atx
                  if (firstv) then
                     if (res.le.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  else if (atpscd.ge.atpmxi) then
                     if (res.gt.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  end if
               end if
               end if 
            end if
         end if
   40 continue

c     ------------------------------------------------------------------
c     see if a feasible and/or infeasible constraint blocks.
c     ------------------------------------------------------------------
      blockf = jhitf .gt. 0
      blocki = jhiti .gt. 0
      unbndd = .not. (blockf .or. blocki)

      if (unbndd) go to 60

      if (blockf) then

c        a constraint is hit which is currently feasible.
c        the corresponding step alfaf is not used, so no need to get it,
c        but we know that alfaf .le. alfap, the step from pass 1.

         jhit = jhitf
         if (jhit.le.n) then
            atp = p(jhit)
         else
            atp = ap(jhit-n)
         end if
         hitlow = atp .lt. 0d0
      end if

c     if there is a choice between alfaf and alfai, it is probably best
c     to take alfai.  however, we can't if alfai is bigger than alfap.

      if (blocki.and.alfai.le.alfap) then

c        an infeasible variable reaches its violated bound.

         jhit = jhiti
         if (jhit.le.n) then
            atp = p(jhit)
         else
            atp = ap(jhit-n)
         end if
         hitlow = atp .gt. 0d0
      end if

      if (jhit.le.n) then
         atx = x(jhit)
      else
         atx = ax(jhit-n)
      end if

c     try to step exactly onto bound, but make sure the exact step
c     is sufficiently positive.  (exact will be alfaf or alfai.)
c     since featol increases by  tolinc  each iteration, we know that
c     a step as large as  stepmn  (below) will not cause any feasible
c     variables to become infeasible (where feasibility is measured
c     by the current featol).

      if (jhit.le.n) then 
      if (hitlow) then
         bound = 0d0
      else
         bound = 1d0
      end if
      else 
         bound = bl(jhit-n)
      end if 

      stepmn = 0d0
      exact = (bound-atx)/atp
      alfa = max(stepmn,exact)
      onbnd = alfa .eq. exact
      move = exact .ge. stepmn
      if (.not. move) ndegen = ndegen + 1

      return

c     unbounded.

   60 alfa = bigalf
      move = .true.
      onbnd = .false.

      end

      subroutine e04nfr(unitq,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,ngq,n,lda,ldq,ldt,kx,condmx,
     *                  a,t,gqm,q,w,c,s)

      implicit none 

c     e04nfr  updates the matrices  z, y, t, r  and  d  associated with
c     factorizations

c              a(free) * q(free)  = ( 0 t)
c                        q(free)  = ( z y)
c                      r' *d * r  =   hz

c     a) the matrices  r  and  t  are upper triangular.
c     b) the arrays  t  and  r  may be the same array.
c     c) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt) of the
c        array  t.   the integer  jt  is always  nz+1.  during regular
c        changes to the working set,  it = 1;  when several constraints
c        are added simultaneously,  it  points to the first row of the
c        existing  t.
c     d) the matrix  r  is stored in the first  nz x nz  rows
c        and columns of the  nfree x nfree  leading principal minor of
c        the array  r.


c     there are three separate cases to consider (although each case
c     shares code with another)...

c     (1) a free variable becomes fixed on one of its bounds when there
c         are already some general constraints in the working set.

c     (2) a free variable becomes fixed on one of its bounds when there
c         are only bound constraints in the working set.

c     (3) a general constraint (corresponding to row  iadd  of  a) is
c         added to the working set.

c     in cases (1) and (2), we assume that  kx(ifix) = jadd.
c     in all cases,  jadd  is the index of the constraint being added.

c     if there are no general constraints in the working set,  the
c     matrix  q = (z y)  is the identity and will not be touched.

c     if  ngq .gt. 0,  the column transformations are applied to the
c     columns of the  (ngq x n)  matrix  gqm'.

      double precision  condmx
      integer           iadd, ifix, inform, it, jadd, lda, ldq,
     *                  ldt, n, nactiv, nfree, ngq, nz
      logical           unitq
      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*),
     *                  s(n), t(ldt,*), w(n)
      integer           kx(n)
      double precision  asize, dtmax, dtmin
      double precision  cond, dtnew, tdtmax, tdtmin
      integer           i, j, jt, k, nanew, npiv
      logical           bound, overfl
      double precision  dnrm2, adivb

      common            /de04nb/asize, dtmax, dtmin

      double precision wmach(9)
      common/ax02za/wmach

      overfl = .false.
      bound = jadd .le. n
      jt = nz + 1

      if (bound) then
c        ===============================================================
c        a simple bound has entered the working set.  iadd is not used.
c        ===============================================================

         nanew = nactiv

         if (unitq) then

c           q is not stored, but  kx  defines an ordering of the columns
c           of the identity matrix that implicitly define q.
c           define the sequence of pairwise interchanges p that moves
c           the newly-fixed variable to position  nfree.
c           reorder  kx  accordingly.

            do i = 1, nfree - 1
               if (i.ge.ifix) then
                  w(i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
            end do 
         else

c           q  is stored explicitly.

c           set  w = the  (ifix)-th  row of  q.
c           move the  (nfree)-th  row of  q  to position ifix.

            call dcopy(nfree,q(ifix,1),ldq,w,1)
            if (ifix.lt.nfree) then
               call dcopy(nfree,q(nfree,1),ldq,q(ifix,1),ldq)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else
c        ===============================================================
c        a general constraint has entered the working set.
c        ifix is not used.
c        ===============================================================

         nanew = nactiv + 1

c        transform the incoming row of a by q'.

         call dcopy(n,a(iadd,1),lda,w,1)
         call e04nbw(8,n,nz,nfree,ldq,unitq,kx,w,q,c)

c        check that the incoming row is not dependent upon those
c        already in the working set.

         dtnew = dnrm2(nz,w,1)
         if (nactiv.eq.0) then

c           this is the only general constraint in the working set.

            cond = adivb(asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else

c           there are already some general constraints in the working
c           set.  update the estimate of the condition number.

            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = adivb(tdtmax,tdtmin,overfl)
         end if

         if (cond.gt.condmx .or. overfl) go to 80

         if (unitq) then

c           first general constraint added.  set  q = i.

            call f06qhf (nfree,nfree,q,ldq)
            unitq = .false.
            it = 0
         end if
      end if

      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if

      if (unitq) then

c        the orthogonal matrix  q  (i.e.,  q) is not stored explicitly.
c        apply  p, the sequence of pairwise interchanges that moves the
c        newly-fixed variable to position  nfree.

         if (ngq.gt.0) call f06qkf(nfree-1,w,ngq,gqm,n)
      else

         call f06fqf (npiv-1,w(npiv),w,c,s)

         if (ngq.gt.0) call f06qxf ('l','f',npiv,ngq,1,npiv,c,s,gqm,n)
         call f06qxf ('r','f',nfree,nfree,1,npiv,c,s,q,ldq)
      end if

      if (.not. unitq) then
         if (bound) then

c           bound constraint added.   the rotations affect columns
c           nz+1  thru  nfree  of  gqm'  and  t.

c           the last row and column of  q  has been transformed to plus
c           or minus the unit vector  e(nfree).  we can reconstitute the
c           column of gqm' corresponding to the new fixed variable.

            if (w(nfree).lt.wmach(3).and.ngq.gt.0) 
     *          call dscal(ngq,-1d0,gqm(nfree,1),n)


            if (nactiv.gt.0) then

               t(it,jt-1) = s(jt-1)*t(it,jt)

               if (dabs(c(jt-1)).lt.wmach(3)) c(jt-1) = 0d0

               t(it,jt) = c(jt-1)*t(it,jt)

               if (nactiv.gt.1) then
                  call f06qvf('right',nactiv,1,nactiv,c(jt),s(jt),
     *                        t(it,jt),ldt)
                  call dcopy(nactiv-1,s(jt),1,t(it+1,jt),ldt+1)
               end if

               jt = jt - 1
               call f06flf(nactiv,t(it,jt),ldt+1,tdtmax,tdtmin)
               cond = adivb(tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint added.  install  w  at the front of  t.
c           if there is no room,  shift all the rows down one position.

            it = it - 1
            if (it.le.0) then
               it = 1
               do 60 k = 1, nactiv
                  j = jt + k - 1
                  do 40 i = k, 1, -1
                     t(i+1,j) = t(i,j)
   40             continue
   60          continue
            end if
            jt = jt - 1
            call dcopy (nanew,w(jt),1,t(it,jt),ldt)
         end if
      end if

c     prepare to exit.  check the magnitude of the condition estimator.

   80 if (nanew.gt.0) then
         if (cond.lt.condmx.and..not. overfl) then

c           the factorization has been successfully updated.

            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin

         else

c           the proposed working set appears to be linearly dependent.

            inform = 1
         end if
      end if

      end


      subroutine e04nfp(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,t,gqm,q,c,s,fail)
      implicit none 

      integer           it, jdel, kdel, lda, ldq, ldt, n, nactiv, nfree,
     *                  ngq, nrz, nz, kactiv(n), kx(n)
      logical           unitq,fail
      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n), 
     *                  t(ldt,*), asize, dtmax, dtmin, cs, sn
      integer           i, ir, itdel, j, jart, jt, k, l, npiv, nrz1,
     *                  nsup, idamax
      common            /de04nb/asize, dtmax, dtmin

      jt = nz + 1
      fail = .false.

      if (jdel.gt.0) then

c        regular constraint or temporary bound deleted.

         if (jdel.le.n) then

c           case 1.  a simple bound has been deleted.
c           =======  columns  nfree+1  and  ir  of gqm' must be swapped.

            ir = nz + kdel

            itdel = nactiv + 1
            nfree = nfree + 1
            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               call dswap(ngq,gqm(nfree,1),n,gqm(ir,1),n)
            end if

            if (.not. unitq) then

c              copy the incoming column of  a(free)  into the end of  t.

               do 20 k = 1, nactiv
                  i = kactiv(k)
                  t(nactiv-k+1,nfree) = a(i,jdel)
   20          continue

c              expand  q  by adding a unit row and column.

               if (nfree.gt.1) then

                  if (nfree.gt.ldq) then 
c                                 dec 07 bug? JADC
c                    call warn (999,zero,nfree,'NLIB')
                     fail = .true.
                     return

                  end if 

                  call f06fbf (nfree-1,0d0,q(nfree,1),ldq)

                  call f06fbf (nfree-1,0d0,q(1,nfree),1)

               end if
               q(nfree,nfree) = 1d0
            end if
         else
c           case 2.  a general constraint has been deleted.

c           delete row  itdel  of  t  and move up the ones below it.
c           t  becomes lower hessenberg.

            itdel = kdel
            do 60 k = itdel, nactiv
               j = jt + k - 1
               do 40 l = itdel, k - 1
                  i = it + l - 1
                  t(i,j) = t(i+1,j)
   40          continue
   60       continue

            do 80 i = nactiv - itdel + 1, nactiv - 1
               kactiv(i) = kactiv(i+1)
   80       continue

            nactiv = nactiv - 1

         end if

         nz = nz + 1

         if (nactiv.eq.0) then

            dtmax = 1d0
            dtmin = 1d0

         else

            nsup = itdel - 1

            if (nsup.gt.0) then

               npiv = jt + itdel - 1

               if (nsup.gt.1) then

                  call dcopy(nsup-1,t(it+1,jt+1),ldt+1,s(jt+1),1)

                  call f06qrf('right',nactiv,1,nsup,c(jt+1),s(jt+1),
     *                        t(it,jt+1),ldt)

               end if

               call f06baf(t(it,jt+1),t(it,jt),cs,sn)

               t(it,jt) = 0d0
               s(jt) = -sn
               c(jt) = cs

               call f06qxf('r','b',nfree,nfree,nz,npiv,c,s,q,ldq)
               call f06qxf('l','b',npiv,ngq,nz,npiv,c,s,gqm,n)

            end if

            jt = jt + 1
            call f06flf(nactiv,t(it,jt),ldt+1,dtmax,dtmin)
         end if
      end if

      nrz1 = nrz + 1

      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamax(nz-nrz1+1,gqm(nrz1,1))
         else
            jart = -jdel
         end if

         if (jart.gt.nrz1) then

            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,q(1,nrz1),1,q(1,jart),1)
            end if

            call dswap(ngq,gqm(nrz1,1),n,gqm(jart,1),n)
         end if
      end if

      nrz = nrz1

      end

      subroutine e04mfl(n,nrz,nz,zerolm,notopt,numinf,trusml,
     *                  smllst,jsmlst,tinyst,jtiny,gq)
      implicit none
      double precision  smllst, tinyst, trusml, zerolm
      integer jsmlst, jtiny, n, notopt, nrz, numinf, nz, j
      double precision  gq(n), rlam

      do 20 j = nrz + 1, nz
         rlam = -abs(gq(j))

         if (rlam.lt.zerolm) then
            if (numinf.eq.0) notopt = notopt + 1

            if (rlam.lt.smllst) then
               trusml = gq(j)
               smllst = rlam
               jsmlst = -j
            end if

         else if (rlam.lt.tinyst) then
            tinyst = rlam
            jtiny = -j
         end if
   20 continue

      end

      subroutine e04mfm(n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,anorms,gq,rlamda,t,wtinf)
c----------------------------------------------------------------------
      implicit none 

      double precision  biggst, smllst, tinyst, trubig, trusml, zerolm
      integer           jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, n, nactiv, nfree, notopt,
     *                  numinf, nz
      double precision  a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer           istate(*), kactiv(n), kx(n)
      double precision  anormj, blam, rlam, scdlam
      integer           i, is, j, k, l, nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree
      jtiny = 0
      jsmlst = 0
      ksmlst = 0
      jbigst = 0
      kbigst = 0

      if (n.gt.nz) call dcopy(n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call dtrsv('t',nactiv,t(1,nz+1),ldt,rlamda)

      do 40 l = 1, nfixed
         j = kx(nfree+l)
         blam = rlamda(nactiv+l)
         do 20 k = 1, nactiv
            i = kactiv(k)
            blam = blam - a(i,j)*rlamda(nactiv-k+1)
   20    continue
         rlamda(nactiv+l) = blam
   40 continue

      do 60 k = 1, n - nz
         if (k.gt.nactiv) then
            j = kx(nz+k)
         else
            j = kactiv(nactiv-k+1) + n
         end if

         is = istate(j)

         i = j - n
         if (j.le.n) anormj = 1d0
         if (j.gt.n) anormj = anorms(i)

         rlam = rlamda(k)

c        change the sign of the estimate if the constraint is in
c        the working set at its upper bound.

         if (is.eq.2) rlam = -rlam
         if (is.eq.3) rlam = dabs(rlam)
         if (is.eq.4) rlam = -dabs(rlam)

         if (is.ne.3) then
            scdlam = rlam*anormj

            if (scdlam.lt.zerolm) then
               if (numinf.eq.0) notopt = notopt + 1

               if (scdlam.lt.smllst) then
                  smllst = scdlam
                  trusml = rlamda(k)
                  jsmlst = j
                  ksmlst = k
               end if
            else if (scdlam.lt.tinyst) then
               tinyst = scdlam
               jtiny = j
            end if
         end if

         scdlam = rlam/wtinf(j)
         if (scdlam.gt.biggst .and.j.gt.jinf) then
            biggst = scdlam
            trubig = rlamda(k)
            jbigst = j
            kbigst = k
         end if
   60 continue

      end

      subroutine e04mfh(n,nclin,lda,istate,numinf,suminf,bl,a,
     *                  featol,cvec,x,wtinf)
c----------------------------------------------------------------------
      implicit none

      double precision  suminf
      integer           lda, n, nclin, numinf

      double precision  a(lda,*), bl(*), cvec(n), featol(*),
     *                  wtinf(*), x(n)
      integer           istate(*)
      double precision  ctx, feasj, s, weight
      integer           j, k
      double precision  ddot
c----------------------------------------------------------------------
      numinf = 0
      suminf = 0d0

      call f06fbf(n,0d0,cvec,1)

      do 40 j = 1, n + nclin

         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot(n,a(k,1),lda,x)
            end if
            istate(j) = 0

c           see if the lower bound is violated.
            if (j.le.n) then      
               s = 0d0 - ctx
               if (s.gt.feasj) then
                  istate(j) = -2
                  weight = -wtinf(j)
                  go to 20
               end if
c           see if the upper bound is violated.

               s = ctx - 1d0
               if (s.le.feasj) go to 40
               istate(j) = -1
               weight = wtinf(j)
            else
               s = bl(j-n) - ctx
               if (s.gt.feasj) then
                  istate(j) = -2
                  weight = -wtinf(j)
                  go to 20
               end if

               if (-s.le.feasj) go to 40
               istate(j) = -1
               weight = wtinf(j)

            end if 

c           add the infeasibility.

   20       numinf = numinf + 1
            suminf = suminf + dabs(weight)*s
            if (j.le.n) then
               cvec(j) = weight
            else
               call daxpy(n,weight,a(k,1),lda,cvec)
            end if
         end if
   40 continue

      end

      subroutine f06baf(a, b, c, s)
c----------------------------------------------------------------------
      double precision   a, b, c, s

      double precision   t
      logical            fail
      double precision   adivb
c----------------------------------------------------------------------
      if (b.eq.0d0) then
         c  = 1d0
         s  = 0d0
      else
         t  = adivb (b,a,fail)
         call f06bcf(t, c, s)
         a  = c*a + s*b
         b  = t
      end if

      end

      subroutine f06fqf (n, alpha, x, c, s)
      implicit none
      double precision   c(*), s(*), x(*), alpha
      integer            n, i, ix

      if (n.gt.0) then

            ix = 1

         do i = 1, n - 1
                  call f06baf(x(ix + 1), x(ix), c(i), s(i))
                  s(i)  = -s(i)
                  x(ix) = -x(ix)
                  ix      =  ix + 1
         end do 
               call f06baf(alpha, x(ix), c(n), s(n))
               s(n)  = -s(n)
               x(ix) = -x(ix)

      end if

      end

      subroutine f06qrf(side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
      implicit none

      integer  k1, k2, lda, n, i, j

      character*1 side

      double precision a(lda,*), c(*), s(*), aij, ctemp, stemp, subh, 
     *                 temp

      double precision wmach(9)
      common/ax02za/wmach
c----------------------------------------------------------------------
      if ((min(n, k1).lt.1).or.(k2.le.k1).or.(k2.gt.n))return

      if (side.eq.'l') then

         do 20 j = k1, n
            aij = a(k1,j)
            do 10 i = k1, min(j, k2) - 1
               temp = a(i + 1,j)
               a(i,j) = s(i)*temp + c(i)*aij
               aij = c(i)*temp - s(i)*aij
   10       continue
            if (j.lt.k2) then

               subh = s(j)
               call f06baf(aij, subh, c(j), s(j))
               a(j,j) = aij
            else
               a(k2,j) = aij
            end if
   20    continue

      else if (side.eq.'r') then

         do 40 j = k2 - 1, k1, -1

            subh = s(j)

            call f06baf(a(j + 1, j + 1), subh, ctemp, stemp)

            stemp = -stemp

            if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

            s(j) = stemp
            c(j) = ctemp

            if ((ctemp.ne.1d0).or.(stemp.ne.0d0)) then

               do 30 i = j, 1, -1
                  temp = a(i, j + 1)
                  a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   30          continue

            end if

   40    continue

      end if

      end

      subroutine f06qvf(side, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, n, i, j

      character*1        side

      double precision   a(lda, *), c(*), s(*), aij, ctemp, stemp, temp

      double precision wmach(9)
      common /ax02za/wmach
c----------------------------------------------------------------------
      if (min(n,k1).lt.1.or.k2.le.k1.or.k2.gt.n) return

      if (side.eq.'l') then

         do 20 j = n, k1, -1

            if (j.ge.k2) then
               aij = a(k2,j)
            else
               aij = c(j)*a(j,j)
               s(j) = -s(j)*a(j,j)
            end if

            do 10 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i + 1,j) = c(i)*aij - s(i)*temp
               aij = s(i)*aij + c(i)*temp
   10       continue
            a(k1,j) = aij
   20    continue

      else if (side.eq.'r') then

         do j = k1, k2 - 1

            if (c(j).ne.1d0.or.s(j).ne.0d0) then
               stemp = s(j)
               ctemp = c(j)

               if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

               do i = 1, j

                  temp = a(i, j + 1)
c                                        could check for underflow.
c                                 added following line 11/06
                  if (dabs(a(i,j)).lt.wmach(3)) a(i,j) = 0d0                            
                  a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)

               end do

               s(j) = stemp*a(j + 1, j + 1)
               a(j + 1, j + 1) = ctemp*a(j + 1, j + 1)

            end if
         end do 
      end if

      end

      subroutine f06qxf(side, direct, m, n, k1, k2, c, s, a, lda)
c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, m, n, i, j

      character*1        direct, side

      double precision   a(lda, *), c(*), s(*), aij, ctemp, stemp, temp

      logical            left, right

      double precision wmach(9)
      common/ax02za/wmach
c----------------------------------------------------------------------
      left = (side.eq.'l')
      right = (side.eq.'r')

      if ((min(m, n, k1).lt.1).or.(k2.le.k1).or.
     $    (left .and.k2.gt.m).or.(right.and.k2.gt.n)) return

      if (left) then

            if (direct.eq.'f') then

               do j = 1, n
                  aij = a(k1,j)
                  do i = k1, k2 - 1
                     temp = a(i + 1,j)

                     if (c(i).lt.wmach(3)) c(i) = 0d0

                     a(i,j) = s(i)*temp + c(i)*aij
                     aij = c(i)*temp - s(i)*aij
                  end do 
                  a(k2,j) = aij
               end do

            else if (direct.eq.'b') then

               do j = 1, n
                  aij = a(k2,j)

                  do i = k2 - 1, k1, -1
                     temp = a(i,j)
                     if (c(i).lt.wmach(3)) c(i) = 0d0
                     a(i + 1,j) = c(i)*aij - s(i)*temp
                     aij = s(i)*aij + c(i)*temp
                  end do

                  a(k1,j) = aij
               end do

            end if

      else if (right) then

            if (direct.eq.'f') then

               do 140 j = k1, k2 - 1

                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then
                     ctemp = c(j)

                     if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

                     stemp = s(j)

                     do 130 i = 1, m
                        temp = a(i, j + 1)
                        a(i, j + 1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
  130                continue

                  end if

  140          continue

            else if (direct.eq.'b') then

               do j = k2 - 1, k1, -1

                  if ((c(j).ne.1d0).or.(s(j).ne.0d0)) then

                     ctemp = c(j)

                     if (dabs(ctemp).lt.wmach(3)) ctemp = 0d0

                     stemp = s(j)

                     do i = m, 1, -1
                        temp = a(i,j + 1)
                        a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                        a(i,j) = stemp*temp + ctemp*a(i,j)
                     end do

                  end if

               end do

            end if

      end if

      end


      subroutine f06qkf (n, perm, k, b, ldb)
c----------------------------------------------------------------------
      implicit none

      integer            k, ldb, n, i, j, l

      double precision   perm(*), b(ldb, *), temp
c----------------------------------------------------------------------
      if (min(n, k).eq.0) return

            do i = 1, n

               l = idint(perm(i))

               if (l.ne.i) then

                  do j = 1, k
                     temp = b(i,j)
                     b(i,j) = b(l,j)
                     b(l,j) = temp
                  end do

               end if

             end do

      end

      subroutine f06frf(n, alpha, x, zeta)
c----------------------------------------------------------------------
      implicit none

      double precision   alpha, zeta, beta, scale, ssq, x(*)

      integer n

      double precision wmach(9)
      common /ax02za/wmach
c----------------------------------------------------------------------
      if (n.lt.1) then

         zeta = 0d0

      else if (n.eq.1.and.x(1).eq.0d0) then

         zeta = 0d0

      else
c        treat case where p is a 2 by 2 matrix specially.

         if (n.eq.1) then

c           deal with cases where  alpha = zero  and
c           dabs(x(1)) .le. max(eps*abs(alpha), tol)  first.

            if (alpha.eq.0d0) then
               zeta   =  1d0
               alpha  =  dabs (x(1))
               x(1) = -dsign(1d0, x(1))
            else if (dabs(x(1)).le.max(wmach(3)*dabs(alpha),0d0)) then
               zeta   =  0d0
            else

               if (dabs(alpha).ge.dabs(x(1))) then
                  beta = dabs(alpha)*dsqrt(1d0 + (x(1)/alpha)**2)
               else
                  beta = dabs(x(1))*dsqrt(1d0 + (alpha/x(1))**2)
               end if

               zeta = dsqrt((dabs(alpha) + beta)/beta)
               if (alpha.ge.0d0) beta = -beta
               x(1) = -x(1)/(zeta*beta)
               alpha  = beta

            end if
         else

c           now p is larger than 2 by 2.

            ssq   = 1d1
            scale = 0d0

            call f06fjf(n, x, 1, scale, ssq)

c           treat cases where  scale = zero,
c           scale .le. max(eps*abs(alpha), tol)  and
c           alpha = zero  specially.
c           note that  scale = max(abs(x(i))).

            if ((scale.eq.0d0).or.
     $          (scale.le.max(wmach(3)*dabs(alpha),0d0))) then

               zeta  = 0d0

            else if (alpha.eq.0d0) then

               zeta  = 1d0
               alpha = scale*dsqrt(ssq)
               call dscal(n, -1d0/alpha, x, 1)

            else

               if (scale.lt.dabs(alpha)) then
                  beta = dabs(alpha)*dsqrt(1d0 + ssq*(scale/alpha)**2)
               else
                  beta = scale       *dsqrt(ssq +   (alpha/scale)**2)
               end if

               zeta = dsqrt((beta + dabs(alpha))/beta)
               if (alpha.gt.0d0) beta = -beta
               call dscal(n,-1d0/(zeta*beta), x, 1)
               alpha = beta

            end if

         end if

      end if

      end

      subroutine f06bcf (t, c, s)
c----------------------------------------------------------------------
      implicit none

      double precision   c, s, t, abst

      double precision wmach(9)
      common /ax02za/wmach
c----------------------------------------------------------------------
      abst = dabs(t)

      if (abst.lt.wmach(4)) then

         c = 1d0
         s = t

      else if (abst*wmach(4).gt.1d0) then

         c = 1d0/abst
         s = sign (1d0,t)

      else

         c = 1d0 / dsqrt(1d0 + abst**2)
         s = c*t

      end if

      end


      subroutine f06flf(n, x, incx, xmax, xmin)
c----------------------------------------------------------------------
      implicit none

      double precision   xmax, xmin, x(*)
      integer            incx, n, ix
c----------------------------------------------------------------------
      if (n.lt.1) then

         xmax = 0d0
         xmin = 0d0

      else

         xmax = dabs(x(1))
         xmin = xmax

         do ix = 1 + incx, 1 + (n - 1)*incx, incx
            xmax = max(xmax,dabs(x(ix)))
            xmin = min(xmin,dabs(x(ix)))
         end do

      end if

      end

      subroutine f06fbf(n, const, x, incx)
c----------------------------------------------------------------------
      implicit none

      double precision   const, x(*)

      integer            incx, n, ix
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (const.ne.0d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = const
            end do

         else

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = 0d0
            end do

         end if

      end if

      end

      subroutine f06fjf(n, x, incx, scale, sumsq)
c----------------------------------------------------------------------
      implicit none

      integer            incx, n, ix

      double precision   scale, sumsq, x(*), absxi
c----------------------------------------------------------------------
      if (n.gt.0) then

         do ix = 1, 1 + (n-1)*incx, incx

            if (x(ix).ne.0d0) then

               absxi = dabs(x(ix))

               if (scale.lt.absxi) then
                  sumsq = 1d0 + sumsq*(scale/absxi)**2
                  scale = absxi
               else
                  sumsq = sumsq + (absxi/scale)**2
               end if

            end if

         end do

      end if

      end

      double precision function adivb (a,b,fail)
c----------------------------------------------------------------------
      implicit none

      logical fail

      double precision   a, b, absb, absa, div

      double precision wmach(9)
      common /ax02za/wmach
c----------------------------------------------------------------------
      absa = dabs(a)
      absb = dabs(b)

      if (absa.le.wmach(3)) then
c                                 a is < machine eps
         div = 0d0

         if (absb.le.wmach(3)) then
c                                 b is also < machine eps
            fail = .true.
         else
            fail = .false.
         end if

      else
c                                 a is > eps
         if (absb.le.wmach(3)) then
c                                 b is < eps, div = sign(a)*huge
            div  =  dsign(wmach(7),a)
            fail = .true.

         else
c                                 b is > eps
            if (absb.ge.1d0) then
c                                 |b| > 1
               fail = .false.

               if (absa.ge.absb*wmach(3)) then
                  div = a/b
               else
                  div = 0d0
               end if

            else
c                                 |b| < 1
               if (absa.le.absb*wmach(7)) then

                  fail = .false.
                  div  =  a/b

               else

                  fail = .true.
                  div  = wmach(7)

                  if (((a.lt.0d0).and.(b.gt.0d0)).or.
     *                ((a.gt.0d0).and.(b.lt.0d0))) div = -div

               end if

            end if

         end if

      end if

      adivb = div

      end


      subroutine f06qhf (m, n, a, lda)
c----------------------------------------------------------------------
      implicit none

      integer i, j, lda, m, n

      double precision a(lda, *)
c----------------------------------------------------------------------
         do j = 1, n
            do i = 1, m
               a(i,j) = 0d0
            end do
         end do

         do i = 1, min(m,n)
            a(i,i) = 1d0
         end do

      end

      double precision function dnrm2 (n, x, incx)
c----------------------------------------------------------------------
      implicit none

      integer incx, n

      double precision x(*), norm, scale, ssq, f06bmf
c----------------------------------------------------------------------
      if (n.lt.1) then

         norm  = 0d0

      else if (n.eq.1) then

         norm  = dabs(x(1))

      else

         scale = 0d0 
         ssq   = 1d0

         call f06fjf(n, x, incx, scale, ssq)

         norm  = f06bmf(scale, ssq)

      end if

      dnrm2 = norm

      end

      subroutine dscal (n, alpha, x, incx)
c----------------------------------------------------------------------
c     x := alpha*x
c----------------------------------------------------------------------
      implicit none

      double precision   alpha, x(*)

      integer            incx, n, ix
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (alpha.eq.0d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = 0d0
            end do

         else if (alpha.ne.1d0) then

            do ix = 1, 1 + (n - 1)*incx, incx
               x(ix) = alpha*x(ix)
            end do 

         end if

      end if

      end


      double precision function ddot  (n, x, incx, y)
c----------------------------------------------------------------------
      implicit none

      integer   i, ix, iy, incx, n

      double precision x(*), y(*), sum
c----------------------------------------------------------------------
      sum = 0d0

      if (n.gt.0) then

         if ((incx.eq.1).and.(incx.gt.0)) then

            do ix = 1, 1 + (n - 1)*incx, incx
               sum = sum + x(ix)*y(ix)
            end do

         else

            iy = 1    
            if (incx.gt.0) then

               do ix = 1, 1 + (n - 1)*incx, incx
                  sum = sum + x(ix)*y(iy)
                  iy  = iy  + 1
               end do

            else

               ix = 1 - (n - 1)*incx

               do i = 1, n
                  sum = sum + x(ix)*y(iy)
                  ix  = ix  + incx
                  iy  = iy  + 1
               end do 

            end if
         end if
      end if

      ddot = sum

      end

      double precision function f06bmf (scale, ssq)
c----------------------------------------------------------------------
      implicit none

      double precision scale, ssq, sqt, norm

      double precision wmach(9)
      common /ax02za/wmach
c----------------------------------------------------------------------
      sqt = dsqrt(ssq)

      if (scale.lt.wmach(7)/sqt) then
         norm = scale*sqt
      else
         norm = wmach(7)
      end if

      f06bmf = norm

      end

      subroutine dgemv (trans, m, n, a, lda, x, beta, y)
c----------------------------------------------------------------------
      implicit none

      integer lda, m, n, i, info, j, jx, kx, leny, m4, n4

      character*1 trans

      double precision a(lda, *), x(*), y(*), beta, temp, temp1, 
     *                 temp2, temp3, temp4      
c----------------------------------------------------------------------

      info = 0

      if (m.lt.0) then
         info = 2
      else if (n.lt.0) then
         info = 3
      else if (lda.lt.max(1, m)) then
         info = 6
      end if
      if (info.ne.0) then
         return
      end if

      if (m.eq.0.or.n.eq.0) return

      if (trans.eq.'n') then
         leny = m
      else
         leny = n
      end if

         kx = 1

      if (beta.ne.1d0) then

         if (beta.eq.0d0) then

            do i = 1, leny
               y(i) = 0d0
            end do

         else

            do i = 1, leny
               y(i) = beta*y(i)
            end do

         end if

      end if

      jx = kx

      if (trans.eq.'n') then

            n4 = 4*(n/4)
            do j = 1, n4, 4
               temp1 = x(jx)
               temp2 = x(jx + 1)
               temp3 = x(jx + 2)
               temp4 = x(jx + 3)

               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     *             temp4.ne.0d0) then

                  do i = 1, m
                     y(i) = ((((y(i) + temp1*a(i,j))
     *                        + temp2*a(i,j + 1))
     *                        + temp3*a(i,j + 2))
     *                        + temp4*a(i,j + 3))
                  end do

               end if

               jx = jx + 4

            end do 

            do j = n4 + 1, n, 1

               temp = x(jx)

               if (temp.ne.0d0) then

                  do i = 1, m
                     y(i) = y(i) + temp*a(i,j)
                  end do 

               end if

               jx = jx + 1

            end do 

      else

            m4 = 4*(m/4)

            do j = 1, m4, 4
               temp1 = x(jx)
               temp2 = x(jx + 1)
               temp3 = x(jx + 2)
               temp4 = x(jx + 3)

               if (temp1.ne.0d0.or.temp2.ne.0d0.or.temp3.ne.0d0.or.
     $             temp4.ne.0d0) then

                  do i = 1, n
                     y(i) = ((((y(i) + temp1*a(j, i))
     $                        + temp2*a(j + 1, i))
     $                        + temp3*a(j + 2, i))
     $                        + temp4*a(j + 3, i))
                  end do

               end if

               jx = jx + 4

            end do
c                                 clean-up
            do j = m4 + 1, m, 1

               temp = x(jx)
               if (temp.ne.0d0) then

                  do i = 1, n
                     y(i) = y(i) + temp*a(j, i)
                  end do

               end if

               jx = jx + 1

            end do

      end if

      end

      subroutine dcopy (n, x, incx, y, incy)
c----------------------------------------------------------------------
      implicit none

      integer i, ix, iy, incx, incy, n

      double precision x(*), y(*)
c----------------------------------------------------------------------
      if (n.gt.0) then
         if (incx.eq.incy.and.incy.gt.0) then
            do iy = 1, 1 + (n - 1)*incy, incy
               y(iy) = x(iy)
            end do
         else
            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incy.gt.0) then
               do iy = 1, 1 + (n - 1)*incy, incy
                  y(iy) = x(ix)
                  ix      = ix      + incx
               end do 
            else
               iy = 1 - (n - 1)*incy
               do i = 1, n
                  y(iy) = x(ix)
                  iy      = iy      + incy
                  ix      = ix      + incx
               end do 
            end if
         end if
      end if

      end

      subroutine daxpy (n, alpha, x, incx, y)
c----------------------------------------------------------------------
      implicit none

      integer incx, n, i, ix, iy

      double precision x(*), y(*), alpha
c----------------------------------------------------------------------
      if (n.gt.0) then

         if (alpha.ne.0d0) then

            if (incx.eq.1.and.incx.gt.0) then

               do ix = 1, 1 + (n - 1)*incx, incx
                  y(ix) = alpha*x(ix) + y(ix)
               end do

            else
              
               iy = 1

               if (incx.gt.0) then

                  do ix = 1, 1 + (n - 1)*incx, incx
                     y(iy) = alpha*x(ix) + y(iy)
                     iy      = iy            + 1
                  end do

               else

                  ix = 1 - (n - 1)*incx

                  do i = 1, n
                     y(iy) = alpha*x(ix) + y(iy)
                     ix      = ix            + incx
                     iy      = iy            + 1
                  end do 

               end if
            end if
         end if
      end if

      end

      subroutine dtrsv (trans, n, a, lda, x)
c----------------------------------------------------------------------
      implicit none

      integer lda, n,  i, j

      character*1 trans

      double precision  a(lda,*), x(*), temp
c----------------------------------------------------------------------
      if (n.eq.0) return

      if (trans.eq.'n') then

               do j = n, 1, -1
                  if (x(j).ne.0d0) then
                     x(j) = x(j)/a(j,j)
                     temp = x(j)
                     do i = j - 1, 1, -1
                        x(i) = x(i) - temp*a(i,j)
                     end do 
                  end if
               end do 

      else
               do j = 1, n
                  temp = x(j)
                  do i = 1, j - 1
                     temp = temp - a(i,j)*x(i)
                  end do 
                  temp = temp/a(j,j)
                  x(j) = temp
               end do 

      end if

      end

      subroutine dger (m, n, x, y, a, lda)
c----------------------------------------------------------------------
      implicit none

      integer lda, m, n, i, j, jy

      double precision a(lda,*), x(*), y(*), alpha, temp
c----------------------------------------------------------------------
c     quick return if possible.

      if (m.eq.0 .or. n.eq.0 ) return

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.


      jy = 1
      alpha = -1d0

      do j = 1, n
            if (y(jy).ne.0d0) then
               temp = alpha*y(jy)
               do i = 1, m
                  a(i,j) = a(i,j) + x(i)*temp
               end do
            end if
            jy = jy + 1
      end do 

      end

      subroutine dswap (n, x, incx, y, incy)
c----------------------------------------------------------------------
      implicit none

      integer incx, incy, n, i, ix, iy

      double precision x(*), y(*), temp 
c----------------------------------------------------------------------

      if (n.gt.0) then

         if ((incx.eq.incy).and.(incy.gt.0)) then

            do iy = 1, 1 + (n - 1)*incy, incy
               temp    = x(iy)
               x(iy) = y(iy)
               y(iy) = temp
            end do

         else

            if (incx.ge.0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if

            if (incy.gt.0) then

               do iy = 1, 1 + (n - 1)*incy, incy
                  temp    = x(ix)
                  x(ix) = y(iy)
                  y(iy) = temp
                  ix      = ix      + incx
               end do 

            else

               iy = 1 - (n - 1)*incy

               do i = 1, n
                  temp    = x(ix)
                  x(ix) = y(iy)
                  y(iy) = temp
                  iy      = iy      + incy
                  ix      = ix      + incx
               end do 

            end if
         end if
      end if

      end

      integer function idamax (n, x)
c----------------------------------------------------------------------
      implicit none

      double precision  x(*), xmax

      integer i, imax, ix, n
c----------------------------------------------------------------------
      if (n.gt.0) then

         imax = 1

         if (n.gt.1) then

            xmax = dabs(x(1))
            ix   = 1

            do i = 2, n

               ix = ix + 1

               if (xmax.lt.dabs(x(ix))) then
                  xmax = dabs(x(ix))
                  imax = i
               end if

            end do

         end if

      else

         imax = 0

      end if

      idamax = imax

      end