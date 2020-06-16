

      subroutine nlpopt(n,nclin,lda,ldr,a,bl,bu,
     *                  objfun,iter,istate,clamda,objf,gradu,r,
     *                  x,iw,leniw,w,lenw,iuser,user,ifail,jprint)

c     n        the number of variables (dimension of  x)
c     nclin    the number of linear constraints (rows of the matrix  a)
c     lda      leading dimension of a
c     ldr      dimension of r
c     a        linear constraint matrix
c     bl       lower bounds on variables, linear constraints, non-linear constraints
c     bu       upper bounds, check bigbnd, -bigbnd
c     objfun   computes objective function
c     iter     output, number of major iterations
c     istate   output (cold start), constraint state (-2 violated lower, 
c              -1 violated upper, 0 inactive, 1 lower, 2 upper, 3 equality)
c              n+nclin
c     clamda   output (cold start), lagrangian multipliers, n+nclin
c     objf     output, value of objfun
c     gradu    output, grad objfun
c     r        output (cold start), upper triangular cholesky factor of r
c     x        input, estimate; output final solution
c     iw       array of dimension leniw = 3*n + nclin + 2*ncnln
c     leniw    as above
c     work     array of dimension lenw
c     lenw     lenw = 20*n ...
c              if nclin > 0, ...+ 2*n^2 + 11*nclin
c              if ncnln > 0, ...+ n*nclin + 2*n*ncnln + 21*ncnln
c     isuer    integer array to be passed to objfun, confun
c     user     real array to be passed to objfun, confun
c     ifail    on entry 0, -1, 1; -1 recommended (see chp p01);
c              on exit 0 is ok, ifail < 0 => mode set < 0; 
c              1-9 error conditions.

c     nlpopt   solves the nonlinear program
c
c            minimize                   f(x)
c
c                                    (   x )
c            subject to    bl  .le.  ( a*x )  .le.  bu
c                                    ( c(x))
c
c     where  f(x)  is a smooth scalar function,  a  is a constant matrix
c     and  c(x)  is a vector of smooth nonlinear functions.  the
c     feasible region is defined by a mixture of linear and nonlinear
c     equality or inequality constraints on  x.
c
c     the dimensions of the problem are...

c     nlpopt   uses a sequential quadratic programming algorithm, with a
c     positive-definite quasi-newton approximation to the transformed
c     hessian  q'hq  of the lagrangian function (which will be stored in
c     the array  r).

      implicit none 

      character*6       srname
      parameter (srname='nlpopt')
      integer mxparm
      parameter (mxparm=30)
      integer lenls
      parameter (lenls=20)
      integer lennp
      parameter (lennp=35)

      double precision zero, point3, point8
      parameter (zero=0.0d+0,point3=3.3d-1,point8=0.8d+0)
      double precision point9, one
      parameter (point9=0.9d+0,one=1.0d+0)
      double precision growth
      parameter (growth=1.0d+2)

      double precision objf
      integer ifail, iter, lda, ldr, leniw, lenw, n,
     *                  nclin, jprint

      double precision a(lda,*), bl(n+nclin), bu(n+nclin),
     *                  clamda(n+nclin),
     *                  gradu(n), r(ldr,*), user(*), w(lenw), x(n)
      integer istate(n+nclin), iuser(*), iw(leniw)

      external objfun

      double precision asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldq, ldt,
     *                  lfdset, lformh, lines1, lines2, lprob,
     *                  lverfy, lvlder, lvldif, lvrfyc, msgls, msgnp,
     *                  nactiv, ncdiff, ncolt, nfdiff, nfree, nlnf,
     *                  nlnj, nlnx, nload, nn, nnclin, nncnln, nout,
     *                  nprob, nsave, nz
      logical           cmdbg, incrun, lsdbg, npdbg, unitq

      double precision rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm), wmach
      integer icmdbg(5), ilsdbg(5), inpdbg(5),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), jverfy(4), locls(lenls),
     *                  locnp(lennp)

      double precision amin, cond, condmx, ctx, epsmch, errmax, fdchk,
     *                  fdnorm, feamax, feamin, obj, rootn, rteps, ssq1,
     *                  suminf, xnorm
      integer i, ianrmj, idbg, idbgsv, ikx, info, inform,
     *                  itmxsv, itns, j, jinf, jmax, lanorm, laqp, lax,
     *                  lcjac, lcjdx, lclam, lcmul, ldaqp, ldcj, ldfju,
     *                  ldx, lfeatl, lgq, lgrad, lhctrl, lhfrwd, liperm,
     *                  litotl, lkactv, lkx, lneedc, lq, lres, lres0,
     *                  lrho, lrlam, lt, lv, lvj, lwrk1, lwrk2, lwrk3,
     *                  lwtinf, lwtotl, m, maxact, maxnz, minact,
     *                  minfxd, mjrdbg, mnrdbg, msgqp, mxfree, nact1,
     *                  nartif, nctotl, nerr, nerror, nfun, ngq, ngrad,
     *                  nlperr, nmajor, nminor, nplin, nrank, nrejtd,
     *                  nres, nstate, numinf, nz1
      logical           cold, linobj, named, needfd, overfl, rowerr,
     *                  vertex
      character*11      title

      double precision rprmls(mxparm), rprmnp(mxparm)
      integer iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(2)

      double precision dnrm2, adivb, dlantr
      integer errmsg
      external          dnrm2, adivb, dlantr, errmsg

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common            /ae04uc/locnp
      common/ cstmch /wmach(9)
      common/ be04nb /ldt, ncolt, ldq
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /ce04uc/lvrfyc, jverfy
      common            /de04nb/asize, dtmax, dtmin
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /fe04nb/icmdbg, cmdbg
      common            /fe04nc/nactiv, nfree, nz, unitq
      common            /fe04uc/inpdbg, npdbg
      common            /ge04uc/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /he04uc/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
      data              title/' ** oink **'/

c     set the machine-dependent constants.

      epsmch = wmach(3)
      rteps = wmach(4)
      nout = 6
      nerr = 6
c
      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9
c
      rhomax = one/epsmch
      rootn = dsqrt(dble(n))
c
c     default names will be provided for variables during printing.
c
      named = .false.
      inform = 0
c                                 print flags
      msgnp = jprint
      msgqp = jprint - 10
c
c     set the default values for the parameters.
c
      call e04ucx(n,nclin,title)
c
      needfd = lvlder .eq. 0.or.lvlder .eq. 2
      cold = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1
c
      nplin = n + nclin
      nctotl = nplin
c
c     assign the dimensions of arrays in the parameter list of e04ucz.
c     economies of storage are possible if the minimum number of active
c     constraints and the minimum number of fixed variables are known in
c     advance.  the expert user should alter minact and minfxd
c     accordingly.
c
      minact = 0
      minfxd = 0
c
      mxfree = n - minfxd
      maxact = max(1,min(n,nclin))
      maxnz = n - (minfxd+minact)
c
      if (nclin.eq.0) then
         ldq = 1
         ldt = 1
         ncolt = 1
      else
         ldq = max(1,mxfree)
         ldt = max(maxnz,maxact)
         ncolt = mxfree
      end if

      m = 1
      ldfju = 2

      ldaqp = max(nclin,1)
      if (nclin.gt.0) ldaqp = lda

c     e04ucp  defines the arrays that contain the locations of various
c     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      call e04ucp(n,nclin,nctotl,litotl,lwtotl)

c     allocate certain addresses that are not allocated in e04ucp.

      lax = lwtotl + 1
      lwtotl = lax + nclin - 1
      lax = min(lax,lwtotl)

c     check input parameters and storage limits.

      call e04nbz(nerror,msgnp,lcrash,leniw,lenw,litotl,lwtotl,n,nclin,
     *            istate,iw,named,names,bigbnd,bl,bu,x,m,lda,ldr,
     *            ldfju,nerr,ifail)

      if (nerror.gt.0) then
         inform = 9
         go to 80
      end if

      lkactv = locls(1)
      lanorm = locls(2)
      lcjdx = locls(3)
      lres = locls(5)
      lres0 = locls(6)
      lgq = locls(9)
      lrlam = locls(10)
      lt = locls(11)
      lq = locls(12)
      lwtinf = locls(13)
      lwrk1 = locls(14)
c
      lkx = locnp(1)
      liperm = locnp(2)
      laqp = locnp(3)
      ldx = locnp(7)
      lfeatl = locnp(10)
      lwrk2 = locnp(12)
c
      lcmul = locnp(16)
      lrho = locnp(20)
      lwrk3 = locnp(21)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)
      lcjac = locnp(27)
      lgrad = locnp(28)
c
      ldcj = 1
c
      tolrnk = zero
      rcndbd = dsqrt(hcndbd)

c                                 load the arrays of feasibility tolerances.
      if (tolfea.gt.zero) call sload (nplin,tolfea,w(lfeatl),1)
c                                 forward finite difference interval, 
c                                 lfdset 0 - evaluate increment
c                                        1 - use specified increment
c                                        2 - use  default
c                                 cdint central difference.
      if (lfdset.eq.0) then
         fdchk = dsqrt(epsrf)
      else if (lfdset.eq.1) then
         fdchk = fdint
         cdint = fdint
      else
         fdchk = w(lhfrwd)
      end if
c                                  optimality tolerance is ftol
      nfun = 0
      ngrad = 0
      nstate = 1

c     if required,  compute the problem functions.
c     if the constraints are nonlinear,  the first call of confun
c     sets up any constant elements in the jacobian matrix.  a copy of
c     the jacobian (with constant elements set) is placed in  cjacu.

      if (lverfy.ge.10) then
         xnorm = dnrm2(n,x,1)
         lvrfyc = lverfy - 10

         call e04ucy(info,msgnp,nstate,lvlder,nfun,ngrad,n,
     *               objfun,iw(lneedc),bigbnd,epsrf,cdint,
     *               fdint,fdchk,fdnorm,objf,xnorm,bl,bu,
     *               w(ldx),w(lgrad),gradu,
     *               w(lhfrwd),w(lhctrl),x,w(lwrk1),w(lwrk2),w,lenw,
     *               iuser,user)

         if (info.ne.0) then
            if (info.gt.0) inform = 7
            if (info.lt.0) inform = info
            go to 80
         end if
         nstate = 0

      end if

      call icopy (5,ilsdbg,1,icmdbg,1)

      if (nclin.gt.0) then
         ianrmj = lanorm
         do 20 j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
   20    continue
         call scond (nclin,w(lanorm),1,asize,amin)
      end if

      call scond (nplin,w(lfeatl),1,feamax,feamin)
      call dcopy (nplin,w(lfeatl),1,w(lwtinf),1)
      call dscal (nplin,(one/feamin),w(lwtinf),1)


c     the input values of x and (optionally)  istate are used by
c     e04ncu  to define an initial working set.

      vertex = .false.
      call e04ncu(cold,vertex,nclin,nplin,nactiv,nartif,nfree,n,lda,
     *            istate,iw(lkactv),bigbnd,tolact,a,w(lax),bl,bu,x,
     *            w(lwrk1),w(lwrk2))

      nres = 0
      ngq = 0
      condmx = one/epspt5

      if (lcrash.le.1) then

c        cold or warm start. the upper-triangular matrix r is the factor
c        of an approximate lagrangian hessian.

         unitq = .true.
         iter = 0

         ikx = lkx
         do 40 i = 1, n
            iw(ikx) = i
            ikx = ikx + 1
   40    continue

         if (cold) then
            call smload('upper-triangular',n,n,zero,one,r,ldr)
            rfrobn = rootn

            nrank = 0

         else

c           r will be updated while finding a feasible x.

            nrank = nlnx
            call sload (nlnx,(zero),w(lres0),1)

         end if
c
         incrun = .true.
         rhonrm = zero
         rhodmp = one
         scale = one
c ncnln
         call sload (0,(zero),w(lrho),1)
c

c        re-order kx so that the free variables come first.
c        if a warm start is required, nrank will be nonzero and the
c        factor r will be updated.

         call e04ncx(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldq,lda,ldr,
     *               ldt,istate,iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)

      end if
c

c     factorize the linear constraints in the initial working set.

      if (nactiv.gt.0) then
         nact1 = nactiv
         nactiv = 0
c
         call e04ncy(unitq,vertex,inform,1,nact1,nactiv,nartif,nz,nfree,
     *               nrank,nrejtd,nres,ngq,n,ldq,lda,ldr,ldt,istate,
     *               iw(lkactv),iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)
      end if
c
      if (lcrash.le.1) then

c        cold or warm start.  move  x  on to the linear constraints and
c        find a feasible point.

         ssq1 = zero
         linobj = .false.
         call e04nch (linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,n,
     *               nplin,ldq,lda,ldr,ldt,istate,iw(lkactv),iw(lkx),
     *               jmax,errmax,ctx,xnorm,a,w(lax),bl,bu,w(lgq),w(lres)
     *               ,w(lres0),w(lfeatl),r,w(lt),x,w(lq),w(lwrk1),
     *               w(lwrk2))

         jinf = 0
         lclam = lwrk2

         idbgsv = idbg
         if (idbg.gt.0) then
            idbg = nminor + 1
         end if

         itmxsv = itmax1
         itmax1 = nminor

         call e04ncz ('fp problem',named,names,linobj,unitq,nlperr,itns,
     *               jinf,nclin,nplin,nactiv,nfree,nrank,nz,nz1,n,lda,
     *               ldr,istate,iw(lkactv),iw(lkx),ctx,obj,ssq1,suminf,
     *               numinf,xnorm,bl,bu,a,w(lclam),w(lax),w(lfeatl),r,x,
     *               w)

         itmax1 = itmxsv

         if (nlperr.gt.0) then
            inform = 2
            go to 80
         else if (msgqp.gt.0) then
            write (rec,fmt=99987)
            call x04bay(iprint,2,rec)
         end if

         idbg = idbgsv
         call icopy(5,inpdbg,1,icmdbg,1)

      end if

      if (lcrash.gt.0) then

c        check for a bad r.

         rfrobn = dlantr('frobenius norm','upper','non-unit diagonal',n,
     *            n,r,ldr,w)
         call scond (n,r,ldr+1,drmax,drmin)
         cond = adivb (drmax,drmin,overfl)

         if (cond.gt.rcndbd.or.rfrobn.gt.rootn*growth*drmax) then
c                                 refactorize the hessian and bound the condition estimator.
            if (msgnp.gt.0) then
               write (rec,fmt=99986)
               call x04baf(iprint,rec(1))
            end if
            call e04udr(unitq,n,nfree,nz,ldq,ldr,iw(liperm),iw(lkx),
     *                  w(lgq),r,w(lq),w(lwrk1),w(lres0))
         end if
      end if
c                                 check the gradients at a feasible x.
      lvrfyc = lverfy
      if (lverfy.ge.10) lvrfyc = -1

      call e04ucy (info,msgnp,nstate,lvlder,nfun,ngrad,n,
     *            objfun,iw(lneedc),bigbnd,epsrf,cdint,
     *            fdint,fdchk,fdnorm,objf,xnorm,bl,bu,
     *            w(ldx),w(lgrad),gradu,
     *            w(lhfrwd),w(lhctrl),x,w(lwrk1),w(lwrk2),w,lenw,iuser,
     *            user)

      if (info.ne.0) then
         if (info.gt.0) inform = 7
         if (info.lt.0) inform = info
         go to 80
      end if

      call dcopy (n,w(lgrad),1,w(lgq),1)
      call cmqmul (6,n,nz,nfree,ldq,unitq,iw(lkx),w(lgq),w(lq),w(lwrk1))
c                                 solve the problem.
      iuser(2) = 1

         call e04ucz (named,names,unitq,inform,iter,n,nclin,nctotl
     *               ,nactiv,nfree,nz,ldaqp,ldr,nfun,ngrad,
     *               istate,iw(lkactv),iw(lkx),objf,fdnorm,xnorm,
     *               objfun,a,w(lax),bl,bu,clamda,
     *               w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw,iuser,user)

c                                if required, form the triangular factor of the hessian.
c                                first,  form the square matrix  r  such that  h = r'r.
c                                compute the  qr  factorization of  r.

      if (lformh.gt.0) then

         write (*,*) 'hessian'
         call errpau

         lv = lwrk2
         do 60 j = 1, n
            if (j.gt.1) call sload (j-1,zero,w(lv),1)
c
            lvj = lv + j - 1
            call dcopy (n-j+1,r(j,j),ldr,w(lvj),1)
            call cmqmul (3,n,nz,nfree,ldq,unitq,iw(lkx),w(lv),w(lq),
     *                  w(lwrk1))
            call dcopy (n,w(lv),1,r(j,1),ldr)
   60    continue

         call sgeqr (n,n,r,ldr,w(lwrk1),info)
      end if

c     print messages if required.

   80 if (msgnp.gt.0) then
         if (inform.lt.0) write (rec,fmt=99999)
         if (inform.eq.0) write (rec,fmt=99998)
         if (inform.eq.1) write (rec,fmt=99997)
         if (inform.eq.2) write (rec,fmt=99996)
         if (inform.eq.3) write (rec,fmt=99995)
         if (inform.eq.4) write (rec,fmt=99994)
         if (inform.eq.6) write (rec,fmt=99993)
         if (inform.eq.7) write (rec,fmt=99992)
         if (inform.eq.9) write (rec,fmt=99991) nerror
         call x04bay(iprint,2,rec)

         if (inform.ge.0.and.inform.lt.7) then
            if (nlperr.eq.0) then
               write (rec,fmt=99990) objf
               call x04bay(iprint,2,rec)
            else
               if (nlperr.eq.3) then
                  write (rec,fmt=99989) suminf
                  call x04bay(iprint,2,rec)
               else
                  write (rec,fmt=99988) suminf
                  call x04bay(iprint,2,rec)
               end if
            end if
         end if
      end if

c     recover the optional parameters set by the user.

      call icopy(mxparm,ipsvls,1,iprmls,1)
      call dcopy(mxparm,rpsvls,1,rprmls,1)
      call icopy(mxparm,ipsvnp,1,iprmnp,1)
      call dcopy(mxparm,rpsvnp,1,rprmnp,1)

      if (inform.lt.9) call dcopy(n,w(lgrad),1,gradu,1)

      if (msgnp.gt.0) then
         if (inform.ne.0.and.(ifail.eq.0.or.ifail.eq.-1)) then
         if (inform.lt.0) write (rec,fmt=99985)
         if (inform.eq.1) write (rec,fmt=99984)
         if (inform.eq.2) write (rec,fmt=99983)
         if (inform.eq.3) write (rec,fmt=99982)
         if (inform.eq.4) write (rec,fmt=99981)
         if (inform.eq.6) write (rec,fmt=99980)
         if (inform.eq.7) write (rec,fmt=99979)
         if (inform.eq.9) write (rec,fmt=99978) nerror
         call x04bay(nerr,2,rec)
         end if
c                                 on input:
c                                 ifail = 1 soft, silent
c                                 ifail = 0 hard, noisy
c                                 ifail = -1, soft noisy
c                                 errmsg sets ifail = inform
         ifail = errmsg (ifail,inform,srname,' ',0,rec)

      end if

c     end of  nlpopt. (npsol)

99999 format (/' exit nlpopt - user requested termination.')
99998 format (/' exit nlpopt - optimal solution found.')
99997 format (/' exit nlpopt - optimal solution found, but requested a',
     *       'ccuracy not achieved.')
99996 format (/' exit nlpopt - no feasible point for the linear constr',
     *       'aints.')
99995 format (/' exit nlpopt - no feasible point for the nonlinear con',
     *       'straints.')
99994 format (/' exit nlpopt - too many major iterations.             ')
99993 format (/' exit nlpopt - current point cannot be improved upon. ')
99992 format (/' exit nlpopt - large errors found in the derivatives. ')
99991 format (/' exit nlpopt - ',i7,' errors found in the input parame',
     *       'ters.  problem abandoned.')
99990 format (/' final objective value =',g16.7)
99989 format (/' minimum sum of infeasibilities =',g16.7)
99988 format (/' final sum of infeasibilities =',g16.7)
99987 format (/' the linear constraints are feasible.')
99986 format (' xxx  bad initial hessian,   r  refactorized.')
99985 format (/' ** user requested termination.')
99984 format (/' ** optimal solution found, but requested accuracy not',
     *       ' achieved.')
99983 format (/' ** no feasible point for the linear constraints.')
99982 format (/' ** no feasible point for the nonlinear constraints.')
99981 format (/' ** too many major iterations.             ')
99980 format (/' ** current point cannot be improved upon. ')
99979 format (/' ** large errors found in the derivatives. ')
99978 format (/' ** ',i7,' errors found in the input parameters.  prob',
     *       'lem abandoned.')
      end


      subroutine e04ucx(n,nclin,title)

c     e04ucx  loads the default values of parameters not set in the
c     options file.

      implicit none

      integer mxparm
      parameter (mxparm=30)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision point3, point8
      parameter (point3=3.3d-1,point8=0.8d+0)
      double precision point9, two
      parameter (point9=0.9d+0,two=2.0d+0)
      double precision tenp6, hundrd
      parameter (tenp6=1.0d+6,hundrd=10.0d+1)
      double precision rdummy
      integer idummy
      parameter (rdummy=-11111.d+0,idummy=-11111)
      double precision gigant
      parameter (gigant=1.0d+20*0.99999d+0)
      double precision wrktol
      parameter (wrktol=1.0d-2)
      integer n, nclin
      character*(*)     title
      double precision bigbnd, bigdx, bndlow, bndupp, cdint, ctol,
     *                  dxlim, epspt3, epspt5, epspt8, epspt9, epsrf,
     *                  eta, fdint, ftol, hcndbd, tolact, tolfea, tolrnk
      integer idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, lfdset,
     *                  lformh, lines1, lines2, lprob, lverfy, lvlder,
     *                  lvldif, lvrfyc, msgls, msgnp, ncdiff, nfdiff,
     *                  nlnf, nlnj, nlnx, nload, nn, nnclin, nncnln,
     *                  nout, nprob, nsave
      logical           cmdbg, lsdbg, newopt, npdbg

      double precision rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm), wmach
      integer icmdbg(5), ilsdbg(5), inpdbg(5),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), jverfy(4)

      double precision condbd, dctol, epsmch
      integer i, idbg, j, k, lent, mjrdbg, mnrdbg, msg1, msg2,
     *                  msgqp, nctotl, nmajor, nminor, nplin
      character*16      key

      double precision rprmls(mxparm), rprmnp(mxparm)
      integer iprmls(mxparm), iprmnp(mxparm)
      character*3       chess(0:1)
      character*4       icrsh(0:2)
      character*80      rec(4)
c     .. intrinsic functions ..
      intrinsic         dble, len, max, mod

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common/ cstmch /wmach(9)
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /ce04uc/lvrfyc, jverfy
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ee04uc/newopt
      common            /fe04nb/icmdbg, cmdbg
      common            /fe04uc/inpdbg, npdbg
      common            /ge04uc/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /he04uc/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
      data              icrsh(0), icrsh(1), icrsh(2)/'cold', 'warm',
     *                  'hot '/
      data              chess(0), chess(1)/' no', 'yes'/

      epsmch = wmach(3)
      nout = 6
c
      condbd = max(one/(hundrd*epsmch*dble(n)),tenp6)
c
      nplin = n + nclin
      nctotl = nplin

c     make a dummy call e04ucq to ensure that the defaults are set.

      call e04ucq(nout,'*',key)
      newopt = .true.
c
c     save the optional parameters set by the user.  the values in
c     iprmls, rprmls, iprmnp and rprmnp may be changed to their
c     default values.
c
      call icopy(mxparm,iprmls,1,ipsvls,1)
      call dcopy (mxparm,rprmls,1,rpsvls,1)
      call icopy(mxparm,iprmnp,1,ipsvnp,1)
      call dcopy (mxparm,rprmnp,1,rpsvnp,1)
c
      if (msgnp.eq.idummy) msgnp = 10
      if (msgqp.eq.idummy) msgqp = 0
      if (iprnt.lt.0) iprnt = nout
      if (isumry.lt.0.or.(msgnp.lt.5.and.msgqp.lt.5)) isumry = -1
      iprint = iprnt
      isumm = isumry
      if (lcrash.lt.0.or.lcrash.gt.2) lcrash = 0
      if (lvlder.lt.0.or.lvlder.gt.3) lvlder = 3
      if (lformh.lt.0.or.lformh.gt.1) lformh = 0
c
      if (nmajor.lt.0) nmajor = max(50,3*nplin)
      if (nminor.lt.1) nminor = max(50,3*nctotl)
      if (mjrdbg.lt.0) mjrdbg = 0
      if (mnrdbg.lt.0) mnrdbg = 0
      if (idbg.lt.0.or.idbg.gt.nmajor) idbg = 0
      if (mjrdbg.eq.0.and.mnrdbg.eq.0) idbg = nmajor + 1
      nlnf = n
      nlnj = n
      nlnx = n
      if (jvrfy2.le.0.or.jvrfy2.gt.n) jvrfy2 = n
      if (jvrfy1.le.0.or.jvrfy1.gt.jvrfy2) jvrfy1 = 1
      if (jvrfy4.le.0.or.jvrfy4.gt.n) jvrfy4 = n
      if (jvrfy3.le.0.or.jvrfy3.gt.jvrfy4) jvrfy3 = 1
      if ((lverfy.lt.-1.or.lverfy.gt.13)
     *   .or.(lverfy.ge.4.and.lverfy.le.9)) lverfy = 0
c
      if (ksave.le.0) ksave = nmajor + 1
      if (nsave.lt.0) nsave = 0
      if (nsave.eq.0) ksave = nmajor + 1
      if (nload.lt.0) nload = 0
      if (lcrash.le.1) nload = 0
      if (nload.eq.0.and.lcrash.eq.2) lcrash = 0
c
      if (tolact.lt.zero.or.tolact.ge.one) tolact = wrktol
      if (tolfea.lt.epsmch.or.tolfea.ge.one) tolfea = epspt5
      if (epsrf.lt.epsmch.or.epsrf.ge.one) epsrf = epspt9
      lfdset = 0
      if (fdint.lt.zero) lfdset = 2
      if (fdint.eq.rdummy) lfdset = 0
      if (fdint.ge.epsmch.and.fdint.lt.one) lfdset = 1
      if (lfdset.eq.1.and.(cdint.lt.epsmch.or.cdint.ge.one))
     *    cdint = epsrf**point3
      if (bigbnd.le.zero) bigbnd = gigant
      if (bigdx.le.zero) bigdx = max(gigant,bigbnd)
      if (dxlim.le.zero) dxlim = two
      if (eta.lt.zero.or.eta.ge.one) eta = point9
c                                  optimality tolerance
      if (ftol.lt.epsrf.or.ftol.ge.one) ftol = epsrf**point8
c
      if (hcndbd.lt.one) hcndbd = condbd
c
      dctol = epspt5
      if (lvlder.lt.2) dctol = epspt3
      if (ctol.lt.epsmch.or.ctol.ge.one) ctol = dctol
c
      itmax1 = max(50,3*(n+nclin))
      jverfy(1) = jvrfy1
      jverfy(2) = jvrfy2
      jverfy(3) = jvrfy3
      jverfy(4) = jvrfy4
c
      npdbg = idbg .eq. 0
      cmdbg = npdbg
      lsdbg = npdbg
c
      k = 1
      msg1 = mjrdbg
      msg2 = mnrdbg
      do 20 i = 1, 5
         inpdbg(i) = mod(msg1/k,10)
         icmdbg(i) = inpdbg(i)
         ilsdbg(i) = mod(msg2/k,10)
         k = k*10
   20 continue
c
      if (msgnp.gt.0) then
c
c        print the title. if no hot start is specified, the parameters
c        are final and can be printed.
c
         lent = len(title)
         write (rec,fmt=99988) (title(j:j),j=1,lent)
         call x04bay(iprint,2,rec)

         if (lcrash.le.1) then
            write (rec,fmt=99999)
            call x04bay(iprint,3,rec)
            write (rec,fmt=99998) nclin, tolfea, n, tolact
            call x04bay(iprint,3,rec)
            write (rec,fmt=99997) bigbnd, icrsh(lcrash), bigdx, epsmch,
     *        dxlim, chess(lformh)
            call x04bay(iprint,4,rec)
            write (rec,fmt=99996) 0, ctol, nlnf, ftol, nlnj, eta
            call x04bay(iprint,4,rec)
            write (rec,fmt=99995) lvlder, epsrf, lverfy, isumry
            call x04bay(iprint,3,rec)
            if (lverfy.gt.0) then
               write (rec,fmt=99994) jvrfy1, jvrfy2
               call x04bay(iprint,2,rec)

            end if
            write (rec,fmt=99992) nmajor, msgnp, nminor, msgqp
            call x04bay(iprint,3,rec)
c
            if (lvlder.lt.3) then
               if (lfdset.eq.0) then
                  write (rec,fmt=99991)
                  call x04bay(iprint,2,rec)
               else if (lfdset.eq.1) then
                  write (rec,fmt=99990) fdint, cdint
                  call x04bay(iprint,2,rec)
               else if (lfdset.eq.2) then
                  write (rec,fmt=99989)
                  call x04bay(iprint,2,rec)
               end if
            end if
c
         end if
      end if
c
      if (msgnp.ge.5.or.msgqp.ge.5) then
         if (isumm.ge.0.and.isumm.ne.iprint) then
            lent = len(title)
            write (rec,fmt=99987) (title(j:j),j=1,lent)
            call x04bay(isumm,2,rec)
         end if
      end if

c     end of  e04ucx. (npdflt)

99999 format (/' parameters',/' ----------')
99998 format (/' linear constraints.....',i10,7x,'linear feasibility..',
     *       '...',1p,d10.2,/' variables..............',i10,7x,'crash ',
     *       'tolerance........',1p,d10.2)
99997 format (/' infinite bound size....',1p,d10.2,7x,a4,' start......',
     *       '.......',/' infinite step size.....',1p,d10.2,7x,'eps (m',
     *       'achine precision)',1p,d10.2,/' step limit.............',
     *       1p,d10.2,7x,'hessian................',7x,a3)
99996 format (/' nonlinear constraints..',i10,7x,'nonlinear feasibilit',
     *       'y..',1p,d10.2,/' nonlinear objectiv vars',i10,7x,'optima',
     *       'lity tolerance...',1p,d10.2,/' nonlinear jacobian vars',
     *       i10,7x,'linesearch tolerance...',1p,d10.2)
99995 format (/' derivative level.......',i10,7x,'function precision..',
     *       '...',1p,d10.2,/' verify level...........',i10,7x,'monito',
     *       'ring file........',i10)
99994 format (/' start obj chk at varble',i10,7x,'stop obj chk at varb',
     *       'le.',i10)
99993 format (' start con chk at varble',i10,7x,'stop con chk at varbl',
     *       'e.',i10)
99992 format (/' major iterations limit.',i10,7x,'major print level...',
     *       '...',i10,/' minor iterations limit.',i10,7x,'minor print',
     *       ' level......',i10)
99991 format (/' difference intervals to be computed.')
99990 format (/' difference interval....',1p,d10.2,7x,'central diffce ',
     *       'interval',1p,d10.2)
99989 format (/' user-supplied difference intervals.')
99988 format (/80a1)
99987 format (/11a1,' monitoring information ')
      end

      subroutine e04ucp(n,nclin,nctotl,litotl,lwtotl)

c     e04ucp   allocates the addresses of the work arrays for e04ucz and
c     e04ncz.

      implicit none

      integer lenls
      parameter (lenls=20)
      integer lennp
      parameter (lennp=35)

      integer litotl, lwtotl, n, nclin, nctotl

      integer ldt, ldzy, ncolt
      logical           npdbg

      integer inpdbg(5), locls(lenls), locnp(lennp)

      integer ladx, lanorm, laqp, lbl, lbu, lc1mul, lcjac,
     *                  lcjdx, lcmul, lcs1, lcs2, ldlam, ldslk, ldx,
     *                  lenaqp, lent, lenzy, lfeatl, lgq, lgq1, lgrad,
     *                  lhctrl, lhfrwd, liperm, lkactv, lkx, lneedc,
     *                  lqpadx, lqpdx, lqpgq, lqphz, lqptol, lrho,
     *                  lrlam, lrpq, lrpq0, lslk, lslk1, lt, lwrk1,
     *                  lwrk2, lwrk3, lwtinf, lx1, lzy, miniw, minw

      common            /ae04nc/locls
      common            /ae04uc/locnp
      common/ be04nb /ldt, ncolt, ldzy
      common            /fe04uc/inpdbg, npdbg

      miniw = litotl + 1
      minw = lwtotl + 1

c     assign array lengths that depend upon the problem dimensions.

      if (nclin.eq.0) then
         lent = 0
         lenzy = 0
      else
         lent = ldt*ncolt
         lenzy = ldzy*ldzy
      end if

      lenaqp = 0

c
      lkactv = miniw
      lkx = lkactv + n
      lneedc = lkx + n
      liperm = lneedc
      miniw = liperm + nctotl
c
      lhfrwd = minw
      lhctrl = lhfrwd + n
      lanorm = lhctrl + n
      lqpgq = lanorm + nclin 
      lgq = lqpgq + n
      lrlam = lgq + n
      lt = lrlam + n
      lzy = lt + lent
      minw = lzy + lenzy
c
      locls(1) = lkactv
      locls(2) = lanorm
      locls(8) = lqpgq
      locls(9) = lgq
      locls(10) = lrlam
      locls(11) = lt
      locls(12) = lzy

c     assign the addresses for the workspace arrays used by  e04ucu.

      lqpadx = minw
      lqpdx = lqpadx + nclin 
      lrpq = lqpdx + n
      lrpq0 = lrpq + n
      lqphz = lrpq0 + n
      lwtinf = lqphz + n
      lwrk1 = lwtinf + nctotl
      lqptol = lwrk1 + nctotl
      minw = lqptol + nctotl

      locls(3) = lqpadx
      locls(4) = lqpdx
      locls(5) = lrpq
      locls(6) = lrpq0
      locls(7) = lqphz
      locls(13) = lwtinf
      locls(14) = lwrk1
      locls(15) = lqptol

c     assign the addresses for arrays used in e04ucz.

      laqp = minw
      ladx = laqp + lenaqp
      lbl = ladx + nclin 
      lbu = lbl + nctotl
      ldx = lbu + nctotl
      lgq1 = ldx + n
      lfeatl = lgq1 + n
      lx1 = lfeatl + nctotl
      lwrk2 = lx1 + n
      minw = lwrk2 + nctotl
c
      locnp(1) = lkx
      locnp(2) = liperm
      locnp(3) = laqp
      locnp(4) = ladx
      locnp(5) = lbl
      locnp(6) = lbu
      locnp(7) = ldx
      locnp(8) = lgq1
      locnp(10) = lfeatl
      locnp(11) = lx1
      locnp(12) = lwrk2
c
      lcs1 = minw
      lcs2 = lcs1
      lc1mul = lcs2 
      lcmul = lc1mul 
      lcjdx = lcmul
      ldlam = lcjdx 
      ldslk = ldlam 
      lrho = ldslk 
      lwrk3 = lrho 
      lslk1 = lwrk3 
      lslk = lslk1 
      minw = lslk 
c
      locnp(13) = lcs1
      locnp(14) = lcs2
      locnp(15) = lc1mul
      locnp(16) = lcmul
      locnp(17) = lcjdx
      locnp(18) = ldlam
      locnp(19) = ldslk
      locnp(20) = lrho
      locnp(21) = lwrk3
      locnp(22) = lslk1
      locnp(23) = lslk
      locnp(24) = lneedc
c
      lcjac = minw
      lgrad = lcjac 
      minw = lgrad + n
c
      locnp(25) = lhfrwd
      locnp(26) = lhctrl
      locnp(27) = lcjac
      locnp(28) = lgrad
c
      litotl = miniw - 1
      lwtotl = minw - 1

c     end of  e04ucp. (nploc)

      end

      subroutine e04nbz(nerror,msglvl,lcrash,liwork,lwork,litotl,lwtotl,
     *                  n,nclin,istate,kx,named,names,bigbnd,bl,
     *                  bu,x,m,lda,ldr,ldfj,nerr,ifail)

c     e04nbz   checks the input data for nlpopt and e04upf.
      implicit none

      double precision bigbnd
      integer ifail, lcrash, lda, ldfj, ldr, litotl,
     *                  liwork, lwork, lwtotl, m, msglvl, n, nclin,
     *                  nerr, nerror
      logical           named
      double precision bl(n+nclin), bu(n+nclin), x(n)
      integer istate(n+nclin), kx(n)
      character*8       names(*)
      integer iprint, isumm, lines1, lines2, nout
      double precision b1, b2
      integer is, j, k, l
      logical           ok
      character*5       id(3)
      character*80      rec(4)

c     .. intrinsic functions ..
      intrinsic         abs, max

      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      data              id(1), id(2), id(3)/'varbl', 'l con', 'n con'/

      nerror = 0

c     check  m.

      if (m.le.0) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99997) m
            call x04bay(nerr,3,rec)
         end if
      end if


c     check  n.

      if (n.le.0) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99996) n
            call x04bay(nerr,3,rec)
         end if
      end if


c     check  nclin and ncnln.



c     check  lda.

      if (lda.lt.max(1,nclin)) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99993) lda, nclin
            call x04bay(nerr,3,rec)
         end if
      end if



c     check  ldfj.

      if (ldfj.lt.m) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99991) ldfj, m
            call x04bay(nerr,3,rec)
         end if
      end if


c     check  ldr.

      if (ldr.lt.n) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99990) ldr, n
            call x04bay(nerr,3,rec)
         end if
      end if


c     check if there is enough workspace to solve the problem.

      ok = litotl .le. liwork.and.lwtotl .le. lwork
      if (.not. ok) then
         nerror = nerror + 1
         if (ifail.eq.0.or.ifail.eq.-1) then
            write (rec,fmt=99999) liwork, lwork, litotl, lwtotl
            call x04bay(nerr,3,rec)
            write (rec,fmt=99998)
            call x04bay(nerr,2,rec)
         end if
      else if (msglvl.gt.0) then
         write (rec,fmt=99999) liwork, lwork, litotl, lwtotl
         call x04bay(iprint,3,rec)
      end if

      if (nerror.eq.0) then


c        check the bounds on all variables and constraints.

         do 20 j = 1, n + nclin
            b1 = bl(j)
            b2 = bu(j)
            ok = b1 .lt. b2.or.(b1.eq.b2.and.abs(b1).lt.bigbnd)
            if (.not. ok) then
               nerror = nerror + 1
               if (j.gt.n+nclin) then
                  k = j - n - nclin
                  l = 3
               else if (j.gt.n) then
                  k = j - n
                  l = 2
               else
                  k = j
                  l = 1
               end if
               if (named) then
                  if (b1.eq.b2) then
                     if (ifail.eq.0.or.ifail.eq.-1) then
                        write (rec,fmt=99989) names(j), j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0.or.ifail.eq.-1) then
                        write (rec,fmt=99988) names(j), j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               else
                  if (b1.eq.b2) then
                     if (ifail.eq.0.or.ifail.eq.-1) then
                        write (rec,fmt=99987) id(l), k, j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0.or.ifail.eq.-1) then
                        write (rec,fmt=99986) id(l), k, j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               end if
            end if
   20    continue


c        if warm start, check  istate.

         if (lcrash.eq.1) then
            do 40 j = 1, n + nclin
               is = istate(j)
               ok = is .ge. (-2).and.is .le. 4
               if (.not. ok) then
                  nerror = nerror + 1
                  if (ifail.eq.0.or.ifail.eq.-1) then
                     write (rec,fmt=99985) j, j, is
                     call x04bay(nerr,3,rec)
                  end if
               end if
   40       continue
         end if
      end if

c     end of  e04nbz. (cmchk)

99999 format (/' workspace provided is     iwork(',i8,'),  work(',i8,
     *       ').',/' to solve problem we need  iwork(',i8,'),  work(',
     *       i8,').')
99998 format (/' ** not enough workspace to solve problem.')
99997 format (/' ** on entry, m.le.0:',/'    m = ',i16)
99996 format (/' ** on entry, n.le.0:',/'    n = ',i16)
99995 format (/' ** on entry, nclin.lt.0:',/'    nclin = ',i16)
99994 format (/' ** on entry, ncnln.lt.0:',/'    ncnln = ',i16)
99993 format (/' ** on entry, lda.lt.max(1,nclin):',/'    lda = ',i16,
     *       '   nclin = ',i16)
99992 format (/' ** on entry, ldcj.lt.max(1,ncnln):',/'    ldcj = ',i16,
     *       '   ncnln = ',i16)
99991 format (/' ** on entry, ldfj.lt.m:',/'    ldfj = ',i16,'   m = ',
     *       i16)
99990 format (/' ** on entry, ldr.lt.n:',/'    ldr = ',i16,'   n = ',
     *       i16)
99989 format (/' ** on entry, the equal bounds on  ',a8,'  are infinit',
     *       'e (because',/'    bl(',i4,').eq.beta and bu(',i4,').eq.b',
     *       'eta, but abs(beta).ge.bigbnd):',/'    beta =',g16.7,' bi',
     *       'gbnd =',g16.7)
99988 format (/' ** on entry, the bounds on  ',a8,'  are inconsistent:',
     *       /'    bl(',i4,') =',g16.7,'   bu(',i4,') =',g16.7)
99987 format (/' ** on entry, the equal bounds on  ',a5,i3,'  are infi',
     *       'nite (because',/'    bl(',i4,').eq.beta and bu(',i4,').e',
     *       'q.beta, but abs(beta).ge.bigbnd):',/'    beta =',g16.7,
     *       ' bigbnd =',g16.7)
99986 format (/' ** on entry, the bounds on  ',a5,i3,'  are inconsiste',
     *       'nt:',/'    bl(',i4,') =',g16.7,'   bu(',i4,') =',g16.7)
99985 format (/' ** on entry with a warm start, istate(',i4,') is out ',
     *       'of range:',/'    istate(',i4,') = ',i16)
      end

      subroutine e04ucq(nout,buffer,key)

c     e04ucq   decodes the option contained in  buffer  in order to set
c     a parameter value in the relevant element of the parameter arrays.

c     input:

c     nout   a unit number for printing error messages.
c            nout  must be a valid unit.

c     output:

c     key    the first keyword contained in buffer.

c     e04ucq  calls e04udx and the subprograms
c                 lookup, scannr, tokens, upcase
c     (now called e04udy, e04udw, e04udv, e04udu)

      integer mxparm
      parameter (mxparm=30)
      integer maxkey, maxtie, maxtok
      parameter (maxkey=41,maxtie=19,maxtok=10)
      integer idummy
      double precision rdummy
      logical           sorted
      double precision zero
      parameter (idummy=-11111,rdummy=-11111.0d+0,sorted=.true.,
     *                  zero=0.0d+0)

      integer nout
      character*16      key
      character*(*)     buffer

      double precision bigbnd, bigdx, bndlow, bndupp, cdint, ctol,
     *                  dxlim, epsrf, eta, fdint, ftol, hcndbd, tolact,
     *                  tolfea, tolrnk
      integer idbgls, idbgnp, iprnt, isumry, itmax1, itmax2,
     *                  itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, ksave,
     *                  lcrash, ldbgls, ldbgnp, lformh, lprob, lverfy,
     *                  lvlder, msgls, msgnp, nlnf, nlnj, nlnx, nload,
     *                  nn, nnclin, nncnln, nprob, nsave

      double precision rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm)

      double precision rvalue
      integer i, idbg, ivalue, lenbuf, loc1, loc2, mjrdbg,
     *                  mnrdbg, msgqp, nmajor, nminor, ntoken
      logical           first, more, number
      character*16      key2, key3, value
      character*80      rec

      double precision rprmls(mxparm), rprmnp(mxparm)
      integer iprmls(mxparm), iprmnp(mxparm)
      character*16      keys(maxkey), ties(maxtie), token(maxtok)

      logical           e04udx
      external          e04udx

c     .. intrinsic functions ..
      intrinsic         index, len

      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ge04uc/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /he04uc/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
     *
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/, first

      data              first/.true./
      data keys/'begin           ', 'central         ',
     *  'cold            ', 'condition       ', 'constraints     ',
     *  'crash           ', 'debug           ', 'defaults        ',
     *  'derivative      ', 'difference      ', 'end             ',
     *  'feasibility     ', 'function        ', 'hessian         ',
     *  'hot             ', 'infinite        ', 'iprmls          ',
     *  'iterations      ', 'iters:iterations', 'itns :iterations',
     *  'linear          ', 'linesearch      ', 'list            ',
     *  'lower           ', 'major           ', 'minor           ',
     *  'monitoring      ', 'nolist          ', 'nonlinear       ',
     *  'optimality      ', 'print           ', 'problem         ',
     *  'row             ', 'rprmls          ', 'start           ',
     *  'step            ', 'stop            ', 'upper           ',
     *  'variables       ', 'verify          ', 'warm            '/
      data              ties/'bound           ', 'constraints     ',
     *                  'debug           ', 'feasibility     ',
     *                  'gradients       ', 'iterations      ',
     *                  'iters:iterations', 'itns :iterations',
     *                  'jacobian        ', 'level           ',
     *                  'no              ', 'no.      :number',
     *                  'number          ', 'objective       ',
     *                  'print           ', 'step            ',
     *                  'tolerance       ', 'variables       ',
     *                  'yes             '/


      if (first) then
         first = .false.
         do 20 i = 1, mxparm
            rprmls(i) = rdummy
            iprmls(i) = idummy
            rprmnp(i) = rdummy
            iprmnp(i) = idummy
   20    continue
      end if

c     eliminate comments and empty lines.
c     a '*' appearing anywhere in buffer terminates the string.

      i = index(buffer,'*')
      if (i.eq.0) then
         lenbuf = len(buffer)
      else
         lenbuf = i - 1
      end if
      if (lenbuf.le.0) then
         key = '*'
         go to 80
      end if


c     extract up to maxtok tokens from the record.
c     ntoken returns how many were actually found.
c     key, key2, key3 are the first tokens if any, otherwise blank.

      call e04udv(buffer(1:lenbuf),maxtok,ntoken,token)
      key = token(1)
      key2 = token(2)
      key3 = token(3)

c     certain keywords require no action.

      if (key.eq.' '.or.key.eq.'begin') go to 80
      if (key.eq.'list'.or.key.eq.'nolist') go to 80
      if (key.eq.'end') go to 80

c     most keywords will have an associated integer or real value,
c     so look for it no matter what the keyword.

      i = 1
      number = .false.

   40 if (i.lt.ntoken.and..not. number) then
         i = i + 1
         value = token(i)
         number = e04udx(value)
         go to 40
      end if
c
      if (number) then
         read (value,fmt='(bn, e16.0)') rvalue
      else
         rvalue = zero
      end if
c
c     convert the keywords to their most fundamental form
c     (upper case, no abbreviations).
c     sorted says whether the dictionaries are in alphabetic order.
c     loci   says where the keywords are in the dictionaries.
c     loci = 0 signals that the keyword wasn't there.
c     loci < 0 signals that the keyword is ambiguous.
c
      call e04udy(maxkey,keys,sorted,key,loc1)
      if (loc1.lt.0) then
         write (rec,fmt=99996) key
         call x04baf(nout,rec)
         return
      end if
      call e04udy(maxtie,ties,sorted,key2,loc2)
c

c     decide what to do about each keyword.
c     the second keyword (if any) might be needed to break ties.
c     some seemingly redundant testing of more is used
c     to avoid compiler limits on the number of consecutive else ifs.

      more = .true.
      if (more) then
         more = .false.
         if (key.eq.'central     ') then
            cdint = rvalue
         else if (key.eq.'cold        ') then
            lcrash = 0
         else if (key.eq.'condition   ') then
            hcndbd = rvalue
         else if (key.eq.'constraints ') then
            nnclin = rvalue
         else if (key.eq.'crash       ') then
            tolact = rvalue
         else if (key.eq.'debug       ') then
            idbg = rvalue
         else if (key.eq.'defaults    ') then
            do 60 i = 1, mxparm
               iprmls(i) = idummy
               rprmls(i) = rdummy
               iprmnp(i) = idummy
               rprmnp(i) = rdummy
   60       continue
         else if (key.eq.'derivative  ') then
            lvlder = rvalue
         else if (key.eq.'difference  ') then
            fdint = rvalue
         else if (key.eq.'feasibility ') then
            tolfea = rvalue
            ctol = rvalue
         else if (key.eq.'function    ') then
            epsrf = rvalue
         else
            more = .true.
         end if
      end if
c
      if (more) then
         more = .false.
         if (key.eq.'hessian     ') then
            lformh = 0
            if (key2.eq.'yes         ') lformh = 1
         else if (key.eq.'hot         ') then
            lcrash = 2
         else if (key.eq.'infinite    ') then
            if (key2.eq.'bound       ') bigbnd = rvalue*0.99999d+0
            if (key2.eq.'step        ') bigdx = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'iprmls      ') then
c           allow things like  iprmls 21 = 100  to set iprmls(21) = 100
            ivalue = rvalue
            if (ivalue.ge.1.and.ivalue.le.mxparm) then
               read (key3,fmt='(bn, i16)') iprmls(ivalue)
            else
               write (rec,fmt=99997) ivalue
               call x04baf(nout,rec)
            end if
         else if (key.eq.'iterations  ') then
            nmajor = rvalue
         else if (key.eq.'linear      ') then
            if (key2.eq.'constraints ') nnclin = rvalue
            if (key2.eq.'feasibility ') tolfea = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'linesearch  ') then
            eta = rvalue
         else if (key.eq.'lower       ') then
            bndlow = rvalue
         else
            more = .true.
         end if
      end if
c
      if (more) then
         more = .false.
         if (key.eq.'major       ') then
            if (key2.eq.'debug       ') mjrdbg = rvalue
            if (key2.eq.'iterations  ') nmajor = rvalue
            if (key2.eq.'print       ') msgnp = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'minor       ') then
            if (key2.eq.'debug       ') mnrdbg = rvalue
            if (key2.eq.'iterations  ') nminor = rvalue
            if (key2.eq.'print       ') msgqp = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'monitoring  ') then
            isumry = rvalue
         else if (key.eq.'optimality  ') then
            ftol = rvalue
         else
            more = .true.
         end if
      end if
c
      if (more) then
         more = .false.
         if (key.eq.'print       ') then
            msgnp = rvalue
         else if (key.eq.'problem     ') then
            if (key2.eq.'number      ') nprob = rvalue
         else if (key.eq.'row         ') then
            if (key2.eq.'tolerance   ') ctol = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'rprmls      ') then
c           allow things like  rprmls 21 = 2  to set rprmls(21) = 2.0
            ivalue = rvalue
            if (ivalue.ge.1.and.ivalue.le.mxparm) then
               read (key3,fmt='(bn, e16.0)') rprmls(ivalue)
            else
               write (rec,fmt=99997) ivalue
               call x04baf(nout,rec)
            end if
         else if (key.eq.'start       ') then
            if (key2.eq.'constraints ') jvrfy3 = rvalue
            if (key2.eq.'objective   ') jvrfy1 = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'step        ') then
            dxlim = rvalue
         else if (key.eq.'stop        ') then
            if (key2.eq.'constraints ') jvrfy4 = rvalue
            if (key2.eq.'objective   ') jvrfy2 = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'upper       ') then
            bndupp = rvalue
         else if (key.eq.'variables   ') then
            nn = rvalue
         else if (key.eq.'verify      ') then
            if (key2.eq.'objective   ') lverfy = 1
            if (key2.eq.'constraints ') lverfy = 2
            if (key2.eq.'no          ') lverfy = -1
            if (key2.eq.'yes         ') lverfy = 3
            if (key2.eq.'gradients   ') lverfy = 3
            if (key2.eq.'level       ') lverfy = rvalue
            if (loc2.eq.0) lverfy = 3
         else if (key.eq.'warm        ') then
            lcrash = 1
         else
            write (rec,fmt=99999) key
            call x04baf(nout,rec)
         end if
      end if
c
   80 return
c
c
c     end of  e04ucq. (npkey)
c
99999 format (' xxx  keyword not recognized:         ',a)
99998 format (' xxx  second keyword not recognized:  ',a)
99997 format (' xxx  the parm subscript is out of range:',i10)
99996 format (' xxx  ambiguous keyword:              ',a)
      end

      subroutine e04udv(string,numin,numout,list)

c     description and usage:
c
c       an aid to parsing input data.  the individual 'tokens' in a
c     character string are isolated, converted to uppercase, and stored
c     in an array.  here, a token is a group of significant, contiguous
c     characters.  the following are non-significant, and hence may
c     serve as separators:  blanks, horizontal tabs, commas, colons,
c     and equal signs.  see e04udw for details.  processing continues
c     until the requested number of tokens have been found or the end
c     of the input string is reached.
c
c
c     parameters:
c
c     name    dimension  type  i/o/s  description
c     string              c    i      input string to be analyzed.
c     numin               i    i/o    number of tokens requested (input)
c     numout                          and found (output).
c     (numin and numout were both called number in the original)
c
c     list    numin       c      o    array of tokens, changed to upper
c                                    case.
c
c
c     external references:
c
c     name    description
c     e04udw  finds positions of first and last significant characters.
c     e04udu  converts a string to uppercase.
c
c
c     environment:  digital vax-11/780 vms fortran (fortran 77).
c               appears to satisfy the ansi fortran 77 standard.
c
c
c     notes:
c
c     (1)  implicit none is non-standard.  (has been commented out.)
c-----------------------------------------------------------------------

      character         blank
      parameter (blank=' ')

      integer numin, numout
      character*(*)     string

      character*(*)     list(numin)

      integer count, first, i, last, mark

c     .. intrinsic functions ..
      intrinsic         len

      first = 1
      last = len(string)
c
      count = 0
   20 continue
c
c        get delimiting indices of next token, if any.
c
      call e04udw(string,first,last,mark)
      if (last.gt.0) then
         count = count + 1
c
c           pass token to output string array, then change case.
c
         list(count) = string(first:mark)
         call e04udu(list(count))
         first = mark + 2
         if (count.lt.numin) go to 20
c
      end if
c
c
c     fill the rest of list with blanks and set number for output.
c
      do 40 i = count + 1, numin
         list(i) = blank
   40 continue
c
      numout = count

c     end of  e04udv. (optokn)
      end

      logical function e04udx(string)

c***********************************************************************
c     description and usage:
c
c        a simple(-minded) test for numeric data is implemented by
c        searching an input string for legitimate characters:
c                digits 0 to 9, d, e, -, + and .
c        insurance is provided by requiring that a numeric string
c        have at least one digit, at most one d, e or .
c        and at most two -s or +s.  note that a few ambiguities remain:
c
c           (a)  a string might have the form of numeric data but be
c                intended as text.  no general test can hope to detect
c                such cases.
c
c           (b)  there is no check for correctness of the data format.
c                for example a meaningless string such as 'e1.+2-'
c                will be accepted as numeric.
c
c        despite these weaknesses, the method should work in the
c        majority of cases.
c
c
c     parameters:
c
c        name    dimension  type  i/o/s  description
c        e04udx              l      o    set .true. if string appears
c                                        to be numerical data.
c        string              c    i      input data to be tested.
c
c
c     environment:  ansi fortran 77.
c
c
c     notes:
c
c        (1)  it is assumed that string is a token extracted by
c             e04udv, which will have converted any lower-case
c             characters to upper-case.
c
c        (2)  e04udv pads string with blanks, so that a genuine
c             number is of the form  '1234        '.
c             hence, the scan of string stops at the first blank.
c
c        (3)  complex data with parentheses will not look numeric.

      character*(*)           string

      integer j, length, ndigit, nexp, nminus, nplus,
     *                        npoint
      logical                 number
      character*1 atom
c     .. intrinsic functions ..
      intrinsic               len, lge, lle

      ndigit = 0
      nexp = 0
      nminus = 0
      nplus = 0
      npoint = 0
      number = .true.
      length = len(string)
      j = 0
c
   20 j = j + 1
      atom = string(j:j)
      if (lge(atom,'0').and.lle(atom,'9')) then
         ndigit = ndigit + 1
      else if (atom.eq.'d'.or.atom.eq.'e'.or.
     *         atom.eq.'D'.or.atom.eq.'E') then
         nexp = nexp + 1
      else if (atom.eq.'-') then
         nminus = nminus + 1
      else if (atom.eq.'+') then
         nplus = nplus + 1
      else if (atom.eq.'.') then
         npoint = npoint + 1
      else if (atom.eq.' ') then
         j = length
      else
         number = .false.
      end if
c
      if (number.and.j.lt.length) go to 20
c
      e04udx = number.and.ndigit .ge. 1.and.nexp .le. 1 .and.
     *         nminus .le. 2.and.nplus .le. 2.and.npoint .le. 1

c     end of  e04udx. (opnumb)
      end

      subroutine e04udy(ndict,dictry,alpha,key,entry)

c     description and usage:
c
c       performs dictionary lookups.  a pointer is returned if a
c     match is found between the input key and the corresponding
c     initial characters of one of the elements of the dictionary.
c     if a 'synonym' has been provided for an entry, the search is
c     continued until a match to a primary dictionary entry is found.
c     cases of no match, or multiple matches, are also provided for.
c
c     dictionary entries must be left-justified, and may be alphabetized
c     for faster searches.  secondary entries, if any, are composed of
c     two words separated by one or more characters such as blank, tab,
c     comma, colon, or equal sign which are treated as non-significant
c     by e04udw.  the first entry of each such pair serves as a synonym
c     for the second, more fundamental keyword.
c
c       the ordered search stops after the section of the dictionary
c     having the same first letters as the key has been checked, or
c     after a specified number of entries have been examined.  a special
c     dictionary entry, the currency symbol '$', will also terminate the
c     search.  this will speed things up if an appropriate dictionary
c     length parameter cannot be determined.  both types of search are
c     sequential.  see 'notes' below for some suggestions if efficiency
c     is an issue.
c
c
c     parameters:
c
c     name    dimension  type  i/o/s  description
c     ndict               i    i     number of dictionary entries to be
c                                    examined.
c     dictry  ndict       c    i     array of dictionary entries,
c                                    left-justified in their fields.
c                                    may be alphabetized for efficiency,
c                                    in which case alpha should be
c                                    .true.  entries with synonyms are
c                                    of the form
c                                    'entry : synonym', where 'synonym'
c                                    is a more fundamental entry in the
c                                    same dictionary.  note: don't build
c                                    'circular' dictionaries.
c     alpha               l    i     indicates whether the dictionary
c                                    is in alphabetical order, in which
c                                    case the search can be terminated
c                                    sooner.
c     key                 c    i/o   string to be compared against the
c                                    dictionary.  abbreviations are ok
c                                    if they correspond to a unique
c                                    entry in the dictionary.  key is
c                                    replaced on termination by its most
c                                    fundamental equivalent dictionary
c                                    entry (uppercase, left-justified)
c                                    if a match was found.
c     entry               i      o   dictionary pointer.  if .gt. 0, it
c                                    indicates which entry matched key.
c                                    in case of trouble, a negative
c                                    value means that a unique match
c                                    was not found - the absolute value
c                                    of entry points to the second
c                                    dictionary entry that matched key.
c                                    zero means that no match could be
c                                    found.  entry always refers to the
c                                    last search performed -
c                                    in searching a chain of synonyms,
c                                    a non-positive value will be
c                                    returned if there is any break,
c                                    even if the original input key
c                                    was found.
c
c
c     external references:
c
c     name    description
c     e04udw  finds first and last significant characters.

c     notes:
c
c     (1)  implicit none is non-standard.  (has been commented out.)
c
c     (2)  we have assumed that the dictionary is not too big.  if
c         many searches are to be done or if the dictionary has more
c         than a dozen or so entries, it may be advantageous to build
c         an index array of pointers to the beginning of the section
c         of the dictionary containing each letter, then pass in the
c         portion of the dictionary beginning with dictry (index).
c         (this won't generally work for dictionaries with synonyms.)
c         for very large problems, a completely different approach may
c         be advisable, e.g. a binary search for ordered dictionaries.
c
c     (3)  e04udy is case sensitive.  in most applications it will be
c         necessary to use an uppercase dictionary, and to convert the
c         input key to uppercase before calling e04udy.  companion
c         routines e04udv and pairs, available from the author, already
c         take care of this.
c
c     (4)  the key need not be left-justified.  any leading (or
c         trailing) characters which are 'non-significant' to e04udw
c         will be ignored.  these include blanks, horizontal tabs,
c         commas, colons, and equal signs.  see e04udw for details.
c
c     (5)  the ascii collating sequence for character data is assumed.
c         (n.b. this means the numerals precede the alphabet, unlike
c         common practice.)  this should not cause trouble on ebcdic
c         machines if dictry just contains alphabetic keywords.
c         otherwise it may be necessary to use the fortran lexical
c         library routines to force use of the ascii sequence.
c
c     (6)  parameter numsig sets a limit on the length of significant
c         dictionary entries.  special applications may require that
c         this be increased.  (it is 16 in the present version.)
c
c     (7)  no protection against 'circular' dictionaries is provided:
c         don't claim that a is b, and that b is a.  all synonym chains
c         must terminate.  other potential errors not checked for
c         include duplicate or mis-ordered entries.
c
c     (8)  the handling of ambiguities introduces some ambiguity:
c
c            alpha = .true.  a potential problem, when one entry
c                            looks like an abbreviation for another
c                            (eg. does 'a' match 'a' or 'ab') was
c                            resolved by dropping out of the search
c                            immediately when an 'exact' match is found.
c
c            alpha = .false. the programmer must ensure that the above
c                            situation does not arise: each dictionary
c                            entry must be recognizable, at least when
c                            specified to full length.  otherwise, the
c                            result of a search will depend on the
c                            order of entries.
c-----------------------------------------------------------------------
      character         blank, curly
      integer numsig
      parameter (blank=' ',curly='$',numsig=16)

      integer entry, ndict
      logical           alpha
      character*(*)     key

      character*(*)     dictry(ndict)

      integer first, i, ifrst, ifrst1, ilast, ilast1, ilen,
     *                  ilen1, ilst, imark, imark1, last, length, mark
      character*16      flag, target, trgt, trgt1
c     .. intrinsic functions ..
      intrinsic         min, len, lle


      entry = 0

c     isolate the significant portion of the input key (if any).
c
      first = 1
      last = min(len(key),numsig)
      call e04udw(key,first,last,mark)
c
      if (mark.gt.0) then
         target = key(first:mark)
c
c        look up target in the dictionary.
c
   20    continue
         length = mark - first + 1
c
c           select search strategy by cunning choice of termination test
c           flag.  the left curly bracket follows all the alphabetic
c           characters in the ascii collating sequence, but precedes the
c           vertical bar.
c
         if (alpha) then
            flag = target
         else
            flag = curly
         end if
c
c
c           perform search.
c           ---------------
c
         i = 0
   40    continue
         i = i + 1
         if (target(1:length).eq.dictry(i)(1:length)) then
            if (entry.eq.0) then
c
c                    first 'hit' - must still guard against ambiguities
c                    by searching until we've gone beyond the key
c                    (ordered dictionary) or until the end-of-dictionary
c                    mark is reached (exhaustive search).
c
               entry = i
c
c                    special handling if match is exact - terminate
c                    search.  we thus avoid confusion if one dictionary
c                    entry looks like an abbreviation of another.
c                    this fix won't generally work for un-ordered
c                    dictionaries.
c
               first = 1
               last = numsig
               call e04udw(dictry(entry),first,last,mark)
               if (mark.eq.length) i = ndict
            else
c                    if two hits check if they are attempting to
c                    indicate the same dictionary entry.
c
c                    extract keyword from first match found
c
               ilst = numsig
               ifrst = mark + 2
               call e04udw(dictry(entry),ifrst,ilst,imark)
               if (imark.gt.0) then
                  trgt = dictry(entry) (ifrst:imark)
                  ilen = imark - ifrst + 1
               else
                  trgt = dictry(entry) (first:mark)
                  ilen = mark - first + 1
               end if
c
c                    extract keyword from next match found
c
               ifrst = 1
               ilast = numsig
               call e04udw(dictry(i),ifrst,ilast,imark)
               ilast1 = numsig
               ifrst1 = imark + 2
               call e04udw(dictry(i),ifrst1,ilast1,imark1)
               if (imark1.gt.0) then
                  trgt1 = dictry(i) (ifrst1:imark1)
                  ilen1 = imark1 - ifrst1 + 1
               else
                  trgt1 = dictry(i) (ifrst:imark)
                  ilen1 = imark - ifrst + 1
               end if
c
c                    if keywords not identical then ambiguity
c
               if (trgt(1:ilen).ne.trgt1(1:ilen1)) then
c
c                       oops - two hits.  abnormal termination.
c                       ---------------------------------------
c
                  entry = -i
                  return
               end if
            end if
         end if
c
c           check whether we've gone past the appropriate section of the
c           dictionary.  the test on the index provides insurance and an
c           optional means for limiting the extent of the search.
c
         if (lle(dictry(i)(1:length),flag).and.i.lt.ndict) go to 40
c
c
c           check for a synonym.
c           --------------------
c
         if (entry.gt.0) then
c
c              look for a second entry 'behind' the first entry.  first
c              and mark were determined above when the hit was detected.
c
            first = mark + 2
            call e04udw(dictry(entry),first,last,mark)
            if (mark.gt.0) then
c
c                 re-set target and dictionary pointer, then repeat the
c                 search for the synonym instead of the original key.
c
               target = dictry(entry) (first:mark)
               entry = 0
               go to 20
c
            end if
         end if
c
      end if
      if (entry.gt.0) key = dictry(entry)

c     end of e04udy.  (cmlook/oplook)
      end

      subroutine e04ucy(inform,msgnp,nstate,lvlder,nfun,ngrad,
     *                  n,objfun,needc,bigbnd,epsrf,
     *                  cdint,fdint,fdchk,fdnorm,objf,xnorm,bl,bu,
     *                  dx,grad,gradu,hforwd,hcntrl,x,
     *                  wrk1,wrk2,w,lenw,iuser,user)


c     e04ucy  performs the following...
c     (1)  computes the objective and constraint values objf and c.
c     (2)  evaluates the user-provided gradients in cjacu and gradu.
c     (3)  counts the missing gradients.
c     (4)  loads the known gradients into grad and cjac.
c     (5)  checks that the known gradients are programmed correctly.
c     (6)  computes the missing gradient elements.

      double precision rdummy
      parameter (rdummy=-11111.0d+0)

      double precision bigbnd, cdint, epsrf, fdchk, fdint, fdnorm,
     *                  objf, xnorm
      integer inform, lenw, lvlder, msgnp, n,
     *                 nfun, ngrad, nstate

      double precision bl(n), bu(n), 
     *                  dx(n), grad(n),
     *                  gradu(n), hcntrl(*), hforwd(*), user(*),
     *                  w(lenw), wrk1(n), wrk2(n), x(n)
      integer iuser(*), needc(*)

      external          objfun

      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  lvrfyc, ncdiff, nfdiff, nout

      integer jverfy(4)

      integer i, infog, infoj, j, mode, ncset
      logical           centrl, needfd

      character*80      rec(3)


      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy

c
      infog = 0
      infoj = 0
      nfdiff = 0
      ncdiff = 0

c     repeat the procedure above for the objective function.

      needfd = lvlder .eq. 0.or.lvlder .eq. 2
c
      if (needfd) call sload (n,rdummy,gradu,1)
c                                 output the initial value
      iuser(2) = 1

      call objfun(mode,n,x,objf,gradu,nstate,iuser,user)
c                                 shut off output for derivative
c                                 evaluation
      iuser(2) = 0 

      call dcopy(n,gradu,1,grad,1)

      if (needfd) then
c
c        count the number of missing gradient elements.
c
         do 60 j = 1, n
            if (gradu(j).eq.rdummy) nfdiff = nfdiff + 1
   60    continue
c
         if (nstate.eq.1) then
            if (nfdiff.eq.0) then
               if (lvlder.eq.0) lvlder = 1
               if (lvlder.eq.2) lvlder = 3
               if (msgnp.gt.0) then
                  write (rec,fmt=99997) lvlder
                  call x04bay(iprint,3,rec)
               end if
            else
               if (msgnp.gt.0) then
                  write (rec,fmt=99996) n - nfdiff, n, nfdiff
                  call x04bay(iprint,3,rec)
               end if
            end if
         end if
      end if
c
      nfun = nfun + 1
      ngrad = ngrad + 1

c     check whatever gradient elements have been provided.

      if (lvrfyc.ge.0) then

         if (needfd) then
            call e04xax(mode,msgnp,n,bigbnd,epsrf,epspt3,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,wrk1,iuser,
     *                  user)

            infog = mode
         end if
      end if

      if (needfd) then

c        compute the missing gradient elements.

         call fdinc (mode,msgnp,lvlder,n,bigbnd,epsrf,
     *               fdnorm,objf,objfun,needc,bl,bu,
     *               grad,gradu,hforwd,hcntrl,x,dx,iuser,user)

         if (lfdset.gt.0) then
            centrl = lvldif .eq. 2
            call evalfd (centrl,mode,n,bigbnd,cdint,
     *                  fdint,fdnorm,objf,objfun,needc,bl,bu,
     *                  grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)

         end if
      end if

      inform = infoj + infog

c     end of  e04ucy. (npchkd)

99999 format (/' all jacobian elements have been set.',/' derivative l',
     *       'evel increased to ',i4)
99998 format (/' the user sets ',i6,'   out of',i6,'   jacobian elemen',
     *       'ts.',/' each iteration, ',i6,'   jacobian elements will ',
     *       'be estimated numerically.')
99997 format (/' all objective gradient elements have been set.',/' de',
     *       'rivative level increased to ',i4)
99996 format (/' the user sets ',i6,'   out of',i6,'   objective gradi',
     *       'ent elements.',/' each iteration, ',i6,'   gradient elem',
     *       'ents will be estimated numerically.')
      end


      subroutine e04udw(string,first,last,mark)

c     description and usage:
c
c       looks for non-blank fields ('tokens') in a string, where the
c     fields are of arbitrary length, separated by blanks, tabs, commas,
c     colons, or equal signs.  the position of the end of the 1st token
c     is also returned, so this routine may be conveniently used within
c     a loop to process an entire line of text.
c
c       the procedure examines a substring, string (first : last), which
c     may of course be the entire string (in which case just call e04udw
c     with first .le. 1 and last .ge. len (string)).  the indices
c     returned are relative to string itself, not the substring.
c
c
c     parameters:
c
c     name    dimension  type  i/o/s  description
c     string              c    i      text string containing data to be
c                                    scanned.
c     first               i    i/o    index of beginning of substring.
c                                    if .le. 1, the search begins with
c                                    1.
c                                    output is index of beginning of
c                                    first non-blank field, or 0 if no
c                                    token was found.
c     last                i    i/o    index of end of substring.
c                                    if .ge. len (string), the search
c                                    begins with len (string).  output
c                                    is index of end of last non-blank
c                                    field, or 0 if no token was found.
c     mark                i      o    points to end of first non-blank
c                                    field in the specified substring.
c                                    set to 0 if no token was found.
c
c
c     environment:  digital vax-11/780 vms fortran (fortran 77).
c               ansi fortran 77, except for the tab character ht.
c
c     notes:
c
c     (1)  implicit none is non-standard.  constant ht (tab) is defined
c         in a non-standard way:  the char function is not permitted
c         in a parameter declaration (ok on vax, though).  for absoft
c         fortran 77 on 68000 machines, use ht = 9.  in other cases, it
c         may be best to declare ht as a variable and assign
c         ht = char(9) on ascii machines, or char(5) for ebcdic.
c
c     (2)  the pseudo-recursive structure was chosen for fun.  it is
c         equivalent to three do loops with embedded go tos in sequence.
c
c     (3)  the variety of separators recognized limits the usefulness of
c         this routine somewhat.  the intent is to facilitate handling
c         such tokens as keywords or numerical values.  in other
c         applications, it may be necessary for all printing characters
c         to be significant.  a simple modification to statement
c         function solid will do the trick.
c
c
c     author:  robert kennelly, informatics general corporation.
c
c
c     development history:
c
c     29 dec. 1984    rak    initial design and coding, (very) loosely
c                           based on scan_string by ralph carmichael.
c     25 feb. 1984    rak    added ':' and '=' to list of separators.
c     16 apr. 1985    rak    defined solid in terms of variable dummy
c                           (previous re-use of string was ambiguous).
c
c-----------------------------------------------------------------------

      character         blank, equal, colon, comma, rparn, lparn
      parameter (blank=' ',equal='=',colon=':',comma=',',
     *                  rparn=')',lparn='(')

      integer first, last, mark
      character*(*)     string

      integer begin, end, length
      character         dummy
c     .. intrinsic functions ..
      intrinsic         max, min, len
c     .. statement functions ..
      logical           solid
c     .. statement function definitions ..
      solid(dummy) = (dummy.ne.blank).and.(dummy.ne.colon)
     *              .and.(dummy.ne.comma).and.(dummy.ne.equal)
     *              .and.(dummy.ne.rparn).and.(dummy.ne.lparn)

      mark = 0
      length = len(string)
      begin = max(first,1)
      end = min(length,last)
c
c     find the first significant character ...
c
      do 60 first = begin, end, +1
         if (solid(string(first:first))) then
c
c           ... then the end of the first token ...
c
            do 40 mark = first, end - 1, +1
               if (.not. solid(string(mark+1:mark+1))) then
c
c                 ... and finally the last significant character.
c
                  do 20 last = end, mark, -1
                     if (solid(string(last:last))) then
                        return
                     end if
   20             continue
c
c                 everything past the first token was a separator.
c
                  last = last + 1
                  return
               end if
   40       continue
c
c           there was nothing past the first token.
c
            last = mark
            return
         end if
   60 continue

c     whoops - the entire substring string (begin : end) was composed of
c     separators .

      first = 0
      mark = 0
      last = 0

c     end of  e04udw. (opscan)
      end

      subroutine e04udu(string)

c     purpose:  this subroutine changes all lower case letters in the
c               character string to upper case.
c
c     method:   each character in string is treated in turn.  the
c               intrinsic function index effectively allows a table
c               lookup, with the local strings low and upp acting as
c               two tables. this method avoids the use of char and
c               ichar, which appear be different on ascii and ebcdic
c               machines.
c
c     arguments
c     arg       dim     type i/o/s description
c     string       *       c   i/o   character string possibly
c                                   containing some lower-case
c                                   letters  on input; strictly
c                                   upper-case letters on output
c                                   with no change to any
c                                   non-alphabetic characters.
c
c     external references:
c     len    - returns the declared length of a character variable.
c     index  - returns the position of second string within first.
c-----------------------------------------------------------------------
      character*(*)     string

      integer i, j
      character*1 c
      character*26      low, upp
c     .. intrinsic functions ..
      intrinsic         index, len, lge, lle

      data              low/'abcdefghijklmnopqrstuvwxyz'/,
     *                  upp/'abcdefghijklmnopqrstuvwxyz'/

c
      do 20 j = 1, len(string)
         c = string(j:j)
         if (lge(c,'a').and.lle(c,'z')) then
c           if (c.ge.'a'.and.c.le.'z') then
            i = index(low,c)
            if (i.gt.0) string(j:j) = upp(i:i)
         end if
   20 continue

c     end of  e04udu. (opuppr)

      end


      subroutine e04xax(inform,msglvl,n,bigbnd,epsrf,oktol,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,y,iuser,user)

c     e04xax  checks if the gradients of the objective function have
c     been coded correctly.
c
c     on input,  the value of the objective function at the point x is
c     stored in objf.  the corresponding gradient is stored in gradu.
c     if any gradient element has not been specified,  it will have a
c     dummy value.  missing values are not checked.
c
c     a cheap test is first undertaken by calculating the directional
c     derivative using two different methods. if this proves
c     satisfactory and no further information is desired, e04xax is
c     terminated. otherwise, the routine fdinc1 is called to give
c     optimal step-sizes and a forward-difference approximation to
c     each element of the gradient for which a test is deemed
c     necessary, either by the program or the user.
c
c     other inputs:
c
c        x         the n-dimensional point at which the
c                  gradient is to be verified.
c        epsrf     the positive bound on the relative error
c                  associated with computing the function at
c                  the point x.
c        oktol     the desired relative accuracy which the
c                  elements of the gradient should satisfy.
c
c     lvrfyc has the following meaning...
c
c     -1        do not perform any check.
c     0        do the cheap test only.
c     1 or 3   do both cheap and full test.


      double precision rdummy
      parameter (rdummy=-11111.0d+0)
      double precision zero, half, point9
      parameter (zero=0.0d+0,half=0.5d+0,point9=0.9d+0)
      double precision one, two, ten
      parameter (one=1.0d+0,two=2.0d+0,ten=1.0d+1)
      character*4       lbad, lgood
      parameter (lbad='bad?',lgood='  ok')

      double precision bigbnd, epsrf, fdchk, objf, oktol, xnorm
      integer inform, msglvl, n

      double precision bl(n), bu(n), dx(n), grad(n), gradu(n), user(*),
     *                  x(n), y(n)
      integer iuser(*)

      external          objfun

      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lines1, lines2, lvrfyc, nout
      logical           npdbg

      integer inpdbg(5), jverfy(4)

      double precision biglow, bigupp, cdest, dxj, dxmult, emax, epsa,
     *                  errbnd, error, f1, f2, fdest, gdiff, gdx, gj,
     *                  gsize, h, hopt, hphi, objf1, sdest, stepbl,
     *                  stepbu, xj
      integer info, iter, itmax, j, j1, j2, jmax, mode,
     *                  ncheck, ngood, nstate, nwrong
      logical           const, debug, done, first, headng, needed, ok
      character*4       key

      character*18      result(0:4)
      character*120     rec(4)

      double precision ddot1
      external          ddot1

c     .. intrinsic functions ..
      intrinsic         abs, max, min, sqrt

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy
      common            /fe04uc/inpdbg, npdbg

      data              result/'                 ', 'constant?      ',
     *                  'linear or odd?   ', 'too nonlinear?',
     *                  'small derivative?'/

c
      inform = 0
      needed = lvrfyc .eq. 0.or.lvrfyc .eq. 1.or.lvrfyc .eq. 3
      if (.not. needed) return
c
      if (msglvl.gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,4,rec)
      end if
      debug = npdbg.and.inpdbg(5) .gt. 0
      nstate = 0
c
      biglow = -bigbnd
      bigupp = bigbnd
      iuser(2) = 0

c     perform the cheap test.

      h = (one+xnorm)*fdchk
c
      if (n.le.100) then
         dxmult = 0.9
      else if (n.le.250) then
         dxmult = 0.99
      else
         dxmult = 0.999
      end if
c
      dxj = one/n
      do 20 j = 1, n
         dx(j) = dxj
         dxj = -dxj*dxmult
   20 continue
c

c     do not perturb x(j) if the  j-th  element is missing.
c     compute the directional derivative.

      ncheck = 0
      do 40 j = 1, n
         if (grad(j).eq.rdummy) then
            dx(j) = zero
         else
            ncheck = ncheck + 1
c
            xj = x(j)
            stepbl = -one
            stepbu = one
            if (bl(j).gt.biglow) stepbl = max(stepbl,bl(j)-xj)
            if (bu(j).lt.bigupp.and.bu(j).gt.bl(j))
     *          stepbu = min(stepbu,bu(j)-xj)
c
            if (half*(stepbl+stepbu).lt.zero) then
               dx(j) = dx(j)*stepbl
            else
               dx(j) = dx(j)*stepbu
            end if
         end if
   40 continue
c
      if (ncheck.eq.0) then
         if (msglvl.gt.0) then
            write (rec,fmt=99989)
            call x04bay(iprint,2,rec)
         end if
         return
      end if
      gdx = ddot1 (n,gradu,1,dx)
c

c     make forward-difference approximation along  p.

      call dcopy (n,x,1,y,1)
      call daxpy (n,h,dx,1,y,1)
c
      mode = 0
      call objfun(mode,n,y,objf1,gradu,nstate,iuser,user)
c
      gdiff = (objf1-objf)/h
      error = abs(gdiff-gdx)/(abs(gdx)+one)
c
      ok = error .le. oktol
c
      if (msglvl.gt.0) then
         if (ok) then
            write (rec,fmt=99998)
            call x04bay(iprint,2,rec)
         else
            write (rec,fmt=99997)
            call x04bay(iprint,2,rec)
         end if
         write (rec,fmt=99996) gdx, gdiff
         call x04bay(iprint,3,rec)
      end if
c
      if (error.ge.point9) inform = 1
c

c     element-wise check.

      if (lvrfyc.eq.1.or.lvrfyc.eq.3) then
         headng = .true.
         itmax = 3
         ncheck = 0
         nwrong = 0
         ngood = 0
         jmax = 0
         emax = zero
         j1 = jverfy(1)
         j2 = jverfy(2)
c

c        loop over each of the elements of  x.

         do 80 j = j1, j2
c
            if (grad(j).ne.rdummy) then
c              ---------------------------------------------------------
c              check this gradient element.
c              ---------------------------------------------------------
               ncheck = ncheck + 1
               gj = grad(j)
               gsize = abs(gj)
               xj = x(j)
c              ---------------------------------------------------------
c              find a finite-difference interval by iteration.
c              ---------------------------------------------------------
               iter = 0
               epsa = epsrf*(one+dabs(objf))
               cdest = zero
               sdest = zero
               first = .true.
c
               stepbl = biglow
               stepbu = bigupp
               if (bl(j).gt.biglow) stepbl = bl(j) - xj
               if (bu(j).lt.bigupp) stepbu = bu(j) - xj
c
               hopt = two*(one+dabs(xj))*dsqrt(epsrf)
               h = ten*hopt
               if (half*(stepbl+stepbu).lt.zero) h = -h
c
c              +             repeat
   60          x(j) = xj + h
               call objfun(mode,n,x,f1,gradu,nstate,iuser,user)

               x(j) = xj + h + h
               call objfun(mode,n,x,f2,gradu,nstate,iuser,user)

               call fdinc1 (done,first,epsa,epsrf,objf,info,iter,
     *                     itmax,cdest,fdest,sdest,errbnd,f1,f2,h,hopt,
     *                     hphi)
c
c              +             until     done
               if (.not. done) go to 60
c
c              ---------------------------------------------------------
c              exit for this variable.
c              ---------------------------------------------------------
               gdiff = cdest
               x(j) = xj
c
               error = abs(gdiff-gj)/(gsize+one)
               if (error.ge.emax) then
                  emax = error
                  jmax = j
               end if
c
               ok = error .le. oktol
               if (ok) then
                  key = lgood
                  ngood = ngood + 1
               else
                  key = lbad
                  nwrong = nwrong + 1
               end if
c
               if (msglvl.gt.0) then
c
c                 zero elements are not printed.
c
                  const = ok.and.info .eq. 1.and.abs(gj) .lt. epspt8
                  if (.not. const) then
                     if (headng) then
                        write (rec,fmt=99995)
                        call x04bay(iprint,4,rec)
                        headng = .false.
                     end if
                     if (ok) then
                        write (rec,fmt=99994) j, xj, hopt, gj, gdiff,
     *                    key, iter
                     else
                        write (rec,fmt=99993) j, xj, hopt, gj, gdiff,
     *                    key, iter, result(info)
                     end if
                     call x04baf(iprint,rec(1))
                  end if
               end if
            end if
   80    continue
c

c        done.

         inform = 0
         if (msglvl.gt.0) then
            if (nwrong.eq.0) then
               write (rec,fmt=99992) ngood, ncheck, j1, j2
               call x04bay(iprint,3,rec)
            else
               write (rec,fmt=99991) nwrong, ncheck, j1, j2
               call x04bay(iprint,3,rec)
            end if
            write (rec,fmt=99990) emax, jmax
            call x04bay(iprint,3,rec)
         end if
         if (error.ge.point9) inform = 1
      end if
c
      call dcopy(n,grad,1,gradu,1)
c
      return
c
  100 inform = mode
      return
c
c
c     end of  e04xax. (chkgrd)
c
99999 format (//' verification of the objective gradients.',/' -------',
     *       '---------------------------------')
99998 format (/' the objective gradients seem to be ok.')
99997 format (/' xxx  the objective gradients seem to be incorrect.')
99996 format (/' directional derivative of the objective',1p,d18.8,
     *       /' difference approximation               ',1p,d18.8)
99995 format (//4x,'j',4x,'x(j)',5x,'dx(j)',11x,'g(j)',11x,'difference',
     *       ' approxn  itns',/)
99994 format (i5,1p,2d10.2,1p,2d18.8,2x,a4,i6)
99993 format (i5,1p,2d10.2,1p,2d18.8,2x,a4,i6,2x,a18)
99992 format (/i7,'  objective gradients out of the',i6,/9x,'set in co',
     *       'ls',i6,'  through',i6,'  seem to be ok.')
99991 format (/' xxx  there seem to be',i6,'  incorrect objective grad',
     *       'ients out of the',i6,/8x,'set in cols',i6,'  through',i6)
99990 format (/' the largest relative error was',1p,d12.2,'   in eleme',
     *       'nt',i6,/)
99989 format (/' no gradient elements assigned.')
      end

      subroutine fdinc (inform,msglvl,lvlder,n,bigbnd,
     *                  epsrf,fdnorm,objf,objfun,needc,bl,bu,
     *                  grad,gradu,hforwd,hcntrl,x,y,
     *                  iuser,user)

c     computes difference intervals for the missing gradients of
c     f(x). intervals are computed using a procedure that
c     usually requires about two function evaluations if the function
c     is well scaled.  central-difference gradients are obtained as a
c     by-product of the algorithm.
c
c     on entry...
c     objf contain the problem functions at the point x.
c     an element of grad not equal to rdummy signifies a known
c     gradient value.  such values are not estimated by differencing.
c     gradu have dummy elements in the same positions as gradu.
c
c     on exit...
c     grad contain central-difference derivative estimates.
c     elements of gradu are unaltered except for those
c     corresponding to constant derivatives, which are given the same
c     values as grad.

      implicit none

      double precision rdummy
      parameter (rdummy=-11111.0d+0)
      double precision factor
      parameter (factor=0.97d+0)
      double precision zero, half, one
      parameter (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision two, four, ten
      parameter (two=2.0d+0,four=4.0d+0,ten=1.0d+1)

      double precision bigbnd, epsrf, fdnorm, objf
      integer inform, ldcj, lvlder, msglvl, n

      double precision bl(n), bu(n), grad(n), gradu(n), hcntrl(*),
     *                  hforwd(*), user(*), x(n), y(n)
      integer iuser(*), needc(*)

      external          objfun

      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  ncdiff, nfdiff, nout
      logical           npdbg

      integer inpdbg(5)

      double precision biglow, bigupp, cdest, cjdiff, d, dx, epsa,
     *                  errbnd, errmax, errmin, f1, f2, fdest, fx,
     *                  gdiff, h, hcd, hfd, hmax, hmin, hopt, hphi,
     *                  objf2, sdest, signh, stepbl, stepbu, sumeps,
     *                  sumsd, test, xj, yj
      integer i, info, irow1, irow2, iter, itmax, j, mode,
     *                  nccnst, ncolj, nfcnst, nstate
      logical           done, first, headng

      character*80      rec(4)

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04uc/inpdbg, npdbg

      inform = 0


      if (lfdset.eq.0) then
         if (msglvl.gt.0) then
            write (rec,fmt=99999)
            call x04bay(iprint,4,rec)
         end if
c
         nstate = 0
         itmax = 3
         mode = 0
c
         nccnst = 0
         nfcnst = 0
         headng = .true.
c
         fdnorm = zero
         iuser(2) = 0

c        for each column of the jacobian augmented by the transpose of
c        the objective gradient, rows irow1 thru irow2 are searched for
c        missing elements.

c
         biglow = -bigbnd
         bigupp = bigbnd
c
         do 60 j = 1, n
            xj = x(j)
            ncolj = 0
            sumsd = zero
            sumeps = zero
            hfd = zero
            hcd = zero
            hmax = zero
            hmin = one/epspt3
            errmax = zero
            errmin = zero
c
            stepbl = biglow
            stepbu = bigupp
            if (bl(j).gt.biglow) stepbl = bl(j) - xj
            if (bu(j).lt.bigupp) stepbu = bu(j) - xj
c
            signh = one
            if (half*(stepbl+stepbu).lt.zero) signh = -one
c
                  test = gradu(j)

               if (test.eq.rdummy) then
c                 ======================================================
c                 get the difference interval for this element.
c                 ======================================================
                  ncolj = ncolj + 1
                     fx = objf
                     epsa = epsrf*(one+abs(fx))
c                 ------------------------------------------------------
c                 find a finite-difference interval by iteration.
c                 ------------------------------------------------------
                  iter = 0
                  hopt = two*(one+dabs(xj))*dsqrt(epsrf)
                  h = signh*ten*hopt
                  cdest = zero
                  sdest = zero
                  first = .true.
c
   20             x(j) = xj + h
                     call objfun(mode,n,x,f1,gradu,nstate,iuser,user)

c
                  x(j) = xj + h + h
                     call objfun(mode,n,x,f2,gradu,nstate,iuser,user)

c
                  call fdinc1 (done,first,epsa,epsrf,fx,info,iter,
     *                        itmax,cdest,fdest,sdest,errbnd,f1,f2,h,
     *                        hopt,hphi)
c
c                 +                until     done
                  if (.not. done) go to 20

                     grad(j) = cdest
                     if (info.eq.1.or.info.eq.2) then
                        nfcnst = nfcnst + 1
                        nfdiff = nfdiff - 1
                        gradu(j) = -rdummy
                     end if

                  sumsd = sumsd + abs(sdest)
                  sumeps = sumeps + epsa
                  if (hopt.gt.hmax) then
                     hmax = hopt
                     errmax = errbnd
                  end if
                  if (hopt.lt.hmin) then
                     hmin = hopt
                     errmin = errbnd
                  end if
c
                  if (info.eq.0) hcd = max(hcd,hphi)
               end if

c
            if (ncolj.gt.0) then
               if (hmin.gt.hmax) then
                  hmin = hmax
                  errmin = errmax
               end if
c
               if (four*sumeps.lt.hmin*hmin*sumsd) then
                  hfd = hmin
                  errmax = errmin
               else if (four*sumeps.gt.hmax*hmax*sumsd) then
                  hfd = hmax
               else
                  hfd = two*dsqrt(sumeps/sumsd)
                  errmax = two*dsqrt(sumeps*sumsd)
               end if
c
               if (hcd.eq.zero) hcd = ten*hfd
c
               if (msglvl.gt.0) then
                  if (headng) then
                     write (rec,fmt=99998)
                     call x04bay(iprint,4,rec)
                  end if
                  write (rec,fmt=99997) j, xj, hfd, hcd, errmax
                  call x04baf(iprint,rec(1))
                  headng = .false.
               end if
               fdnorm = max(fdnorm,hfd)
               hforwd(j) = hfd/(one+abs(xj))
               hcntrl(j) = hcd/(one+abs(xj))
            end if
            x(j) = xj
   60    continue
c
         if (nccnst+nfcnst.gt.0) then
c
c           check that the constants have been set properly by
c           evaluating the gradients at a strange (but feasible) point.
c
            d = one/n
c
            do 80 j = 1, n
               xj = x(j)
               stepbl = -one
               stepbu = one
               if (bl(j).gt.biglow) stepbl = max(stepbl,bl(j)-xj)
               if (bu(j).lt.bigupp.and.bu(j).gt.bl(j))
     *             stepbu = min(stepbu,bu(j)-xj)
c
               if (half*(stepbl+stepbu).lt.zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if
               d = factor*d
   80       continue

            call objfun(mode,n,y,objf2,gradu,nstate,iuser,user)

c           loop over each of the elements of  x.

            do 140 j = 1, n
               yj = y(j)
               dx = half*(x(j)-yj)
               y(j) = yj + dx

c              now check the objective gradient element.

               if (gradu(j).eq.-rdummy) then

                  call objfun(mode,n,y,f1,gradu,nstate,iuser,user)

                  gdiff = (f1-objf2)/dx
                  if (gdiff.eq.grad(j)) then
                     gradu(j) = gdiff
                  else
                     gradu(j) = rdummy
                     nfdiff = nfdiff + 1
                     nfcnst = nfcnst - 1
                  end if
               end if

               y(j) = yj
  140       continue
c
            if (msglvl.gt.0) then
               if (lvlder.lt.2.and.nccnst.gt.0) then
                  write (rec,fmt=99996) nccnst
                  call x04bay(iprint,2,rec)
               end if
               if (lvlder.ne.1.and.nfcnst.gt.0) then
                  write (rec,fmt=99995) nfcnst
                  call x04bay(iprint,2,rec)
               end if
            end if
c
            if (ncdiff.eq.0.and.lvlder.lt.2) then
               if (lvlder.eq.0) lvlder = 2
               if (lvlder.eq.1) lvlder = 3
               if (msglvl.gt.0) then
                  write (rec,fmt=99994) lvlder
                  call x04bay(iprint,4,rec)
               end if
            end if
c
            if (nfdiff.eq.0.and.lvlder.ne.1) then
               if (lvlder.eq.0) lvlder = 1
               if (lvlder.eq.2) lvlder = 3
               if (msglvl.gt.0) then
                  write (rec,fmt=99993) lvlder
                  call x04bay(iprint,4,rec)
               end if
            end if
         end if
      else if (lfdset.eq.2) then
c
c        the user has supplied hforwd and hcntrl.
c        check for wild values.
c
         do 160 j = 1, n
            if (hforwd(j).le.zero) then
               if (msglvl.gt.0) then
                  write (rec,fmt=99992) j, hforwd(j), epspt5
                  call x04baf(iprint,rec(1))
               end if
               hforwd(j) = epspt5
            end if
  160    continue
         do 180 j = 1, n
            if (hcntrl(j).le.zero) then
               if (msglvl.gt.0) then
                  write (rec,fmt=99991) j, hcntrl(j), epspt3
                  call x04baf(iprint,rec(1))
               end if
               hcntrl(j) = epspt3
            end if
  180    continue
      end if
c
      return

c     end of  e04xay. (chfd)

99999 format (//' computation of the finite-difference intervals',/' -',
     *       '---------------------------------------------')
99998 format (//'    j      x(j)   forward dx(j)   central dx(j)      ',
     *       'error est.',/)
99997 format (i5,1p,d10.2,1p,d16.6,1p,2d16.6)
99996 format (/i5,'  constant constraint gradient elements assigned.')
99995 format (/i5,'  constant  objective gradient elements assigned.')
99994 format (//' all missing jacobian elements are constants.',/' der',
     *       'ivative level increased to ',i4)
99993 format (//' all missing objective gradients are constants.',/' d',
     *       'erivative level increased to ',i4)
99992 format (' xxx  ',i4,'-th difference interval ',1p,d10.2,' replac',
     *       'ed by ',1p,d10.2)
99991 format (' xxx  ',i4,'-th central-difference interval ',1p,d10.2,
     *       ' replaced by ',1p,d10.2)

      end


      subroutine fdinc1 (done,first,epsa,epsr,fx,inform,iter,itmax,
     *                   cdest,fdest,sdest,errbnd,f1,f2,h,hopt,hphi)

c     implements algorithm  fd, the method described in
c     gill, p.e., murray, w., saunders, m.a., and wright, m. h.,
c     computing forward-difference intervals for numerical optimization,
c     siam journal on scientific and statistical computing, vol. 4,
c     pp. 310-321, june 1983.
c
c     the procedure is based on finding an interval (hphi) that
c     produces an acceptable estimate of the second derivative, and
c     then using that estimate to compute an interval that should
c     produce a reasonable forward-difference approximation.
c
c     one-sided difference estimates are used to ensure feasibility with
c     respect to an upper or lower bound on x. if x is close to an upper
c     bound, the trial intervals will be negative. the final interval is
c     always positive.
c
c     fdinc has a reverse communication
c     control structure, i.e., all evaluations of the function occur
c     outside this routine. the calling routine repeatedly calls fdinc1
c     after computing the indicated function values.

      double precision bndlo, bndup
      parameter (bndlo=1.0d-3,bndup=1.0d-1)
      double precision zero, sixth, fourth
      parameter (zero=0.0d+0,sixth=1.6d-1,fourth=2.5d-1)
      double precision half, two
      parameter (half=5.0d-1,two=2.0d+0)
      double precision three, four, ten
      parameter (three=3.0d+0,four=4.0d+0,ten=1.0d+1)
      double precision cdest, epsa, epsr, errbnd, f1, f2, fdest, fx, h,
     *                  hopt, hphi, sdest
      integer inform, iter, itmax
      logical done, first
      integer iprint, isumm, lines1, lines2, nout
      double precision afdmin, cdsave, err1, err2, fdcerr, fdest2,
     *                  fdsave, hsave, oldcd, oldh, oldsd, rho, sdcerr,
     *                  sdsave
      logical           ce1big, ce2big, overfl, te2big

      double precision adivb 
      external          adivb 
      intrinsic         abs, max, min, sqrt

      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      save              cdsave, fdsave, hsave, oldh, rho, sdsave,
     *                  ce1big, ce2big, te2big

c     bndlo, bndup, and rho control the logic of the routine.
c     bndlo and bndup are the lower and upper bounds that define an
c     acceptable value of the bound on the relative condition error in
c     the second derivative estimate.
c
c     the scalar rho is the factor by which the interval is multiplied
c     or divided, and also the multiple of the well-scaled interval
c     that is used as the initial trial interval.

      iter = iter + 1
c
c     compute the forward-,  backward-,  central-  and second-order
c     difference estimates.
c
      fdest = adivb (f1-fx,h,overfl)
      fdest2 = adivb (f2-fx,two*h,overfl)
c
      oldcd = cdest
      cdest = adivb (four*f1-three*fx-f2,two*h,overfl)
c
      oldsd = sdest
      sdest = adivb (fx-two*f1+f2,h*h,overfl)
c
c     compute  fdcerr  and  sdcerr,  bounds on the relative condition
c     errors in the first and second derivative estimates.
c
      afdmin = min(abs(fdest),abs(fdest2))
      fdcerr = adivb (epsa,half*abs(h)*afdmin,overfl)
      sdcerr = adivb (epsa,fourth*abs(sdest)*h*h,overfl)

c     select the correct case.

      if (first) then

c        first time through.
c        check whether sdcerr lies in the acceptable range.
c        ------------------------------------------------------------
         first = .false.
         done = sdcerr .ge. bndlo.and.sdcerr .le. bndup
         te2big = sdcerr .lt. bndlo
         ce2big = sdcerr .gt. bndup
         ce1big = fdcerr .gt. bndup

         if (.not. ce1big) then
            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if

         rho = epsr**(-sixth)/four
         if (te2big) then

c           the truncation error may be too big  (same as saying
c           sdcerr is too small).  decrease the trial interval.

            rho = ten*rho
            oldh = h
            h = h/rho
         else if (ce2big) then
c
c           sdcerr is too large.  increase the trial interval.
c
            oldh = h
            h = h*rho
         end if
      else if (ce2big) then

c        during the last iteration,  the trial interval was
c        increased in order to decrease sdcerr.

         if (ce1big.and.fdcerr.le.bndup) then
            ce1big = .false.
            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if
c
c        if sdcerr is small enough, accept h.  otherwise,
c        increase h again.
c
         done = sdcerr .le. bndup
         if (.not. done) then
            oldh = h
            h = h*rho
         end if
      else if (te2big) then

c        during the last iteration,  the interval was decreased in order
c        to reduce the truncation error.

         done = sdcerr .gt. bndup
         if (done) then
c
c           sdcerr has jumped from being too small to being too
c           large.  accept the previous value of h.

            h = oldh
            sdest = oldsd
            cdest = oldcd
         else

c           test whether fdcerr is sufficiently small.

            if (fdcerr.le.bndup) then
               ce1big = .false.
               hsave = h
               fdsave = fdest
               cdsave = cdest
               sdsave = sdest
            end if
c
c           check whether sdcerr is in range.
c
            done = sdcerr .ge. bndlo
c
            if (.not. done) then
c
c              sdcerr is still too small, decrease h again.
c
               oldh = h
               h = h/rho
            end if
         end if
c
      end if

c     we have either finished or have a new estimate of h.

      if (done) then
c
c        sufficiently good second-derivative estimate found.
c        compute the optimal interval.
c
         hphi = dabs(h)
         hopt = two*dsqrt(epsa)/dsqrt(dabs(sdest))
c
c        err1 is the error bound on the forward-difference estimate
c        with the final value of h.  err2 is the difference of fdest
c        and the central-difference estimate with hphi.
c
         err1 = hopt*dabs(sdest)
         err2 = dabs(fdest-cdest)
         errbnd = max(err1,err2)
c
c        set inform = 4  if the forward- and central-difference
c        estimates are not close.
c
         inform = 0
         if (errbnd.gt.half*abs(fdest)) inform = 4
      else

c        check whether the maximum number of iterations has been
c        exceeded.  if not, exit.

         done = iter .ge. itmax
         if (done) then
            if (ce1big) then
c
c              fdcerr was never small.  probably a constant function.
c
               inform = 1
               hphi = hopt
               fdest = zero
               cdest = zero
               sdest = zero
               errbnd = zero
            else if (ce2big) then
c
c              fdcerr was small,  but sdcerr was never small.
c              probably a linear or odd function.
c
               inform = 2
               hphi = abs(hsave)
               hopt = hphi
               fdest = fdsave
               cdest = cdsave
               sdest = zero
               errbnd = two*epsa/hopt
            else
c
c              the only remaining case occurs when the second
c              derivative is changing too rapidly for an adequate
c              interval to be found (sdcerr remained small even
c              though h was decreased itmax times).
c
               inform = 3
               hphi = abs(hsave)
               hopt = hphi
               fdest = fdsave
               cdest = cdsave
               sdest = sdsave
               errbnd = hopt*abs(sdest)/two + two*epsa/hopt
            end if
         end if
      end if

99999 format (/' //e04xaz//  itn ',i3,' fx     h',11x,1p,2d16.6,/' //e',
     *       '04xaz//  f1      fdest',14x,1p,2d16.6,/' //e04xaz//  f2 ',
     *       '     fdest2',13x,1p,2d16.6,/' //e04xaz//  cdest   sdest',
     *       14x,1p,2d16.6,/' //e04xaz//  fdcerr  sdcerr',13x,1p,2d16.6)
99998 format (' //e04xaz//  ce1big  ce2big  te2big',5x,3l2)
99997 format (' //e04xaz//  inform  hopt    errbnd',i5,1p,2d16.6)
      end

      subroutine evalfd (centrl,inform,n,bigbnd,cdint,
     *                  fdint,fdnorm,objf,objfun,needc,bl,bu,
     *                  grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)

c evaluates any missing gradients.

      implicit none

      double precision rdummy
      parameter (rdummy=-11111.0d+0)
      double precision zero, half, one
      parameter (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision three, four
      parameter (three=3.0d+0,four=4.0d+0)
      double precision bigbnd, cdint, fdint, fdnorm, objf
      integer inform, lenw, n
      logical           centrl
      double precision bl(n), bu(n), 
     *                  grad(n), gradu(n), hcntrl(n),
     *                  hforwd(n), user(*), w(lenw), x(n)
      integer iuser(*), needc(*)
      external          objfun
      double precision epspt3, epspt5, epspt8, epspt9
      integer lfdset, lvldif, ncdiff, nfdiff
      logical           npdbg
      integer inpdbg(5)
      double precision biglow, bigupp, delta, objf1, objf2, stepbl,
     *                  stepbu, xj
      integer i, j, mode, nstate
      intrinsic         abs, max
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04uc/inpdbg, npdbg

      inform = 0
      iuser(2) = 0

c     use the pre-assigned difference intervals to approximate the
c     derivatives.

c     use either the same interval for each element (lfdset = 1),
c     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).

      nstate = 0
      mode = 0
c
      biglow = -bigbnd
      bigupp = bigbnd
c
      fdnorm = zero
c
      do 80 j = 1, n
         xj = x(j)

         if (gradu(j).eq.rdummy) then
            stepbl = biglow
            stepbu = bigupp
            if (bl(j).gt.biglow) stepbl = bl(j) - xj
            if (bu(j).lt.bigupp) stepbu = bu(j) - xj
c
            if (centrl) then
               if (lfdset.eq.1) then
                  delta = cdint
               else
                  delta = hcntrl(j)
               end if
            else
               if (lfdset.eq.1) then
                  delta = fdint
               else
                  delta = hforwd(j)
               end if
            end if
c
            delta = delta*(one+abs(xj))
            fdnorm = max(fdnorm,delta)
            if (half*(stepbl+stepbu).lt.zero) delta = -delta
c
            x(j) = xj + delta

            if (gradu(j).eq.rdummy) then
               call objfun(mode,n,x,objf1,gradu,nstate,iuser,user)
            end if
c
            if (centrl) then
c              ---------------------------------------------------------
c              central differences.
c              ---------------------------------------------------------
               x(j) = xj + delta + delta

               if (gradu(j).eq.rdummy) then
                  call objfun(mode,n,x,objf2,gradu,nstate,iuser,user)
                  grad(j) = (four*objf1-three*objf-objf2)/(delta+delta)
               end if
            else
c              ---------------------------------------------------------
c              forward differences.
c              ---------------------------------------------------------
               if (gradu(j).eq.rdummy) grad(j) = (objf1-objf)/delta
c
            end if
         end if
         x(j) = xj
c
   80 continue

c     end of  e04uds. (npfd)

      end

      subroutine e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)

c     e04uch  finds steps palfa1, palfa2 such that
c        x + palfa1*p  reaches a linear constraint that is currently not
c                      in the working set but is satisfied.
c        x + palfa2*p  reaches a linear constraint that is currently not
c                      in the working set but is violated.
c     the constraints are perturbed by an amount featol, so that palfa1
c     is slightly larger than it should be,  and palfa2 is slightly
c     smaller than it should be.  this gives some leeway later when the
c     exact steps are computed by e04ucg.
c
c     constraints in the working set are ignored  (istate(j) .ge. 1).
c
c     if negstp is true, the search direction will be taken to be  - p.
c
c
c     values of istate(j)....
c
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     the values  -2  and  -1  do not occur once a feasible point has
c     been found.

      integer lcmdbg
      parameter (lcmdbg=5)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision bigalf, bigbnd, palfa1, palfa2, pnorm
      integer jadd1, jadd2, n, nctotl
      logical           firstv, negstp
      double precision anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer istate(nctotl)
      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lines1, lines2, nout
      logical           cmdbg
      integer icmdbg(lcmdbg)
      double precision absatp, atp, atx, res, rownrm
      integer i, j, js
      logical           lastv
      character*120     rec(3)
      intrinsic         abs
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg

      if (cmdbg.and.icmdbg(3).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,3,rec)
      end if
      lastv = .not. firstv
      jadd1 = 0
      jadd2 = 0
      palfa1 = bigalf
c
      palfa2 = zero
      if (firstv) palfa2 = bigalf
c
      do 20 j = 1, nctotl
         js = istate(j)
         if (js.le.0) then
            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               rownrm = one
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               rownrm = one + anorm(i)
            end if
            if (negstp) atp = -atp
c
            if (abs(atp).le.epspt9*rownrm*pnorm) then
c
c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.
c
               res = -one
            else if (atp.le.zero.and.js.ne.-2) then
c              ---------------------------------------------------------
c              a'x  is decreasing and the lower bound is not violated.
c              ---------------------------------------------------------
c              first test for smaller palfa1.
c
               absatp = -atp
               if (bl(j).gt.(-bigbnd)) then
                  res = atx - bl(j) + featol(j)
                  if (bigalf*absatp.gt.abs(res)) then
                     if (palfa1*absatp.gt.res) then
                        palfa1 = res/absatp
                        jadd1 = j
                     end if
                  end if
               end if
c
               if (js.eq.-1) then
c
c                 the upper bound is violated.  test for either larger
c                 or smaller palfa2, depending on the value of firstv.
c
                  res = atx - bu(j) - featol(j)
                  if (bigalf*absatp.gt.abs(res)) then
                     if (firstv.and.palfa2*absatp.gt.res .or.
     *                   lastv.and.palfa2*absatp.lt.res) then
                        palfa2 = res/absatp
                        jadd2 = j
                     end if
                  end if
               end if
            else if (atp.gt.zero.and.js.ne.-1) then
c              ---------------------------------------------------------
c              a'x  is increasing and the upper bound is not violated.
c              ---------------------------------------------------------
c              test for smaller palfa1.
c
               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx + featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (palfa1*atp.gt.res) then
                        palfa1 = res/atp
                        jadd1 = j
                     end if
                  end if
               end if
c
               if (js.eq.-2) then
c
c                 the lower bound is violated.  test for a new palfa2.
c
                  res = bl(j) - atx - featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (firstv.and.palfa2*atp.gt.res.or.lastv .and.
     *                   palfa2*atp.lt.res) then
                        palfa2 = res/atp
                        jadd2 = j
                     end if
                  end if
               end if
            end if
c
            if (cmdbg.and.icmdbg(3).gt.0) then
               write (rec,fmt=99998) j, js, featol(j), res, atp, jadd1,
     *           palfa1, jadd2, palfa2
               call x04baf(iprint,rec(1))
            end if
         end if
   20 continue

c     end of  e04uch. (cmalf1)

99999 format (/'    j  js         featol        res             ap    ',
     *       ' jadd1       palfa1     jadd2       palfa2',/)
99998 format (i5,i4,3g15.5,2(i6,g17.7))
      end


      subroutine e04nbv(n,nrank,ldr,lenv,lenw,r,u,v,w,c,s)

c     e04nbv  modifies the  nrank*n  upper-triangular matrix  r  so that
c     q*(r + v*w')  is upper triangular,  where  q  is orthogonal,
c     v  and  w  are vectors, and the modified  r  overwrites the old.
c     q  is the product of two sweeps of plane rotations (not stored).
c     if required,  the rotations are applied to the nu columns of
c     the matrix  u.
c
c     the matrix v*w' is an (lenv) by (lenw) matrix.
c     the vector v is overwritten.

      implicit none

      integer ldr, lenv, lenw, n, nrank
      double precision c(n), r(ldr,*), s(n), u(n,*), v(n), w(n)
      integer j
      intrinsic         min

      j = min(lenv,nrank)
      if (nrank.gt.0) then

c        reduce  v to beta*e(j)  using a backward sweep of rotations
c        in planes (j-1, j), (j-2, j), ..., (1, j).

         call ssrotg ('fixed','backwards',j-1,v(j),v,1,c,s)

c        apply the sequence of rotations to r. this generates a spike in
c        the j-th row of r, which is stored in s.

         call sutsrs ('left',n,1,j,c,s,r,ldr)
c

c        form  beta*e(j)*w' + r.  this a spiked matrix, with a row
c        spike in row j.

         call daxpy(min(j-1,lenw),v(j),w,1,s,1)
         call daxpy(lenw-j+1,v(j),w(j),1,r(j,j),ldr)
c

c        eliminate the spike using a forward sweep of rotations in
c        planes (1, j), (2, j), ..., (j-1, j).

         call susqr('left',n,1,j,c,s,r,ldr)

      end if
c
c     end of  e04nbv. (cmr1md)
c
      end


      subroutine e04ucl(lsumry,unitq,n,nfree,nz,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,
     *                  cs1,cs2,gq1,gq2,hpq,rpq,qpmul,r,omega,zy,
     *                  wrk1,wrk2)

c     e04ucl  computes the bfgs update for the approximate hessian of
c     the lagrangian.  if the approximate curvature of the lagrangian
c     function is negative,  a nonnegative penalty vector omega(i) of
c     minimum two norm is computed such that the approximate curvature
c     of the augmented lagrangian will be positive. if no finite penalty
c     vector exists,  the bfgs update is performed with the approximate
c     curvature modified to be a small positive value.
c
c     on entry,  gq1 and gq2 contain the transformed objective gradients
c     at x1 and x2,  hpq contains  r'r(pq), the transformed hessian
c     times the transformed search direction.  the vectors gq1 and hpq
c     are not saved.  if the regular bfgs quasi-newton update could not
c     be performed, the first character of lsumry is loaded with 'm'.

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision tolg
      parameter (tolg=1.0d-1)
      double precision alfa, glf1, glf2, qpcurv
      integer ldcj1, ldcj2, ldr, ldzy, n, ncnln, nfree, nz
      logical           unitq
      character*5       lsumry
      double precision cs1(*), cs2(*), gq1(n), gq2(n),
     *                  hpq(n), omega(*), qpmul(*), r(ldr,*), rpq(n),
     *                  wrk1(n), wrk2(n), zy(ldzy,*)
      integer kx(n)
      double precision drmax, drmin, rcndbd, rfrobn, rhodmp, rhomax,
     *                  rhonrm, scale
      integer iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg
      integer inpdbg(5)
      double precision beta, curvl, eta, qi, qmax, qnorm, rtgtp, rtyts,
     *                  test, tinycl, trace1, trace2
      integer i, imax, j
      logical           overfl, ssbfgs
      character*80      rec(3)
      double precision dnrm2, adivb 
      integer idamx1 
      external          dnrm2, adivb , idamx1 
      intrinsic         abs, max, min, sqrt
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /fe04uc/inpdbg, npdbg

c     set curvl = (g2 - g1)'dx,  the approximate curvature along dx of
c     the (augmented) lagrangian.  at first, the curvature is not scaled
c     by the steplength alfa.

      curvl = glf2 - glf1
      tinycl = qpcurv*tolg
      ssbfgs = curvl .le. alfa*tinycl

      if (npdbg.and.inpdbg(1).gt.0) then
         write (rec,fmt=99999) ssbfgs, tinycl, curvl
         call x04bay(iprint,3,rec)
      end if


c     test if curvl is sufficiently positive.  if there are no nonlinear
c     constraints,  no update can be performed.

      if (curvl.lt.tinycl) then
         lsumry(1:1) = 'modified bfgs'
      end if
c
      if (npdbg.and.inpdbg(1).gt.0) then
         write (rec,fmt=99998) alfa, curvl
         call x04bay(iprint,3,rec)
      end if
c
      if (curvl.lt.tinycl) curvl = tinycl
c
      do 120 j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
  120 continue
c
      rtgtp = dsqrt(qpcurv)
      rtyts = dsqrt(alfa*curvl)
      eta = one
      if (ssbfgs) eta = rtyts/(rtgtp*alfa)
c
      trace1 = dnrm2(n,hpq,1)/rtgtp
      trace2 = dnrm2(n,wrk2,1)/(rtyts*eta)
      rfrobn = eta*dsqrt(dabs((rfrobn-trace1)*
     *                  (rfrobn+trace1)+trace2**2))
c

c     update the cholesky factor of  q'hq.

c     normalize the vector  rpq (= r(pq)).

      call dscal(n,(one/rtgtp),rpq,1)

c     do the self-scaled or regular bfgs update.
c     form the vector wrk1 = gamma * (gq2 - gq1) - beta * r'r*pq,
c     where  gamma = 1/sqrt(curv) = 1/sqrt((gq2 - gq1)'sq)

      call dscal (n,(one/rtgtp),hpq,1)

      if (ssbfgs) then
         do 140 j = 1, n
            call dscal (j,eta,r(1,j),1)
            wrk1(j) = wrk2(j)/rtyts - eta*hpq(j)
  140    continue
      else
         do 160 j = 1, n
            wrk1(j) = wrk2(j)/rtyts - hpq(j)
  160    continue
      end if
c
c     perform the update to  r = r + rpq*wrk1'.
c     rpq is overwritten. arrays gq1 and hpq are used to store the
c     sines and cosines defined by the plane rotations.
c
      call e04nbv(n,n,ldr,n,n,r,hpq,rpq,wrk1,gq1,hpq)

c     end of  e04ucl. (npupdt)

99999 format (/' //e04ucl// ssbfgs    min. curvl         curvl ',/' //',
     *       'e04ucl//   ',l4,1p,2d14.2)
99998 format (/' //e04ucl//          alfa         curvl ',
     *       /' //e04ucl//',1p,2d14.2)
99997 format (/' //e04ucl//   omega(imax)',/' //e04ucl//',1p,d14.2)
99996 format (/' //e04ucl//  penalty parameters = ')
99995 format (1p,5d15.6)
      end

      subroutine e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,bl,bu,rpq,rpq0,dx,gq,r,t,zy,work)

c     e04ucm   defines a point which lies on the initial working set for
c     the qp subproblem.  this routine is similar to e04nch except
c     that advantage is taken of the fact that the initial estimate of
c     the solution of the least-squares subproblem is zero.

      implicit none 

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision dxnorm, gdx
      integer ldaqp, ldr, ldt, ldzy, n, nactiv, ncqp, nctotl,
     *                  nfree, nlnx, nz
      logical           unitq
      double precision adx(*), aqp(ldaqp,*), bl(nctotl), bu(nctotl),
     *                  dx(n), gq(n), r(ldr,*), rpq(nlnx), rpq0(nlnx),
     *                  t(ldt,*), work(n), zy(ldzy,*)
      integer istate(nctotl), kactiv(n), kx(n)
      integer iprint, isumm, lines1, lines2, nout
      logical           npdbg
      integer inpdbg(5)
      double precision bnd
      integer i, j, k, nfixed, nr
      character*80      rec(2)
      double precision ddot1, dnrm2
      external          ddot1, dnrm2
      intrinsic         min
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04uc/inpdbg, npdbg

      nfixed = n - nfree
      gdx = zero
      call sload (n,zero,dx,1)
      call sload (nlnx,zero,rpq,1)
      call sload (nlnx,zero,rpq0,1)
c
      if (nactiv+nfixed.gt.0) then
c
c        set  work = residuals for constraints in the working set.
c        solve for  dx,  the smallest correction to  x  that gives a
c        point on the constraints in the working set.
c        set the fixed variables on their bounds,  solve the triangular
c        system  t*(dxy) = residuals,  and define  dx = y*(dxy).
c        use  (dxy)  to update  d(=pr)  as  d = d - r'( 0 ).
c                                                     (dxy)
c
         do 20 i = 1, nfixed
            j = kx(nfree+i)
            if (istate(j).le.3) then
               bnd = bl(j)
               if (istate(j).eq.2) bnd = bu(j)
               dx(j) = bnd
               work(nfree+i) = bnd
            else
               work(nfree+i) = zero
            end if
   20    continue
c
         do 40 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(nz+i) = bnd - ddot1 (n,aqp(k,1),ldaqp,dx)
   40    continue
c
         if (nactiv.gt.0) call e04nbt(1,ldt,nactiv,t(1,nz+1),work(nz+1))
         call dcopy(nactiv+nfixed,work(nz+1),1,dx(nz+1),1)
         if (nz.gt.0) call sload (nz,zero,dx,1)
c
         gdx = ddot1 (nactiv+nfixed,gq(nz+1),1,dx(nz+1))
c
         if (nz.lt.n) then
            call dgemv ('n',nz,n-nz,-one,r(1,nz+1),ldr,dx(nz+1),one,
     *                 rpq)
            if (nz.lt.nlnx) then
               nr = ldr
               if (nz+1.eq.n) nr = 1
               call dcopy (nlnx-nz,dx(nz+1),1,rpq(nz+1),1)
               call dscal (nlnx-nz,(-one),rpq(nz+1),1)
               call dtrmv ('u','n','n',nlnx-nz,r(nz+1,nz+1),nr,
     *                    rpq(nz+1),1)
               if (nlnx.lt.n) then
                  nr = ldr
                  if (nlnx+1.eq.n) nr = n - nz
                  call dgemv ('n',nlnx-nz,n-nlnx,-one,r(nz+1,nlnx+1),nr,
     *                       dx(nlnx+1),one,rpq(nz+1))
               end if
            end if
         end if

         call cmqmul (2,n,nz,nfree,ldzy,unitq,kx,dx,zy,work)
      end if


c     compute the 2-norm of  dx.
c     initialize  a*dx.

      dxnorm = dnrm2(n,dx,1)
      if (ncqp.gt.0) call dgemv('n',ncqp,n,one,aqp,ldaqp,dx,zero,adx)
c
      if (npdbg.and.inpdbg(2).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,2,rec)
         do 60 i = 1, n, 5
            write (rec,fmt=99998) (dx(j),j=i,min(i+4,n))
            call x04baf(iprint,rec(1))
   60    continue
      end if

c     end of  e04ucm. (npsetx)

99999 format (/' //e04ucm// variables after e04ucm ... ')
99998 format (5g12.3)
      end


      subroutine e04nbx(msglvl,nfree,nrowa,n,nclin,nctotl,bigbnd,named,
     *                  names,nactiv,istate,kactiv,kx,a,bl,bu,clamda,
     *                  rlamda,x)

c     e04nbx   creates the expanded lagrange multiplier vector clamda.
c     if msglvl .eq 1 or msglvl .ge. 10,  e04nbx prints  x,  a*x,
c     c(x),  their bounds, the multipliers, and the residuals (distance
c     to the nearer bound).
c
c     e04nbx is called by e04ncz, e04ucz and e04upz just before exiting.

      integer lcmdbg
      parameter (lcmdbg=5)
      double precision zero
      parameter (zero=0.0d+0)

      double precision bigbnd
      integer msglvl, n, nactiv, nclin, nctotl, nfree, nrowa
      logical           named

      double precision a(nrowa,*), bl(nctotl), bu(nctotl),
     *                  clamda(nctotl), rlamda(n), x(n)
      integer istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)

      integer iprint, isumm, lines1, lines2, nout
      logical           cmdbg
      integer icmdbg(lcmdbg)
      double precision b1, b2, res, res2, v, wlam
      integer ip, is, j, k, nfixed, nplin, nz
      character*2       ls
      character*5       id3
      character*8       id4

      character*2       lstate(7)
      character*5       id(3)
      character*80      rec(4)

      double precision ddot1
      external          ddot1
      intrinsic         abs
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
      data              id(1)/'varbl'/
      data              id(2)/'l con'/
      data              id(3)/'n con'/
      data              lstate(1)/'--'/, lstate(2)/'++'/
      data              lstate(3)/'fr'/, lstate(4)/'ll'/
      data              lstate(5)/'ul'/, lstate(6)/'eq'/
      data              lstate(7)/'tf'/

      nplin = n + nclin
      nz = nfree - nactiv

c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      call sload (nctotl,zero,clamda,1)
      nfixed = n - nfree
      do 20 k = 1, nactiv + nfixed
         if (k.le.nactiv) j = kactiv(k) + n
         if (k.gt.nactiv) j = kx(nz+k)
         clamda(j) = rlamda(k)
   20 continue

      if (msglvl.lt.10.and.msglvl.ne.1) return

      write (rec,fmt=99999)
      call x04bay(iprint,4,rec)
      id3 = id(1)

      do 40 j = 1, nctotl
         b1 = bl(j)
         b2 = bu(j)
         wlam = clamda(j)
         is = istate(j)
         ls = lstate(is+3)
         if (j.le.n) then

c           section 1 -- the variables  x.
c           ------------------------------
            k = j
            v = x(j)

         else if (j.le.nplin) then

c           section 2 -- the linear constraints  a*x.
c           -----------------------------------------
            if (j.eq.n+1) then
               write (rec,fmt=99998)
               call x04bay(iprint,4,rec)
               id3 = id(2)
            end if
c
            k = j - n
            v = ddot1 (n,a(k,1),nrowa,x)

         end if

c        print a line for the j-th variable or constraint.
c        -------------------------------------------------
         res = v - b1
         res2 = b2 - v
         if (abs(res).gt.abs(res2)) res = res2
         ip = 1
         if (b1.le.(-bigbnd)) ip = 2
         if (b2.ge.bigbnd) ip = ip + 2
         if (named) then

            id4 = names(j)
            if (ip.eq.1) then
               write (rec,fmt=99996) id4, ls, v, b1, b2, wlam, res
            else if (ip.eq.2) then
               write (rec,fmt=99995) id4, ls, v, b2, wlam, res
            else if (ip.eq.3) then
               write (rec,fmt=99994) id4, ls, v, b1, wlam, res
            else
               write (rec,fmt=99993) id4, ls, v, wlam, res
            end if
            call x04baf(iprint,rec(1))

         else

            if (ip.eq.1) then
               write (rec,fmt=99992) id3, k, ls, v, b1, b2, wlam, res
            else if (ip.eq.2) then
               write (rec,fmt=99991) id3, k, ls, v, b2, wlam, res
            else if (ip.eq.3) then
               write (rec,fmt=99990) id3, k, ls, v, b1, wlam, res
            else
               write (rec,fmt=99989) id3, k, ls, v, wlam, res
            end if
            call x04baf(iprint,rec(1))
         end if
   40 continue

c     end of  e04nbx. (cmprt)

99999 format (//1x,'varbl',1x,'state',5x,'value',5x,'lower bound',3x,
     *       'upper bound',4x,'lagr mult',3x,'residual',/)
99998 format (//1x,'l con',1x,'state',5x,'value',5x,'lower bound',3x,
     *       'upper bound',4x,'lagr mult',3x,'residual',/)
99997 format (//1x,'n con',1x,'state',5x,'value',5x,'lower bound',3x,
     *       'upper bound',4x,'lagr mult',3x,'residual',/)
99996 format (1x,a4,4x,a2,1x,1p,3g14.6,1p,2g12.4)
99995 format (1x,a4,4x,a2,1x,1p,g14.6,5x,'none',5x,1p,g14.6,1p,2g12.4)
99994 format (1x,a4,4x,a2,1x,1p,2g14.6,5x,'none',5x,1p,2g12.4)
99993 format (1x,a4,4x,a2,1x,1p,g14.6,5x,'none',10x,'none',5x,1p,2g12.4)
99992 format (1x,a1,i3,4x,a2,1x,1p,3g14.6,1p,2g12.4)
99991 format (1x,a1,i3,4x,a2,1x,1p,g14.6,5x,'none',5x,1p,g14.6,1p,
     *       2g12.4)
99990 format (1x,a1,i3,4x,a2,1x,1p,2g14.6,5x,'none',5x,2g12.4)
99989 format (1x,a1,i3,4x,a2,1x,1p,g14.6,5x,'none',10x,'none',5x,1p,
     *       2g12.4)
      end


      subroutine e04ncr(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *                  ax,bl,bu,featol,x,work)

c     e04ncr  computes the following...
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.

      implicit none

      double precision zero
      parameter (zero=0.0d+0)

      double precision bigbnd, cvnorm, errmax
      integer jmax, n, nclin, nviol

      double precision ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), work(n+nclin), x(n)
      integer istate(n+nclin)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg
      integer ilsdbg(5)
      double precision biglow, bigupp, con, feasj, res, tolj
      integer i, is, j
      character*80      rec(2)
      double precision dnrm2
      integer idamx1 
      external          dnrm2, idamx1 
      intrinsic         abs
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

      biglow = -bigbnd
      bigupp = bigbnd

c     compute nviol,  the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint violations and
c     residuals of the constraints in the working set.

      nviol = 0
c
      do 40 j = 1, n + nclin
         feasj = featol(j)
         is = istate(j)
         res = zero
c
         if (is.ge.0.and.is.lt.4) then
            if (j.le.n) then
               con = x(j)
            else
               i = j - n
               con = ax(i)
            end if
c
            tolj = feasj
c
c           check for constraint violations.
c
            if (bl(j).gt.biglow) then
               res = bl(j) - con
               if (res.gt.feasj) nviol = nviol + 1
               if (res.gt.tolj) go to 20
            end if
c
            if (bu(j).lt.bigupp) then
               res = bu(j) - con
               if (res.lt.(-feasj)) nviol = nviol + 1
               if (res.lt.(-tolj)) go to 20
            end if
c
c           this constraint is satisfied,  but count the residual as a
c           violation if the constraint is in the working set.
c
            if (is.le.0) res = zero
            if (is.eq.1) res = bl(j) - con
            if (is.ge.2) res = bu(j) - con
            if (abs(res).gt.feasj) nviol = nviol + 1
         end if
   20    work(j) = res
   40 continue
c
      jmax = idamx1 (n+nclin,work)
      errmax = abs(work(jmax))
c
      if (lsdbg.and.ilsdbg(1).gt.0) then
         write (rec,fmt=99999) errmax, jmax
         call x04bay(iprint,2,rec)
      end if
c
      cvnorm = dnrm2(n+nclin,work,1)

c     end of  e04ncr. (lsfeas)

99999 format (/' //e04ncr//  the maximum violation is ',1p,d14.2,' in ',
     *       'constraint',i5)
      end


      subroutine e04ncl(hitcon,hitlow,linobj,unitgz,nclin,nrank,nrz,n,
     *                  ldr,jadd,numinf,alfa,ctp,ctx,xnorm,ap,ax,bl,bu,
     *                  gq,hz,p,res,r,x,work)

c     e04ncl  changes x to x + alfa*p and updates ctx, ax, res and gq
c     accordingly.
c
c     if a bound was added to the working set,  move x exactly on to it,
c     except when a negative step was taken (e04ucg may have had to move
c     to some other closer constraint.)

      implicit none 

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision alfa, ctp, ctx, xnorm
      integer jadd, ldr, n, nclin, nrank, nrz, numinf
      logical           hitcon, hitlow, linobj, unitgz

      double precision ap(*), ax(*), bl(*), bu(*), gq(*), hz(*), p(n),
     *                  r(ldr,*), res(*), work(*), x(n)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision bnd

      double precision dnrm2
      external          dnrm2

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

      call daxpy (n,alfa,p,1,x,1)
      if (linobj) ctx = ctx + alfa*ctp
c
      if (hitcon.and.jadd.le.n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa.ge.zero) x(jadd) = bnd
      end if
      xnorm = dnrm2(n,x,1)
c
      if (nclin.gt.0) call daxpy (nclin,alfa,ap,1,ax,1)
c
      if (nrz.le.nrank) then
         if (unitgz) then
            res(nrz) = res(nrz) - alfa*hz(nrz)
         else
            call daxpy (nrz,(-alfa),hz,1,res,1)
         end if
c
         if (numinf.eq.0) then
c
c           update the transformed gradient gq so that
c           gq = gq + alfa*r'(hz).
c                            (0 )
c
            if (unitgz) then
               call daxpy (n-nrz+1,alfa*hz(nrz),r(nrz,nrz),ldr,gq(nrz)
     *                                                        ,1)
            else
               call dcopy(nrz,hz,1,work,1)
               call dtrmv('u','t','n',nrz,r,ldr,work,1)
               if (nrz.lt.n) call dgemv('t',nrz,n-nrz,one,r(1,nrz+1),
     *                                  ldr,hz,zero,work(nrz+1))
c
               call daxpy (n,alfa,work,1,gq,1)
            end if
         end if
      end if

c     end of  e04ncl. (lsmove)

      end


      subroutine e04ucg(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfa,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  anorm,ap,ax,bl,bu,featol,p,x)

c     e04ucg finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).  two possible
c     steps are defined as follows...
c
c     alfa1   is the maximum step that can be taken without violating
c             one of the linear constraints that is currently satisfied.
c     alfa2   reaches a linear constraint that is currently violated.
c             usually this will be the furthest such constraint along p,
c             but if firstv = .true. it will be the first one along p.
c             this is used only when the problem has been determined to
c             be infeasible, and the sum of infeasibilities are being
c             minimized.  (alfa2  is not defined if numinf = 0.)
c
c     alfa will usually be the minimum of alfa1 and alfa2.
c     alfa could be negative (since we allow inactive constraints
c     to be violated by as much as featol).  in such cases, a
c     third possible step is computed, to find the nearest satisfied
c     constraint (perturbed by featol) along the direction  - p.
c     alfa  will be reset to this step if it is shorter.  this is the
c     only case for which the final step  alfa  does not move x exactly
c     onto a constraint (the one denoted by jadd).
c
c     constraints in the working set are ignored  (istate(j) ge 1).
c
c     jadd    denotes which linear constraint is reached.
c
c     hitlow  indicates whether it is the lower or upper bound that
c             has restricted alfa.
c
c     values of istate(j)....
c
c     - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     the values -2 and -1 do not occur once a feasible point has been
c     found.

      integer lcmdbg
      parameter (lcmdbg=5)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision alfa, atphit, bigalf, bigbnd, palfa, pnorm
      integer inform, jadd, n, nctotl, numinf
      logical           firstv, hitlow

      double precision anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer istate(nctotl)

      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lines1, lines2, nout
      logical           cmdbg

      integer icmdbg(lcmdbg)

      double precision absatp, alfa1, alfa2, apmax1, apmax2, atp, atp1,
     *                  atp2, atx, palfa1, palfa2, res, rownrm
      integer i, j, jadd1, jadd2, js, jsave1, jsave2
      logical           hlow1, hlow2, lastv, negstp, step2

      character*120     rec(4)

c     .. intrinsic functions ..
      intrinsic         abs, min

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg

c
      inform = 0
c

c     first pass -- find steps to perturbed constraints, so that
c     palfa1 will be slightly larger than the true step, and
c     palfa2 will be slightly smaller than it should be.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

c
      negstp = .false.
      call e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,palfa1,
     *            palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,featol,p,x)
c
      jsave1 = jadd1
      jsave2 = jadd2
c

c     second pass -- recompute step-lengths without perturbation.
c     amongst constraints that are less than the perturbed steps,
c     choose the one (of each type) that makes the largest angle
c     with the search direction.

      if (cmdbg.and.icmdbg(3).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,4,rec)
      end if
      alfa1 = bigalf
      alfa2 = zero
      if (firstv) alfa2 = bigalf
c
      apmax1 = zero
      apmax2 = zero
      atp1 = zero
      atp2 = zero
      hlow1 = .false.
      hlow2 = .false.
      lastv = .not. firstv
c
      do 20 j = 1, nctotl
         js = istate(j)
         if (js.le.0) then
            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               rownrm = one
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               rownrm = anorm(i) + one
            end if
c
            if (abs(atp).le.epspt9*rownrm*pnorm) then
c
c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.
c
               res = -one
            else if (atp.le.zero.and.js.ne.-2) then
c              ---------------------------------------------------------
c              a'x  is decreasing.
c              ---------------------------------------------------------
c              the lower bound is satisfied.  test for smaller alfa1.
c
               absatp = -atp
               if (bl(j).gt.(-bigbnd)) then
                  res = atx - bl(j)
                  if (palfa1*absatp.ge.res.or.j.eq.jsave1) then
                     if (apmax1*rownrm*pnorm.lt.absatp) then
                        apmax1 = absatp/(rownrm*pnorm)
                        alfa1 = res/absatp
                        jadd1 = j
                        atp1 = atp
                        hlow1 = .true.
                     end if
                  end if
               end if
c
               if (js.eq.-1) then
c
c                 the upper bound is violated.  test for either a bigger
c                 or smaller alfa2,  depending on the value of firstv.
c
                  res = atx - bu(j)
                  if ((firstv.and.palfa2*absatp.ge.res .or.
     *                lastv.and.palfa2*absatp.le.res)
     *               .or.j.eq.jsave2) then
                     if (apmax2*rownrm*pnorm.lt.absatp) then
                        apmax2 = absatp/(rownrm*pnorm)
                        if (absatp.ge.one) then
                           alfa2 = res/absatp
                        else if (res.lt.bigalf*absatp) then
                           alfa2 = res/absatp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2 = j
                        atp2 = atp
                        hlow2 = .false.
                     end if
                  end if
               end if
            else if (atp.gt.zero.and.js.ne.-1) then
c              ---------------------------------------------------------
c              a'x  is increasing and the upper bound is not violated.
c              ---------------------------------------------------------
c              test for smaller alfa1.
c
               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx
                  if (palfa1*atp.ge.res.or.j.eq.jsave1) then
                     if (apmax1*rownrm*pnorm.lt.atp) then
                        apmax1 = atp/(rownrm*pnorm)
                        alfa1 = res/atp
                        jadd1 = j
                        atp1 = atp
                        hlow1 = .false.
                     end if
                  end if
               end if
c
               if (js.eq.-2) then
c
c                 the lower bound is violated.  test for a new alfa2.
c
                  res = bl(j) - atx
                  if ((firstv.and.palfa2*atp.ge.res.or.lastv .and.
     *                palfa2*atp.le.res).or.j.eq.jsave2) then
                     if (apmax2*rownrm*pnorm.lt.atp) then
                        apmax2 = atp/(rownrm*pnorm)
                        if (atp.ge.one) then
                           alfa2 = res/atp
                        else if (res.lt.bigalf*atp) then
                           alfa2 = res/atp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2 = j
                        atp2 = atp
                        hlow2 = .true.
                     end if
                  end if
               end if
            end if
c
            if (cmdbg.and.icmdbg(3).gt.0) then
               write (rec,fmt=99998) j, js, featol(j), res, atp, jadd1,
     *           alfa1, jadd2, alfa2
               call x04baf(iprint,rec(1))
            end if
         end if
   20 continue
c

c     determine alfa, the step to be taken.

c     in the infeasible case, check whether to take the step alfa2
c     rather than alfa1...
c
      step2 = numinf .gt. 0.and.jadd2 .gt. 0
c
c     we do so if alfa2 is less than alfa1 or (if firstv is false)
c     lies in the range  (alfa1, palfa1)  and has a smaller value of
c     atp.
c
      step2 = step2.and.(alfa2.lt.alfa1.or.lastv.and.alfa2.le.
     *        palfa1.and.apmax2.ge.apmax1)
c
      if (step2) then
         alfa = alfa2
         palfa = palfa2
         jadd = jadd2
         atphit = atp2
         hitlow = hlow2
      else
         alfa = alfa1
         palfa = palfa1
         jadd = jadd1
         atphit = atp1
         hitlow = hlow1
c
c        if alfa1 is negative, the constraint to be added (jadd)
c        remains unchanged, but alfa may be shortened to the step
c        to the nearest perturbed satisfied constraint along  - p.
c
         negstp = alfa .lt. zero
         if (negstp) then
            call e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)
c
            if (cmdbg.and.icmdbg(1).gt.0) then
               write (rec,fmt=99997) alfa, palfa1
               call x04bay(iprint,4,rec)
            end if
c
            alfa = -min(abs(alfa),palfa1)
         end if
      end if
c
c     test for undefined or infinite step.
c
      if (jadd.eq.0) then
         alfa = bigalf
         palfa = bigalf
      end if
c
      if (alfa.ge.bigalf) inform = 3
      if (cmdbg.and.icmdbg(1).gt.0.and.inform.gt.0) then
         write (rec,fmt=99996) jadd, alfa
         call x04bay(iprint,4,rec)
      end if

c     end of  e04ucg. (cmalf)
c
99999 format (/' e04ucg  entered',/'    j  js         featol        re',
     *       's             ap     jadd1        alfa1     jadd2       ',
     *       ' alfa2 ',/)
99998 format (i5,i4,3g15.5,2(i6,g17.7))
99997 format (/' //e04ucg //  negative step',/' //e04ucg //           ',
     *       'alfa          palfa',/' //e04ucg //',2g15.4)
99996 format (/' //e04ucg //  unbounded step.',/' //e04ucg //  jadd   ',
     *       '        alfa',/' //e04ucg //  ',i4,g15.4)
      end


      subroutine e04ncq(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,ap,res,
     *                  hz,p,gq,cq,r,zy,work)

c     e04ncq  computes the following quantities for  e04ncz.
c     (1) the vector  (hz1) = (rz1)(pz1).
c         if x is not yet feasible,  the product is computed directly.
c         if  rz1 is singular,  hz1  is zero.  otherwise  hz1  satisfies
c         the equations
c                        rz1'hz1 = -gz1,
c         where  g  is the total gradient.  if there is no linear term
c         in the objective,  hz1  is set to  dz1  directly.
c     (2) the search direction p (and its 2-norm).  the vector p is
c         defined as  z*(pz1), where  (pz1)  depends upon whether or
c         not x is feasible and the nonsingularity of  (rz1).
c         if  numinf .gt. 0,  (pz1)  is the steepest-descent direction.
c         otherwise,  x  is the solution of the  nrz*nrz  triangular
c         system   (rz1)*(pz1) = (hz1).
c     (3) the vector ap,  where a is the matrix of linear constraints.

      implicit none 

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision ctp, pnorm
      integer lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf
      logical           linobj, singlr, unitgz, unitq

      double precision a(lda,*), ap(*), cq(*), gq(n), hz(*), p(n),
     *                  r(ldr,*), res(*), work(n), zy(ldzy,*)
      integer kx(n)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision gtp
      integer i, j

      character*80      rec(2)

      double precision ddot1, dnrm2
      external          ddot1, dnrm2

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

      if (singlr) then

c        the triangular factor for the current objective function is
c        singular,  i.e., the objective is linear along the last column
c        of z1.  this can only occur when unitgz is true.

         if (nrz.gt.1) then
            call dcopy(nrz-1,r(1,nrz),1,p,1)
            call dtrsv('u','n','n',nrz-1,r,ldr,p,1)
         end if
         p(nrz) = -one
c
         gtp = ddot1 (nrz,gq,1,p)
         if (gtp.gt.zero) call dscal(nrz,(-one),p,1)
c
         if (nrz.le.nrank) then
            if (numinf.eq.0) then
               if (unitgz) then
                  hz(nrz) = r(nrz,nrz)*p(nrz)
               else
                  call sload (nrz,(zero),hz,1)
               end if
            else
               hz(1) = r(1,1)*p(1)
            end if
         end if
      else

c        the objective is quadratic in the space spanned by z1.

         if (linobj) then
            if (unitgz) then
               if (nrz.gt.1) call sload (nrz-1,(zero),hz,1)
               hz(nrz) = -gq(nrz)/r(nrz,nrz)
            else
               call dcopy(nrz,gq,1,hz,1)
               call dscal(nrz,(-one),hz,1)
               call dtrsv('u','t','n',nrz,r,ldr,hz,1)
            end if
         else
            call dcopy(nrz,res,1,hz,1)
         end if
c
c        solve  rz1*pz1 = hz1.
c
         call dcopy(nrz,hz,1,p,1)
         call dtrsv('u','n','n',nrz,r,ldr,p,1)
c
      end if
c
c     compute  p = z1*pz1  and its norm.
c
      if (linobj) ctp = ddot1 (nrz,cq,1,p)
      pnorm = dnrm2(nrz,p,1)
c
      call cmqmul (1,n,nrz,nfree,ldzy,unitq,kx,p,zy,work)
c
      if (lsdbg.and.ilsdbg(2).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,2,rec)
         do 20 i = 1, n, 5
            write (rec,fmt=99997) (p(j),j=i,min(i+4,n))
            call x04baf(iprint,rec(1))
   20    continue
      end if
c
c     compute  ap.
c
      if (nclin.gt.0) then
         call dgemv ('n',nclin,n,one,a,lda,p,zero,ap)
         if (lsdbg.and.ilsdbg(2).gt.0) then
            write (rec,fmt=99998)
            call x04bay(iprint,2,rec)
            do 40 i = 1, n, 5
               write (rec,fmt=99997) (ap(j),j=i,min(i+4,n))
               call x04baf(iprint,rec(1))
   40       continue
         end if
      end if

c     end of  e04ncq. (lsgetp)

99999 format (/' //e04ncq//   p ... ')
99998 format (/' //e04ncq//  ap ... ')
99997 format (1p,5d15.5)
      end


      subroutine e04ncj(prbtyp,isdel,iter,jadd,jdel,msglvl,nactiv,nfree,
     *                  n,nclin,nrank,ldr,ldt,nz,nrz,istate,alfa,condrz,
     *                  condt,gfnorm,gzrnrm,numinf,suminf,ctx,ssq,ax,r,
     *                  t,x,work)

c     e04ncj  prints various levels of output for  e04ncz.
c
c           msg        cumulative result
c           ---        -----------------
c
c        le   0        no output.
c
c        eq   1        nothing now (but full output later).
c
c        eq   5        one terse line of output.
c
c        ge  10        same as 5 (but full output later).
c
c        ge  20        constraint status,  x  and  ax.
c
c        ge  30        diagonals of  t  and  r.
c
c
c     debug printing is performed depending on the logical variable
c     lsdbg.  lsdbg  is set true when  idbg  major iterations have
c     been performed. at this point,  printing is done according to
c     a string of binary digits of the form  svt  (stored in the
c     integer array  ilsdbg).
c
c     s  set 'on'  gives information from the maximum step routine
c                  e04ucg.
c     v  set 'on'  gives various vectors in  e04ncz  and its
c                  auxiliaries.
c     t  set 'on'  gives a trace of which routine was called and an
c                  indication of the progress of the run.

      implicit none 

      integer mline1, mline2
      parameter (mline1=50000,mline2=50000)

      double precision alfa, condrz, condt, ctx, gfnorm, gzrnrm, ssq,
     *                  suminf
      integer isdel, iter, jadd, jdel, ldr, ldt, msglvl, n,
     *                  nactiv, nclin, nfree, nrank, nrz, numinf, nz
      character*2       prbtyp

      double precision ax(*), r(ldr,*), t(ldt,*), work(n), x(n)
      integer istate(*)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision obj
      integer i, itn, j, k, kadd, kdel, nart, ndf
      logical           first, linobj, newset, prthdr
      character*2       ladd, ldel

      character*2       lstate(0:5)
      character*120     rec(4)

c     .. intrinsic functions ..
      intrinsic         min, mod

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

      data              lstate(0), lstate(1), lstate(2)/'  ', 'l ',
     *                  'u '/
      data              lstate(3), lstate(4), lstate(5)/'e ', 'f ',
     *                  'a '/

      if (msglvl.ge.15) then
         if (isumm.ge.0) then
            write (rec,fmt=99999) prbtyp, iter
            call x04bay(isumm,4,rec)
         end if
      end if
c
      if (msglvl.ge.5) then
c
         first = iter .eq. 0
         linobj = nrank .eq. 0
c
         itn = mod(iter,10000)
         ndf = mod(nrz,10000)
c
         nart = nz - nrz
c
         if (jdel.gt.0) then
            kdel = isdel
         else if (jdel.lt.0) then
            jdel = nart + 1
            kdel = 5
         else
            kdel = 0
         end if
c
         if (jadd.gt.0) then
            kadd = istate(jadd)
         else
            kadd = 0
         end if
c
         ldel = lstate(kdel)
         ladd = lstate(kadd)
c
         if (numinf.gt.0) then
            obj = suminf
         else
            obj = ssq + ctx
         end if

c        if necessary, print a header.
c        print a single line of information.

         if (isumm.ge.0) then
c           -----------------------------------
c           terse line for the monitoring file.
c           -----------------------------------
            newset = lines1 .ge. mline1
            prthdr = msglvl .ge. 15.or.first.or.newset
c
            if (prthdr) then
               if (linobj) then
                  write (rec,fmt=99998)
                  call x04bay(isumm,3,rec)
               else
                  write (rec,fmt=99996)
                  call x04bay(isumm,3,rec)
               end if
               lines1 = 0
            end if
c
            if (linobj) then
               write (rec,fmt=99994) itn, jdel, ldel, jadd, ladd, alfa,
     *           numinf, obj, n - nfree, nactiv, nart, ndf, gzrnrm,
     *           gfnorm, condt
               call x04baf(isumm,rec(1))
            else
               write (rec,fmt=99994) itn, jdel, ldel, jadd, ladd, alfa,
     *           numinf, obj, n - nfree, nactiv, nart, ndf, gzrnrm,
     *           gfnorm, condt, condrz
               call x04baf(isumm,rec(1))
            end if
            lines1 = lines1 + 1
         end if
c
         if (iprint.ge.0.and.isumm.ne.iprint) then
c           ------------------------------
c           terse line for the print file.
c           ------------------------------
            newset = lines2 .ge. mline2
            prthdr = first.or.newset
c
            if (prthdr) then
               write (rec,fmt=99997)
               call x04bay(iprint,3,rec)
               lines2 = 0
            end if
c
            write (rec,fmt=99995) itn, alfa, numinf, obj, gzrnrm
            call x04baf(iprint,rec(1))
            lines2 = lines2 + 1
         end if
c
         if (msglvl.ge.20) then
            if (isumm.ge.0) then
               write (rec,fmt=99993) prbtyp
               call x04bay(isumm,3,rec)
               write (rec,fmt=99992)
               call x04bay(isumm,2,rec)
               do 20 i = 1, n, 5
                  write (rec,fmt=99987) (x(j),istate(j),j=i,min(i+4,n))
                  call x04baf(isumm,rec(1))
   20          continue
               if (nclin.gt.0) then
                  write (rec,fmt=99991)
                  call x04bay(isumm,2,rec)
                  do 40 i = 1, nclin, 5
                     write (rec,fmt=99987) (ax(k),istate(n+k),k=i,
     *                 min(i+4,nclin))
                     call x04baf(isumm,rec(1))
   40             continue
               end if
c
               if (msglvl.ge.30) then
c                 ------------------------------------------------------
c                 print the diagonals of  t  and  r.
c                 ------------------------------------------------------
                  if (nactiv.gt.0) then
                     call dcopy(nactiv,t(nactiv,nz+1),ldt-1,work,1)
                     write (rec,fmt=99990) prbtyp
                     call x04bay(isumm,2,rec)
                     do 60 i = 1, nactiv, 5
                        write (rec,fmt=99986) (work(j),j=i,
     *                    min(i+4,nactiv))
                        call x04baf(isumm,rec(1))
   60                continue
                  end if
                  if (nrank.gt.0) then
                     write (rec,fmt=99989) prbtyp
                     call x04bay(isumm,2,rec)
                     do 80 i = 1, nrank, 5
                        write (rec,fmt=99986) (r(j,j),j=i,min(i+4,nrank)
     *                   )
                        call x04baf(isumm,rec(1))
   80                continue
                  end if
               end if
               write (rec,fmt=99988)
               call x04bay(isumm,3,rec)
            end if
         end if
      end if

c     end of  e04ncj. (lsprt)

99999 format (//' ',a2,' iteration',i5,/' =================')
99998 format (//' itn jdel  jadd      step ninf  sinf/objective  bnd  ',
     *       'lin  art   zr  norm gz  norm gf   cond t')
99997 format (//' itn     step ninf sinf/objective  norm gz')
99996 format (//' itn jdel  jadd      step ninf  sinf/objective  bnd  ',
     *       'lin  art   zr  norm gz  norm gf   cond t  cond rz')
99995 format (i4,1p,d9.1,i5,d15.6,d9.1)
99994 format (i4,i5,a1,i5,a1,1p,d9.1,i5,d16.8,4i5,4d9.1)
99993 format (/' values and status of the ',a2,' constraints',/' -----',
     *       '----------------------------------')
99992 format (/' variables...')
99991 format (/' general linear constraints...')
99990 format (/' diagonals of ',a2,' working set factor t')
99989 format (/' diagonals of ',a2,' triangle r         ')
99988 format (//' ----------------------------------------------------',
     *       '-------------------------------------------')
99987 format (1x,5(1p,d15.6,i5))
99986 format (1p,5d15.6)
      end



      subroutine e04ncp(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,res,featol,gq,cq,r,
     *                  x,wtinf,zy,wrk)

c     e04ncp  finds the number and weighted sum of infeasibilities for
c     the bounds and linear constraints.   an appropriate transformed
c     gradient vector is returned in  gq.
c
c     positive values of  istate(j)  will not be altered.  these mean
c     the following...
c
c               1             2           3
c           a'x = bl      a'x = bu     bl = bu
c
c     other values of  istate(j)  will be reset as follows...
c           a'x lt bl     a'x gt bu     a'x free
c              - 2           - 1           0
c
c     if  x  is feasible,  e04ncp computes the vector q(free)'g(free),
c     where  g  is the gradient of the the sum of squares plus the
c     linear term.  the matrix q is of the form
c                    (q(free)  0      ),
c                    (  0      i(fixed))
c     where  q(free)  is the orthogonal factor of  a(free)  and  a  is
c     the matrix of constraints in the working set.  the transformed
c     gradients are stored in gq.

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision bigbnd, suminf, tolrnk
      integer lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf, nz
      logical           linobj, singlr, unitgz, unitq
      character*2       prbtyp

      double precision a(lda,*), bl(*), bu(*), cq(*), featol(*), gq(n),
     *                  r(ldr,*), res(*), wrk(n), wtinf(*), x(n),
     *                  zy(ldzy,*)
      integer istate(*), kx(n)

      double precision biglow, bigupp, ctx, feasj, rownrm, s, weight
      integer j, k

      double precision ddot1, dnrm2
      integer isrank
      external          ddot1, dnrm2, isrank
c     .. intrinsic functions ..
      intrinsic         abs, min

      bigupp = bigbnd
      biglow = -bigbnd
c
      numinf = 0
      suminf = zero
      call sload (n,zero,gq,1)
c
      do 40 j = 1, n + nclin
         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot1 (n,a(k,1),lda,x)
            end if
            istate(j) = 0
c
c           see if the lower bound is violated.
c
            if (bl(j).gt.biglow) then
               s = bl(j) - ctx
               if (s.gt.feasj) then
                  istate(j) = -2
                  weight = -wtinf(j)
                  go to 20
               end if
            end if
c
c           see if the upper bound is violated.
c
            if (bu(j).ge.bigupp) go to 40
            s = ctx - bu(j)
            if (s.le.feasj) go to 40
            istate(j) = -1
            weight = wtinf(j)
c
c           add the infeasibility.
c
   20       numinf = numinf + 1
            suminf = suminf + abs(weight)*s
            if (j.le.n) then
               gq(j) = weight
            else
               call daxpy (n,weight,a(k,1),lda,gq,1)
            end if
         end if
   40 continue
c

c     install  gq,  the transformed gradient.

      singlr = .false.
      unitgz = .true.
c
      if (numinf.gt.0) then
         call cmqmul (6,n,nz,nfree,ldzy,unitq,kx,gq,zy,wrk)
         unitgz = .true.
      else if (numinf.eq.0.and.prbtyp.eq.'fp') then
         call sload (n,zero,gq,1)
      else
c
c        ready for the optimality phase.
c        set nrz so that rz1 is nonsingular.
c
         if (nrank.eq.0) then
            if (linobj) then
               call dcopy(n,cq,1,gq,1)
            else
               call sload (n,zero,gq,1)
            end if
            nrz = 0
         else
c
c           compute gq = - r' * (transformed residual)
c
            call sscmv (nrank,-one,res,1,gq,1)
            call dtrmv('u','t','n',nrank,r,ldr,gq,1)
            if (nrank.lt.n) call dgemv ('t',nrank,n-nrank,-one,
     *                                 r(1,nrank+1),ldr,res,zero,
     *                                 gq(nrank+1))
            if (linobj) call daxpy (n,one,cq,1,gq,1)
            unitgz = .false.
            rownrm = dnrm2(n,r(1,1),ldr)
            if (rownrm.le.tolrnk.or.abs(r(1,1)).le.rownrm*tolrnk) then
               nrz = 0
            else
               nrz = isrank(min(nrank,nz),r,ldr+1,tolrnk)
            end if
         end if
         singlr = .false.
      end if

c     end of  e04ncp. (lsgset)
c
      end


      subroutine e04uck(first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *                  toltny,alfa,alfbst,fbest,gbest)

c     e04uck  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     e04uck  requires both  f(alpha)  and  f'(alpha) to be evaluated at
c     points in the interval.  estimates of the minimizer are computed
c     using safeguarded cubic interpolation.
c
c     reverse communication is used to allow the calling program to
c     evaluate f and f'.  some of the parameters must be set or tested
c     by the calling program.  the remainder would ordinarily be local
c     variables.
c
c     input parameters (relevant to the calling program)
c     --------------------------------------------------
c
c     first         must be .true. on the first entry.
c                   it is subsequently altered by e04uck.
c
c     debug         specifies whether detailed output is wanted.
c
c     maxf          is an upper limit on the number of times e04uck is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).
c
c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by e04uck (see below).
c
c     alfmax        is the upper limit of the interval to be searched.
c
c     epsaf         is an estimate of the absolute precision in the
c                   computed value of f(0).
c
c     ftry, gtry    are the values of f, f'  at the new point
c                   alfa = alfbst + xtry.
c
c     g0            is the value of f'(0).  g0 must be negative.
c
c     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
c                   such that if f has already been evaluated at alfa,
c                   it will not be evaluated closer than tol(alfa).
c                   these values may be reduced by e04uck.
c
c     targtg        is the target value of abs(f'(alfa)). the search
c                   is terminated when
c                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
c
c     toltny        is the smallest value that tolabs is allowed to be
c                   reduced to.
c
c     output parameters (relevant to the calling program)
c     ---------------------------------------------------
c
c     imprvd        is .true. if the previous alfa was the best point so
c                   far.  any related quantities should be saved by the
c                   calling program (e.g., gradient arrays) before
c                   paying attention to the variable done.
c
c     done = .false.  means the calling program should evaluate
c                      ftry = f(alfa),  gtry = f'(alfa)
c                   for the new trial alfa, and re-enter e04uck.
c
c     done = .true.   means that no new alfa was calculated.  the value
c                   of inform gives the result of the search as follows
c
c                   inform = 1 means the search has terminated
c                              successfully with alfbst < alfmax.
c
c                   inform = 2 means the search has terminated
c                              successfully with alfbst = alfmax.
c
c                   inform = 3 means that the search failed to find a
c                              point of sufficient decrease in maxf
c                              functions, but a lower point was found.
c
c                   inform = 4 means alfmax is so small that a search
c                              should not have been attempted.
c
c                   inform = 5 is never set by e04uck.
c
c                   inform = 6 means the search has failed to find a
c                              useful step.  the interval of uncertainty
c                              is [0,b] with b < 2*tolabs. a minimizer
c                              lies very close to alfa = 0, or f'(0) is
c                              not sufficiently accurate.
c
c                   inform = 7 if no better point could be found after
c                              maxf  function calls.
c
c                   inform = 8 means the input parameters were bad.
c                              alfmax le toltny  or g0 ge zero.
c                              no function evaluations were made.
c
c     numf          counts the number of times e04uck has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).
c
c     alfa          is the point at which the next function ftry and
c                   derivative gtry must be computed.
c
c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever e04uck returns
c                   inform = 1 or 2 (and possibly 3).
c
c     fbest, gbest  will be the corresponding values of f, f'.
c
c
c     the following parameters retain information between entries
c     -----------------------------------------------------------
c
c     braktd        is .false. if f and f' have not been evaluated at
c                   the far end of the interval of uncertainty.  in this
c                   case, the point b will be at alfmax + tol(alfmax).
c
c     crampd        is .true. if alfmax is very small (le tolabs).  if
c                   the search fails, this indicates that a zero step
c                   should be taken.
c
c     extrap        is .true. if xw lies outside the interval of
c                   uncertainty.  in this case, extra safeguards are
c                   applied to allow for instability in the polynomial
c                   fit.
c
c     moved         is .true. if a better point has been found, i.e.,
c                   alfbst gt 0.
c
c     wset          records whether a second-best point has been
c                   determined it will always be .true. when convergence
c                   is tested.
c
c     nsamea        is the number of consecutive times that the
c                   left-hand end point of the interval of uncertainty
c                   has remained the same.
c
c     nsameb        similarly for the right-hand end.
c
c     a, b, alfbst  define the current interval of uncertainty.
c                   a minimizer lies somewhere in the interval
c                   [alfbst + a, alfbst + b].
c
c     alfbst        is the best point so far.  it is always at one end
c                   of the interval of uncertainty.  hence we have
c                   either  a lt 0,  b = 0  or  a = 0,  b gt 0.
c
c     fbest, gbest  are the values of f, f' at the point alfbst.
c
c     factor        controls the rate at which extrapolated estimates
c                   of alfa may expand into the interval of uncertainty.
c                   factor is not used if a minimizer has been bracketed
c                   (i.e., when the variable braktd is .true.).
c
c     fw, gw        are the values of f, f' at the point alfbst + xw.
c                   they are not defined until wset is .true..
c
c     xtry          is the trial point in the shifted interval (a, b).
c
c     xw            is such that  alfbst + xw  is the second-best point.
c                   it is not defined until  wset  is .true..
c                   in some cases,  xw  will replace a previous  xw
c                   that has a lower function but has just been excluded
c                   from the interval of uncertainty.

      double precision zero, point1, half
      parameter (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision one, three, five
      parameter (one=1.0d+0,three=3.0d+0,five=5.0d+0)
      double precision ten, eleven
      parameter (ten=1.0d+1,eleven=1.1d+1)

      double precision alfa, alfbst, alfmax, epsaf, fbest, ftry, g0,
     *                  gbest, gtry, targtg, tolabs, tolrel, toltny
      integer inform, maxf, nout, numf
      logical           debug, done, first, imprvd

      double precision a, absr, artifa, artifb, b, daux, dtry, factor,
     *                  fw, gw, q, r, s, scale, tol, tolmax, truea,
     *                  trueb, xmidpt, xtry, xw
      integer nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, fitok,
     *                  found, moved, quitf, quiti, setxw, wset

      character*120     rec(7)

c     .. intrinsic functions ..
      intrinsic         abs, sqrt
     *
      save              braktd, crampd, extrap, moved, wset, nsamea,
     *                  nsameb, a, b, factor, xtry, xw, fw, gw, tolmax

c

c     local variables
c     ===============
c
c     closef     is .true. if the new function ftry is within epsaf of
c                fbest (up or down).
c
c     found      is .true. if the sufficient decrease conditions hold at
c                alfbst.
c
c     quitf      is .true. when  maxf  function calls have been made.
c
c     quiti      is .true. when the interval of uncertainty is less than
c                2*tol.

c
      badfun = .false.
      quitf = .false.
      quiti = .false.
      imprvd = .false.
c
      if (first) then

c        first entry.  initialize various quantities, check input data
c        and prepare to evaluate the function at the initial alfa.

         first = .false.
         numf = 0
         alfbst = zero
         badfun = alfmax .le. toltny.or.g0 .ge. zero
         done = badfun
         moved = .false.
c
         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            wset = .false.
            nsamea = 0
            nsameb = 0
c
            tolmax = tolabs + tolrel*alfmax
            a = zero
            b = alfmax + tolmax
            factor = five
            tol = tolabs
            xtry = alfa
            if (debug) then
               write (rec,fmt=99999) g0, tolabs, alfmax, targtg, tolrel,
     *           epsaf, crampd
               call x04bay(nout,4,rec)
            end if
         end if
      else

c        subsequent entries. the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry and gtry.

         if (debug) then
            write (rec,fmt=99998) alfa, ftry, gtry
            call x04bay(nout,2,rec)
         end if
c
         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1
c
         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if
c
c        see if the new step is better.  if alfa is large enough that
c        ftry can be distinguished numerically from zero,  the function
c        is required to be sufficiently negative.
c
         closef = abs(ftry-fbest) .le. epsaf
         if (closef) then
            imprvd = abs(gtry) .le. abs(gbest)
         else
            imprvd = ftry .lt. fbest
         end if
c
         if (imprvd) then
c
c           we seem to have an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.
c
            fw = fbest
            fbest = ftry
            gw = gbest
            gbest = gtry
            alfbst = alfa
            moved = .true.
c
            a = a - xtry
            b = b - xtry
            xw = zero - xtry
            wset = .true.
            extrap = xw .lt. zero.and.gbest .lt. zero.or.xw .gt.
     *               zero.and.gbest .gt. zero
c
c           decrease the length of the interval of uncertainty.
c
            if (gtry.le.zero) then
               a = zero
               nsamea = 0
            else
               b = zero
               nsameb = 0
               braktd = .true.
            end if
         else
c
c           the new function value is not better than the best point so
c           far.  the origin remains unchanged but the new point may
c           qualify as xw.  xtry must be a new bound on the best point.
c
            if (xtry.le.zero) then
               a = xtry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if
c
c           if xw has not been set or ftry is better than fw, update the
c           points accordingly.
c
            if (wset) then
               setxw = ftry .lt. fw.or..not. extrap
            else
               setxw = .true.
            end if
c
            if (setxw) then
               xw = xtry
               fw = ftry
               gw = gtry
               wset = .true.
               extrap = .false.
            end if
         end if
c

c        check the termination criteria.  wset will always be .true..

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b
c
         found = abs(gbest) .le. targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol
c
         if (quiti.and..not. moved) then
c
c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.
c
            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf.or.tol .le. toltny
         end if
c
         done = quitf.or.quiti.or.found
c
         if (debug) then
            write (rec,fmt=99997) truea, trueb, b - a, tol, nsamea,
     *        nsameb, numf, braktd, extrap, closef, imprvd, found,
     *        quiti, alfbst, fbest, gbest, alfbst + xw, fw, gw
            call x04bay(nout,7,rec)
         end if
c

c        proceed with the computation of a trial steplength.
c        the choices are...
c        1. parabolic fit using derivatives only, if the f values are
c           close.
c        2. cubic fit for a minimizer, using both f and f'.
c        3. damped cubic or parabolic fit if the regular fit appears to
c           be consistently overestimating the distance to a minimizer.
c        4. bisection, geometric bisection, or a step of  tol  if
c           choices 2 or 3 are unsatisfactory.

         if (.not. done) then
            xmidpt = half*(a+b)
            s = zero
            q = zero
c
            if (closef) then
c              ---------------------------------------------------------
c              fit a parabola to the two best gradient values.
c              ---------------------------------------------------------
               s = gbest
               q = gbest - gw
               if (debug) then
                  write (rec,fmt=99995)
                  call x04baf(nout,rec(1))
               end if
            else
c              ---------------------------------------------------------
c              fit cubic through  fbest  and  fw.
c              ---------------------------------------------------------
               if (debug) then
                  write (rec,fmt=99996)
                  call x04baf(nout,rec(1))
               end if
               fitok = .true.
               r = three*(fbest-fw)/xw + gbest + gw
               absr = dabs(r)
               s = dsqrt(dabs(gbest))*dsqrt(dabs(gw))
c
c              compute  q =  the square root of  r*r - gbest*gw.
c              the method avoids unnecessary underflow and overflow.
c
               if ((gw.lt.zero.and.gbest.gt.zero)
     *            .or.(gw.gt.zero.and.gbest.lt.zero)) then
                  scale = absr + s
                  if (scale.eq.zero) then
                     q = zero
                  else
                     q = scale*dsqrt((absr/scale)**2+(s/scale)**2)
                  end if
               else if (absr.ge.s) then
                  q = dsqrt(absr+s)*dsqrt(absr-s)
               else
                  fitok = .false.
               end if
c
               if (fitok) then
c
c                 compute a minimizer of the fitted cubic.
c
                  if (xw.lt.zero) q = -q
                  s = gbest - r - q
                  q = gbest - gw - q - q
               end if
            end if

c           construct an artificial interval  (artifa, artifb)  in which
c           the new estimate of a minimizer must lie.  set a default
c           value of xtry that will be used if the polynomial fit fails.

            artifa = a
            artifb = b
            if (.not. braktd) then
c
c              a minimizer has not been bracketed.  set an artificial
c              upper bound by expanding the interval  xw  by a suitable
c              factor.
c
               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = five*factor
c
            else if (extrap) then
c
c              the points are configured for an extrapolation.
c              set a default value of  xtry  in the interval  (a, b)
c              that will be used if the polynomial fit is rejected.  in
c              the following,  dtry  and  daux  denote the lengths of
c              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
c              appropriate).  the value of  xtry is the point at which
c              the exponents of  dtry  and  daux  are approximately
c              bisected.
c
               daux = abs(xw)
               dtry = b - a
               if (daux.ge.dtry) then
                  xtry = five*dtry*(point1+dtry/daux)/eleven
               else
                  xtry = half*dsqrt(daux)*dsqrt(dtry)
               end if
               if (xw.gt.zero) xtry = -xtry
               if (debug) then
                  write (rec,fmt=99993) xtry, daux, dtry
                  call x04baf(nout,rec(1))
               end if
c
c              reset the artificial bounds.  if the point computed by
c              extrapolation is rejected,  xtry will remain at the
c              relevant artificial bound.
c
               if (xtry.le.zero) artifa = xtry
               if (xtry.gt.zero) artifb = xtry
            else
c
c              the points are configured for an interpolation.  the
c              default value xtry bisects the interval of uncertainty.
c              the artificial interval is just (a, b).
c
               xtry = xmidpt
               if (debug) then
                  write (rec,fmt=99994) xtry
                  call x04baf(nout,rec(1))
               end if
               if (nsamea.ge.3.or.nsameb.ge.3) then
c
c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation.
c
                  factor = factor/five
                  s = factor*s
               else
                  factor = one
               end if
            end if

c           the polynomial fits give  (s/q)*xw  as the new step.
c           reject this step if it lies outside  (artifa, artifb).

            if (q.ne.zero) then
               if (q.lt.zero) s = -s
               if (q.lt.zero) q = -q
               if (s*xw.ge.q*artifa.and.s*xw.le.q*artifb) then
c
c                 accept the polynomial fit.
c
                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) then
                     write (rec,fmt=99992) xtry
                     call x04baf(nout,rec(1))
                  end if
               end if
            end if
         end if
      end if
c

c
      if (.not. done) then
         alfa = alfbst + xtry
         if (braktd.or.alfa.lt.alfmax-tolmax) then
c
c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)
c
            if (xtry.le.a+tol.or.xtry.ge.b-tol) then
               if (half*(a+b).le.zero) then
                  xtry = -tol
               else
                  xtry = tol
               end if
               alfa = alfbst + xtry
            end if
         else
c
c           the step is close to, or larger than alfmax, replace it by
c           alfmax to force evaluation of  f  at the boundary.
c
            braktd = .true.
            xtry = alfmax - alfbst
            alfa = alfmax
         end if
      end if
c

c     exit.

      if (done) then
         if (badfun) then
            inform = 8
         else if (found) then
            if (alfbst.lt.alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved) then
            inform = 3
         else if (quitf) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if
c
      if (debug) then
         write (rec,fmt=99991)
         call x04bay(nout,2,rec)
      end if

c     end of e04uck. (srchc)

99999 format (/'     g0  tolabs  alfmax        ',1p,2d22.14,d16.8,/' t',
     *       'argtg  tolrel   epsaf        ',1p,2d22.14,d16.8,/' cramp',
     *       'd                        ',l3)
99998 format (/' alfa    ftry    gtry          ',1p,2d22.14,d16.8)
99997 format (/' a       b       b - a   tol   ',1p,2d22.14,2d16.8,
     *       /' nsamea  nsameb  numf          ',3i3,/' braktd  extrap ',
     *       ' closef  imprvd',4l3,/' found   quiti                 ',
     *       2l3,/' alfbst  fbest   gbest         ',1p,3d22.14,/' alfa',
     *       'w   fw      gw            ',1p,3d22.14)
99996 format (' cubic.   ')
99995 format (' parabola.')
99994 format (' bisection.              xmidpt',1p,d22.14)
99993 format (' geo. bisection. xtry,daux,dtry',1p,3d22.14)
99992 format (' polynomial fit accepted.  xtry',1p,d22.14)
99991 format (' ----------------------------------------------------',/)
      end


      subroutine e04ucj(first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,
     *                  tolrel,toltny,alfa,alfbst,fbest)

c     e04ucj  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     e04ucj  requires  f(alpha) (but not f'(alpha)) to be evaluated
c     in the interval.  new estimates of a minimizer are computed using
c     safeguarded quadratic interpolation.
c
c     reverse communication is used to allow the calling program to
c     evaluate f.  some of the parameters must be set or tested by the
c     calling program.  the remainder would ordinarily be local
c     variables.
c
c     input parameters (relevant to the calling program)
c     --------------------------------------------------
c
c     first         must be .true. on the first entry.
c                   it is subsequently altered by e04ucj.
c
c     debug         specifies whether detailed output is wanted.
c
c     maxf          is an upper limit on the number of times e04ucj is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).
c
c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by e04ucj (see below).
c
c     alfmax        is the upper limit of the interval to be searched.
c
c     alfsml        is intended to prevent inefficiency when a minimizer
c                   is very small, for cases where the calling program
c                   would prefer to redefine f'(alfa).  alfsml is
c                   allowed to be zero.  early termination will occur if
c                   e04ucj determines that a minimizer lies somewhere in
c                   the interval [0, alfsml) (but not if alfmax is
c                   smaller that alfsml).
c
c     epsaf         is an estimate of the absolute precision in the
c                   computed value of f(0).
c
c     ftry          the value of f at the new point
c                   alfa = alfbst + xtry.
c
c     g0            is the value of f'(0).  g0 must be negative.
c
c     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
c                   such that if f has already been evaluated at alfa,
c                   it will not be evaluated closer than tol(alfa).
c                   these values may be reduced by e04uck.
c
c     targtg        is the target value of abs(f'(alfa)). the search
c                   is terminated when
c                    abs(f'(alfa)) le targtg and f(alfa) lt 0.
c
c     toltny        is the smallest value that tolabs is allowed to be
c                   reduced to.
c
c     output parameters (relevant to the calling program)
c     ---------------------------------------------------
c
c     imprvd        is .true. if the previous alfa was the best point so
c                   far.  any related quantities should be saved by the
c                   calling program (e.g., arrays) before paying
c                   attention to the variable done.
c
c     done = .false.  means the calling program should evaluate ftry
c                   for the new trial step alfa, and reenter e04ucj.
c
c     done = .true.   means that no new alfa was calculated.  the value
c                   of inform gives the result of the search as follows
c
c                   inform = 1 means the search has terminated
c                              successfully with alfbst < alfmax.
c
c                   inform = 2 means the search has terminated
c                              successfully with alfbst = alfmax.
c
c                   inform = 3 means that the search failed to find a
c                              point of sufficient decrease in maxf
c                              functions, but a lower point was found.
c
c                   inform = 4 means alfmax is so small that a search
c                              should not have been attempted.
c
c                   inform = 5 means that the search was terminated
c                              because of alfsml (see above).
c
c                   inform = 6 means the search has failed to find a
c                              useful step.  the interval of uncertainty
c                              is [0,b] with b < 2*tolabs. a minimizer
c                              lies very close to alfa = 0, or f'(0) is
c                              not sufficiently accurate.
c
c                   inform = 7 if no better point could be found after
c                              maxf  function calls.
c
c                   inform = 8 means the input parameters were bad.
c                              alfmax le toltny  or  g0 ge zero.
c                              no function evaluations were made.
c
c     numf          counts the number of times e04ucj has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).
c
c     alfa          is the point at which the next function ftry must
c                   be computed.
c
c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever e04ucj returns
c                   inform = 1, 2 or 3.
c
c     fbest         will be the corresponding value of f.
c
c     the following parameters retain information between entries
c     -----------------------------------------------------------
c
c     braktd        is .false. if f has not been evaluated at the far
c                   end of the interval of uncertainty.  in this case,
c                   the point b will be at alfmax + tol(alfmax).
c
c     crampd        is .true. if alfmax is very small (le tolabs).  if
c                   the search fails, this indicates that a zero step
c                   should be taken.
c
c     extrap        is .true. if alfbst has moved at least once and xv
c                   lies outside the interval of uncertainty.  in this
c                   case, extra safeguards are applied to allow for
c                   instability in the polynomial fit.
c
c     moved         is .true. if a better point has been found, i.e.,
c                   alfbst gt 0.
c
c     vset          records whether a third-best point has been defined.
c
c     wset          records whether a second-best point has been
c                   defined.  it will always be .true. by the time the
c                   convergence test is applied.
c
c     nsamea        is the number of consecutive times that the
c                   left-hand end point of the interval of uncertainty
c                   has remained the same.
c
c     nsameb        similarly for the right-hand end.
c
c     a, b, alfbst  define the current interval of uncertainty.
c                   a minimizer lies somewhere in the  interval
c                   [alfbst + a, alfbst + b].
c
c     alfbst        is the best point so far.  it lies strictly within
c                   [atrue,btrue]  (except when alfbst has not been
c                   moved, in which case it lies at the left-hand end
c                   point).  hence we have a .le. 0 and b .gt. 0.
c
c     fbest         is the value of f at the point alfbst.
c
c     fa            is the value of f at the point alfbst + a.
c
c     factor        controls the rate at which extrapolated estimates of
c                   alfa  may expand into the interval of uncertainty.
c                   factor is not used if a minimizer has been bracketed
c                   (i.e., when the variable braktd is .true.).
c
c     fv, fw        are the values of f at the points alfbst + xv  and
c                   alfbst + xw.  they are not defined until  vset  or
c                   wset  are .true..
c
c     xtry          is the trial point within the shifted interval
c                   (a, b).  the new trial function value must be
c                   computed at the point alfa = alfbst + xtry.
c
c     xv            is such that alfbst + xv is the third-best point.
c                   it is not defined until vset is .true..
c
c     xw            is such that alfbst + xw is the second-best point.
c                   it is not defined until wset is .true..  in some
c                   cases,  xw will replace a previous xw that has a
c                   lower function but has just been excluded from
c                   (a,b).

      double precision zero, point1, half
      parameter (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision one, two, five
      parameter (one=1.0d+0,two=2.0d+0,five=5.0d+0)
      double precision ten, eleven
      parameter (ten=1.0d+1,eleven=1.1d+1)

      double precision alfa, alfbst, alfmax, alfsml, epsaf, fbest,
     *                  ftry, g0, targtg, tolabs, tolrel, toltny
      integer inform, maxf, nout, numf
      logical           debug, done, first, imprvd

      double precision a, artifa, artifb, b, daux, dtry, endpnt, fa,
     *                  factor, fv, fw, gv, gw, q, s, tol, tolmax,
     *                  truea, trueb, xmidpt, xtry, xv, xw
      integer nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, found,
     *                  moved, quitf, quitfz, quiti, quits, setxv, vset,
     *                  wset, xinxw

      character*120     rec(7)

c     .. intrinsic functions ..
      intrinsic         abs, dsqrt
      save              braktd, crampd, extrap, moved, vset, wset,
     *                  nsamea, nsameb, a, b, fa, factor, xtry, xw, fw,
     *                  xv, fv, tolmax


c     local variables
c     ===============
c
c     closef     is .true. if the worst function fv is within epsaf of
c                fbest (up or down).
c
c     found      is .true. if the sufficient decrease conditions holds
c                at alfbst.
c
c     quitf      is .true. when  maxf  function calls have been made.
c
c     quitfz     is .true. when the three best function values are
c                within epsaf of each other, and the new point satisfies
c                fbest le ftry le fbest+epsaf.
c
c     quiti      is .true. when the interval of uncertainty is less than
c                2*tol.
c
c     quits      is .true. as soon as alfa is too small to be useful;
c                i.e., btrue le alfsml.
c
c     xinxw      is .true. if xtry is in (xw,0) or (0,xw).

c
      imprvd = .false.
      badfun = .false.
      quitf = .false.
      quitfz = .false.
      quits = .false.
      quiti = .false.
c
      if (first) then

c        first entry.  initialize various quantities, check input data
c        and prepare to evaluate the function at the initial step alfa.

         first = .false.
         numf = 0
         alfbst = zero
         badfun = alfmax .le. toltny.or.g0 .ge. zero
         done = badfun
         moved = .false.
c
         if (.not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            vset = .false.
            wset = .false.
            nsamea = 0
            nsameb = 0
c
            tolmax = tolrel*alfmax + tolabs
            a = zero
            b = alfmax + tolmax
            fa = zero
            factor = five
            tol = tolabs
            xtry = alfa
            if (debug) then
               write (rec,fmt=99999) g0, tolabs, alfmax, targtg, tolrel,
     *           epsaf, crampd
               call x04bay(nout,4,rec)
            end if
         end if
      else

c        subsequent entries.  the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry.

         if (debug) then
            write (rec,fmt=99998) alfa, ftry
            call x04bay(nout,2,rec)
         end if
c
         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1
c
         if (.not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if
c
c        check if xtry is in the interval (xw,0) or (0,xw).
c
         if (wset) then
            xinxw = zero .lt. xtry.and.xtry .le. xw.or.xw .le.
     *              xtry.and.xtry .lt. zero
         else
            xinxw = .false.
         end if
c
         imprvd = ftry .lt. fbest
         if (vset) then
            closef = abs(fbest-fv) .le. epsaf
         else
            closef = .false.
         end if
c
         if (imprvd) then
c
c           we seem to have an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.
c
            if (wset) then
               xv = xw - xtry
               fv = fw
               vset = .true.
            end if
c
            xw = zero - xtry
            fw = fbest
            wset = .true.
            fbest = ftry
            alfbst = alfa
            moved = .true.
c
            a = a - xtry
            b = b - xtry
            extrap = .not. xinxw
c
c           decrease the length of (a,b).
c
            if (xtry.ge.zero) then
               a = xw
               fa = fw
               nsamea = 0
            else
               b = xw
               nsameb = 0
               braktd = .true.
            end if
         else if (closef.and.ftry-fbest.lt.epsaf) then
c
c           quit if there has been no progress and ftry, fbest, fw
c           and fv are all within epsaf of each other.
c
            quitfz = .true.
         else
c
c           the new function value is no better than the current best
c           point.  xtry must an end point of the new (a,b).
c
            if (xtry.lt.zero) then
               a = xtry
               fa = ftry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if
c
c           the origin remains unchanged but xtry may qualify as xw.
c
            if (wset) then
               if (ftry.lt.fw) then
                  xv = xw
                  fv = fw
                  vset = .true.
c
                  xw = xtry
                  fw = ftry
                  if (moved) extrap = xinxw
               else if (moved) then
                  if (vset) then
                     setxv = ftry .lt. fv.or..not. extrap
                  else
                     setxv = .true.
                  end if
c
                  if (setxv) then
                     if (vset.and.xinxw) then
                        xw = xv
                        fw = fv
                     end if
                     xv = xtry
                     fv = ftry
                     vset = .true.
                  end if
               else
                  xw = xtry
                  fw = ftry
               end if
            else
               xw = xtry
               fw = ftry
               wset = .true.
            end if
         end if
c

c        check the termination criteria.

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b
c
         found = moved.and.abs(fa-fbest) .le. -a*targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol
         quits = trueb .le. alfsml
c
         if (quiti.and..not. moved) then
c
c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.
c
            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf.or.tol .le. toltny
         end if
c
         done = quitf.or.quitfz.or.quits.or.quiti.or.found
c
         if (debug) then
            write (rec,fmt=99997) truea, trueb, b - a, tol, nsamea,
     *        nsameb, numf, braktd, extrap, closef, imprvd, found,
     *        quiti, quitfz, quits, alfbst, fbest, alfbst + xw, fw
            call x04bay(nout,7,rec)
            if (vset) then
               write (rec,fmt=99996) alfbst + xv, fv
               call x04bay(nout,2,rec)
            end if
         end if
c

c        proceed with the computation of an estimate of a minimizer.
c        the choices are...
c        1. parabolic fit using function values only.
c        2. damped parabolic fit if the regular fit appears to be
c           consistently overestimating the distance to a minimizer.
c        3. bisection, geometric bisection, or a step of tol if the
c           parabolic fit is unsatisfactory.

         if (.not. done) then
            xmidpt = half*(a+b)
            s = zero
            q = zero
c
c           ============================================================
c           fit a parabola.
c           ============================================================
c           see if there are two or three points for the parabolic fit.
c
            gw = (fw-fbest)/xw
            if (vset.and.moved) then
c
c              three points available.  use fbest, fw and fv.
c
               gv = (fv-fbest)/xv
               s = gv - (xv/xw)*gw
               q = two*(gv-gw)
               if (debug) then
                  write (rec,fmt=99994)
                  call x04baf(nout,rec(1))
               end if
            else
c
c              only two points available.  use fbest, fw and g0.
c
               if (moved) then
                  s = g0 - two*gw
               else
                  s = g0
               end if
               q = two*(g0-gw)
               if (debug) then
                  write (rec,fmt=99995)
                  call x04baf(nout,rec(1))
               end if
            end if
c

c           construct an artificial interval (artifa, artifb) in which
c           the new estimate of the steplength must lie.  set a default
c           value of  xtry  that will be used if the polynomial fit is
c           rejected. in the following, the interval (a,b) is considered
c           the sum of two intervals of lengths  dtry  and  daux, with
c           common end point the best point (zero).  dtry is the length
c           of the interval into which the default xtry will be placed
c           and endpnt denotes its non-zero end point.  the magnitude of
c           xtry is computed so that the exponents of dtry and daux are
c           approximately bisected.

            artifa = a
            artifb = b
            if (.not. braktd) then
c
c              a minimizer has not yet been bracketed.
c              set an artificial upper bound by expanding the interval
c              xw  by a suitable factor.
c
               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = five*factor
            else if (vset.and.moved) then
c
c              three points exist in the interval of uncertainty.
c              check if the points are configured for an extrapolation
c              or an interpolation.
c
               if (extrap) then
c
c                 the points are configured for an extrapolation.
c
                  if (xw.lt.zero) endpnt = b
                  if (xw.gt.zero) endpnt = a
               else
c
c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation step.
c
                  if (nsamea.ge.3.or.nsameb.ge.3) then
                     factor = factor/five
                     s = factor*s
                  else
                     factor = one
                  end if
c
c                 the points are configured for an interpolation.  the
c                 artificial interval will be just (a,b).  set endpnt so
c                 that xtry lies in the larger of the intervals (a,b)
c                 and  (0,b).
c
                  if (xmidpt.gt.zero) then
                     endpnt = b
                  else
                     endpnt = a
                  end if
c
c                 if a bound has remained the same for three iterations,
c                 set endpnt so that  xtry  is likely to replace the
c                 offending bound.
c
                  if (nsamea.ge.3) endpnt = a
                  if (nsameb.ge.3) endpnt = b
               end if
c
c              compute the default value of  xtry.
c
               dtry = abs(endpnt)
               daux = b - a - dtry
               if (daux.ge.dtry) then
                  xtry = five*dtry*(point1+dtry/daux)/eleven
               else
                  xtry = half*dsqrt(daux)*dsqrt(dtry)
               end if
               if (endpnt.lt.zero) xtry = -xtry
               if (debug) then
                  write (rec,fmt=99992) xtry, daux, dtry
                  call x04baf(nout,rec(1))
               end if
c
c              if the points are configured for an extrapolation set the
c              artificial bounds so that the artificial interval lies
c              within (a,b).  if the polynomial fit is rejected,  xtry
c              will remain at the relevant artificial bound.
c
               if (extrap) then
                  if (xtry.le.zero) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else
c
c              the gradient at the origin is being used for the
c              polynomial fit.  set the default xtry to one tenth xw.
c
               if (extrap) then
                  xtry = -xw
               else
                  xtry = xw/ten
               end if
               if (debug) then
                  write (rec,fmt=99993) xtry
                  call x04baf(nout,rec(1))
               end if
            end if
c

c           the polynomial fits give (s/q)*xw as the new step.  reject
c           this step if it lies outside (artifa, artifb).

            if (q.ne.zero) then
               if (q.lt.zero) s = -s
               if (q.lt.zero) q = -q
               if (s*xw.ge.q*artifa.and.s*xw.le.q*artifb) then
c
c                 accept the polynomial fit.
c
                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (debug) then
                     write (rec,fmt=99991) xtry
                     call x04baf(nout,rec(1))
                  end if
               end if
            end if
         end if
      end if

c
      if (.not. done) then
         alfa = alfbst + xtry
         if (braktd.or.alfa.lt.alfmax-tolmax) then
c
c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)
c
            xmidpt = half*(a+b)
            if (xtry.le.a+tol.or.xtry.ge.b-tol) then
               if (xmidpt.le.zero) then
                  xtry = -tol
               else
                  xtry = tol
               end if
            end if
c
            if (abs(xtry).lt.tol) then
               if (xmidpt.le.zero) then
                  xtry = -tol
               else
                  xtry = tol
               end if
            end if
            alfa = alfbst + xtry
         else
c
c           the step is close to or larger than alfmax, replace it by
c           alfmax to force evaluation of the function at the boundary.
c
            braktd = .true.
            xtry = alfmax - alfbst
            alfa = alfmax
         end if
      end if

c     exit.

      if (done) then
         if (badfun) then
            inform = 8
         else if (quits) then
            inform = 5
         else if (found) then
            if (alfbst.lt.alfmax) then
               inform = 1
            else
               inform = 2
            end if
         else if (moved) then
            inform = 3
         else if (quitf) then
            inform = 7
         else if (crampd) then
            inform = 4
         else
            inform = 6
         end if
      end if
c
      if (debug) then
         write (rec,fmt=99990)
         call x04bay(nout,2,rec)
      end if

c     end of  e04ucj. (srchq)
c
99999 format (/'     g0  tolabs  alfmax        ',1p,2d22.14,d16.8,/' t',
     *       'argtg  tolrel   epsaf        ',1p,2d22.14,d16.8,/' cramp',
     *       'd                        ',l3)
99998 format (/' alfa    ftry                  ',1p,2d22.14)
99997 format (/' a       b       b - a   tol   ',1p,2d22.14,2d16.8,
     *       /' nsamea  nsameb  numf          ',3i3,/' braktd  extrap ',
     *       ' closef  imprvd',4l3,/' found   quiti   quitfz  quits ',
     *       4l3,/' alfbst  fbest                 ',1p,2d22.14,/' alfa',
     *       'w   fw                    ',1p,2d22.14)
99996 format (' alfav   fv                    ',1p,2d22.14,/)
99995 format (' parabolic fit,    two points. ')
99994 format (' parabolic fit,  three points. ')
99993 format (' exponent reduced.  trial point',1p,d22.14)
99992 format (' geo. bisection. xtry,daux,dtry',1p,3d22.14)
99991 format (' polynomial fit accepted.  xtry',1p,d22.14)
99990 format (' ----------------------------------------------------',/)
      end


      subroutine e04nbu(n,nu,nrank,ldr,i,j,r,u,c,s)

c     e04nbu  interchanges the  i-th  and  j-th  (i .lt. j)  columns of
c     an  nrank*n  upper-trapezoidal matrix  r   and restores the
c     resulting matrix to upper-trapezoidal form using two sweeps of
c     plane rotations applied on the left.  r is overwritten.
c
c     if nu .gt. 0,  the rotations are applied to the  nu  columns of
c     the matrix  u.
      implicit none

      double precision zero
      parameter (zero=0.0d+0)

      integer i, j, ldr, n, nrank, nu

      double precision c(n), r(ldr,*), s(n), u(n,*)

      integer lenj

c     swap the elements of the i-th and j-th columns of r on, or above,
c     the main diagonal.

      call dswap(min(i,nrank),r(1,i),1,r(1,j),1)
      lenj = min(j,nrank)
c
      if (lenj.gt.i) then

c        reduce elements  r(i+1,j), ..., r(lenj,j)  to  beta*e(lenj)
c        using a backward sweep in planes
c        (lenj-1,lenj), (lenj-2,lenj), ..., (i+1,lenj).
c        if required, apply the sequence of rotations to u.

         call ssrotg ('fixed','backwards',lenj-i-1,r(lenj,j),r(i+1,j),1,
     *               c(i+1),s(i+1))

         if (nu.gt.0) call sgesrc ('left','bottom','backwards',n,nu,i+1,
     *                            lenj,c,s,u,n)

c        put zeros into the j-th column of r in positions corresponding
c        to the sub-diagonals of the i-th column.

         s(i) = r(lenj,j)
         call sload (lenj-i,zero,r(i+1,j),1)

c        apply the sequence of rotations to r.  this generates a spike
c        in the lenj-th row of r, which is stored in s.

         call sutsrs ('left',n,i+1,lenj,c,s,r,ldr)

c        eliminate the spike using a forward sweep in planes
c        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
c        if necessary, apply the sequence of rotations to u.

         call susqr('left',n,i,lenj,c,s,r,ldr)

         if (nu.gt.0) call sgesrc ('left','bottom','forwards',lenj,nu,i,
     *                            lenj,c,s,u,n)
      end if

      end

      subroutine e04nct(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,res,r,t,gq,
     *                  zy,c,s)
c     e04nct  updates the least-squares factor r and the factorization
c     a(free) (z y) = (0 t) when a regular, temporary or artificial
c     constraint is deleted from the working set.
      implicit none

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      integer jdel, kdel, lda, ldr, ldt, ldzy, n, nactiv,
     *                  nfree, ngq, nrank, nres, nrz, nz
      logical           unitq

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), zy(ldzy,*)
      integer kactiv(n), kx(n)

      double precision asize, dtmax, dtmin
      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision cs, sn
      integer i, ir, itdel, jart, k, ka, ld, npiv, nrz1, nsup,
     *                  nt

      character*80      rec(4)

      integer idamx1 
      external          idamx1 
c     .. intrinsic functions ..
      intrinsic         max, min

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin

      if (jdel.gt.0) then

c        regular constraint or temporary bound deleted.

c
         if (jdel.le.n) then
c
c           case 1.  a simple bound has been deleted.
c           =======  columns nfree+1 and ir of r must be swapped.
c
            ir = nz + kdel
            if (lsdbg.and.ilsdbg(1).gt.0) then
               write (rec,fmt=99998) nactiv, nz, nfree, ir, jdel, unitq
               call x04bay(iprint,4,rec)
            end if
c
            itdel = 1
            nfree = nfree + 1
c
            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,nfree,ir,r,
     *                                     res,c,s)
               call dswap(ngq,gq(nfree,1),n,gq(ir,1),n)
            end if
c
            if (.not. unitq) then
c
c              copy the incoming column of  a(free)  into the end of t.
c
               do 20 ka = 1, nactiv
                  i = kactiv(ka)
                  t(ka,nfree) = a(i,jdel)
   20          continue
c
c              expand q by adding a unit row and column.
c
               if (nfree.gt.1) then
                  call sload (nfree-1,zero,zy(nfree,1),ldzy)
                  call sload (nfree-1,zero,zy(1,nfree),1)
               end if
               zy(nfree,nfree) = one
            end if
         else
c
c           case 2.  a general constraint has been deleted.
c           =======
c
            if (lsdbg.and.ilsdbg(1).gt.0) then
               write (rec,fmt=99997) nactiv, nz, nfree, kdel, jdel,
     *           unitq
               call x04bay(iprint,4,rec)
            end if
c
            itdel = kdel
            nactiv = nactiv - 1
c
c           delete row  kdel  of t and move up the ones below it.
c           t becomes reverse lower hessenberg.
c
            do 40 i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld = nfree - i
               call dcopy(i+1,t(i+1,ld),ldt,t(i,ld),ldt)
   40       continue
         end if
c
         nz = nz + 1
c
         if (nactiv.eq.0) then
            dtmax = one
            dtmin = one
         else

c           restore the nactiv x (nactiv+1) reverse-hessenberg matrix  t
c           to reverse-triangular form.  the last nactiv super-diagonal
c           elements are removed using a backward sweep of plane
c           rotations.  the rotation for the singleton in the first
c           column is generated separately.

            nsup = nactiv - itdel + 1
c
            if (nsup.gt.0) then
               npiv = nfree - itdel + 1
               if (nsup.gt.1) then
                  call dcopy(nsup-1,t(nactiv-1,nz+1),ldt-1,s(nz+1),1)
                  call f06qzz('remove',nactiv,1,nsup,c(nz+1),s(nz+1),
     *                        t(1,nz+1),ldt)
               end if
c
               call srotg1 (t(nactiv,nz+1),t(nactiv,nz),cs,sn)
               t(nactiv,nz) = zero
               s(nz) = -sn
               c(nz) = cs
c
               call sgesrc ('right','variable','backwards',nfree,nfree,
     *                     nz,npiv,c,s,zy,ldzy)
               call sgesrc ('left ','variable','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gq,n)
c
               nt = min(nrank,npiv)
c
               if (nt.lt.npiv.and.nt.gt.0) then
c
c                 r is upper trapezoidal, pretend r is (nt x n) and
c                 apply the rotations in columns  max(nt,nz)  thru npiv.
c
                  call sgesrc ('right','variable','backwards',nt,n,
     *                        max(nt,nz),npiv,c,s,r,ldr)
               end if
c
c              apply the column transformations to the triangular part
c              of r.  a subdiagonal element is generated that must be
c              eliminated by a row rotation before the next column
c              transformation can be applied.
c
               if (nz.lt.nt) call sutsqr ('right',nt,nz,nt,c,s,r,ldr)

c              apply the row rotations to the remaining rows of r.
c
               call sgesrc ('left','variable','backwards',nt,n-nt,nz,nt,
     *                     c,s,r(1,min(nt+1,n)),ldr)
c
               if (nres.gt.0) call sgesrc ('left','variable','backwards'
     *                                     ,nt,nres,nz,nt,c,s,res,n)
c
            end if
            call scond (nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
         end if
      end if
c
      nrz1 = nrz + 1
c
      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamx1 (nz-nrz1+1,gq(nrz1,1))
         else
            jart = -jdel
         end if
c
         if (lsdbg.and.ilsdbg(1).gt.0) then
            write (rec,fmt=99999) nz, nrz1, jart
            call x04bay(iprint,4,rec)
         end if
c
         if (jart.gt.nrz1) then
c
c           swap columns nrz1 and jart of r.
c
            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,zy(1,nrz1),1,zy(1,jart),1)
            end if
c
            call dswap(ngq,gq(nrz1,1),n,gq(jart,1),n)
            if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,nrz1,jart,r,
     *                                  res,c,s)
         end if
      end if
c
      nrz = nrz1

c     end of  e04nct (lsdel).

99999 format (/' //e04nct //  artificial constraint deleted.      ',
     *       /' //e04nct //      nz   nrz   jart                 ',
     *       /' //e04nct //  ',3i6)
99998 format (/' //e04nct //  simple bound deleted.               ',
     *       /' //e04nct //  nactiv    nz nfree    ir  jdel unitq',
     *       /' //e04nct //  ',5i6,l6)
99997 format (/' //e04nct //  general constraint deleted.         ',
     *       /' //e04nct //  nactiv    nz nfree  kdel  jdel unitq',
     *       /' //e04nct //  ',5i6,l6)
      end


      subroutine e04nck(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,
     *                  nrz,istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,
     *                  jtiny,jbigst,kbigst,trulam,a,anorms,gq,rlamda,t,
     *                  wtinf)

c     e04nck  first computes the lagrange multiplier estimates for the
c     given working set.  it then determines the values and indices of
c     certain significant multipliers.  in this process, the multipliers
c     for inequalities at their upper bounds are adjusted so that a
c     negative multiplier for an inequality constraint indicates non-
c     optimality.  all adjusted multipliers are scaled by the 2-norm
c     of the associated constraint row.  in the following, the term
c     minimum refers to the ordering of numbers on the real line,  and
c     not to their magnitude.
c
c     jsmlst  is the index of the minimum of the set of adjusted
c             multipliers with values less than  - dinky.  a negative
c             jsmlst defines the index in q'g of the artificial
c             constraint to be deleted.
c     ksmlst  marks the position of general constraint jsmlst in kactiv.
c
c     jbigst  is the index of the largest of the set of adjusted
c             multipliers with values greater than (1 + dinky).
c     kbigst  marks its position in kactiv.
c
c     on exit,  elements 1 thru nactiv of rlamda contain the unadjusted
c     multipliers for the general constraints.  elements nactiv onwards
c     of rlamda contain the unadjusted multipliers for the bounds.

      implicit none

      double precision one
      parameter (one=1.0d+0)

      double precision dinky, trulam
      integer jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, msglvl, n, nactiv, nfree, nrz, numinf,
     *                  nz
      character*2       prbtyp

      double precision a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer istate(*), kactiv(n), kx(n)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision anormj, biggst, blam, rlam, scdlam, smllst,
     *                  tinylm
      integer i, is, j, k, l, nfixed

      character*80      rec(80)

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

      nfixed = n - nfree
c
      jsmlst = 0
      ksmlst = 0
      smllst = -dinky
c
      tinylm = dinky
      jtiny = 0
c
      jbigst = 0
      kbigst = 0
      biggst = one + dinky
c
      if (nrz.lt.nz) then

c        compute jsmlst for the artificial constraints.

         do 20 j = nrz + 1, nz
            rlam = -abs(gq(j))
            if (rlam.lt.smllst) then
               smllst = rlam
               jsmlst = -j
            else if (rlam.lt.tinylm) then
               tinylm = rlam
               jtiny = j
            end if
   20    continue
c
         if (msglvl.ge.20) then
            if (isumm.ge.0) then
               write (rec,fmt=99999)
               call x04bay(isumm,2,rec)
               do 40 j = nrz + 1, nz, 4
                  write (rec,fmt=99993) (gq(k),k=j,min(j+3,nz))
                  call x04baf(isumm,rec(1))
   40          continue
            end if
         end if
      end if
c     ---------------------------------------------------------------
c     compute jsmlst for regular constraints and temporary bounds.
c     ---------------------------------------------------------------
c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.
c
      if (n.gt.nz) call dcopy(n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call e04nbt(2,ldt,nactiv,t(1,nz+1),rlamda)
c     --------------------------------------------------------------
c     now set elements nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.
c     --------------------------------------------------------------
      do 80 l = 1, nfixed
         j = kx(nfree+l)
         blam = rlamda(nactiv+l)
         do 60 k = 1, nactiv
            i = kactiv(k)
            blam = blam - a(i,j)*rlamda(k)
   60    continue
         rlamda(nactiv+l) = blam
   80 continue
c
c     --------------------------------------------------------------
c     find jsmlst and ksmlst.
c     --------------------------------------------------------------
      do 100 k = 1, n - nz
         if (k.gt.nactiv) then
            j = kx(nz+k)
         else
            j = kactiv(k) + n
         end if
c
         is = istate(j)
c
         i = j - n
         if (j.le.n) anormj = one
         if (j.gt.n) anormj = anorms(i)
c
         rlam = rlamda(k)
c
c        change the sign of the estimate if the constraint is in
c        the working set at its upper bound.
c
         if (is.eq.2) rlam = -rlam
         if (is.eq.3) rlam = abs(rlam)
         if (is.eq.4) rlam = -abs(rlam)
c
         if (is.ne.3) then
            scdlam = rlam*anormj
            if (scdlam.lt.smllst) then
               smllst = scdlam
               jsmlst = j
               ksmlst = k
            else if (scdlam.lt.tinylm) then
               tinylm = scdlam
               jtiny = j
            end if
         end if
c
         if (numinf.gt.0.and.j.gt.jinf) then
            scdlam = rlam/wtinf(j)
            if (scdlam.gt.biggst) then
               biggst = scdlam
               trulam = rlamda(k)
               jbigst = j
               kbigst = k
            end if
         end if
  100 continue
c
c     --------------------------------------------------------------
c     if required, print the multipliers.
c     --------------------------------------------------------------
      if (msglvl.ge.20) then
         if (isumm.ge.0) then
            if (nfixed.gt.0) then
               write (rec,fmt=99998) prbtyp
               call x04bay(isumm,2,rec)
               do 120 j = 1, nfixed, 4
                  write (rec,fmt=99992) (kx(nfree+k),rlamda(nactiv+k),
     *              k=j,min(j+3,nfixed))
                  call x04baf(isumm,rec(1))
  120          continue
            end if
            if (nactiv.gt.0) then
               write (rec,fmt=99997) prbtyp
               call x04bay(isumm,2,rec)
               do 140 j = 1, nactiv, 4
                  write (rec,fmt=99992) (kactiv(k),rlamda(k),k=j,
     *              min(j+3,nactiv))
                  call x04baf(isumm,rec(1))
  140          continue
            end if
         end if
      end if
c
      if (lsdbg.and.ilsdbg(1).gt.0) then
         write (rec,fmt=99996) jsmlst, smllst, ksmlst
         call x04bay(iprint,3,rec)
         write (rec,fmt=99995) jbigst, biggst, kbigst
         call x04bay(iprint,2,rec)
         write (rec,fmt=99994) jtiny, tinylm
         call x04bay(iprint,2,rec)
      end if

c     end of  e04nck. (lsmuls)

99999 format (/' multipliers for the artificial constraints        ')
99998 format (/' multipliers for the ',a2,' bound  constraints   ')
99997 format (/' multipliers for the ',a2,' linear constraints   ')
99996 format (/' //e04nck//  jsmlst     smllst     ksmlst (scaled) ',
     *       /' //e04nck//  ',i6,1p,d11.2,5x,i6)
99995 format (' //e04nck//  jbigst     biggst     kbigst (scaled) ',
     *       /' //e04nck//  ',i6,1p,d11.2,5x,i6)
99994 format (' //e04nck//   jtiny     tinylm                     ',
     *       /' //e04nck//  ',i6,1p,d11.2)
99993 format (4(5x,1p,d11.2))
99992 format (4(i5,1p,d11.2))
      end



      subroutine e04ucr(needfd,inform,n,nfun,ngrad,
     *                  needc,objfun,alfa,alfbnd,alfmax,alfsml,
     *                  dxnorm,epsrf,eta,gdx,grdalf,glf1,glf,objf,
     *                  objalf,qpcurv,xnorm,
     *                  cmul1,cmul,cs1,cs,dx,dlam,dslk,grad,gradu,qpmul,
     *                  rho,slk1,slk,x1,x,work,w,lenw,iuser,user)

c     e04ucr finds the steplength alfa that gives sufficient decrease in
c     the augmented lagrangian merit function.
c
c     on exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
c     with an associated merit function value  objalf  which is lower
c     than that at the base point. if  inform = 4, 5, 6, 7 or 8,  alfa
c     is zero and  objalf will be the merit value at the base point.

      implicit none

      double precision zero, half, one
      parameter (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision two
      parameter (two=2.0d+0)
      double precision tolg
      parameter (tolg=1.0d-1)
      double precision rmu
      parameter (rmu=1.0d-4)

      double precision alfa, alfbnd, alfmax, alfsml, dxnorm, epsrf,
     *                  eta, gdx, glf, glf1, grdalf, objalf, objf,
     *                  qpcurv, xnorm
      integer inform, ldcj, lenw, n, nfun, ngrad
      logical           needfd

      double precision cmul(*), cmul1(*), cs(*),
     *                  cs1(*), dlam(*), dslk(*), dx(n), grad(n),
     *                  gradu(n), qpmul(*), rho(*), slk(*), slk1(*),
     *                  user(*), w(lenw), work(*), x(n), x1(n)
      integer iuser(*), needc(*)

      external          objfun

      double precision epspt3, epspt5, epspt8, epspt9, rhodmp, rhomax,
     *                  rhonrm, scale
      integer iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  ncdiff, nfdiff, nout
      logical           incrun, npdbg

      double precision wmach
      integer inpdbg(5)

      double precision alfbst, cs1jdx, csjdx, curvc, curvlf, epsaf,
     *                  epsmch, fbest, fterm, ftry, g0, gbest, gtry,
     *                  oldf, oldg, q, rhobfs, s, t, targtg, tgdx, tglf,
     *                  tobj, tobjm, tolabs, tolax, tolrel, tolrx,
     *                  toltny
      integer j, maxf, mode, nstate, numf
      logical           debug, done, first, imprvd

      character*80      rec(3)

      double precision ddot1
      external          ddot1

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common/ cstmch /wmach(9)
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /npdebg/inpdbg, npdbg

      epsmch = wmach(3)
c

c     set the input parameters and tolerances for e04uck and e04ucj.
c
c     tolrx   is the tolerance on relative changes in dx resulting from
c             changes in alfa.
c
c     tolax   is the tolerance on absolute changes in dx resulting from
c             changes in alfa.
c
c     tolabs  is the tolerance on absolute changes in alfa.
c
c     tolrel  is the tolerance on relative changes in alfa.
c
c     toltny  is the magnitude of the smallest allowable value of alfa.
c             if  m(tolabs) - m(0) .gt. epsaf,  the linesearch tries
c             steps in the range  toltny .le. alfa .le. tolabs.

      nstate = 0
      debug = npdbg.and.inpdbg(4) .gt. 0
c
      if (needfd) then
         maxf = 15
      else
         maxf = 10
      end if
c
      epsaf = epsrf*(one+abs(objalf))
      tolax = epspt8
      tolrx = epspt8
c
      if (tolrx*xnorm+tolax.lt.dxnorm*alfmax) then
         tolabs = (tolrx*xnorm+tolax)/dxnorm
      else
         tolabs = alfmax
      end if
      tolrel = max(tolrx,epsmch)
c
      t = zero
      do 20 j = 1, n
         s = abs(dx(j))
         q = abs(x(j))*tolrx + tolax
         if (s.gt.t*q) t = s/q
   20 continue
c
      if (t*tolabs.gt.one) then
         toltny = one/t
      else
         toltny = tolabs
      end if
c
      oldf = objalf
      oldg = grdalf
      alfbst = zero
      fbest = zero
      gbest = (one-rmu)*oldg
      targtg = (rmu-eta)*oldg
      g0 = gbest

      if (needfd) then
         mode = 0
      else
         mode = 2
      end if
c
      first = .true.
c

c     commence main loop, entering e04uck or e04ucj two or more times.
c     first = .true. for the first entry, .false. for subsequent entries
c     done  = .true. indicates termination, in which case the value of
c     inform gives the result of the search.
c     inform = 1 if the search is successful and alfa < alfmax.
c            = 2 if the search is successful and alfa = alfmax.
c            = 3 if a better point was found but too many functions
c                were needed (not sufficient decrease).
c            = 4 if alfmax < tolabs (too small to do a search).
c            = 5 if alfa < alfsml (e04ucj only -- maybe want to switch
c                to central differences to get a better direction).
c            = 6 if the search found that there is no useful step.
c                the interval of uncertainty is less than 2*tolabs.
c                the minimizer is very close to alfa = zero
c                or the gradients are not sufficiently accurate.
c            = 7 if there were too many function calls.
c            = 8 if the input parameters were bad
c                (alfmax le toltny  or  oldg ge 0).

c     +    repeat
   40 if (needfd) then
         iuser(2) = 0
         call e04ucj(first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest)
      else
         call e04uck(first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest,gbest)
      end if
c
      if (imprvd) then
         objf = tobj
         objalf = tobjm
c
         if (.not. needfd) then
            call dcopy(n,gradu,1,grad,1)
            gdx = tgdx
            glf = tglf

         end if
      end if
c
c     ---------------------------------------------------------------
c     if done = .false.,  the problem functions must be computed for
c     the next entry to e04uck or e04ucj.
c     if done = .true.,   this is the last time through.
c     ---------------------------------------------------------------
      if (.not. done) then
c
         call dcopy (n,x1,1,x,1)
         call daxpy (n,alfa,dx,1,x,1)
c        ------------------------------------------------------------
c        compute the value and gradient of the objective function.
c        ------------------------------------------------------------
c        if (done) then 
            iuser(2) = 1
c         else 
c           iuser(2) = 0
c         end if 
         call objfun(mode,n,x,tobj,gradu,nstate,iuser,user)

            tobjm = tobj

         ftry = tobjm - oldf - rmu*oldg*alfa
c
c        ------------------------------------------------------------
c        compute auxiliary gradient information.
c        ------------------------------------------------------------
         if (.not. needfd) then
            gtry = ddot1 (n,gradu,1,dx)
            tgdx = gtry
            tglf = gtry
            gtry = gtry - rmu*oldg
c
         end if
      end if
c     +    until (     done)
      if (.not. done) go to 40
c
      nfun = nfun + numf
      if (.not. needfd) ngrad = ngrad + numf
      alfa = alfbst
c
      if (.not. imprvd) then
         call dcopy (n,x1,1,x,1)
         call daxpy (n,alfa,dx,1,x,1)

      end if
c
      if (npdbg.and.inpdbg(1).gt.0) then
         write (rec,fmt=99999) inform
         call x04bay(iprint,2,rec)
      end if
c
      return
c
c     the user wants to stop.  who am i to object?
c
   60 inform = mode

c     end of  e04ucr. (npsrch)
c
99999 format (/' //e04ucr// inform  = ',i4)
99998 format (/' //e04ucr//        alfbnd          alfa        alfmax',
     *       /' //e04ucr//',1p,3d14.2)
99997 format (/' //e04ucr//         csjdx        cs1jdx        curvlf',
     *       /' //e04ucr//',1p,3d14.2)
      end


      subroutine e04udt(inform,n,nclin,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,anorm,adx,ax,bl,bu,dslk,dx,slk,x)

c     e04udt  finds a step alfa such that the point x + alfa*p reaches
c     one of the slacks or linear constraints.  the step alfa is the
c     maximum step that can be taken without violating one of the slacks
c     or linear constraints that is currently satisfied.

      implicit none

      integer lcmdbg
      parameter (lcmdbg=5)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision alfa, alfmax, alfmin, bigbnd, dxnorm
      integer inform, n, nclin

      double precision adx(*), anorm(*), ax(*), bl(*), bu(*), dslk(*),
     *                  dx(n), slk(*), x(n)

      double precision epspt3, epspt5, epspt8, epspt9
      integer iprint, isumm, lines1, lines2, nout
      logical           cmdbg

      integer icmdbg(lcmdbg)

      double precision adxi, axi, res, rownrm
      integer i, j

      character*80      rec(4)

c     .. intrinsic functions ..
      intrinsic         max

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg

c
      if (cmdbg.and.icmdbg(3).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,4,rec)
      end if
c
      alfa = alfmax
      j = 1
c
c     +    while (j .le. n+nclin.and.alfa .gt. alfmin) do
   20 if (j.le.n+nclin.and.alfa.gt.alfmin) then
c
         if (j.le.n) then
            axi = x(j)
            adxi = dx(j)
            rownrm = one
         else if (j.le.n+nclin) then
c
            i = j - n
            axi = ax(i)
            adxi = adx(i)
            rownrm = anorm(i) + one
         else
c
            i = j - n - nclin
            axi = slk(i)
            adxi = dslk(i)
            rownrm = one
         end if
c
         res = -one
         if (adxi.le.-epspt9*rownrm*dxnorm) then
c
c           constraint decreasing.
c
            adxi = -adxi
            if (bl(j).gt.-bigbnd) res = axi - bl(j)
         else if (adxi.gt.epspt9*rownrm*dxnorm) then
c
c           constraint increasing.
c
            if (bu(j).lt.bigbnd) res = bu(j) - axi
         end if
c
         if (res.gt.zero.and.alfa*adxi.gt.res) alfa = res/adxi
c
         if (cmdbg.and.icmdbg(3).gt.0) then
            write (rec,fmt=99998) j, res, adxi, alfa
            call x04baf(iprint,rec(1))
         end if
c
         j = j + 1
         go to 20
c        +    end while
      end if
c

c     determine alfa, the bound on the step to be taken.

      alfa = max(alfa,alfmin)
c
      inform = 0
      if (alfa.ge.alfmax) inform = 1
c
      if (cmdbg.and.icmdbg(1).gt.0.and.inform.gt.0) then
         write (rec,fmt=99997) alfa
         call x04bay(iprint,4,rec)
      end if

c     end of  e04udt. (npalf)

99999 format (/' e04udt entered',/'    j            res             ap',
     *       '           alfa ',/)
99998 format (i5,3g15.5)
99997 format (/' //e04udt//  no finite step.',/' //e04udt//           ',
     *       '  alfa',/' //e04udt//  ',g15.4)
      end


      subroutine e04uct(ktcond,convrg,lsumry,msgnp,msgqp,ldr,ldt,n,
     *                  nclin,nctotl,nactiv,linact,nlnact,nz,
     *                  nfree,majit0,majits,minits,istate,alfa,nfun,
     *                  condhz,condh,condt,objalf,objf,gfnorm,gznorm,
     *                  cvnorm,ax,r,t,violn,x,work)

c     e04uct  prints various levels of output for e04ucz and e04upz.
c
c           msg        cumulative result
c           ---        -----------------
c
c        le   0        no output.
c
c        eq   1        nothing now (but full output later).
c
c        eq   5        one terse line of output.
c
c        ge  10        same as 5 (but full output later).
c
c        ge  20        objective function,  x,  ax  and  c.
c
c        ge  30        diagonals of  t  and  r.
c
c     debug print is performed depending on the logical variable npdbg.
c     npdbg is set true when idbg major iterations have been performed.
c     at this point,  printing is done according to a string of binary
c     digits of the form clsvt (stored in the integer array inpdbg).
c
c     c  set 'on' gives detailed information from the checking routines.
c     l  set 'on' gives information from the linesearch.
c     s  set 'on' gives information from the maximum step routine e04udt
c     v  set 'on' gives various vectors in  e04ucz  and its auxiliaries.
c     t  set 'on' gives a trace of which routine was called and an
c                 indication of the progress of the run.
c     for example, `major debug level 11000' gives much output from
c                  the checking routines and the linesearch.

      implicit none 

      double precision alfa, condh, condhz, condt, cvnorm, gfnorm,
     *                  gznorm, objalf, objf
      integer ldr, ldt, linact, majit0, majits, minits, msgnp,
     *                  msgqp, n, nactiv, nclin, nctotl, nfree,
     *                  nfun, nlnact, nz
      logical           convrg
      character*5       lsumry

      double precision ax(*), r(ldr,*), t(ldt,*), violn(*),
     *                  work(n), x(n)
      integer istate(nctotl)
      logical           ktcond(2)

      double precision rhodmp, rhomax, rhonrm, scale
      integer iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg

      integer inpdbg(5)

      double precision cviols
      integer i, inct, j, k, mjr, mnr, ndf, neval
      logical           first, newset, nlncon, prthdr

      character*132     rec(4)

      double precision dnrm2
      external          dnrm2
      intrinsic         min, mod

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /fe04uc/inpdbg, npdbg

      if (msgnp.ge.20) then
         if (isumm.ge.0) then
            write (rec,fmt=99999) majits
            call x04bay(isumm,4,rec)
         end if
      end if
c
      if (msgnp.ge.5) then
c
         mjr = mod(majits,1000)
         mnr = mod(minits,1000)
         neval = mod(nfun,1000)
         ndf = mod(nz,1000)
         nlncon = .false.
         first = majits .eq. majit0
c

c        if necessary, print a header.
c        print a single line of information.

         if (isumm.ge.0) then
c           -----------------------------------
c           terse line for the monitoring file.
c           -----------------------------------
            newset = lines1 .ge. 50000
            prthdr = msgqp .gt. 0.or.first.or.msgnp .ge. 20 .or.
     *               newset
c
            if (prthdr) then

                  write (rec,fmt=99996)
                  call x04bay(isumm,3,rec)

               lines1 = 0
            end if
c

               write (rec,fmt=99995) mjr, mnr, alfa, neval, objalf,
     *           gznorm, ndf, n - nfree, linact, gfnorm, condh, condhz,
     *           condt, convrg, ktcond(1), ktcond(2), lsumry
               call x04baf(isumm,rec(1))
            end if
            lines1 = lines1 + 1

c
         if (iprint.ge.0.and.isumm.ne.iprint) then
c           ------------------------------
c           terse line for the print file.
c           ------------------------------
            newset = lines2 .ge. 50000
            prthdr = msgqp .gt. 0.or.first.or.newset
c
            if (prthdr) then

                  write (rec,fmt=99992)
                  call x04bay(iprint,3,rec)

               lines2 = 0
            end if
c
            if (nlncon) then
               write (rec,fmt=99993) mjr, mnr, alfa, objalf, cvnorm,
     *           gznorm, condhz, lsumry
               call x04baf(iprint,rec(1))
            else
               write (rec,fmt=99991) mjr, mnr, alfa, objalf, gznorm,
     *           condhz, lsumry
               call x04baf(iprint,rec(1))
            end if
            lines2 = lines2 + 1
         end if
c
         if (msgnp.ge.20) then
            if (isumm.ge.0) then
 
                  write (rec,fmt=99990) objf
                  call x04bay(isumm,2,rec)

c
c              ---------------------------------------------------------
c              print the constraint values.
c              ---------------------------------------------------------
               write (rec,fmt=99988)
               call x04bay(isumm,3,rec)
               write (rec,fmt=99987)
               call x04bay(isumm,2,rec)
               do 20 i = 1, n, 5
                  write (rec,fmt=99982) (x(j),istate(j),j=i,min(i+4,n))
                  call x04baf(isumm,rec(1))
   20          continue
               if (nclin.gt.0) then
                  write (rec,fmt=99986)
                  call x04bay(isumm,2,rec)
                  do 40 i = 1, nclin, 5
                     write (rec,fmt=99982) (ax(k),istate(n+k),k=i,
     *                 min(i+4,nclin))
                     call x04baf(isumm,rec(1))
   40             continue
               end if

c
               if (msgnp.ge.30) then
c                 ------------------------------------------------------
c                 print the diagonals of  t  and  r.
c                 ------------------------------------------------------
                  inct = ldt - 1
                  if (nactiv.gt.0) then
                     call dcopy(nactiv,t(nactiv,nz+1),inct,work,1)
                     write (rec,fmt=99984)
                     call x04bay(isumm,2,rec)
                     do 80 i = 1, nactiv, 5
                        write (rec,fmt=99981) (work(j),j=i,
     *                    min(i+4,nactiv))
                        call x04baf(isumm,rec(1))
   80                continue
                  end if
                  write (rec,fmt=99983)
                  call x04bay(isumm,2,rec)
                  do 100 i = 1, n, 5
                     write (rec,fmt=99981) (r(j,j),j=i,min(i+4,n))
                     call x04baf(isumm,rec(1))
  100             continue
               end if
            end if
         end if
      end if
c
      if (msgnp.ge.20) then
         if (isumm.ge.0) then
            write (rec,fmt=99980)
            call x04bay(isumm,3,rec)
         end if
      end if
c
      lsumry(1:2) = '  '
      lsumry(4:5) = '  '

c     end of e04uct. (npprt)

99999 format (//' major iteration',i5,/' ====================')
99998 format (//'  maj  mnr    step nfun  merit function  violtn norm ',
     *       'gz   nz  bnd  lin  nln penalty norm gf  cond h cond hz  ',
     *       'cond t conv')
99997 format (2i5,1p,d8.1,i5,d16.8,2d8.1,4i5,5d8.1,1x,l1,1x,2l1,a5)
99996 format (//'  maj  mnr    step nfun       objective norm gz   nz ',
     *       ' bnd  lin norm gf  cond h cond hz  cond t conv')
99995 format (2i5,1p,d8.1,i5,d16.8,d8.1,3i5,4d8.1,1x,l1,1x,2l1,a5)
99994 format (//'  maj  mnr    step merit function  violtn norm gz con',
     *       'd hz')
99993 format (2i5,1p,d8.1,d15.6,3d8.1,2x,a5)
99992 format (//'  maj  mnr    step      objective norm gz cond hz')
99991 format (2i5,1p,d8.1,d15.6,2d8.1,2x,a5)
99990 format (/' nonlinear objective value = ',1p,d15.6)
99989 format (/' nonlinear objective value = ',1p,d15.6,'   norm of th',
     *       'e nonlinear constraint violations = ',d15.6)
99988 format (/' values of the constraints and their predicted status',
     *       /' ----------------------------------------------------')
99987 format (/' variables                  ')
99986 format (/' general linear constraints ')
99985 format (/' nonlinear constraints      ')
99984 format (/' diagonals of  t  =         ')
99983 format (/' diagonals of  r  =         ')
99982 format (1x,5(1p,d15.6,i4))
99981 format (1p,5d15.6)
99980 format (' ======================================================',
     *       '===================================',//)
      end



      subroutine e04ucw(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,
     *                  nviol,ax,bl,bu,featol,x,work)

c     e04ucw  computes the following...
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.

      implicit none 

      double precision zero
      parameter (zero=0.0d+0)

      double precision bigbnd, cvnorm, errmax
      integer jmax, n, nclin, nviol

      double precision ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin),
     *                  work(n+nclin), x(n)
      integer istate(n+nclin)

      integer iprint, isumm, lines1, lines2, nout
      logical           npdbg

      integer inpdbg(5)

      double precision biglow, bigupp, con, feasj, res, tolj
      integer is, j

      character*80      rec(2)

      double precision dnrm2
      integer idamx1 
      external          dnrm2, idamx1 

c     .. intrinsic functions ..
      intrinsic         abs

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04uc/inpdbg, npdbg

c
      biglow = -bigbnd
      bigupp = bigbnd
c

c     compute nviol, the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint
c     violations and residuals of the constraints in the qp working set.

      nviol = 0
c
      do 40 j = 1, n + nclin
         feasj = featol(j)
         res = zero
c
         if (j.le.n+nclin) then
c
c           bound or general linear constraint.
c
            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if
c
            tolj = feasj

         end if
c
c        check for constraint violations.
c
         if (bl(j).gt.biglow) then
            res = bl(j) - con
            if (res.gt.feasj) nviol = nviol + 1
            if (res.gt.tolj) go to 20
         end if
c
         if (bu(j).lt.bigupp) then
            res = bu(j) - con
            if (res.lt.(-feasj)) nviol = nviol + 1
            if (res.lt.(-tolj)) go to 20
         end if
c
c        this constraint is satisfied,  but count the residual as a
c        violation if the constraint is in the working set.
c
         is = istate(j)
c
         if (is.eq.0) then
            res = zero
         else if (is.eq.1.or.is.le.-2) then
            res = bl(j) - con
         else if (is.ge.2.or.is.eq.-1) then
            res = bu(j) - con
         end if
c
         if (abs(res).gt.feasj) nviol = nviol + 1
c
c        set the array of violations.
c
   20    work(j) = res
   40 continue
c
      jmax = idamx1 (n+nclin,work)
      errmax = abs(work(jmax))
c
      if (npdbg.and.inpdbg(1).gt.0) then
         write (rec,fmt=99999) errmax, jmax
         call x04bay(iprint,2,rec)
      end if
c
      cvnorm = dnrm2(n+nclin,work,1)

c     end of  e04ucw. (npfeas)

99999 format (/' the maximum violation is ',1p,d14.2,' in ',
     *       'constraint',i5)
      end

      subroutine e04ucu(feasqp,unitq,nqperr,majits,minits,n,nclin,
     *                  ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,
     *                  numinf,istate,kactiv,kx,dxnorm,gdx,qpcurv,aqp,
     *                  adx,anorm,ax,bl,bu,clamda,cmul,cs,dlam,
     *                  dslk,dx,qpbl,qpbu,qptol,r,rho,slk,violn,x,wtinf,
     *                  w)

c     (1)  generate the upper and lower bounds for the qp  subproblem.
c
c     (2)  compute the  tq  factors of the rows of  aqp  specified by
c          the array  istate.  the part of the factorization defined by
c          the first contiguous group of linear constraints does not
c          need to be recomputed.  the remaining rows (which could be
c          comprised of both linear and nonlinear constraints) are
c          included as new rows of the  tq  factorization stored in
c          t and zy.  note that if there are no nonlinear constraints,
c          no factorization is required.
c
c     (3)  solve the  qp  subproblem.
c                 minimize     1/2 (w p - d)'(wp - d) + g'p
c
c                 subject to   qpbl .le. ( p) .le. qpbu,
c                                        (ap)
c
c          where  w  is a matrix (not stored) such that  w'w = h  and
c          wq = r,  d  is the zero vector,  and  g  is the gradient.
c          if the subproblem is infeasible, compute the point which
c          minimizes the sum of infeasibilities.
c
c     (4)   find the value of each slack variable for which the merit
c          function is minimized.
c
c     (5)   compute  dslk,  dlam  and  dx,  the search directions for
c          the slack variables, the multipliers and the variables.

      implicit none

      integer lenls
      parameter (lenls=20)

      integer mxparm
      parameter (mxparm=30)
      logical           qpnamd, vertex
      parameter (qpnamd=.false.,vertex=.false.)
      double precision zero, one, two
      parameter (zero=0.0d+0,one=1.0d+0,two=2.0d+0)
      double precision hundrd
      parameter (hundrd=1.0d+2)

      double precision dxnorm, gdx, qpcurv
      integer ldaqp, ldcj, ldr, linact, majits, minits, n,
     *                  nactiv, nclin, nfree, nlnact, nqperr,
     *                  numinf, nz
      logical           feasqp, unitq

      double precision adx(*), anorm(*), aqp(ldaqp,*), ax(*), bl(*),
     *                  bu(*), clamda(*), cmul(*),
     *                  cs(*), dlam(*), dslk(*), dx(n), qpbl(*),
     *                  qpbu(*), qptol(*), r(ldr,*), rho(*), slk(*),
     *                  violn(*), w(*), wtinf(*), x(n)
      integer istate(*), kactiv(n), kx(n)

      double precision asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldt,
     *                  ldzy, lformh, lines1, lines2, lprob,
     *                  lverfy, lvlder, msgls, msgnp, ncolt, nlnf, nlnj,
     *                  nlnx, nload, nn, nnclin, nncnln, nout, nprob,
     *                  nsave
      logical           cmdbg, incrun, lsdbg, npdbg

      double precision rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer icmdbg(5), ilsdbg(5), inpdbg(5),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), locls(lenls)

      double precision biglow, bigupp, blj, buj, con, condmx,
     *                  quotnt, ssq, ssq1, suminf, viol, weight, wscale,
     *                  wtmax, wtmin
      integer i, idbg, idbgsv, inform, iswap, j, jinf, k, k1,
     *                  k2, kviol, l, lgq, lhpq, lrlam, lrpq, lrpq0, lt,
     *                  lwrk1, lzy, mjrdbg, mnrdbg, msgqp, nartif, ncqp,
     *                  nctotl, ngq, nmajor, nminor, nplin, nrank,
     *                  nrejtd, nrpq, ntry, nviol, nz1
      logical           linobj

      double precision rprmls(mxparm), rprmnp(mxparm)
      integer iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(3)

      double precision ddot1, dnrm2, adivb 
      external          ddot1, dnrm2, adivb 
c     .. intrinsic functions ..
      intrinsic         abs, max, min

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common/ be04nb /ldt, ncolt, ldzy
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /fe04nb/icmdbg, cmdbg
      common            /fe04uc/inpdbg, npdbg
      common            /ge04uc/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /he04uc/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/

      idbgsv = idbg
      if (npdbg) then
         idbg = 0
      else
         idbg = nminor + 1
      end if
      lsdbg = npdbg
      cmdbg = npdbg
      call icopy(5,ilsdbg,1,icmdbg,1)
c
      lrpq = locls(5)
      lrpq0 = locls(6)
      lhpq = locls(8)
      lgq = locls(9)
      lrlam = locls(10)
      lt = locls(11)
      lzy = locls(12)
      lwrk1 = locls(14)
c
      nrpq = 0
      ngq = 1
c
      feasqp = .true.
      linobj = .true.
c
      biglow = -bigbnd
      bigupp = bigbnd
      ssq1 = zero
c
      nplin = n + nclin
      nctotl = nplin
      ncqp = nclin 
      nrank = n
      nrejtd = 0
c
      if (msgqp.gt.0) then
         write (rec,fmt=99999) majits
         call x04bay(iprint,3,rec)
      end if
c

c     generate the upper and lower bounds upon the search direction, the
c     weights on the sum of infeasibilities and the nonlinear constraint
c     violations.

      wscale = -one
      do 40 j = 1, nctotl
c
         if (j.le.n) then
            con = x(j)
         else if (j.le.nplin) then
            con = ax(j-n)
         end if
c
         blj = bl(j)
         buj = bu(j)
         if (blj.gt.biglow) blj = blj - con
         if (buj.lt.bigupp) buj = buj - con
c
         weight = one
         if (j.le.nplin) then
            if (abs(blj).le.qptol(j)) blj = zero
            if (abs(buj).le.qptol(j)) buj = zero
         else
            i = j - nplin
            viol = zero
            if (bl(j).gt.biglow) then
               if (blj.gt.zero) then
                  viol = blj
                  if (rho(i).gt.zero) then
                     weight = viol*rho(i)
                  else
                     weight = viol
                  end if
                  wscale = max(wscale,weight)
                  go to 20
               end if
            end if
c
            if (bu(j).lt.bigupp) then
               if (buj.lt.zero) then
                  viol = buj
                  if (rho(i).gt.zero) then
                     weight = -viol*rho(i)
                  else
                     weight = -viol
                  end if
                  wscale = max(wscale,weight)
               end if
            end if
c
c           set the vector of nonlinear constraint violations.
c
   20       violn(i) = viol
         end if
c
         wtinf(j) = weight
         qpbl(j) = blj
         qpbu(j) = buj
c
   40 continue
c
      if (wscale.gt.zero) then
         wscale = one/wscale
         call dscal(nctotl,(wscale),wtinf,1)
      end if
c
      call scond (nctotl,wtinf,1,wtmax,wtmin)
      wtmin = epspt9*wtmax
      do 60 j = 1, nctotl
         wtinf(j) = max(wtinf(j),wtmin)
   60 continue
c
c     set the maximum allowable condition estimator of the constraints
c     in the working set.  note that a relatively well-conditioned
c     working set is used to start the qp iterations.
c
      condmx = max(one/epspt3,hundrd)

c     solve for dx, the vector of minimum two-norm that satisfies the
c     constraints in the working set.

      call e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,ldaqp,
     *            ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,adx,qpbl,qpbu,
     *            w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt),w(lzy),w(lwrk1))
c

c     solve a quadratic program for the search direction  dx  and
c     multiplier estimates  clamda.

c     if there is no feasible point for the subproblem,  the sum of
c     infeasibilities is minimized subject to the linear constraints
c     (1  thru  jinf)  being satisfied.
c
      jinf = n + nclin
c
      ntry = 1
c     +    repeat
  180 call e04ncz('qp subproblem',qpnamd,names,linobj,unitq,nqperr,
     *            minits,jinf,ncqp,nctotl,nactiv,nfree,nrank,nz,nz1,n,
     *            ldaqp,ldr,istate,kactiv,kx,gdx,ssq,ssq1,suminf,numinf,
     *            dxnorm,qpbl,qpbu,aqp,clamda,adx,qptol,r,dx,w)
c
      if (npdbg.and.inpdbg(1).gt.0) then
         write (rec,fmt=99998) nqperr
         call x04bay(iprint,3,rec)
      end if
c
      nviol = 0
      if (numinf.gt.0) then
c
c           count the violated linear constraints.
c
         do 200 j = 1, nplin
            if (istate(j).lt.0) nviol = nviol + 1
  200    continue
c
         if (nviol.gt.0) then
            ntry = ntry + 1
            unitq = .true.
            nactiv = 0
            nfree = n
            nz = n
            call iload(nctotl,(0),istate,1)
c
            call e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,qpbl,qpbu,w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt)
     *                  ,w(lzy),w(lwrk1))
         end if
      end if
      if (.not. (nviol.eq.0.or.ntry.gt.2)) go to 180
c     +    until (   nviol .eq. 0 .or. ntry .gt. 2)
c

c     count the number of nonlinear constraint gradients in the  qp
c     working set.  make sure that all small  qp  multipliers associated
c     with nonlinear inequality constraints have the correct sign.

      linact = nactiv
c

c     extract various useful quantities from the qp solution.

c     compute  hpq = r'r(pq)  from the transformed gradient of the qp
c     objective function and  r(pq)  from the transformed residual.
c
      call dscal (n,(-one),w(lrpq),1)
      call daxpy (n,(-one),w(lgq),1,w(lhpq),1)
      qpcurv = two*ssq

      call icopy(5,inpdbg,1,icmdbg,1)
      idbg = idbgsv

c     end of  e04ucu. (npiqp)

99999 format (/1x,79('-'),/' start of major itn',i6)
99998 format (/' //e04ucu // nqperr',/' //e04ucu // ',i6)
99997 format (/' //e04ucu // dx recomputed with null space portion...')
99996 format (/' //e04ucu // violations = ')
99995 format (/' //e04ucu // slacks     = ')
99994 format (1p,5d15.6)
99993 format (5g12.3)
      end


      subroutine e04ucz(named,names,unitq,inform,majits,n,nclin,
     *                  nctotl,nactiv,nfree,nz,ldaqp,ldr,
     *                  nfun,ngrad,istate,kactiv,kx,objf,fdnorm,xnorm,
     *                  objfun,aqp,ax,bl,bu,clamda,
     *                  featol,grad,gradu,r,x,iw,w,lenw,iuser,user)

c     e04ucz  is the core routine for  nlpopt,  a sequential quadratic
c     programming (sqp) method for nonlinearly constrained optimization.

      implicit none 

      integer lenls
      parameter (lenls=20)
      integer lennp
      parameter (lennp=35)
      integer mxparm
      parameter (mxparm=30)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)
      double precision growth
      parameter (growth=1.0d+2)

      double precision fdnorm, objf, xnorm
      integer inform, ldaqp, ldcj, ldr, lenw, majits,
     *                  n, nactiv, nclin, nctotl, nfree, nfun,
     *                  ngrad, nz
      logical           named, unitq

      double precision aqp(ldaqp,*), ax(*), bl(nctotl), bu(nctotl),
     *                  clamda(nctotl), featol(nctotl), grad(n),
     *                  gradu(n), r(ldr,*), user(*), w(lenw), x(n)
      integer istate(*), iuser(*), iw(*), kactiv(n), kx(n)
      character*8       names(*)

      external          objfun

      double precision asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldt,
     *                  ldzy, lfdset, lformh, lines1, lines2,
     *                  lprob, lverfy, lvlder, lvldif, lvrfyc, msgls,
     *                  msgnp, ncdiff, ncolt, nfdiff, nlnf, nlnj, nlnx,
     *                  nload, nn, nnclin, nncnln, nout, nprob, nsave
      logical           cmdbg, incrun, npdbg

      double precision rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm), wmach
      integer icmdbg(5), inpdbg(5), ipadls(18),
     *                  ipadnp(12), ipsvls(mxparm), ipsvnp(mxparm),
     *                  jverfy(4), locls(lenls), locnp(lennp)

      double precision alfa, alfbnd, alfdx, alflim, alfmax, alfmin,
     *                  alfsml, cnorm, cond, condh, condhz, condt,
     *                  cvnorm, dinky, drzmax, drzmin, dxnorm, errmax,
     *                  flmax, gdx, gfnorm, glf1, glf2, glnorm, gltest,
     *                  grdalf, gtest, gznorm, obj, objalf, objsiz,
     *                  qpcurv, rootn, rtftol, rtmax, xsize
      integer idbg, info, jmax, ladx, lanorm, laqp, lbl, lbu,
     *                  lc1mul, lcjac1, lcjdx, lcjdx1, lcmul, lcs1,
     *                  lcs2, ldcj1, ldlam, ldslk, ldx, lgq, lgq1,
     *                  lhctrl, lhfrwd, lhpq, linact, liperm, lneedc,
     *                  lqptol, lqrwrk, lrho, lrlam, lrpq, lslk, lslk1,
     *                  lt, lvioln, lwrk1, lwrk2, lwrk3, lwtinf, lx1,
     *                  lzy, majit0, minits, mjrdbg, mnr, mnrdbg,
     *                  mnrsum, mode, msgqp, msgsv1, msgsv2, ncqp, nl,
     *                  nlnact, nlserr, nmajor, nminor, nplin, nqperr,
     *                  nqpinf, nstate, numinf, nviol
      logical           centrl, convpt, convrg, done, error, feasqp,
     *                  goodgq, infeas, needfd, newgq, optiml, overfl
      character*5       mjrmsg

      double precision rprmls(mxparm), rprmnp(mxparm)
      integer iprmls(mxparm), iprmnp(mxparm)
      logical           ktcond(2)
      character*80      rec(5)

      double precision ddot1, dnrm2, adivb 
      external          ddot1, dnrm2, adivb 

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common            /ae04uc/locnp
      common/ cstmch /wmach(9)
      common/ be04nb /ldt, ncolt, ldzy
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy
      common            /de04nb/asize, dtmax, dtmin
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /fe04nb/icmdbg, cmdbg
      common            /fe04uc/inpdbg, npdbg
      common            /ge04uc/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /he04uc/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)

      flmax = wmach(7)
      rtmax = wmach(8)

      lanorm = locls(2)
      lrpq = locls(5)
      lqrwrk = locls(6)
      lhpq = locls(8)
      lgq = locls(9)
      lrlam = locls(10)
      lt = locls(11)
      lzy = locls(12)
      lwtinf = locls(13)
      lwrk1 = locls(14)
      lqptol = locls(15)
c
      liperm = locnp(2)
      laqp = locnp(3)
      ladx = locnp(4)
      lbl = locnp(5)
      lbu = locnp(6)
      ldx = locnp(7)
      lgq1 = locnp(8)
      lx1 = locnp(11)
      lwrk2 = locnp(12)
      lcs1 = locnp(13)
      lcs2 = locnp(14)
      lc1mul = locnp(15)
      lcmul = locnp(16)
      lcjdx1 = locnp(17)
      ldlam = locnp(18)
      ldslk = locnp(19)
      lrho = locnp(20)
      lwrk3 = locnp(21)
      lslk1 = locnp(22)
      lslk = locnp(23)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)

      lcjac1 = laqp + nclin
      lcjdx = ladx + nclin
      lvioln = lwrk3

      mjrmsg = '     '
      nqpinf = 0
      mnrsum = 0

      majit0 = majits
      nplin = n + nclin
      ncqp = nclin
      nl = min(nplin+1,nctotl)

      ldcj1 = max(ncqp,1)

      needfd = lvlder .eq. 0.or.lvlder .eq. 2

      alfa = zero
      alfdx = zero
      rtftol = dsqrt(ftol)
      rootn = dsqrt(dble(n))
c
c     if debug printing is required,  turn off any extensive printing
c     until iteration  idbg.
c
      msgsv1 = msgnp
      msgsv2 = msgqp
      if (idbg.le.nmajor.and.idbg.gt.0) then
         msgnp = 0
         if (msgsv1.ge.5) msgnp = 5
         msgqp = 0
         if (msgsv2.ge.5) msgqp = 5
      end if
c     information from the feasibility phase will be used to generate a
c     hot start for the first qp subproblem.

      call dcopy(nctotl,featol,1,w(lqptol),1)

      nstate = 0
c
      objalf = objf
      newgq = .false.

   20 minits = 0
c
   40 centrl = lvldif .eq. 2
c
      if (newgq) then
         if (needfd) then
c           ------------------------------------------------------
c           compute any missing gradient elements and the
c           transformed gradient of the objective.
c           ------------------------------------------------------
            call evalfd (centrl,mode,n,bigbnd,cdint,
     *                  fdint,fdnorm,objf,objfun,iw(lneedc),bl,
     *                  bu,grad,gradu,
     *                  w(lhfrwd),w(lhctrl),x,w,lenw,iuser,user)
            inform = mode
            if (mode.lt.0) go to 60
c
         end if
c
         call dcopy(n,grad,1,w(lgq),1)
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),w(lwrk1))
         newgq = .false.
      end if
c
c     ============================================================
c     (1) solve an inequality quadratic program (iqp) for the
c         search direction and multiplier estimates.
c     (2) for each nonlinear inequality constraint,  compute
c         the slack variable for which the merit function is
c         minimized.
c     (3) compute the search direction for the slack variables
c         and multipliers.
c
c     note that the array violn is wrk3.
c     ============================================================
      call e04ucu(feasqp,unitq,nqperr,majits,mnr,n,nclin,
     *            ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,numinf,istate,
     *            kactiv,kx,dxnorm,gdx,qpcurv,aqp,w(ladx),w(lanorm),ax,
     *            bl,bu,clamda,w(lcmul),w(lcs1),w(ldlam),w(ldslk)
     *            ,w(ldx),w(lbl),w(lbu),w(lqptol),r,w(lrho),w(lslk),
     *            w(lvioln),x,w(lwtinf),w)
c
      minits = minits + mnr
      mnrsum = mnrsum + mnr
c
      if (feasqp) then
         nqpinf = 0
      else
         nqpinf = nqpinf + 1
         mjrmsg(2:2) = 'infeasible subproblem'
      end if
c
c     ============================================================
c     compute quantities needed for the convergence test.
c     ============================================================
c     compute the norms of the projected gradient and the
c     gradient with respect to the free variables.
c
      gznorm = zero
      if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)
      gfnorm = gznorm
      if (nfree.gt.0.and.nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),1)
c
c     if the forward-difference estimate of the transformed
c     gradient of the lagrangian function is small,  switch to
c     central differences, recompute the derivatives and re-solve
c     the qp.
c
      goodgq = .true.
      if (needfd.and..not. centrl) then
         glnorm = dnrm2(n,w(lhpq),1)

            cnorm = zero
c
         gltest = (one+abs(objf)+abs(cnorm))*epsrf/fdnorm
         if (glnorm.le.gltest) then
            goodgq = .false.
            mjrmsg(3:3) = 'central differences'
            lvldif = 2
            newgq = .true.
            if (msgnp.ge.5) then
               if (minits.gt.0) then
                  write (rec,fmt=99999) minits
                  call x04baf(iprint,rec(1))
               end if
            end if
         end if
c
      end if
c
c     +       until     (goodgq)
      if (.not. goodgq) go to 40
c
c     ===============================================================
c     (1) compute the number of constraints that are violated by more
c         than featol.
c     (2) compute the 2-norm of the residuals of the constraints in
c         the qp working set.
c     ===============================================================
      call e04ucw(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *            ax,bl,bu,featol,x,w(lwrk2))
c
c     define small quantities that reflect the magnitude of objf and
c     the norm of grad(free).
c
      objsiz = one + abs(objf)
      xsize = one + xnorm
      gtest = max(objsiz,gfnorm)
      dinky = rtftol*gtest
c
      if (nactiv.eq.0) then
         condt = zero
      else if (nactiv.eq.1) then
         condt = dtmin
      else
         condt = adivb (dtmax,dtmin,overfl)
      end if
c
      call scond (n,r,ldr+1,drmax,drmin)
c
      condh = adivb (drmax,drmin,overfl)
      if (condh.lt.rtmax) then
         condh = condh*condh
      else
         condh = flmax
      end if
c
      if (nz.eq.0) then
         condhz = one
      else if (nz.eq.n) then
         condhz = condh
      else
         call scond (nz,r,ldr+1,drzmax,drzmin)
         condhz = adivb (drzmax,drzmin,overfl)
         if (condhz.lt.rtmax) then
            condhz = condhz*condhz
         else
            condhz = flmax
         end if
      end if
c
c     ---------------------------------------------------------------
c     test for convergence.
c     the point test convpt checks for a k-t point at the initial
c     point or after a large change in x.
c     ---------------------------------------------------------------
      convpt = dxnorm .le. epspt8*gtest.and.nviol .eq. 0 .and.
     *         nqperr .le. 1
c
      ktcond(1) = gznorm .lt. dinky
      ktcond(2) = nviol .eq. 0
      optiml = ktcond(1).and.ktcond(2)
c
      convrg = majits .gt. 0.and.alfdx .le. rtftol*xsize

c     write (*,*) majits, alfdx,rtftol*xsize
c
      infeas = convrg.and..not. feasqp.or.nqpinf .gt. 7
c
      done = convpt.or.(convrg.and.optiml).or.infeas
c
      objalf = objf
      grdalf = gdx
      glf1 = gdx

c     print the details of this iteration.

      call e04uct(ktcond,convrg,mjrmsg,msgnp,msgqp,ldr,ldt,n,nclin,
     *            nctotl,nactiv,linact,nlnact,nz,nfree,majit0,
     *            majits,minits,istate,alfa,nfun,condhz,condh,condt,
     *            objalf,objf,gfnorm,gznorm,cvnorm,ax,r,w(lt),
     *            w(lvioln),x,w(lwrk1))
c
      alfa = zero
      error = majits .ge. nmajor
c
      if (.not. (done.or.error)) then
         majits = majits + 1
c
         if (majits.eq.idbg) then
            npdbg = .true.
            cmdbg = npdbg
            msgnp = msgsv1
            msgqp = msgsv2
         end if

c        make copies of information needed for the bfgs update.

         call dcopy(n,x,1,w(lx1),1)
         call dcopy(n,w(lgq),1,w(lgq1),1)

c        ============================================================
c        compute the parameters for the linesearch.
c        ============================================================
c        alfmin is the smallest allowable step predicted by the qp
c        subproblem.

         alfmin = one
         if (.not. feasqp) alfmin = zero

c        ------------------------------------------------------------
c        alfmax is the largest feasible steplength subject to a user-
c        defined limit alflim on the change in x.
c        ------------------------------------------------------------

            alfmax = adivb (bigdx,dxnorm,overfl)
            call e04udt(info,n,nclin,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,w(lanorm),w(ladx),ax,bl,bu,w(ldslk),
     *                  w(ldx),w(lslk),x)
            alfmax = alfa
            if (alfmax.lt.one+epspt3.and.feasqp) alfmax = one

c        ------------------------------------------------------------
c        alfbnd is a tentative upper bound on the steplength.  if the
c        merit function is decreasing at alfbnd and certain
c        conditions hold,  alfbnd will be increased in multiples of
c        two (subject to not being greater than alfmax).
c        ------------------------------------------------------------

            alfbnd = alfmax
c        ------------------------------------------------------------
c        alfsml trips the computation of central differences.  if a
c        trial steplength falls below alfsml, the linesearch is
c        terminated.
c        ------------------------------------------------------------
         alfsml = zero
         if (needfd.and..not. centrl) then
            alfsml = adivb (fdnorm,dxnorm,overfl)
            alfsml = min(alfsml,alfmax)
         end if

c        ============================================================
c        compute the steplength using safeguarded interpolation.
c        ============================================================
         alflim = adivb ((one+xnorm)*dxlim,dxnorm,overfl)
         alfa = min(alflim,one)

         call e04ucr(needfd,nlserr,n,nfun,ngrad,
     *               iw(lneedc),objfun,alfa,alfbnd,alfmax,alfsml,
     *               dxnorm,epsrf,eta,gdx,grdalf,glf1,glf2,objf,objalf,
     *               qpcurv,xnorm,
     *               w(lc1mul),w(lcmul),w(lcs1),w(lcs2),w(ldx),
     *               w(ldlam),w(ldslk),grad,gradu,clamda(nl),w(lrho),
     *               w(lslk1),w(lslk),w(lx1),x,w(lwrk2),w,lenw,iuser,
     *               user)

c           e04ucr  sets nlserr to the following values...

c           < 0  if the user wants to stop.
c             1  if the search is successful and alfa < alfmax.
c             2  if the search is successful and alfa = alfmax.
c             3  if a better point was found but too many functions
c                were needed (not sufficient decrease).

c           values of nlserr occurring with a nonzero value of alfa.
c             4  if alfmax < tolabs (too small to do a search).
c             5  if alfa  < alfsml (e04ucj only -- maybe want to switch
c                to central differences to get a better direction).
c             6  if the search found that there is no useful step.
c                the interval of uncertainty is less than 2*tolabs.
c                the minimizer is very close to alfa = zero
c                or the gradients are not sufficiently accurate.
c             7  if there were too many function calls.
c             8  if the input parameters were bad
c                (alfmax le toltny  or  uphill).

         if (nlserr.lt.0) then
            inform = nlserr
            go to 60
         end if

         if (alfa.gt.alflim) mjrmsg(4:4) = 'l'

         error = nlserr .ge. 4
         if (error) then
c           ---------------------------------------------------------
c           the linesearch failed to find a better point.
c           if exact gradients or central differences are being used,
c           or the kt conditions are satisfied, stop.  otherwise,
c           switch to central differences and solve the qp again.
c           ---------------------------------------------------------
            if (needfd.and..not. centrl) then
               if (.not. optiml) then
                  error = .false.
                  mjrmsg(3:3) = 'central differences'
                  lvldif = 2
                  newgq = .true.
                  if (msgnp.ge.5) then
                     if (minits.gt.0) then
                        write (rec,fmt=99999) minits
                        call x04baf(iprint,rec(1))
                     end if
                  end if
               end if
            end if
         else
            if (needfd) then
c              ======================================================
c              compute the missing gradients.
c              ======================================================
               mode = 1
               ngrad = ngrad + 1

               iuser(2) = 0

               call objfun(mode,n,x,obj,gradu,nstate,iuser,user)
c
               call dcopy(n,gradu,1,grad,1)
c
               call evalfd (centrl,mode,n,bigbnd,cdint,
     *                     fdint,fdnorm,objf,objfun,iw(lneedc),
     *                     bl,bu,grad,
     *                     gradu,w(lhfrwd),w(lhctrl),x,w,lenw,iuser,
     *                     user)
c
               inform = mode
               if (mode.lt.0) go to 60
c
               gdx = ddot1 (n,grad,1,w(ldx))
               glf2 = gdx

            end if
c
            call dcopy(n,grad,1,w(lgq),1)
            call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),
     *                  w(lwrk1))
c
            xnorm = dnrm2(n,x,1)

            if (nclin.gt.0) call daxpy (nclin,alfa,w(ladx),1,ax,1)
            alfdx = alfa*dxnorm
c
c           =========================================================
c           update the factors of the approximate hessian of the
c           lagrangian function.
c           =========================================================
            call e04ucl(mjrmsg,unitq,n,nfree,nz,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,
     *                  w(lcs1),w(lcs2),w(lgq1),
     *                  w(lgq),w(lhpq),w(lrpq),clamda(nl),r,w(lwrk3),
     *                  w(lzy),w(lwrk2),w(lwrk1))
 
            call scond (n,r,ldr+1,drmax,drmin)
            cond = adivb (drmax,drmin,overfl)
 
            if (cond.gt.rcndbd.or.rfrobn.gt.rootn*growth*drmax) then
c              ------------------------------------------------------
c              reset the condition estimator and range-space
c              partition of q'hq.
c              ------------------------------------------------------
               if (npdbg.and.inpdbg(1).gt.0) then
                  write (rec,fmt=99997) rfrobn, drmax, drmin, cond,
     *              rcndbd
                  call x04bay(iprint,5,rec)
               end if
c
               mjrmsg(5:5) = 'refactorize hessian'
c
               call e04udr(unitq,n,nfree,nz,ldzy,ldr,iw(liperm),kx,
     *                     w(lgq),r,w(lzy),w(lwrk1),w(lqrwrk))
            end if
         end if
      end if
c
c     +    until     (done .or. error)
      if (.not. (done.or.error)) go to 20
c
c     ======================end of main loop============================
c
      if (done) then
         if (convrg.and.optiml) then
            inform = 0
         else if (convpt) then
            inform = 1
         else if (infeas) then
            inform = 3
         end if
      else if (error) then
         if (majits.ge.nmajor) then
            inform = 4
         else if (optiml) then
            inform = 1
         else
            inform = 6
         end if
      end if


c     set  clamda.  print the full solution.

   60 msgnp = msgsv1
      msgqp = msgsv2
      if (msgnp.gt.0) then
         write (rec,fmt=99998) majits, mnrsum
         call x04bay(iprint,3,rec)
      end if
c
      call e04nbx(msgnp,nfree,ldaqp,n,nclin,nctotl,bigbnd,named,names,
     *            nactiv,istate,kactiv,kx,aqp,bl,bu,clamda,w(lrlam),x)


c     end of  e04ucz. (npcore)

99999 format (' mnr itn ',i4,' -- re-solve qp subproblem.')
99998 format (/' exit from np problem after ',i5,' major iterations,',
     *       /'                            ',i5,' minor iterations.')
99997 format (/' //e04ucz//        rfrobn         drmax         drmin',
     *       /' //e04ucz//',1p,3d14.2,/' //e04ucz//          cond     ',
     *       '   rcndbd',/' //e04ucz//',1p,2d14.2)
      end


      subroutine e04udr(unitq,n,nfree,nz,nq,nrowr,iperm,kx,gq,r,zy,work,
     *                  qrwork)

c     e04udr  bounds the condition estimator of the transformed hessian.
c     on exit, r is of the form
c                  (drz   0    )
c                  ( 0  sigma*i)
c     where d is a diagonal matrix such that drz has a bounded condition
c     number,  i is the identity matrix and sigma  is the geometric mean
c     of the largest and smallest elements of drz. the qr factorization
c     with interchanges is used to give diagonals of drz that are
c     decreasing in modulus.

      implicit none

      double precision zero, half, one
      parameter (zero=0.0d+0,half=0.5d+0,one=1.0d+0)

      integer n, nfree, nq, nrowr, nz
      logical           unitq

      double precision gq(n), qrwork(2*n), r(nrowr,*), work(n),
     *                  zy(nq,*)
      integer iperm(n), kx(n)

      double precision drmax, drmin, rcndbd, rfrobn
      logical           npdbg

      integer inpdbg(5)

      double precision drgm, drgs, gjmax, scle, sumsq
      integer info, j, jmax, jsave, nrank

      double precision snorm
      integer isrank
      external          snorm, isrank
c     .. intrinsic functions ..
      intrinsic         abs, dble, sqrt

      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /fe04uc/inpdbg, npdbg

c

c     bound the condition estimator of q'hq.

      if (nz.gt.1) then

c        refactorize rz.  interchanges are used to give diagonals
c        of decreasing magnitude.

         do 20 j = 1, nz - 1
            call sload (nz-j,zero,r(j+1,j),1)
   20    continue
c
         call sgeqrp('column iterchanges',nz,nz,r,nrowr,work,iperm,
     *               qrwork,info)
c
         do 40 j = 1, nz
            jmax = iperm(j)
            if (jmax.gt.j) then
               if (unitq) then
                  jsave = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j) = jsave
               else
                  call dswap(nfree,zy(1,jmax),1,zy(1,j),1)
               end if
c
               gjmax = gq(jmax)
               gq(jmax) = gq(j)
               gq(j) = gjmax
            end if
   40    continue
      end if
c
      drgm = one
c
      if (nz.gt.0) then
         nrank = isrank(nz,r,nrowr+1,one/rcndbd)
         drgm = half*dsqrt(dabs(r(1,1)*r(nrank,nrank)))
         drgs = dabs(r(1,1))/rcndbd
c
         if (nz.gt.nrank) then
            do 60 j = nrank + 1, nz
               call sload (j-1,zero,r(1,j),1)
   60       continue
            call sload (nz-nrank,drgs,r(nrank+1,nrank+1),nrowr+1)
         end if
      end if
c

c     reset the range-space partition of the hessian.

      if (nz.lt.n) then
         do 80 j = nz + 1, n
            call sload (j,zero,r(1,j),1)
   80    continue
         call sload (n-nz,drgm,r(nz+1,nz+1),nrowr+1)
      end if
c
c     recompute the frobenius norm of r.
c
      scle = dsqrt(dble(n-nz))*drgm
      sumsq = one
      do 100 j = 1, nz
         call sssq (j,r(1,j),1,scle,sumsq)
  100 continue

      rfrobn = snorm (scle,sumsq)

c     end of  e04udr. (nprset)
c
      end



      subroutine e04ncz(prbtyp,named,names,linobj,unitq,inform,iter,
     *                  jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nrz,n,
     *                  lda,ldr,istate,kactiv,kx,ctx,ssq,ssq1,suminf,
     *                  numinf,xnorm,bl,bu,a,clamda,ax,featol,r,x,w)

c     e04ncz  is a subroutine for linearly constrained linear-least
c     squares.  on entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.
c     the arrays istate, kactiv and kx will have been set accordingly
c     and the arrays t and zy will contain the tq factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  the tq factorization of the resulting (nactiv by nfree)
c     matrix is  a(free)*q = (0 t),  where q is (nfree by nfree) and t
c     is reverse-triangular.
c
c     values of istate(j) for the linear constraints.......
c
c     istate(j)
c     ---------
c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.
c
c     constraint j may be violated by as much as featol(j).

      implicit none 

      integer lenls
      parameter (lenls=20)

      integer mxparm
      parameter (mxparm=30)
      double precision zero, half, one
      parameter (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      integer mstall, mrefn
      parameter (mstall=50,mrefn=1)

      double precision ctx, ssq, ssq1, suminf, xnorm
      integer inform, iter, jinf, lda, ldr, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nrz, numinf, nz
      logical           linobj, named, unitq
      character*2       prbtyp

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  clamda(nctotl), featol(nctotl), r(ldr,*), w(*),
     *                  x(n)
      integer istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)

      double precision asize, bigbnd, bigdx, bndlow, bndupp, dtmax,
     *                  dtmin, epspt3, epspt5, epspt8, epspt9, tolact,
     *                  tolfea, tolrnk
      integer idbgls, iprint, iprnt, isumm, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, ldt, ldzy,
     *                  lines1, lines2, lprob, msgls, ncolt, nn, nnclin,
     *                  nout, nprob
      logical           cmdbg, lsdbg

      double precision rpadls(23), rpsvls(mxparm), wmach
      integer icmdbg(5), ilsdbg(5), ipadls(18),
     *                  ipsvls(mxparm), locls(lenls)

      double precision absrzz, alfa, alfhit, atphit, bigalf, cnorm,
     *                  condmx, condrz, condt, ctp, dinky, drzmax,
     *                  drzmin, err1, err2, flmax, gfnorm, grznrm,
     *                  gznorm, objsiz, palfa, pnorm, resnrm, rownrm,
     *                  trulam, wssize
      integer iadd, idbg, ifix, irefn, is, isdel, itmax, jadd,
     *                  jbigst, jdel, jmax1, jsmlst, jtiny, kbigst,
     *                  kdel, ksmlst, lanorm, lap, lcq, lgq, lhz, lpx,
     *                  lres, lres0, lrlam, lt, lwrk, lwtinf, lzy,
     *                  msgdbg, msglvl, msgsvd, ngq, nphase, nres,
     *                  nstall, nviol
      logical           convrg, cyclin, error, firstv, hitcon, hitlow,
     *                  needfg, overfl, prnt, rowerr, singlr, stall,
     *                  statpt, unbndd, uncon, unitgz, weak

      double precision rprmls(mxparm)
      integer iprmls(mxparm)
      character*80      rec(3)

      double precision dnrm2, adivb 
      external          dnrm2, adivb 
c     .. intrinsic functions ..
      intrinsic         abs, max

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common/ cstmch /wmach(9)
      common/ be04nb /ldt, ncolt, ldzy
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /fe04nb/icmdbg, cmdbg

      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (msgls,msglvl), (idbgls,idbg), (ldbgls,msgdbg)

      flmax = wmach(7)
      lanorm = locls(2)
      lap = locls(3)
      lpx = locls(4)
      lres = locls(5)
      lres0 = locls(6)
      lhz = locls(7)
      lgq = locls(8)
      lcq = locls(9)
      lrlam = locls(10)
      lt = locls(11)
      lzy = locls(12)
      lwtinf = locls(13)
      lwrk = locls(14)
c
c     set up the adresses of the contiguous arrays  (res0, res)
c     and  (gq, cq).
c
      nres = 0
      if (nrank.gt.0) nres = 2
      ngq = 1
      if (linobj) ngq = 2
c
c     initialize.
c
      irefn = 0
      iter = 0
      jadd = 0
      jdel = 0
      nphase = 1
      nstall = 0
      numinf = -1
      nrz = 0
c
      if (prbtyp.eq.'fp') then
         itmax = itmax2
      else
         itmax = itmax1
      end if
c
      alfa = zero
      condmx = flmax
      drzmax = one
      drzmin = one
      ssq = zero
c
      cyclin = .false.
      error = .false.
      firstv = .false.
      prnt = .true.
      needfg = .true.
      stall = .true.
      uncon = .false.
      unitgz = .true.
      unbndd = .false.
c
c     if debug output is required,  print nothing until iteration idbg.
c
      msgsvd = msglvl
      if (idbg.gt.0.and.idbg.le.itmax) then
         msglvl = 0
      end if
c
c     =================== start of the main loop =======================
c
c      cyclin = false
c      unbndd = false
c      error  = false
c      k      = 0
c
c      repeat
c            repeat
c                  compute z'g,  print details of this iteration
c                  stat pt = (z'g .eq. 0)
c                  if (not stat pt) then
c                     error =  k .ge. itmax
c                     if (not error) then
c                        compute p, alfa
c                        error = unbndd  or  cyclin
c                        if (not error) then
c                           k = k + 1
c                           x = x + alfa p
c                           if (feasible) update z'g
c                           if necessary, add a constraint
c                        end if
c                     end if
c                  end if
c            until  stat pt  or  error
c
c            compute lam1, lam2, smllst
c            optmul =  smllst .gt. 0
c            if (not (optmul.or.error)) then
c                  delete an artificial or regular constraint
c            end if
c      until optmul  or  error
c

c
c     repeat
c        repeat
   20 if (needfg) then
         if (nrank.gt.0) then
            resnrm = dnrm2(nrank,w(lres),1)
            ssq = half*(ssq1**2+resnrm**2)
         end if
c
         if (numinf.ne.0) then
c
c           compute the transformed gradient of either the sum of
c           infeasibilities or the objective.  initialize
c           singlr and unitgz.
c
            call e04ncp(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,w(lres),featol,
     *                  w(lgq),w(lcq),r,x,w(lwtinf),w(lzy),w(lwrk))
            if (numinf.eq.0.and.prbtyp.ne.'fp'.and.nphase.eq.1) then
               itmax = iter + itmax2
               nphase = 2
            end if
         end if
      end if
c
      gznorm = zero
      if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)
c
      if (nrz.eq.nz) then
         grznrm = gznorm
      else
         grznrm = zero
         if (nrz.gt.0) grznrm = dnrm2(nrz,w(lgq),1)
      end if
c
      gfnorm = gznorm
      if (nfree.gt.0.and.nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),1)
c
c        ------------------------------------------------------------
c        print the details of this iteration.
c        ------------------------------------------------------------
c        define small quantities that reflect the magnitude of  x,
c        r  and  the matrix of constraints in the working set.
c        use the largest and smallest diagonals of  r  to estimate
c        the condition number of  rz1.
c
      if (nrz.eq.0) then
         singlr = .false.
      else
         if (numinf.gt.0.or.nrz.gt.nrank) then
            absrzz = zero
            singlr = .true.
         else
            call scond (nrz,r,ldr+1,drzmax,drzmin)
            absrzz = abs(r(nrz,nrz))
            rownrm = dnrm2(n,r(1,1),ldr)
            singlr = absrzz .le. drzmax*tolrnk.or.rownrm .le.
     *               tolrnk.or.abs(r(1,1)) .le. rownrm*tolrnk
         end if
c
         if (lsdbg.and.ilsdbg(1).gt.0) then
            write (rec,fmt=99995) singlr, absrzz, drzmax, drzmin
            call x04bay(iprint,3,rec)
         end if
c
      end if
c
      condrz = adivb (drzmax,drzmin,overfl)
c
      condt = one
      if (nactiv.gt.0) condt = adivb (dtmax,dtmin,overfl)
c
      if (prnt) then
         call e04ncj(prbtyp,isdel,iter,jadd,jdel,msglvl,nactiv,nfree,n,
     *               nclin,nrank,ldr,ldt,nz,nrz,istate,alfa,condrz,
     *               condt,gfnorm,grznrm,numinf,suminf,ctx,ssq,ax,r,
     *               w(lt),x,w(lwrk))
c
         jdel = 0
         jadd = 0
         alfa = zero
      end if
c
      if (numinf.gt.0) then
         dinky = zero
      else
         objsiz = one + abs(ssq+ctx)
         wssize = zero
         if (nactiv.gt.0) wssize = dtmax
         dinky = epspt8*max(wssize,objsiz,gfnorm)
         if (uncon) then
            unitgz = grznrm .le. dinky
         end if
      end if
c
      if (lsdbg.and.ilsdbg(1).gt.0) then
         write (rec,fmt=99996) unitgz, irefn, grznrm, dinky
         call x04bay(iprint,3,rec)
      end if
c
c     if the projected gradient  z'g  is small and rz is of full
c     rank, x is a minimum on the working set.  an additional
c     refinement step is allowed to take care of an inaccurate
c     value of dinky.
c
      statpt = .not. singlr.and.grznrm .le. dinky.or.irefn .gt.
     *         mrefn
c
      if (.not. statpt) then
c        ---------------------------------------------------------
c        compute a search direction.
c        ---------------------------------------------------------
         prnt = .true.
c
         error = iter .ge. itmax
         if (.not. error) then
c
            irefn = irefn + 1
            iter = iter + 1
c
            if (iter.eq.idbg) then
               lsdbg = .true.
               cmdbg = lsdbg
               msglvl = msgsvd
            end if
c
            call e04ncq(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,w(lap),
     *                  w(lres),w(lhz),w(lpx),w(lgq),w(lcq),r,w(lzy),
     *                  w(lwrk))
c
c           ------------------------------------------------------
c           find the constraint we bump into along p.
c           update x and ax if the step alfa is nonzero.
c           ------------------------------------------------------
c           alfhit is initialized to bigalf.  if it remains
c           that way after the call to e04ucg, it will be
c           regarded as infinite.
c
            bigalf = adivb (bigdx,pnorm,overfl)
c
            call e04ucg(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfhit,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  w(lanorm),w(lap),ax,bl,bu,featol,w(lpx),x)
c
c           if  rz1  is nonsingular,  alfa = 1.0  will be the
c           step to the least-squares minimizer on the
c           current subspace. if the unit step does not violate
c           the nearest constraint by more than featol,  the
c           constraint is not added to the working set.
c
            hitcon = singlr.or.palfa .le. one
            uncon = .not. hitcon
c
            if (hitcon) then
               alfa = alfhit
            else
               jadd = 0
               alfa = one
            end if
c
c           check for an unbounded solution or negligible step.
c
            unbndd = alfa .ge. bigalf
            stall = abs(alfa*pnorm) .le. epspt9*xnorm
            if (stall) then
               nstall = nstall + 1
               cyclin = nstall .gt. mstall
            else
               nstall = 0
            end if
c
            error = unbndd.or.cyclin
            if (.not. error) then
c              ---------------------------------------------------
c              set x = x + alfa*p.  update ax, gq, res and ctx.
c              ---------------------------------------------------
               if (alfa.ne.zero) call e04ncl(hitcon,hitlow,linobj,
     *                                unitgz,nclin,nrank,nrz,n,ldr,jadd,
     *                                numinf,alfa,ctp,ctx,xnorm,w(lap),
     *                                ax,bl,bu,w(lgq),w(lhz),w(lpx),
     *                                w(lres),r,x,w(lwrk))
c
               if (hitcon) then
c                 ------------------------------------------------
c                 add a constraint to the working set.
c                 update the tq factors of the working set.
c                 use p as temporary work space.
c                 ------------------------------------------------
c                 update  istate.
c
                  if (bl(jadd).eq.bu(jadd)) then
                     istate(jadd) = 3
                  else if (hitlow) then
                     istate(jadd) = 1
                  else
                     istate(jadd) = 2
                  end if
                  iadd = jadd - n
                  if (jadd.le.n) then
c
                     do 40 ifix = 1, nfree
                        if (kx(ifix).eq.jadd) go to 60
   40                continue
                  end if
   60             continue
c
                  call e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,
     *                        nfree,nrank,nres,ngq,n,lda,ldzy,ldr,ldt,
     *                        kx,condmx,a,r,w(lt),w(lres),w(lgq),w(lzy),
     *                        w(lwrk),w(lrlam),w(lpx),msglvl)
c
                  nrz = nrz - 1
                  nz = nz - 1
c
                  if (jadd.le.n) then
c
c                    a simple bound has been added.
c
                     nfree = nfree - 1
                  else
c
c                    a general constraint has been added.
c
                     nactiv = nactiv + 1
                     kactiv(nactiv) = iadd
                  end if
c
                  irefn = 0
               end if
c
c              ---------------------------------------------------
c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  refine the current x and set inform so
c              that feasibility is checked in e04ncp.
c              ---------------------------------------------------
               call e04ncr(n,nclin,istate,bigbnd,cnorm,err1,jmax1,nviol,
     *                     ax,bl,bu,featol,x,w(lwrk))
c
               if (err1.gt.featol(jmax1)) then
                  call e04nch(linobj,rowerr,unitq,nclin,nactiv,nfree,
     *                        nrank,nz,n,nctotl,ldzy,lda,ldr,ldt,istate,
     *                        kactiv,kx,jmax1,err2,ctx,xnorm,a,ax,bl,bu,
     *                        w(lcq),w(lres),w(lres0),featol,r,w(lt),x,
     *                        w(lzy),w(lpx),w(lwrk))
c
                  if (lsdbg.and.ilsdbg(1).gt.0) then
                     write (rec,fmt=99998) err1, err2
                     call x04bay(iprint,2,rec)
                  end if
                  if (rowerr) then
                     if (msglvl.gt.0) then
                        write (rec,fmt=99997)
                        call x04baf(iprint,rec(1))
                     end if
                     numinf = 1
                     error = .true.
                  else
                     numinf = -1
                     uncon = .false.
                     irefn = 0
                  end if
               end if
               needfg = alfa .ne. zero
            end if
         end if
      end if
c
c        until      statpt .or. error
      if (.not. (statpt.or.error)) go to 20
c
c     ===============================================================
c     try and find the index jdel of a constraint to drop from
c     the working set.
c     ===============================================================
      jdel = 0
c
      if (numinf.eq.0.and.prbtyp.eq.'fp') then
         if (n.gt.nz) call sload (n-nz,(zero),w(lrlam),1)
         jtiny = 0
         jsmlst = 0
         jbigst = 0
      else
c
         call e04nck(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,nrz,
     *               istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,jtiny,
     *               jbigst,kbigst,trulam,a,w(lanorm),w(lgq),w(lrlam),
     *               w(lt),w(lwtinf))
c
      end if
c
      if (.not. error) then
         if (jsmlst.gt.0) then
c
c           e04nck found a regular constraint with multiplier less
c           than (-dinky).
c
            jdel = jsmlst
            kdel = ksmlst
            isdel = istate(jdel)
            istate(jdel) = 0
c
         else if (jsmlst.lt.0) then
c
            jdel = jsmlst
c
         else if (numinf.gt.0.and.jbigst.gt.0) then
c
c           no feasible point exists for the constraints but the
c           sum of the constraint violations may be reduced by
c           moving off constraints with multipliers greater than 1.
c
            jdel = jbigst
            kdel = kbigst
            isdel = istate(jdel)
            if (trulam.le.zero) is = -1
            if (trulam.gt.zero) is = -2
            istate(jdel) = is
            firstv = .true.
            numinf = numinf + 1
         end if
c
         if (jdel.ne.0.and.singlr) then
c
c           cannot delete a constraint when rz is singular.
c           probably a weak minimum.
c
            jdel = 0
         else if (jdel.ne.0) then
c
c           constraint jdel has been deleted.
c           update the matrix factorizations.
c
            call e04nct(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,w(lres),r,
     *                  w(lt),w(lgq),w(lzy),w(lwrk),w(lpx))
         end if
      end if
c
      irefn = 0
      convrg = jdel .eq. 0
c
      prnt = .false.
      uncon = .false.
      needfg = .false.
c
c     until       convrg .or. error
      if (.not. (convrg.or.error)) go to 20
c
c     .....................end of main loop........................
c
      weak = jtiny .gt. 0.or.singlr
c
      if (error) then
         if (unbndd) then
            inform = 2
            if (numinf.gt.0) inform = 3
         else if (iter.ge.itmax) then
            inform = 4
         else if (cyclin) then
            inform = 5
         end if
      else if (convrg) then
         inform = 0
         if (numinf.gt.0) then
            inform = 3
         else if (prbtyp.ne.'fp'.and.weak) then
            inform = 1
         end if
      end if
c

c     set   clamda.  print the full solution.

      msglvl = msgsvd
      if (msglvl.gt.0) then
         write (rec,fmt=99999) prbtyp, iter
         call x04bay(iprint,2,rec)
      end if
c
      call e04nbx(msglvl,nfree,lda,n,nclin,nctotl,bigbnd,named,names,
     *            nactiv,istate,kactiv,kx,a,bl,bu,clamda,w(lrlam),x)

c     end of  e04ncz. (lscore)

99999 format (/' exit from ',a2,' problem after ',i5,' iterations.')
99998 format (' xxx  iterative refinement.  maximum errors before and ',
     *       'after refinement are',/6x,1p,2d14.2)
99997 format (' xxx  warning.  cannot satisfy the constraints to the a',
     *       'ccuracy requested.')
99996 format (/' //e04ncz//  unitgz irefn     grznrm      dinky   ',
     *       /' //e04ncz//  ',l6,i6,1p,2d11.2)
99995 format (/' //e04ncz//  singlr   abs(rzz1)      drzmax      drzmin'
     *       ,/' //e04ncz//  ',l6,1p,3d12.4)
      end


      subroutine e04ncx(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldzy,lda,
     *                  ldr,ldt,istate,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)

c     e04ncx updates the factor r as kx is reordered to reflect the
c     status of the bound constraints given by istate.  kx is reordered
c     so that the fixed variables come last.  one of two alternative
c     are used to reorder kx. one method needs fewer accesses to kx, the
c     other gives a matrix rz with more rows and columns.

      double precision condmx
      integer inform, lda, ldr, ldt, ldzy, msglvl, n, nfree,
     *                  ngq, nrank, nres, nz
      logical           unitq

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer istate(*), kx(n)

      integer iadd, ifix, j, j2, jadd, k, l, lstart, nactv,
     *                  nfixed

      nfixed = n - nfree
c
      if (nrank.lt.n.and.nrank.gt.0) then

c        r is specified but singular.  try and keep the dimension of rz
c        as large as possible.

         nactv = 0
         nfree = n
         nz = n
c
         j = n
c        +       while (j .gt. 0 .and. n-nfree .lt. nfixed) do
   20    if (j.gt.0.and.n-nfree.lt.nfixed) then
            if (istate(j).gt.0) then
               jadd = j
               do 40 ifix = nfree, 1, -1
                  if (kx(ifix).eq.jadd) go to 60
   40          continue
c
c              add bound jadd.
c
   60          call e04ncv(unitq,inform,ifix,iadd,jadd,nactv,nz,nfree,
     *                     nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,
     *                     a,r,t,res,gq,zy,w,c,s,msglvl)
c
               nfree = nfree - 1
               nz = nz - 1
            end if
            j = j - 1
            go to 20
c           +       end while
         end if
      else

c        r is of full rank,  or is not specified.

         if (nfixed.gt.0) then
c
c           order kx so that the free variables come first.
c
            lstart = nfree + 1
            do 120 k = 1, nfree
               j = kx(k)
               if (istate(j).gt.0) then
                  do 80 l = lstart, n
                     j2 = kx(l)
                     if (istate(j2).eq.0) go to 100
   80             continue
c
  100             kx(k) = j2
                  kx(l) = j
                  lstart = l + 1
c
                  if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,k,l,r,
     *                                 res,c,s)
               end if
  120       continue
c
         end if
         nz = nfree
      end if

c     end of  e04ncx. (lsbnds)

      end


      subroutine e04ncu(cold,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *                  lda,istate,kactiv,bigbnd,tolact,a,ax,bl,bu,x,wx,
     *                  work)

c     e04ncu  computes the quantities  istate (optionally), kactiv,
c     nactiv, nz and nfree  associated with the working set at x.
c     the computation depends upon the value of the input parameter
c     cold,  as follows...
c
c     cold = true.  an initial working set will be selected. first,
c                   nearly-satisfied or violated bounds are added.
c                   next,  general linear constraints are added that
c                   have small residuals.
c
c     cold = false. the quantities kactiv, nactiv, nz and nfree are
c                   computed from istate,  specified by the user.
c
c     values of istate(j)....
c
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu

      implicit none

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision bigbnd, tolact
      integer lda, n, nactiv, nartif, nclin, nctotl, nfree
      logical           cold, vertex

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  work(n), wx(n), x(n)
      integer istate(nctotl), kactiv(n)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      double precision wmach
      integer ilsdbg(5)

      double precision b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, toobig
      integer i, imin, is, j, jmin, k, nfixed

      character*80      rec(4)

      double precision ddot1
      external          ddot1

c     .. intrinsic functions ..
      intrinsic         abs, min

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common/ cstmch /wmach(9)
      common            /ce04nc/ilsdbg, lsdbg

      flmax = wmach(7)
      biglow = -bigbnd
      bigupp = bigbnd
c

c     move the variables inside their bounds.

c
      do 20 j = 1, n
         b1 = bl(j)
         b2 = bu(j)
c
         if (b1.gt.biglow) then
            if (x(j).lt.b1) x(j) = b1
         end if
c
         if (b2.lt.bigupp) then
            if (x(j).gt.b2) x(j) = b2
         end if
   20 continue
c
      call dcopy(n,x,1,wx,1)
c
      if (lsdbg) then
         if (ilsdbg(1).gt.0) then
            write (rec,fmt=99999) cold, nclin, nctotl
            call x04bay(iprint,3,rec)
         end if
         if (ilsdbg(2).gt.0) then
            write (rec,fmt=99998)
            call x04bay(iprint,2,rec)
            do 40 i = 1, n, 5
               write (rec,fmt=99995) (wx(j),j=i,min(i+4,n))
               call x04baf(iprint,rec(1))
   40       continue
         end if
      end if
c
      nfixed = 0
      nactiv = 0
      nartif = 0
c
c     if a cold start is being made, initialize  istate.
c     if  bl(j) = bu(j),  set  istate(j)=3  for all variables and linear
c     constraints.
c
      if (cold) then
         do 60 j = 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
   60    continue
      else
         do 80 j = 1, nctotl
            if (istate(j).gt.3.or.istate(j).lt.0) istate(j) = 0
            if (bl(j).ne.bu(j).and.istate(j).eq.3) istate(j) = 0
   80    continue
      end if
c
c     initialize nfixed, nfree and kactiv.
c     ensure that the number of bounds and general constraints in the
c     working set does not exceed n.
c
      do 100 j = 1, nctotl
         if (nfixed+nactiv.eq.n) istate(j) = 0
         if (istate(j).gt.0) then
            if (j.le.n) then
               nfixed = nfixed + 1
               if (istate(j).eq.1) wx(j) = bl(j)
               if (istate(j).ge.2) wx(j) = bu(j)
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
  100 continue
c

c     if a cold start is required,  attempt to add as many
c     constraints as possible to the working set.

      if (cold) then
c
c        see if any bounds are violated or nearly satisfied.
c        if so,  add these bounds to the working set and set the
c        variables exactly on their bounds.
c
         j = n
c        +       while (j .ge. 1 .and. nfixed + nactiv .lt. n) do
  120    if (j.ge.1.and.nfixed+nactiv.lt.n) then
            if (istate(j).eq.0) then
               b1 = bl(j)
               b2 = bu(j)
               is = 0
               if (b1.gt.biglow) then
                  if (wx(j)-b1.le.(one+abs(b1))*tolact) is = 1
               end if
               if (b2.lt.bigupp) then
                  if (b2-wx(j).le.(one+abs(b2))*tolact) is = 2
               end if
               if (is.gt.0) then
                  istate(j) = is
                  if (is.eq.1) wx(j) = b1
                  if (is.eq.2) wx(j) = b2
                  nfixed = nfixed + 1
               end if
            end if
            j = j - 1
            go to 120
c           +       end while
         end if
c

c        the following loop finds the linear constraint (if any) with
c        smallest residual less than or equal to tolact  and adds it
c        to the working set.  this is repeated until the working set
c        is complete or all the remaining residuals are too large.

c        first, compute the residuals for all the constraints not in the
c        working set.
c
         if (nclin.gt.0.and.nactiv+nfixed.lt.n) then
            do 140 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot1 (n,a(i,1),lda,wx)
  140       continue
c
            is = 1
            toobig = tolact + tolact
c
c           + while (is .gt. 0 .and. nfixed + nactiv .lt. n) do
  160       if (is.gt.0.and.nfixed+nactiv.lt.n) then
               is = 0
               resmin = tolact
c
               do 180 i = 1, nclin
                  j = n + i
                  if (istate(j).eq.0) then
                     b1 = bl(j)
                     b2 = bu(j)
                     resl = toobig
                     resu = toobig
                     if (b1.gt.biglow) resl = abs(ax(i)-b1)/(one+abs(b1)
     *                                       )
                     if (b2.lt.bigupp) resu = abs(ax(i)-b2)/(one+abs(b2)
     *                                       )
                     residl = min(resl,resu)
                     if (residl.lt.resmin) then
                        resmin = residl
                        imin = i
                        is = 1
                        if (resl.gt.resu) is = 2
                     end if
                  end if
  180          continue
c
               if (is.gt.0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j = n + imin
                  istate(j) = is
               end if
               go to 160
c              +          end while
            end if
         end if
      end if
c
      if (vertex.and.nactiv+nfixed.lt.n) then

c        find an initial vertex by temporarily fixing some variables.

c        compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).

         do 220 j = 1, n
            if (istate(j).eq.0) then
               colsiz = zero
               do 200 k = 1, nclin
                  if (istate(n+k).gt.0) colsiz = colsiz + abs(a(k,j))
  200          continue
               work(j) = colsiz
            end if
  220    continue

c        find the  nartif  smallest such columns.
c        this is an expensive loop.  later we can replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).
c        (this comment written in 1980).

c        +       while (nfixed + nactiv .lt. n) do
  240    if (nfixed+nactiv.lt.n) then
            colmin = flmax
            do 260 j = 1, n
               if (istate(j).eq.0) then
                  if (nclin.eq.0) go to 280
                  colsiz = work(j)
                  if (colmin.gt.colsiz) then
                     colmin = colsiz
                     jmin = j
                  end if
               end if
  260       continue
            j = jmin
  280       istate(j) = 4
            nartif = nartif + 1
            nfixed = nfixed + 1
            go to 240
c           +       end while
         end if
      end if
c
      nfree = n - nfixed
c
      if (lsdbg) then
         if (ilsdbg(1).gt.0) then
            write (rec,fmt=99996) nfixed, nactiv, nartif
            call x04bay(iprint,4,rec)
         end if
         if (ilsdbg(2).gt.0) then
            write (rec,fmt=99997)
            call x04bay(iprint,2,rec)
            do 300 i = 1, n, 5
               write (rec,fmt=99995) (wx(j),j=i,min(i+4,n))
               call x04baf(iprint,rec(1))
  300       continue
         end if
      end if

c     end of  e04ncu. (lscrsh)
c
99999 format (/' //e04ncu// cold nclin nctotl',/' //e04ncu// ',l4,i6,i7)
99998 format (/' //e04ncu// variables before crash... ')
99997 format (/' //e04ncu// variables after  crash... ')
99996 format (/' //e04ncu// working set selected ...             ',
     *       /' //e04ncu// nfixed nactiv nartif      ',/' //e04ncu// ',
     *       i6,2i7)
99995 format (5g12.3)
      end


      subroutine e04nbt(mode,nrowt,n,t,y)

c     e04nbt  solves equations involving a reverse-triangular matrix  t
c     and a right-hand-side vector  y,  returning the solution in  y.

      double precision zero
      parameter (zero=0.0d+0)

      integer mode, n, nrowt

      double precision t(nrowt,*), y(n)

      double precision yj
      integer j, jj, l, n1

      n1 = n + 1
      if (mode.eq.1) then

c        mode = 1  ---  solve  t * y(new) = y(old).

         do 20 j = 1, n
            jj = n1 - j
            yj = y(j)/t(j,jj)
            y(j) = yj
            l = jj - 1
            if (l.gt.0.and.yj.ne.zero) call daxpy (l,(-yj),t(j+1,jj),1
     *          ,y(j+1),1)
   20    continue

      else
c
c        mode = 2  ---  solve  t' y(new) = y(old).
c
         do 40 j = 1, n
            jj = n1 - j
            yj = y(j)/t(jj,j)
            y(j) = yj
            l = jj - 1
            if (l.gt.0.and.yj.ne.zero) call daxpy (l,(-yj),t(jj,j+1),
     *          nrowt,y(j+1),1)
   40    continue
      end if
c
c     reverse the solution vector.
c
      if (n.gt.1) then
         l = n/2
         do 60 j = 1, l
            jj = n1 - j
            yj = y(j)
            y(j) = y(jj)
            y(jj) = yj
   60    continue
      end if

c     end of  e04nbt. (cmtsol)

      end


      subroutine e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s,msglvl)

c     e04ncv  updates the factorization,  a(free) * (z y) = (0 t),  when
c     a constraint is added to the working set.  if  nrank .gt. 0, the
c     factorization  (r) = pcq  is also updated,  where  c  is the
c                    (0)
c     least squares matrix,  r  is upper-triangular,  and  p  is an
c     orthogonal matrix.  the matrices  c  and  p  are not stored.
c
c     there are three separate cases to consider (although each case
c     shares code with another)...
c
c     (1) a free variable becomes fixed on one of its bounds when there
c         are already some general constraints in the working set.
c
c     (2) a free variable becomes fixed on one of its bounds when there
c         are only bound constraints in the working set.
c
c     (3) a general constraint (corresponding to row  iadd  of  a) is
c         added to the working set.
c
c     in cases (1) and (2), we assume that  kx(ifix) = jadd.
c     in all cases,  jadd  is the index of the constraint being added.
c
c     if there are no general constraints in the working set,  the
c     matrix  q = (z y)  is the identity and will not be touched.
c
c     if  nres .gt. 0,  the row transformations are applied to the rows
c     of the  (n by nres)  matrix  res.
c     if  ngq .gt. 0,  the column transformations are applied to the
c     columns of the  (ngq by n)  matrix  gq'.

      implicit none
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision condmx
      integer iadd, ifix, inform, jadd, lda, ldr, ldt, ldzy,
     *                  msglvl, n, nactiv, nfree, ngq, nrank, nres, nz
      logical           unitq

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer kx(n)

      double precision asize, dtmax, dtmin, epspt3, epspt5, epspt8,
     *                  epspt9
      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision cond, condbd, dtnew, tdtmax, tdtmin
      integer i, nanew, nfmin, npiv, nt
      logical           bound, overfl

      character*80      rec(5)

      double precision dnrm2, adivb 
      external          dnrm2, adivb 

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin

c     if the condition estimator of the updated factors is greater than
c     condbd,  a warning message is printed.

      condbd = one/epspt9
      overfl = .false.
      bound = jadd .le. n

      if (bound) then

c        a simple bound has entered the working set.  iadd  is not used.

         if (lsdbg.and.ilsdbg(1).gt.0) then
            write (rec,fmt=99999) nactiv, nz, nfree, ifix, jadd, unitq
            call x04bay(iprint,4,rec)
         end if
         nanew = nactiv
c
         if (unitq) then
c
c           q  is not stored, but kx defines an ordering of the columns
c           of the identity matrix that implicitly define  q.
c           define the sequence of pairwise interchanges p that moves
c           the newly-fixed variable to position nfree.
c           reorder kx accordingly.
c
            do 20 i = 1, nfree - 1
               if (i.ge.ifix) then
                  w(i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
   20       continue
c
         else

c           q  is stored explicitly.

c           set  w = the  (ifix)-th  row of  q.
c           move the  (nfree)-th  row of  q  to position  ifix.
c
            call dcopy(nfree,zy(ifix,1),ldzy,w,1)
            if (ifix.lt.nfree) then
               call dcopy(nfree,zy(nfree,1),ldzy,zy(ifix,1),ldzy)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else

c        a general constraint has entered the working set.
c        ifix  is not used.

         if (lsdbg.and.ilsdbg(1).gt.0) then
            write (rec,fmt=99998) nactiv, nz, nfree, iadd, jadd, unitq
            call x04bay (iprint,4,rec)
         end if
c
         nanew = nactiv + 1
c
c        transform the incoming row of  a  by  q'.  use c as workspace.
c
         call dcopy (n,a(iadd,1),lda,w,1)
         call cmqmul (8,n,nz,nfree,ldzy,unitq,kx,w,zy,c)
c
c        check that the incoming row is not dependent upon those
c        already in the working set.
c
         dtnew = dnrm2 (nz,w,1)
         if (nactiv.eq.0) then
c
c           this is the only general constraint in the working set.
c
            cond = adivb (asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else
c
c           there are already some general constraints in the working
c           set. update the estimate of the condition number.
c
            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = adivb (tdtmax,tdtmin,overfl)
         end if
c
         if (cond.gt.condmx.or.overfl) go to 60
c
         if (unitq) then
c
c           first general constraint added.  set  q = i.
c
            call smload('general',nfree,nfree,zero,one,zy,ldzy)
            unitq = .false.
         end if
      end if
c
      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if
c
      nt = min(nrank,npiv)
c
      if (unitq) then

c        q (i.e., zy) is not stored explicitly.
c        apply the sequence of pairwise interchanges p that moves the
c        newly-fixed variable to position nfree.

         if (ngq.gt.0) call sgeapr ('left','transpose',nfree-1,w,ngq,gq,
     *                             n)
c
         if (nrank.gt.0) then
c
c           apply the pairwise interchanges to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1), s(2), ..., s(nt-1).
c
            call sutsr1('right',n,ifix,nt,s,r,ldr)
c
            if (nt.lt.npiv) then
c
c              r is upper trapezoidal.  apply the interchanges in
c              columns  nt  thru  npiv.
c
               do 40 i = ifix, nt - 1
                  w(i) = i
   40          continue
c
               call sgeapr ('right','normal',nfree-1,w,nt,r,ldr)
            end if
c
c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.
c
            call suhqr ('left ',n,ifix,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc ('left','variable','forwards',nt,
     *                                 nres,ifix,nt,c,s,res,n)
         end if
      else

c        full matrix q.  define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1,2), (2,3), ...,
c        (npiv-1,npiv).  the rotations must be applied to zy, r, t
c        and gq'.

         call ssrotg ('varble','forwrds',npiv-1,w(npiv),w,1,c,s)
c
         if (bound.and.nactiv.gt.0) then
c
            call dcopy(nactiv,s(nz),1,w(nz),1)
c
            s(nz) = s(nz)*t(nactiv,nz+1)
            t(nactiv,nz+1) = c(nz)*t(nactiv,nz+1)
c
            call f06qzz('create',nactiv,1,nactiv,c(nz+1),s(nz+1),
     *                  t(1,nz+1),ldt)
            call dcopy(nactiv,s(nz),1,t(nactiv,nz),ldt-1)
c
            call dcopy(nactiv,w(nz),1,s(nz),1)
         end if
c
         if (ngq.gt.0) call sgesrc ('left ','variable','forwards',npiv,
     *                             ngq,1,npiv,c,s,gq,n)
         call sgesrc ('right','variable','forwards',nfree,nfree,1,npiv,c
     *               ,s,zy,ldzy)
c
         if (nrank.gt.0) then
c
c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nt-1).
c
            nt = min(nrank,npiv)
            call sutsrh ('right',n,1,nt,c,s,r,ldr)
c
            if (nt.lt.npiv) then
c
c              r is upper trapezoidal.  pretend r is (nt x n) and
c              apply the rotations in columns  nt  thru  npiv.
c
               call sgesrc ('right','variable','forwards',nt,n,nt,npiv,c
     *                     ,s,r,ldr)
            end if
c
c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.
c
            call suhqr ('left ',n,1,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc ('left','variable','forwards',nt,
     *                                 nres,1,nt,c,s,res,n)
         end if
c
         if (bound) then
c
c           the last row and column of zy has been transformed to plus
c           or minus the unit vector e(nfree).  we can reconstitute the
c           columns of gq and r corresponding to the new fixed variable.
c
            if (w(nfree).lt.zero) then
               nfmin = min(nrank,nfree)
               if (nfmin.gt.0) call dscal(nfmin,-one,r(1,nfree),1)
               if (ngq.gt.0) call dscal(ngq,-one,gq(nfree,1),n)
            end if
c

c           the diagonals of t have been altered.  recompute the
c           largest and smallest values.

            if (nactiv.gt.0) then
               call scond (nactiv,t(nactiv,nz),ldt-1,tdtmax,tdtmin)
               cond = adivb (tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint.  install the new row of t.

            call dcopy(nanew,w(nz),1,t(nanew,nz),ldt)
         end if
      end if
c

c     prepare to exit.  check the magnitude of the condition estimator.

   60 if (nanew.gt.0) then
         if (cond.lt.condmx.and..not. overfl) then
c
c           the factorization has been successfully updated.
c
            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin
            if (cond.ge.condbd) then
               if (msglvl.gt.0) then
                  write (rec,fmt=99997) jadd
                  call x04bay(iprint,5,rec)
               end if
            end if
         else
c
c           the proposed working set appears to be linearly dependent.
c
            inform = 1
            if (lsdbg.and.ilsdbg(1).gt.0) then
               write (rec,fmt=99996)
               call x04bay(iprint,2,rec)
               if (bound) then
                  write (rec,fmt=99995) asize, dtmax, dtmin
                  call x04bay(iprint,3,rec)
               else
                  if (nactiv.gt.0) then
                     write (rec,fmt=99994) asize, dtmax, dtmin, dtnew
                     call x04bay(iprint,3,rec)
                  else
                     write (rec,fmt=99993) asize, dtnew
                     call x04bay(iprint,3,rec)
                  end if
               end if
            end if
         end if
      end if

c     end of  e04ncv. (lsadd)

99999 format (/' //e04ncv //  simple bound added.',/' //e04ncv //  nac',
     *       'tiv    nz nfree  ifix  jadd unitq',/' //e04ncv //  ',5i6,
     *       l6)
99998 format (/' //e04ncv //  general constraint added.           ',
     *       /' //e04ncv //  nactiv    nz nfree  iadd  jadd unitq',
     *       /' //e04ncv //  ',5i6,l6)
99997 format (/' xxx  serious ill-conditioning in the working set afte',
     *       'r adding constraint ',i5,/' xxx  overflow may occur in s',
     *       'ubsequent iterations.',//)
99996 format (/' //e04ncv //  dependent constraint rejected.')
99995 format (/' //e04ncv //     asize     dtmax     dtmin        ',
     *       /' //e04ncv //',1p,3d10.2)
99994 format (/' //e04ncv //     asize     dtmax     dtmin     dtnew',
     *       /' //e04ncv //',1p,4d10.2)
99993 format (/' //e04ncv //     asize     dtnew',/' //e04ncv //',1p,
     *       2d10.2)
      end


      subroutine e04ncy(unitq,vertex,inform,k1,k2,nactiv,nartif,nz,
     *                  nfree,nrank,nrejtd,nres,ngq,n,ldzy,lda,ldr,ldt,
     *                  istate,kactiv,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)

c     e04ncy  includes general constraints k1 thru k2 as new rows of
c     the tq factorization stored in t, zy.  if nrank is nonzero, the
c     changes in q are reflected in nrank by n triangular factor r such
c     that
c                         c  =  p (r) q,
c                                 (0)
c     where  p  is orthogonal.

      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision condmx
      integer inform, k1, k2, lda, ldr, ldt, ldzy, msglvl, n,
     *                  nactiv, nartif, nfree, ngq, nrank, nrejtd, nres,
     *                  nz
      logical           unitq, vertex

      double precision a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer istate(*), kactiv(n), kx(n)

      double precision asize, dtmax, dtmin

      double precision wmach

      double precision cndmax, rnorm, rowmax, rtmax
      integer i, iadd, iartif, ifix, iswap, jadd, k, l, nzadd
      double precision dnrm2
      external          dnrm2

      common/ cstmch /wmach(9)
      common            /de04nb/asize, dtmax, dtmin

      rtmax = wmach(8)

c     estimate the condition number of the constraints that are not
c     to be refactorized.

      if (nactiv.eq.0) then
         dtmax = zero
         dtmin = one
      else
         call scond (nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
      end if
c
      do 20 k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv.lt.nfree) then
c
            call e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s,msglvl)
c
            if (inform.eq.0) then
               nactiv = nactiv + 1
               nz = nz - 1
            else
               istate(jadd) = 0
               kactiv(k) = -kactiv(k)
            end if
         end if
   20 continue

      if (nactiv.lt.k2) then

c        some of the constraints were classed as dependent and not
c        included in the factorization.  re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.

         l = k1 - 1
         do 40 k = k1, k2
            i = kactiv(k)
            if (i.ge.0) then
               l = l + 1
               if (l.ne.k) then
                  iswap = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
   40    continue

c        if a vertex is required, add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.

         if (vertex) then
            cndmax = rtmax
            nzadd = nz
            do 80 iartif = 1, nzadd
               if (unitq) then
                  ifix = nfree
                  jadd = kx(ifix)
               else
                  rowmax = zero
                  do 60 i = 1, nfree
                     rnorm = dnrm2(nz,zy(i,1),ldzy)
                     if (rowmax.lt.rnorm) then
                        rowmax = rnorm
                        ifix = i
                     end if
   60             continue
                  jadd = kx(ifix)
c
                  call e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,
     *                        nfree,nrank,nres,ngq,n,lda,ldzy,ldr,ldt,
     *                        kx,cndmax,a,r,t,res,gq,zy,w,c,s,msglvl)
c
               end if
               nfree = nfree - 1
               nz = nz - 1
               nartif = nartif + 1
               istate(jadd) = 4
   80       continue
         end if
      end if
c
      nrejtd = k2 - nactiv
c
      return
c
c     end of  e04ncy. (lsadds)
c
      end

      subroutine e04nch(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,
     *                  n,nctotl,ldzy,lda,ldr,ldt,istate,kactiv,kx,jmax,
     *                  errmax,ctx,xnorm,a,ax,bl,bu,cq,res,res0,featol,
     *                  r,t,x,zy,p,work)

c     e04nch  computes the point on a working set that is closest to the
c     input vector  x  (in the least-squares sense).  the norm of  x,
c     the transformed residual vector  pr - rq'x,  and the constraint
c     values ax  are also initialized.
c
c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable  rowerr
c     is set.

      implicit none

      integer ntry
      parameter (ntry=5)
      double precision zero, one
      parameter (zero=0.0d+0,one=1.0d+0)

      double precision ctx, errmax, xnorm
      integer jmax, lda, ldr, ldt, ldzy, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nz
      logical           linobj, rowerr, unitq

      double precision a(lda,*), ax(*), bl(nctotl), bu(nctotl), cq(*),
     *                  featol(nctotl), p(n), r(ldr,*), res(*), res0(*),
     *                  t(ldt,*), work(nctotl), x(n), zy(ldzy,*)
      integer istate(nctotl), kactiv(n), kx(n)

      integer iprint, isumm, lines1, lines2, nout
      logical           lsdbg

      integer ilsdbg(5)

      double precision bnd
      integer i, is, j, k, ktry

      character*80      rec(2)

      double precision ddot1, dnrm2
      integer idamx1 
      external          ddot1, dnrm2, idamx1 

      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg

c     move  x  onto the simple bounds in the working set.

      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = bl(j)
         if (is.ge.2) bnd = bu(j)
         if (is.ne.4) x(j) = bnd
   20 continue
c

c     move  x  onto the general constraints in the working set.
c     we shall make  ntry  tries at getting acceptable row errors.

      ktry = 1
      jmax = 1
      errmax = zero
c
c     repeat
   40 if (nactiv.gt.0) then
c
c        set  work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.
c
         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(i) = bnd - ddot1 (n,a(k,1),lda,x)
   60    continue
c
         call e04nbt(1,ldt,nactiv,t(1,nz+1),work)
         call sload (n,zero,p,1)
         call dcopy(nactiv,work,1,p(nz+1),1)
c
         call cmqmul(2,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
         call daxpy (n,one,p,1,x,1)
      end if
c
c     ---------------------------------------------------------------
c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.
c     ---------------------------------------------------------------
      xnorm = dnrm2(n,x,1)
      if (nclin.gt.0) call dgemv('n',nclin,n,one,a,lda,x,zero,ax)
c
c     ---------------------------------------------------------------
c     check the row residuals.
c     ---------------------------------------------------------------
      if (nactiv.gt.0) then
         do 80 k = 1, nactiv
            i = kactiv(k)
            j = n + i
            is = istate(j)
            if (is.eq.1) work(k) = bl(j) - ax(i)
            if (is.ge.2) work(k) = bu(j) - ax(i)
   80    continue
c
         jmax = idamx1 (nactiv,work)
         errmax = abs(work(jmax))
      end if
c
      ktry = ktry + 1
c     until    (errmax .le. featol(jmax).or.ktry .gt. ntry
      if (.not. (errmax.le.featol(jmax).or.ktry.gt.ntry)) go to 40
c
      rowerr = errmax .gt. featol(jmax)
c

c     compute the linear objective value  c'x  and the transformed
c     residual  pr  -  rq'x = res0  -  rq'x.

      if (nrank.gt.0.or.linobj) then
         call dcopy(n,x,1,p,1)
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
      end if
c
      ctx = zero
      if (linobj) ctx = ddot1 (n,cq,1,p)
c
      if (nrank.gt.0) then
         call dtrmv('u','n','n',nrank,r,ldr,p,1)
         if (nrank.lt.n) call dgemv('n',nrank,n-nrank,one,r(1,nrank+1),
     *                              ldr,p(nrank+1),one,p)
         call dcopy(nrank,res0,1,res,1)
         call daxpy (nrank,-one,p,1,res,1)
      end if
c
      if (lsdbg.and.ilsdbg(2).gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,2,rec)
         do 100 i = 1, n, 5
            write (rec,fmt=99998) (x(j),j=i,min(i+4,n))
            call x04baf(iprint,rec(1))
  100    continue
      end if
c
      return
c
c
c     end of  e04nch. (lssetx)
c
99999 format (/' //e04nch// variables after refinement ... ')
99998 format (5g12.3)
      end


      SUBROUTINE E04UEF(STRING)

C***********************************************************************
C     E04UEF  loads the option supplied in STRING into the relevant
C     element of IPRMLS, RPRMLS, IPRMNP or RPRMNP.
C***********************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER*(*)     STRING
C     .. Scalars in Common ..
      LOGICAL           NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH
C     .. Local Scalars ..
      INTEGER           NOUT
      LOGICAL           FIRST, PRNT
      CHARACTER*16      KEY
      CHARACTER*72      BUFFER
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Subroutines ..
      EXTERNAL          E04UCQ
C     .. Common blocks ..
      common/ cstmch /wmach(9)
      COMMON            /EE04UC/NEWOPT
C     .. Save statement ..
      save prnt
      SAVE              /EE04UC/
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
C
C     If first time in, set NOUT.
C     NEWOPT is true first time into E04UDF or E04UEF
C     and just after a call to an optimization routine.
C     PRNT is set to true whenever NEWOPT is true.
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         NEWOPT = .TRUE.
      END IF

      NOUT = 6
      BUFFER = STRING
C
C     Call E04UCQ to decode the option and set the parameter value.
C     If NEWOPT is true, reset PRNT and test specially for NOLIST.
C
      IF (NEWOPT) THEN
         NEWOPT = .FALSE.
         PRNT = .TRUE.
         CALL E04UCQ(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'nolist') THEN
            PRNT = .FALSE.
         ELSE
            WRITE (REC,FMT='(// A / A /)') ' Calls to E04UEF',
     *        ' ---------------'
            CALL X04BAY(NOUT,5,REC)
            WRITE (REC,FMT='( 6X, A )') BUFFER
            CALL X04BAF(NOUT,REC(1))
         END IF
      ELSE
         IF (PRNT) THEN
            WRITE (REC,FMT='( 6X, A )') BUFFER
            CALL X04BAF(NOUT,REC(1))
         END IF
         CALL E04UCQ(NOUT,BUFFER,KEY)
C
         IF (KEY.EQ.'list') PRNT = .TRUE.
         IF (KEY.EQ.'nolist') PRNT = .FALSE.
      END IF
C
      RETURN
C
C     End of  E04UEF. (NPOPTN)
C
      END


      subroutine x04bay(nout,nrec,rec)

c     if nrec is 0 then no records are output.

      integer nout, nrec

      character*(*)     rec(*)

      integer i

      do 20 i = 1, nrec
         call x04baf(nout,rec(i))
   20 continue
      end
