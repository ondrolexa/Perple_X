c this file contains:

c 1) blas level 2 fortran subroutines. modified subroutines are named by 
c    appending a digit to the original name or in the case of 6 char
c    names by replacing the final character by a digit (1).
c    blas is a freely-available software package at www.netlib.org/blas/

c 2) lpsol a fortran subroutine, and any non-blas subroutines it invokes,
c    to solve linear programming problems; and nlpsol a fortran subroutine,
c    and any non-blas subroutines it invokes, to solve non-linear programming 
c    problems by succesive quadratic programming from gill p e, hammarling s,
c    murray w, saunders m a and wright m h (1986) user’s guide for lssol
c    (version 1.0) report sol 86–1 department of operations research, 
c    stanford university.

      subroutine lpsol (n,nclin,a,lda,bl,bu,cvec,istate,x,iter,obj,ax,
     *                  clamda,iw,leniw,w,lenw,ifail,
     *                  nsglvl,istart,jtmax2,tol,lpprob)
c----------------------------------------------------------------------
c     lpsol  solves the linear programming problem
c
c           minimize               c' x
c              x
c                                 (  x )
c           subject to    bl  .le.(    ).ge.  bu,
c                                 ( ax )
c
c     where  a  is a constant  nclin by n  matrix.
c     the feasible region is defined by a mixture of linear equality or
c     inequality constraints on  x.
c
c     n  is the number of variables (dimension of x).
c        (n must be positive.)
c
c     nclin  is the number of general linear constraints (rows of  a).
c        (nclin may be zero.)
c
c     the first  n  elements of  bl  and   bu  are lower and upper
c     bounds on the variables.  the next  nclin  elements are
c     lower and upper bounds on the general linear constraints.
c
c     the matrix  a  of coefficients in the general linear constraints
c     is entered as the two-dimensional array  a  (of dimension
c     lda by n).  if nclin = 0, a is not referenced.
c
c     the vector  x  must contain an initial estimate of the solution,
c     and will contain the computed solution on output.
c----------------------------------------------------------------------
      implicit none

      integer nsglvl, jtmax2, istart, lpprob
      double precision tol

      integer           lenlc
      parameter         (lenlc=20)
      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, point3
      parameter         (zero=0.0d+0,point3=3.3d-1)
      double precision  point8, point9, one
      parameter         (point8=0.8d+0,point9=0.9d+0,one=1.0d+0)
      double precision  hundrd
      parameter         (hundrd=1.0d+2)

      double precision  obj
      integer           ifail, iter, lda, leniw, lenw, n, nclin

      double precision  a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  clamda(n+nclin), cvec(*), w(lenw), x(n)
      integer           istate(n+nclin), iw(leniw)
c     .. scalars in common ..
      double precision  alfa, asize, bigbnd, bigdx, bndlow, bndupp,
     *                  dtmax, dtmin, epspt3, epspt5, epspt8, epspt9,
     *                  tolact, tolfea, tolinc, tolrnk, tolx0, trulam
      integer           idbglc, iprint, iprnt, isdel, isumm, isumry,
     *                  itmax1, itmax2, itnfix, jadd, jdel, kchk,
     *                  kcycle, kdegen, lcrash, ldbglc, ldq, ldt,
     *                  lennam, lines1, lines2, lprob, maxact, maxnz,
     *                  mm, msglc, mxfree, ncolt, ndegen, nn, nnclin,
     *                  nout, nprob
      logical           header, lcdbg, prnt
c     .. arrays in common ..
      double precision  rpadlc(23), rpsvlc(mxparm)
      integer           ilcdbg(ldbg), ipadlc(14), ipsvlc(mxparm),
     *                  loclc(lenlc), nfix(2)

      double precision  amin, condmx, epsmch, errmax, feamax, feamin,
     *                  rteps, xnorm
      integer           ianrmj, idbg, inform, it, itmax, j, jinf, jmax,
     *                  lanorm, lcq, ld, ldh, ldr, lfeatu, lgq, litotl,
     *                  lkactv, lkx, llptyp, lq, lr, lrlam, lt, lwrk,
     *                  lwtinf, lwtotl, minact, minfxd, msgdbg, msglvl,
     *                  nact1, nactiv, nartif, ncnln, nctotl, nerr,
     *                  nerror, nfree, ngq, nmoved, nrejtd, nrz, numinf,
     *                  nviol, nz
      logical           cold, cset, done, found, halted, hot, named,
     *                  rowerr, rset, unitq, vertex, warm
      character*2       prbtyp
      character*4       start
      character*6       msg
      character*11      title

      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*8       names(1)
      character*80      errrec(2), rec(2)
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. external subroutines ..
      external          dcopy, dscal, cmsetx, cmprnt, cmdgen,
     *                  cmcrsh, lploc, lpcore, cmqmul,
     *                  rzadds, icopy , sload, scond


      common            /ngg003/loclc
      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg004/lennam, ldt, ncolt, ldq
      common            /ngg005/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg007/alfa, trulam, isdel, jdel, jadd, header,
     *                  prnt
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg009/ilcdbg, lcdbg
      common            /ngg010/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /ngg011/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)
c----------------------------------------------------------------------
c                                 iprint - print level
c                                 istart - 0 - cold start, 1 - warm start, 2 - hot (no benefit)
c                                 jtmax2 - maximum number of iterations, l6
c                                 tol    - feasibility tolerance
c                                 lpprob - problem type
      epsmch = wmach(3)
      rteps = wmach(4)
      msglvl = nsglvl
      itmax2 = jtmax2
      lcrash = istart

      nout = 6
      nerr = 6
c                                 defaults, assigned to iprmlc and 
c                                 rprmlc by equivalence
c                                 ----------------------------------
c                                 feasibility tolerance, sqrt(eps)
      tolfea = tol
c                                 problem type 1 - fp, 2 - lp
      lprob = lpprob

      iprnt = nout
      isumry = -1
      iprint = iprnt
      isumm = isumry
      kchk = 50
      kcycle = 10000
      kdegen = kcycle
      itmax1 = max(50,5*(n+nclin))
      maxact = max(1,min(n,nclin))
      maxnz = n
      mxfree = n

      if (nclin.lt.n) then
         mxfree = nclin + 1
         maxnz = mxfree
      end if

      msgdbg = 0
      idbg = itmax1 + itmax2 + 1
      tolact = 1d-2
      bigbnd = 1.0d+20*0.99999d+0
      bigdx = bigbnd

      lcdbg = .false.
c                                     iprmlc and rprmlc may be changed(?) so
c                                     copy into ipsvlc and rpsvlc
      call icopy (mxparm,iprmlc,1,ipsvlc,1)
      call dcopy(mxparm,rprmlc,1,rpsvlc,1)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      inform = 0
      iter = 0
      header = .true.
      prnt = .true.
      condmx = max(one/epspt5,hundrd)
      llptyp = lprob
      nctotl = n + nclin

c     set all parameters determined by the problem type.

      if (llptyp.eq.1) then
         prbtyp = 'fp'
         cset = .false.
      else if (llptyp.eq.2) then
         prbtyp = 'lp'
         cset = .true.
      end if

c     assign the dimensions of arrays in the parameter list of lpcore.
c     economies of storage are possible if the minimum number of active
c     constraints and the minimum number of fixed variables are known in
c     advance.  the expert user should alter minact and minfxd
c     accordingly.
c     if a linear program is being solved and the matrix of general
c     constraints has fewer rows than columns, i.e.,  nclin .lt. n,
c     a non-zero value is known for minfxd.  note that in this case,
c     vertex must be set  .true..

      vertex = nclin .lt. n

      minfxd = n - mxfree
      minact = mxfree - maxnz

      ldt = max(maxnz,maxact)
      ncolt = mxfree
      if (nclin.eq.0) then
         ldq = 1
      else
         ldq = max(1,mxfree)
      end if
      ldr = ldt

      ncnln = 0
      lennam = 1
      ldh = 1
      mm = 0

c     cold start:  only  x  is provided.
c     warm start:  initial working set is specified in  istate.
c     hot  start:  the work arrays  iw  and  w  are assumed to have been
c                  initialized during a previous run.
c                  the first three elements of  iw  contain details
c                  on the dimension of the initial working set.

      if (lcrash.eq.0) then
         start = 'cold'
      else if (lcrash.eq.1) then
         start = 'warm'
         lcrash = 2
         start  = 'hot '
      else if (lcrash.eq.2) then
         start = 'hot '
      end if

      cold = lcrash .eq. 0
      warm = lcrash .eq. 1
      hot = lcrash .eq. 2

c     allocate remaining work arrays.

      litotl = 3
      lwtotl = 0
      call lploc(cset,n,nclin,litotl,lwtotl)

      lkactv = loclc(1)
      lkx = loclc(2)

      lfeatu = loclc(3)
      lanorm = loclc(4)
      ld = loclc(7)
      lgq = loclc(8)
      lcq = loclc(9)
      lrlam = loclc(10)
      lr = loclc(11)
      lt = loclc(12)
      lq = loclc(13)
      lwtinf = loclc(14)
      lwrk = loclc(15)

c     define the initial feasibility tolerances in clamda.

      if (tolfea.gt.zero) call sload(n+nclin,(tolfea),w(lfeatu),1)

      call cmdgen('initialize anti-cycling variables',msglvl,n,nclin,
     *            nmoved,iter,numinf,istate,bigbnd,ax,bl,bu,clamda,
     *            w(lfeatu),x)

      if (cold .or. warm) then

c        cold or warm start.  just about everything must be initialized.
c        the only exception is istate during a warm start.

         ianrmj = lanorm
         do 20 j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
   20    continue
         if (nclin.gt.0) call scond (nclin,w(lanorm),1,asize,amin)

         call scond (nctotl,w(lfeatu),1,feamax,feamin)
         call dcopy(nctotl,w(lfeatu),1,w(lwtinf),1)
         call dscal(nctotl,(one/feamin),w(lwtinf),1)

c        define the initial working set.
c               nfree ,  nactiv,  kactiv, kx,
c               istate (if start  = 'cold')
c               nartif (if vertex = 'true')

         call cmcrsh(start,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *               lda,istate,iw(lkactv),iw(lkx),bigbnd,tolact,a,ax,
     *               bl,bu,clamda,x,w(lgq),w(lwrk))

c        compute the tq factorization of the working set matrix.

         unitq = .true.
         nz = nfree

         if (nactiv.gt.0) then
            it = nactiv + 1
            nact1 = nactiv
            nactiv = 0
            ngq = 0

            call rzadds(unitq,vertex,1,nact1,it,nactiv,nartif,nz,nfree,
     *                  nrejtd,ngq,n,ldq,lda,ldt,istate,iw(lkactv),
     *                  iw(lkx),condmx,a,w(lt),w(lgq),w(lq),w(lwrk),
     *                  w(ld),w(lrlam),msglvl)
         end if
      else if (hot) then

c        arrays  iw  and  w  have been defined in a previous run.
c        the first three elements of  iw  are  unitq,  nfree and nactiv.

         unitq = iw(1) .eq. 1
         nfree = iw(2)
         nactiv = iw(3)
         nz = nfree - nactiv
      end if

      if (cset) then

c        install the transformed linear term in cq.

         call dcopy(n,cvec,1,w(lcq),1)
         call cmqmul(6,n,nz,nfree,ldq,unitq,iw(lkx),w(lcq),w(lq),w(lwrk)
     *               )
      end if

      rset = .false.
      itmax = itmax2
      jinf = 0

c     +    take your pick when minimizing the sum of infeasibilities:
c     +    nrz    =  nz  implies steepest-descent in the two-norm.
c     +    nrz    =  0   implies steepest-descent in the infinity norm.
      nrz = 0

c     repeat               (until working set residuals are acceptable)

c     move x onto the constraints in the working set.

   40 call cmsetx(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,ldt,
     *            istate,iw(lkactv),iw(lkx),jmax,errmax,xnorm,a,ax,bl,
     *            bu,w(lfeatu),w(lt),x,w(lq),w(ld),w(lwrk))

      if (rowerr) then
         msg = 'infeas'
         numinf = 1
         obj = errmax
         go to 60
      end if

      call lpcore(prbtyp,msg,cset,named,names,rset,unitq,iter,itmax,
     *            jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,nz,istate,
     *            iw(lkactv),iw(lkx),obj,numinf,xnorm,a,ax,bl,bu,
     *            cvec,clamda,w(lfeatu),x,iw,w)

      found = msg .eq. 'feasbl' .or. msg .eq. 'optiml' .or. msg .eq.
     *        'weak  ' .or. msg .eq. 'unbndd' .or. msg .eq. 'infeas'
      halted = msg .eq. 'itnlim'

      if (found) then
         call cmdgen('optimal',msglvl,n,nclin,nmoved,iter,numinf,istate,
     *               bigbnd,ax,bl,bu,clamda,w(lfeatu),x)
      end if

      done = found .and. nviol .eq. 0 .and. nmoved .eq. 0

c     until      done  .or.  halted
      if ( .not. (done .or. halted)) go to 40

c     set   clamda.
c     clean up.  save values for a subsequent hot start.

      call cmprnt(msglvl,nfree,lda,n,nclin,ncnln,nctotl,bigbnd,named,
     *            names,lennam,nactiv,istate,iw(lkactv),iw(lkx),a,bl,bu,
     *            x,clamda,w(lrlam),x)

      iw(1) = 0
      if (unitq) iw(1) = 1
      iw(2) = nfree
      iw(3) = nactiv

   60 if (msg.eq.'optiml') then
         inform = 0
      else if (msg.eq.'feasbl') then
         inform = 0
      else if (msg.eq.'weak  ') then
         inform = 1
      else if (msg.eq.'unbndd') then
         inform = 2
      else if (msg.eq.'infeas') then
         inform = 3
      else if (msg.eq.'itnlim') then
         inform = 4
      else if (msg.eq.'errors') then
         inform = 6
      else if (msg.eq.'noprob') then
         inform = 7
      end if

      if (inform.ge.0) then

         ifail = inform

         if (ifail.lt.3) ifail = 0

         if (ifail.lt.4) then
            istart = 1
         else
            istart = 0
         end if

      else

         call errdbg ('wanola')

      end if
c                                 end of lpsol
      end

      subroutine nlpsol (n,nclin,ncnln,lda,ldcju,ldr,a,bl,bu,confun,
     *                  objfun,iter,istate,c,cjacu,clamda,objf,gradu,r,
     *                  x,iw,leniw,w,lenw,iuser,user,ifail)
c----------------------------------------------------------------------
c     nlpsol   solves the nonlinear program

c            minimize                   f(x)

c                                    (    x  )
c            subject to    bl  .le.  (  a*x  )  .le.  bu
c                                    (  c(x) )

c     where  f(x)  is a smooth scalar function,  a  is a constant matrix
c     and  c(x)  is a vector of smooth nonlinear functions.  the
c     feasible region is defined by a mixture of linear and nonlinear
c     equality or inequality constraints on  x.

c     the dimensions of the problem are...

c     n        the number of variables (dimension of  x),

c     nclin    the number of linear constraints (rows of the matrix  a),

c     ncnln    the number of nonlinear constraints (dimension of  c(x)),

c     nlpsol   uses a sequential quadratic programming algorithm, with a
c     positive-definite quasi-newton approximation to the transformed
c     hessian  q'hq  of the lagrangian function (which will be stored in
c     the array  r).

c     this material is based upon work partially supported by the
c     national science foundation under grants mcs-7926009 and
c     ecs-8312142; the department of energy contract am03-76sf00326,
c     pa no. de-at03-76er72018; the army research office contract
c     daa29-84-k-0156; and the office of naval research grant
c     n00014-75-c-0267.
c----------------------------------------------------------------------
      integer           mxparm
      parameter         (mxparm=30)
      integer           lenls
      parameter         (lenls=20)
      integer           lennp
      parameter         (lennp=35)
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, point3, point8
      parameter         (zero=0.0d+0,point3=3.3d-1,point8=0.8d+0)
      double precision  point9, one
      parameter         (point9=0.9d+0,one=1.0d+0)
      double precision  growth
      parameter         (growth=1.0d+2)

      double precision  objf
      integer           ifail, iter, lda, ldcju, ldr, leniw, lenw, n,
     *                  nclin, ncnln

      double precision  a(lda,*), bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     *                  c(*), cjacu(ldcju,*), clamda(n+nclin+ncnln),
     *                  gradu(n), r(ldr,*), user(*), w(lenw), x(n)
      integer           istate(n+nclin+ncnln), iuser(*), iw(leniw)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer           idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldq, ldt,
     *                  lennam, lfdset, lformh, lines1, lines2, lprob,
     *                  lverfy, lvlder, lvldif, lvrfyc, msgls, msgnp,
     *                  nactiv, ncdiff, ncolt, nfdiff, nfree, nlnf,
     *                  nlnj, nlnx, nload, nn, nnclin, nncnln, nout,
     *                  nprob, nsave, nz
      logical           cmdbg, incrun, lsdbg, npdbg, unitq
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           icmdbg(ldbg), ilsdbg(ldbg), inpdbg(ldbg),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), jverfy(4), locls(lenls),
     *                  locnp(lennp)

      double precision  amin, cond, condmx, ctx, epsmch, errmax, fdchk,
     *                  fdnorm, feamax, feamin, obj, rootn, rteps, ssq1,
     *                  suminf, xnorm
      integer           i, ianrmj, idbg, idbgsv, ikx, info, inform,
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

      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(2)
c     .. external functions ..
      double precision  dnrm2, sdiv , dlantr

      external          dnrm2, sdiv , dlantr
c     .. external subroutines ..
      external          dcopy, dscal, cmqmul, lssetx, lscrsh,
     *                  lsbnds, lsadds, lscore, nploc, npcrsh, npdflt,
     *                  npchkd, npcore, nprset, sgeqr , icopy , sload,
     *                  scond , smcopy, smload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg012/locls
      common            /ngg013/locnp

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg004/lennam, ldt, ncolt, ldq
      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg015/lvrfyc, jverfy
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg016/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ngg018/rcndbd, rfrobn, drmax, drmin
      common            /ngg019/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ngg4fb/icmdbg, cmdbg
      common            /ngg001/nactiv, nfree, nz, unitq
      common            /ngg002/inpdbg, npdbg
      common            /ngg020/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /ngg021/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c----------------------------------------------------------------------
      epsmch = wmach(3)
      rteps = wmach(4)
      nout = 6
      nerr = 6

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      rhomax = one/epsmch
      rootn = sqrt(dble(n))

      inform = 0
c                                 set default parameter values
      call npdflt(n,nclin,ncnln,title)
c                                 set values from user/iuser
c                                 epsrf, function precision
      epsrf = user(1)
c                                 ftol,.optimality tolerance
      ftol = user(2)
c                                 ctol and tolfea,feasibility tolerance
      ctol = user(3)
      tolfea = ctol
c                                 eta, step limit < nopt(5) leads to bad results, coincidence?
      eta = user(5)
c                                 dxlim, linesearch tolerance, low values -> more accurate search -> more function calls
c                                 0.05-.4 seem best
      dxlim = user(4)
c                                 fdint, finite difference interval, forward.
      fdint = user(6)
c                                 cdint, centered finite difference interval
      cdint = fdint**(0.67d0)
      cdint = epsrf**(0.33d0)

      if (fdint.gt.0d0) then
         lfdset = 1
      else 
         lfdset = 0
      end if 
c                                 lverfy, verify level, default off, may be reset. 0 - off, 1 - on
      lverfy = iuser(11)
c                                 lvlder, derivative level, 3 - all available, 1 - some, 0 - none
      lvlder = iuser(13)
c                                 copy the equivalenced parameter arrays iprmls, rprmls, iprmnp and 
c                                 rprmnp to ipsvls, rpsvls, ipsvnp, and rpsvnp so the parameters
c                                 can be reset to their original values (i don't think this is 
c                                 necessary.
      call icopy (mxparm,iprmls,1,ipsvls,1)
      call dcopy(mxparm,rprmls,1,rpsvls,1)
      call icopy (mxparm,iprmnp,1,ipsvnp,1)
      call dcopy(mxparm,rprmnp,1,rpsvnp,1)

      needfd = lvlder .eq. 0 .or. lvlder .eq. 2 .or.
     *         (lvlder.eq.1 .and. ncnln.gt.0)
      cold = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1

      nplin = n + nclin
      nctotl = nplin + ncnln

c     assign the dimensions of arrays in the parameter list of npcore.
c     economies of storage are possible if the minimum number of active
c     constraints and the minimum number of fixed variables are known in
c     advance.  the expert user should alter minact and minfxd
c     accordingly.

      minact = 0
      minfxd = 0

      mxfree = n - minfxd
      maxact = max(1,min(n,nclin))
      maxnz = n - (minfxd+minact)

      if (nclin+ncnln.eq.0) then
         ldq = 1
         ldt = 1
         ncolt = 1
      else
         ldq = max(1,mxfree)
         ldt = max(maxnz,maxact)
         ncolt = mxfree
      end if

      lennam = 1
      m = 1
      ldfju = 2

      ldaqp = max(nclin+ncnln,1)
      if (ncnln.eq.0 .and. nclin.gt.0) ldaqp = lda

c     nploc  defines the arrays that contain the locations of various
c     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      call nploc(n,nclin,ncnln,nctotl,litotl,lwtotl)

c     allocate certain addresses that are not allocated in nploc.

      lax = lwtotl + 1
      lwtotl = lax + nclin - 1
      lax = min(lax,lwtotl)

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

      lkx = locnp(1)
      liperm = locnp(2)
      laqp = locnp(3)
      ldx = locnp(7)
      lfeatl = locnp(10)
      lwrk2 = locnp(12)

      lcmul = locnp(16)
      lrho = locnp(20)
      lwrk3 = locnp(21)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)
      lcjac = locnp(27)
      lgrad = locnp(28)

      ldcj = max(ncnln,1)

      tolrnk = zero
      rcndbd = sqrt(hcndbd)

c     load the arrays of feasibility tolerances.

      if (tolfea.gt.zero) call sload(nplin,tolfea,w(lfeatl),1)

      if (ncnln.gt.0 .and. ctol.gt.zero) call sload(ncnln,ctol,
     *    w(lfeatl+nplin),1)

      if (lfdset.eq.0) then
         fdchk = sqrt(epsrf)
      else if (lfdset.eq.1) then
         fdchk = fdint
      else
         fdchk = w(lhfrwd)
      end if

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

         call npchkd(info,msgnp,nstate,lvlder,nfun,ngrad,ldcj,ldcju,n,
     *               ncnln,confun,objfun,iw(lneedc),bigbnd,epsrf,cdint,
     *               fdint,fdchk,fdnorm,objf,xnorm,bl,bu,c,w(lwrk3),
     *               w(lcjac),cjacu,w(lcjdx),w(ldx),w(lgrad),gradu,
     *               w(lhfrwd),w(lhctrl),x,w(lwrk1),w(lwrk2),w,lenw,
     *               iuser,user)

         if (info.ne.0) then
            if (info.gt.0) inform = 7
            if (info.lt.0) inform = info
            go to 80
         end if
         nstate = 0
      end if

      call icopy (ldbg,ilsdbg,1,icmdbg,1)

      if (nclin.gt.0) then
         ianrmj = lanorm
         do 20 j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
   20    continue
         call scond (nclin,w(lanorm),1,asize,amin)
      end if

      call scond (nplin,w(lfeatl),1,feamax,feamin)
      call dcopy(nplin,w(lfeatl),1,w(lwtinf),1)
      call dscal(nplin,(one/feamin),w(lwtinf),1)

c     the input values of x and (optionally)  istate are used by
c     lscrsh  to define an initial working set.

      vertex = .false.
      call lscrsh(cold,vertex,nclin,nplin,nactiv,nartif,nfree,n,lda,
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
            if (ncnln.gt.0) call sload(ncnln,(zero),w(lcmul),1)
         else

c           r will be updated while finding a feasible x.

            nrank = nlnx
            call sload(nlnx,(zero),w(lres0),1)
            if (ncnln.gt.0) call dcopy(ncnln,clamda(nplin+1),1,w(lcmul),
     *                                 1)

         end if

         incrun = .true.
         rhonrm = zero
         rhodmp = one
         scale = one
         call sload(ncnln,(zero),w(lrho),1)

c        re-order kx so that the free variables come first.
c        if a warm start is required, nrank will be nonzero and the
c        factor r will be updated.

         call lsbnds(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldq,lda,ldr,
     *               ldt,istate,iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)

      end if

c     factorize the linear constraints in the initial working set.

      if (nactiv.gt.0) then
         nact1 = nactiv
         nactiv = 0

         call lsadds(unitq,vertex,inform,1,nact1,nactiv,nartif,nz,nfree,
     *               nrank,nrejtd,nres,ngq,n,ldq,lda,ldr,ldt,istate,
     *               iw(lkactv),iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)
      end if

      if (lcrash.le.1) then

c        cold or warm start.  move  x  on to the linear constraints and
c        find a feasible point.

         ssq1 = zero
         linobj = .false.
         call lssetx(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,n,
     *               nplin,ldq,lda,ldr,ldt,istate,iw(lkactv),iw(lkx),
     *               jmax,errmax,ctx,xnorm,a,w(lax),bl,bu,w(lgq),w(lres)
     *               ,w(lres0),w(lfeatl),r,w(lt),x,w(lq),w(lwrk1),
     *               w(lwrk2))

c        call  lscore  to find a feasible  x.
c        use  work2  as the multiplier vector.

         jinf = 0
         lclam = lwrk2

         itmxsv = itmax1
         itmax1 = nminor

         call lscore('fp problem',named,names,linobj,unitq,nlperr,itns,
     *               jinf,nclin,nplin,nactiv,nfree,nrank,nz,nz1,n,lda,
     *               ldr,istate,iw(lkactv),iw(lkx),ctx,obj,ssq1,suminf,
     *               numinf,xnorm,bl,bu,a,w(lclam),w(lax),w(lfeatl),r,x,
     *               w)

         itmax1 = itmxsv

         if (nlperr.gt.0) then
            inform = 2
            go to 80
         end if

         call icopy (ldbg,inpdbg,1,icmdbg,1)

      end if

      if (lcrash.gt.0) then

c        check for a bad r.

         rfrobn = dlantr('frobenius norm','upper','non-unit diagonal',n,
     *            n,r,ldr,w)
         call scond (n,r,ldr+1,drmax,drmin)
         cond = sdiv (drmax,drmin,overfl)

         if (cond.gt.rcndbd .or. rfrobn.gt.rootn*growth*drmax) then

c           refactorize the hessian and bound the condition estimator.

            call nprset(unitq,n,nfree,nz,ldq,ldr,iw(liperm),iw(lkx),
     *                  w(lgq),r,w(lq),w(lwrk1),w(lres0))
         end if
      end if

c     check the gradients at a feasible x.

      lvrfyc = lverfy
      if (lverfy.ge.10) lvrfyc = -1

      call npchkd(info,msgnp,nstate,lvlder,nfun,ngrad,ldcj,ldcju,n,
     *            ncnln,confun,objfun,iw(lneedc),bigbnd,epsrf,cdint,
     *            fdint,fdchk,fdnorm,objf,xnorm,bl,bu,c,w(lwrk3),
     *            w(lcjac),cjacu,w(lcjdx),w(ldx),w(lgrad),gradu,
     *            w(lhfrwd),w(lhctrl),x,w(lwrk1),w(lwrk2),w,lenw,iuser,
     *            user)

      if (info.ne.0) then
         if (info.gt.0) inform = 7
         if (info.lt.0) inform = info
         go to 80
      end if

      call dcopy(n,w(lgrad),1,w(lgq),1)
      call cmqmul(6,n,nz,nfree,ldq,unitq,iw(lkx),w(lgq),w(lq),w(lwrk1))
c debug 691
      iuser(2) = 1

c     solve the problem.

      if (ncnln.eq.0) then

c        the problem has only linear constraints and bounds.

         call npcore(named,names,unitq,inform,iter,n,nclin,ncnln,nctotl,
     *               nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,nfun,ngrad,
     *               istate,iw(lkactv),iw(lkx),objf,fdnorm,xnorm,confun,
     *               objfun,a,w(lax),bl,bu,c,w(lcjac),cjacu,clamda,
     *               w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw,iuser,user)
      else

c        the problem has some nonlinear constraints.

         if (nclin.gt.0) call smcopy('general',nclin,n,a,lda,w(laqp),
     *                               ldaqp)

c        try and add some nonlinear constraint indices to kactiv.

         call npcrsh(cold,n,nclin,ncnln,nctotl,nactiv,nfree,nz,istate,
     *               iw(lkactv),bigbnd,tolact,bl,bu,c)

         call npcore(named,names,unitq,inform,iter,n,nclin,ncnln,nctotl,
     *               nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,nfun,ngrad,
     *               istate,iw(lkactv),iw(lkx),objf,fdnorm,xnorm,confun,
     *               objfun,w(laqp),w(lax),bl,bu,c,w(lcjac),cjacu,
     *               clamda,w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw,
     *               iuser,user)
c
      end if


c     if required, form the triangular factor of the hessian.

c     first,  form the square matrix  r  such that  h = r'r.
c     compute the  qr  factorization of  r.

      if (lformh.gt.0) then
         lv = lwrk2
         do 60 j = 1, n
            if (j.gt.1) call sload(j-1,zero,w(lv),1)
            lvj = lv + j - 1
            call dcopy(n-j+1,r(j,j),ldr,w(lvj),1)
            call cmqmul(3,n,nz,nfree,ldq,unitq,iw(lkx),w(lv),w(lq),
     *                  w(lwrk1))
            call dcopy(n,w(lv),1,r(j,1),ldr)
   60    continue

         call sgeqr (n,n,r,ldr,w(lwrk1),info)

      end if

c     recover the optional parameters set by the user.

80    call icopy (mxparm,ipsvls,1,iprmls,1)
      call dcopy(mxparm,rpsvls,1,rprmls,1)
      call icopy (mxparm,ipsvnp,1,iprmnp,1)
      call dcopy(mxparm,rpsvnp,1,rprmnp,1)

      if (inform.lt.9) then
         if (ncnln.gt.0) call smcopy('general',ncnln,n,w(lcjac),ldcj,
     *                               cjacu,ldcju)
         call dcopy(n,w(lgrad),1,gradu,1)
      end if
c                                 the diagnositics are not so hot,
c                                 here let the calling routine decide
c                                 what to do. this could be checked 
c                                 again
      ifail = 0
c                                 end of nlpsol
      end

      subroutine lsmove(hitcon,hitlow,linobj,unitgz,nclin,nrank,nrz,n,
     *                  ldr,jadd,numinf,alfa,ctp,ctx,xnorm,ap,ax,bl,bu,
     *                  gq,hz,p,res,r,x,work)
c----------------------------------------------------------------------
c     lsmove  changes x to x + alfa*p and updates ctx, ax, res and gq
c     accordingly.
c     if a bound was added to the working set,  move x exactly on to it,
c     except when a negative step was taken (cmalf may have had to move
c     to some other closer constraint.)
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  alfa, ctp, ctx, xnorm
      integer           jadd, ldr, n, nclin, nrank, nrz, numinf
      logical           hitcon, hitlow, linobj, unitgz

      double precision  ap(*), ax(*), bl(*), bu(*), gq(*), hz(*), p(n),
     *                  r(ldr,*), res(*), work(*), x(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  bnd
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dtrmv

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
      call daxpy(n,alfa,p,1,x,1)
      if (linobj) ctx = ctx + alfa*ctp

      if (hitcon .and. jadd.le.n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa.ge.zero) x(jadd) = bnd
      end if
      xnorm = dnrm2(n,x,1)

      if (nclin.gt.0) call daxpy(nclin,alfa,ap,1,ax,1)

      if (nrz.le.nrank) then
         if (unitgz) then
            res(nrz) = res(nrz) - alfa*hz(nrz)
         else
            call daxpy(nrz,(-alfa),hz,1,res,1)
         end if

         if (numinf.eq.0) then

c           update the transformed gradient gq so that
c           gq = gq + alfa*r'( hz ).
c                            ( 0  )

            if (unitgz) then
               call daxpy(n-nrz+1,alfa*hz(nrz),r(nrz,nrz),ldr,gq(nrz),1)
            else
               call dcopy(nrz,hz,1,work,1)
               call dtrmv('u','t','n',nrz,r,ldr,work,1)
               if (nrz.lt.n) call dgemv('t',nrz,n-nrz,one,r(1,nrz+1),
     *                                  ldr,hz,1,zero,work(nrz+1),1)

               call daxpy(n,alfa,work,1,gq,1)
            end if
         end if
      end if
c                                 end of lsmove
      end

      subroutine chcore(debug,done,first,epsa,epsr,fx,inform,iter,itmax,
     *                  cdest,fdest,sdest,errbnd,f1,f2,h,hopt,hphi)
c----------------------------------------------------------------------
c     chcore  implements algorithm  fd, the method described in
c     gill, p.e., murray, w., saunders, m.a., and wright, m. h.,
c     computing forward-difference intervals for numerical optimization,
c     siam journal on scientific and statistical computing, vol. 4,
c     pp. 310-321, june 1983.

c     the procedure is based on finding an interval (hphi) that
c     produces an acceptable estimate of the second derivative, and
c     then using that estimate to compute an interval that should
c     produce a reasonable forward-difference approximation.

c     one-sided difference estimates are used to ensure feasibility with
c     respect to an upper or lower bound on x. if x is close to an upper
c     bound, the trial intervals will be negative. the final interval is
c     always positive.

c     chcore has been designed to use a reverse communication
c     control structure, i.e., all evaluations of the function occur
c     outside this routine. the calling routine repeatedly calls  chcore
c     after computing the indicated function values.
c----------------------------------------------------------------------
      double precision  bndlo, bndup
      parameter         (bndlo=1.0d-3,bndup=1.0d-1)
      double precision  zero, sixth, fourth
      parameter         (zero=0.0d+0,sixth=1.6d-1,fourth=2.5d-1)
      double precision  half, two
      parameter         (half=5.0d-1,two=2.0d+0)
      double precision  three, four, ten
      parameter         (three=3.0d+0,four=4.0d+0,ten=1.0d+1)

      double precision  cdest, epsa, epsr, errbnd, f1, f2, fdest, fx, h,
     *                  hopt, hphi, sdest
      integer           inform, iter, itmax
      logical           debug, done, first
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout

      double precision  afdmin, cdsave, err1, err2, fdcerr, fdest2,
     *                  fdsave, hsave, oldcd, oldh, oldsd, rho, sdcerr,
     *                  sdsave
      logical           ce1big, ce2big, overfl, te2big

      character*80      rec(6)
c     .. external functions ..
      double precision  sdiv 
      external          sdiv 
c     .. external subroutines ..
      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      save              cdsave, fdsave, hsave, oldh, rho, sdsave,
     *                  ce1big, ce2big, te2big
c----------------------------------------------------------------------
c     bndlo, bndup, and rho control the logic of the routine.
c     bndlo and bndup are the lower and upper bounds that define an
c     acceptable value of the bound on the relative condition error in
c     the second derivative estimate.

c     the scalar rho is the factor by which the interval is multiplied
c     or divided, and also the multiple of the well-scaled interval
c     that is used as the initial trial interval.

      iter = iter + 1

c     compute the forward-,  backward-,  central-  and second-order
c     difference estimates.

      fdest = sdiv (f1-fx,h,overfl)
      fdest2 = sdiv (f2-fx,two*h,overfl)

      oldcd = cdest
      cdest = sdiv (four*f1-three*fx-f2,two*h,overfl)

      oldsd = sdest
      sdest = sdiv (fx-two*f1+f2,h*h,overfl)

c     compute  fdcerr  and  sdcerr,  bounds on the relative condition
c     errors in the first and second derivative estimates.

      afdmin = min(abs(fdest),abs(fdest2))
      fdcerr = sdiv (epsa,half*abs(h)*afdmin,overfl)
      sdcerr = sdiv (epsa,fourth*abs(sdest)*h*h,overfl)

c     select the correct case.

      if (first) then

c        first time through.
c        check whether sdcerr lies in the acceptable range.
c        ------------------------------------------------------------
         first = .false.
         done = sdcerr .ge. bndlo .and. sdcerr .le. bndup
         te2big = sdcerr .lt. bndlo
         ce2big = sdcerr .gt. bndup
         ce1big = fdcerr .gt. bndup

         if ( .not. ce1big) then
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

c           sdcerr is too large.  increase the trial interval.

            oldh = h
            h = h*rho
         end if
      else if (ce2big) then

c        during the last iteration,  the trial interval was
c        increased in order to decrease sdcerr.

         if (ce1big .and. fdcerr.le.bndup) then
            ce1big = .false.
            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if

c        if sdcerr is small enough, accept h.  otherwise,
c        increase h again.

         done = sdcerr .le. bndup
         if ( .not. done) then
            oldh = h
            h = h*rho
         end if
      else if (te2big) then

c        during the last iteration,  the interval was decreased in order
c        to reduce the truncation error.

         done = sdcerr .gt. bndup
         if (done) then

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

c           check whether sdcerr is in range.

            done = sdcerr .ge. bndlo

            if ( .not. done) then

c              sdcerr is still too small, decrease h again.

               oldh = h
               h = h/rho
            end if
         end if

      end if

c     we have either finished or have a new estimate of h.

      if (done) then

c        sufficiently good second-derivative estimate found.
c        compute the optimal interval.

         hphi = abs(h)
         hopt = two*sqrt(epsa)/sqrt(abs(sdest))

c        err1 is the error bound on the forward-difference estimate
c        with the final value of h.  err2 is the difference of fdest
c        and the central-difference estimate with hphi.

         err1 = hopt*abs(sdest)
         err2 = abs(fdest-cdest)
         errbnd = max(err1,err2)

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

c              fdcerr was never small.  probably a constant function.

               inform = 1
               hphi = hopt
               fdest = zero
               cdest = zero
               sdest = zero
               errbnd = zero
            else if (ce2big) then

c              fdcerr was small,  but sdcerr was never small.
c              probably a linear or odd function.

               inform = 2
               hphi = abs(hsave)
               hopt = hphi
               fdest = fdsave
               cdest = cdsave
               sdest = zero
               errbnd = two*epsa/hopt
            else

c              the only remaining case occurs when the second
c              derivative is changing too rapidly for an adequate
c              interval to be found (sdcerr remained small even
c              though h was decreased itmax times).

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
c                                 end of chcore
      end

      subroutine cmalf1(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)
c----------------------------------------------------------------------
c     cmalf1  finds steps palfa1, palfa2 such that
c        x + palfa1*p  reaches a linear constraint that is currently not
c                      in the working set but is satisfied.
c        x + palfa2*p  reaches a linear constraint that is currently not
c                      in the working set but is violated.
c     the constraints are perturbed by an amount featol, so that palfa1
c     is slightly larger than it should be,  and palfa2 is slightly
c     smaller than it should be.  this gives some leeway later when the
c     exact steps are computed by cmalf.
c
c     constraints in the working set are ignored  (istate(j) .ge. 1).
c
c     if negstp is true, the search direction will be taken to be  - p.

c     values of istate(j)....

c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     the values  -2  and  -1  do not occur once a feasible point has
c     been found.

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  bigalf, bigbnd, palfa1, palfa2, pnorm
      integer           jadd1, jadd2, n, nctotl
      logical           firstv, negstp

      double precision  anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer           istate(nctotl)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)

      double precision  absatp, atp, atx, res, rownrm
      integer           i, j, js
      logical           lastv

      character*120     rec(3)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      lastv = .not. firstv
      jadd1 = 0
      jadd2 = 0
      palfa1 = bigalf

      palfa2 = zero
      if (firstv) palfa2 = bigalf

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

            if (abs(atp).le.epspt9*rownrm*pnorm) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

               res = -one
            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

c              first test for smaller palfa1.

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

               if (js.eq.-1) then

c                 the upper bound is violated.  test for either larger
c                 or smaller palfa2, depending on the value of firstv.

                  res = atx - bu(j) - featol(j)
                  if (bigalf*absatp.gt.abs(res)) then
                     if (firstv .and. palfa2*absatp.gt.res .or.
     *                   lastv .and. palfa2*absatp.lt.res) then
                        palfa2 = res/absatp
                        jadd2 = j
                     end if
                  end if
               end if
            else if (atp.gt.zero .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.
c              test for smaller palfa1.

               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx + featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (palfa1*atp.gt.res) then
                        palfa1 = res/atp
                        jadd1 = j
                     end if
                  end if
               end if

               if (js.eq.-2) then

c                 the lower bound is violated.  test for a new palfa2.

                  res = bl(j) - atx - featol(j)
                  if (bigalf*atp.gt.abs(res)) then
                     if (firstv .and. palfa2*atp.gt.res .or. lastv .and.
     *                   palfa2*atp.lt.res) then
                        palfa2 = res/atp
                        jadd2 = j
                     end if
                  end if
               end if
            end if

         end if
   20 continue
c                                 end of cmalf1
      end

      subroutine cmr1md(n,nu,nrank,ldr,lenv,lenw,r,u,v,w,c,s)
c----------------------------------------------------------------------
c     cmr1md  modifies the  nrank*n  upper-triangular matrix  r  so that
c     q*(r + v*w')  is upper triangular,  where  q  is orthogonal,
c     v  and  w  are vectors, and the modified  r  overwrites the old.
c     q  is the product of two sweeps of plane rotations (not stored).
c     if required,  the rotations are applied to the nu columns of
c     the matrix  u.
c
c     the matrix v*w' is an (lenv) by (lenw) matrix.
c     the vector v is overwritten.

      integer           ldr, lenv, lenw, n, nrank, nu

      double precision  c(n), r(ldr,*), s(n), u(n,*), v(n), w(n)

      integer           j
c     .. external subroutines ..
      external          daxpy, ssrotg, susqr , sutsrs, sgesrc
c----------------------------------------------------------------------
      j = min(lenv,nrank)
      if (nrank.gt.0) then

c        reduce  v to beta*e( j )  using a backward sweep of rotations
c        in planes (j-1, j), (j-2, j), ..., (1, j).

         call ssrotg('fixed','backwards',j-1,v(j),v,1,c,s)

c        apply the sequence of rotations to u.

         if (nu.gt.0) call sgesrc('left','bottom','backwards',j,nu,1,j,
     *                            c,s,u,n)

c        apply the sequence of rotations to r. this generates a spike in
c        the j-th row of r, which is stored in s.

         call sutsrs('left',n,1,j,c,s,r,ldr)

c        form  beta*e(j)*w' + r.  this a spiked matrix, with a row
c        spike in row j.

         call daxpy(min(j-1,lenw),v(j),w,1,s,1)
         call daxpy(lenw-j+1,v(j),w(j),1,r(j,j),ldr)
c

c        eliminate the spike using a forward sweep of rotations in
c        planes (1, j), (2, j), ..., (j-1, j).

         call susqr ('left',n,1,j,c,s,r,ldr)

c        apply the rotations to u.

         if (nu.gt.0) call sgesrc('left','bottom','forwards',j,nu,1,j,c,
     *                            s,u,n)
      end if
c                                 end of cmr1md
      end

      subroutine lpcolr(nrz,ldr,r,rzz)
c----------------------------------------------------------------------
c     lpcolr  loads the last column of the  nrz x nrz  triangular factor
c     rz  with the multiple  rzz  of the  nrz-th unit vector.
c----------------------------------------------------------------------
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  rzz
      integer           ldr, nrz

      double precision  r(ldr,*)
c     .. external subroutines ..
      external          sload
c----------------------------------------------------------------------
      if (nrz.eq.0) return

      call sload(nrz-1,zero,r(1,nrz),1)
      r(nrz,nrz) = rzz
c                                 end of lpcolr
      end

      subroutine npalf(inform,n,nclin,ncnln,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,anorm,adx,ax,bl,bu,dslk,dx,slk,x)
c----------------------------------------------------------------------
c     npalf  finds a step alfa such that the point x + alfa*p reaches
c     one of the slacks or linear constraints.  the step alfa is the
c     maximum step that can be taken without violating one of the slacks
c     or linear constraints that is currently satisfied.
c----------------------------------------------------------------------
      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  alfa, alfmax, alfmin, bigbnd, dxnorm
      integer           inform, n, nclin, ncnln

      double precision  adx(*), anorm(*), ax(*), bl(*), bu(*), dslk(*),
     *                  dx(n), slk(*), x(n)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)

      double precision  adxi, axi, res, rownrm
      integer           i, j

      character*80      rec(4)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      alfa = alfmax
      j = 1

c     +    while (j .le. n+nclin+ncnln .and. alfa .gt. alfmin) do
   20 if (j.le.n+nclin+ncnln .and. alfa.gt.alfmin) then

         if (j.le.n) then
            axi = x(j)
            adxi = dx(j)
            rownrm = one
         else if (j.le.n+nclin) then

            i = j - n
            axi = ax(i)
            adxi = adx(i)
            rownrm = anorm(i) + one
         else

            i = j - n - nclin
            axi = slk(i)
            adxi = dslk(i)
            rownrm = one
         end if

         res = -one
         if (adxi.le.-epspt9*rownrm*dxnorm) then

c           constraint decreasing.

            adxi = -adxi
            if (bl(j).gt.-bigbnd) res = axi - bl(j)
         else if (adxi.gt.epspt9*rownrm*dxnorm) then

c           constraint increasing.

            if (bu(j).lt.bigbnd) res = bu(j) - axi
         end if

         if (res.gt.zero .and. alfa*adxi.gt.res) alfa = res/adxi

         j = j + 1
         go to 20
c        +    end while
      end if

c     determine alfa, the bound on the step to be taken.

      alfa = max(alfa,alfmin)

      inform = 0
      if (alfa.ge.alfmax) inform = 1
c                                 end of npalf
      end

      subroutine cmtsol(mode,nrowt,n,t,y)
c----------------------------------------------------------------------
c     cmtsol  solves equations involving a reverse-triangular matrix  t
c     and a right-hand-side vector  y,  returning the solution in  y.
c----------------------------------------------------------------------
      double precision  zero
      parameter         (zero=0.0d+0)

      integer           mode, n, nrowt

      double precision  t(nrowt,*), y(n)

      double precision  yj
      integer           j, jj, l, n1
c     .. external subroutines ..
      external          daxpy
c----------------------------------------------------------------------
      n1 = n + 1
      if (mode.eq.1) then

c        mode = 1  ---  solve  t * y(new) = y(old).

         do 20 j = 1, n
            jj = n1 - j
            yj = y(j)/t(j,jj)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.zero) call daxpy(l,(-yj),t(j+1,jj),1,
     *          y(j+1),1)
   20    continue
      else

c        mode = 2  ---  solve  t' y(new) = y(old).

         do 40 j = 1, n
            jj = n1 - j
            yj = y(j)/t(jj,j)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.zero) call daxpy(l,(-yj),t(jj,j+1),
     *          nrowt,y(j+1),1)
   40    continue
      end if

c     reverse the solution vector.

      if (n.gt.1) then
         l = n/2
         do 60 j = 1, l
            jj = n1 - j
            yj = y(j)
            y(j) = y(jj)
            y(jj) = yj
   60    continue
      end if
c                                 end of cmtsol
      end

      subroutine nggnbu(n,nu,nrank,ldr,i,j,r,u,c,s)
c----------------------------------------------------------------------
c     nggnbu  interchanges the  i-th  and  j-th  (i .lt. j)  columns of
c     an  nrank*n  upper-trapezoidal matrix  r   and restores the
c     resulting matrix to upper-trapezoidal form using two sweeps of
c     plane rotations applied on the left.  r is overwritten.

c     if nu .gt. 0,  the rotations are applied to the  nu  columns of
c     the matrix  u.
c----------------------------------------------------------------------
      double precision  zero
      parameter         (zero=0.0d+0)

      integer           i, j, ldr, n, nrank, nu

      double precision  c(n), r(ldr,*), s(n), u(n,*)

      integer           lenj
c     .. external subroutines ..
      external          sload, ssrotg, susqr , sutsrs, sgesrc, dswap
c----------------------------------------------------------------------
c     swap the elements of the i-th and j-th columns of r on, or above,
c     the main diagonal.

      call dswap(min(i,nrank),r(1,i),1,r(1,j),1)
      lenj = min(j,nrank)

      if (lenj.gt.i) then

c        reduce elements  r(i+1,j), ..., r(lenj,j)  to  beta*e(lenj)
c        using a backward sweep in planes
c        (lenj-1,lenj), (lenj-2,lenj), ..., (i+1,lenj).
c        if required, apply the sequence of rotations to u.

         call ssrotg('fixed','backwards',lenj-i-1,r(lenj,j),r(i+1,j),1,
     *               c(i+1),s(i+1))

         if (nu.gt.0) call sgesrc('left','bottom','backwards',n,nu,i+1,
     *                            lenj,c,s,u,n)

c        put zeros into the j-th column of r in positions corresponding
c        to the sub-diagonals of the i-th column.

         s(i) = r(lenj,j)
         call sload(lenj-i,zero,r(i+1,j),1)

c        apply the sequence of rotations to r.  this generates a spike
c        in the lenj-th row of r, which is stored in s.

         call sutsrs('left',n,i+1,lenj,c,s,r,ldr)

c        eliminate the spike using a forward sweep in planes
c        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
c        if necessary, apply the sequence of rotations to u.

         call susqr ('left',n,i,lenj,c,s,r,ldr)
c
         if (nu.gt.0) call sgesrc('left','bottom','forwards',lenj,nu,i,
     *                            lenj,c,s,u,n)
      end if
c                                 end of nggnbu
      end

      subroutine lsdel(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,res,r,t,gq,
     *                  zy,c,s)
c----------------------------------------------------------------------
c     lsdel  updates the least-squares factor r and the factorization
c     a(free) (z y) = (0 t) when a regular, temporary or artificial
c     constraint is deleted from the working set.
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      integer           jdel, kdel, lda, ldr, ldt, ldzy, n, nactiv,
     *                  nfree, ngq, nrank, nres, nrz, nz
      logical           unitq

      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), zy(ldzy,*)
      integer           kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  cs, sn
      integer           i, ir, itdel, jart, k, ka, ld, npiv, nrz1, nsup,
     *                  nt

      character*80      rec(4)
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. external subroutines ..
      external          dcopy, dswap, nggnbu, srotgc, sload, scond ,
     *                  sutsqr, sgesrc, nggqzz

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg008/asize, dtmax, dtmin
c----------------------------------------------------------------------
      if (jdel.gt.0) then

c        regular constraint or temporary bound deleted.

         if (jdel.le.n) then

c           case 1.  a simple bound has been deleted.
c           =======  columns nfree+1 and ir of r must be swapped.

            ir = nz + kdel

            itdel = 1
            nfree = nfree + 1
c
            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,nfree,ir,r,
     *                                     res,c,s)
               call dswap(ngq,gq(nfree,1),n,gq(ir,1),n)
            end if

            if ( .not. unitq) then

c              copy the incoming column of  a(free)  into the end of t.

               do 20 ka = 1, nactiv
                  i = kactiv(ka)
                  t(ka,nfree) = a(i,jdel)
   20          continue

c              expand q by adding a unit row and column.

               if (nfree.gt.1) then
                  call sload(nfree-1,zero,zy(nfree,1),ldzy)
                  call sload(nfree-1,zero,zy(1,nfree),1)
               end if
               zy(nfree,nfree) = one
            end if
         else

c           case 2.  a general constraint has been deleted.
c           =======

            itdel = kdel
            nactiv = nactiv - 1

c           delete row  kdel  of t and move up the ones below it.
c           t becomes reverse lower hessenberg.

            do 40 i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld = nfree - i
               call dcopy(i+1,t(i+1,ld),ldt,t(i,ld),ldt)
   40       continue
         end if

         nz = nz + 1

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

            if (nsup.gt.0) then
               npiv = nfree - itdel + 1
               if (nsup.gt.1) then
                  call dcopy(nsup-1,t(nactiv-1,nz+1),ldt-1,s(nz+1),1)
                  call nggqzz('remove',nactiv,1,nsup,c(nz+1),s(nz+1),
     *                        t(1,nz+1),ldt)
               end if

               call srotgc(t(nactiv,nz+1),t(nactiv,nz),cs,sn)
               t(nactiv,nz) = zero
               s(nz) = -sn
               c(nz) = cs

               call sgesrc('right','variable','backwards',nfree,nfree,
     *                     nz,npiv,c,s,zy,ldzy)
               call sgesrc('left ','variable','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gq,n)

               nt = min(nrank,npiv)

               if (nt.lt.npiv .and. nt.gt.0) then

c                 r is upper trapezoidal, pretend r is (nt x n) and
c                 apply the rotations in columns  max(nt,nz)  thru npiv.

                  call sgesrc('right','variable','backwards',nt,n,
     *                        max(nt,nz),npiv,c,s,r,ldr)
               end if

c              apply the column transformations to the triangular part
c              of r.  a subdiagonal element is generated that must be
c              eliminated by a row rotation before the next column
c              transformation can be applied.

               if (nz.lt.nt) then
                  call sutsqr('right',nt,nz,nt,c,s,r,ldr)
               end if

c              apply the row rotations to the remaining rows of r.

               call sgesrc('left','variable','backwards',nt,n-nt,nz,nt,
     *                     c,s,r(1,min(nt+1,n)),ldr)

               if (nres.gt.0) call sgesrc('left','variable','backwards',
     *                                    nt,nres,nz,nt,c,s,res,n)

            end if
            call scond (nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
         end if
      end if

      nrz1 = nrz + 1

      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamax(nz-nrz1+1,gq(nrz1,1),1)
         else
            jart = -jdel
         end if

         if (jart.gt.nrz1) then

c           swap columns nrz1 and jart of r.

            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,zy(1,nrz1),1,zy(1,jart),1)
            end if

            call dswap(ngq,gq(nrz1,1),n,gq(jart,1),n)
            if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,nrz1,jart,r,
     *                                  res,c,s)
         end if
      end if

      nrz = nrz1
c                                 end of lsdel
      end

      subroutine npsetx(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,bl,bu,rpq,rpq0,dx,gq,r,t,zy,work)
c----------------------------------------------------------------------
c     npsetx   defines a point which lies on the initial working set for
c     the qp subproblem.  this routine is similar to lssetx except
c     that advantage is taken of the fact that the initial estimate of
c     the solution of the least-squares subproblem is zero.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  dxnorm, gdx
      integer           ldaqp, ldr, ldt, ldzy, n, nactiv, ncqp, nctotl,
     *                  nfree, nlnx, nz
      logical           unitq

      double precision  adx(*), aqp(ldaqp,*), bl(nctotl), bu(nctotl),
     *                  dx(n), gq(n), r(ldr,*), rpq(nlnx), rpq0(nlnx),
     *                  t(ldt,*), work(n), zy(ldzy,*)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  bnd
      integer           i, j, k, nfixed, nr

      character*80      rec(2)
c     .. external functions ..
      double precision  ddot, dnrm2
      external          ddot, dnrm2
c     .. external subroutines ..
      external          dcopy, dgemv, dscal, dtrmv, cmtsol, cmqmul,
     *                  sload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
      nfixed = n - nfree

      gdx = zero

      call sload(n,zero,dx,1)
      call sload(nlnx,zero,rpq,1)
      call sload(nlnx,zero,rpq0,1)

      if (nactiv+nfixed.gt.0) then

c        set  work = residuals for constraints in the working set.
c        solve for  dx,  the smallest correction to  x  that gives a
c        point on the constraints in the working set.
c        set the fixed variables on their bounds,  solve the triangular
c        system  t*(dxy) = residuals,  and define  dx = y*(dxy).
c        use  (dxy)  to update  d(=pr)  as  d = d - r'(  0  ).
c                                                     ( dxy )
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
            work(nz+i) = bnd - ddot(n,aqp(k,1),ldaqp,dx,1)
   40    continue
c
         if (nactiv.gt.0) call cmtsol(1,ldt,nactiv,t(1,nz+1),work(nz+1))
         call dcopy(nactiv+nfixed,work(nz+1),1,dx(nz+1),1)
         if (nz.gt.0) call sload(nz,zero,dx,1)
c
         gdx = ddot(nactiv+nfixed,gq(nz+1),1,dx(nz+1),1)
c
         if (nz.lt.n) then
            call dgemv('n',nz,n-nz,-one,r(1,nz+1),ldr,dx(nz+1),1,one,
     *                 rpq,1)
            if (nz.lt.nlnx) then
               nr = ldr
               if (nz+1.eq.n) nr = 1
               call dcopy(nlnx-nz,dx(nz+1),1,rpq(nz+1),1)
               call dscal(nlnx-nz,(-one),rpq(nz+1),1)
               call dtrmv('u','n','n',nlnx-nz,r(nz+1,nz+1),nr,rpq(nz+1),
     *                    1)
               if (nlnx.lt.n) then
                  nr = ldr
                  if (nlnx+1.eq.n) nr = n - nz
                  call dgemv('n',nlnx-nz,n-nlnx,-one,r(nz+1,nlnx+1),nr,
     *                       dx(nlnx+1),1,one,rpq(nz+1),1)
               end if
            end if
         end if

         call cmqmul(2,n,nz,nfree,ldzy,unitq,kx,dx,zy,work)
      end if

c     compute the 2-norm of  dx.
c     initialize  a*dx.

      dxnorm = dnrm2(n,dx,1)
      if (ncqp.gt.0) call dgemv('n',ncqp,n,one,aqp,ldaqp,dx,1,zero,adx,
     *                          1)
c                                 end of npsetx
      end

      subroutine lsmuls(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,
     *                  nrz,istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,
     *                  jtiny,jbigst,kbigst,trulam,a,anorms,gq,rlamda,t,
     *                  wtinf)
c----------------------------------------------------------------------
c     lsmuls  first computes the lagrange multiplier estimates for the
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

      integer           ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)

      double precision  dinky, trulam
      integer           jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, msglvl, n, nactiv, nfree, nrz, numinf,
     *                  nz
      character*2       prbtyp

      double precision  a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  anormj, biggst, blam, rlam, scdlam, smllst,
     *                  tinylm
      integer           i, is, j, k, l, nfixed

      character*80      rec(80)
c     .. external subroutines ..
      external          dcopy, cmtsol

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
      nfixed = n - nfree

      jsmlst = 0
      ksmlst = 0
      smllst = -dinky

      tinylm = dinky
      jtiny = 0

      jbigst = 0
      kbigst = 0
      biggst = one + dinky

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

      end if


c     compute jsmlst for regular constraints and temporary bounds.
c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.

      if (n.gt.nz) call dcopy(n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call cmtsol(2,ldt,nactiv,t(1,nz+1),rlamda)
c     set elements nactiv, nactiv+1,... of  rlamda  equal to
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
         if (numinf.gt.0 .and. j.gt.jinf) then
            scdlam = rlam/wtinf(j)
            if (scdlam.gt.biggst) then
               biggst = scdlam
               trulam = rlamda(k)
               jbigst = j
               kbigst = k
            end if
         end if
  100 continue

c                                 end of lsmuls

      end

      subroutine cmfeas(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,bl,
     *                  bu,featol,x)


c     cmfeas  checks the residuals of the constraints that are believed
c     to be feasible.  the number of constraints violated by more than
c     featol is computed, along with the maximum constraint violation.

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd, errmax
      integer           jmax, n, nclin, nviol

      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)

      double precision  biglow, bigupp, con, feasj, res
      integer           is, j

      character*80      rec(2)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd


c     compute the number of constraints (nviol) violated by more than
c     featol and  the maximum constraint violation (errmax).
c     (the residual of a constraint in the working set is treated as if
c     it were an equality constraint fixed at that bound.)

      nviol = 0
      jmax = 0
      errmax = zero
c
      do 40 j = 1, n + nclin
         is = istate(j)
c
         if (is.ge.0) then
            feasj = featol(j)
c
            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if
c
c           check for constraint violations.
c
            if (bl(j).gt.biglow) then
               res = bl(j) - con
               if (res.gt.feasj) then
                  nviol = nviol + 1
                  go to 20
               end if
            end if
c
            if (bu(j).lt.bigupp) then
               res = bu(j) - con
               if (res.lt.(-feasj)) then
                  nviol = nviol + 1
                  res = -res
                  go to 20
               end if
            end if
c
c           this constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.
c
            res = zero
c
            if (is.eq.1) then
               res = abs(bl(j)-con)
c
            else if (is.eq.2) then
               res = abs(bu(j)-con)
c
            else if (is.eq.3) then
               res = abs(bu(j)-con)
            end if
c
            if (res.gt.feasj) nviol = nviol + 1
c
   20       if (res.gt.errmax) then
               jmax = j
               errmax = res
            end if
         end if
   40 continue

c                                 end of cmfeas

      end

      subroutine srchc (first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *                  toltny,alfa,alfbst,fbest,gbest)
c----------------------------------------------------------------------
c     srchc   finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     srchc   requires both  f(alpha)  and  f'(alpha) to be evaluated at
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
c                   it is subsequently altered by srchc .
c
c     debug         specifies whether detailed output is wanted.
c
c     maxf          is an upper limit on the number of times srchc  is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).
c
c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by srchc  (see below).
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
c                   these values may be reduced by srchc .
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
c                   for the new trial alfa, and re-enter srchc .
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
c                   inform = 5 is never set by srchc .
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
c     numf          counts the number of times srchc  has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).
c
c     alfa          is the point at which the next function ftry and
c                   derivative gtry must be computed.
c
c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever srchc  returns
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

      double precision  zero, point1, half
      parameter         (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision  one, three, five
      parameter         (one=1.0d+0,three=3.0d+0,five=5.0d+0)
      double precision  ten, eleven
      parameter         (ten=1.0d+1,eleven=1.1d+1)

      double precision  alfa, alfbst, alfmax, epsaf, fbest, ftry, g0,
     *                  gbest, gtry, targtg, tolabs, tolrel, toltny
      integer           inform, maxf, nout, numf
      logical           debug, done, first, imprvd

      double precision  a, absr, artifa, artifb, b, daux, dtry, factor,
     *                  fw, gw, q, r, s, scale, tol, tolmax, truea,
     *                  trueb, xmidpt, xtry, xw
      integer           nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, fitok,
     *                  found, moved, quitf, quiti, setxw, wset

      character*120     rec(7)

      save              braktd, crampd, extrap, moved, wset, nsamea,
     *                  nsameb, a, b, factor, xtry, xw, fw, gw, tolmax
c----------------------------------------------------------------------
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
         badfun = alfmax .le. toltny .or. g0 .ge. zero
         done = badfun
         moved = .false.
c
         if ( .not. done) then
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
         end if
      else

c        subsequent entries. the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry and gtry.

         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1
c
         if ( .not. braktd) then
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
            extrap = xw .lt. zero .and. gbest .lt. zero .or. xw .gt.
     *               zero .and. gbest .gt. zero
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
               setxw = ftry .lt. fw .or. .not. extrap
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
         if (quiti .and. .not. moved) then
c
c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.
c
            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if
c
         done = quitf .or. quiti .or. found


c        proceed with the computation of a trial steplength.
c        the choices are...
c        1. parabolic fit using derivatives only, if the f values are
c           close.
c        2. cubic fit for a minimizer, using both f and f'.
c        3. damped cubic or parabolic fit if the regular fit appears to
c           be consistently overestimating the distance to a minimizer.
c        4. bisection, geometric bisection, or a step of  tol  if
c           choices 2 or 3 are unsatisfactory.

         if ( .not. done) then
            xmidpt = half*(a+b)
            s = zero
            q = zero
c
            if (closef) then

c              fit a parabola to the two best gradient values.

               s = gbest
               q = gbest - gw

            else

c              fit cubic through  fbest  and  fw.


               fitok = .true.
               r = three*(fbest-fw)/xw + gbest + gw
               absr = abs(r)
               s = sqrt(abs(gbest))*sqrt(abs(gw))
c
c              compute  q =  the square root of  r*r - gbest*gw.
c              the method avoids unnecessary underflow and overflow.
c
               if ((gw.lt.zero .and. gbest.gt.zero)
     *             .or. (gw.gt.zero .and. gbest.lt.zero)) then
                  scale = absr + s
                  if (scale.eq.zero) then
                     q = zero
                  else
                     q = scale*sqrt((absr/scale)**2+(s/scale)**2)
                  end if
               else if (absr.ge.s) then
                  q = sqrt(absr+s)*sqrt(absr-s)
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
            if ( .not. braktd) then
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

               daux = abs(xw)
               dtry = b - a
               if (daux.ge.dtry) then
                  xtry = five*dtry*(point1+dtry/daux)/eleven
               else
                  xtry = half*sqrt(daux)*sqrt(dtry)
               end if
               if (xw.gt.zero) xtry = -xtry

c              reset the artificial bounds.  if the point computed by
c              extrapolation is rejected,  xtry will remain at the
c              relevant artificial bound.

               if (xtry.le.zero) artifa = xtry
               if (xtry.gt.zero) artifb = xtry
            else
c
c              the points are configured for an interpolation.  the
c              default value xtry bisects the interval of uncertainty.
c              the artificial interval is just (a, b).
c
               xtry = xmidpt

               if (nsamea.ge.3 .or. nsameb.ge.3) then
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
               if (s*xw.ge.q*artifa .and. s*xw.le.q*artifb) then
c
c                 accept the polynomial fit.
c
                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if

               end if
            end if
         end if
      end if
c

c
      if ( .not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then
c
c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)
c
            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
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

c                                 end of srchc

      end

      subroutine npfd (centrl,inform,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,c2,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)
c----------------------------------------------------------------------
c     npfd  evaluates any missing gradients.
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  three, four
      parameter         (three=3.0d+0,four=4.0d+0)

      double precision  bigbnd, cdint, fdint, fdnorm, objf
      integer           inform, ldcj, ldcju, lenw, n, ncnln
      logical           centrl

      double precision  bl(n), bu(n), c(*), c1(*), c2(*), cjac(ldcj,*),
     *                  cjacu(ldcju,*), grad(n), gradu(n), hcntrl(n),
     *                  hforwd(n), user(*), w(lenw), x(n)
      integer           iuser(*), needc(*)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           lfdset, lvldif, ncdiff, nfdiff
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  biglow, bigupp, delta, objf1, objf2, stepbl,
     *                  stepbu, xj
      integer           i, j, mode, ncolj, nstate

      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------

      inform = 0

c     use the pre-assigned difference intervals to approximate the
c     derivatives.

c     use either the same interval for each element (lfdset = 1),
c     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).

      nstate = 0
      mode = 0

      biglow = -bigbnd
      bigupp = bigbnd
c
      fdnorm = zero
c
      do 80 j = 1, n
         xj = x(j)
         ncolj = 0
         if (ncdiff.gt.0) then
            do 20 i = 1, ncnln
               if (cjacu(i,j).eq.rdummy) then
                  needc(i) = 1
                  ncolj = ncolj + 1
               else
                  needc(i) = 0
               end if
   20       continue
         end if
c
         if (ncolj.gt.0 .or. gradu(j).eq.rdummy) then
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

c           if (x(j).lt.0d0.or.x(j).gt.1d0) then 
c              xj = xj*one
c              write (*,*) x(j)
c           end if

            if (ncolj.gt.0) then
               call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,nstate,
     *                     iuser,user)
               if (mode.lt.0) go to 100
            end if
c
            if (gradu(j).eq.rdummy) then
               call objfun(mode,n,x,objf1,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100
            end if
c
            if (centrl) then

c              central differences.

               x(j) = xj + delta + delta

c           if (x(j).lt.0d0.or.x(j).gt.1d0) then
c              xj = xj*one
c              write (*,*) x(j)
c           end if
c
               if (ncolj.gt.0) then
                  call confun(mode,ncnln,n,ldcju,needc,x,c2,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 100
c
                  do 40 i = 1, ncnln
                     if (needc(i).eq.1) cjac(i,j) = (four*c1(i)
     *                   -three*c(i)-c2(i))/(delta+delta)
   40             continue
               end if
c
               if (gradu(j).eq.rdummy) then
                  call objfun(mode,n,x,objf2,gradu,nstate,iuser,user)
                  if (mode.lt.0) go to 100
c
                  grad(j) = (four*objf1-three*objf-objf2)/(delta+delta)
c
               end if
            else

c              forward differences.

               if (ncolj.gt.0) then
                  do 60 i = 1, ncnln
                     if (needc(i).eq.1) cjac(i,j) = (c1(i)-c(i))/delta
   60             continue
               end if
c
               if (gradu(j).eq.rdummy) grad(j) = (objf1-objf)/delta
c
            end if
         end if
         x(j) = xj
c
   80 continue
c
      return
c
  100 inform = mode

c                                 end of npfd

      end

      subroutine lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s,msglvl)
c----------------------------------------------------------------------
c     lsadd   updates the factorization,  a(free) * (z y) = (0 t),  when
c     a constraint is added to the working set.  if  nrank .gt. 0, the
c     factorization  ( r ) = pcq  is also updated,  where  c  is the
c                    ( 0 )
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

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  condmx
      integer           iadd, ifix, inform, jadd, lda, ldr, ldt, ldzy,
     *                  msglvl, n, nactiv, nfree, ngq, nrank, nres, nz
      logical           unitq

      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer           kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin, epspt3, epspt5, epspt8,
     *                  epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  cond, condbd, dtnew, tdtmax, tdtmin
      integer           i, nanew, nfmin, npiv, nt
      logical           bound, overfl

      character*80      rec(5)
c     .. external functions ..
      double precision  dnrm2, sdiv 
      external          dnrm2, sdiv 
c     .. external subroutines ..
      external          dcopy, dscal, cmqmul, scond , ssrotg, smload,
     *                  sgeapr, sutsr1, suhqr, sutsrh, sgesrc, nggqzz

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg008/asize, dtmax, dtmin
c----------------------------------------------------------------------
c
c     if the condition estimator of the updated factors is greater than
c     condbd,  a warning message is printed.
c
      condbd = one/epspt9
c
      overfl = .false.
      bound = jadd .le. n
c
      if (bound) then

c        a simple bound has entered the working set.  iadd  is not used.


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


c
         nanew = nactiv + 1
c
c        transform the incoming row of  a  by  q'.  use c as workspace.
c
         call dcopy(n,a(iadd,1),lda,w,1)
         call cmqmul(8,n,nz,nfree,ldzy,unitq,kx,w,zy,c)
c
c        check that the incoming row is not dependent upon those
c        already in the working set.
c
         dtnew = dnrm2(nz,w,1)
         if (nactiv.eq.0) then
c
c           this is the only general constraint in the working set.
c
            cond = sdiv (asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else
c
c           there are already some general constraints in the working
c           set. update the estimate of the condition number.
c
            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = sdiv (tdtmax,tdtmin,overfl)
         end if
c
         if (cond.gt.condmx .or. overfl) go to 60
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

         if (ngq.gt.0) call sgeapr('left','transpose',nfree-1,w,ngq,gq,
     *                             n)
c
         if (nrank.gt.0) then
c
c           apply the pairwise interchanges to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1), s(2), ..., s(nt-1).
c
            call sutsr1 ('right',n,ifix,nt,s,r,ldr)
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
               call sgeapr('right','normal',nfree-1,w,nt,r,ldr)
            end if
c
c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.
c
            call suhqr('left ',n,ifix,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc('left','variable','forwards',nt,
     *                                 nres,ifix,nt,c,s,res,n)
         end if
      else

c        full matrix q.  define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1,2), (2,3), ...,
c        (npiv-1,npiv).  the rotations must be applied to zy, r, t
c        and gq'.

         call ssrotg('varble','forwrds',npiv-1,w(npiv),w,1,c,s)
c
         if (bound .and. nactiv.gt.0) then
c
            call dcopy(nactiv,s(nz),1,w(nz),1)
c
            s(nz) = s(nz)*t(nactiv,nz+1)
            t(nactiv,nz+1) = c(nz)*t(nactiv,nz+1)
c
            call nggqzz('create',nactiv,1,nactiv,c(nz+1),s(nz+1),
     *                  t(1,nz+1),ldt)
            call dcopy(nactiv,s(nz),1,t(nactiv,nz),ldt-1)
c
            call dcopy(nactiv,w(nz),1,s(nz),1)
         end if
c
         if (ngq.gt.0) call sgesrc('left ','variable','forwards',npiv,
     *                             ngq,1,npiv,c,s,gq,n)
         call sgesrc('right','variable','forwards',nfree,nfree,1,npiv,c,
     *               s,zy,ldzy)
c
         if (nrank.gt.0) then
c
c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nt-1).
c
            nt = min(nrank,npiv)
            call sutsrh('right',n,1,nt,c,s,r,ldr)
c
            if (nt.lt.npiv) then
c
c              r is upper trapezoidal.  pretend r is (nt x n) and
c              apply the rotations in columns  nt  thru  npiv.
c
               call sgesrc('right','variable','forwards',nt,n,nt,npiv,c,
     *                     s,r,ldr)
            end if
c
c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.
c
            call suhqr('left ',n,1,nt,c,s,r,ldr)
            if (nres.gt.0) call sgesrc('left','variable','forwards',nt,
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
               cond = sdiv (tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint.  install the new row of t.

            call dcopy(nanew,w(nz),1,t(nanew,nz),ldt)
         end if
      end if
c

c     prepare to exit.  check the magnitude of the condition estimator.

   60 if (nanew.gt.0) then
         if (cond.lt.condmx .and. .not. overfl) then
c
c           the factorization has been successfully updated.
c
            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin

         else
c
c           the proposed working set appears to be linearly dependent.
c
            inform = 1

         end if
      end if

c                                 end of lsadd

      end

      subroutine chcjac(inform,lvlder,msglvl,ncset,n,ncnln,ldcj,ldcju,
     *                  bigbnd,epsrf,oktol,fdchk,xnorm,confun,needc,bl,
     *                  bu,c,c1,cjac,cjacu,cjdx,dx,err,x,y,iuser,user)
c----------------------------------------------------------------------

c     chcjac  checks if the gradients of the constraints have been coded
c     correctly.
c
c     on input,  the values of the constraints at the point x are stored
c     in c.  their corresponding gradients are stored in cjacu.  if any
c     jacobian element has not been specified,  it will have a dummy
c     value.  missing values are not checked.
c
c     a cheap test is first undertaken by calculating the directional
c     derivative using two different methods.  if this proves
c     satisfactory and no further information is desired, chcjac is
c     terminated. otherwise, chcore is called to give optimal
c     step-sizes and a central-difference approximation to each
c     element of the jacobian for which a test is deemed necessary,
c     either by the program or the user.
c
c     lvrfyc has the following meaning...
c
c     -1        do not do any check.
c     0        do the cheap test only.
c     2 or 3   do both cheap and full test.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, point9
      parameter         (zero=0.0d+0,half=0.5d+0,point9=0.9d+0)
      double precision  one, two, ten
      parameter         (one=1.0d+0,two=2.0d+0,ten=1.0d+1)
      character*4       lbad, lgood
      parameter         (lbad='bad?',lgood='  ok')

      double precision  bigbnd, epsrf, fdchk, oktol, xnorm
      integer           inform, ldcj, ldcju, lvlder, msglvl, n, ncnln,
     *                  ncset

      double precision  bl(n), bu(n), c(*), c1(*), cjac(ldcj,*),
     *                  cjacu(ldcju,*), cjdx(*), dx(n), err(*), user(*),
     *                  x(n), y(n)
      integer           iuser(*), needc(*)
c     .. subroutine arguments ..
      external          confun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, lvrfyc, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg), jverfy(4)

      double precision  biglow, bigupp, cdest, cij, cjdiff, cjsize,
     *                  colmax, dxj, dxmult, emax, epsaci, errbnd, f1,
     *                  f2, fdest, h, hopt, hphi, sdest, signh, stepbl,
     *                  stepbu, xj
      integer           i, imax, info, irow, iter, itmax, j, j3, j4,
     *                  jcol, mode, ncheck, ncolj, ngood, nstate, nwrong
      logical           const, debug, done, first, headng, needed, ok
      character*4       key

      character*18      result(0:4)
      character*120     rec(4)
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, chcore, iload , sload,
     *                  smcopy, smload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg015/lvrfyc, jverfy
      common            /ngg002/inpdbg, npdbg

      data              result/'                 ', 'constant?      ',
     *                  'linear or odd?   ', 'too nonlinear?',
     *                  'small derivative?'/
c----------------------------------------------------------------------
c
      inform = 0
      needed = ncnln .gt. 0 .and. lvrfyc .eq. 0 .or. lvrfyc .eq. 2 .or.
     *         lvrfyc .eq. 3
      if ( .not. needed) return

      debug = npdbg .and. inpdbg(5) .gt. 0
      nstate = 0

      biglow = -bigbnd
      bigupp = bigbnd

c     do the cheap test.

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

c     do not perturb  x(j)  if the  j-th  column contains any
c     unknown elements.  compute the directional derivative for each
c     constraint gradient.

      ncheck = 0
      do 60 j = 1, n
         do 40 i = 1, ncnln
            if (cjac(i,j).eq.rdummy) then
               dx(j) = zero
               go to 60
            end if
   40    continue
         ncheck = ncheck + 1
c
         xj = x(j)
         stepbl = -one
         stepbu = one
         if (bl(j).gt.biglow) stepbl = max(stepbl,bl(j)-xj)
         if (bu(j).lt.bigupp .and. bu(j).gt.bl(j)) stepbu = min(stepbu,
     *       bu(j)-xj)
c
         if (half*(stepbl+stepbu).lt.zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
   60 continue
c
      if (ncheck.eq.0) then

      else
c
c        compute  (jacobian)*dx.
c
         call dgemv('normal',ncnln,n,one,cjacu,ldcju,dx,1,zero,cjdx,1)
c

c        make forward-difference approximation along dx.

         call dcopy(n,x,1,y,1)
         call daxpy(n,h,dx,1,y,1)
c
         call iload (ncnln,(1),needc,1)
c
         mode = 0
         call confun(mode,ncnln,n,ldcju,needc,y,c1,cjacu,nstate,iuser,
     *               user)
         if (mode.lt.0) go to 160
c
c        set  err = (c1 - c)/h  - jacobian*dx.  this should be small.
c
         do 80 i = 1, ncnln
            err(i) = (c1(i)-c(i))/h - cjdx(i)
   80    continue
         imax = idamax(ncnln,err,1)
         emax = abs(err(imax))/(abs(cjdx(imax))+one)
c

         if (emax.ge.point9) inform = 1
      end if
c

c     element-wise check.

      if (lvrfyc.ge.2) then
         if (lvlder.eq.3) then
c
c           recompute the jacobian to find the non-constant elements.
c
            call smload('general',ncnln,n,rdummy,rdummy,cjacu,ldcju)
c
            call iload (ncnln,(1),needc,1)
            nstate = 0
            mode = 2
c
            call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,nstate,
     *                  iuser,user)
            if (mode.lt.0) go to 160
c
         end if
c
         call iload (ncnln,(0),needc,1)
c
         itmax = 3
         ncheck = 0
         nwrong = 0
         ngood = 0
         colmax = -one
         jcol = 0
         irow = 0
         mode = 0
         j3 = jverfy(3)
         j4 = jverfy(4)
c

c        loop over each column.

         do 140 j = j3, j4
c
            call sload(ncnln,zero,err,1)
            ncolj = 0
            headng = .true.
            xj = x(j)
c
            stepbl = biglow
            stepbu = bigupp
            if (bl(j).gt.biglow) stepbl = bl(j) - xj
            if (bu(j).lt.bigupp) stepbu = bu(j) - xj
c
            signh = one
            if (half*(stepbl+stepbu).lt.zero) signh = -one
c
            do 120 i = 1, ncnln
               epsaci = epsrf*(one+abs(c(i)))
c
               if (cjacu(i,j).ne.rdummy) then
c                 ------------------------------------------------------
c                 check this jacobian element.
c                 ------------------------------------------------------
                  ncheck = ncheck + 1
                  ncolj = ncolj + 1
                  needc(i) = 1
c
                  cij = cjac(i,j)
                  cjsize = abs(cij)
c                 ------------------------------------------------------
c                 find a finite-difference interval by iteration.
c                 ------------------------------------------------------
                  iter = 0
                  hopt = two*(one+abs(xj))*sqrt(epsrf)
                  h = ten*hopt*signh
                  cdest = zero
                  sdest = zero
                  first = .true.
c
c                 +                repeat
  100             x(j) = xj + h
                  call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 160
                  f1 = c1(i)
c
                  x(j) = xj + h + h
                  call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 160
                  f2 = c1(i)
c
                  call chcore(debug,done,first,epsaci,epsrf,c(i),info,
     *                        iter,itmax,cdest,fdest,sdest,errbnd,f1,f2,
     *                        h,hopt,hphi)
c
c                 +                until     done
                  if ( .not. done) go to 100
c
c                 ------------------------------------------------------
c                 exit for this element.
c                 ------------------------------------------------------
                  cjdiff = cdest
                  err(i) = abs(cjdiff-cij)/(cjsize+one)
c
                  ok = err(i) .le. oktol
                  if (ok) then
                     key = lgood
                     ngood = ngood + 1
                  else
                     key = lbad
                     nwrong = nwrong + 1
                  end if
c
                  needc(i) = 0
               end if
  120       continue
c

c           finished with this column.

            if (ncolj.gt.0) then
               imax = idamax(ncnln,err,1)
               emax = abs(err(imax))
c
               if (emax.ge.colmax) then
                  irow = imax
                  jcol = j
                  colmax = emax
               end if
            end if
            x(j) = xj
c
  140    continue
c
         inform = 0
         if (colmax.ge.point9) inform = 1

c
      end if
c
c     copy  ( constants + gradients + dummy values )  back into cjacu.
c
      call smcopy('general',ncnln,n,cjac,ldcj,cjacu,ldcju)
c
      return

  160 inform = mode

c                                 end of chcjac

      end

      subroutine lsgetp(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,ap,res,
     *                  hz,p,gq,cq,r,zy,work)
c----------------------------------------------------------------------
c     lsgetp  computes the following quantities for  lscore.
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

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  ctp, pnorm
      integer           lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf
      logical           linobj, singlr, unitgz, unitq

      double precision  a(lda,*), ap(*), cq(*), gq(n), hz(*), p(n),
     *                  r(ldr,*), res(*), work(n), zy(ldzy,*)
      integer           kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  gtp
      integer           i, j

      character*80      rec(2)
c     .. external functions ..
      double precision  ddot, dnrm2
      external          ddot, dnrm2
c     .. external subroutines ..
      external          dcopy, dgemv, dscal, dtrsv, cmqmul, sload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
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
         gtp = ddot(nrz,gq,1,p,1)
         if (gtp.gt.zero) call dscal(nrz,(-one),p,1)
c
         if (nrz.le.nrank) then
            if (numinf.eq.0) then
               if (unitgz) then
                  hz(nrz) = r(nrz,nrz)*p(nrz)
               else
                  call sload(nrz,(zero),hz,1)
               end if
            else
               hz(1) = r(1,1)*p(1)
            end if
         end if
      else

c        the objective is quadratic in the space spanned by z1.

         if (linobj) then
            if (unitgz) then
               if (nrz.gt.1) call sload(nrz-1,(zero),hz,1)
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
      if (linobj) ctp = ddot(nrz,cq,1,p,1)
      pnorm = dnrm2(nrz,p,1)
c
      call cmqmul(1,n,nrz,nfree,ldzy,unitq,kx,p,zy,work)

c     compute  ap.

      if (nclin.gt.0) then
         call dgemv('no transpose',nclin,n,one,a,lda,p,1,zero,ap,1)
      end if

c                                 end of lsgetp

      end

      subroutine lsgset(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,res,featol,gq,cq,r,
     *                  x,wtinf,zy,wrk)
c----------------------------------------------------------------------
c     lsgset  finds the number and weighted sum of infeasibilities for
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
c     if  x  is feasible,  lsgset computes the vector q(free)'g(free),
c     where  g  is the gradient of the the sum of squares plus the
c     linear term.  the matrix q is of the form
c                    ( q(free)  0       ),
c                    (   0      i(fixed))
c     where  q(free)  is the orthogonal factor of  a(free)  and  a  is
c     the matrix of constraints in the working set.  the transformed
c     gradients are stored in gq.

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  bigbnd, suminf, tolrnk
      integer           lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf, nz
      logical           linobj, singlr, unitgz, unitq
      character*2       prbtyp

      double precision  a(lda,*), bl(*), bu(*), cq(*), featol(*), gq(n),
     *                  r(ldr,*), res(*), wrk(n), wtinf(*), x(n),
     *                  zy(ldzy,*)
      integer           istate(*), kx(n)

      double precision  biglow, bigupp, ctx, feasj, rownrm, s, weight
      integer           j, k
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           isrank
      external          ddot, dnrm2, isrank
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dtrmv, cmqmul, sload,
     *                  sscmv 
c     .. intrinsic functions ..
      intrinsic         abs, min
c----------------------------------------------------------------------
      bigupp = bigbnd
      biglow = -bigbnd
c
      numinf = 0
      suminf = zero
      call sload(n,zero,gq,1)
c
      do 40 j = 1, n + nclin
         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot(n,a(k,1),lda,x,1)
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
               call daxpy(n,weight,a(k,1),lda,gq,1)
            end if
         end if
   40 continue
c

c     install  gq,  the transformed gradient.

      singlr = .false.
      unitgz = .true.
c
      if (numinf.gt.0) then
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,gq,zy,wrk)
         unitgz = .true.
      else if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         call sload(n,zero,gq,1)
      else
c
c        ready for the optimality phase.
c        set nrz so that rz1 is nonsingular.
c
         if (nrank.eq.0) then
            if (linobj) then
               call dcopy(n,cq,1,gq,1)
            else
               call sload(n,zero,gq,1)
            end if
            nrz = 0
         else
c
c           compute gq = - r' * (transformed residual)
c
            call sscmv (nrank,-one,res,1,gq,1)
            call dtrmv('u','t','n',nrank,r,ldr,gq,1)
            if (nrank.lt.n) call dgemv('t',nrank,n-nrank,-one,
     *                                 r(1,nrank+1),ldr,res,1,zero,
     *                                 gq(nrank+1),1)
            if (linobj) call daxpy(n,one,cq,1,gq,1)
            unitgz = .false.
            rownrm = dnrm2(n,r(1,1),ldr)
            if (rownrm.le.tolrnk .or. abs(r(1,1)).le.rownrm*tolrnk) then
               nrz = 0
            else
               nrz = isrank(min(nrank,nz),r,ldr+1,tolrnk)
            end if
         end if
         singlr = .false.
      end if

c                                 end of lsgset

      end

      subroutine lsfeas(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *                  ax,bl,bu,featol,x,work)
c----------------------------------------------------------------------

c     lsfeas  computes 
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd, cvnorm, errmax
      integer           jmax, n, nclin, nviol

      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), work(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  biglow, bigupp, con, feasj, res, tolj
      integer           i, is, j

      character*80      rec(2)
c     .. external functions ..
      double precision  dnrm2
      integer           idamax
      external          dnrm2, idamax

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
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
         if (is.ge.0 .and. is.lt.4) then
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

      jmax = idamax(n+nclin,work,1)
      errmax = abs(work(jmax))

      cvnorm = dnrm2(n+nclin,work,1)

c                                 end of lsfeas

      end

      subroutine cmmul2(msglvl,n,nrz,nz,zerolm,notopt,numinf,trusml,
     *                  smllst,jsmlst,tinyst,jtiny,gq)

c     cmmul2  updates jsmlst and smllst when there are artificial
c     constraints.
c
c     on input,  jsmlst  is the index of the minimum of the set of
c     adjusted multipliers.
c     on output, a negative jsmlst defines the index in q'g of the
c     artificial constraint to be deleted.

      integer           ldbg
      parameter         (ldbg=5)

      double precision  smllst, tinyst, trusml, zerolm
      integer           jsmlst, jtiny, msglvl, n, notopt, nrz, numinf,
     *                  nz

      double precision  gq(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)

      double precision  rlam
      integer           j, k, kk, length

      character*80      rec(3)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
c
      do 20 j = nrz + 1, nz
         rlam = -abs(gq(j))
c
         if (rlam.lt.zerolm) then
            if (numinf.eq.0) notopt = notopt + 1
c
            if (rlam.lt.smllst) then
               trusml = gq(j)
               smllst = rlam
               jsmlst = -j
            end if
c
         else if (rlam.lt.tinyst) then
            tinyst = rlam
            jtiny = -j
         end if
   20 continue

c                                 end of cmmul2

      end

      subroutine nggnfm(side,n,k1,k2,s,a,lda)
c----------------------------------------------------------------------
c  nggnfm applies a  sequence  of  pairwise interchanges to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the interchanges are
c  applied in planes k1 up to k2.
c  based on f06qnf.
c
c  the upper hessenberg matrix, h, is formed as
c
c     h = p*u,    when   side = 'l' or 'l',  (  left-hand side )
c
c  where p is a permutation matrix of the form
c
c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 )
c
c  and is formed as
c
c     h = u*p',   when   side = 'r' or 'r',  ( right-hand side )
c
c  where p is a permutation matrix of the form
c
c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c  p( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
c  the  two by two
c  interchange part of p( k ), r( k ), is assumed to have the form
c
c     r( k ) = ( 0  1 ).
c              ( 1  0 )
c
c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of  h.
c
c  the  sub-diagonal elements of  h, h( k + 1, k ),  are returned in the
c  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision  zero
      parameter         (zero=0.0d+0)

      integer           k1, k2, lda, n
      character*1       side

      double precision  a(lda,*), s(*)

      double precision  aij, temp
      integer           i, j
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return
      if (side.eq.'l') then
c
c        apply the permutations to columns n back to k1.
c
         do 40 j = n, k1, -1
            if (j.ge.k2) then
               aij = a(k2,j)
            else
c
c              form  the  additional sub-diagonal element  h( j + 1, j )
c              and store it in s( j ).
c
               aij = zero
               s(j) = a(j,j)
            end if
            do 20 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i+1,j) = temp
               aij = aij
   20       continue
            a(k1,j) = aij
   40    continue
      else if (side.eq.'r') then
c
c        apply  the  plane interchanges to  columns  k1  up to
c        ( k2 - 1 ) and  form   the   additional  sub-diagonal
c        elements,   storing  h( j + 1, j ) in s( j ).
c
         do 80 j = k1, k2 - 1
            do 60 i = 1, j
               temp = a(i,j+1)
               a(i,j+1) = a(i,j)
               a(i,j) = temp
   60       continue
            s(j) = a(j+1,j+1)
            a(j+1,j+1) = zero
   80    continue
      end if

c                                 end of nggnfm.

      end

      subroutine srchq(first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,
     *                  tolrel,toltny,alfa,alfbst,fbest)
c----------------------------------------------------------------------
c     srchq  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     srchq  requires  f(alpha) (but not f'(alpha)) to be evaluated
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
c                   it is subsequently altered by srchq.
c
c     debug         specifies whether detailed output is wanted.
c
c     maxf          is an upper limit on the number of times srchq is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).
c
c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by srchq (see below).
c
c     alfmax        is the upper limit of the interval to be searched.
c
c     alfsml        is intended to prevent inefficiency when a minimizer
c                   is very small, for cases where the calling program
c                   would prefer to redefine f'(alfa).  alfsml is
c                   allowed to be zero.  early termination will occur if
c                   srchq determines that a minimizer lies somewhere in
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
c                   these values may be reduced by srchc .
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
c                   for the new trial step alfa, and reenter srchq.
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
c     numf          counts the number of times srchq has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).
c
c     alfa          is the point at which the next function ftry must
c                   be computed.
c
c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever srchq returns
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

      double precision  zero, point1, half
      parameter         (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision  one, two, five
      parameter         (one=1.0d+0,two=2.0d+0,five=5.0d+0)
      double precision  ten, eleven
      parameter         (ten=1.0d+1,eleven=1.1d+1)

      double precision  alfa, alfbst, alfmax, alfsml, epsaf, fbest,
     *                  ftry, g0, targtg, tolabs, tolrel, toltny
      integer           inform, maxf, nout, numf
      logical           debug, done, first, imprvd

      double precision  a, artifa, artifb, b, daux, dtry, endpnt, fa,
     *                  factor, fv, fw, gv, gw, q, s, tol, tolmax,
     *                  truea, trueb, xmidpt, xtry, xv, xw
      integer           nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, found,
     *                  moved, quitf, quitfz, quiti, quits, setxv, vset,
     *                  wset, xinxw

      character*120     rec(7)

      save              braktd, crampd, extrap, moved, vset, wset,
     *                  nsamea, nsameb, a, b, fa, factor, xtry, xw, fw,
     *                  xv, fv, tolmax
c----------------------------------------------------------------------
c

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
         badfun = alfmax .le. toltny .or. g0 .ge. zero
         done = badfun
         moved = .false.
c
         if ( .not. done) then
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
         end if
      else

c        subsequent entries.  the function has just been evaluated at
c        alfa = alfbst + xtry,  giving ftry.


         numf = numf + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if ( .not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if
c
c        check if xtry is in the interval (xw,0) or (0,xw).
c
         if (wset) then
            xinxw = zero .lt. xtry .and. xtry .le. xw .or. xw .le.
     *              xtry .and. xtry .lt. zero
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
         else if (closef .and. ftry-fbest.lt.epsaf) then
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
                     setxv = ftry .lt. fv .or. .not. extrap
                  else
                     setxv = .true.
                  end if
c
                  if (setxv) then
                     if (vset .and. xinxw) then
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
         found = moved .and. abs(fa-fbest) .le. -a*targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol
         quits = trueb .le. alfsml
c
         if (quiti .and. .not. moved) then
c
c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.
c
            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if
c
         done = quitf .or. quitfz .or. quits .or. quiti .or. found


c        proceed with the computation of an estimate of a minimizer.
c        the choices are...
c        1. parabolic fit using function values only.
c        2. damped parabolic fit if the regular fit appears to be
c           consistently overestimating the distance to a minimizer.
c        3. bisection, geometric bisection, or a step of tol if the
c           parabolic fit is unsatisfactory.

         if ( .not. done) then
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
            if (vset .and. moved) then
c
c              three points available.  use fbest, fw and fv.
c
               gv = (fv-fbest)/xv
               s = gv - (xv/xw)*gw
               q = two*(gv-gw)

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

            end if

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
            if ( .not. braktd) then
c
c              a minimizer has not yet been bracketed.
c              set an artificial upper bound by expanding the interval
c              xw  by a suitable factor.
c
               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = five*factor
            else if (vset .and. moved) then
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
                  if (nsamea.ge.3 .or. nsameb.ge.3) then
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
                  xtry = half*sqrt(daux)*sqrt(dtry)
               end if
               if (endpnt.lt.zero) xtry = -xtry

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

            end if


c           the polynomial fits give (s/q)*xw as the new step.  reject
c           this step if it lies outside (artifa, artifb).

            if (q.ne.zero) then
               if (q.lt.zero) s = -s
               if (q.lt.zero) q = -q
               if (s*xw.ge.q*artifa .and. s*xw.le.q*artifb) then
c
c                 accept the polynomial fit.
c
                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if

               end if
            end if
         end if
      end if

c
      if ( .not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then
c
c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)
c
            xmidpt = half*(a+b)
            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
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

c                                 end of srchq

      end

      subroutine cmmul1(prbtyp,msglvl,n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,anorms,gq,rlamda,t,wtinf)

c     cmmul1  first computes the lagrange multiplier estimates for the
c     given working set.  it then determines the values and indices of
c     certain significant multipliers.  in this process, the multipliers
c     for inequalities at their upper bounds are adjusted so that a
c     negative multiplier for an inequality constraint indicates non-
c     optimality.  all adjusted multipliers are scaled by the 2-norm
c     of the associated constraint row.  in the following, the term
c     minimum refers to the ordering of numbers on the real line,  and
c     not to their magnitude.
c
c     jsmlst          is the index of the constraint whose multiplier is
c                     the minimum of the set of adjusted multipliers
c                     with values less than  small.
c     rlamda(ksmlst)  is the associated multiplier.
c
c     jbigst          is the index of the constraint whose multiplier is
c                     the largest of the set of adjusted multipliers
c                     with values greater than (1 + small).
c     rlamda(kbigst)  is the associated multiplier.
c
c     on exit,  elements  1  thru  nactiv  of  rlamda  contain the
c     unadjusted multipliers for the general constraints.  elements
c     nactiv  onwards of  rlamda  contain the unadjusted multipliers
c     for the bounds.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)

      double precision  biggst, smllst, tinyst, trubig, trusml, zerolm
      integer           jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, msglvl, n, nactiv, nfree, notopt,
     *                  numinf, nz
      character*2       prbtyp

      double precision  a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lcdbg
c     .. arrays in common ..
      integer           ilcdbg(ldbg)

      double precision  anormj, blam, rlam, scdlam
      integer           i, is, j, k, kk, l, nfixed

      character*80      rec(3)
c     .. external subroutines ..
      external          dcopy, dtrsv

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg009/ilcdbg, lcdbg
c----------------------------------------------------------------------
c
      nfixed = n - nfree
c
      jtiny = 0
      jsmlst = 0
      ksmlst = 0
c
      jbigst = 0
      kbigst = 0
c

c     compute  jsmlst  for regular constraints and temporary bounds.

c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.
c
      if (n.gt.nz) call dcopy(n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call dtrsv('u','t','n',nactiv,t(1,nz+1),ldt,
     *                            rlamda,1)

c     set elements  nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.

      do 40 l = 1, nfixed
         j = kx(nfree+l)
         blam = rlamda(nactiv+l)
         do 20 k = 1, nactiv
            i = kactiv(k)
            blam = blam - a(i,j)*rlamda(nactiv-k+1)
   20    continue
         rlamda(nactiv+l) = blam
   40 continue

c     find  jsmlst  and  ksmlst.

      do 60 k = 1, n - nz
         if (k.gt.nactiv) then
            j = kx(nz+k)
         else
            j = kactiv(nactiv-k+1) + n
         end if

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
c
            if (scdlam.lt.zerolm) then
               if (numinf.eq.0) notopt = notopt + 1
c
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
c
         scdlam = rlam/wtinf(j)
         if (scdlam.gt.biggst .and. j.gt.jinf) then
            biggst = scdlam
            trubig = rlamda(k)
            jbigst = j
            kbigst = k
         end if
   60 continue
c                                 end of cmmul1
      end

      subroutine chkgrd(inform,msglvl,n,bigbnd,epsrf,oktol,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,y,iuser,user)
c----------------------------------------------------------------------
c     chkgrd  checks if the gradients of the objective function have
c     been coded correctly.
c
c     on input,  the value of the objective function at the point x is
c     stored in objf.  the corresponding gradient is stored in gradu.
c     if any gradient element has not been specified,  it will have a
c     dummy value.  missing values are not checked.
c
c     a cheap test is first undertaken by calculating the directional
c     derivative using two different methods. if this proves
c     satisfactory and no further information is desired, chkgrd is
c     terminated. otherwise, the routine chcore is called to give
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
c     -1        do not do any check.
c     0        do the cheap test only.
c     1 or 3   do both cheap and full test.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, point9
      parameter         (zero=0.0d+0,half=0.5d+0,point9=0.9d+0)
      double precision  one, two, ten
      parameter         (one=1.0d+0,two=2.0d+0,ten=1.0d+1)
      character*4       lbad, lgood
      parameter         (lbad='bad?',lgood='  ok')

      double precision  bigbnd, epsrf, fdchk, objf, oktol, xnorm
      integer           inform, msglvl, n

      double precision  bl(n), bu(n), dx(n), grad(n), gradu(n), user(*),
     *                  x(n), y(n)
      integer           iuser(*)
c     .. subroutine arguments ..
      external          objfun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, lvrfyc, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg), jverfy(4)

      double precision  biglow, bigupp, cdest, dxj, dxmult, emax, epsa,
     *                  errbnd, error, f1, f2, fdest, gdiff, gdx, gj,
     *                  gsize, h, hopt, hphi, objf1, sdest, stepbl,
     *                  stepbu, xj
      integer           info, iter, itmax, j, j1, j2, jmax, mode,
     *                  ncheck, ngood, nstate, nwrong
      logical           const, debug, done, first, headng, needed, ok
      character*4       key

      character*18      result(0:4)
      character*120     rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. external subroutines ..
      external          daxpy, dcopy, chcore

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg015/lvrfyc, jverfy
      common            /ngg002/inpdbg, npdbg

      data              result/'                 ', 'constant?      ',
     *                  'linear or odd?   ', 'too nonlinear?',
     *                  'small derivative?'/
c----------------------------------------------------------------------
c
      inform = 0
      needed = lvrfyc .eq. 0 .or. lvrfyc .eq. 1 .or. lvrfyc .eq. 3
      if ( .not. needed) return
c
      debug = npdbg .and. inpdbg(5) .gt. 0
      nstate = 0
c
      biglow = -bigbnd
      bigupp = bigbnd
c

c     do the cheap test.

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
            if (bu(j).lt.bigupp .and. bu(j).gt.bl(j))
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

      gdx = ddot(n,gradu,1,dx,1)
c

c     make forward-difference approximation along  p.

      call dcopy(n,x,1,y,1)
      call daxpy(n,h,dx,1,y,1)
c
      mode = 0
      call objfun(mode,n,y,objf1,gradu,nstate,iuser,user)
      if (mode.lt.0) go to 100
c
      gdiff = (objf1-objf)/h
      error = abs(gdiff-gdx)/(abs(gdx)+one)
c
      ok = error .le. oktol

c
      if (error.ge.point9) inform = 1
c

c     element-wise check.

      if (lvrfyc.eq.1 .or. lvrfyc.eq.3) then
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

c              check this gradient element.

               ncheck = ncheck + 1
               gj = grad(j)
               gsize = abs(gj)
               xj = x(j)

c              find a finite-difference interval by iteration.

               iter = 0
               epsa = epsrf*(one+abs(objf))
               cdest = zero
               sdest = zero
               first = .true.
c
               stepbl = biglow
               stepbu = bigupp
               if (bl(j).gt.biglow) stepbl = bl(j) - xj
               if (bu(j).lt.bigupp) stepbu = bu(j) - xj
c
               hopt = two*(one+abs(xj))*sqrt(epsrf)
               h = ten*hopt
               if (half*(stepbl+stepbu).lt.zero) h = -h
c
c              +             repeat
   60          x(j) = xj + h
               call objfun(mode,n,x,f1,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100
c
               x(j) = xj + h + h
               call objfun(mode,n,x,f2,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100
c
               call chcore(debug,done,first,epsa,epsrf,objf,info,iter,
     *                     itmax,cdest,fdest,sdest,errbnd,f1,f2,h,hopt,
     *                     hphi)
c
c              +             until     done
               if ( .not. done) go to 60
c

c              exit for this variable.

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
            end if
   80    continue
c

c        done.

         inform = 0
         if (error.ge.point9) inform = 1
      end if
c
      call dcopy(n,grad,1,gradu,1)
c
      return
c
  100 inform = mode

c                                 end of chkgrd

      end

      subroutine npiqp (feasqp,unitq,nqperr,majits,minits,n,nclin,ncnln,
     *                  ldcj,ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,
     *                  numinf,istate,kactiv,kx,dxnorm,gdx,qpcurv,aqp,
     *                  adx,anorm,ax,bl,bu,c,cjac,clamda,cmul,cs,dlam,
     *                  dslk,dx,qpbl,qpbu,qptol,r,rho,slk,violn,x,wtinf,
     *                  w)
c----------------------------------------------------------------------
c     npiqp    does the following:
c
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
c                 subject to   qpbl .le. (  p ) .le. qpbu,
c                                        ( ap )
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

      integer           lenls
      parameter         (lenls=20)
      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      logical           qpnamd, vertex
      parameter         (qpnamd=.false.,vertex=.false.)
      double precision  zero, one, two
      parameter         (zero=0.0d+0,one=1.0d+0,two=2.0d+0)
      double precision  hundrd
      parameter         (hundrd=1.0d+2)

      double precision  dxnorm, gdx, qpcurv
      integer           ldaqp, ldcj, ldr, linact, majits, minits, n,
     *                  nactiv, nclin, ncnln, nfree, nlnact, nqperr,
     *                  numinf, nz
      logical           feasqp, unitq

      double precision  adx(*), anorm(*), aqp(ldaqp,*), ax(*), bl(*),
     *                  bu(*), c(*), cjac(ldcj,*), clamda(*), cmul(*),
     *                  cs(*), dlam(*), dslk(*), dx(n), qpbl(*),
     *                  qpbu(*), qptol(*), r(ldr,*), rho(*), slk(*),
     *                  violn(*), w(*), wtinf(*), x(n)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer           idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldt,
     *                  ldzy, lennam, lformh, lines1, lines2, lprob,
     *                  lverfy, lvlder, msgls, msgnp, ncolt, nlnf, nlnj,
     *                  nlnx, nload, nn, nnclin, nncnln, nout, nprob,
     *                  nsave
      logical           cmdbg, incrun, lsdbg, npdbg
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           icmdbg(ldbg), ilsdbg(ldbg), inpdbg(ldbg),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), locls(lenls)

      double precision  amin, biglow, bigupp, blj, buj, con, condmx,
     *                  quotnt, ssq, ssq1, suminf, viol, weight, wscale,
     *                  wtmax, wtmin
      integer           i, idbg, idbgsv, inform, iswap, j, jinf, k, k1,
     *                  k2, kviol, l, lgq, lhpq, lrlam, lrpq, lrpq0, lt,
     *                  lwrk1, lzy, mjrdbg, mnrdbg, msgqp, nartif, ncqp,
     *                  nctotl, ngq, nmajor, nminor, nplin, nrank,
     *                  nrejtd, nrpq, ntry, nviol, nz1
      logical           linobj, overfl

      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(3)
c     .. external functions ..
      double precision  ddot, dnrm2, sdiv 
      external          ddot, dnrm2, sdiv 
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dscal, dtrmv, dtrsv,
     *                  cmqmul, lsadds, lscore, npsetx, iload , icopy ,
     *                  sload, scond , smcopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg012/locls
      common            /ngg004/lennam, ldt, ncolt, ldzy
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg016/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ngg018/rcndbd, rfrobn, drmax, drmin
      common            /ngg019/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ngg4fb/icmdbg, cmdbg
      common            /ngg002/inpdbg, npdbg
      common            /ngg020/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /ngg021/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c----------------------------------------------------------------------

      call icopy (ldbg,ilsdbg,1,icmdbg,1)
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
      nctotl = nplin + ncnln
      ncqp = nclin + ncnln
      nrank = n
      nrejtd = 0


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
         else
            con = c(j-nplin)
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
c
      if (ncnln.gt.0) then

c        refactorize part of the  qp  constraint matrix.

c        load the new jacobian into the  qp  matrix  a.  compute the
c        2-norms of the rows of the jacobian.
c
         call smcopy('general',ncnln,n,cjac,ldcj,aqp(nclin+1,1),ldaqp)
c
         do 80 j = nclin + 1, ncqp
            anorm(j) = dnrm2(n,aqp(j,1),ldaqp)
   80    continue
c
c        count the number of linear constraints in the working set and
c        move them to the front of kactiv.  compute the norm of the
c        matrix of constraints in the working set.
c        let k1  point to the first nonlinear constraint.  constraints
c        with indices kactiv(k1),..., kactiv(nactiv)  must be
c        refactorized.
c
         asize = zero
         linact = 0
         k1 = nactiv + 1
         do 100 k = 1, nactiv
            i = kactiv(k)
            asize = max(asize,anorm(i))
c
            if (i.le.nclin) then
               linact = linact + 1
               if (linact.ne.k) then
                  iswap = kactiv(linact)
                  kactiv(linact) = i
                  kactiv(k) = iswap
               end if
            else
c
c              record the old position of the 1st. nonlinear constraint.
c
               if (k1.gt.nactiv) k1 = k
            end if
  100    continue
c
         if (nactiv.le.1) call scond (ncqp,anorm,1,asize,amin)
c
c        compute the absolute values of the nonlinear constraints in
c        the working set.  use dx as workspace.
c
         do 120 k = linact + 1, nactiv
            j = n + kactiv(k)
            if (istate(j).eq.1) dx(k) = abs(qpbl(j))
            if (istate(j).ge.2) dx(k) = abs(qpbu(j))
  120    continue
c
c        sort the elements of kactiv corresponding to nonlinear
c        constraints in descending order of violation (i.e.,
c        the first element of kactiv for a nonlinear constraint
c        is associated with the most violated constraint.)
c        in this way, the rows of the jacobian corresponding
c        to the more violated constraints tend to be included
c        in the  tq  factorization.
c
c        the sorting procedure is taken from the simple insertion
c        sort in d. knuth, acp volume 3, sorting and searching,
c        page 81.  it should be replaced by a faster sort if the
c        number of active nonlinear constraints becomes large.
c
         do 160 k = linact + 2, nactiv
            l = k
            viol = dx(l)
            kviol = kactiv(l)
c           while (l .gt. linact+1  .and.  dx(l-1) .lt. viol) do
  140       if (l.gt.linact+1) then
               if (dx(l-1).lt.viol) then
                  dx(l) = dx(l-1)
                  kactiv(l) = kactiv(l-1)
                  l = l - 1
                  go to 140
               end if
c              end while
            end if
            dx(l) = viol
            kactiv(l) = kviol
  160    continue
c
         k2 = nactiv
         nactiv = k1 - 1
         nz = nfree - nactiv
c
c        update the factors  r,  t  and  q  to include constraints
c        k1  through  k2.
c
         if (k1.le.k2) call lsadds(unitq,vertex,inform,k1,k2,nactiv,
     *                             nartif,nz,nfree,nrank,nrejtd,nrpq,
     *                             ngq,n,ldzy,ldaqp,ldr,ldt,istate,
     *                             kactiv,kx,condmx,aqp,r,w(lt),w(lrpq),
     *                             w(lgq),w(lzy),w(lwrk1),dx,w(lrlam),
     *                             msgqp)
      end if
c

c     solve for dx, the vector of minimum two-norm that satisfies the
c     constraints in the working set.

      call npsetx(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,ldaqp,
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
  180 call lscore('qp subproblem',qpnamd,names,linobj,unitq,nqperr,
     *            minits,jinf,ncqp,nctotl,nactiv,nfree,nrank,nz,nz1,n,
     *            ldaqp,ldr,istate,kactiv,kx,gdx,ssq,ssq1,suminf,numinf,
     *            dxnorm,qpbl,qpbu,aqp,clamda,adx,qptol,r,dx,w)

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
            call iload (nctotl,(0),istate,1)
c
            call npsetx(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,qpbl,qpbu,w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt)
     *                  ,w(lzy),w(lwrk1))
         end if
      end if
      if ( .not. (nviol.eq.0 .or. ntry.gt.2)) go to 180
c     +    until (    nviol .eq. 0  .or.  ntry .gt. 2)
c

c     count the number of nonlinear constraint gradients in the  qp
c     working set.  make sure that all small  qp  multipliers associated
c     with nonlinear inequality constraints have the correct sign.

      nlnact = 0
      if (nactiv.gt.0 .and. ncnln.gt.0) then
         do 220 k = 1, nactiv
            l = kactiv(k)
            if (l.gt.nclin) then
               nlnact = nlnact + 1
               j = n + l
               if (istate(j).eq.1) clamda(j) = max(zero,clamda(j))
               if (istate(j).eq.2) clamda(j) = min(zero,clamda(j))
            end if
  220    continue
      end if

      linact = nactiv - nlnact


c     extract various useful quantities from the qp solution.

c     compute  hpq = r'r(pq)  from the transformed gradient of the qp
c     objective function and  r(pq)  from the transformed residual.

      call dscal(n,(-one),w(lrpq),1)
      call daxpy(n,(-one),w(lgq),1,w(lhpq),1)
      qpcurv = two*ssq

      if (ncnln.gt.0) then
         if (numinf.gt.0) then
            feasqp = .false.
            call sload(nctotl,(zero),clamda,1)
c
            if (nz.gt.0) then

c              compute a null space element for the search direction
c              as the solution of  z'hz(pz) = -z'g - z'hy(py).

c              overwrite dx with the transformed search direction
c              q'(dx).  the first nz elements of dx are zero.
c
               call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,dx,w(lzy),w(lwrk1)
     *                     )
c
c              overwrite the first nz elements of dx with the solution
c              of  (rz)u = -(v + w),  where  (rz)'w = z'g  and  v  is
c              vector of first nz elements of  r(pq).
c
               call dcopy(nz,w(lgq),1,dx,1)
               call dtrsv('u','t','n',nz,r,ldr,dx,1)
c
               call daxpy(nz,(one),w(lrpq),1,dx,1)
c
               call dtrsv('u','n','n',nz,r,ldr,dx,1)
               call dscal(nz,(-one),dx,1)
c
c              recompute rpq, hpq, gdx and qpcurv.
c
               call dcopy(nlnx,dx,1,w(lrpq),1)
               call dtrmv('u','n','n',nlnx,r,ldr,w(lrpq),1)
               if (nlnx.lt.n) call dgemv('n',nlnx,n-nlnx,one,r(1,nlnx+1)
     *                                   ,ldr,dx(nlnx+1),1,one,w(lrpq),
     *                                   1)
c
               gdx = ddot(n,w(lgq),1,dx,1)
               qpcurv = ddot(n,w(lrpq),1,w(lrpq),1)
c
               call cmqmul(3,n,nz,nfree,ldzy,unitq,kx,dx,w(lzy),w(lwrk1)
     *                     )
c

c              recompute adx and the 2-norm of dx.

               dxnorm = dnrm2(n,dx,1)
               if (ncqp.gt.0) call dgemv('n',ncqp,n,one,aqp,ldaqp,dx,1,
     *                                   zero,adx,1)

            end if
c
            call dcopy(nlnx,w(lrpq),1,w(lhpq),1)
            call dtrmv('u','t','n',nlnx,r,ldr,w(lhpq),1)
            if (nlnx.lt.n) call dgemv('t',nlnx,n-nlnx,one,r(1,nlnx+1),
     *                                ldr,w(lrpq),1,zero,w(lhpq+nlnx),1)
         end if
c

c        for given values of the objective function and constraints,
c        attempt to minimize the merit function with respect to each
c        slack variable.

         do 280 i = 1, ncnln
            j = nplin + i
            con = c(i)
c
            if ( .not. feasqp .and. violn(i).ne.zero .and. rho(i)
     *          .le.zero) rho(i) = one
c
            quotnt = sdiv (cmul(i),scale*rho(i),overfl)
c
c           define the slack variable to be  con - mult / rho.
c           force each slack to lie within its upper and lower bounds.
c
            if (bl(j).gt.biglow) then
               if (qpbl(j).ge.-quotnt) then
                  slk(i) = bl(j)
                  go to 260
               end if
            end if
c
            if (bu(j).lt.bigupp) then
               if (qpbu(j).le.-quotnt) then
                  slk(i) = bu(j)
                  go to 260
               end if
            end if
c
            slk(i) = con - quotnt
c
c           the slack has been set within its bounds.
c
  260       cs(i) = con - slk(i)
c

c           compute the search direction for the slacks and multipliers.

            dslk(i) = adx(nclin+i) + cs(i)
c
            if (feasqp) then
c
c              if any constraint is such that  (dlam)*(c - s)  is
c              positive,  the merit function may be reduced immediately
c              by substituting the qp multiplier.
c
               dlam(i) = clamda(j) - cmul(i)
               if (dlam(i)*cs(i).ge.zero) then
                  cmul(i) = clamda(j)
                  dlam(i) = zero
               end if
            else
c
c              the  qp  subproblem was infeasible.
c
               dlam(i) = zero
c
               if (istate(j).lt.0 .or. violn(i).ne.zero) dslk(i) = zero
c
            end if
  280    continue
c
         if ( .not. feasqp) rhonrm = dnrm2(ncnln,rho,1)

      end if

      call icopy (ldbg,inpdbg,1,icmdbg,1)

c                                 end of npiqp

      end


      subroutine cmsinf(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,a,
     *                  featol,cvec,x,wtinf)

c     cmsinf  finds the number and weighted sum of infeasibilities for
c     the bounds and linear constraints.   an appropriate gradient
c     is returned in cvec.
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

      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd, suminf
      integer           lda, n, nclin, numinf

      double precision  a(lda,*), bl(*), bu(*), cvec(n), featol(*),
     *                  wtinf(*), x(n)
      integer           istate(*)

      double precision  biglow, bigupp, ctx, feasj, s, weight
      integer           j, k
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. external subroutines ..
      external          daxpy, sload
c     .. intrinsic functions ..
      intrinsic         abs
c----------------------------------------------------------------------
c
      bigupp = bigbnd
      biglow = -bigbnd
c
      numinf = 0
      suminf = zero
      call sload(n,(zero),cvec,1)
c
      do 40 j = 1, n + nclin
         if (istate(j).le.0) then
            feasj = featol(j)
            if (j.le.n) then
               ctx = x(j)
            else
               k = j - n
               ctx = ddot(n,a(k,1),lda,x,1)
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
               cvec(j) = weight
            else
               call daxpy(n,weight,a(k,1),lda,cvec,1)
            end if
         end if
   40 continue

c                                 end of cmsinf

      end

      subroutine rzdel(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,t,gqm,q,work,c,s)
c----------------------------------------------------------------------
c     rzdel   updates the matrices  z, y, t, r  and  d  associated with
c     factorizations
c
c              a(free) * q(free)  = (  0 t )
c                        q(free)  = (  z y )
c
c     when a regular, temporary or artificial constraint is deleted
c     from the working set.
c
c     the  nactiv x nactiv  upper-triangular matrix  t  is stored
c     with its (1,1) element in position  (it,jt)  of the array  t.
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      integer           it, jdel, kdel, lda, ldq, ldt, n, nactiv, nfree,
     *                  ngq, nrz, nz
      logical           unitq

      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n),
     *                  t(ldt,*), work(n)
      integer           kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)

      double precision  cs, sn
      integer           i, ir, itdel, j, jart, jt, k, l, npiv, nrz1,
     *                  nsup

      character*80      rec(4)
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. external subroutines ..
      external          dcopy, dswap, srotgc, sload, scond , suhqr,
     *                  sgesrc

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      jt = nz + 1

      if (jdel.gt.0) then
c
c        regular constraint or temporary bound deleted.
c
         if (jdel.le.n) then
c
c           case 1.  a simple bound has been deleted.
c           =======  columns  nfree+1  and  ir  of gqm' must be swapped.
c
            ir = nz + kdel
c
            itdel = nactiv + 1
            nfree = nfree + 1
            if (nfree.lt.ir) then
               kx(ir) = kx(nfree)
               kx(nfree) = jdel
               call dswap(ngq,gqm(nfree,1),n,gqm(ir,1),n)
            end if
c
            if ( .not. unitq) then
c
c              copy the incoming column of  a(free)  into the end of  t.
c
               do 20 k = 1, nactiv
                  i = kactiv(k)
                  t(nactiv-k+1,nfree) = a(i,jdel)
   20          continue
c
c              expand  q  by adding a unit row and column.
c
               if (nfree.gt.ldq) then
c debug debug 691
                  write (*,*) 'wtf nfree > ldq we are gonna crash'
               else
                  if (nfree.gt.1) then
                     call sload(nfree-1,zero,q(nfree,1),ldq)
                     call sload(nfree-1,zero,q(1,nfree),1)
                  end if
                  q(nfree,nfree) = one
               end if
            end if
         else
c
c           case 2.  a general constraint has been deleted.
c           =======
c
c           delete row  itdel  of  t  and move up the ones below it.
c           t  becomes lower hessenberg.
c
            itdel = kdel
            do 60 k = itdel, nactiv
               j = jt + k - 1
               do 40 l = itdel, k - 1
                  i = it + l - 1
                  t(i,j) = t(i+1,j)
   40          continue
   60       continue
c
            do 80 i = nactiv - itdel + 1, nactiv - 1
               kactiv(i) = kactiv(i+1)
   80       continue
            nactiv = nactiv - 1
         end if
c
         nz = nz + 1
c
         if (nactiv.eq.0) then
            dtmax = one
            dtmin = one
         else

c           restore the nactiv x (nactiv+1) upper-hessenberg matrix  t
c           to upper-triangular form.  the  nsup  super-diagonal
c           elements are removed by a backward sweep of rotations.
c           the rotation for the  (1,1)-th  element of  t  is generated
c           separately.

            nsup = itdel - 1
c
            if (nsup.gt.0) then
               npiv = jt + itdel - 1
               if (nsup.gt.1) then
                  call dcopy(nsup-1,t(it+1,jt+1),ldt+1,s(jt+1),1)
                  call suhqr('right',nactiv,1,nsup,c(jt+1),s(jt+1),
     *                        t(it,jt+1),ldt)
               end if
c
               call srotgc(t(it,jt+1),t(it,jt),cs,sn)
               t(it,jt) = zero
               s(jt) = -sn
               c(jt) = cs
               call sgesrc('right','variable','backwards',nfree,nfree,
     *                     nz,npiv,c,s,q,ldq)
               call sgesrc('left ','variable','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gqm,n)
            end if
c
            jt = jt + 1
            call scond (nactiv,t(it,jt),ldt+1,dtmax,dtmin)
         end if
      end if
c
      nrz1 = nrz + 1
c
      if (nz.gt.nrz) then
         if (jdel.gt.0) then
            jart = nrz1 - 1 + idamax(nz-nrz1+1,gqm(nrz1,1),1)
         else
            jart = -jdel
         end if
c
         if (jart.gt.nrz1) then
c
c           swap columns  nrz1  and  jart  of  q  and  gqm.
c
            if (unitq) then
               k = kx(nrz1)
               kx(nrz1) = kx(jart)
               kx(jart) = k
            else
               call dswap(nfree,q(1,nrz1),1,q(1,jart),1)
            end if
c
            call dswap(ngq,gqm(nrz1,1),n,gqm(jart,1),n)
         end if
      end if
c
      nrz = nrz1

c                                 end of rzdel

      end

      subroutine cmchzr(firstv,n,nclin,istate,bigalf,bigbnd,pnorm,
     *                  hitlow,move,onbnd,unbndd,alfa,alfap,jhit,anorm,
     *                  ap,ax,bl,bu,featol,featlu,p,x)
c----------------------------------------------------------------------
c     cmchzr  finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).
c
c     in this version of cmchzr, when x is infeasible, the number of
c     infeasibilities will never increase.  if the number stays the
c     same, the sum of infeasibilities will decrease.  if the number
c     decreases by one or more,  the sum of infeasibilities will usually
c     decrease also, but occasionally it will increase after the step
c     alfa  is taken.  (convergence is still assured because the number
c     has decreased.)
c
c     three possible steps are computed as follows:
c
c     alfaf = the maximum step that can be taken without violating
c              one of the constraints that are currently satisfied.
c
c     alfai = reaches a linear constraint that is currently violated.
c              usually this will be the furthest such constraint along
c              p, subject to the angle between the constraint normal and
c              p being reasonably close to the maximum value among
c              infeasible constraints,  but if firstv = .true. it will
c              be the first one along p.  the latter case applies only
c              when the problem has been determined to be infeasible,
c              and the sum of infeasibilities are being minimized.
c              (alfai is not defined when x is feasible.)
c
c     alfai is needed occasionally when infeasible, to prevent
c     going unnecessarily far when alfaf is quite large.  it will
c     always come into effect when x is about to become feasible.
c     (the sum of infeasibilities will decrease initially as alfa
c     increases from zero, but may start increasing for larger steps.
c     choosing a large alfai allows several elements of  x  to
c     become feasible at the same time.
c
c     in the end, we take  alfa = alfaf  if x is feasible, or if
c     alfai > alfap (where  alfap  is the perturbed step from pass 1).
c     otherwise,  we take  alfa = alfai.
c
c     input parameters
c     ----------------
c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            if alfa reaches a bound as large as bigbnd, it is
c            classed as an unbounded step.
c     featol is the array of current feasibility tolerances used by
c            cmsinf.  typically in the range 0.5*tolx to 0.99*tolx,
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
c     bu     the upper bounds on ditto.
c     x      the values of       ditto.
c     p      the search direction.
c
c
c     output parameters
c     -----------------
c     hitlow  = true  if a lower bound restricted alfa.
c             = false otherwise.
c     move    = true  if  exact ge stepmn  (defined at end of code).
c     onbnd   = true  if  alfa = exact.  this means that the step  alfa
c                     moves x  exactly onto one of its constraints,
c                     namely  bound.
c             = false if the exact step would be too small
c                     ( exact .lt. stepmn ).
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
c               otherwise,          alfa = max( stepmn, exact ).
c
c
c     cmchzr is based on minos 5.2 routine m5chzr, which implements the
c     expand procedure to deal with degeneracy. the step alfaf is
c     chosen as in the two-pass approach of paula harris (1973), except
c     that this version insists on returning a positive step, alfa.
c     two features make this possible:
c
c        1. featol increases slightly each iteration.
c
c        2. the blocking constraint, when added to the working set,
c           retains the value ax(jhit) + alfa * ap(jhit),
c           even if this is not exactly on the blocking bound.
c
c     for infeasible variables moving towards their bound, we require
c     the rate of change of the chosen constraint to be at least gamma
c     times as large as the biggest available.  this still gives us
c     freedom in pass 2.
c     gamma = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
c     gamma = 0.001 seems to be safe.

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      double precision  gamma
      parameter         (gamma=1.0d-3)

      double precision  alfa, alfap, bigalf, bigbnd, pnorm
      integer           jhit, n, nclin
      logical           firstv, hitlow, move, onbnd, unbndd

      double precision  anorm(*), ap(*), ax(*), bl(n+nclin),
     *                  bu(n+nclin), featlu(n+nclin), featol(n+nclin),
     *                  p(n), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9, tolinc, tolx0
      integer           iprint, isumm, itnfix, kdegen, lines1, lines2,
     *                  ndegen, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg), nfix(2)

      double precision  alfai, atp, atpabs, atpmxf, atpmxi, atpscd, atx,
     *                  biglow, bigupp, bound, delta, exact, res,
     *                  stepmn, tolpiv
      integer           i, j, jhitf, jhiti, js
      logical           blockf, blocki

      character*80      rec(4)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg005/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------

c     tolpiv is a tolerance to exclude negligible elements of a'p.

      biglow = -bigbnd
      bigupp = bigbnd
      tolpiv = epspt9*pnorm

c     first pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

      atpmxi = zero
      alfap = bigalf
c
      do 20 j = 1, n + nclin
         js = istate(j)
c
         if (js.le.0) then
            delta = featol(j)
c
            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = abs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = abs(atp)
               atpscd = atpabs/(one+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

               res = -one

            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

c              test for smaller alfap.
c              if the upper bound is violated. test for bigger atp.
c
               if (bl(j).gt.biglow) then
                  res = atx - bl(j) + delta
                  if (res.lt.alfap*atpabs) alfap = res/atpabs
               end if

               if (js.eq.-1) atpmxi = max(atpmxi,atpscd)

            else if (atp.gt.zero .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfap.
c              if the lower bound is violated. test for bigger atp.
c
               if (bu(j).lt.bigupp) then
                  res = bu(j) - atx + delta

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
         alfai = zero
      end if

      atpmxf = zero
      atpmxi = gamma*atpmxi
      jhitf = 0
      jhiti = 0

      do 40 j = 1, n + nclin
         js = istate(j)

         if (js.le.0) then

            if (j.le.n) then
               atx = x(j)
               atp = p(j)
               atpabs = abs(atp)
               atpscd = atpabs
            else
               i = j - n
               atx = ax(i)
               atp = ap(i)
               atpabs = abs(atp)
               atpscd = atpabs/(one+anorm(i))
            end if

            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

               res = -one

            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing.
c              test for bigger a'p if the lower bound is satisfied.
c              test for smaller alfaf.

               if (atpscd.gt.atpmxf) then

                  if (bl(j).gt.biglow) then
                     res = atx - bl(j)

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

                  if (firstv) then
                     res = atx - bu(j)

                     if (res.le.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if

                  else if (atpscd.ge.atpmxi) then
                     res = atx - bu(j)

                     if (res.gt.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
                  end if
               end if

            else if (atp.gt.zero .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfap.
c
               if (atpscd.gt.atpmxf) then
c
                  if (bu(j).lt.bigupp) then
                     res = bu(j) - atx
c
                     if (res.le.alfap*atp) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                  end if
               end if
c
               if (js.eq.-2) then
c
c                 the lower bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.
c
                  if (firstv) then
                     res = bl(j) - atx
c
                     if (res.le.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  else if (atpscd.ge.atpmxi) then
                     res = bl(j) - atx
c
                     if (res.gt.alfai*atp) then
                        alfai = res/atp
                        jhiti = j
                     end if
                  end if
               end if
            end if
         end if
   40 continue

c     see if a feasible and/or infeasible constraint blocks.

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
         hitlow = atp .lt. zero
      end if

c     if there is a choice between alfaf and alfai, it is probably best
c     to take alfai.  however, we can't if alfai is bigger than alfap.

      if (blocki .and. alfai.le.alfap) then

c        an infeasible variable reaches its violated bound.

         jhit = jhiti
         if (jhit.le.n) then
            atp = p(jhit)
         else
            atp = ap(jhit-n)
         end if
         hitlow = atp .gt. zero
      end if
c
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

      if (hitlow) then
         bound = bl(jhit)
      else
         bound = bu(jhit)
      end if
c
      unbndd = abs(bound) .ge. bigbnd
      if (unbndd) go to 60
c
      stepmn = tolinc*featlu(jhit)/abs(atp)
      exact = (bound-atx)/atp
      alfa = max(stepmn,exact)
      onbnd = alfa .eq. exact
      move = exact .ge. stepmn
      if ( .not. move) ndegen = ndegen + 1
c
      return

c     unbounded.

   60 alfa = bigalf
      move = .true.
      onbnd = .false.
c                                 end of cmchzr
      end

      subroutine cmprt(msglvl,nfree,nrowa,n,nclin,nctotl,bigbnd,named,
     *                  names,nactiv,istate,kactiv,kx,a,bl,bu,c,clamda,
     *                  rlamda,x)
c----------------------------------------------------------------------
c     cmprt   creates the expanded lagrange multiplier vector clamda.
c----------------------------------------------------------------------
      implicit none

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd
      integer           msglvl, n, nactiv, nclin, nctotl, nfree, nrowa
      logical           named

      double precision  a(nrowa,*), bl(nctotl), bu(nctotl), c(*),
     *                  clamda(nctotl), rlamda(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
      integer           ip, is, j, k, nfixed, nplin, nz
      character*8       names(*)

      external          sload
c----------------------------------------------------------------------
      nplin = n + nclin
      nz = nfree - nactiv

c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      call sload(nctotl,zero,clamda,1)
      nfixed = n - nfree
      do 20 k = 1, nactiv + nfixed
         if (k.le.nactiv) j = kactiv(k) + n
         if (k.gt.nactiv) j = kx(nz+k)
         clamda(j) = rlamda(k)
   20 continue
c                                 end of cmprt
      end

      subroutine rzadd(unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,nrz,ngq,n,lda,ldq,ldr,ldt,kx,condmx,drzz,
     *                  a,r,t,gqm,q,w,c,s,msglvl)
c----------------------------------------------------------------------
c     rzadd  updates the matrices  z, y, t, r  and  d  associated with
c     factorizations
c
c              a(free) * q(free)  = (  0 t )
c                        q(free)  = (  z y )
c                      r' *d * r  =   hz
c
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
c     e) if  rset  is  false,   r  is not touched.
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
c     if  ngq .gt. 0,  the column transformations are applied to the
c     columns of the  (ngq x n)  matrix  gqm'.
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  condmx, drzz
      integer           iadd, ifix, inform, it, jadd, lda, ldq, ldr,
     *                  ldt, msglvl, n, nactiv, nfree, ngq, nrz, nz
      logical           rset, unitq

      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), r(ldr,*),
     *                  s(n), t(ldt,*), w(n)
      integer           kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin, epspt3, epspt5, epspt8,
     *                  epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)

      double precision  cond, condbd, dtnew, tdtmax, tdtmin
      integer           i, j, jt, k, nanew, npiv, nsup
      logical           bound, overfl

      character*80      rec(5)
c     .. external functions ..
      double precision  dnrm2, sdiv 
      external          dnrm2, sdiv 
c     .. external subroutines ..
      external          dcopy, dscal, cmqmul, nggnfm, scond , ssrotg,
     *                  smload, sgeapr, suhqr, sutsrh, sgesrc

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      condbd = one/epspt9
c
      overfl = .false.
      bound = jadd .le. n
      jt = nz + 1
c
      if (bound) then

c        a simple bound has entered the working set.  iadd is not used.

         nanew = nactiv
c
         if (unitq) then

c           q is not stored, but  kx  defines an ordering of the columns
c           of the identity matrix that implicitly define q.
c           define the sequence of pairwise interchanges p that moves
c           the newly-fixed variable to position  nfree.
c           reorder  kx  accordingly.

            do 20 i = 1, nfree - 1
               if (i.ge.ifix) then
                  w(i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
   20       continue
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

c        a general constraint has entered the working set.
c        ifix is not used.

         nanew = nactiv + 1

c        transform the incoming row of a by q'.

         call dcopy(n,a(iadd,1),lda,w,1)
         call cmqmul(8,n,nz,nfree,ldq,unitq,kx,w,q,c)

c        check that the incoming row is not dependent upon those
c        already in the working set.

         dtnew = dnrm2(nz,w,1)
         if (nactiv.eq.0) then

c          this is the only general constraint in the working set.

            cond = sdiv (asize,dtnew,overfl)
            tdtmax = dtnew
            tdtmin = dtnew
         else
c
c           there are already some general constraints in the working
c           set.  update the estimate of the condition number.
c
            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = sdiv (tdtmax,tdtmin,overfl)
         end if

         if (cond.gt.condmx .or. overfl) go to 80

         if (unitq) then

c           first general constraint added.  set  q = i.

            call smload('general',nfree,nfree,zero,one,q,ldq)
            unitq = .false.
            it = 0
         end if
      end if
c
      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if
c
      if (unitq) then

c        the orthogonal matrix  q  (i.e.,  q) is not stored explicitly.
c        apply  p, the sequence of pairwise interchanges that moves the
c        newly-fixed variable to position  nfree.

         if (ngq.gt.0) call sgeapr('left','transpose',nfree-1,w,ngq,gqm,
     *                             n)
c
         if (rset) then
c
c           apply the pairwise interchanges to  rz.
c           the subdiagonal elements generated by this process are
c           stored in  s(ifix), s(2), ..., s(nrz-1).
c
            nsup = nrz - ifix
            call nggnfm('right',nrz,ifix,nrz,s,r,ldr)
         end if
      else

c        the matrix  q  is stored explicitly.
c        define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1, 2), (2, 3), ...,
c        (npiv-1, npiv).  the rotations must be applied to q, gqm', r
c        and t.

         call ssrotg('varble','forwrds',npiv-1,w(npiv),w,1,c,s)
c
         if (ngq.gt.0) call sgesrc('left ','variable','forwards',npiv,
     *                             ngq,1,npiv,c,s,gqm,n)
         call sgesrc('right','variable','forwards',nfree,nfree,1,npiv,c,
     *               s,q,ldq)
c
         if (rset) then
c
c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nrz-1).
c
            nsup = nrz - 1
            call sutsrh('right',nrz,1,nrz,c,s,r,ldr)
         end if
      end if
c
      if (rset) then

c        eliminate the  nsup  subdiagonal elements of  r  stored in
c        s(nrz-nsup), ..., s(nrz-1)  with a left-hand sweep of rotations
c        in planes (nrz-nsup, nrz-nsup+1), ..., (nrz-1, nrz).

         call suhqr('left ',nrz,nrz-nsup,nrz,c,s,r,ldr)
c
         if (nsup.gt.0 .and. drzz.ne.one) then
            drzz = c(nrz-1)**2 + drzz*s(nrz-1)**2
         end if
      end if
c
      if ( .not. unitq) then
         if (bound) then

c           bound constraint added.   the rotations affect columns
c           nz+1  thru  nfree  of  gqm'  and  t.

c           the last row and column of  q  has been transformed to plus
c           or minus the unit vector  e(nfree).  we can reconstitute the
c           column of gqm' corresponding to the new fixed variable.
c
            if (w(nfree).lt.zero) then
               if (ngq.gt.0) call dscal(ngq,-one,gqm(nfree,1),n)
            end if
c
            if (nactiv.gt.0) then
               t(it,jt-1) = s(jt-1)*t(it,jt)
               t(it,jt) = c(jt-1)*t(it,jt)
c
               if (nactiv.gt.1) then
                  call sutsrh('right',nactiv,1,nactiv,c(jt),s(jt),
     *                        t(it,jt),ldt)
                  call dcopy(nactiv-1,s(jt),1,t(it+1,jt),ldt+1)
               end if
c
               jt = jt - 1
               call scond (nactiv,t(it,jt),ldt+1,tdtmax,tdtmin)
               cond = sdiv (tdtmax,tdtmin,overfl)
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
            call dcopy(nanew,w(jt),1,t(it,jt),ldt)
         end if
      end if
c

c     prepare to exit.  check the magnitude of the condition estimator.

   80 if (nanew.gt.0) then
         if (cond.lt.condmx .and. .not. overfl) then
c
c           the factorization has been successfully updated.
c
            inform = 0
            dtmax = tdtmax
            dtmin = tdtmin

         else
c
c           the proposed working set appears to be linearly dependent.
c
            inform = 1

         end if
      end if

c                                 end of rzadd

      end

      subroutine npmrt (feasqp,n,nclin,ncnln,objalf,grdalf,qpcurv,
     *                  istate,cjdx,cmul,cs,dlam,rho,violn,work1,work2)
c----------------------------------------------------------------------
c     npmrt    computes the value and directional derivative of the
c     augmented lagrangian merit function.  the penalty parameters
c     rho(j) are boosted if the directional derivative of the resulting
c     augmented lagrangian function is not sufficiently negative.  if
c     rho needs to be increased,  the perturbation with minimum two-norm
c     is found that gives a directional derivative equal to  - p'hp.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two
      parameter         (two=2.0d+0)

      double precision  grdalf, objalf, qpcurv
      integer           n, nclin, ncnln
      logical           feasqp

      double precision  cjdx(*), cmul(*), cs(*), dlam(*), rho(*),
     *                  violn(*), work1(*), work2(*)
      integer           istate(*)
c     .. scalars in common ..
      double precision  rhodmp, rhomax, rhonrm, scale
      integer           iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  pterm, pterm2, qnorm, rho1, rhoi, rhomin,
     *                  rhonew, rtmin, tscl
      integer           i, l, nplin
      logical           boost, overfl

      character*80      rec(3)
c     .. external functions ..
      double precision  ddot, dnrm2, sdiv 
      external          ddot, dnrm2, sdiv 
c     .. external subroutines ..
      external          dcopy, sdscl 

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
      if (ncnln.eq.0) return
c
      rtmin = wmach(6)
c
      objalf = objalf - ddot(ncnln,cmul,1,cs,1)
      grdalf = grdalf - ddot(ncnln,dlam,1,cs,1)
c
      call dcopy(ncnln,cs,1,work1,1)
c
      if ( .not. feasqp) then
         nplin = n + nclin
c
         do 20 i = 1, ncnln
            if (istate(nplin+i).lt.0 .or. violn(i).ne.zero) work1(i)
     *          = -cjdx(i)
   20    continue
      end if
c
      grdalf = grdalf + ddot(ncnln,work1,1,cmul,1)

      if (feasqp) then
c
c        find the quantities that define  rhomin, the vector of minimum
c        two-norm such that the directional derivative is one half of
c        approximate curvature   - (dx)'h(dx).
c
         do 40 i = 1, ncnln
            if (abs(cs(i)).le.rtmin) then
               work2(i) = zero
            else
               work2(i) = cs(i)**2
            end if
   40    continue
c
         qnorm = dnrm2(ncnln,work2,1)
         tscl = sdiv (grdalf+half*qpcurv,qnorm,overfl)
         if (abs(tscl).le.rhomax .and. .not. overfl) then

c           bounded  rhomin  found.  the final value of  rho(j)  will
c           never be less than  rhomin(j).  if the  qp  was feasible,  a
c           trial value  rhonew  is computed that is equal to the
c           geometric mean of the previous  rho  and a damped value of
c           rhomin.  the new  rho  is defined as  rhonew  if it is less
c           than half the previous  rho  and greater than  rhomin.

            scale = one
            do 60 i = 1, ncnln
               rhomin = max((work2(i)/qnorm)*tscl,zero)
               rhoi = rho(i)
c
               rhonew = sqrt(rhoi*(rhodmp+rhomin))
               if (rhonew.lt.half*rhoi) rhoi = rhonew
               if (rhoi.lt.rhomin) rhoi = rhomin
               rho(i) = rhoi
   60       continue
c
            rho1 = rhonrm
            rhonrm = dnrm2(ncnln,rho,1)
c

c           if  incrun = .true.,  there has been a run of iterations in
c           which the norm of  rho  has not decreased.  conversely,
c           incrun = false  implies that there has been a run of
c           iterations in which the norm of rho has not increased.  if
c           incrun changes during this iteration the damping parameter
c           rhodmp is increased by a factor of two.  this ensures that
c           rho(j) will oscillate only a finite number of times.

            boost = .false.
            if (incrun .and. rhonrm.lt.rho1) boost = .true.
            if ( .not. incrun .and. rhonrm.gt.rho1) boost = .true.
            if (boost) then
               rhodmp = two*rhodmp
               incrun = .not. incrun
            end if
         end if

      else
c
c        the  qp  was infeasible.  do not alter the penalty parameters,
c        but compute the scale factor so that the constraint violations
c        are reduced.
c
         call sdscl (ncnln,rho,1,work1,1)
         pterm2 = ddot(ncnln,work1,1,cs,1)
c
         scale = rhomax
         tscl = sdiv (grdalf,pterm2,overfl)
         if (tscl.gt.scale .and. tscl.le.rhomax/(one+rhonrm)
     *       .and. .not. overfl) scale = tscl
c
         call dcopy(ncnln,cs,1,work1,1)
      end if
c

c     compute the new value and directional derivative of the
c     merit function.

      call sdscl (ncnln,rho,1,work1,1)
c
      pterm = ddot(ncnln,work1,1,cs,1)
      objalf = objalf + half*scale*pterm
c
      if (feasqp) pterm2 = pterm
c
      grdalf = grdalf - scale*pterm2

c                                 end of npmrt

      end

      subroutine chfd(inform,msglvl,lvlder,n,ncnln,ldcj,ldcju,bigbnd,
     *                  epsrf,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,c2,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,y,
     *                  iuser,user)
c----------------------------------------------------------------------
c     chfd  computes difference intervals for the missing gradients of
c     f(x) and c(x). intervals are computed using a procedure that
c     usually requires about two function evaluations if the function
c     is well scaled.  central-difference gradients are obtained as a
c     by-product of the algorithm.
c
c     on entry...
c     objf and c contain the problem functions at the point x.
c     an element of cjac or grad not equal to rdummy signifies a known
c     gradient value.  such values are not estimated by differencing.
c     cjacu and gradu have dummy elements in the same positions as
c     cjac and gradu.
c
c     on exit...
c     cjac and grad contain central-difference derivative estimates.
c     elements of cjacu and gradu are unaltered except for those
c     corresponding to constant derivatives, which are given the same
c     values as cjac or grad.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  factor
      parameter         (factor=0.97d+0)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two, four, ten
      parameter         (two=2.0d+0,four=4.0d+0,ten=1.0d+1)

      double precision  bigbnd, epsrf, fdnorm, objf
      integer           inform, ldcj, ldcju, lvlder, msglvl, n, ncnln

      double precision  bl(n), bu(n), c(*), c1(*), c2(*), cjac(ldcj,*),
     *                  cjacu(ldcju,*), grad(n), gradu(n), hcntrl(*),
     *                  hforwd(*), user(*), x(n), y(n)
      integer           iuser(*), needc(*)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  ncdiff, nfdiff, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  biglow, bigupp, cdest, cjdiff, d, dx, epsa,
     *                  errbnd, errmax, errmin, f1, f2, fdest, fx,
     *                  gdiff, h, hcd, hfd, hmax, hmin, hopt, hphi,
     *                  objf2, sdest, signh, stepbl, stepbu, sumeps,
     *                  sumsd, test, xj, yj
      integer           i, info, irow1, irow2, iter, itmax, j, mode,
     *                  nccnst, ncolj, nfcnst, nstate
      logical           debug, done, first, headng, needed

      character*80      rec(4)
c     .. external subroutines ..
      external          chcore, iload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
      inform = 0
      needed = lvlder .eq. 0 .or. lvlder .eq. 2 .or. lvlder .eq. 1 .and.
     *         ncnln .gt. 0
      if ( .not. needed) return
c
      debug = npdbg .and. inpdbg(5) .gt. 0
      if (lfdset.eq.0) then
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
c

c        for each column of the jacobian augmented by the transpose of
c        the objective gradient, rows irow1 thru irow2 are searched for
c        missing elements.

         irow1 = 1
         irow2 = ncnln + 1
         if (lvlder.eq.1) irow2 = ncnln
         if (lvlder.eq.2) irow1 = ncnln + 1
c
         biglow = -bigbnd
         bigupp = bigbnd
c
         if (ncnln.gt.0) call iload (ncnln,(0),needc,1)
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
            do 40 i = irow1, irow2
c
               if (i.le.ncnln) then
                  test = cjacu(i,j)
               else
                  test = gradu(j)
               end if
c
               if (test.eq.rdummy) then
c                 ======================================================
c                 get the difference interval for this element.
c                 ======================================================
                  ncolj = ncolj + 1
c
                  if (i.le.ncnln) then
                     needc(i) = 1
                     fx = c(i)
                     epsa = epsrf*(one+abs(c(i)))
                  else
                     fx = objf
                     epsa = epsrf*(one+abs(fx))
                  end if
c
c                 ------------------------------------------------------
c                 find a finite-difference interval by iteration.
c                 ------------------------------------------------------
                  iter = 0
                  hopt = two*(one+abs(xj))*sqrt(epsrf)
                  h = signh*ten*hopt
                  cdest = zero
                  sdest = zero
                  first = .true.
c
c                 +                repeat
   20             x(j) = xj + h
                  if (i.le.ncnln) then
                     call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                           nstate,iuser,user)
                     if (mode.lt.0) go to 200
                     f1 = c1(i)
                  else
                     call objfun(mode,n,x,f1,gradu,nstate,iuser,user)
                     if (mode.lt.0) go to 200
                  end if
c
                  x(j) = xj + h + h
                  if (i.le.ncnln) then
                     call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                           nstate,iuser,user)
                     if (mode.lt.0) go to 200
                     f2 = c1(i)
                  else
                     call objfun(mode,n,x,f2,gradu,nstate,iuser,user)
                     if (mode.lt.0) go to 200
                  end if
c
                  call chcore(debug,done,first,epsa,epsrf,fx,info,iter,
     *                        itmax,cdest,fdest,sdest,errbnd,f1,f2,h,
     *                        hopt,hphi)
c
c                 +                until     done
                  if ( .not. done) go to 20
c
                  if (i.le.ncnln) then
                     cjac(i,j) = cdest
                     if (info.eq.1 .or. info.eq.2) then
                        nccnst = nccnst + 1
                        ncdiff = ncdiff - 1
                        cjacu(i,j) = -rdummy
                     end if
                  else
                     grad(j) = cdest
                     if (info.eq.1 .or. info.eq.2) then
                        nfcnst = nfcnst + 1
                        nfdiff = nfdiff - 1
                        gradu(j) = -rdummy
                     end if
                  end if
c
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
   40       continue
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
                  hfd = two*sqrt(sumeps/sumsd)
                  errmax = two*sqrt(sumeps*sumsd)
               end if
c
               if (hcd.eq.zero) hcd = ten*hfd
c

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
               if (bu(j).lt.bigupp .and. bu(j).gt.bl(j))
     *             stepbu = min(stepbu,bu(j)-xj)
c
               if (half*(stepbl+stepbu).lt.zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if
c
               d = factor*d
   80       continue
c
            if (ncnln.gt.0) then
               call iload (ncnln,(1),needc,1)
               call confun(mode,ncnln,n,ldcju,needc,y,c2,cjacu,nstate,
     *                     iuser,user)
               if (mode.lt.0) go to 200
            end if
c
            call objfun(mode,n,y,objf2,gradu,nstate,iuser,user)
            if (mode.lt.0) go to 200
c

c           loop over each of the elements of  x.

            do 140 j = 1, n
               yj = y(j)
               dx = half*(x(j)-yj)
               y(j) = yj + dx
c
               if (ncnln.gt.0) then
                  ncolj = 0
                  do 100 i = 1, ncnln
                     if (cjacu(i,j).eq.-rdummy) then
                        needc(i) = 1
                        ncolj = ncolj + 1
                     else
                        needc(i) = 0
                     end if
  100             continue
c
                  if (ncolj.gt.0) then
                     call confun(mode,ncnln,n,ldcju,needc,y,c1,cjacu,
     *                           nstate,iuser,user)
                     if (mode.lt.0) go to 200
c
                     do 120 i = 1, ncnln
                        if (needc(i).eq.1) then
                           cjdiff = (c1(i)-c2(i))/dx
                           if (cjdiff.eq.cjac(i,j)) then
                              cjacu(i,j) = cjdiff
                           else
                              cjacu(i,j) = rdummy
                              nccnst = nccnst - 1
                              ncdiff = ncdiff + 1
                           end if
                        end if
  120                continue
                  end if
               end if
c
c              check the objective gradient element.
c
               if (gradu(j).eq.-rdummy) then
c
                  call objfun(mode,n,y,f1,gradu,nstate,iuser,user)
                  if (mode.lt.0) go to 200
c
                  gdiff = (f1-objf2)/dx
                  if (gdiff.eq.grad(j)) then
                     gradu(j) = gdiff
                  else
                     gradu(j) = rdummy
                     nfdiff = nfdiff + 1
                     nfcnst = nfcnst - 1
                  end if
               end if
c
               y(j) = yj
  140       continue

c
            if (ncdiff.eq.0 .and. lvlder.lt.2) then
               if (lvlder.eq.0) lvlder = 2
               if (lvlder.eq.1) lvlder = 3

            end if
c
            if (nfdiff.eq.0 .and. lvlder.ne.1) then
               if (lvlder.eq.0) lvlder = 1
               if (lvlder.eq.2) lvlder = 3

            end if
         end if
      else if (lfdset.eq.2) then
c
c        the user has supplied hforwd and hcntrl.
c        check for wild values.
c
         do 160 j = 1, n
            if (hforwd(j).le.zero) then
               hforwd(j) = epspt5
            end if
  160    continue
         do 180 j = 1, n
            if (hcntrl(j).le.zero) then
               hcntrl(j) = epspt3
            end if
  180    continue
      end if
c
      return
c
  200 inform = mode

c                                 end of chfd

      end

      subroutine npfeas(n,nclin,ncnln,istate,bigbnd,cvnorm,errmax,jmax,
     *                  nviol,ax,bl,bu,c,featol,x,work)
c----------------------------------------------------------------------
c     npfeas  computes the following...
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd, cvnorm, errmax
      integer           jmax, n, nclin, ncnln, nviol

      double precision  ax(*), bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     *                  c(*), featol(n+nclin+ncnln),
     *                  work(n+nclin+ncnln), x(n)
      integer           istate(n+nclin+ncnln)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  biglow, bigupp, con, feasj, res, tolj
      integer           is, j

      character*80      rec(2)
c     .. external functions ..
      double precision  dnrm2
      integer           idamax
      external          dnrm2, idamax

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
c
      biglow = -bigbnd
      bigupp = bigbnd
c

c     compute nviol, the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint
c     violations and residuals of the constraints in the qp working set.

      nviol = 0
c
      do 40 j = 1, n + nclin + ncnln
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
         else
c
c           nonlinear constraint.
c
            con = c(j-n-nclin)
            tolj = zero
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
         else if (is.eq.1 .or. is.le.-2) then
            res = bl(j) - con
         else if (is.ge.2 .or. is.eq.-1) then
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
      jmax = idamax(n+nclin+ncnln,work,1)
      errmax = abs(work(jmax))

      cvnorm = dnrm2(n+nclin+ncnln,work,1)

c                                 end of npfeas

      end

      subroutine npsrch(needfd,inform,n,ncnln,ldcj,ldcju,nfun,ngrad,
     *                  needc,confun,objfun,alfa,alfbnd,alfmax,alfsml,
     *                  dxnorm,epsrf,eta,gdx,grdalf,glf1,glf,objf,
     *                  objalf,qpcurv,xnorm,c,c2,cjac,cjacu,cjdx,cjdx2,
     *                  cmul1,cmul,cs1,cs,dx,dlam,dslk,grad,gradu,qpmul,
     *                  rho,slk1,slk,x1,x,work,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     npsrch finds the steplength alfa that gives sufficient decrease in
c     the augmented lagrangian merit function.
c
c     on exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
c     with an associated merit function value  objalf  which is lower
c     than that at the base point. if  inform = 4, 5, 6, 7 or 8,  alfa
c     is zero and  objalf will be the merit value at the base point.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two
      parameter         (two=2.0d+0)
      double precision  tolg
      parameter         (tolg=1.0d-1)
      double precision  rmu
      parameter         (rmu=1.0d-4)

      double precision  alfa, alfbnd, alfmax, alfsml, dxnorm, epsrf,
     *                  eta, gdx, glf, glf1, grdalf, objalf, objf,
     *                  qpcurv, xnorm
      integer           inform, ldcj, ldcju, lenw, n, ncnln, nfun, ngrad
      logical           needfd

      double precision  c(*), c2(*), cjac(ldcj,*), cjacu(ldcju,*),
     *                  cjdx(*), cjdx2(*), cmul(*), cmul1(*), cs(*),
     *                  cs1(*), dlam(*), dslk(*), dx(n), grad(n),
     *                  gradu(n), qpmul(*), rho(*), slk(*), slk1(*),
     *                  user(*), w(lenw), work(*), x(n), x1(n)
      integer           iuser(*), needc(*)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9, rhodmp, rhomax,
     *                  rhonrm, scale
      integer           iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  ncdiff, nfdiff, nout
      logical           incrun, npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  alfbst, cs1jdx, csjdx, curvc, curvlf, epsaf,
     *                  epsmch, fbest, fterm, ftry, g0, gbest, gtry,
     *                  oldf, oldg, q, rhobfs, s, t, targtg, tgdx, tglf,
     *                  tobj, tobjm, tolabs, tolax, tolrel, tolrx,
     *                  toltny
      integer           j, maxf, mode, nstate, numf
      logical           debug, done, first, imprvd

      character*80      rec(3)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, srchq, srchc , iload ,
     *                  sdscl , smcopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /npdebg/inpdbg, npdbg

      epsmch = wmach(3)

      if ( .not. needfd .and. ncnln.gt.0) cs1jdx = ddot(ncnln,cs1,1,
     *    cjdx,1)


c     set the input parameters and tolerances for srchc  and srchq.

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
      debug = npdbg .and. inpdbg(4) .gt. 0
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
c
      if (ncnln.gt.0) call iload (ncnln,(1),needc,1)
c
      if (needfd) then
         mode = 0
      else
         mode = 2
      end if
c
      first = .true.
c

c     commence main loop, entering srchc  or srchq two or more times.
c     first = .true. for the first entry, .false. for subsequent entries
c     done  = .true. indicates termination, in which case the value of
c     inform gives the result of the search.
c     inform = 1 if the search is successful and alfa < alfmax.
c            = 2 if the search is successful and alfa = alfmax.
c            = 3 if a better point was found but too many functions
c                were needed (not sufficient decrease).
c            = 4 if alfmax < tolabs (too small to do a search).
c            = 5 if alfa < alfsml (srchq only -- maybe want to switch
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
         call srchq(first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest)
      else
         call srchc (first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest,gbest)
      end if
c
      if (imprvd) then
         objf = tobj
         objalf = tobjm
c
         if (ncnln.gt.0) call dcopy(ncnln,c2,1,c,1)
c
         if ( .not. needfd) then
            call dcopy(n,gradu,1,grad,1)
            gdx = tgdx
            glf = tglf
c
            if (ncnln.gt.0) then
               call dcopy(ncnln,cjdx2,1,cjdx,1)
               call smcopy('general',ncnln,n,cjacu,ldcju,cjac,ldcj)
            end if
         end if
      end if
c

c     if done = .false.,  the problem functions must be computed for
c     the next entry to srchc  or srchq.
c     if done = .true.,   this is the last time through.

      if ( .not. done) then
c
         call dcopy(n,x1,1,x,1)
         call daxpy(n,alfa,dx,1,x,1)
c
         if (ncnln.gt.0) then
c
c           compute the new estimates of the multipliers and slacks.
c           if the step length is greater than one,  the multipliers
c           are fixed as the qp-multipliers.
c
            if (alfa.le.one) then
               call dcopy(ncnln,cmul1,1,cmul,1)
               call daxpy(ncnln,alfa,dlam,1,cmul,1)
            end if
            call dcopy(ncnln,slk1,1,slk,1)
            call daxpy(ncnln,alfa,dslk,1,slk,1)
c

c           compute the new constraint vector and jacobian.

            call confun(mode,ncnln,n,ldcju,needc,x,c2,cjacu,nstate,
     *                  iuser,user)
            if (mode.lt.0) go to 60
c
            call dcopy(ncnln,c2,1,cs,1)
            call daxpy(ncnln,(-one),slk,1,cs,1)
c
            call dcopy(ncnln,cs,1,work,1)
            call sdscl (ncnln,rho,1,work,1)
c
            fterm = ddot(ncnln,cmul,1,cs,1) - half*scale*ddot(ncnln,
     *              work,1,cs,1)
         end if
c
c        ------------------------------------------------------------
c        compute the value and gradient of the objective function.
c        ------------------------------------------------------------
         call objfun(mode,n,x,tobj,gradu,nstate,iuser,user)
         if (mode.lt.0) go to 60
c
         if (ncnln.gt.0) then
            tobjm = tobj - fterm
         else
            tobjm = tobj
         end if
c
         ftry = tobjm - oldf - rmu*oldg*alfa
c
c        ------------------------------------------------------------
c        compute auxiliary gradient information.
c        ------------------------------------------------------------
         if ( .not. needfd) then
            gtry = ddot(n,gradu,1,dx,1)
            tgdx = gtry
            tglf = gtry
            if (ncnln.gt.0) then
c
c              compute the jacobian times the search direction.
c
               call dgemv('n',ncnln,n,one,cjacu,ldcju,dx,1,zero,cjdx2,1)
c
               call dcopy(ncnln,cjdx2,1,work,1)
               call daxpy(ncnln,(-one),dslk,1,work,1)
c
               gtry = gtry - ddot(ncnln,cmul,1,work,1)
               if (alfa.le.one) gtry = gtry - ddot(ncnln,dlam,1,cs,1)
c
               call sdscl (ncnln,rho,1,work,1)
               gtry = gtry + scale*ddot(ncnln,work,1,cs,1)
c
               tglf = tgdx - ddot(ncnln,cjdx2,1,qpmul,1)
c
c              ------------------------------------------------------
c              if alfbnd .le. alfa .lt. alfmax and the norm of the
c              quasi-newton update is bounded, set alfmax to be alfa.
c              this will cause the linesearch to stop if the merit
c              function is decreasing at the boundary.
c              ------------------------------------------------------
               if (alfbnd.le.alfa .and. alfa.lt.alfmax) then
c
                  csjdx = ddot(ncnln,cs,1,cjdx2,1)

                  curvlf = tglf - glf1
                  curvc = abs(csjdx-cs1jdx)
                  rhobfs = max(qpcurv*tolg-curvlf,zero)
                  if (rhobfs.le.curvc*rhomax) then
                     alfmax = alfa
                  else
                     alfbnd = min(two*alfa,alfmax)
                  end if

               end if
            end if
c
            gtry = gtry - rmu*oldg
c
         end if
      end if
c     +    until (      done)
      if ( .not. done) go to 40
c
      nfun = nfun + numf
      if ( .not. needfd) ngrad = ngrad + numf
      alfa = alfbst
c
      if ( .not. imprvd) then
         call dcopy(n,x1,1,x,1)
         call daxpy(n,alfa,dx,1,x,1)
         if (ncnln.gt.0) then
            if (alfa.le.one) then
               call dcopy(ncnln,cmul1,1,cmul,1)
               call daxpy(ncnln,alfa,dlam,1,cmul,1)
            end if
            call dcopy(ncnln,slk1,1,slk,1)
            call daxpy(ncnln,alfa,dslk,1,slk,1)
            call dcopy(ncnln,c,1,cs,1)
            call daxpy(ncnln,(-one),slk,1,cs,1)
         end if
      end if

      return

c     the user wants to stop.  who am i to object?

   60 inform = mode

c                                 end of npsrch

      end

      subroutine npupdt(lsumry,unitq,n,ncnln,nfree,nz,ldcj1,ldcj2,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,cjac1,cjac2,cjdx1,
     *                  cjdx2,cs1,cs2,gq1,gq2,hpq,rpq,qpmul,r,omega,zy,
     *                  wrk1,wrk2)
c----------------------------------------------------------------------
c     npupdt  computes the bfgs update for the approximate hessian of
c     the lagrangian.  if the approximate curvature of the lagrangian
c     function is negative,  a nonnegative penalty vector omega(i) of
c     minimum two norm is computed such that the approximate curvature
c     of the augmented lagrangian will be positive. if no finite penalty
c     vector exists,  the bfgs update is doed with the approximate
c     curvature modified to be a small positive value.
c
c     on entry,  gq1 and gq2 contain the transformed objective gradients
c     at x1 and x2,  hpq contains  r'r(pq), the transformed hessian
c     times the transformed search direction.  the vectors gq1 and hpq
c     are not saved.  if the regular bfgs quasi-newton update could not
c     be done, the first character of lsumry is loaded with 'm'.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      double precision  tolg
      parameter         (tolg=1.0d-1)

      double precision  alfa, glf1, glf2, qpcurv
      integer           ldcj1, ldcj2, ldr, ldzy, n, ncnln, nfree, nz
      logical           unitq
      character*5       lsumry

      double precision  cjac1(ldcj1,*), cjac2(ldcj2,*), cjdx1(*),
     *                  cjdx2(*), cs1(*), cs2(*), gq1(n), gq2(n),
     *                  hpq(n), omega(*), qpmul(*), r(ldr,*), rpq(n),
     *                  wrk1(n+ncnln), wrk2(n), zy(ldzy,*)
      integer           kx(n)
c     .. scalars in common ..
      double precision  drmax, drmin, rcndbd, rfrobn, rhodmp, rhomax,
     *                  rhonrm, scale
      integer           iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  beta, curvl, eta, qi, qmax, qnorm, rtgtp, rtyts,
     *                  test, tinycl, trace1, trace2
      integer           i, imax, j
      logical           overfl, ssbfgs

      character*80      rec(3)
c     .. external functions ..
      double precision  dnrm2, sdiv 
      integer           idamax
      external          dnrm2, sdiv , idamax
c     .. external subroutines ..
      external          daxpy, dgemv, dscal, cmr1md, cmqmul, sload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ngg018/rcndbd, rfrobn, drmax, drmin
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
      if (ncnln.gt.0) call sload(ncnln,zero,omega,1)

c     set curvl = (g2 - g1)'dx,  the approximate curvature along dx of
c     the (augmented) lagrangian.  at first, the curvature is not scaled
c     by the steplength alfa.

      curvl = glf2 - glf1
      tinycl = qpcurv*tolg
      ssbfgs = curvl .le. alfa*tinycl

c     test if curvl is sufficiently positive.  if there are no nonlinear
c     constraints,  no update can be doed.

      if (curvl.lt.tinycl) then
         lsumry(1:1) = 'modified bfgs'
         if (ncnln.gt.0) then
            qmax = zero
            do 20 i = 1, ncnln
               qi = cjdx2(i)*cs2(i) - cjdx1(i)*cs1(i)
               qmax = max(qmax,qi)
               if (qi.le.zero) wrk1(i) = zero
               if (qi.gt.zero) wrk1(i) = qi
   20       continue
c
            qnorm = dnrm2(ncnln,wrk1,1)
c
            test = max(tinycl-curvl,zero)
            beta = sdiv (qmax*test,qnorm*qnorm,overfl)
            if (beta.lt.rhomax .and. .not. overfl) then
               lsumry(1:1) = ' '
               beta = test/(qnorm*qnorm)
               do 40 i = 1, ncnln
                  qi = wrk1(i)
                  omega(i) = beta*qi
                  curvl = curvl + beta*qi*qi
   40          continue
            end if
         end if
      end if
c

c     compute the difference in the augmented lagrangian gradient.

c     update gq1 to include the augmented lagrangian terms.
c
      if (ncnln.gt.0) then
c
         do 80 i = 1, ncnln
            wrk1(i) = -qpmul(i) + omega(i)*cs1(i)
   80    continue
         call dgemv('t',ncnln,n,one,cjac1,ldcj1,wrk1,1,zero,wrk2,1)
c
         do 100 i = 1, ncnln
            wrk1(i) = qpmul(i) - omega(i)*cs2(i)
  100    continue
         call dgemv('t',ncnln,n,one,cjac2,ldcj2,wrk1,1,one,wrk2,1)
c
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,wrk2,zy,wrk1)
         call daxpy(n,one,wrk2,1,gq1,1)
      end if

      if (curvl.lt.tinycl) curvl = tinycl
c
      do 120 j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
  120 continue
c
      rtgtp = sqrt(qpcurv)
      rtyts = sqrt(alfa*curvl)
      eta = one
      if (ssbfgs) eta = rtyts/(rtgtp*alfa)
c
      trace1 = dnrm2(n,hpq,1)/rtgtp
      trace2 = dnrm2(n,wrk2,1)/(rtyts*eta)
      rfrobn = eta*sqrt(abs((rfrobn-trace1)*(rfrobn+trace1)+trace2**2))
c

c     update the cholesky factor of  q'hq.

c     normalize the vector  rpq ( = r(pq) ).
c
      call dscal(n,(one/rtgtp),rpq,1)
c
c     do the self-scaled or regular bfgs update.
c     form the vector wrk1 = gamma * (gq2 - gq1) - beta * r'r*pq,
c     where  gamma = 1/sqrt( curv ) = 1/sqrt( (gq2 - gq1)'sq )
c
      call dscal(n,(one/rtgtp),hpq,1)
c
      if (ssbfgs) then
         do 140 j = 1, n
            call dscal(j,eta,r(1,j),1)
            wrk1(j) = wrk2(j)/rtyts - eta*hpq(j)
  140    continue
      else
         do 160 j = 1, n
            wrk1(j) = wrk2(j)/rtyts - hpq(j)
  160    continue
      end if

c     do the update to  r = r + rpq*wrk1'.
c     rpq is overwritten. arrays gq1 and hpq are used to store the
c     sines and cosines defined by the plane rotations.

      call cmr1md(n,0,n,ldr,n,n,r,hpq,rpq,wrk1,gq1,hpq)
c                                 end of npupdt
      end

      subroutine cmalf(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfa,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  anorm,ap,ax,bl,bu,featol,p,x)
c----------------------------------------------------------------------
c     cmalf finds a step alfa such that the point x + alfa*p reaches
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

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  alfa, atphit, bigalf, bigbnd, palfa, pnorm
      integer           inform, jadd, n, nctotl, numinf
      logical           firstv, hitlow

      double precision  anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer           istate(nctotl)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)

      double precision  absatp, alfa1, alfa2, apmax1, apmax2, atp, atp1,
     *                  atp2, atx, palfa1, palfa2, res, rownrm
      integer           i, j, jadd1, jadd2, js, jsave1, jsave2
      logical           hlow1, hlow2, lastv, negstp, step2

      character*120     rec(4)
c     .. external subroutines ..
      external          cmalf1

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------

      inform = 0


c     first pass -- find steps to perturbed constraints, so that
c     palfa1 will be slightly larger than the true step, and
c     palfa2 will be slightly smaller than it should be.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

c
      negstp = .false.
      call cmalf1(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,palfa1,
     *            palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,featol,p,x)
c
      jsave1 = jadd1
      jsave2 = jadd2
c

c     second pass -- recompute step-lengths without perturbation.
c     amongst constraints that are less than the perturbed steps,
c     choose the one (of each type) that makes the largest angle
c     with the search direction.


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
            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing.

c              the lower bound is satisfied.  test for smaller alfa1.
c
               absatp = -atp
               if (bl(j).gt.(-bigbnd)) then
                  res = atx - bl(j)
                  if (palfa1*absatp.ge.res .or. j.eq.jsave1) then
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
                  if ((firstv .and. palfa2*absatp.ge.res .or.
     *                lastv .and. palfa2*absatp.le.res)
     *                .or. j.eq.jsave2) then
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
            else if (atp.gt.zero .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfa1.
c
               if (bu(j).lt.bigbnd) then
                  res = bu(j) - atx
                  if (palfa1*atp.ge.res .or. j.eq.jsave1) then
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
                  if ((firstv .and. palfa2*atp.ge.res .or. lastv .and.
     *                palfa2*atp.le.res) .or. j.eq.jsave2) then
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

         end if
   20 continue
c

c     determine alfa, the step to be taken.

c     in the infeasible case, check whether to take the step alfa2
c     rather than alfa1...
c
      step2 = numinf .gt. 0 .and. jadd2 .gt. 0
c
c     we do so if alfa2 is less than alfa1 or (if firstv is false)
c     lies in the range  (alfa1, palfa1)  and has a smaller value of
c     atp.
c
      step2 = step2 .and. (alfa2.lt.alfa1 .or. lastv .and. alfa2.le.
     *        palfa1 .and. apmax2.ge.apmax1)
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
            call cmalf1(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)

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

      if (alfa.ge.bigalf) inform = 3
c                                 end of cmalf
      end

      subroutine lssetx(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,
     *                  n,nctotl,ldzy,lda,ldr,ldt,istate,kactiv,kx,jmax,
     *                  errmax,ctx,xnorm,a,ax,bl,bu,cq,res,res0,featol,
     *                  r,t,x,zy,p,work)
c----------------------------------------------------------------------
c     lssetx  computes the point on a working set that is closest to the
c     input vector  x  (in the least-squares sense).  the norm of  x,
c     the transformed residual vector  pr - rq'x,  and the constraint
c     values ax  are also initialized.
c
c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable  rowerr
c     is set.

      integer           ldbg
      parameter         (ldbg=5)
      integer           ntry
      parameter         (ntry=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  ctx, errmax, xnorm
      integer           jmax, lda, ldr, ldt, ldzy, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nz
      logical           linobj, rowerr, unitq

      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl), cq(*),
     *                  featol(nctotl), p(n), r(ldr,*), res(*), res0(*),
     *                  t(ldt,*), work(nctotl), x(n), zy(ldzy,*)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  bnd
      integer           i, is, j, k, ktry

      character*80      rec(2)
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           idamax
      external          ddot, dnrm2, idamax
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dtrmv, cmtsol, cmqmul,
     *                  sload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
c

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
            work(i) = bnd - ddot(n,a(k,1),lda,x,1)
   60    continue
c
         call cmtsol(1,ldt,nactiv,t(1,nz+1),work)
         call sload(n,zero,p,1)
         call dcopy(nactiv,work,1,p(nz+1),1)
c
         call cmqmul(2,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
         call daxpy(n,one,p,1,x,1)
      end if
c

c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.

      xnorm = dnrm2(n,x,1)
      if (nclin.gt.0) call dgemv('n',nclin,n,one,a,lda,x,1,zero,ax,1)
c

c     check the row residuals.

      if (nactiv.gt.0) then
         do 80 k = 1, nactiv
            i = kactiv(k)
            j = n + i
            is = istate(j)
            if (is.eq.1) work(k) = bl(j) - ax(i)
            if (is.ge.2) work(k) = bu(j) - ax(i)
   80    continue
c
         jmax = idamax(nactiv,work,1)
         errmax = abs(work(jmax))
      end if
c
      ktry = ktry + 1
c     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      if ( .not. (errmax.le.featol(jmax) .or. ktry.gt.ntry)) go to 40
c
      rowerr = errmax .gt. featol(jmax)
c

c     compute the linear objective value  c'x  and the transformed
c     residual  pr  -  rq'x = res0  -  rq'x.

      if (nrank.gt.0 .or. linobj) then
         call dcopy(n,x,1,p,1)
         call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
      end if
c
      ctx = zero
      if (linobj) ctx = ddot(n,cq,1,p,1)
c
      if (nrank.gt.0) then
         call dtrmv('u','n','n',nrank,r,ldr,p,1)
         if (nrank.lt.n) call dgemv('n',nrank,n-nrank,one,r(1,nrank+1),
     *                              ldr,p(nrank+1),1,one,p,1)
         call dcopy(nrank,res0,1,res,1)
         call daxpy(nrank,-one,p,1,res,1)
      end if
c                                 end of lssetx
      end

      subroutine lsadds(unitq,vertex,inform,k1,k2,nactiv,nartif,nz,
     *                  nfree,nrank,nrejtd,nres,ngq,n,ldzy,lda,ldr,ldt,
     *                  istate,kactiv,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)
c----------------------------------------------------------------------
c     lsadds  includes general constraints k1 thru k2 as new rows of
c     the tq factorization stored in t, zy.  if nrank is nonzero, the
c     changes in q are reflected in nrank by n triangular factor r such
c     that
c                         c  =  p ( r ) q,
c                                 ( 0 )
c     where  p  is orthogonal.

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  condmx
      integer           inform, k1, k2, lda, ldr, ldt, ldzy, msglvl, n,
     *                  nactiv, nartif, nfree, ngq, nrank, nrejtd, nres,
     *                  nz
      logical           unitq, vertex

      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin

      double precision  cndmax, rnorm, rowmax, rtmax
      integer           i, iadd, iartif, ifix, iswap, jadd, k, l, nzadd
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. external subroutines ..
      external          lsadd , scond 

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg008/asize, dtmax, dtmin

      rtmax = wmach(8)
c
c     estimate the condition number of the constraints that are not
c     to be refactorized.
c
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
            call lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
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
c
      if (nactiv.lt.k2) then
c
c        some of the constraints were classed as dependent and not
c        included in the factorization.  re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.
c
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
c
c        if a vertex is required, add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.
c
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
                  call lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,
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

      nrejtd = k2 - nactiv
c                                 end of lsadds
      end

      subroutine lpcore(prbtyp,msg,cset,named,names,rset,unitq,iter,
     *                  itmax,jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,
     *                  nz,istate,kactiv,kx,obj,numinf,xnorm,a,
     *                  ax,bl,bu,cvec,featol,featlu,x,iw,w)
c----------------------------------------------------------------------
c     lpcore  is a subroutine for linear programming.
c     on entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.  the arrays  istate,
c     kactiv  and  kx  will have been set accordingly
c     and the arrays  t  and  q  will contain the tq factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  the tq factorization of the resulting (nactiv by nfree)
c     matrix is  a(free)*q = (0 t),  where q is (nfree by nfree) and t
c     is upper-triangular.
c
c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     after kdegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.
c
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
c----------------------------------------------------------------------

      integer           lenlc
      parameter         (lenlc=20)
      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      character*6       empty
      parameter         (empty='      ')

      double precision  obj, xnorm
      integer           iter, itmax, jinf, lda, n, nactiv, nclin, nfree,
     *                  nrz, numinf, nviol, nz
      logical           cset, named, rset, unitq
      character*2       prbtyp
      character*6       msg

      double precision  a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  cvec(*), featlu(n+nclin), featol(n+nclin), w(*),
     *                  x(n)
      integer           istate(n+nclin), iw(*), kactiv(n), kx(n)
      character*8       names(*)

c     .. scalars in common ..
      double precision  alfa, asize, bigbnd, bigdx, bndlow, bndupp,
     *                  dtmax, dtmin, epspt3, epspt5, epspt8, epspt9,
     *                  tolact, tolfea, tolinc, tolrnk, tolx0, trulam
      integer           idbglc, iprint, iprnt, isdel, isumm, isumry,
     *                  itmax1, itmax2, itnfix, jadd, jdel, kchk,
     *                  kcycle, kdegen, lcrash, ldbglc, ldq, ldt,
     *                  lennam, lines1, lines2, lprob, maxact, maxnz,
     *                  mm, msglc, mxfree, ncolt, ndegen, nn, nnclin,
     *                  nout, nprob
      logical           cmdbg, header, lcdbg, prnt
c     .. arrays in common ..
      double precision  rpadlc(23), rpsvlc(mxparm)
      integer           icmdbg(ldbg), ilcdbg(ldbg), ipadlc(14),
     *                  ipsvlc(mxparm), loclc(lenlc), nfix(2)

      double precision  alfap, alfhit, bigalf, biggst, condmx, condrz,
     *                  condt, dinky, dnorm, dzz, epsmch, errmax, flmax,
     *                  gfnorm, grznrm, gznorm, objsiz, rtmax, smllst,
     *                  suminf, tinyst, trubig, trusml, wssize, zerolm
      integer           iadd, idbg, ifix, inform, is, it, j, jbigst,
     *                  jmax, jsmlst, jtiny, kbigst, kdel, ksmlst, lad,
     *                  lanorm, lcq, ld, ldr, lgq, lq, lr, lrlam, lt,
     *                  ltmp, lwrk, lwtinf, msgdbg, msglvl, msgsvd,
     *                  nctotl, nfixed, ngq, nmoved, notopt, ntfixd
      logical           firstv, fp, hitlow, lp, move, onbnd, overfl,
     *                  unbndd

      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*80      rec(3)
c     .. external functions ..
      double precision  ddot, dnrm2, sdiv 
      external          ddot, dnrm2, sdiv 
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dscal, cmsinf, cmmul2,
     *                  cmmul1, cmfeas, cmdgen, cmchzr, lpcolr,
     *                  cmqmul, rzdel, rzadd, sload

      common            /ngg003/loclc
      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg004/lennam, ldt, ncolt, ldq
      common            /ngg005/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg007/alfa, trulam, isdel, jdel, jadd, header,
     *                  prnt
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg009/ilcdbg, lcdbg
      common            /ngg010/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /ngg4fb/icmdbg, cmdbg
      common            /ngg011/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)

      save              firstv
c----------------------------------------------------------------------
c     specify the machine-dependent parameters.

      epsmch = wmach(3)
      flmax = wmach(7)
      rtmax = wmach(8)

      if (cset) then
         ngq = 2
      else
         ngq = 1
      end if

      lp = prbtyp .eq. 'lp' .or. prbtyp .eq. 'lp'
      fp = .not. lp
c
      ldr = ldt
      it = 1
c
      lanorm = loclc(4)
      lad = loclc(5)
c
      ld = loclc(7)
      lgq = loclc(8)
      lcq = loclc(9)
      lrlam = loclc(10)
c
      lr = loclc(11)
      lt = loclc(12)
      lq = loclc(13)
      lwtinf = loclc(14)
      lwrk = loclc(15)

c     we need a temporary array when changing the active set.
c     use the multiplier array.

      ltmp = lrlam

      if (iter.eq.0) then
c        -------------------------
c        first entry.  initialize.
c        -------------------------
         jadd = 0
         jdel = 0
         isdel = 0
         firstv = .false.

         alfa = zero
         dzz = one
      end if

      nctotl = n + nclin
      nviol = 0

      condmx = flmax

      call cmsinf(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,a,
     *            featol,w(lgq),x,w(lwtinf))
c
      if (numinf.gt.0) then
         call cmqmul(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),w(lwrk))
      else if (lp) then
         call dcopy(n,w(lcq),1,w(lgq),1)
      end if
c
      if (numinf.eq.0 .and. lp) then
         obj = ddot(n,cvec,1,x,1)
      else
         obj = suminf
      end if
c
      msg = empty
c
c*    ======================start of main loop==========================
c     +    do while (msg .eq. empty)
   20 if (msg.eq.empty) then
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
         if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),
     *       1)
c

c        print the details of this iteration.

c        define small quantities that reflect the size of x, r and
c        the constraints in the working set.
c
         if (prnt) then
            condt = one
            if (nactiv.gt.0) condt = sdiv (dtmax,dtmin,overfl)

            header = .false.
            jdel = 0
            jadd = 0
            alfa = zero
         end if
c
         if (numinf.gt.0) then
            dinky = epspt8*abs(suminf)
         else
            objsiz = one + abs(obj)
            wssize = zero
            if (nactiv.gt.0) wssize = dtmax
            dinky = epspt8*max(wssize,objsiz,gfnorm)
         end if

c        if the reduced gradient z'g is small enough,
c        lagrange multipliers will be computed.
c
         if (numinf.eq.0 .and. fp) then
            msg = 'feasbl'
            nfixed = n - nfree
            call sload(nactiv+nfixed,zero,w(lrlam),1)
            go to 20
         end if
c
         if (grznrm.le.dinky) then
c           =========================================================
c           the point  x  is a constrained stationary point.
c           compute lagrange multipliers.
c           =========================================================
c           define what we mean by 'tiny' and non-optimal multipliers.
c
            notopt = 0
            jdel = 0
            zerolm = -dinky
            smllst = -dinky
            biggst = dinky + one
            tinyst = dinky
c
            call cmmul1(prbtyp,msglvl,n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,w(lanorm),w(lgq),w(lrlam),w(lt),
     *                  w(lwtinf))
c
            if (nrz.lt.nz) call cmmul2(msglvl,n,nrz,nz,zerolm,notopt,
     *                                 numinf,trusml,smllst,jsmlst,
     *                                 tinyst,jtiny,w(lgq))
c
            if (abs(jsmlst).gt.0) then
c              ------------------------------------------------------
c              delete a constraint.
c              ------------------------------------------------------
c              cmmul1  or  cmmul2  found a non-optimal multiplier.
c
               trulam = trusml
               jdel = jsmlst
c
               if (jsmlst.gt.0) then
c
c                 regular constraint.
c
                  kdel = ksmlst
                  isdel = istate(jdel)
                  istate(jdel) = 0
               end if
            else
               if (numinf.gt.0 .and. jbigst.gt.0) then
c
c                 no feasible point exists for the constraints but
c                 the sum of the constraint violations can be reduced
c                 by moving off constraints with multipliers greater
c                 than 1.
c
                  jdel = jbigst
                  kdel = kbigst
                  isdel = istate(jdel)
                  if (trubig.le.zero) is = -1
                  if (trubig.gt.zero) is = -2
                  istate(jdel) = is
                  trulam = trubig
                  firstv = .true.
                  numinf = numinf + 1
               end if
            end if
c
            if (jdel.eq.0) then
               if (numinf.gt.0) then
                  msg = 'infeas'
               else
                  msg = 'optiml'
               end if
               go to 20
            else 
c debug debug
               if (jdel.gt.0.and.nfree.eq.ldq) then 
c                 write (*,*) 'bugwandita!'
                  msg = 'infeas'
                  goto 20
               end if

            end if
c
c           constraint  jdel  has been deleted.
c           update the  tq  factorization.
c
            call rzdel(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,w(lt),w(lgq),w(lq),w(lwrk)
     *                  ,w(ld),w(lrlam))
            if (rset) call lpcolr(nrz,ldr,w(lr),one)
c
            prnt = .false.
         else
c           ============================================================
c           compute a search direction.
c           ============================================================
            if (iter.ge.itmax) then
               msg = 'itnlim'
               go to 20
            end if
c
            prnt = .true.
            iter = iter + 1

            call dcopy(nrz,w(lgq),1,w(ld),1)
            call dscal(nrz,(-one),w(ld),1)

            dnorm = dnrm2(nrz,w(ld),1)

            call cmqmul(1,n,nrz,nfree,ldq,unitq,kx,w(ld),w(lq),w(lwrk))
            call dgemv('no transpose',nclin,n,one,a,lda,w(ld),1,zero,
     *                 w(lad),1)

c           find the constraint we bump into along d.
c           update  x  and  ax  if the step alfa is nonzero.

c           alfhit is initialized to bigalf. if it remains that value
c           after the call to  cmchzr, it is regarded as infinite.

            bigalf = sdiv (bigdx,dnorm,overfl)

            call cmchzr(firstv,n,nclin,istate,bigalf,bigbnd,dnorm,
     *                  hitlow,move,onbnd,unbndd,alfhit,alfap,jadd,
     *                  w(lanorm),w(lad),ax,bl,bu,featol,featlu,w(ld),x)
c
            if (unbndd) then
               msg = 'unbndd'
               go to 20
            end if

            alfa = alfhit
            call daxpy(n,alfa,w(ld),1,x,1)

            if (nclin.gt.0) call daxpy(nclin,alfa,w(lad),1,ax,1)
            xnorm = dnrm2(n,x,1)

c           add a constraint to the working set.
c           update the  tq  factors of the working set.
c           use  d  as temporary work space.

            if (bl(jadd).eq.bu(jadd)) then
               istate(jadd) = 3
            else if (hitlow) then
               istate(jadd) = 1
            else
               istate(jadd) = 2
            end if

            if (jadd.gt.n) then
               iadd = jadd - n
            else
               if (alfa.ge.zero) then
                  if (hitlow) then
                     x(jadd) = bl(jadd)
                  else
                     x(jadd) = bu(jadd)
                  end if
               end if
               do 40 ifix = 1, nfree
                  if (kx(ifix).eq.jadd) go to 60
   40          continue
   60       end if

            call rzadd(unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,nrz,ngq,n,lda,ldq,ldr,ldt,kx,condmx,dzz,a,
     *                  w(lr),w(lt),w(lgq),w(lq),w(lwrk),w(lrlam),w(ld),
     *                  msglvl)

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

            call daxpy(nctotl,tolinc,featlu,1,featol,1)

            if (mod(iter,kchk).eq.0) then
c              ------------------------------------------------------
c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  set inform to force iterative
c              refinement and a switch to phase 1.
c              ------------------------------------------------------
               call cmfeas(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,
     *                     bl,bu,featol,x)

            end if
c
            if (mod(iter,kdegen).eq.0) then
c
c              every  kdegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.
c
               call cmdgen('end of cycle',msglvl,n,nclin,nmoved,iter,
     *                     numinf,istate,bigbnd,ax,bl,bu,featol,featlu,
     *                     x)
c
               nviol = nviol + nmoved
            end if
c
            if (nviol.gt.0) then
               msg = 'resetx'
               go to 20
            end if
c
            if (numinf.ne.0) then
               call cmsinf(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,
     *                     bu,a,featol,w(lgq),x,w(lwtinf))
c
               if (numinf.gt.0) then
                  call cmqmul(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),
     *                        w(lwrk))
               else if (lp) then
                  call dcopy(n,w(lcq),1,w(lgq),1)
               end if
            end if
c
            if (numinf.eq.0 .and. lp) then
               obj = ddot(n,cvec,1,x,1)
            else
               obj = suminf
            end if
         end if
         go to 20
c        +    end while
      end if
c     ======================end of main loop============================
c
      if (msg.eq.'optiml') then
         if (lp) then
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
         end if
      else if (msg.eq.'unbndd' .and. numinf.gt.0) then
         msg = 'infeas'
      end if
c                                 end of lpcore
      end

      subroutine nploc(n,nclin,ncnln,nctotl,litotl,lwtotl)
c----------------------------------------------------------------------
c     nploc   allocates the addresses of the work arrays for npcore and
c     lscore.
c----------------------------------------------------------------------
      integer           lenls
      parameter         (lenls=20)
      integer           lennp
      parameter         (lennp=35)
      integer           ldbg
      parameter         (ldbg=5)

      integer           litotl, lwtotl, n, nclin, ncnln, nctotl
c     .. scalars in common ..
      integer           ldt, ldzy, lennam, ncolt
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg), locls(lenls), locnp(lennp)

      integer           ladx, lanorm, laqp, lbl, lbu, lc1mul, lcjac,
     *                  lcjdx, lcmul, lcs1, lcs2, ldlam, ldslk, ldx,
     *                  lenaqp, lent, lenzy, lfeatl, lgq, lgq1, lgrad,
     *                  lhctrl, lhfrwd, liperm, lkactv, lkx, lneedc,
     *                  lqpadx, lqpdx, lqpgq, lqphz, lqptol, lrho,
     *                  lrlam, lrpq, lrpq0, lslk, lslk1, lt, lwrk1,
     *                  lwrk2, lwrk3, lwtinf, lx1, lzy, miniw, minw

      common            /ngg012/locls
      common            /ngg013/locnp
      common            /ngg004/lennam, ldt, ncolt, ldzy
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
c
      miniw = litotl + 1
      minw = lwtotl + 1
c
c     assign array lengths that depend upon the problem dimensions.
c
      if (nclin+ncnln.eq.0) then
         lent = 0
         lenzy = 0
      else
         lent = ldt*ncolt
         lenzy = ldzy*ldzy
      end if
c
      if (ncnln.eq.0) then
         lenaqp = 0
      else
         lenaqp = (nclin+ncnln)*n
      end if
c
      lkactv = miniw
      lkx = lkactv + n
      lneedc = lkx + n
      liperm = lneedc + ncnln
      miniw = liperm + nctotl
c
      lhfrwd = minw
      lhctrl = lhfrwd + n
      lanorm = lhctrl + n
      lqpgq = lanorm + nclin + ncnln
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
c
c     assign the addresses for the workspace arrays used by  npiqp .
c
      lqpadx = minw
      lqpdx = lqpadx + nclin + ncnln
      lrpq = lqpdx + n
      lrpq0 = lrpq + n
      lqphz = lrpq0 + n
      lwtinf = lqphz + n
      lwrk1 = lwtinf + nctotl
      lqptol = lwrk1 + nctotl
      minw = lqptol + nctotl
c
      locls(3) = lqpadx
      locls(4) = lqpdx
      locls(5) = lrpq
      locls(6) = lrpq0
      locls(7) = lqphz
      locls(13) = lwtinf
      locls(14) = lwrk1
      locls(15) = lqptol
c
c     assign the addresses for arrays used in npcore.
c
      laqp = minw
      ladx = laqp + lenaqp
      lbl = ladx + nclin + ncnln
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
      lcs2 = lcs1 + ncnln
      lc1mul = lcs2 + ncnln
      lcmul = lc1mul + ncnln
      lcjdx = lcmul + ncnln
      ldlam = lcjdx + ncnln
      ldslk = ldlam + ncnln
      lrho = ldslk + ncnln
      lwrk3 = lrho + ncnln
      lslk1 = lwrk3 + ncnln
      lslk = lslk1 + ncnln
      minw = lslk + ncnln
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
      lgrad = lcjac + ncnln*n
      minw = lgrad + n
c
      locnp(25) = lhfrwd
      locnp(26) = lhctrl
      locnp(27) = lcjac
      locnp(28) = lgrad
c
      litotl = miniw - 1
      lwtotl = minw - 1
c                                 end of nploc
      end

      subroutine npcore(named,names,unitq,inform,majits,n,nclin,ncnln,
     *                  nctotl,nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,
     *                  nfun,ngrad,istate,kactiv,kx,objf,fdnorm,xnorm,
     *                  confun,objfun,aqp,ax,bl,bu,c,cjac,cjacu,clamda,
     *                  featol,grad,gradu,r,x,iw,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     npcore  is the core routine for  e04ucf,  a sequential quadratic
c     programming (sqp) method for nonlinearly constrained optimization.

      integer           lenls
      parameter         (lenls=20)
      integer           lennp
      parameter         (lennp=35)
      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      double precision  growth
      parameter         (growth=1.0d+2)

      double precision  fdnorm, objf, xnorm
      integer           inform, ldaqp, ldcj, ldcju, ldr, lenw, majits,
     *                  n, nactiv, nclin, ncnln, nctotl, nfree, nfun,
     *                  ngrad, nz
      logical           named, unitq

      double precision  aqp(ldaqp,*), ax(*), bl(nctotl), bu(nctotl),
     *                  c(*), cjac(ldcj,*), cjacu(ldcju,*),
     *                  clamda(nctotl), featol(nctotl), grad(n),
     *                  gradu(n), r(ldr,*), user(*), w(lenw), x(n)
      integer           istate(*), iuser(*), iw(*), kactiv(n), kx(n)
      character*8       names(*)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  asize, bigbnd, bigdx, bndlow, bndupp, cdint,
     *                  ctol, drmax, drmin, dtmax, dtmin, dxlim, epspt3,
     *                  epspt5, epspt8, epspt9, epsrf, eta, fdint, ftol,
     *                  hcndbd, rcndbd, rfrobn, rhodmp, rhomax, rhonrm,
     *                  scale, tolact, tolfea, tolrnk
      integer           idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, ldt,
     *                  ldzy, lennam, lfdset, lformh, lines1, lines2,
     *                  lprob, lverfy, lvlder, lvldif, lvrfyc, msgls,
     *                  msgnp, ncdiff, ncolt, nfdiff, nlnf, nlnj, nlnx,
     *                  nload, nn, nnclin, nncnln, nout, nprob, nsave
      logical           cmdbg, incrun, npdbg
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           icmdbg(ldbg), inpdbg(ldbg), ipadls(18),
     *                  ipadnp(12), ipsvls(mxparm), ipsvnp(mxparm),
     *                  jverfy(4), locls(lenls), locnp(lennp)

      double precision  alfa, alfbnd, alfdx, alflim, alfmax, alfmin,
     *                  alfsml, cnorm, cond, condh, condhz, condt,
     *                  cvnorm, dinky, drzmax, drzmin, dxnorm, errmax,
     *                  flmax, gdx, gfnorm, glf1, glf2, glnorm, gltest,
     *                  grdalf, gtest, gznorm, obj, objalf, objsiz,
     *                  qpcurv, rootn, rtftol, rtmax, xsize
      integer           idbg, info, jmax, ladx, lanorm, laqp, lbl, lbu,
     *                  lc1mul, lcjac1, lcjdx, lcjdx1, lcmul, lcs1,
     *                  lcs2, ldcj1, ldlam, ldslk, ldx, lgq, lgq1,
     *                  lhctrl, lhfrwd, lhpq, linact, liperm, lneedc,
     *                  lqptol, lqrwrk, lrho, lrlam, lrpq, lslk, lslk1,
     *                  lt, lvioln, lwrk1, lwrk2, lwrk3, lwtinf, lx1,
     *                  lzy, majit0, minits, mjrdbg, mnr, mnrdbg,
     *                  mnrsum, mode, msgqp, ncqp, nl,
     *                  nlnact, nlserr, nmajor, nminor, nplin, nqperr,
     *                  nqpinf, nstate, numinf, nviol
      logical           centrl, convpt, convrg, done, error, feasqp,
     *                  goodgq, infeas, needfd, newgq, optiml, overfl
      character*5       mjrmsg

      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      logical           ktcond(2)
      character*80      rec(5)
c     .. external functions ..
      double precision  ddot, dnrm2, sdiv 
      external          ddot, dnrm2, sdiv 
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, cmqmul, cmprt, npupdt,
     *                  npmrt , npsrch, npiqp , npfeas, nprset,
     *                  npfd , npalf, iload , scond , smcopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg012/locls
      common            /ngg013/locnp

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg004/lennam, ldt, ncolt, ldzy
      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg015/lvrfyc, jverfy
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg016/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ngg017/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ngg018/rcndbd, rfrobn, drmax, drmin
      common            /ngg019/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ngg4fb/icmdbg, cmdbg
      common            /ngg002/inpdbg, npdbg
      common            /ngg020/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /ngg021/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c----------------------------------------------------------------------
c
c     specify machine-dependent parameters.
c
      flmax = wmach(7)
      rtmax = wmach(8)
c
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
c
      lcjac1 = laqp + nclin
      lcjdx = ladx + nclin
      lvioln = lwrk3
c
c     initialize
c
      mjrmsg = '     '
      nqpinf = 0
      mnrsum = 0
c
      majit0 = majits
      nplin = n + nclin
      ncqp = nclin + ncnln
      nl = min(nplin+1,nctotl)
c
      ldcj1 = max(ncqp,1)
c
      needfd = lvlder .eq. 0 .or. lvlder .eq. 2 .or.
     *         (lvlder.eq.1 .and. ncnln.gt.0)
c
      alfa = zero
      alfdx = zero
      rtftol = sqrt(ftol)
      rootn = sqrt(dble(n))

c     information from the feasibility phase will be used to generate a
c     hot start for the first qp subproblem.

      call dcopy(nctotl,featol,1,w(lqptol),1)
c
      nstate = 0
c
      objalf = objf
      if (ncnln.gt.0) then
         objalf = objalf - ddot(ncnln,w(lcmul),1,c,1)
      end if
c
      newgq = .false.
c
c*    ==================================================================
c+    repeat                             (until converged or error exit)
c

c     see if we want to save the details of this iteration.

c     20 if (mod(majits,ksave).eq.0 .and. majits.ne.majit0) then
c        call npsavr(unitq,n,nclin,ncnln,ldr,ldq,nfree,nsave,majits,
c    *               istate,kx,w(lhfrwd),w(lhctrl),w(lcmul),r,w(lrho),x,
c    *               x,w(lq))
c     end if
c
c   *    ===============================================================
c   +    repeat                         (until a good gradient is found)
c
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
            call npfd (centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,iw(lneedc),bl,
     *                  bu,c,w(lwrk2),w(lwrk3),cjac,cjacu,grad,gradu,
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

c     (1) solve an inequality quadratic program (iqp) for the
c         search direction and multiplier estimates.
c     (2) for each nonlinear inequality constraint,  compute
c         the slack variable for which the merit function is
c         minimized.
c     (3) compute the search direction for the slack variables
c         and multipliers.

c     note that the array violn is wrk3.

      call npiqp (feasqp,unitq,nqperr,majits,mnr,n,nclin,ncnln,ldcj,
     *            ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,numinf,istate,
     *            kactiv,kx,dxnorm,gdx,qpcurv,aqp,w(ladx),w(lanorm),ax,
     *            bl,bu,c,cjac,clamda,w(lcmul),w(lcs1),w(ldlam),w(ldslk)
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

c     compute quantities needed for the convergence test.

c     compute the norms of the projected gradient and the
c     gradient with respect to the free variables.

      gznorm = zero
      if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)
      gfnorm = gznorm
      if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),1)

c     if the forward-difference estimate of the transformed
c     gradient of the lagrangian function is small,  switch to
c     central differences, recompute the derivatives and re-solve
c     the qp.

      goodgq = .true.
      if (needfd .and. .not. centrl) then
         glnorm = dnrm2(n,w(lhpq),1)
         if (ncnln.eq.0) then
            cnorm = zero
         else
            cnorm = dnrm2(ncnln,c,1)
         end if
c
         gltest = (one+abs(objf)+abs(cnorm))*epsrf/fdnorm
         if (glnorm.le.gltest) then
            goodgq = .false.
            mjrmsg(3:3) = 'central differences'
            lvldif = 2
            newgq = .true.

         end if
c
      end if
c
c     +       until     (goodgq)
      if ( .not. goodgq) go to 40
c

c     (1) compute the number of constraints that are violated by more
c         than featol.
c     (2) compute the 2-norm of the residuals of the constraints in
c         the qp working set.

      call npfeas(n,nclin,ncnln,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *            ax,bl,bu,c,featol,x,w(lwrk2))
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
         condt = sdiv (dtmax,dtmin,overfl)
      end if
c
      call scond (n,r,ldr+1,drmax,drmin)
c
      condh = sdiv (drmax,drmin,overfl)
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
         condhz = sdiv (drzmax,drzmin,overfl)
         if (condhz.lt.rtmax) then
            condhz = condhz*condhz
         else
            condhz = flmax
         end if
      end if
c

c     test for convergence.
c     the point test convpt checks for a k-t point at the initial
c     point or after a large change in x.

      convpt = dxnorm .le. epspt8*gtest .and. nviol .eq. 0 .and.
     *         nqperr .le. 1
c
      ktcond(1) = gznorm .lt. dinky
      ktcond(2) = nviol .eq. 0
      optiml = ktcond(1) .and. ktcond(2)
c
      convrg = majits .gt. 0 .and. alfdx .le. rtftol*xsize
c
      infeas = convrg .and. .not. feasqp .or. nqpinf .gt. 7
c
      done = convpt .or. (convrg .and. optiml) .or. infeas
c
      objalf = objf
      grdalf = gdx
      glf1 = gdx
      if (ncnln.gt.0) then
         glf1 = glf1 - ddot(ncnln,w(lcjdx),1,clamda(nl),1)
c
c        compute the value and directional derivative of the
c        augmented lagrangian merit function.
c        the penalty parameters may be increased or decreased.
c
         call npmrt (feasqp,n,nclin,ncnln,objalf,grdalf,qpcurv,istate,
     *               w(lcjdx),w(lcmul),w(lcs1),w(ldlam),w(lrho),
     *               w(lvioln),w(lwrk1),w(lwrk2))
      end if

      alfa = zero
      error = majits .ge. nmajor
c
      if ( .not. (done .or. error)) then
         majits = majits + 1

c        make copies of information needed for the bfgs update.

         call dcopy(n,x,1,w(lx1),1)
         call dcopy(n,w(lgq),1,w(lgq1),1)
c
         if (ncnln.gt.0) then
            call dcopy(ncnln,w(lcjdx),1,w(lcjdx1),1)
            call dcopy(ncnln,w(lcmul),1,w(lc1mul),1)
            call dcopy(ncnln,w(lslk),1,w(lslk1),1)
         end if

c        ============================================================
c        compute the parameters for the linesearch.
c        ============================================================
c        alfmin is the smallest allowable step predicted by the qp
c        subproblem.

         alfmin = one
         if ( .not. feasqp) alfmin = zero

c        ------------------------------------------------------------
c        alfmax is the largest feasible steplength subject to a user-
c        defined limit alflim on the change in x.
c        ------------------------------------------------------------
         if (ncnln.gt.0 .and. needfd) then
            alfmax = one
         else
            alfmax = sdiv (bigdx,dxnorm,overfl)
            call npalf(info,n,nclin,ncnln,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,w(lanorm),w(ladx),ax,bl,bu,w(ldslk),
     *                  w(ldx),w(lslk),x)
            alfmax = alfa
            if (alfmax.lt.one+epspt3 .and. feasqp) alfmax = one
         end if
c
c        ------------------------------------------------------------
c        alfbnd is a tentative upper bound on the steplength.  if the
c        merit function is decreasing at alfbnd and certain
c        conditions hold,  alfbnd will be increased in multiples of
c        two (subject to not being greater than alfmax).
c        ------------------------------------------------------------
         if (ncnln.eq.0) then
            alfbnd = alfmax
         else
            alfbnd = min(one,alfmax)
         end if
c
c        ------------------------------------------------------------
c        alfsml trips the computation of central differences.  if a
c        trial steplength falls below alfsml, the linesearch is
c        terminated.
c        ------------------------------------------------------------
         alfsml = zero
         if (needfd .and. .not. centrl) then
            alfsml = sdiv (fdnorm,dxnorm,overfl)
            alfsml = min(alfsml,alfmax)
         end if
c
c        ============================================================
c        compute the steplength using safeguarded interpolation.
c        ============================================================
         alflim = sdiv ((one+xnorm)*dxlim,dxnorm,overfl)
         alfa = min(alflim,one)
c
         call npsrch(needfd,nlserr,n,ncnln,ldcj,ldcju,nfun,ngrad,
     *               iw(lneedc),confun,objfun,alfa,alfbnd,alfmax,alfsml,
     *               dxnorm,epsrf,eta,gdx,grdalf,glf1,glf2,objf,objalf,
     *               qpcurv,xnorm,c,w(lwrk1),cjac,cjacu,w(lcjdx),
     *               w(lwrk3),w(lc1mul),w(lcmul),w(lcs1),w(lcs2),w(ldx),
     *               w(ldlam),w(ldslk),grad,gradu,clamda(nl),w(lrho),
     *               w(lslk1),w(lslk),w(lx1),x,w(lwrk2),w,lenw,iuser,
     *               user)
c

c           npsrch  sets nlserr to the following values...
c
c           < 0  if the user wants to stop.
c             1  if the search is successful and alfa < alfmax.
c             2  if the search is successful and alfa = alfmax.
c             3  if a better point was found but too many functions
c                were needed (not sufficient decrease).
c
c           values of nlserr occurring with a nonzero value of alfa.
c             4  if alfmax < tolabs (too small to do a search).
c             5  if alfa  < alfsml (srchq only -- maybe want to switch
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
c
         if (alfa.gt.alflim) mjrmsg(4:4) = 'l'
c
         error = nlserr .ge. 4
         if (error) then

c           the linesearch failed to find a better point.
c           if exact gradients or central differences are being used,
c           or the kt conditions are satisfied, stop.  otherwise,
c           switch to central differences and solve the qp again.

            if (needfd .and. .not. centrl) then
               if ( .not. optiml) then
                  error = .false.
                  mjrmsg(3:3) = 'central differences'
                  lvldif = 2
                  newgq = .true.
               end if
            end if
         else
            if (needfd) then
c              ======================================================
c              compute the missing gradients.
c              ======================================================
               mode = 1
               ngrad = ngrad + 1
c
               if (ncnln.gt.0) then
                  call iload (ncnln,(1),iw(lneedc),1)
c
                  call confun(mode,ncnln,n,ldcju,iw(lneedc),x,w(lwrk1),
     *                        cjacu,nstate,iuser,user)
                  inform = mode
                  if (mode.lt.0) go to 60
c
                  call smcopy('general',ncnln,n,cjacu,ldcju,cjac,ldcj)
               end if
c
               call objfun(mode,n,x,obj,gradu,nstate,iuser,user)
               inform = mode
               if (mode.lt.0) go to 60
c
               call dcopy(n,gradu,1,grad,1)
c
               call npfd (centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                     fdint,fdnorm,objf,confun,objfun,iw(lneedc),
     *                     bl,bu,c,w(lwrk2),w(lwrk3),cjac,cjacu,grad,
     *                     gradu,w(lhfrwd),w(lhctrl),x,w,lenw,iuser,
     *                     user)
c
               inform = mode
               if (mode.lt.0) go to 60
c
               gdx = ddot(n,grad,1,w(ldx),1)
               glf2 = gdx
               if (ncnln.gt.0) then
                  call dgemv('n',ncnln,n,one,cjac,ldcj,w(ldx),1,zero,
     *                       w(lcjdx),1)
                  glf2 = glf2 - ddot(ncnln,w(lcjdx),1,clamda(nl),1)
               end if
            end if
c
            call dcopy(n,grad,1,w(lgq),1)
            call cmqmul(6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),
     *                  w(lwrk1))
c
            xnorm = dnrm2(n,x,1)
c
            if (ncnln.gt.0 .and. alfa.ge.one) call dcopy(ncnln,
     *          clamda(nl),1,w(lcmul),1)
c
            if (nclin.gt.0) call daxpy(nclin,alfa,w(ladx),1,ax,1)
            alfdx = alfa*dxnorm
c
c           =========================================================
c           update the factors of the approximate hessian of the
c           lagrangian function.
c           =========================================================
            call npupdt(mjrmsg,unitq,n,ncnln,nfree,nz,ldcj1,ldcj,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,w(lcjac1),cjac,
     *                  w(lcjdx1),w(lcjdx),w(lcs1),w(lcs2),w(lgq1),
     *                  w(lgq),w(lhpq),w(lrpq),clamda(nl),r,w(lwrk3),
     *                  w(lzy),w(lwrk2),w(lwrk1))
c
            call scond (n,r,ldr+1,drmax,drmin)
            cond = sdiv (drmax,drmin,overfl)
c
            if (cond.gt.rcndbd .or. rfrobn.gt.rootn*growth*drmax) then
c              ------------------------------------------------------
c              reset the condition estimator and range-space
c              partition of q'hq.
c              ------------------------------------------------------
               mjrmsg(5:5) = 'refactorize hessian'
c
               call nprset(unitq,n,nfree,nz,ldzy,ldr,iw(liperm),kx,
     *                     w(lgq),r,w(lzy),w(lwrk1),w(lqrwrk))
            end if
         end if
      end if
c
c     +    until     (done  .or.  error)
      if ( .not. (done .or. error)) go to 20
c
c     ======================end of main loop============================
c
      if (done) then
         if (convrg .and. optiml) then
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

c     set  clamda

   60 call cmprt(msgnp,nfree,ldaqp,n,nclin,nctotl,bigbnd,named,names,
     *            nactiv,istate,kactiv,kx,aqp,bl,bu,c,clamda,w(lrlam),x)
      if (ncnln.gt.0) call dcopy(ncnln,w(lcmul),1,clamda(n+nclin+1),1)

c                                 end of npcore
      end

      subroutine cmsetx(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,
     *                  ldt,istate,kactiv,kx,jmax,errmax,xnorm,a,ax,bl,
     *                  bu,featol,t,x,q,p,work)
c----------------------------------------------------------------------
c     cmsetx  computes the point on a working set that is closest in the
c     least-squares sense to the input vector x.
c
c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable rowerr
c     is set.
c----------------------------------------------------------------------


      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      integer           ntry
      parameter         (ntry=5)

      double precision  errmax, xnorm
      integer           jmax, lda, ldq, ldt, n, nactiv, nclin, nfree, nz
      logical           rowerr, unitq

      double precision  a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), p(n), q(ldq,*), t(ldt,*),
     *                  work(n), x(n)
      integer           istate(n+nclin), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)

      double precision  bnd
      integer           i, is, j, k, ktry
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           idamax
      external          ddot, dnrm2, idamax
c     .. external subroutines ..
      external          daxpy, dcopy, dgemv, dtrsv, cmqmul,
     *                  sload
c     .. intrinsic functions ..
      intrinsic         abs

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
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
c     ntry  attempts are made to get acceptable row errors.

      ktry = 1
      jmax = 1
      errmax = zero
c
c     repeat
   40 if (nactiv.gt.0) then
c
c        set work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.
c
         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(nactiv-i+1) = bnd - ddot(n,a(k,1),lda,x,1)
   60    continue
c
         call dtrsv('u','n','n',nactiv,t(1,nz+1),ldt,work,1)
         call sload(n,zero,p,1)
         call dcopy(nactiv,work,1,p(nz+1),1)
c
         call cmqmul(2,n,nz,nfree,ldq,unitq,kx,p,q,work)
         call daxpy(n,one,p,1,x,1)
      end if
c

c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.

      xnorm = dnrm2(n,x,1)
      if (nclin.gt.0) call dgemv('n',nclin,n,one,a,lda,x,1,zero,ax,1)
c

c     check the row residuals.

      if (nactiv.gt.0) then
         do 80 k = 1, nactiv
            i = kactiv(k)
            j = n + i
            is = istate(j)
            if (is.eq.1) work(k) = bl(j) - ax(i)
            if (is.ge.2) work(k) = bu(j) - ax(i)
   80    continue
c
         jmax = idamax(nactiv,work,1)
         errmax = abs(work(jmax))
      end if
c
      ktry = ktry + 1
c     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      if ( .not. (errmax.le.featol(jmax) .or. ktry.gt.ntry)) go to 40
c
      rowerr = errmax .gt. featol(jmax)
c                                 end of cmsetx
      end

      subroutine lploc(cset,n,nclin,litotl,lwtotl)
c----------------------------------------------------------------------
c     lploc   allocates the addresses of the work arrays for lpcore.

c     note that the arrays ( gq, cq ) lie in contiguous areas of
c     workspace.
c----------------------------------------------------------------------
      integer           lenlc
      parameter         (lenlc=20)
      integer           ldbg
      parameter         (ldbg=5)

      integer           litotl, lwtotl, n, nclin
      logical           cset
c     .. scalars in common ..
      integer           iprint, isumm, ldq, ldt, lennam, lines1, lines2,
     *                  ncolt, nout
      logical           lcdbg
c     .. arrays in common ..
      integer           ilcdbg(ldbg), loclc(lenlc)

      integer           lad, lanorm, lcq, ld, lencq, lenq, lenrt,
     *                  lfeatu, lgq, lkactv, lkx, lq, lr, lrlam, lt,
     *                  lwrk, lwtinf, miniw, minw

      common            /ngg003/loclc
      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg004/lennam, ldt, ncolt, ldq
      common            /ngg009/ilcdbg, lcdbg
c----------------------------------------------------------------------
c     refer to the first free space in the work arrays.

      miniw = litotl + 1
      minw = lwtotl + 1

c     integer workspace.

      lkactv = miniw
      lkx = lkactv + n
      miniw = lkx + n
c

c     real workspace.
c     assign array lengths that depend upon the problem dimensions.

      lenrt = ldt*ncolt
      if (nclin.eq.0) then
         lenq = 0
      else
         lenq = ldq*ldq
      end if
c
      if (cset) then
         lencq = n
      else
         lencq = 0
      end if
c

c     we start with arrays that can be preloaded by smart users.

      lfeatu = minw
      minw = lfeatu + nclin + n
c
c     next comes stuff used by  lpcore  and  e04nfz.
c
      lanorm = minw
      lad = lanorm + nclin
      ld = lad + nclin
      lgq = ld + n
      lcq = lgq + n
      lrlam = lcq + lencq
      lr = lrlam + n
      lt = lr
      lq = lt + lenrt
      lwtinf = lq + lenq
      lwrk = lwtinf + n + nclin
      minw = lwrk + n + nclin
c
c     load the addresses in loclc.
c
      loclc(1) = lkactv
      loclc(2) = lkx
c
      loclc(3) = lfeatu
      loclc(4) = lanorm
      loclc(5) = lad
c
      loclc(7) = ld
      loclc(8) = lgq
      loclc(9) = lcq
      loclc(10) = lrlam
      loclc(11) = lr
      loclc(12) = lt
      loclc(13) = lq
      loclc(14) = lwtinf
      loclc(15) = lwrk

      litotl = miniw - 1
      lwtotl = minw - 1
c                                 end of lploc
      end

      subroutine rzadds(unitq,vertex,k1,k2,it,nactiv,nartif,nz,nfree,
     *                  nrejtd,ngq,n,ldq,lda,ldt,istate,kactiv,kx,
     *                  condmx,a,t,gqm,q,w,c,s,msglvl)
c----------------------------------------------------------------------
c     rzadds  includes general constraints  k1  thru  k2  as new rows of
c     the  tq  factorization:
c              a(free) * q(free)  = (  0 t )
c                        q(free)  = (  z y )
c
c     a) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt)  of the array  t.
c----------------------------------------------------------------------
      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  condmx
      integer           it, k1, k2, lda, ldq, ldt, msglvl, n, nactiv,
     *                  nartif, nfree, ngq, nrejtd, nz
      logical           unitq, vertex

      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n),
     *                  t(ldt,*), w(n)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin, epspt3, epspt5, epspt8,
     *                  epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)

      double precision  cndmax, cond, delta, drzz, dtnew, rnorm, rowmax,
     *                  rtmax, tdtmax, tdtmin
      integer           i, iadd, iartif, ifix, inform, iswap, j, jadd,
     *                  jt, k, l, nzadd
      logical           overfl, rset

      character*80      rec(4)
c     .. external functions ..
      double precision  dnrm2, sdiv 
      external          dnrm2, sdiv 
c     .. external subroutines ..
      external          dcopy, dgemv, dger, cmqmul, rzadd, scond ,
     *                  sgrfg, smload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      rtmax = wmach(8)

      jt = nz + 1

c     estimate the condition number of the constraints already
c     factorized.

      if (nactiv.eq.0) then
         dtmax = zero
         dtmin = one
         if (unitq) then

c           first general constraint added.  set  q = i.

            call smload('general',nfree,nfree,zero,one,q,ldq)
            unitq = .false.
         end if
      else
         call scond (nactiv,t(it,jt),ldt+1,dtmax,dtmin)
      end if

      do 20 k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv.lt.nfree) then

            overfl = .false.

c           transform the incoming row of  a  by  q'.

            call dcopy(n,a(iadd,1),lda,w,1)
            call cmqmul(8,n,nz,nfree,ldq,unitq,kx,w,q,s)

c           check that the incoming row is not dependent upon those
c           already in the working set.

            dtnew = dnrm2(nz,w,1)
            if (nactiv.eq.0) then

c              this is the first general constraint in the working set.

               cond = sdiv (asize,dtnew,overfl)
               tdtmax = dtnew
               tdtmin = dtnew
            else
c
c              there are already some general constraints in the working
c              set. update the estimate of the condition number.
c
               tdtmax = max(dtnew,dtmax)
               tdtmin = min(dtnew,dtmin)
               cond = sdiv (tdtmax,tdtmin,overfl)
            end if
c
            if (cond.ge.condmx .or. overfl) then

c              this constraint appears to be dependent on those already
c              in the working set.  skip it.


c
               istate(jadd) = 0
               kactiv(k) = -kactiv(k)
            else
               if (nz.gt.1) then
c                 ------------------------------------------------------
c                 use a single column transformation to reduce the first
c                 nz-1  elements of  w  to zero.
c                 ------------------------------------------------------
c                 apply the householder reflection  i  -  w w'.
c                 the reflection is applied to  z  and gqm so that
c                    y  =    z  * w,   z    =  z    -  y w'  and
c                    y  =  gqm' * w,   gqm  =  gqm  -  w y',
c                 where  w = wrk1 (from householder),
c                 and    y = wrk2 (workspace).
c
c                 note that delta  has to be stored after the reflection
c                 is used.
c
                  delta = w(nz)
                  call sgrfg(nz-1,delta,w,1,zero,w(nz))
                  if (w(nz).gt.zero) then
c
                     call dgemv('n',nfree,nz,one,q,ldq,w,1,zero,s,1)
                     call dger(nfree,nz,(-one),s,1,w,1,q,ldq)
c
                     if (ngq.gt.0) then
                        call dgemv('t',nz,ngq,one,gqm,n,w,1,zero,s,1)
                        call dger(nz,ngq,(-one),w,1,s,1,gqm,n)
                     end if
                  end if
c
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
   20 continue
c
      if (nactiv.lt.k2) then
c
c        some of the constraints were classed as dependent and not
c        included in the factorization.  re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.
c
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
c
c        if a vertex is required,  add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.
c
         if (vertex) then
            rset = .false.
            cndmax = rtmax
            drzz = one
            nzadd = nz
            do 80 iartif = 1, nzadd
               if (unitq) then
                  ifix = nfree
                  jadd = kx(ifix)
               else
                  rowmax = zero
                  do 60 i = 1, nfree
                     rnorm = dnrm2(nz,q(i,1),ldq)
                     if (rowmax.lt.rnorm) then
                        rowmax = rnorm
                        ifix = i
                     end if
   60             continue
                  jadd = kx(ifix)
c
                  call rzadd(unitq,rset,inform,ifix,iadd,jadd,it,
     *                        nactiv,nz,nfree,nz,ngq,n,lda,ldq,ldt,ldt,
     *                        kx,cndmax,drzz,a,t,t,gqm,q,w,c,s,msglvl)
               end if
               nfree = nfree - 1
               nz = nz - 1
               nartif = nartif + 1
               istate(jadd) = 4
   80       continue
         end if
c
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

c                                 end of rzadds

      end

      subroutine npdflt(n,nclin,ncnln,title)
c----------------------------------------------------------------------
c     npdflt  loads the default values of parameters not set in the
c     options file.

      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      double precision  point3, point8
      parameter         (point3=3.3d-1,point8=0.8d+0)
      double precision  point9, two
      parameter         (point9=0.9d+0,two=2.0d+0)
      double precision  tenp6, hundrd
      parameter         (tenp6=1.0d+6,hundrd=10.0d+1)
      double precision  rdummy
      integer           idummy
      parameter         (rdummy=-11111.d+0,idummy=-11111)
      double precision  gigant
      parameter         (gigant=1.0d+20*0.99999d+0)
      double precision  wrktol
      parameter         (wrktol=1.0d-2)

      integer           n, nclin, ncnln
      character*(*)     title
c     .. scalars in common ..
      double precision  bigbnd, bigdx, bndlow, bndupp, cdint, ctol,
     *                  dxlim, epspt3, epspt5, epspt8, epspt9, epsrf,
     *                  eta, fdint, ftol, hcndbd, tolact, tolfea, tolrnk
      integer           idbgls, idbgnp, iprint, iprnt, isumm, isumry,
     *                  itmax1, itmax2, itmxnp, jvrfy1, jvrfy2, jvrfy3,
     *                  jvrfy4, ksave, lcrash, ldbgls, ldbgnp, lfdset,
     *                  lformh, lines1, lines2, lprob, lverfy, lvlder,
     *                  lvldif, lvrfyc, msgls, msgnp, ncdiff, nfdiff,
     *                  nlnf, nlnj, nlnx, nload, nn, nnclin, nncnln,
     *                  nout, nprob, nsave
      logical           cmdbg, lsdbg, npdbg
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           icmdbg(ldbg), ilsdbg(ldbg), inpdbg(ldbg),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), jverfy(4)

      double precision  condbd, dctol, epsmch
      integer           i, idbg, j, k, lent, mjrdbg, mnrdbg, msg1, msg2,
     *                  msgqp, nctotl, nmajor, nminor, nplin
      character*16      key

      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*3       chess(0:1)
      character*4       icrsh(0:2)
      character*80      rec(4)
c     .. external subroutines ..
      external          dcopy, icopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg015/lvrfyc, jverfy
      common            /ngg016/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ngg019/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls

      common            /ngg4fb/icmdbg, cmdbg
      common            /ngg002/inpdbg, npdbg
      common            /ngg020/ipsvnp, idbgnp, itmxnp, jvrfy1, jvrfy2,
     *                  jvrfy3, jvrfy4, ldbgnp, lformh, lvlder, lverfy,
     *                  msgnp, nlnf, nlnj, nlnx, nncnln, nsave, nload,
     *                  ksave, ipadnp
      common            /ngg021/rpsvnp, cdint, ctol, dxlim, epsrf, eta,
     *                  fdint, ftol, hcndbd, rpadnp
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)

      data              icrsh(0), icrsh(1), icrsh(2)/'cold', 'warm',
     *                  'hot '/
      data              chess(0), chess(1)/' no', 'yes'/
c----------------------------------------------------------------------
      epsmch = wmach(3)
      nout = 6

      condbd = max(one/(hundrd*epsmch*dble(n)),tenp6)

      nplin = n + nclin
      nctotl = nplin + ncnln

      msgnp = 0
      msgqp = 0
      iprnt = nout
      isumry = -1
      iprint = iprnt
      isumm = isumry
      lcrash = 0
      lvlder = 3
      lformh = 0

      nmajor = max(50,3*nplin+10*ncnln)
      nminor = max(50,3*nctotl)
      mjrdbg = 0
      mnrdbg = 0
      idbg = nmajor + 1
      nlnf = n
      nlnj = n
      nlnx = n

      jvrfy2 = n
      jvrfy1 = 1
      jvrfy4 = n
      jvrfy3 = 1
      lverfy = 0

      nsave = 0
      ksave = nmajor + 1
      nload = 0

      tolact = wrktol
      tolfea = epspt5
      epsrf = epspt9
      lfdset = 0
      bigbnd = gigant
      bigdx = max(gigant,bigbnd)
      dxlim = two
      eta = point9
      ftol = epsrf**point8

      hcndbd = condbd
      dctol = epspt5
      ctol = dctol

      itmax1 = max(50,3*(n+nclin+ncnln))
      jverfy(1) = jvrfy1
      jverfy(2) = jvrfy2
      jverfy(3) = jvrfy3
      jverfy(4) = jvrfy4

      npdbg = .false.
      cmdbg = .false.
      lsdbg = .false.
c                                 end of npdflt
      end

      subroutine lscrsh(cold,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *                  lda,istate,kactiv,bigbnd,tolact,a,ax,bl,bu,x,wx,
     *                  work)
c----------------------------------------------------------------------
c     lscrsh  computes the quantities  istate (optionally), kactiv,
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

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  bigbnd, tolact
      integer           lda, n, nactiv, nartif, nclin, nctotl, nfree
      logical           cold, vertex

      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  work(n), wx(n), x(n)
      integer           istate(nctotl), kactiv(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)

      double precision  b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, toobig
      integer           i, imin, is, j, jmin, k, nfixed

      character*80      rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. external subroutines ..
      external          dcopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg4nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
      flmax = wmach(7)
      biglow = -bigbnd
      bigupp = bigbnd

c     move the variables inside their bounds.


      do 20 j = 1, n
         b1 = bl(j)
         b2 = bu(j)

         if (b1.gt.biglow) then
            if (x(j).lt.b1) x(j) = b1
         end if

         if (b2.lt.bigupp) then
            if (x(j).gt.b2) x(j) = b2
         end if
   20 continue

      call dcopy(n,x,1,wx,1)

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
            if (istate(j).gt.3 .or. istate(j).lt.0) istate(j) = 0
            if (bl(j).ne.bu(j) .and. istate(j).eq.3) istate(j) = 0
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
c        +       while (j .ge. 1  .and.  nfixed + nactiv .lt. n) do
  120    if (j.ge.1 .and. nfixed+nactiv.lt.n) then
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
         if (nclin.gt.0 .and. nactiv+nfixed.lt.n) then
            do 140 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot(n,a(i,1),lda,wx,1)
  140       continue
c
            is = 1
            toobig = tolact + tolact
c
c           + while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  160       if (is.gt.0 .and. nfixed+nactiv.lt.n) then
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
     *                                        )
                     if (b2.lt.bigupp) resu = abs(ax(i)-b2)/(one+abs(b2)
     *                                        )
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
      if (vertex .and. nactiv+nfixed.lt.n) then

c        find an initial vertex by temporarily fixing some variables.

c        compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).
c
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

c                                 end of lscrsh

      end

      subroutine cmqmul(mode,n,nz,nfree,nq,unitq,kx,v,zy,wrk)
c----------------------------------------------------------------------
c     cmqmul  transforms the vector  v  in various ways using the
c     matrix  q = ( z  y )  defined by the input parameters.
c
c        mode               result
c        ----               ------
c
c          1                v = z v
c          2                v = y v
c          3                v = q v
c
c     on input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
c     on output, v  is a full n-vector.
c
c
c          4                v = z'v
c          5                v = y'v
c          6                v = q'v
c
c     on input,  v  is a full n-vector.
c     on output, v  is ordered as  ( v(free)  v(fixed) ).
c
c          7                v = y'v
c          8                v = q'v
c
c     on input,  v  is a full n-vector.
c     on output, v  is as in modes 5 and 6 except that v(fixed) is not
c     set.
c
c     modes  1, 4, 7 and 8  do not involve  v(fixed).

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      integer           mode, n, nfree, nq, nz
      logical           unitq

      double precision  v(n), wrk(n), zy(nq,*)
      integer           kx(n)

      integer           j, j1, j2, k, l, lenv, nfixed
c     .. external subroutines ..
      external          sload, dcopy, dgemv
c----------------------------------------------------------------------
c
      nfixed = n - nfree
      j1 = 1
      j2 = nfree
      if (mode.eq.1 .or. mode.eq.4) j2 = nz
      if (mode.eq.2 .or. mode.eq.5 .or. mode.eq.7) j1 = nz + 1
      lenv = j2 - j1 + 1
      if (mode.le.3) then

c        mode = 1, 2  or  3.

c
         if (nfree.gt.0) call sload(nfree,zero,wrk,1)
c
c        copy  v(fixed)  into the end of  wrk.
c
         if (mode.ge.2 .and. nfixed.gt.0) call dcopy(nfixed,v(nfree+1),
     *       1,wrk(nfree+1),1)
c
c        set  wrk  =  relevant part of  zy * v.
c
         if (lenv.gt.0) then
            if (unitq) then
               call dcopy(lenv,v(j1),1,wrk(j1),1)
            else
               call dgemv('n',nfree,j2-j1+1,one,zy(1,j1),nq,v(j1),1,one,
     *                    wrk,1)
            end if
         end if
c
c        expand  wrk  into  v  as a full n-vector.
c
         call sload(n,zero,v,1)
         do 20 k = 1, nfree
            j = kx(k)
            v(j) = wrk(k)
   20    continue
c
c        copy  wrk(fixed)  into the appropriate parts of  v.
c
         if (mode.gt.1) then
            do 40 l = 1, nfixed
               j = kx(nfree+l)
               v(j) = wrk(nfree+l)
   40       continue
         end if

      else

c        mode = 4, 5, 6, 7  or  8.

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
                  call dgemv('t',nfree,j2-j1+1,one,zy(1,j1),nq,wrk,1,
     *                       zero,v(j1),1)
               end if
            end if
         end if

c        copy the fixed components of  wrk  into the end of  v.

         if (nfixed.gt.0 .and. (mode.eq.5 .or. mode.eq.6))
     *       call dcopy(nfixed,wrk(nfree+1),1,v(nfree+1),1)
      end if

c                                 end of cmqmul

      end

      subroutine lscore(prbtyp,named,names,linobj,unitq,inform,iter,
     *                  jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nrz,n,
     *                  lda,ldr,istate,kactiv,kx,ctx,ssq,ssq1,suminf,
     *                  numinf,xnorm,bl,bu,a,clamda,ax,featol,r,x,w)
c----------------------------------------------------------------------
c     lscore  is a subroutine for linearly constrained linear-least
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

c     this material may be reproduced by or for the u.s. government
c     pursuant to the copyright license under dar clause 7-104.9(a)
c     (1979 mar).
c
c     this material is based upon work partially supported by the
c     national science foundation under grants mcs-7926009 and
c     ecs-8012974; the department of energy contract am03-76sf00326, pa
c     no. de-at03-76er72018; and the army research office contract daa29
c     -79-c-0110.

      integer           lenls
      parameter         (lenls=20)
      integer           ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      integer           mstall, mrefn
      parameter         (mstall=50,mrefn=1)

      double precision  ctx, ssq, ssq1, suminf, xnorm
      integer           inform, iter, jinf, lda, ldr, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nrz, numinf, nz
      logical           linobj, named, unitq
      character*2       prbtyp

      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  clamda(nctotl), featol(nctotl), r(ldr,*), w(*),
     *                  x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)
c     .. scalars in common ..
      double precision  asize, bigbnd, bigdx, bndlow, bndupp, dtmax,
     *                  dtmin, epspt3, epspt5, epspt8, epspt9, tolact,
     *                  tolfea, tolrnk
      integer           idbgls, iprint, iprnt, isumm, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, ldt, ldzy, lennam,
     *                  lines1, lines2, lprob, msgls, ncolt, nn, nnclin,
     *                  nout, nprob
      logical           cmdbg, lsdbg
c     .. arrays in common ..
      double precision  rpadls(23), rpsvls(mxparm)
      integer           icmdbg(ldbg), ilsdbg(ldbg), ipadls(18),
     *                  ipsvls(mxparm), locls(lenls)

      double precision  absrzz, alfa, alfhit, atphit, bigalf, cnorm,
     *                  condmx, condrz, condt, ctp, dinky, drzmax,
     *                  drzmin, err1, err2, flmax, gfnorm, grznrm,
     *                  gznorm, objsiz, palfa, pnorm, resnrm, rownrm,
     *                  trulam, wssize
      integer           iadd, idbg, ifix, irefn, is, isdel, itmax, jadd,
     *                  jbigst, jdel, jmax1, jsmlst, jtiny, kbigst,
     *                  kdel, ksmlst, lanorm, lap, lcq, lgq, lhz, lpx,
     *                  lres, lres0, lrlam, lt, lwrk, lwtinf, lzy,
     *                  msgdbg, msglvl, msgsvd, ngq, nphase, nres,
     *                  nstall, nviol
      logical           convrg, cyclin, error, firstv, hitcon, hitlow,
     *                  needfg, overfl, prnt, rowerr, singlr, stall,
     *                  statpt, unbndd, uncon, unitgz, weak

      double precision  rprmls(mxparm)
      integer           iprmls(mxparm)
      character*80      rec(3)
c     .. external functions ..
      double precision  dnrm2, sdiv 
      external          dnrm2, sdiv 
c     .. external subroutines ..
      external          cmprt, lssetx, lsmuls, lsmove, lsgset,
     *                  lsgetp, lsfeas, lsdel, lsadd , cmalf, sload,
     *                  scond

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg012/locls

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg004/lennam, ldt, ncolt, ldzy
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg4nc/ilsdbg, lsdbg
      common            /ngg008/asize, dtmax, dtmin
      common            /ngg016/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ngg019/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /ngg4fb/icmdbg, cmdbg
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (msgls,msglvl), (idbgls,idbg), (ldbgls,msgdbg)
c----------------------------------------------------------------------
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

c     set up the adresses of the contiguous arrays  ( res0, res )
c     and  ( gq, cq ).

      nres = 0
      if (nrank.gt.0) nres = 2
      ngq = 1
      if (linobj) ngq = 2

c     initialize.

      irefn = 0
      iter = 0
      jadd = 0
      jdel = 0
      nphase = 1
      nstall = 0
      numinf = -1
      nrz = 0

      if (prbtyp.eq.'fp') then
         itmax = itmax2
      else
         itmax = itmax1
      end if

      alfa = zero
      condmx = flmax
      drzmax = one
      drzmin = one
      ssq = zero

      cyclin = .false.
      error = .false.
      firstv = .false.
      prnt = .true.
      needfg = .true.
      stall = .true.
      uncon = .false.
      unitgz = .true.
      unbndd = .false.
      msglvl = 0
c     =================== start of the main loop =======================
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
            call lsgset(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,w(lres),featol,
     *                  w(lgq),w(lcq),r,x,w(lwtinf),w(lzy),w(lwrk))
            if (numinf.eq.0 .and. prbtyp.ne.'fp' .and. nphase.eq.1) then
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
      if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),1)

c        define small quantities that reflect the magnitude of  x,
c        r  and  the matrix of constraints in the working set.
c        use the largest and smallest diagonals of  r  to estimate
c        the condition number of  rz1.

      if (nrz.eq.0) then
         singlr = .false.
      else
         if (numinf.gt.0 .or. nrz.gt.nrank) then
            absrzz = zero
            singlr = .true.
         else
            call scond (nrz,r,ldr+1,drzmax,drzmin)
            absrzz = abs(r(nrz,nrz))
            rownrm = dnrm2(n,r(1,1),ldr)
            singlr = absrzz .le. drzmax*tolrnk .or. rownrm .le.
     *               tolrnk .or. abs(r(1,1)) .le. rownrm*tolrnk
         end if

      end if

      condrz = sdiv (drzmax,drzmin,overfl)

      condt = one
      if (nactiv.gt.0) condt = sdiv (dtmax,dtmin,overfl)

      if (prnt) then

         jdel = 0
         jadd = 0
         alfa = zero
      end if

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

c     if the projected gradient  z'g  is small and rz is of full
c     rank, x is a minimum on the working set.  an additional
c     refinement step is allowed to take care of an inaccurate
c     value of dinky.

      statpt = .not. singlr .and. grznrm .le. dinky .or. irefn .gt.
     *         mrefn

      if ( .not. statpt) then
c        ---------------------------------------------------------
c        compute a search direction.
c        ---------------------------------------------------------
         prnt = .true.

         error = iter .ge. itmax
         if ( .not. error) then

            irefn = irefn + 1
            iter = iter + 1

            call lsgetp(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,w(lap),
     *                  w(lres),w(lhz),w(lpx),w(lgq),w(lcq),r,w(lzy),
     *                  w(lwrk))

c           ------------------------------------------------------
c           find the constraint we bump into along p.
c           update x and ax if the step alfa is nonzero.
c           ------------------------------------------------------
c           alfhit is initialized to bigalf.  if it remains
c           that way after the call to cmalf, it will be
c           regarded as infinite.

            bigalf = sdiv (bigdx,pnorm,overfl)

            call cmalf(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfhit,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  w(lanorm),w(lap),ax,bl,bu,featol,w(lpx),x)

c           if  rz1  is nonsingular,  alfa = 1.0  will be the
c           step to the least-squares minimizer on the
c           current subspace. if the unit step does not violate
c           the nearest constraint by more than featol,  the
c           constraint is not added to the working set.

            hitcon = singlr .or. palfa .le. one
            uncon = .not. hitcon
c
            if (hitcon) then
               alfa = alfhit
            else
               jadd = 0
               alfa = one
            end if

c           check for an unbounded solution or negligible step.

            unbndd = alfa .ge. bigalf
            stall = abs(alfa*pnorm) .le. epspt9*xnorm
            if (stall) then
               nstall = nstall + 1
               cyclin = nstall .gt. mstall
            else
               nstall = 0
            end if

            error = unbndd .or. cyclin
            if ( .not. error) then
c              ---------------------------------------------------
c              set x = x + alfa*p.  update ax, gq, res and ctx.
c              ---------------------------------------------------
               if (alfa.ne.zero) call lsmove(hitcon,hitlow,linobj,
     *                                unitgz,nclin,nrank,nrz,n,ldr,jadd,
     *                                numinf,alfa,ctp,ctx,xnorm,w(lap),
     *                                ax,bl,bu,w(lgq),w(lhz),w(lpx),
     *                                w(lres),r,x,w(lwrk))

               if (hitcon) then
c                 ------------------------------------------------
c                 add a constraint to the working set.
c                 update the tq factors of the working set.
c                 use p as temporary work space.
c                 ------------------------------------------------
c                 update  istate.

                  if (bl(jadd).eq.bu(jadd)) then
                     istate(jadd) = 3
                  else if (hitlow) then
                     istate(jadd) = 1
                  else
                     istate(jadd) = 2
                  end if
                  iadd = jadd - n
                  if (jadd.le.n) then

                     do 40 ifix = 1, nfree
                        if (kx(ifix).eq.jadd) go to 60
   40                continue
                  end if
   60             continue

                  call lsadd (unitq,inform,ifix,iadd,jadd,nactiv,nz,
     *                        nfree,nrank,nres,ngq,n,lda,ldzy,ldr,ldt,
     *                        kx,condmx,a,r,w(lt),w(lres),w(lgq),w(lzy),
     *                        w(lwrk),w(lrlam),w(lpx),msglvl)

                  nrz = nrz - 1
                  nz = nz - 1

                  if (jadd.le.n) then

c                    a simple bound has been added.

                     nfree = nfree - 1
                  else

c                    a general constraint has been added.

                     nactiv = nactiv + 1
                     kactiv(nactiv) = iadd
                  end if
c
                  irefn = 0
               end if
c              ---------------------------------------------------
c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  refine the current x and set inform so
c              that feasibility is checked in lsgset.
c              ---------------------------------------------------
               call lsfeas(n,nclin,istate,bigbnd,cnorm,err1,jmax1,nviol,
     *                     ax,bl,bu,featol,x,w(lwrk))
c
               if (err1.gt.featol(jmax1)) then
                  call lssetx(linobj,rowerr,unitq,nclin,nactiv,nfree,
     *                        nrank,nz,n,nctotl,ldzy,lda,ldr,ldt,istate,
     *                        kactiv,kx,jmax1,err2,ctx,xnorm,a,ax,bl,bu,
     *                        w(lcq),w(lres),w(lres0),featol,r,w(lt),x,
     *                        w(lzy),w(lpx),w(lwrk))

                  if (rowerr) then

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

c        until      statpt  .or.  error
      if ( .not. (statpt .or. error)) go to 20

c     try and find the index jdel of a constraint to drop from
c     the working set.

      jdel = 0
c
      if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         if (n.gt.nz) call sload(n-nz,(zero),w(lrlam),1)
         jtiny = 0
         jsmlst = 0
         jbigst = 0
      else

         call lsmuls(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,nrz,
     *               istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,jtiny,
     *               jbigst,kbigst,trulam,a,w(lanorm),w(lgq),w(lrlam),
     *               w(lt),w(lwtinf))

      end if

      if ( .not. error) then
         if (jsmlst.gt.0) then

c           lsmuls found a regular constraint with multiplier less
c           than (-dinky).

            jdel = jsmlst
            kdel = ksmlst
            isdel = istate(jdel)
            istate(jdel) = 0

         else if (jsmlst.lt.0) then

            jdel = jsmlst

         else if (numinf.gt.0 .and. jbigst.gt.0) then

c           no feasible point exists for the constraints but the
c           sum of the constraint violations may be reduced by
c           moving off constraints with multipliers greater than 1.

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
         if (jdel.ne.0 .and. singlr) then
c
c           cannot delete a constraint when rz is singular.
c           probably a weak minimum.
c
            jdel = 0
         else if (jdel.ne.0) then

c           constraint jdel has been deleted.
c           update the matrix factorizations.

            call lsdel(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,w(lres),r,
     *                  w(lt),w(lgq),w(lzy),w(lwrk),w(lpx))
         end if
      end if

      irefn = 0
      convrg = jdel .eq. 0

      prnt = .false.
      uncon = .false.
      needfg = .false.

c     until       convrg  .or.  error
      if ( .not. (convrg .or. error)) go to 20

c     .....................end of main loop........................
      weak = jtiny .gt. 0 .or. singlr
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
         else if (prbtyp.ne.'fp' .and. weak) then
            inform = 1
         end if
      end if
c                                 set   clamda

      call cmprt(msglvl,nfree,lda,n,nclin,nctotl,bigbnd,named,names,
     *           nactiv,istate,kactiv,kx,a,bl,bu,x,clamda,w(lrlam),x)

c                                 end of lscore
      end

      subroutine npchkd(inform,msgnp,nstate,lvlder,nfun,ngrad,ldcj,
     *                  ldcju,n,ncnln,confun,objfun,needc,bigbnd,epsrf,
     *                  cdint,fdint,fdchk,fdnorm,objf,xnorm,bl,bu,c,c1,
     *                  cjac,cjacu,cjdx,dx,grad,gradu,hforwd,hcntrl,x,
     *                  wrk1,wrk2,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     npchkd  does the following...
c     (1)  computes the objective and constraint values objf and c.
c     (2)  evaluates the user-provided gradients in cjacu and gradu.
c     (3)  counts the missing gradients.
c     (4)  loads the known gradients into grad and cjac.
c     (5)  checks that the known gradients are programmed correctly.
c     (6)  computes the missing gradient elements.
c----------------------------------------------------------------------
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)

      double precision  bigbnd, cdint, epsrf, fdchk, fdint, fdnorm,
     *                  objf, xnorm
      integer           inform, ldcj, ldcju, lenw, lvlder, msgnp, n,
     *                  ncnln, nfun, ngrad, nstate

      double precision  bl(n), bu(n), c(*), c1(*), cjac(ldcj,*),
     *                  cjacu(ldcju,*), cjdx(*), dx(n), grad(n),
     *                  gradu(n), hcntrl(*), hforwd(*), user(*),
     *                  w(lenw), wrk1(n+ncnln), wrk2(n+ncnln), x(n)
      integer           iuser(*), needc(*)
c     .. subroutine arguments ..
      external          confun, objfun
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lfdset, lines1, lines2, lvldif,
     *                  lvrfyc, ncdiff, nfdiff, nout
c     .. arrays in common ..
      integer           jverfy(4)

      integer           i, infog, infoj, j, mode, ncset
      logical           centrl, needfd

      character*80      rec(3)
c     .. external subroutines ..
      external          dcopy, npfd , chcjac, chkgrd, chfd, iload ,
     *                  sload, smcopy, smload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg014/lvldif, ncdiff, nfdiff, lfdset
      common            /ngg006/epspt3, epspt5, epspt8, epspt9
      common            /ngg015/lvrfyc, jverfy
c----------------------------------------------------------------------
c
      infog = 0
      infoj = 0
      nfdiff = 0
      ncdiff = 0
      ncset = n*ncnln

      if (ncnln.gt.0) then

c        compute the constraints and jacobian matrix.

c        if some derivatives are missing, load the jacobian with dummy
c        values.  any elements left unaltered after the call to confun
c        must be estimated.  a record of the missing jacobian elements
c        is stored in  cjacu.
c
         needfd = lvlder .eq. 0 .or. lvlder .eq. 1
c
         if (needfd) call smload('general',ncnln,n,rdummy,rdummy,cjacu,
     *                           ldcju)

         call iload (ncnln,(1),needc,1)

         mode = 2
         call confun(mode,ncnln,n,ldcju,needc,x,c,cjacu,nstate,iuser,
     *               user)
         if (mode.lt.0) go to 80

         call smcopy('general',ncnln,n,cjacu,ldcju,cjac,ldcj)

         if (needfd) then

c           count the number of missing jacobian elements.

            do 40 j = 1, n
               do 20 i = 1, ncnln
                  if (cjacu(i,j).eq.rdummy) ncdiff = ncdiff + 1
   20          continue
   40       continue

            ncset = ncset - ncdiff
            if (nstate.eq.1) then
               if (ncdiff.eq.0) then
                  if (lvlder.eq.0) lvlder = 2
                  if (lvlder.eq.1) lvlder = 3
               end if
            end if
         end if
      end if

c     repeat the procedure above for the objective function.

      needfd = lvlder .eq. 0 .or. lvlder .eq. 2

      if (needfd) call sload(n,rdummy,gradu,1)
c                                 output the initial value
      iuser(2) = 1

      mode = 2
      call objfun(mode,n,x,objf,gradu,nstate,iuser,user)
      if (mode.lt.0) go to 80

      call dcopy(n,gradu,1,grad,1)
c                                 shut off output for derivative
c                                 evaluation
      iuser(2) = 0

      if (needfd) then

c        count the number of missing gradient elements.

         do 60 j = 1, n
            if (gradu(j).eq.rdummy) nfdiff = nfdiff + 1
   60    continue

         if (nstate.eq.1) then
            if (nfdiff.eq.0) then
               if (lvlder.eq.0) lvlder = 1
               if (lvlder.eq.2) lvlder = 3
            end if
         end if
      end if

      nfun = nfun + 1
      ngrad = ngrad + 1

c     check whatever gradient elements have been provided.

      if (lvrfyc.ge.0) then
         if (ncset.gt.0) then
            call chcjac(mode,lvlder,msgnp,ncset,n,ncnln,ldcj,ldcju,
     *                  bigbnd,epsrf,epspt3,fdchk,xnorm,confun,needc,bl,
     *                  bu,c,c1,cjac,cjacu,cjdx,dx,wrk2,x,wrk1,iuser,
     *                  user)
            if (mode.lt.0) go to 80
            infoj = mode
         end if

         if (nfdiff.lt.n) then
            call chkgrd(mode,msgnp,n,bigbnd,epsrf,epspt3,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,wrk1,iuser,
     *                  user)
            if (mode.lt.0) go to 80
            infog = mode
         end if
      end if

      needfd = ncdiff .gt. 0 .or. nfdiff .gt. 0
      if (needfd) then

c        compute the missing gradient elements.

         call chfd(mode,msgnp,lvlder,n,ncnln,ldcj,ldcju,bigbnd,epsrf,
     *               fdnorm,objf,confun,objfun,needc,bl,bu,c,c1,cjdx,
     *               cjac,cjacu,grad,gradu,hforwd,hcntrl,x,dx,iuser,
     *               user)

         if (mode.lt.0) go to 80

         if (lfdset.gt.0) then
            centrl = lvldif .eq. 2
            call npfd (centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,cjdx,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)
c
            if (mode.lt.0) go to 80
         end if
      end if

      inform = infoj + infog
      return

c     the user requested termination.

   80 inform = mode
c                                 end of npchkd
      end

      subroutine cmprnt(msglvl,nfree,lda,n,nclin,ncnln,nctotl,bigbnd,
     *                  named,names,lennam,nactiv,istate,kactiv,kx,a,bl,
     *                  bu,c,clamda,rlamda,x)
c----------------------------------------------------------------------
c     cmprnt   creates the expanded lagrange multiplier vector clamda.
c----------------------------------------------------------------------
      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)

      double precision  bigbnd
      integer           lda, lennam, msglvl, n, nactiv, nclin, ncnln,
     *                  nctotl, nfree
      logical           named

      double precision  a(lda,*), bl(nctotl), bu(nctotl), c(*),
     *                  clamda(nctotl), rlamda(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)

      double precision  b1, b2, res, res2, rlam, v, wlam
      integer           ip, is, j, k, nfixed, nplin, nz
c     .. external subroutines ..
      external          sload

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2
      common            /ngg4fb/icmdbg, cmdbg
c----------------------------------------------------------------------
      nplin = n + nclin
      nz = nfree - nactiv

c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      call sload(nctotl,zero,clamda,1)
      nfixed = n - nfree
      do 20 k = 1, nactiv + nfixed
         if (k.le.nactiv) then
            j = kactiv(k) + n
            rlam = rlamda(nactiv-k+1)
         else
            j = kx(nz+k)
            rlam = rlamda(k)
         end if
         clamda(j) = rlam
   20 continue
c                                 end of cmprnt
      end

      subroutine lsbnds(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldzy,lda,
     *                  ldr,ldt,istate,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)
c----------------------------------------------------------------------
c     lsbnds updates the factor r as kx is reordered to reflect the
c     status of the bound constraints given by istate.  kx is reordered
c     so that the fixed variables come last.  one of two alternative
c     are used to reorder kx. one method needs fewer accesses to kx, the
c     other gives a matrix rz with more rows and columns.
c----------------------------------------------------------------------
      double precision  condmx
      integer           inform, lda, ldr, ldt, ldzy, msglvl, n, nfree,
     *                  ngq, nrank, nres, nz
      logical           unitq

      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer           istate(*), kx(n)

      integer           iadd, ifix, j, j2, jadd, k, l, lstart, nactv,
     *                  nfixed
c     .. external subroutines ..
      external          nggnbu, lsadd 
c----------------------------------------------------------------------
      nfixed = n - nfree

      if (nrank.lt.n .and. nrank.gt.0) then

c        r is specified but singular.  try and keep the dimension of rz
c        as large as possible.

         nactv = 0
         nfree = n
         nz = n
c
         j = n
c        +       while (j .gt. 0  .and.  n-nfree .lt. nfixed) do
   20    if (j.gt.0 .and. n-nfree.lt.nfixed) then
            if (istate(j).gt.0) then
               jadd = j
               do 40 ifix = nfree, 1, -1
                  if (kx(ifix).eq.jadd) go to 60
   40          continue
c
c              add bound jadd.
c
   60          call lsadd (unitq,inform,ifix,iadd,jadd,nactv,nz,nfree,
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
                  if (nrank.gt.0) call nggnbu(n,nres,nrank,ldr,k,l,r,
     *                                 res,c,s)
               end if
  120       continue
c
         end if
         nz = nfree
      end if
c                                 end of lsbnds
      end

      subroutine cmdgen(job,msglvl,n,nclin,nmoved,iter,numinf,istate,
     *                  bigbnd,ax,bl,bu,featol,featlu,x)
c----------------------------------------------------------------------
c     cmdgen does most of the manoeuvres associated with degeneracy.
c     the degeneracy-resolving strategy operates in the following way.
c
c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     after kdegen consecutive iterations, variables within
c     featol of their bounds are set exactly on their bounds and x is
c     recomputed to satisfy the general constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.
c
c     featlu  is the array of user-supplied feasibility tolerances.
c     featol  is the array of current feasibility tolerances.
c
c     if job = 'i', cmdgen initializes the parameters in
c     common block ngg005:
c
c     tolx0   is the minimum (scaled) feasibility tolerance.
c     tolx1   is the maximum (scaled) feasibility tolerance.
c     tolinc  is the scaled increment to the current featol.
c     kdegen  is the expand frequency (specified by the user).
c             it is the frequency of resetting featol to (scaled) tolx0.
c     ndegen  counts the number of degenerate steps (incremented
c             by cmchzr).
c     itnfix  is the last iteration at which a job = 'e' or 'o' entry
c             caused an x to be put on a constraint.
c     nfix(j) counts the number of times a job = 'o' entry has
c             caused the variables to be placed on the working set,
c             where j=1 if infeasible, j=2 if feasible.
c
c     tolx0*featlu and tolx1*featlu are both close to the feasibility
c     tolerance featlu specified by the user.  (they must both be less
c     than featlu.)
c
c
c     if job = 'e',  cmdgen has been called after a cycle of kdegen
c     iterations.  constraints in the working set are examined to see if
c     any are off their bounds by an amount approaching featol.  nmoved
c     returns how many.  if nmoved is positive,  x  is moved onto the
c     constraints in the working set.  it is assumed that the calling
c     routine will then continue iterations.
c
c
c     if job = 'o',  cmdgen is being called after a subproblem has been
c     judged optimal, infeasible or unbounded.  constraint violations
c     are examined as above.
c
c     cmdgen is based on minos 5.2 routine m5dgen,
c----------------------------------------------------------------------

      double precision  zero, point6
      parameter         (zero=0.0d+0,point6=0.6d+0)

      double precision  bigbnd
      integer           iter, msglvl, n, nclin, nmoved, numinf
      character         job

      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featlu(n+nclin), featol(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      double precision  tolinc, tolx0
      integer           iprint, isumm, itnfix, kdegen, lines1, lines2,
     *                  ndegen, nout

      integer           nfix(2)

      double precision  d, epsmch, tolx1, tolz
      integer           is, j, maxfix
      character*80      rec

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg005/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix

      save              tolz
c----------------------------------------------------------------------
c
      nmoved = 0
      if (job.eq.'i' .or. job.eq.'i') then

c        job = 'initialize'.
c        initialize at the start of each linear problem.
c        kdegen  is the expand frequency      and
c        featlu  are the user-supplied feasibility tolerances.
c        they are not changed.

         epsmch = wmach(3)
c
         ndegen = 0
         itnfix = 0
         nfix(1) = 0
         nfix(2) = 0
         tolx0 = 0.5d+0
         tolx1 = 0.99d+0
         tolz = epsmch**point6
c
         if (kdegen.lt.9999999) then
            tolinc = (tolx1-tolx0)/kdegen
         else
            tolinc = zero
         end if
c
         do 20 j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
   20    continue
      else

c        job = 'end of cycle' or 'optimal'.
c        initialize local variables maxfix and tolz.

         maxfix = 2
c
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
c
            if (nfix(j).ge.maxfix) return
            nfix(j) = nfix(j) + 1
         end if
c
c        reset featol to its minimum value.
c
         do 40 j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
   40    continue
c
c        count the number of times a variable is moved a nontrivial
c        distance onto its bound.
c
         itnfix = iter
c
         do 60 j = 1, n
            is = istate(j)
            if (is.gt.0 .and. is.lt.4) then
               if (is.eq.1) then
                  d = abs(x(j)-bl(j))
               else
                  d = abs(x(j)-bu(j))
               end if
c
               if (d.gt.tolz) nmoved = nmoved + 1
            end if
   60    continue

      end if

c                                 end of cmdgen

      end

      subroutine npcrsh(cold,n,nclin,ncnln,nctotl,nactiv,nfree,nz,
     *                  istate,kactiv,bigbnd,tolact,bl,bu,c)

c     npcrsh  adds indices of nonlinear constraints to the initial
c     working set.

      integer           ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)

      double precision  bigbnd, tolact
      integer           n, nactiv, nclin, ncnln, nctotl, nfree, nz
      logical           cold

      double precision  bl(nctotl), bu(nctotl), c(*)
      integer           istate(nctotl), kactiv(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg

      integer           inpdbg(ldbg)

      double precision  b1, b2, biglow, bigupp, cmin, res, resl, resu,
     *                  toobig
      integer           i, imin, is, j, linact, nfixed, nlnact, nplin

      character*80      rec(4)

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg002/inpdbg, npdbg

      nfixed = n - nfree
      linact = nactiv
      nplin = n + nclin
c
c     if a cold start is being made, initialize the status of the qp
c     working set.  first,  if  bl(j) = bu(j),  set istate(j)=3.
c
      if (cold) then
         do 20 j = nplin + 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
   20    continue
      end if
c
c     increment nactiv and kactiv.
c     ensure that the number of bounds and general constraints in the
c     qp  working set does not exceed n.
c
      do 40 j = nplin + 1, nctotl
         if (nfixed+nactiv.eq.n) istate(j) = 0
         if (istate(j).gt.0) then
            nactiv = nactiv + 1
            kactiv(nactiv) = j - n
         end if
   40 continue
c
      if (cold) then
c

c        if a cold start is required, an attempt is made to add as many
c        nonlinear constraints as possible to the working set.

c        the following loop finds the most violated constraint.  if
c        there is room in kactiv, it will be added to the working set
c        and the process will be repeated.
c
c
         is = 1
         biglow = -bigbnd
         bigupp = bigbnd
         toobig = tolact + tolact
c
c        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
   60    if (is.gt.0 .and. nfixed+nactiv.lt.n) then
            is = 0
            cmin = tolact
c
            do 80 i = 1, ncnln
               j = nplin + i
               if (istate(j).eq.0) then
                  b1 = bl(j)
                  b2 = bu(j)
                  resl = toobig
                  resu = toobig
                  if (b1.gt.biglow) resl = abs(c(i)-b1)/(one+abs(b1))
                  if (b2.lt.bigupp) resu = abs(c(i)-b2)/(one+abs(b2))
                  res = min(resl,resu)
                  if (res.lt.cmin) then
                     cmin = res
                     imin = i
                     is = 1
                     if (resl.gt.resu) is = 2
                  end if
               end if
   80       continue
c
            if (is.gt.0) then
               nactiv = nactiv + 1
               kactiv(nactiv) = nclin + imin
               j = nplin + imin
               istate(j) = is
            end if
            go to 60
c           end while
         end if
      end if

c     an initial working set has been selected.

      nlnact = nactiv - linact
      nz = nfree - nactiv

c                                 end of npcrsh

      end

      subroutine nprset(unitq,n,nfree,nz,nq,nrowr,iperm,kx,gq,r,zy,work,
     *                  qrwork)

c     nprset  bounds the condition estimator of the transformed hessian.
c     on exit, r is of the form
c                  ( drz   0     )
c                  (  0  sigma*i )
c     where d is a diagonal matrix such that drz has a bounded condition
c     number,  i is the identity matrix and sigma  is the geometric mean
c     of the largest and smallest elements of drz. the qr factorization
c     with interchanges is used to give diagonals of drz that are
c     decreasing in modulus.
c

c

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)

      integer           n, nfree, nq, nrowr, nz
      logical           unitq

      double precision  gq(n), qrwork(2*n), r(nrowr,*), work(n),
     *                  zy(nq,*)
      integer           iperm(n), kx(n)
c     .. scalars in common ..
      double precision  drmax, drmin, rcndbd, rfrobn
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)

      double precision  drgm, drgs, gjmax, scle, sumsq
      integer           info, j, jmax, jsave, nrank
c     .. external functions ..
      double precision  snorm
      integer           isrank
      external          snorm, isrank
c     .. external subroutines ..
      external          dswap, sgeqrp, sload, sssq 
c     .. intrinsic functions ..
      intrinsic         abs, dble, sqrt

      common            /ngg018/rcndbd, rfrobn, drmax, drmin
      common            /ngg002/inpdbg, npdbg
c----------------------------------------------------------------------
c

c     bound the condition estimator of q'hq.

      if (nz.gt.1) then

c        refactorize rz.  interchanges are used to give diagonals
c        of decreasing magnitude.

         do 20 j = 1, nz - 1
            call sload(nz-j,zero,r(j+1,j),1)
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
         drgm = half*sqrt(abs(r(1,1)*r(nrank,nrank)))
         drgs = abs(r(1,1))/rcndbd
c
         if (nz.gt.nrank) then
            do 60 j = nrank + 1, nz
               call sload(j-1,zero,r(1,j),1)
   60       continue
            call sload(nz-nrank,drgs,r(nrank+1,nrank+1),nrowr+1)
         end if
      end if
c

c     reset the range-space partition of the hessian.

      if (nz.lt.n) then
         do 80 j = nz + 1, n
            call sload(j,zero,r(1,j),1)
   80    continue
         call sload(n-nz,drgm,r(nz+1,nz+1),nrowr+1)
      end if
c
c     recompute the frobenius norm of r.
c
      scle = sqrt(dble(n-nz))*drgm
      sumsq = one
      do 100 j = 1, nz
         call sssq (j,r(1,j),1,scle,sumsq)
  100 continue
      rfrobn = snorm(scle,sumsq)

c                                 end of nprset

      end

      subroutine cmcrsh(start,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *                  lda,istate,kactiv,kx,bigbnd,tolact,a,ax,bl,bu,
     *                  featol,x,wx,work)

c     cmcrsh  computes the quantities  istate (optionally),  kactiv,
c     nactiv,  nz  and  nfree  associated with the working set at x.
c
c     the computation depends upon the value of the input parameter
c     start,  as follows...
c
c     start = 'cold'  an initial working set will be selected. first,
c                     nearly-satisfied or violated bounds are added.
c                     next,  general linear constraints are added that
c                     have small residuals.
c
c     start = 'warm'  the quantities kactiv, nactiv and nfree are
c                     initialized from istate,  specified by the user.
c
c     if vertex is true, an artificial vertex is defined by fixing some
c     variables on their bounds.  infeasible variables selected for the
c     artificial vertex are fixed at their nearest bound.  otherwise,
c     the variables are unchanged.
c
c     values of istate(j)....
c
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c

c

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      double precision  bigbnd, tolact
      integer           lda, n, nactiv, nartif, nclin, nctotl, nfree
      logical           vertex
      character*4       start

      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), work(n), wx(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg

      integer           icmdbg(ldbg)

      double precision  b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, tol, toobig
      integer           i, imin, is, j, jfix, jfree, jmin, k

      character*80      rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. external subroutines ..
      external          dcopy

      common            /ngg4nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ngg4fb/icmdbg, cmdbg

      flmax = wmach(7)
      biglow = -bigbnd
      bigupp = bigbnd


c     move the variables inside their bounds.


      do 20 j = 1, n
         b1 = bl(j)
         b2 = bu(j)
         tol = featol(j)

         if (b1.gt.biglow) then
            if (x(j).lt.b1-tol) x(j) = b1
         end if

         if (b2.lt.bigupp) then
            if (x(j).gt.b2+tol) x(j) = b2
         end if
   20 continue

      call dcopy(n,x,1,wx,1)

      nfree = n
      nactiv = 0
      nartif = 0

      if (start.eq.'cold') then
         do 40 j = 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
   40    continue
c
      else if (start.eq.'warm') then
         do 60 j = 1, nctotl
            if (istate(j).gt.3 .or. istate(j).lt.0) istate(j) = 0
            if (bl(j).ne.bu(j) .and. istate(j).eq.3) istate(j) = 0
   60    continue
      end if
c
c     define nfree and kactiv.
c     ensure that the number of bounds and general constraints in the
c     working set does not exceed n.
c
      do 80 j = 1, nctotl
         if (nactiv.eq.nfree) istate(j) = 0
c
         if (istate(j).gt.0) then
            if (j.le.n) then
               nfree = nfree - 1
c
               if (istate(j).eq.1) then
                  wx(j) = bl(j)
               else if (istate(j).ge.2) then
                  wx(j) = bu(j)
               end if
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
   80 continue
c

c     if a cold start is required,  attempt to add as many
c     constraints as possible to the working set.

      if (start.eq.'cold') then
c
c        see if any bounds are violated or nearly satisfied.
c        if so,  add these bounds to the working set and set the
c        variables exactly on their bounds.
c
         j = n
c        +       while (j .ge. 1  .and.  nactiv .lt. nfree) do
  100    if (j.ge.1 .and. nactiv.lt.nfree) then
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
                  nfree = nfree - 1
               end if
            end if
            j = j - 1
            go to 100
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
         if (nclin.gt.0 .and. nactiv.lt.nfree) then
            do 120 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot(n,a(i,1),lda,wx,1)
  120       continue
c
            is = 1
            toobig = tolact + tolact
c
c           +          while (is .gt. 0  .and.  nactiv .lt. nfree) do
  140       if (is.gt.0 .and. nactiv.lt.nfree) then
               is = 0
               resmin = tolact
c
               do 160 i = 1, nclin
                  j = n + i
                  if (istate(j).eq.0) then
                     b1 = bl(j)
                     b2 = bu(j)
                     resl = toobig
                     resu = toobig
                     if (b1.gt.biglow) resl = abs(ax(i)-b1)/(one+abs(b1)
     *                                        )
                     if (b2.lt.bigupp) resu = abs(ax(i)-b2)/(one+abs(b2)
     *                                        )
                     residl = min(resl,resu)
                     if (residl.lt.resmin) then
                        resmin = residl
                        imin = i
                        is = 1
                        if (resl.gt.resu) is = 2
                     end if
                  end if
  160          continue
c
               if (is.gt.0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j = n + imin
                  istate(j) = is
               end if
               go to 140
c              +          end while
            end if
         end if
      end if
c
      if (vertex .and. nactiv.lt.nfree) then

c        find an initial vertex by temporarily fixing some variables.

c        compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).
c
         do 200 j = 1, n
            if (istate(j).eq.0) then
               colsiz = zero
               do 180 k = 1, nclin
                  if (istate(n+k).gt.0) colsiz = colsiz + abs(a(k,j))
  180          continue
               work(j) = colsiz
            end if
  200    continue
c
c        find the  nartif  smallest such columns.
c        this is an expensive loop.  later we can replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).

c        +       while (nactiv .lt. nfree) do
  220    if (nactiv.lt.nfree) then
            colmin = flmax
            do 240 j = 1, n
               if (istate(j).eq.0) then
                  if (nclin.eq.0) go to 260
                  colsiz = work(j)
                  if (colmin.gt.colsiz) then
                     colmin = colsiz
                     jmin = j
                  end if
               end if
  240       continue
            j = jmin
c
c           fix x(j) at its current value.
c
  260       istate(j) = 4
            nartif = nartif + 1
            nfree = nfree - 1
            go to 220
c           +       end while
         end if
      end if
c
      jfree = 1
      jfix = nfree + 1
      do 280 j = 1, n
         if (istate(j).le.0) then
            kx(jfree) = j
            jfree = jfree + 1
         else
            kx(jfix) = j
            jfix = jfix + 1
         end if
  280 continue

c                                 end of cmcrsh

      end

      subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )

      integer            incx, lda, n
      character*1        diag, trans, uplo

      double precision   a( lda, * ), x( * )

c  dtrmv  does one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,
c
c  where x is n element vector and a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.

c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be doed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := a'*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.

c  level 2 blas routine.

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit

      if( n.eq.0 ) return
c
      nounit = (diag.eq.'n')

c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.

      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if( trans.eq.'n' )then

c        form  x := a*x.

         if( uplo.eq.'u' )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     *                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     *                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     *                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     *                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x.
c
         if( uplo.eq.'u' )then
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  if( nounit )
     *               temp = temp*a( j, j )
                  do 90, i = j - 1, 1, -1
                     temp = temp + a( i, j )*x( i )
   90             continue
                  x( j ) = temp
  100          continue
            else
               jx = kx + ( n - 1 )*incx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     *               temp = temp*a( j, j )
                  do 110, i = j - 1, 1, -1
                     ix   = ix   - incx
                     temp = temp + a( i, j )*x( ix )
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  if( nounit ) temp = temp*a( j, j )
                  do 130, i = j + 1, n
                     temp = temp + a( i, j )*x( i )
  130             continue
                  x( j ) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( nounit ) temp = temp*a( j, j )
                  do 150, i = j + 1, n
                     ix   = ix   + incx
                     temp = temp + a( i, j )*x( ix )
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if
c                                 end of dtrmv.
      end

      integer function idamax( n, x, incx )


      integer                  incx, n

      double precision         x( * )

c  idamax returns the smallest value of i such that

c     abs( x( i ) ) = max( abs( x( j ) ) )


      double precision         xmax
      integer                  i, imax, ix
c----------------------------------------------------------------------
      if( n.gt.0 )then
         imax = 1
         if( n.gt.1 )then
            xmax = abs( x( 1 ) )
            ix   = 1
            do 10, i = 2, n
               ix = ix + incx
               if( xmax.lt.abs( x( ix ) ) )then
                  xmax = abs( x( ix ) )
                  imax = i
               end if
   10       continue
         end if
      else
         imax = 0
      end if

      idamax = imax

c                                 end of idamax

      end

      subroutine daxpy ( n, alpha, x, incx, y, incy )


      double precision   alpha
      integer            incx, incy, n

      double precision   x( * ), y( * )

c  daxpy  does the operation
c
c     y := alpha*x + y
c

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      integer            i, ix, iy

c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( alpha.ne.zero )then
            if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
               do 10, ix = 1, 1 + ( n - 1 )*incx, incx
                  y( ix ) = alpha*x( ix ) + y( ix )
   10          continue
            else
               if( incy.ge.0 )then
                  iy = 1
               else
                  iy = 1 - ( n - 1 )*incy
               end if
               if( incx.gt.0 )then
                  do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                     y( iy ) = alpha*x( ix ) + y( iy )
                     iy      = iy            + incy
   20             continue
               else
                  ix = 1 - ( n - 1 )*incx
                  do 30, i = 1, n
                     y( iy ) = alpha*x( ix ) + y( iy )
                     ix      = ix            + incx
                     iy      = iy            + incy
   30             continue
               end if
            end if
         end if
      end if

c                                 end of daxpy

      end

      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,
     *                   beta, y, incy )


      double precision   alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans

      double precision   a( lda, * ), x( * ), y( * )

c  dgemv  does one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be doed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - double precision.
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - double precision.
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - double precision array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.


      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   temp, temp1, temp2, temp3, temp4
      integer            i, iy, j, jx, kx, ky, lenx, leny, m4, n4
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     *    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     *   return
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( trans.eq.'n' )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the inner loops are all
c     equivalent to axpy operations.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero ) return
      jx = kx
      if( trans.eq.'n' )then

c        form  y := alpha*a*x + y.

         if( incy.eq.1 )then
c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*( n/4 )
            do 60, j = 1, n4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     *             temp4.ne.zero )then
                  do 50, i = 1, m
                     y( i ) = ( ( ( ( y( i ) + temp1*a( i, j ) )
     *                        + temp2*a( i, j + 1 ) )
     *                        + temp3*a( i, j + 2 ) )
     *                        + temp4*a( i, j + 3 ) )
   50             continue
               end if
               jx = jx + 4*incx
   60       continue
c**** clean-up loop ****************************************************
            do 80, j = n4 + 1, n, 1
               temp = alpha*x( jx )
               if( temp.ne.zero )then
                  do 70, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   70             continue
               end if
               jx = jx + incx
   80       continue
         else
c**** u n r o l l   t o   d e p t h   4 ********************************
            n4 = 4*( n/4 )
            do 100, j = 1, n4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     *             temp4.ne.zero )then
                  iy = ky
                  do 90, i = 1, m
                     y( iy ) = ( ( ( ( y( iy ) + temp1*a( i, j ) )
     *                         + temp2*a( i, j + 1 ) )
     *                         + temp3*a( i, j + 2 ) )
     *                         + temp4*a( i, j + 3 ) )
                     iy = iy + incy
   90             continue
               end if
               jx = jx + 4*incx
  100       continue
c**** clean-up loop ****************************************************
            do 120, j = n4 + 1, n, 1
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy = ky
                  do 110, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy = iy + incy
  110             continue
               end if
               jx = jx + incx
  120       continue
         end if
      else
c
c        form  y := alpha*a'*x + y.
c
         if( incy.eq.1 )then
c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*( m/4 )
            do 140, j = 1, m4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     *             temp4.ne.zero )then
                  do 130, i = 1, n
                     y( i ) = ( ( ( ( y( i ) + temp1*a( j, i ) )
     *                        + temp2*a( j + 1, i ) )
     *                        + temp3*a( j + 2, i ) )
     *                        + temp4*a( j + 3, i ) )
  130             continue
               end if
               jx = jx + 4*incx
  140       continue
c**** clean-up loop ****************************************************
            do 160, j = m4 + 1, m, 1
               temp = alpha*x( jx )
               if( temp.ne.zero )then
                  do 150, i = 1, n
                     y( i ) = y( i ) + temp*a( j, i )
  150             continue
               end if
               jx = jx + incx
  160       continue
         else
c**** u n r o l l   t o   d e p t h   4 ********************************
            m4 = 4*( m/4 )
            do 180, j = 1, m4, 4
               temp1 = alpha*x( jx )
               temp2 = alpha*x( jx + incx )
               temp3 = alpha*x( jx + 2*incx )
               temp4 = alpha*x( jx + 3*incx )
               if( temp1.ne.zero.or.temp2.ne.zero.or.temp3.ne.zero.or.
     *             temp4.ne.zero )then
                  iy = ky
                  do 170, i = 1, n
                     y( iy ) = ( ( ( ( y( iy ) + temp1*a( j, i ) )
     *                         + temp2*a( j + 1, i ) )
     *                         + temp3*a( j + 2, i ) )
     *                         + temp4*a( j + 3, i ) )
                     iy = iy + incy
  170             continue
               end if
               jx = jx + 4*incx
  180       continue
c**** clean-up loop ****************************************************
            do 200, j = m4 + 1, m, 1
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy = ky
                  do 190, i = 1, n
                     y( iy ) = y( iy ) + temp*a( j, i )
                     iy = iy + incy
  190             continue
               end if
               jx = jx + incx
  200       continue
         end if
      end if

c                                 end of dgemv.

      end

      subroutine dtrsv ( uplo, trans, diag, n, a, lda, x, incx )

      integer            incx, lda, n
      character*1        diag, trans, uplo

      double precision   a( lda, * ), x( * )

c  dtrsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be doed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.

c  level 2 blas routine.

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, ix, j, jx, kx
      logical            nounit

      if( n.eq.0 ) return
c
      nounit = diag.eq.'n'
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if( trans.eq.'n' )then

c        form  x := inv( a )*x.

         if( uplo.eq.'u' )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit ) x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit ) x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit ) x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit ) x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else

c        form  x := inv( a' )*x.

         if( uplo.eq.'u' )then
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  do 90, i = 1, j - 1
                     temp = temp - a( i, j )*x( i )
   90             continue
                  if( nounit ) temp = temp/a( j, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  do 110, i = 1, j - 1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit ) temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  do 130, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( i )
  130             continue
                  if( nounit ) temp = temp/a( j, j )
                  x( j ) = temp
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  do 150, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit ) temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   - incx
  160          continue
            end if
         end if
      end if
c                                 end of dtrsv.
      end

      subroutine dswap ( n, x, incx, y, incy )


      integer            incx, incy, n

      double precision   x( * ), y( * )
c     ..
c
c dswap does the operations
c
c     temp := x,   x := y,   y := temp.

      double precision   temp
      integer            i, ix, iy
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               temp    = x( iy )
               x( iy ) = y( iy )
               y( iy ) = temp
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  temp    = x( ix )
                  x( ix ) = y( iy )
                  y( iy ) = temp
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  temp    = x( ix )
                  x( ix ) = y( iy )
                  y( iy ) = temp
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if

c                                 end of dswap

      end

      double precision function ddot ( n, x, incx, y, incy )

      integer                           incx, incy, n

      double precision                  x( * ), y( * )

c  ddot returns the value
c
c     ddot = x'y

      double precision      zero
      parameter           ( zero = 0.0d+0 )

      double precision      sum
      integer               i, ix, iy
c----------------------------------------------------------------------
      sum = zero
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               sum = sum + x( ix )*y( ix )
   10       continue
         else
            if( incy.ge.0 )then
               iy = 1
            else
               iy = 1 - ( n - 1 )*incy
            end if
            if( incx.gt.0 )then
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                  sum = sum + x( ix )*y( iy )
                  iy  = iy  + incy
   20          continue
            else
               ix = 1 - ( n - 1 )*incx
               do 30, i = 1, n
                  sum = sum + x( ix )*y( iy )
                  ix  = ix  + incx
                  iy  = iy  + incy
   30          continue
            end if
         end if
      end if
c
      ddot = sum

c                                 end of ddot

      end

      subroutine dger ( m, n, alpha, x, incx, y, incy, a, lda )

      double precision   alpha
      integer            incx, incy, lda, m, n

      double precision   a( lda, * ), x( * ), y( * )

c  dger   does the rank 1 operation
c
c     a := alpha*x*y' + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c
c  parameters
c  ==========
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - double precision.
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
c           ( 1 + ( m - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the m
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - double precision array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients. on exit, a is
c           overwritten by the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.

c  level 2 blas routine.

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, ix, j, jy, kx

      if(m.eq.0.or.n.eq.0.or.( alpha.eq.zero ) ) return

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if

c                                 end of dger.

      end

      subroutine dcopy ( n, x, incx, y, incy )

      integer            incx, incy, n

      double precision   x( * ), y( * )

c  dcopy does the operation

c     y := x

      integer            i, ix, iy

      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( iy )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  y( iy ) = x( ix )
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  y( iy ) = x( ix )
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if

c                                 end of dcopy

      end


      double precision function dnrm2 ( n, x, incx )

      integer                           incx, n

      double precision                  x( * )

c  dnrm2 returns the euclidean norm of a vector via the function
c  name, so that
c
c     dnrm2  := sqrt( x'*x )

      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      double precision      norm, scale, ssq
c     .. external functions ..
      double precision      snorm
      external              snorm
c     .. external subroutines ..
      external              sssq 
c----------------------------------------------------------------------
      if( n.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
         call sssq ( n, x, incx, scale, ssq )
         norm  = snorm( scale, ssq )
      end if
c
      dnrm2  = norm

c                                 end of dnrm2

      end

      subroutine dscal ( n, alpha, x, incx )

      double precision   alpha
      integer            incx, n

      double precision   x( * )

c  dscal does the operation

c     x := alpha*x

      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      integer            ix
c     ..
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( alpha.eq.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   10       continue
         else if( alpha.eq.( -one ) )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = -x( ix )
   20       continue
         else if( alpha.ne.one )then
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = alpha*x( ix )
   30       continue
         end if
      end if

c                                 end of dscal

      end

      subroutine sgeqrp(pivot,m,n,a,lda,zeta,perm,work,ifail)

c  1. purpose
c     =======
c
c  sgeqrp  finds  a  qr factorization  of  the  real  m by n  matrix  a,
c  incorporating  column interchanges,  so that  a  is reduced to  upper
c  triangular form  by means of  orthogonal transformations  and  column
c  permutations.
c
c  2. description
c     ===========
c
c  the m by n matrix a is factorized as
c
c     a = q*( r )*p'      when   m.gt.n,
c           ( 0 )
c
c     a = q*r*p'          when   m = n,
c
c     a = q*( r  x )*p'   when   m.lt.n,
c
c  where  q  is an  m by m  orthogonal matrix, r  is  a  min( m, n )  by
c  min( m, n )  upper triangular matrix and  p is an  n by n permutation
c  matrix.
c
c  the  factorization  is  obtained  by  householder's  method. the  kth
c  transformation matrix, q( k ), which is used  to introduce zeros into
c  the kth column of a is given in the form
c
c     q( k ) = ( i     0   ),
c              ( 0  t( k ) )
c
c  where
c
c     t( k ) = i - u( k )*u( k )',
c
c     u( k ) = ( zeta( k ) ),
c              (    z( k ) )
c
c  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
c  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
c  triangular part of  a.
c
c  the vector  u( k ) is returned in the kth element of  zeta and in the
c  kth column of a, such that zeta( k ) is in zeta( k ) and the elements
c  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  the elements of  r
c  are returned in the upper triangular part of  a.
c
c  q is given by
c
c     q = ( q( p )*q( p - 1 )*...*q( 1 ) )',   p = min( m, n ).
c
c  two options are available for the column permutations. in either case
c  the column for which the  sub-diagonal elements are to be annihilated
c  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
c  the  particular column chosen as the pivot column is either that  for
c  which  the  unreduced  part  ( elements k onwards )  has the  largest
c  euclidean  length, or  is that for  which the ratio of the  euclidean
c  length  of the  unreduced part  to the  euclidean length of the whole
c  column is a maximum.
c
c  3. parameters
c     ==========
c
c  pivot  - character*1.
c
c           on  entry, pivot  specifies  the  pivoting  strategy  to  be
c           doed as follows.
c
c           pivot = 'c' or 'c'   ( column interchanges )
c
c              column  interchanges  are  to be  incorporated  into  the
c              factorization, such that the  column whose unreduced part
c              has  maximum  euclidean  length  is chosen  as the  pivot
c              column at each step.
c
c           pivot = 's' or 's'   ( scaled column interchanges )
c
c              scaled  column interchanges  are to be  incorporated into
c              the  factorization, such  that the  column for which  the
c              ratio  of the  euclidean  length of the unreduced part of
c              the column to the original euclidean length of the column
c              is a maximum is chosen as the  pivot column at each step.
c
c           unchanged on exit.
c
c  m      - integer.
c
c           on entry, m  must specify the number of rows of a. m must be
c           at  least  zero. when  m = 0  then  an  immediate return  is
c           effected.
c
c           unchanged on exit.
c
c  n      - integer.
c
c           on entry, n  must specify the number of columns of a. n must
c           be  at least zero. when  n = 0  then an immediate return  is
c           effected.
c
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c
c           before entry, the leading  m by n  part of the array  a must
c           contain the matrix to be factorized.
c
c           on  exit, the  min( m, n ) by min( m, n )  upper  triangular
c           part of a will contain the upper triangular matrix r and the
c           m by min( m, n )  strictly lower triangular part of  a  will
c           contain details  of the  factorization  as  described above.
c           when m.lt.n then the remaining m by ( n - m ) part of a will
c           contain the matrix x.
c
c  lda    - integer.
c
c           on entry, lda  must  specify  the  leading dimension of  the
c           array  a  as declared in the calling (sub) program. lda must
c           be at least  m.
c
c           unchanged on exit.
c
c  zeta   - real             array of dimension at least ( n ).
c
c           on exit,  zeta( k )  contains the scalar  zeta( k )  for the
c           kth  transformation.  if  t( k ) = i  then  zeta( k ) = 0.0,
c           otherwise  zeta( k )  contains  zeta( k ) as described above
c           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
c           when n.gt.m the elements  zeta( m + 1 ), zeta( m + 2 ), ...,
c           zeta( n )  are used as internal workspace.
c
c  perm   - integer array of dimension at least  min( m, n ).
c
c           on exit, perm  contains details of the permutation matrix p,
c           such  that  perm( k ) = k  if no  column interchange occured
c           at  the  kth  step  and  perm( k ) = j, ( k .lt. j .le. n ),
c           if columns  k  and  j  were  interchanged  at the  kth step.
c           note that there are  min( m, n ) permutations.
c
c  work   - real array of dimension at least ( 2*n ).
c
c           used as internal workspace.
c
c           on exit, work( j ), j = 1, 2, ..., n, contains the euclidean
c           length of the jth column of the permuted matrix a*p'.
c
c  ifail  - integer.
c
c           before entry,  ifail  must contain one of the values -1 or 0
c           or 1 to specify noisy soft failure or noisy hard failure  or
c           silent soft failure.
c
c           on  successful exit, ifail  will be  zero,  otherwise  ifail
c           will  be set to   -1  indicating that an input parameter has
c           been  incorrectly supplied. see the next section for further
c           details.
c
c  4. diagnostic information
c     ======================
c
c  ifail = -1
c
c     one or more of the following conditions holds:
c
c        pivot .ne. 'c' or 'c' or 's' or 's'
c        m     .lt. 0
c        n     .lt. 0
c        lda   .lt. m

      double precision  lamda, one, zero
      parameter         (lamda=1.0d-2,one=1.0d+0,zero=0.0d+0)


      integer           ifail, lda, m, n
      character*1       pivot

      double precision  a(lda,*), work(*), zeta(*)
      integer           perm(*)

      double precision  eps, maxnrm, norm, temp, tol
      integer           ierr, j, jmax, k, la

      character*46      rec(1)
c     .. external functions ..
      double precision  dnrm2, x02ajf

      external          dnrm2, x02ajf
c     .. external subroutines ..
      external          dgemv, dger, dswap, sgrfg

      double precision wmach
      common/ cstmch /wmach(10)

c
c     compute eps and the initial column norms.
c
      if (min(m,n).eq.0) call errdbg ('sgeqrp')

      eps = wmach(3)
      do 20 j = 1, n
         work(j) = dnrm2(m,a(1,j),1)
         work(j+n) = work(j)
   20 continue
c
c     do the factorization. tol is the tolerance for sgrfg.
c
      la = lda
      do 120 k = 1, min(m,n)
c
c        find the pivot column.
c
         maxnrm = zero
         jmax = k
         if ( pivot.eq.'c' ) then
            do 40 j = k, n
               if (work(j+n).gt.maxnrm) then
                  maxnrm = work(j+n)
                  jmax = j
               end if
   40       continue
         else
            do 60 j = k, n
               if (work(j).gt.zero) then
                  if (k.le.1) then
                     jmax = j
                     go to 80
                  else if ((work(j+n)/work(j)).gt.maxnrm) then
                     maxnrm = work(j+n)/work(j)
                     jmax = j
                  end if
               end if
   60       continue
   80       continue
         end if
         perm(k) = jmax
         if (jmax.gt.k) then
            call dswap(m,a(1,k),1,a(1,jmax),1)
            temp = work(k)
            work(k) = work(jmax)
            work(jmax) = temp
            work(jmax+n) = work(k+n)
         end if
         tol = eps*work(k)
         if (k.lt.m) then
c
c           use a householder reflection to zero the kth column of a.
c           first set up the reflection.
c
            call sgrfg (m-k,a(k,k),a(k+1,k),1,tol,zeta(k))
            if (k.lt.n) then
               if (zeta(k).gt.zero) then
                  if ((k+1).eq.n) la = m - k + 1
c
c                 temporarily store beta and put zeta( k ) in a( k, k ).
c
                  temp = a(k,k)
                  a(k,k) = zeta(k)
c
c                 do the operation  a := q( k )*a.
c
c                 let  b  denote  the bottom  ( m - k + 1 ) by ( n - k )
c                 part of a.

c                 first  form   work = b'*u.  ( work  is  stored  in the
c                 elements zeta( k + 1 ), ..., zeta( n ). )

                  call dgemv('transpose',m-k+1,n-k,one,a(k,k+1),la,
     *                       a(k,k),1,zero,zeta(k+1),1)

c                 form  b := b - u*work'.

                  call dger(m-k+1,n-k,-one,a(k,k),1,zeta(k+1),1,a(k,k+1)
     *                      ,la)

c                 restore beta.

                  a(k,k) = temp
               end if
c
c              update  the  unreduced  column  norms.  use  the  linpack
c              criterion for when to recompute the norms, except that we
c              retain  the original column lengths throughout  and use a
c              smaller lamda.
c
               do 100 j = k + 1, n
                  if (work(j+n).gt.zero) then
                     temp = abs(a(k,j))/work(j+n)
                     temp = max((one+temp)*(one-temp),zero)
                     norm = temp
                     temp = one + lamda*temp*(work(j+n)/work(j))**2
                     if (temp.gt.one) then
                        work(j+n) = work(j+n)*sqrt(norm)
                     else
                        work(j+n) = dnrm2(m-k,a(k+1,j),1)
                     end if
                  end if
  100          continue
            end if
         end if
  120 continue
c
c     set the final  zeta  when  m.le.n.
c
      if (m.le.n) zeta(m) = zero

c                                 end of sgeqrp

      end

      subroutine sscmv ( n, alpha, x, incx, y, incy )

      double precision   alpha
      integer            incx, incy, n

      double precision   x( * ), y( * )
c
c  sscmv  does the operation
c
c     y := alpha*x

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      integer            i, ix, iy
c     .. external subroutines ..
      external           sload
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( ( alpha.eq.zero ).and.( incy.ne.0 ) )then
            call sload( n, zero, y, abs( incy ) )
         else
            if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
               do 10, ix = 1, 1 + ( n - 1 )*incx, incx
                  y( ix ) = alpha*x( ix )
   10          continue
            else
               if( incy.ge.0 )then
                  iy = 1
               else
                  iy = 1 - ( n - 1 )*incy
               end if
               if( incx.gt.0 )then
                  do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                     y( iy ) = alpha*x( ix )
                     iy      = iy            + incy
   20             continue
               else
                  ix = 1 - ( n - 1 )*incx
                  do 30, i = 1, n
                     y( iy ) = alpha*x( ix )
                     ix      = ix            + incx
                     iy      = iy            + incy
   30             continue
               end if
            end if
         end if
      end if

c                                 end of sscmv

      end

      subroutine sgeqr (m,n,a,lda,zeta,ifail)

c  1. purpose
c     =======
c
c  sgeqr   finds  the  qr factorization  of the real  m by n,  m .ge. n,
c  matrix a,  so that  a is reduced to upper triangular form by means of
c  orthogonal transformations.
c
c  2. description
c     ===========
c
c  the m by n matrix a is factorized as
c
c     a = q*( r )   when   m.gt.n,
c           ( 0 )
c
c     a = q*r       when   m = n,
c
c  where  q  is an  m by m orthogonal matrix and  r  is an  n by n upper
c  triangular matrix.
c
c  the  factorization  is  obtained  by  householder's  method. the  kth
c  transformation matrix, q( k ), which is used  to introduce zeros into
c  the kth column of a is given in the form
c
c     q( k ) = ( i     0   ),
c              ( 0  t( k ) )
c
c  where
c
c     t( k ) = i - u( k )*u( k )',
c
c     u( k ) = ( zeta( k ) ),
c              (    z( k ) )
c
c  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
c  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
c  triangular part of  a.
c
c  the vector  u( k ) is returned in the kth element of  zeta and in the
c  kth column of a, such that zeta( k ) is in zeta( k ) and the elements
c  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  the elements of  r
c  are returned in the upper triangular part of  a.
c
c  q is given by
c
c     q = ( q( n )*q( n - 1 )*...*q( 1 ) )'.
c
c  3. parameters
c     ==========
c
c  m      - integer.
c
c           on entry, m must specify the number of rows of  a. m must be
c           at least  n.
c
c           unchanged on exit.
c
c  n      - integer.
c
c           on entry, n must specify the number of columns of  a. n must
c           be  at  least zero. when  n = 0  then an immediate return is
c           effected.
c
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c
c           before entry, the leading  m by n  part of the array  a must
c           contain the matrix to be factorized.
c
c           on exit, the  n by n upper triangular part of a will contain
c           the upper triangular matrix r and the  m by n strictly lower
c           triangular  part   of   a   will  contain  details   of  the
c           factorization as described above.
c
c  lda    - integer.
c
c           on entry, lda  must  specify  the  leading dimension of  the
c           array  a  as declared in the calling (sub) program. lda must
c           be at least  m.
c
c           unchanged on exit.
c
c  zeta   - real             array of dimension at least ( n ).
c
c           on exit,  zeta( k )  contains the scalar  zeta( k )  for the
c           kth  transformation.  if  t( k ) = i  then  zeta( k ) = 0.0,
c           otherwise  zeta( k )  contains  zeta( k ) as described above
c           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
c
c  ifail  - integer.
c
c           before entry,  ifail  must contain one of the values -1 or 0
c           or 1 to specify noisy soft failure or noisy hard failure  or
c           silent soft failure.
c
c           on successful  exit  ifail  will be  zero,  otherwise  ifail
c           will  be set to  -1  indicating that an  input parameter has
c           been  incorrectly  set. see  the  next section  for  further
c           details.
c
c  4. diagnostic information
c     ======================
c
c  ifail = -1
c
c     one or more of the following conditions holds:
c
c        m   .lt. n
c        n   .lt. 0
c        lda .lt. m

      double precision  one, zero
      parameter         (one=1.0d+0,zero=0.0d+0)

      integer           ifail, lda, m, n

      double precision  a(lda,*), zeta(*)

      double precision  temp
      integer           ierr, k, la

      character*46      rec(1)
c     .. external subroutines ..
      external          dgemv, dger, sgrfg

c     do the factorization.

      if (n.eq.0) call errdbg ('sgeqr')

      la = lda
      do 20 k = 1, min(m-1,n)
c
c        use a  householder reflection  to  zero the  kth column  of  a.
c        first set up the reflection.
c
         call sgrfg(m-k,a(k,k),a(k+1,k),1,zero,zeta(k))
         if ((zeta(k).gt.zero) .and. (k.lt.n)) then
            if ((k+1).eq.n) la = m - k + 1
c
c           temporarily  store  beta and  put  zeta( k )  in  a( k, k ).
c
            temp = a(k,k)
            a(k,k) = zeta(k)
c
c           do the operation  a := q( k )*a.
c
c           let  b  denote  the bottom  ( m - k + 1 ) by ( n - k )  part
c           of  a.
c
c           first form   work = b'*u.  ( work  is stored in the elements
c           zeta( k + 1 ), ..., zeta( n ). )
c
            call dgemv('transpose',m-k+1,n-k,one,a(k,k+1),la,a(k,k),1,
     *                 zero,zeta(k+1),1)
c
c           form  b := b - u*work'.
c
            call dger(m-k+1,n-k,-one,a(k,k),1,zeta(k+1),1,a(k,k+1),la)
c
c           restore beta.
c
            a(k,k) = temp
         end if
   20 continue
c
c     set the final  zeta  when  m.eq.n.
c
      if (m.eq.n) zeta(n) = zero

c                                 end of sgeqr

      end

      double precision function dlantr(norm,uplo,diag,m,n,a,lda,work)

c  purpose
c  =======
c
c  dlantr  returns the value of the one norm,  or the frobenius norm, or
c  the  infinity norm,  or the  element of  largest absolute value  of a
c  trapezoidal or triangular matrix a.
c
c  description
c  ===========
c
c  dlantr returns the value
c
c     dlantr = ( max(abs(a(i,j))), norm = 'm' or 'm'
c              (
c              ( norm1(a),         norm = '1', 'o' or 'o'
c              (
c              ( normi(a),         norm = 'i' or 'i'
c              (
c              ( normf(a),         norm = 'f', 'f', 'e' or 'e'
c
c  where  norm1  denotes the  one norm of a matrix (maximum column sum),
c  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
c  normf  denotes the  frobenius norm of a matrix (square root of sum of
c  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.
c
c  arguments
c  =========
c
c  norm    (input) character*1
c          specifies the value to be returned in dlantr as described
c          above.
c
c  uplo    (input) character*1
c          specifies whether the matrix a is upper or lower trapezoidal.
c          = 'u':  upper trapezoidal
c          = 'l':  lower trapezoidal
c          note that a is triangular instead of trapezoidal if m = n.
c
c  diag    (input) character*1
c          specifies whether or not the matrix a has unit diagonal.
c          = 'n':  non-unit diagonal
c          = 'u':  unit diagonal
c
c  m       (input) integer
c          the number of rows of the matrix a.  m >= 0, and if
c          uplo = 'u', m <= n.  when m = 0, dlantr is set to zero.
c
c  n       (input) integer
c          the number of columns of the matrix a.  n >= 0, and if
c          uplo = 'l', n <= m.  when n = 0, dlantr is set to zero.
c
c  a       (input) real array, dimension (lda,n)
c          the trapezoidal matrix a (a is triangular if m = n).
c          if uplo = 'u', the leading m by n upper trapezoidal part of
c          the array a contains the upper trapezoidal matrix, and the
c          strictly lower triangular part of a is not referenced.
c          if uplo = 'l', the leading m by n lower trapezoidal part of
c          the array a contains the lower trapezoidal matrix, and the
c          strictly upper triangular part of a is not referenced.  note
c          that when diag = 'u', the diagonal elements of a are not
c          referenced and are assumed to be one.
c
c  lda     (input) integer
c          the leading dimension of the array a.  lda >= max(m,1).
c
c  work    (workspace) real array, dimension (lwork),
c          where lwork >= m when norm = 'i'; otherwise, work is not
c          referenced.

c  -- lapack auxiliary routine

      double precision                 one, zero
      parameter                        (one=1.0d+0,zero=0.0d+0)

      integer                          lda, m, n
      character                        diag, norm, uplo

      double precision                 a(lda,*), work(*)

      double precision                 scale, sum, value
      integer                          i, j
      logical                          udiag
c     .. external subroutines ..
      external                         sssq 
c----------------------------------------------------------------------
      if (min(m,n).eq.0) then
         value = zero
      else if ( norm.eq.'m') then

c        find max(abs(a(i,j))).

         if ( diag.eq.'u' ) then
            value = one
            if ( uplo.eq.'u') then
               do 40 j = 1, n
                  do 20 i = 1, min(m,j-1)
                     value = max(value,abs(a(i,j)))
   20             continue
   40          continue
            else
               do 80 j = 1, n
                  do 60 i = j + 1, m
                     value = max(value,abs(a(i,j)))
   60             continue
   80          continue
            end if
         else
            value = zero
            if ( uplo.eq.'u' ) then
               do 120 j = 1, n
                  do 100 i = 1, min(m,j)
                     value = max(value,abs(a(i,j)))
  100             continue
  120          continue
            else
               do 160 j = 1, n
                  do 140 i = j, m
                     value = max(value,abs(a(i,j)))
  140             continue
  160          continue
            end if
         end if
      else if (norm.eq.'o' .or. norm.eq.'1') then
c
c        find norm1(a).
c
         value = zero
         udiag = diag.eq.'u'
         if ( uplo.eq.'u' ) then
            do 220 j = 1, n
               if ((udiag) .and. (j.le.m)) then
                  sum = one
                  do 180 i = 1, j - 1
                     sum = sum + abs(a(i,j))
  180             continue
               else
                  sum = zero
                  do 200 i = 1, min(m,j)
                     sum = sum + abs(a(i,j))
  200             continue
               end if
               value = max(value,sum)
  220       continue
         else
            do 280 j = 1, n
               if (udiag) then
                  sum = one
                  do 240 i = j + 1, m
                     sum = sum + abs(a(i,j))
  240             continue
               else
                  sum = zero
                  do 260 i = j, m
                     sum = sum + abs(a(i,j))
  260             continue
               end if
               value = max(value,sum)
  280       continue
         end if
      else if ( norm.eq.'i' ) then

c        find normi(a).

         if ( uplo.eq.'u' ) then
            if ( diag.eq.'u' ) then
               do 300 i = 1, m
                  work(i) = one
  300          continue
               do 340 j = 1, n
                  do 320 i = 1, min(m,j-1)
                     work(i) = work(i) + abs(a(i,j))
  320             continue
  340          continue
            else
               do 360 i = 1, m
                  work(i) = zero
  360          continue
               do 400 j = 1, n
                  do 380 i = 1, min(m,j)
                     work(i) = work(i) + abs(a(i,j))
  380             continue
  400          continue
            end if
         else
            if ( diag.eq.'u' ) then
               do 420 i = 1, n
                  work(i) = one
  420          continue
               do 440 i = n + 1, m
                  work(i) = zero
  440          continue
               do 480 j = 1, n
                  do 460 i = j + 1, m
                     work(i) = work(i) + abs(a(i,j))
  460             continue
  480          continue
            else
               do 500 i = 1, m
                  work(i) = zero
  500          continue
               do 540 j = 1, n
                  do 520 i = j, m
                     work(i) = work(i) + abs(a(i,j))
  520             continue
  540          continue
            end if
         end if
         value = zero
         do 560 i = 1, m
            value = max(value,work(i))
  560    continue
      else if (norm.eq.'f' .or. norm.eq.'e') then

c        find normf(a).

         if ( uplo.eq.'u' ) then
            if ( diag.eq.'u' ) then
               scale = one
               sum = min(m,n)
               do 580 j = 2, n
                  call sssq (min(m,j-1),a(1,j),1,scale,sum)
  580          continue
            else
               scale = zero
               sum = one
               do 600 j = 1, n
                  call sssq (min(m,j),a(1,j),1,scale,sum)
  600          continue
            end if
         else
            if ( diag.eq.'u' ) then
               scale = one
               sum = min(m,n)
               do 620 j = 1, n
                  call sssq (m-j,a(min(m,j+1),j),1,scale,sum)
  620          continue
            else
               scale = zero
               sum = one
               do 640 j = 1, n
                  call sssq (m-j+1,a(j,j),1,scale,sum)
  640          continue
            end if
         end if
         value = scale*sqrt(sum)
      end if

      dlantr = value
c                                 end of dlantr
      end

      subroutine scsg ( t, c, s )

      double precision   c, s, t

c  scsg  returns values c and s such that
c
c     c = cos( theta ),   s = sin( theta )
c
c  for a given value of
c
c     t = tan( theta ).
c
c  c is always non-negative and s has the same sign as t, so that
c
c     c = 1.0/sqrt( 1.0 + t**2 ),   s = t/sqrt( 1.0 + t**2 ).

      double precision   one
      parameter        ( one = 1.0d+0 )

      double precision   abst, eps, reps, rrteps, rteps
      logical            first
c     .. external functions ..
      double precision   x02ajf
      external           x02ajf

      save               first, eps, reps, rteps, rrteps

      data               first/ .true. /

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      if( first )then
         first  = .false.
         eps    =  wmach(3)
         reps   =  1/eps
         rteps  =  sqrt( eps )
         rrteps =  1/rteps
      end if
c
      abst = abs( t )
      if( abst.lt.rteps )then
         c = one
         s = t
      else if( abst.gt.rrteps )then
         c = 1/abst
         s = sign( one, t )
      else
         c = 1/sqrt( 1 + abst**2 )
         s = c*t
      end if
c                                 end of scsg
      end

      subroutine sssq ( n, x, incx, scale, sumsq )

      double precision   scale, sumsq
      integer            incx, n

      double precision   x( * )
c
c  sssq  returns the values scl and smsq such that
c
c     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
c
c  where x( i ) = x( 1 + ( i - 1 )*incx ). the value of sumsq is assumed
c  to be at least unity and the value of smsq will then satisfy
c
c     1.0 .le. smsq .le. ( sumsq + n ) .
c
c  scale is assumed to be non-negative and scl returns the value
c
c     scl = max( scale, abs( x( i ) ) ) .
c
c  scale and sumsq must be supplied in scale and sumsq respectively.
c  scl and smsq are overwritten on scale and sumsq respectively.
c
c  the routine makes only one pass through the vector x.

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   absxi
      integer            ix
c----------------------------------------------------------------------
      if( n.gt.0 )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  sumsq = 1     + sumsq*( scale/absxi )**2
                  scale = absxi
               else
                  sumsq = sumsq +       ( absxi/scale )**2
               end if
            end if
   10    continue
      end if
c                                 end of sssq
      end

      subroutine scond ( n, x, incx, xmax, xmin )

      double precision   xmax, xmin
      integer            incx, n

      double precision   x( * )

c  scond  returns the values xmax and xmin given by

c     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
c             i                              i

c  if n is less than unity then xmax and xmin are returned as zero.

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      integer            ix
c----------------------------------------------------------------------
      if( n.lt.1 )then
         xmax = zero
         xmin = zero
      else
         xmax = abs( x( 1 ) )
         xmin = xmax
         do 10 ix = 1 + incx, 1 + ( n - 1 )*incx, incx
            xmax = max( xmax, abs( x( ix ) ) )
            xmin = min( xmin, abs( x( ix ) ) )
   10    continue
      end if
c                                 end of scond
      end

      subroutine iload ( n, const, x, incx )

      integer            const, incx, n

      integer            x( * )

c  iload  does the operation

c     x = const*e,   e' = ( 1  1 ... 1 ).

      integer            ix
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( const.ne.0 )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = const
   10       continue
         else
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = 0
   20       continue
         end if
      end if
c                                 end of iload
      end

      integer function isrank( n, x, incx, tol )

      double precision         tol
      integer                  incx, n

      double precision         x( * )

c  isrank finds the first element of the n element vector x for which
c
c     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
c
c  and returns the value ( k - 1 ) in the function name isrank. if no
c  such k exists then isrank is returned as n.
c
c  if tol is supplied as less than zero then the value epsmch, where
c  epsmch is the relative machine precision, is used in place of tol.

      double precision         zero
      parameter              ( zero = 0.0d+0 )

      double precision         tl, xmax
      integer                  ix, k
c     .. external functions ..
      double precision         x02ajf
      external                 x02ajf

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      k = 0
      if( n.ge.1 )then
         ix = 1
         if( tol.lt.zero )then
            tl = wmach(3)
         else
            tl = tol
         end if
         xmax = abs( x( ix ) )

c+       while( k.lt.n )loop
   10    if   ( k.lt.n )then
            if( abs( x( ix ) ).le.tl*xmax ) go to 20
            xmax = max( xmax, abs( x( ix ) ) )
            k    = k  + 1
            ix   = ix + incx
            go to 10
         end if
c+       end while

      end if

   20 isrank = k
c                                 end of isrank
      end

      subroutine sutsqr( side, n, k1, k2, c, s, a, lda )

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )

c  sutsqr does the transformation

c     r := p*u*q'  when  side = 'l' or 'l'  (  left-hand side )

c     r := q*u*p'  when  side = 'r' or 'r'  ( right-hand side ),

c  where  u and r  are  n by n  upper  triangular  matrices,   p  is  an
c  orthogonal matrix,  consisting of a given sequence of plane rotations
c  to be  applied  in  planes  k1 to k2,  and  q  is  a  unitary  matrix
c  consisting of a sequence of plane rotations, applied in planes  k1 to
c  k2,  chosen to make  r  upper triangular.

c  when  side = 'l' or 'l'  then  p  is  given  as a  sequence of  plane
c  rotation matrices
c
c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c  where  p( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
c  in this case the matrix q is given as
c
c     q = q( k2 - 1 )*...*q( k1 + 1 )*q( k1 ),
c
c  where  q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
c
c  when  side = 'r' or 'r'  then  p  is  given  as a  sequence of  plane
c  rotation matrices
c
c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c  where  p( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
c  in this case the matrix q is given as
c
c     q = q( k1 )*q( k1 + 1 )*...*q( k2 - 1 ),
c
c  where  q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
c
c  the  upper  triangular  matrix  u  must  be  supplied  in the  n by n
c  leading upper triangular part of  a,  and this  is overwritten by the
c  upper triangular matrix  r.  the cosine  and  sine  that  define  the
c  plane rotation matrix  p( k )  must be supplied in  c( k ) and s( k )
c  respectively,  and  the two by two rotation part of  p( k ),  t( k ),
c  is assumed to be of the form
c
c     t( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  the cosine  and  sine that define  q( k )  are overwritten on  c( k )
c  and  s( k )  respectively and the two by two rotation part of  q( k )
c  will have the form of  t( k )  above.
c
c  if  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
c  k2  is greater than  n  then an immediate return is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, fill, stemp, temp
      integer            i, i1, j
c     .. external subroutines ..
      external           srotgc
c     .. intrinsic functions ..
      intrinsic          min
c     ..
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return
      if( side.eq.'l' )then

c        apply the left-hand transformations,  column by column,  to the
c        triangular part of  u,  but not to  anywhere  that would  cause
c        fill.

         do 20 j = k1 + 1, n

c           apply  p( k1 ) ... p( j - 1 )  to column j.

            aij = a( k1, j )
            do 10 i = k1, min( j - 1, k2 - 1 )
               a( i, j ) = s( i )*a( i + 1, j ) + c( i )*aij
               aij = c( i )*a( i + 1, j ) - s( i )*aij
   10       continue
            a( i, j ) = aij
   20    continue

c           apply each  left-hand tranformation  to form the fill-in
c           elements and apply a  right-hand transformation to eliminate
c           the fill-in element.

         do 40 j = k1, k2 - 1

c           apply  p( j )  to the jth diagonal element  and the  fill-in
c           position.

            fill = -s( j )*a( j, j )
            a( j, j ) = c( j )*a( j, j )

c            set up  the rotation  q( j )  to eliminate the  fill-in
c           element,  and  apply  q( j )  to  the  jth  and  ( j + 1 )th
c           columns.

            call srotgc( a( j + 1, j + 1 ), fill, ctemp, stemp )
            c( j ) = ctemp
            s( j ) = -stemp
            if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
               stemp = -stemp
               do 30 i = 1, j
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
   30          continue
            end if
   40    continue
      else if ( side.eq.'r' ) then

c        we intermingle the  left and right hand transformations so that
c        at the kth step we form

c           a := q( k )*a*p( k )'.

c        first  apply  the  transformations  in  columns  k2 back to k1.

         do 60 j = k2 - 1, k1, -1

c           first apply  p( j ).

            if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
               ctemp = c( j )
               stemp = s( j )
               do 50 i = 1, j
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
   50          continue

c              next form the fill-in element  a( j + 1, j )  by applying
c              p( j ).

               fill = s( j )*a( j + 1, j + 1 )
               a( j + 1, j + 1 ) = c( j )*a( j + 1, j + 1 )

c              set up the rotation  q( j )  to eliminate the fill-in
c              element.

               call srotgc( a( j, j ), fill, c( j ), s( j ) )
            end if
   60    continue

c        finally  apply  q( k2 - 1 ) ... q( k1 )  to columns  n  back to
c        ( k1 + 1 ).

         do 80 j = n, k1 + 1, -1
            i1 = min( k2, j )
            aij = a( i1, j )
            do 70 i = i1 - 1, k1, -1
               temp = a( i, j )
               a( i + 1, j ) = c( i )*aij - s( i )*temp
               aij = s( i )*aij + c( i )*temp
   70       continue
            a( k1, j ) = aij
   80    continue
      end if
c                                 end of sutsqr
      end

      subroutine sdscl ( n, d, incd, x, incx )

      integer            incd, incx, n

      double precision   d( * ), x( * )

c  sdscl  does the operation
c
c     x := diag( d )*x

      integer            i, id, ix
c     .. external subroutines ..
      external           dscal

      if( n.gt.0 )then
         if( ( incd.eq.0 ).and.( incx.ne.0 ) )then
            call dscal( n, d( 1 ), x, abs( incx ) )
         else if( ( incd.eq.incx ).and.( incd.gt.0 ) )then
            do 10, id = 1, 1 + ( n - 1 )*incd, incd
               x( id ) = d( id )*x( id )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incd.gt.0 )then
               do 20, id = 1, 1 + ( n - 1 )*incd, incd
                  x( ix ) = d( id )*x( ix )
                  ix      = ix              + incx
   20          continue
            else
               id = 1 - ( n - 1 )*incd
               do 30, i = 1, n
                  x( ix ) = d( id )*x( ix )
                  id      = id              + incd
                  ix      = ix              + incx
   30          continue
            end if
         end if
      end if
c                                 end of sdscl
      end

      subroutine ssrotg( pivot, direct, n, alpha, x, incx, c, s )

      double precision   alpha
      integer            incx, n
      character*1        direct, pivot

      double precision   c( * ), s( * ), x( * )

c  ssrotg generates the parameters of an orthogonal matrix p such that

c     when   pivot = 'f' or 'f'   and   direct = 'f' or 'f'
c     or     pivot = 'v' or 'v'   and   direct = 'b' or 'b'
c
c        p*( alpha ) = ( beta ),
c          (   x   )   (   0  )
c
c     when   pivot = 'f' or 'f'   and   direct = 'b' or 'b'
c     or     pivot = 'v' or 'v'   and   direct = 'f' or 'f'
c
c        p*(   x   ) = (   0  ),
c          ( alpha ) = ( beta )
c
c  where alpha is a scalar and x is an n element vector.
c
c  when  pivot = 'f' or 'f'  ( fixed pivot )
c  and  direct = 'f' or 'f'  ( forward sequence ) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p( n )*p( n - 1 )*...*p( 1 )
c
c     where p( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'v' or 'v'  ( variable pivot )
c  and  direct = 'b' or 'b'  ( backward sequence ) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p( 1 )*p( 2 )*...*p( n )
c
c     where p( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'f' or 'f'  ( fixed pivot )
c  and  direct = 'b' or 'b'  ( backward sequence ) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p( 1 )*p( 2 )*...*p( n )
c
c     where p( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
c     designed to annihilate the kth element of x.
c
c  when  pivot = 'v' or 'v'  ( variable pivot )
c  and  direct = 'f' or 'f'  ( forward sequence ) then
c
c     p is given as the sequence of plane rotation matrices
c
c        p = p( n )*p( n - 1 )*...*p( 1 )
c
c     where p( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
c     designed to annihilate the kth element of x.
c
c  the routine returns the cosine, c( k ), and sine, s( k ) that define
c  the matrix p( k ), such that the two by two rotation part of p( k ),
c  r( k ), has the form
c
c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  on entry, alpha must contain  the scalar alpha and on exit, alpha is
c  overwritten by beta. the cosines and sines are returned in the arrays
c  c and s and the vector x is overwritten by the tangents of the plane
c  rotations ( t( k ) = s( k )/c( k ) ).

      integer            i, ix
c     .. external subroutines ..
      external           srotgc
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( direct.eq.'b' )then
            ix = 1 + ( n - 1 )*incx
            if( pivot.eq.'v' )then
               do 10, i = n, 2, -1
                  call srotgc( x( ix - incx ), x( ix ), c( i ), s( i ) )
                  ix = ix - incx
   10          continue
               call srotgc( alpha, x( ix ), c( 1 ), s( 1 ) )
            else if( pivot.eq.'f' )then

c              here we choose c and s so that

c                 ( alpha ) := (  c  s )*( alpha  )
c                 (   0   )    ( -s  c ) ( x( i ) )

c              which is equivalent to

c                 (   0   ) := ( c  -s )*( x( i ) )
c                 ( alpha )    ( s   c ) ( alpha  )
c
c              and so we need to return  s( i ) = -s  in order to make
c              r( i ) look like
c
c                 r( i ) = (  c( i )  s( i ) ).
c                          ( -s( i )  c( i ) )
c
               do 20, i = n, 1, -1
                  call srotgc( alpha, x( ix ), c( i ), s( i ) )
                  s( i )  = -s( i )
                  x( ix ) = -x( ix )
                  ix      =  ix      - incx
   20          continue
            end if
         else if( direct.eq.'f' )then
            ix = 1
            if( pivot.eq.'v'  )then

c              here we choose c and s so that
c
c                 ( x( i + 1 ) ) := (  c  s )*( x( i + 1 ) )
c                 (    0       )    ( -s  c ) ( x( i )     )
c
c              which is equivalent to
c
c                 (    0       ) := ( c  -s )*( x( i )     )
c                 ( x( i + 1 ) )    ( s   c ) ( x( i + 1 ) )
c
c              and so we need to return  s( i ) = -s  in order to make
c              r( i ) look like
c
c                 r( i ) = (  c( i )  s( i ) ).
c                          ( -s( i )  c( i ) )
c
               do 30, i = 1, n - 1
                  call srotgc( x( ix + incx ), x( ix ), c( i ), s( i ) )
                  s( i )  = -s( i )
                  x( ix ) = -x( ix )
                  ix      =  ix      + incx
   30          continue
               call srotgc( alpha, x( ix ), c( n ), s( n ) )
               s( n )  = -s( n )
               x( ix ) = -x( ix )
            else if( pivot.eq.'f' )then
               do 40, i = 1, n
                  call srotgc( alpha, x( ix ), c( i ), s( i ) )
                  ix = ix + incx
   40          continue
            end if
         end if
      end if
c                                 end of ssrotg
      end

      double precision function sdiv ( a, b, fail )

      double precision                  a, b
      logical                           fail

c  sdiv  returns the value div given by

c     div = ( a/b                 if a/b does not overflow,
c           (
c           ( 0.0                 if a .eq. 0.0,
c           (
c           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,

c  where  flmax  is a large value, via the function name. in addition if
c  a/b would overflow then  fail is returned as true, otherwise  fail is
c  returned as false.

c  note that when  a and b  are both zero, fail is returned as true, but
c  div  is returned as  0.0. in all other cases of overflow  div is such
c  that  abs( div ) = flmax.

c  when  b = 0  then  sign( a/b )  is taken as  sign( a ).

      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      double precision      absb, div, flmax, flmin
      logical               first

      double precision wmach
      common/ cstmch /wmach(10)

      save                  first, flmin, flmax

      data                  first/ .true. /
c----------------------------------------------------------------------
      if( a.eq.zero )then
         div = zero
         if( b.eq.zero )then
            fail = .true.
         else
            fail = .false.
         end if
      else
c
         if( first )then
            first = .false.
            flmin =  wmach(10)
            flmax =  1/flmin
         end if
c
         if( b.eq.zero )then
            div  =  sign( flmax, a )
            fail = .true.
         else
            absb = abs( b )
            if( absb.ge.one )then
               fail = .false.
               if( abs( a ).ge.absb*flmin )then
                  div = a/b
               else
                  div = zero
               end if
            else
               if( abs( a ).le.absb*flmax )then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if( ( ( a.lt.zero ).and.( b.gt.zero ) ).or.
     *                ( ( a.gt.zero ).and.( b.lt.zero ) )     )
     *               div = -div
               end if
            end if
         end if
      end if

      sdiv  = div
c                                 end of sdiv
      end

      subroutine susqr ( side, n, k1, k2, c, s, a, lda )

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )
c     ..
c
c  susqr  restores an upper spiked matrix  h to upper triangular form by
c  applying a sequence of plane rotations, in planes  k1 up to k2,  from
c  either the left, or the right.
c
c  the matrix  h is assumed to have non-zero elements only in the spiked
c  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
c  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
c  elements s( k ).
c
c  when  side = 'l' or 'l'  (  left-hand side )
c
c     h  is  assumed  to have a  row spike  and is restored to the upper
c     triangular matrix  r as
c
c        r = p*h,
c
c     where p is an orthogonal matrix of the form
c
c        p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c     p( k )  being a  plane rotation  matrix for the  ( k, k2 )  plane.
c
c  when  side = 'r' or 'r'  ( right-hand side )
c
c     h  is assumed to have a  column spike and is restored to the upper
c     triangular matrix r as
c
c        r = h*p',
c
c     where p is an orthogonal matrix of the form
c
c        p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c     p( k ) being a plane rotation matrix for the  ( k1, k + 1 ) plane.
c
c  the  two by two  rotation  part of  p( k ),  q( k ),  is of  the form
c
c     q( k ) = (  c( k )  s( k ) )
c              ( -s( k )  c( k ) )
c
c  and  c( k ) and s( k ) are returned in the kth elements of the arrays
c  c and s respectively.
c
c  the upper triangular part of the matrix  h must be supplied in the  n
c  by n  leading upper triangular part of  a, and this is overwritten by
c  the upper triangular matrix r.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, spike, stemp, temp
      integer            i, j
c     .. external subroutines ..
      external           srotgc
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return
      if( side.eq.'l' )then

c        restore h to upper triangular form by annihilating the elements
c        in  the  spike  of  h.  the  jth rotation  is  chosen  so  that

c        ( h( j, j ) ) := (  c  s )*( h( j , j ) ).
c        (     0     )    ( -s  c ) ( h( k2, j ) )

c        apply the rotations in columns k1 up to ( k2 - 1 ).

         do 20 j = k1, k2 - 1
            spike = s( j )
            do 10 i = k1, j - 1
               aij = a( i, j )
               a( i, j ) = s( i )*spike + c( i )*aij
               spike = c( i )*spike - s( i )*aij
   10       continue

c           set up the rotation.

            call srotgc( a( j, j ), spike, c( j ), s( j ) )
   20    continue

c        apply the rotations to columns k2 up to n.

         do 40 j = k2, n
            temp = a( k2, j )
            do 30 i = k1, k2 - 1
               aij = a( i, j )
               a( i, j ) = s( i )*temp + c( i )*aij
               temp = c( i )*temp - s( i )*aij
   30       continue
            a( k2, j ) = temp
   40    continue
      else if( side.eq.'r' )then

c        restore h to upper triangular form by annihilating the spike of
c        h. the jth rotation is chosen so that

c           ( h( j, j ) ) := (  c  s )*( h( j, j )  ),
c           (     0     )    ( -s  c ) ( h( j, k1 ) )

c        which can be expressed as

c           ( 0  h( j, j ) ) := ( h( j, k1 )  h( j, j ) )*(  c  s ).
c                                                         ( -s  c )

c        thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
c        rotation matrix look like

c           q( j ) = (  c( j )  s( j ) ).
c                    ( -s( j )  c( j ) )

         do 70 j = k2, k1 + 1, -1
            call srotgc( a( j, j ), s( j - 1 ), ctemp, stemp )
            stemp = -stemp
            s( j - 1 ) = stemp
            c( j - 1 ) = ctemp
            if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
               do 50 i = j - 1, k1 + 1, -1
                  spike = s( i - 1 )
                  s( i - 1 ) = stemp*a( i, j ) + ctemp*spike
                  a( i, j ) = ctemp*a( i, j ) - stemp*spike
   50          continue
               do 60 i = k1, 1, -1
                  temp = a( i, k1 )
                  a( i, k1 ) = stemp*a( i, j ) + ctemp*temp
                  a( i, j ) = ctemp*a( i, j ) - stemp*temp
   60          continue
            end if
   70    continue
      end if

c                                 end of susqr

      end

      subroutine sgesrc( side, pivot, direct, m, n, k1, k2, c, s, a,
     *                   lda )

      integer            k1, k2, lda, m, n
      character*1        direct, pivot, side

      double precision   a( lda, * ), c( * ), s( * )

c  sgesrc  does the transformation
c
c     a := p*a,   when   side = 'l' or 'l'  (  left-hand side )
c
c     a := a*p',  when   side = 'r' or 'r'  ( right-hand side )
c
c  where a is an m by n matrix and p is an orthogonal matrix, consisting
c  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
c  determined by the parameters pivot and direct as follows:
c
c     when  pivot  = 'v' or 'v'  ( variable pivot )
c     and   direct = 'f' or 'f'  ( forward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c        where  p( k )  is a plane rotation matrix for the  ( k, k + 1 )
c        plane.
c
c     when  pivot  = 'v' or 'v'  ( variable pivot )
c     and   direct = 'b' or 'b'  ( backward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c        where  p( k )  is a plane rotation matrix for the  ( k, k + 1 )
c        plane.
c
c     when  pivot  = 't' or 't'  ( top pivot )
c     and   direct = 'f' or 'f'  ( forward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k2 - 1 )*p( k2 - 2 )*...*p( k1 ),
c
c        where  p( k )  is a plane rotation matrix for the ( k1, k + 1 )
c        plane.
c
c     when  pivot  = 't' or 't'  ( top pivot )
c     and   direct = 'b' or 'b'  ( backward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c        where  p( k )  is a plane rotation matrix for the ( k1, k + 1 )
c        plane.
c
c     when  pivot  = 'b' or 'b'  ( bottom pivot )
c     and   direct = 'f' or 'f'  ( forward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k2 - 1 )*p( k2 - 2 )*...*p( k1 ),
c
c        where  p( k )  is a  plane rotation  matrix  for the  ( k, k2 )
c        plane.
c
c     when  pivot  = 'b' or 'b'  ( bottom pivot )
c     and   direct = 'b' or 'b'  ( backward sequence ) then
c
c        p is given as a sequence of plane rotation matrices
c
c           p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c        where  p( k )  is a  plane rotation  matrix  for the  ( k, k2 )
c        plane.
c
c  c( k ) and s( k )  must contain the  cosine and sine  that define the
c  matrix  p( k ).  the  two by two  plane rotation  part of the  matrix
c  p( k ), r( k ), is assumed to be of the form
c
c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  if m, n or k1 are less than unity,  or k2 is not greater than k1,  or
c  side = 'l' or 'l'  and  k2  is greater than  m, or  side = 'r' or 'r'
c  and  k2  is greater than  n,  then an  immediate return  is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, stemp, temp
      integer            i, j
      logical            left, right
c----------------------------------------------------------------------
      left = ( side.eq.'l' )
      right = ( side.eq.'r' )
      if( ( min( m, n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *    ( ( left ).and.( k2.gt.m ) ).or.
     *    ( ( right ).and.( k2.gt.n ) ) )return
      if( left )then
         if( pivot.eq.'v' )then
            if( direct.eq.'f'  )then
               do 20 j = 1, n
                  aij = a( k1, j )
                  do 10 i = k1, k2 - 1
                     temp = a( i + 1, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     aij = c( i )*temp - s( i )*aij
   10             continue
                  a( k2, j ) = aij
   20          continue
            else if( direct.eq.'b' )then
               do 40 j = 1, n
                  aij = a( k2, j )
                  do 30 i = k2 - 1, k1, -1
                     temp = a( i, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     aij = s( i )*aij + c( i )*temp
   30             continue
                  a( k1, j ) = aij
   40          continue
            end if
         else if( pivot.eq.'t' )then
            if( direct.eq.'f' )then
               do 60 j = 1, n
                  temp = a( k1, j )
                  do 50 i = k1, k2 - 1
                     aij = a( i + 1, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     temp = s( i )*aij + c( i )*temp
   50             continue
                  a( k1, j ) = temp
   60          continue
            else if( direct.eq.'b' )then
               do 80 j = 1, n
                  temp = a( k1, j )
                  do 70 i = k2 - 1, k1, -1
                     aij = a( i + 1, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     temp = s( i )*aij + c( i )*temp
   70             continue
                  a( k1, j ) = temp
   80          continue
            end if
         else if( pivot.eq.'b' )then
            if( direct.eq.'f' )then
               do 100 j = 1, n
                  temp = a( k2, j )
                  do 90 i = k1, k2 - 1
                     aij = a( i, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     temp = c( i )*temp - s( i )*aij
   90             continue
                  a( k2, j ) = temp
  100          continue
            else if( direct.eq.'b' )then
               do 120 j = 1, n
                  temp = a( k2, j )
                  do 110 i = k2 - 1, k1, -1
                     aij = a( i, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     temp = c( i )*temp - s( i )*aij
  110             continue
                  a( k2, j ) = temp
  120          continue
            end if
         end if
      else if( right )then
         if( pivot.eq.'v' )then
            if( direct.eq.'f' )then
               do 140 j = k1, k2 - 1
                  if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
                     ctemp = c( j )
                     stemp = s( j )
                     do 130 i = 1, m
                        temp = a( i, j + 1 )
                        a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                        a( i, j ) = stemp*temp + ctemp*a( i, j )
  130                continue
                  end if
  140          continue
            else if( direct.eq.'b' )then
               do 160 j = k2 - 1, k1, -1
                  if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
                     ctemp = c( j )
                     stemp = s( j )
                     do 150 i = m, 1, -1
                        temp = a( i, j + 1 )
                        a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                        a( i, j ) = stemp*temp + ctemp*a( i, j )
  150                continue
                  end if
  160          continue
            end if
         else if( pivot.eq.'t' )then
            if( direct.eq.'f' )then
               do 180 j = k1 + 1, k2
                  ctemp = c( j - 1 )
                  stemp = s( j - 1 )
                  if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
                     do 170 i = 1, m
                        temp = a( i, j )
                        a( i, j ) = ctemp*temp - stemp*a( i, k1 )
                        a( i, k1 ) = stemp*temp + ctemp*a( i, k1 )
  170                continue
                  end if
  180          continue
            else if( direct.eq.'b' )then
               do 200 j = k2, k1 + 1, -1
                  ctemp = c( j - 1 )
                  stemp = s( j - 1 )
                  if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
                     do 190 i = m, 1, -1
                        temp = a( i, j )
                        a( i, j ) = ctemp*temp - stemp*a( i, k1 )
                        a( i, k1 ) = stemp*temp + ctemp*a( i, k1 )
  190                continue
                  end if
  200          continue
            end if
         else if( pivot.eq.'b' )then
            if( direct.eq.'f' )then
               do 220 j = k1, k2 - 1
                  if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
                     ctemp = c( j )
                     stemp = s( j )
                     do 210 i = 1, m
                        temp = a( i, j )
                        a( i, j ) = stemp*a( i, k2 ) + ctemp*temp
                        a( i, k2 ) = ctemp*a( i, k2 ) - stemp*temp
  210                continue
                  end if
  220          continue
            else if( direct.eq.'b' )then
               do 240 j = k2 - 1, k1, -1
                  if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
                     ctemp = c( j )
                     stemp = s( j )
                     do 230 i = m, 1, -1
                        temp = a( i, j )
                        a( i, j ) = stemp*a( i, k2 ) + ctemp*temp
                        a( i, k2 ) = ctemp*a( i, k2 ) - stemp*temp
  230                continue
                  end if
  240          continue
            end if
         end if
      end if
c                                 end of sgesrc
      end

      subroutine icopy ( n, x, incx, y, incy )

      integer            incx, incy, n
      integer            x( * ), y( * )

c  icopy  does the operation

c     y := x

      integer            i, ix, iy
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( iy )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  y( iy ) = x( ix )
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  y( iy ) = x( ix )
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if

c                                 end of icopy

      end

      double precision function snorm( scale, ssq )

      double precision                  scale, ssq

c  snorm returns the value norm given by
c
c     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
c            (
c            ( flmax,             scale*sqrt( ssq ) .ge. flmax
c
c  via the function name.

      double precision      flmax, flmin, norm, sqt
      logical               first

      double precision wmach
      common/ cstmch /wmach(10)

      save                  first, flmax

      data                  first/ .true. /
c----------------------------------------------------------------------
      if( first )then
         first = .false.
         flmin =  wmach(10)
         flmax =  1/flmin
      end if
c
      sqt = sqrt( ssq )
      if( scale.lt.flmax/sqt )then
         norm = scale*sqt
      else
         norm = flmax
      end if
c
      snorm = norm

c                                 end of snorm.

      end

      subroutine sutsrh( side, n, k1, k2, c, s, a, lda )

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )

c  sutsrh applies a  given sequence  of  plane rotations  to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the rotations are applied
c  in planes k1 up to k2.
c
c  the upper hessenberg matrix, h, is formed as

c     h = p*u,    when   side = 'l' or 'l',  (  left-hand side )

c  where p is an orthogonal matrix of the form

c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 )

c  and is formed as

c     h = u*p',   when   side = 'r' or 'r',  ( right-hand side )

c  where p is an orthogonal matrix of the form

c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),

c  p( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. the
c  cosine and sine that define p( k ), k = k1, k1 + 1, ..., k2 - 1, must
c  be  supplied  in  c( k )  and  s( k )  respectively.  the  two by two
c  rotation part of p( k ), r( k ), is assumed to have the form

c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of  h.

c  the  sub-diagonal elements of  h, h( k + 1, k ),  are returned in the
c  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.

c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, stemp, temp
      integer            i, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return
      if( side.eq.'l' )then

c        apply the plane rotations to columns n back to k1.

         do 20 j = n, k1, -1
            if( j.ge.k2 )then
               aij = a( k2, j )
            else

c              form  the  additional sub-diagonal element  h( j + 1, j )
c              and store it in s( j ).

               aij = c( j )*a( j, j )
               s( j ) = -s( j )*a( j, j )
            end if
            do 10 i = min( k2, j ) - 1, k1, -1
               temp = a( i, j )
               a( i + 1, j ) = c( i )*aij - s( i )*temp
               aij = s( i )*aij + c( i )*temp
   10       continue
            a( k1, j ) = aij
   20    continue
      else if( side.eq.'r' )then

c        apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
c        and  form   the   additional  sub-diagonal  elements,   storing
c        h( j + 1, j ) in s( j ).

         do 40 j = k1, k2 - 1
            if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
               stemp = s( j )
               ctemp = c( j )
               do 30 i = 1, j
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
   30          continue
               s( j ) = stemp*a( j + 1, j + 1 )
               a( j + 1, j + 1 ) = ctemp*a( j + 1, j + 1 )
            end if
   40    continue
      end if
c                                 end of sutsrh.
      end

      subroutine sgeapr( side, trans, n, perm, k, b, ldb )

      integer            k, ldb, n
      character*1        side, trans

      double precision   perm( * ), b( ldb, * )

c  sgeapr does one of the transformations

c     b := p'*b   or   b := p*b,   where b is an m by k matrix,

c  or

c     b := b*p'   or   b := b*p,   where b is a k by m matrix,

c  p being an m by m permutation matrix of the form

c     p = p( 1, index( 1 ) )*p( 2, index( 2 ) )*...*p( n, index( n ) ),

c  where  p( i, index( i ) ) is the permutation matrix that interchanges
c  items i and index( i ). that is p( i, index( i ) ) is the unit matrix
c  with rows and columns  i and  index( i )  interchanged. of course, if
c  index( i ) = i  then  p( i, index( i ) ) = i.

c  parameters
c  ==========

c  side   - character*1.
c  trans
c           on entry,  side  ( left-hand side, or right-hand side )  and
c           trans  ( transpose, or no transpose )  specify the operation
c           to be doed as follows.

c           side = 'l' or 'l'   and   trans = 't' or 't'

c              do the operation   b := p'*b.

c           side = 'l' or 'l'   and   trans = 'n' or 'n'

c              do the operation   b := p*b.

c           side = 'r' or 'r'   and   trans = 't' or 't'

c              do the operation   b := b*p'.
c
c           side = 'r' or 'r'   and   trans = 'n' or 'n'

c              do the operation   b := b*p.

c           unchanged on exit.

c  n      - integer.

c           on entry, n must specify the value of n.  n must be at least
c           zero.  when  n = 0  then an  immediate  return  is effected.

c           unchanged on exit.

c  perm   - real             array of dimension at least ( n ).

c           before  entry,  perm  must  contain  the  n indices  for the
c           permutation matrices. index( i ) must satisfy

c              1 .le. index( i ) .le. m.

c           it is usual for index( i ) to be at least i, but this is not
c           necessary for this routine. it is assumed that the statement
c           index = perm( i )  returns the correct integer in  index, so
c           that,  if necessary,  perm( i )  should contain a real value
c           slightly larger than  index.

c           unchanged on exit.

c  k      - integer.

c           on entry with  side = 'l' or 'l',  k must specify the number
c           of columns of b and on entry with  side = 'r' or 'r', k must
c           specify the number of rows of  b.  k must be at least  zero.
c           when  k = 0  then an immediate return is effected.

c           unchanged on exit.

c  b      - real  array  of  dimension ( ldb, ncolb ),  where  ncolb = k
c           when  side = 'l' or 'l'  and  ncolb = m  when  side = 'r' or
c           'r'.

c           before entry  with  side = 'l' or 'l',  the  leading  m by k
c           part  of  the  array   b  must  contain  the  matrix  to  be
c           transformed  and before  entry with  side = 'r' or 'r',  the
c           leading  k by m part of the array  b must contain the matrix
c           to  be  transformed.  on exit,   b  is  overwritten  by  the
c           transformed matrix.

c  ldb    - integer.

c           on entry,  ldb  must specify  the  leading dimension  of the
c           array  b  as declared  in the  calling  (sub) program.  when
c           side = 'l' or 'l'   then  ldb  must  be  at  least  m,  when
c           side = 'r' or 'r'   then  ldb  must  be  at  least  k.
c           unchanged on exit.

      logical            left, null, right, trnsp
      integer            i, j, l
      double precision   temp
c----------------------------------------------------------------------
      if( min( n, k ).eq.0 ) return
      left = ( side.eq.'l' )
      right = ( side.eq.'r' )
      null = ( trans.eq.'n' )
      trnsp = ( trans.eq.'t' )
      if( left )then
         if( trnsp )then
            do 20 i = 1, n
               l = perm( i )
               if( l.ne.i )then
                  do 10 j = 1, k
                     temp = b( i, j )
                     b( i, j ) = b( l, j )
                     b( l, j ) = temp
   10             continue
               end if
   20       continue
         else if( null )then
            do 40 i = n, 1, -1
               l = perm( i )
               if( l.ne.i )then
                  do 30 j = 1, k
                     temp = b( l, j )
                     b( l, j ) = b( i, j )
                     b( i, j ) = temp
   30             continue
               end if
   40       continue
         end if
      else if( right )then
         if( trnsp )then
            do 60 j = n, 1, -1
               l = perm( j )
               if( l.ne.j )then
                  do 50 i = 1, k
                     temp = b( i, j )
                     b( i, j ) = b( i, l )
                     b( i, l ) = temp
   50             continue
               end if
   60       continue
         else if( null )then
            do 80 j = 1, n
               l = perm( j )
               if( l.ne.j )then
                  do 70 i = 1, k
                     temp = b( i, l )
                     b( i, l ) = b( i, j )
                     b( i, j ) = temp
   70             continue
               end if
   80       continue
         end if
      end if
c                                 end of sgeapr.
      end

      subroutine sutsr1 (side,n,k1,k2,s,a,lda)

c  sutsr1 applies a  sequence  of  pairwise interchanges to either  the
c  left,  or the right,  of the  n by n  upper triangular matrix  u,  to
c  transform u to an  upper hessenberg matrix. the interchanges are
c  applied in planes k1 up to k2.
c
c  the upper hessenberg matrix, h, is formed as
c
c     h = p*u,    when   side = 'l' or 'l',  (  left-hand side )
c
c  where p is a permutation matrix of the form
c
c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 )
c
c  and is formed as
c
c     h = u*p',   when   side = 'r' or 'r',  ( right-hand side )
c
c  where p is a permutation matrix of the form
c
c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c  p( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
c  the  two by two
c  interchange part of p( k ), r( k ), is assumed to have the form
c
c     r( k ) = ( 0  1 ).
c              ( 1  0 )
c
c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of  h.
c
c  the  sub-diagonal elements of  h, h( k + 1, k ),  are returned in the
c  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision  zero
      parameter         (zero=0.0d+0)

      integer           k1, k2, lda, n
      character*1       side

      double precision  a(lda,*), s(*)

      double precision  aij, temp
      integer           i, j
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return
      if (side.eq.'l') then

c        apply the permutations to columns n back to k1.

         do 40 j = n, k1, -1
            if (j.ge.k2) then
               aij = a(k2,j)
            else

c              form  the  additional sub-diagonal element  h( j + 1, j )
c              and store it in s( j ).

               aij = zero
               s(j) = a(j,j)
            end if
            do 20 i = min(k2,j) - 1, k1, -1
               temp = a(i,j)
               a(i+1,j) = temp
               aij = aij
   20       continue
            a(k1,j) = aij
   40    continue
      else if (side.eq.'r') then

c        apply  the  plane interchanges to  columns  k1  up to
c        ( k2 - 1 ) and  form   the   additional  sub-diagonal
c        elements,   storing  h( j + 1, j ) in s( j ).

         do 80 j = k1, k2 - 1
            do 60 i = 1, j
               temp = a(i,j+1)
               a(i,j+1) = a(i,j)
               a(i,j) = temp
   60       continue
            s(j) = a(j+1,j+1)
            a(j+1,j+1) = zero
   80    continue
      end if
c                                 end of sutsrh.
      end

      subroutine sgrfg( n, alpha, x, incx, tol, zeta )
c----------------------------------------------------------------------
      double precision   alpha, tol, zeta
      integer            incx, n

      double precision   x( * )

c  sgrfg generates details of a generalized householder reflection such
c  that

c     p*( alpha ) = ( beta ),   p'*p = i.
c       (   x   )   (   0  )

c  p is given in the form

c     p = i - ( zeta )*( zeta  z' ),
c             (   z  )

c  where z is an n element vector and zeta is a scalar that satisfies

c     1.0 .le. zeta .le. sqrt( 2.0 ).

c  zeta is returned in zeta unless x is such that

c     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )

c  where eps is the relative machine precision and tol is the user
c  supplied value tol, in which case zeta is returned as 0.0 and p can
c  be taken to be the unit matrix.

c  beta is overwritten on alpha and z is overwritten on x.
c  the routine may be called with  n = 0  and advantage is taken of the
c  case where  n = 1.

      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   beta, eps, scale, ssq
      logical            first
c     .. external functions ..
      double precision   x02ajf
      external           x02ajf
c     .. external subroutines ..
      external           sssq , dscal

      double precision wmach
      common/ cstmch /wmach(10)

      save               eps, first

      data               first/ .true. /
c----------------------------------------------------------------------
      if( n.lt.1 )then
         zeta = zero
      else if( ( n.eq.1 ).and.( x( 1 ).eq.zero ) )then
         zeta = zero
      else
c
         if( first )then
            first = .false.
            eps   =  wmach(3)
         end if
c
c        treat case where p is a 2 by 2 matrix specially.
c
         if( n.eq.1 )then
c
c           deal with cases where  alpha = zero  and
c           abs( x( 1 ) ) .le. max( eps*abs( alpha ), tol )  first.
c
            if( alpha.eq.zero )then
               zeta   =  one
               alpha  =  abs ( x( 1 ) )
               x( 1 ) = -sign( one, x( 1 ) )
            else if( abs( x( 1 ) ).le.max( eps*abs( alpha ), tol ) )then
               zeta   =  zero
            else
               if( abs( alpha ).ge.abs( x( 1 ) ) )then
                  beta = abs( alpha ) *sqrt( 1 + ( x( 1 )/alpha )**2 )
               else
                  beta = abs( x( 1 ) )*sqrt( 1 + ( alpha/x( 1 ) )**2 )
               end if
               zeta = sqrt( ( abs( alpha ) + beta )/beta )
               if( alpha.ge.zero ) beta = -beta
               x( 1 ) = -x( 1 )/( zeta*beta )
               alpha  = beta
            end if
         else

c           p is larger than 2 by 2.

            ssq   = one
            scale = zero
            call sssq ( n, x, incx, scale, ssq )

c           treat cases where  scale = zero,
c           scale .le. max( eps*abs( alpha ), tol )  and
c           alpha = zero  specially.
c           note that  scale = max( abs( x( i ) ) ).

            if( ( scale.eq.zero ).or.
     *          ( scale.le.max( eps*abs( alpha ), tol ) ) )then
               zeta  = zero
            else if( alpha.eq.zero )then
               zeta  = one
               alpha = scale*sqrt( ssq )
               call dscal( n, -1/alpha, x, incx )
            else
               if( scale.lt.abs( alpha ) )then
                  beta = abs( alpha )*sqrt( 1 + ssq*( scale/alpha )**2 )
               else
                  beta = scale       *sqrt( ssq +   ( alpha/scale )**2 )
               end if
               zeta = sqrt( ( beta + abs( alpha ) )/beta )
               if( alpha.gt.zero ) beta = -beta
               call dscal( n, -1/( zeta*beta ), x, incx )
               alpha = beta
            end if
         end if
      end if
c                                 end of sgrfg
      end

      subroutine nggqzz(hess,n,k1,k2,c,s,a,lda)

c  nggqzz  either applies a  given sequence  of  plane rotations  to the
c  right of the n by n reverse lower triangular matrix t, to transform t
c  to a  reverse lower hessenberg matrix  h, or restores a reverse lower
c  hessenberg matrix h to reverse lower triangular form t, by applying a
c  sequence of plane rotations from the right.
c
c  the rotations are applied  in planes k1 up to k2.
c
c  when   hess = 'c' or 'c',   ( create ),  then   the   reverse   lower
c  hessenberg matrix, h, is formed as
c
c     h = t*p',
c
c  where p is an orthogonal matrix of the form
c
c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c  p( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. the
c  cosine and sine that define p( k ), k = k1, k1 + 1, ..., k2 - 1, must
c  be  supplied  in  c( k )  and  s( k )  respectively.  the  two by two
c  rotation part of p( k ), r( k ), is assumed to have the form
c
c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  the matrix  t must be supplied in the n by n reverse lower triangular
c  part  of the array  a,  and this is overwritten by the  reverse lower
c  triangular part of  h.
c
c  the super-diagonal elements of  h, h( n - k, k ), are returned in the
c  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c  when   hess = 'r' or 'r',   ( remove ),  then   the   reverse   lower
c  hessenberg matrix  h  is  assumed  to  have  non-zero  super-diagonal
c  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
c  only and  h( n - k, k ) must be supplied in  s( k ). h is restored to
c  the reverse lower triangular matrix t as

c     t = h*p',

c  where p is an orthogonal matrix of the form

c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),

c  p( k ) being a plane rotation for the  ( k, k + 1 ) plane. the cosine
c  and  sine  that  define  p( k )  are  returned  in  c( k ) and s( k )
c  respectively.  the  two by two  rotation part of  p( k ),  r( k ), is
c  of the form

c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  the reverse lower triangular part of the matrix h must be supplied in
c  the  n by n  reverse  lower  triangular  part  of  a,   and  this  is
c  overwritten by the reverse triangular matrix t.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c  when   n = 7, k1 = 2 and k2 = 5   then  t  and  h  are  of  the  form
c
c     t = ( 0  0  0  0  0  0  x ),   h = ( 0  0  0  0  0  0  x ).
c         ( 0  0  0  0  0  x  x )        ( 0  0  0  0  x  x  x )
c         ( 0  0  0  0  x  x  x )        ( 0  0  0  x  x  x  x )
c         ( 0  0  0  x  x  x  x )        ( 0  0  x  x  x  x  x )
c         ( 0  0  x  x  x  x  x )        ( 0  x  x  x  x  x  x )
c         ( 0  x  x  x  x  x  x )        ( 0  x  x  x  x  x  x )
c         ( x  x  x  x  x  x  x )        ( x  x  x  x  x  x  x )
c

      double precision  one, zero
      parameter         (one=1.0d+0,zero=0.0d+0)

      integer           k1, k2, lda, n
      character*1       hess

      double precision  a(lda,*), c(*), s(*)

      double precision  ctemp, stemp, suph, temp
      integer           i, j
c     .. external subroutines ..
      external          srotgc
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return
      if (hess.eq.'c') then

c        apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
c        and  form   the  additional  super-diagonal  elements,  storing
c        h( n - j, j ) in s( j ).

         do 40 j = k1, k2 - 1
            if ((c(j).ne.one) .or. (s(j).ne.zero)) then
               stemp = s(j)
               ctemp = c(j)
               s(j) = stemp*a(n-j,j+1)
               a(n-j,j+1) = ctemp*a(n-j,j+1)
               do 20 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   20          continue
            end if
   40    continue
      else if (hess.eq.'r') then

c        restore  h to reverse lower triangular form by annihilating the
c        super-diagonal elements of  h.  the  jth rotation  is chosen so
c        that

c          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
c          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )

c        which can be expressed as

c           ( 0  h( n - j, n - j ) ) :=

c               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
c                                                            ( -s  c )

c        thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
c        rotation matrix look like

c           r( j ) = (  c( j )  s( j ) ).
c                    ( -s( j )  c( j ) )

         do 80 j = k2 - 1, k1, -1
            suph = s(j)
            call srotgc(a(n-j,j+1),suph,ctemp,stemp)
            stemp = -stemp
            s(j) = stemp
            c(j) = ctemp
            if ((ctemp.ne.one) .or. (stemp.ne.zero)) then
               do 60 i = n - j + 1, n
                  temp = a(i,j+1)
                  a(i,j+1) = ctemp*temp - stemp*a(i,j)
                  a(i,j) = stemp*temp + ctemp*a(i,j)
   60          continue
            end if
   80    continue
      end if
c                                 end of nggqzz.
      end

      subroutine sutsrs( side, n, k1, k2, c, s, a, lda )

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )

c  sutsrs applies a  given sequence  of  plane rotations  to either  the
c  left,  or the right,  of the  n by n  upper triangular  matrix  u  to
c  transform  u  to an upper spiked matrix. the rotations are applied in
c  planes k1 up to k2.

c  the upper spiked matrix, h, is formed as

c     h = p*u,   when   side = 'l' or 'l',  ( left-hand side )

c  where p is an orthogonal matrix of the form

c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),

c  p( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
c  formed as

c     h = u*p',   when   side = 'r' or 'r',  ( right-hand side )

c  where p is an orthogonal matrix of the form

c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),

c  p( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.

c  the cosine and sine that define  p( k ), k = k1, k1 + 1, ..., k2 - 1,
c  must be  supplied  in  c( k ) and s( k ) respectively. the two by two
c  rotation part of p( k ), r( k ), is assumed to have the form

c     r( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )

c  the matrix  u must be supplied in the n by n leading upper triangular
c  part of the array  a, and this is overwritten by the upper triangular
c  part of h.

c  when  side = 'l' or 'l'  then a  row spike  is  generated  in  h  and
c  when  side = 'r' or 'r'  then a  column spike is generated. for a row
c  spike the elements  h( k2, k )  and for a  column spike  the elements
c  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).

c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, spike, stemp, temp
      integer            i, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return
      if( side.eq.'l' )then

c        apply the plane rotations to columns n back to k2.

         do 20 j = n, k2, -1
            temp = a( k2, j )
            do 10 i = k2 - 1, k1, -1
               aij = a( i, j )
               a( i, j ) = s( i )*temp + c( i )*aij
               temp = c( i )*temp - s( i )*aij
   10       continue
            a( k2, j ) = temp
   20    continue

c        form  the spike  and apply the rotations in columns  ( k2 - 1 )
c        back to k1.

         do 40 j = k2 - 1, k1, -1
            spike = -s( j )*a( j, j )
            a( j, j ) = c( j )*a( j, j )
            do 30 i = j - 1, k1, -1
               aij = a( i, j )
               a( i, j ) = s( i )*spike + c( i )*aij
               spike = c( i )*spike - s( i )*aij
   30       continue
            s( j ) = spike
   40    continue
      else if( side.eq.'r' )then

c        apply the  plane rotations to columns  ( k1 + 1 ) up to k2  and
c        form the spike.

         do 70 j = k1 + 1, k2
            ctemp = c( j - 1 )
            stemp = s( j - 1 )
            if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
               do 50 i = 1, k1
                  temp = a( i, k1 )
                  a( i, k1 ) = stemp*a( i, j ) + ctemp*temp
                  a( i, j ) = ctemp*a( i, j ) - stemp*temp
   50          continue
               do 60 i = k1 + 1, j - 1
                  spike = s( i - 1 )
                  s( i - 1 ) = stemp*a( i, j ) + ctemp*spike
                  a( i, j ) = ctemp*a( i, j ) - stemp*spike
   60          continue
               s( j - 1 ) = stemp*a( j, j )
               a( j, j ) = ctemp*a( j, j )
            end if
   70    continue
      end if
c                                 end of sutsrs
      end

      subroutine smload( matrix, m, n, const, diag, a, lda )

      character*1        matrix
      double precision   const, diag
      integer            lda, m, n

      double precision   a( lda, * )

c  smload forms the m by n matrix a given by
c
c     a( i, j ) = (  diag  i.eq.j,
c                 (
c                 ( const  i.ne.j.
c
c  if   matrix = 'g' or 'g'   then  a  is regarded  as a general matrix,
c  if   matrix = 'u' or 'u'   then  a  is regarded  as upper triangular,
c                             and only  elements  for which  i.le.j  are
c                             referenced,
c  if   matrix = 'l' or 'l'   then  a  is regarded  as lower triangular,
c                             and only  elements  for which  i.ge.j  are
c                             referenced.

      integer            i, j

      if( matrix.eq.'g' )then
         do 20 j = 1, n
            do 10 i = 1, m
               a( i, j ) = const
   10       continue
   20    continue
         if( const.ne.diag )then
            do 30 i = 1, min( m, n )
               a( i, i ) = diag
   30       continue
         end if
      else if( matrix.eq.'u' )then
         do 50 j = 1, n
            do 40 i = 1, min( m, j )
               a( i, j ) = const
   40       continue
   50    continue
         if( const.ne.diag )then
            do 60 i = 1, min( m, n )
               a( i, i ) = diag
   60       continue
         end if
      else if( matrix.eq.'l' )then
         do 80 j = 1, min( m, n )
            do 70 i = j, m
               a( i, j ) = const
   70       continue
   80    continue
         if( const.ne.diag )then
            do 90 i = 1, min( m, n )
               a( i, i ) = diag
   90       continue
         end if
      end if
c                                 end of smload
      end

      subroutine suhqr( side, n, k1, k2, c, s, a, lda )

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )

c  suhqr restores an upper hessenberg matrix h to upper triangular form
c  by  applying a sequence of  plane rotations  from either the left, or
c  the right.  the matrix  h  is assumed to have  non-zero  sub-diagonal
c  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
c  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).
c
c  h is restored to the upper triangular matrix r either as
c
c     r = p*h,   when   side = 'l' or 'l'  (  left-hand side )
c
c  where p is an orthogonal matrix of the form
c
c     p = p( k2 - 1 )*...*p( k1 + 1 )*p( k1 ),
c
c  or as
c
c     r = h*p',  when   side = 'r' or 'r'  ( right-hand side )
c
c  where p is an orthogonal matrix of the form
c
c     p = p( k1 )*p( k1 + 1 )*...*p( k2 - 1 ),
c
c  in both cases  p( k )  being a  plane rotation  for the  ( k, k + 1 )
c  plane.  the cosine and sine that define p( k ) are returned in c( k )
c  and  s( k )  respectively.  the two by two  rotation part of  p( k ),
c  q( k ), is of the form
c
c     q( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  the upper triangular part of the matrix  h  must be supplied in the n
c  by n  leading upper triangular part of  a, and this is overwritten by
c  the upper triangular matrix r.
c
c  if n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, stemp, subh, temp
      integer            i, j
c     .. external subroutines ..
      external           srotgc
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return
      if( side.eq.'l' )then
c
c        restore   h  to  upper  triangular  form  by  annihilating  the
c        sub-diagonal elements of h.  the jth rotation is chosen so that
c
c           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
c           (     0     )    ( -s  c ) ( h( j + 1, j ) )
c
c        apply the rotations in columns k1 up to n.
c
         do 20 j = k1, n
            aij = a( k1, j )
            do 10 i = k1, min( j, k2 ) - 1
               temp = a( i + 1, j )
               a( i, j ) = s( i )*temp + c( i )*aij
               aij = c( i )*temp - s( i )*aij
   10       continue
            if( j.lt.k2 )then

c              set up the rotation.

               subh = s( j )
               call srotgc( aij, subh, c( j ), s( j ) )
               a( j, j ) = aij
            else
               a( k2, j ) = aij
            end if
   20    continue
      else if( side.eq.'r' )then

c        restore   h  to  upper  triangular  form  by  annihilating  the
c        sub-diagonal elements of h.  the jth rotation is chosen so that

c           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
c           (         0         )    ( -s  c ) ( h( j + 1, j )     )

c        which can be expressed as

c           ( 0  h( j + 1, j + 1 ) ) :=

c               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
c                                                    ( -s  c )

c        thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
c        rotation matrix look like

c           q( j ) = (  c( j )  s( j ) ).
c                    ( -s( j )  c( j ) )

         do 40 j = k2 - 1, k1, -1
            subh = s( j )
            call srotgc( a( j + 1, j + 1 ), subh, ctemp, stemp )
            stemp = -stemp
            s( j ) = stemp
            c( j ) = ctemp
            if( ( ctemp.ne.one ).or.( stemp.ne.zero ) )then
               do 30 i = j, 1, -1
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
   30          continue
            end if
   40    continue
      end if
c                                 end of suhqr.
      end

      subroutine smcopy( matrix, m, n, a, lda, b, ldb )

      character*1        matrix
      integer            m, n, lda, ldb

      double precision   a( lda, * ), b( ldb, * )

c  smcopy  copies  the  m by n  matrix  a  into  the  m by n  matrix  b.
c
c  if   matrix = 'g' or 'g'   then  a  and  b  are  regarded as  general
c                             matrices,
c  if   matrix = 'u' or 'u'   then  a  and  b  are  regarded  as   upper
c                             triangular,  and only  elements  for which
c                             i.le.j  are referenced,
c  if   matrix = 'l' or 'l'   then  a  and  b  are  regarded  as   lower
c                             triangular,  and only  elements  for which
c                             i.ge.j  are referenced.
      integer            i, j
c----------------------------------------------------------------------
      if( matrix.eq.'g' )then
         do 20 j = 1, n
            do 10 i = 1, m
               b( i, j ) = a( i, j )
   10       continue
   20    continue
      else if( matrix.eq.'u' )then
         do 40 j = 1, n
            do 30 i = 1, min( m, j )
               b( i, j ) = a( i, j )
   30       continue
   40    continue
      else if( matrix.eq.'l' )then
         do 60 j = 1, min( m, n )
            do 50 i = j, m
               b( i, j ) = a( i, j )
   50       continue
   60    continue
      end if
c                                 end of smcopy
      end

      subroutine srotgc( a, b, c, s )
c----------------------------------------------------------------------
      double precision   a, b, c, s

c  srotg blas except that c is always
c  returned as non-negative and  b  is overwritten by the tangent of the
c  angle that defines the plane rotation.
c
c  c and s are given as
c
c     c = 1.0/sqrt( 1.0 + t**2 ),   s = c*t   where   t = b/a.
c
c  when  abs( b ) .le. eps*abs( a ),  where  eps is the relative machine
c  precision as  returned by routine  x02ajf,  then  c and s  are always
c  returned as
c
c     c = 1.0  and  s = 0.0
c
c  and when  abs( a ) .le. eps*abs( b ) then c and s are always returned
c  as
c
c     c = 0.0  and  s = sign( t ).
c
c  note that t is always returned as  b/a, unless this would overflow in
c  which  case the value  sign( t )*flmax  is returned,  where  flmax is
c  the value given by  1/x02amf( ).
c
c  c and s  can be reconstructed from the tangent,  t,  by a call to the
c  basic linear algebra routine scsg .

      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   t
      logical            fail
c     .. external functions ..
      double precision   sdiv 
      external           sdiv 
c     .. external subroutines ..
      external           scsg 
c----------------------------------------------------------------------
      if( b.eq.zero )then
         c  = one
         s  = zero
      else
         t  = sdiv ( b, a, fail )
         call scsg ( t, c, s )
         a  = c*a + s*b
         b  = t
      end if
c                                 end of srotgc
      end

      subroutine sload( n, const, x, incx )
c----------------------------------------------------------------------

      double precision   const
      integer            incx, n

      double precision   x( * )
c  sload does the operation

c     x = const*e,   e' = ( 1  1 ... 1 ).
      double precision   zero
      parameter        ( zero = 0.0d+0 )

      integer            ix
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( const.ne.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = const
   10       continue
         else
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   20       continue
         end if
      end if
c                                 end of sload
      end
