      subroutine lpsol (n,nclin,a,lda,bl,bu,cvec,istate,x,iter,obj,ax,
     *                  clamda,iw,leniw,w,lenw,ifail)
c----------------------------------------------------------------------
c lpsol is the linear programming algorithm of 

c Gill P E, Murray W, Saunders M A and Wright M H (1984) Procedures for 
c optimization problems with a mixture of bounds and general linear 
c constraints ACM Trans. Math. Software 10 282–298

c     a  is a constant  nclin by n  matrix.
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
c-----------------------------------------------------------------------
      implicit none

      character*6       srname
      parameter         (srname='lpsol')
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
c     .. scalar arguments ..
      double precision  obj
      integer           ifail, iter, lda, leniw, lenw, n, nclin
c     .. array arguments ..
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
      logical           header, lcdbg, newopt, prnt
c     .. arrays in common ..
      double precision  rpadlc(23), rpsvlc(mxparm)
      integer           ilcdbg(ldbg), ipadlc(14), ipsvlc(mxparm),
     *                  loclc(lenlc), nfix(2)
c     .. local scalars ..
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
c     .. local arrays ..
      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*8       names(1)
      character*80      errrec(2), rec(2)
c     .. external functions ..
      double precision  dnrm2
      integer           p01abf
      external          dnrm2, p01abf
c                                 external subroutine
      external          e04mfu
c     .. common blocks ..
      common            /ae04mf/loclc
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04mf/newopt
      common            /be04nb/lennam, ldt, ncolt, ldq
      common            /ce04mf/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04mf/alfa, trulam, isdel, jdel, jadd, header,
     *                  prnt
      common            /de04nb/asize, dtmax, dtmin
      common            /ee04mf/ilcdbg, lcdbg
      common            /fe04mf/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /ge04mf/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)
c     .. save statement ..
      save              /ae04mf/, /be04mf/, /fe04mf/, /ge04mf/
c     .. data statements ..
      data              title/' *** lpsol'/
c----------------------------------------------------------------------
      epsmch = wmach(3)
      rteps = wmach(4)
      nout = 6
      nerr = 6

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      named = .false.

      msglvl = 0
      inform = 0
      iter = 0

      header = .true.
      prnt = .true.

      condmx = max(one/epspt5,hundrd)

c     set the default values of the parameters.

      call e04mfw (n,nclin,title)

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

c     assign the dimensions of arrays in the parameter list of e04mfz.
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
      else if (lcrash.eq.2) then
         start = 'hot '
      end if

      cold = lcrash .eq. 0
      warm = lcrash .eq. 1
      hot = lcrash .eq. 2

c     allocate remaining work arrays.

      litotl = 3
      lwtotl = 0
      call e04mfv(cset,n,nclin,litotl,lwtotl)

c     check input parameters and storage limits.

      call e04mfp(nerror,msglvl,start,leniw,lenw,litotl,lwtotl,n,nclin,
     *            ncnln,istate,named,names,bigbnd,bl,bu,x,lda,ldh,mm,
     *            nerr,llptyp,ifail)

      if (nerror.gt.0) then
         msg = 'errors'
         go to 60
      end if

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

      if (tolfea.gt.zero) call sload (n+nclin,(tolfea),w(lfeatu),1)
c
      call e04mfr('initialize anti-cycling variables',msglvl,n,nclin,
     *            nmoved,iter,numinf,istate,bigbnd,ax,bl,bu,clamda,
     *            w(lfeatu),x)

      if (cold .or. warm) then
c        cold or warm start.  just about everything must be initialized.
c        the only exception is istate during a warm start.

         ianrmj = lanorm

         do j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
         end do

         if (nclin.gt.0) call f06flf(nclin,w(lanorm),1,asize,amin)

         call f06flf(nctotl,w(lfeatu),1,feamax,feamin)
         call dcopy (nctotl,w(lfeatu),1,w(lwtinf),1)
         call dscal (nctotl,(one/feamin),w(lwtinf),1)

c        define the initial working set.
c               nfree ,  nactiv,  kactiv, kx,
c               istate (if start  = 'cold')
c               nartif (if vertex = 'true')

         call e04mft(start,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
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

            call e04nfq(unitq,vertex,1,nact1,it,nactiv,nartif,nz,nfree,
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

         call dcopy (n,cvec,1,w(lcq),1)
         call e04nbw(6,n,nz,nfree,ldq,unitq,iw(lkx),w(lcq),w(lq),w(lwrk)
     *               )
      end if

      rset = .false.
      itmax = itmax2
      jinf = 0

c     +    take your pick when minimizing the sum of infeasibilities:
c     +    nrz    =  nz  implies steepest-descent in the two-norm.
c     +    nrz    =  0   implies steepest-descent in the infinity norm.
      nrz = 0

c     move x onto the constraints in the working set.

   40 call e04mfj(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,ldt,
     *            istate,iw(lkactv),iw(lkx),jmax,errmax,xnorm,a,ax,bl,
     *            bu,w(lfeatu),w(lt),x,w(lq),w(ld),w(lwrk))

      if (rowerr) then
         if (msglvl.gt.0) then
            write (rec,fmt=99981)
            call x04baf(iprint,rec(1))
         end if
         msg = 'infeas'
         numinf = 1
         obj = errmax
         go to 60
      end if

      call e04mfz(prbtyp,msg,cset,named,names,rset,unitq,iter,itmax,
     *            jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,nz,istate,
     *            iw(lkactv),iw(lkx),e04mfu,obj,numinf,xnorm,a,ax,bl,bu,
     *            cvec,clamda,w(lfeatu),x,iw,w)

      found = msg .eq. 'feasbl' .or. msg .eq. 'optiml' .or. msg .eq.
     *        'weak  ' .or. msg .eq. 'unbndd' .or. msg .eq. 'infeas'
      halted = msg .eq. 'itnlim'

      if (found) then
         call e04mfr('optimal',msglvl,n,nclin,nmoved,iter,numinf,istate,
     *               bigbnd,ax,bl,bu,clamda,w(lfeatu),x)
      end if

      done = found .and. nviol .eq. 0 .and. nmoved .eq. 0

      if ( .not. (done .or. halted)) go to 40

c     set   clamda.  print the full solution.
c     clean up.  save values for a subsequent hot start.

      call e04mfk(msglvl,nfree,lda,n,nclin,ncnln,nctotl,bigbnd,named,
     *            names,lennam,nactiv,istate,iw(lkactv),iw(lkx),a,bl,bu,
     *            x,clamda,w(lrlam),x)

      iw(1) = 0
      if (unitq) iw(1) = 1
      iw(2) = nfree
      iw(3) = nactiv

c     print messages if required.
c     recover the optional parameters set by the user.

   60 if (msg.eq.'optiml') then
         inform = 0
         if (msglvl.gt.0) then
            write (rec,fmt=99997) prbtyp
            call x04bay(iprint,2,rec)
         end if
c
      else if (msg.eq.'feasbl') then
         inform = 0
         if (msglvl.gt.0) then
            write (rec,fmt=99998)
            call x04bay(iprint,2,rec)
         end if
c
      else if (msg.eq.'weak  ') then
         inform = 1
         if (msglvl.gt.0) then
            write (rec,fmt=99996) prbtyp
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99995)
     *       prbtyp
c
      else if (msg.eq.'unbndd') then
         inform = 2
         if (msglvl.gt.0) then
            write (rec,fmt=99994) prbtyp
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99993)
     *       prbtyp
c
      else if (msg.eq.'infeas') then
         inform = 3
         if (msglvl.gt.0) then
            write (rec,fmt=99992)
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99991)
c
      else if (msg.eq.'itnlim') then
         inform = 4
         if (msglvl.gt.0) then
            write (rec,fmt=99990)
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99989)
c
      else if (msg.eq.'errors') then
         inform = 6
         if (msglvl.gt.0) then
            write (rec,fmt=99988) nerror
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99987)
     *       nerror
c
      else if (msg.eq.'noprob') then
         inform = 7
         if (msglvl.gt.0) then
            write (rec,fmt=99986)
            call x04bay(iprint,2,rec)
         end if
         if (ifail.eq.0 .or. ifail.eq.-1) write (errrec,fmt=99985)
      end if
c
      if (msglvl.gt.0) then
c
         if (inform.lt.6) then
            if (numinf.eq.0) then
               if (prbtyp.ne.'fp') then
                  write (rec,fmt=99984) prbtyp, obj
                  call x04bay(iprint,2,rec)
               end if
            else if (inform.eq.3) then
               write (rec,fmt=99983) obj
               call x04bay(iprint,2,rec)
            else
               write (rec,fmt=99982) obj
               call x04bay(iprint,2,rec)
            end if
         end if
      end if

      if (inform.lt.6) then
         if (msglvl.gt.0) then
            write (rec,fmt=99999) prbtyp, iter
            call x04bay(iprint,2,rec)
         end if
      end if

      call f06dff(mxparm,ipsvlc,1,iprmlc,1)
      call dcopy (mxparm,rpsvlc,1,rprmlc,1)

      if ((inform.ge.1 .and. inform.le.7)
     *    .and. (ifail.eq.0 .or. ifail.eq.-1)) call x04bay(nerr,2,
     *    errrec)
      ifail = p01abf(ifail,inform,srname,0,rec)

      if (ifail.lt.4) ifail = 0

99999 format (/' exit from ',a2,' problem after ',i5,' iterations.')
99998 format (/' exit lpsol - feasible point found.     ')
99997 format (/' exit lpsol - optimal ',a2,' solution.')
99996 format (/' exit lpsol - weak ',a2,' solution.')
99995 format (/' ** weak ',a2,' solution.')
99994 format (/' exit lpsol - ',a2,' solution is unbounded.')
99993 format (/' ** ',a2,' solution is unbounded.')
99992 format (/' exit lpsol - no feasible point for the linear constr',
     *       'aints.')
99991 format (/' ** no feasible point for the linear constraints.')
99990 format (/' exit lpsol - too many iterations.')
99989 format (/' ** too many iterations.')
99988 format (/' exit lpsol - ',i7,' errors found in the input parame',
     *       'ters.  problem abandoned.')
99987 format (/' ** ',i7,' errors found in the input parameters.  prob',
     *       'lem abandoned.')
99986 format (/' exit lpsol - problem type not recognized.  problem a',
     *       'bandoned.')
99985 format (/' ** problem type not recognized.  problem abandoned.')
99984 format (/' final ',a2,' objective value =',g16.7)
99983 format (/' minimum sum of infeasibilities =',g16.7)
99982 format (/' final sum of infeasibilities =',g16.7)
99981 format (' xxx  cannot satisfy the working set constraints to the',
     *       ' accuracy requested.')
      end

      subroutine nlpsol (n,nclin,ncnln,lda,ldcju,ldr,a,bl,bu,confun,
     *                  objfun,iter,istate,c,cjacu,clamda,objf,gradu,r,
     *                  x,iw,leniw,w,lenw,iuser,user,ifail,jprint)
c-----------------------------------------------------------------------
c nlpsol is the non-linear programming algorithm of 

c Gill P E, Murray W, Saunders M A and Wright M H (1984) Procedures for 
c optimization problems with a mixture of bounds and general linear 
c constraints ACM Trans. Math. Software 10 282–298

c     it solves the nonlinear program
c
c            minimize                   f(x)
c
c                                    (    x  )
c            subject to    bl  .le.  (  a*x  )  .le.  bu
c                                    (  c(x) )
c
c     where  f(x)  is a smooth scalar function,  a  is a constant matrix
c     and  c(x)  is a vector of smooth nonlinear functions.  the
c     feasible region is defined by a mixture of linear and nonlinear
c     equality or inequality constraints on  x.

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

c     nlpsol uses a sequential quadratic programming algorithm, with a
c     positive-definite quasi-newton approximation to the transformed
c     hessian  q'hq  of the lagrangian function (stored in r).

c     this material is based upon work partially supported by the
c     national science foundation under grants mcs-7926009 and
c     ecs-8312142; the department of energy contract am03-76sf00326,
c     pa no. de-at03-76er72018; the army research office contract
c     daa29-84-k-0156; and the office of naval research grant
c     n00014-75-c-0267.
c-----------------------------------------------------------------------
      implicit none

      character*6       srname
      parameter         (srname='nlpsol')
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
c     .. scalar arguments ..
      double precision  objf
      integer           ifail, iter, lda, ldcju, ldr, leniw, lenw, n,
     *                  nclin, ncnln, jprint
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  amin, cond, condmx, ctx, epsmch, errmax, fdchk,
     *                  fdnorm, feamax, feamin, obj, rootn, rteps, ssq1,
     *                  suminf, xnorm
      integer           i, ianrmj, idbg, ikx, info, inform,
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
c     .. local arrays ..
      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(2)
c     .. external functions ..
      double precision  dnrm2, adivb, f06rjf
      integer           p01abf
      external          dnrm2, adivb, f06rjf, p01abf
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common            /ae04uc/locnp

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04nb/lennam, ldt, ncolt, ldq
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
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c     .. save statement ..
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/,
     *                  /fe04nc/
c     .. data statements ..
      data              title/' *** nlpsol'/
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

c     default names will be provided for variables during printing.

      named = .false.
      inform = 0
c                                 print flags
      msgnp = jprint
      msgqp = jprint - 10

c     set the default values for the parameters.

      call e04ucx(n,nclin,ncnln,title)

      needfd = lvlder .eq. 0 .or. lvlder .eq. 2 .or.
     *         (lvlder.eq.1 .and. ncnln.gt.0)
      cold = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1

      nplin = n + nclin
      nctotl = nplin + ncnln

c     assign the dimensions of arrays in the parameter list of e04ucz.
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

c     e04ucp  defines the arrays that contain the locations of various
c     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      call e04ucp(n,nclin,ncnln,nctotl,litotl,lwtotl)

c     allocate certain addresses that are not allocated in e04ucp.

      lax = lwtotl + 1
      lwtotl = lax + nclin - 1
      lax = min(lax,lwtotl)

c     check input parameters and storage limits.

      call e04nbz(nerror,msgnp,lcrash,leniw,lenw,litotl,lwtotl,n,nclin,
     *            ncnln,istate,iw,named,names,bigbnd,bl,bu,x,m,lda,ldr,
     *            ldcju,ldfju,nerr,ifail)

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

      if (tolfea.gt.zero) call sload (nplin,tolfea,w(lfeatl),1)
c
      if (ncnln.gt.0 .and. ctol.gt.zero) call sload (ncnln,ctol,
     *    w(lfeatl+nplin),1)
c
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

         call e04ucy(info,msgnp,nstate,lvlder,nfun,ngrad,ldcj,ldcju,n,
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

      call f06dff(ldbg,ilsdbg,1,icmdbg,1)

      if (nclin.gt.0) then
         ianrmj = lanorm
         do 20 j = 1, nclin
            w(ianrmj) = dnrm2(n,a(j,1),lda)
            ianrmj = ianrmj + 1
   20    continue
         call f06flf(nclin,w(lanorm),1,asize,amin)
      end if

      call f06flf(nplin,w(lfeatl),1,feamax,feamin)
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
            call f06qhf('upper-triangular',n,n,zero,one,r,ldr)
            rfrobn = rootn

            nrank = 0
            if (ncnln.gt.0) call sload (ncnln,(zero),w(lcmul),1)
         else

c           r will be updated while finding a feasible x.

            nrank = nlnx
            call sload (nlnx,(zero),w(lres0),1)
            if (ncnln.gt.0) call dcopy (ncnln,clamda(nplin+1),1,w(lcmul)
     *                                 ,1)

         end if

         incrun = .true.
         rhonrm = zero
         rhodmp = one
         scale = one
         call sload (ncnln,(zero),w(lrho),1)

c        re-order kx so that the free variables come first.
c        if a warm start is required, nrank will be nonzero and the
c        factor r will be updated.

         call e04ncx(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldq,lda,ldr,
     *               ldt,istate,iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)

      end if

c     factorize the linear constraints in the initial working set.

      if (nactiv.gt.0) then
         nact1 = nactiv
         nactiv = 0

         call e04ncy(unitq,vertex,inform,1,nact1,nactiv,nartif,nz,nfree,
     *               nrank,nrejtd,nres,ngq,n,ldq,lda,ldr,ldt,istate,
     *               iw(lkactv),iw(lkx),condmx,a,r,w(lt),w(lres0),w(lgq)
     *               ,w(lq),w(lwrk1),w(lwrk2),w(lrlam),msgnp)
      end if

      if (lcrash.le.1) then

c        cold or warm start.  move  x  on to the linear constraints and
c        find a feasible point.

         ssq1 = zero
         linobj = .false.
         call e04nch(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,n,
     *               nplin,ldq,lda,ldr,ldt,istate,iw(lkactv),iw(lkx),
     *               jmax,errmax,ctx,xnorm,a,w(lax),bl,bu,w(lgq),w(lres)
     *               ,w(lres0),w(lfeatl),r,w(lt),x,w(lq),w(lwrk1),
     *               w(lwrk2))

c        call  e04ncz  to find a feasible  x.
c        use  work2  as the multiplier vector.

         jinf = 0
         lclam = lwrk2

         itmxsv = itmax1
         itmax1 = nminor

         call e04ncz('fp problem',named,names,linobj,unitq,nlperr,itns,
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

      end if

      if (lcrash.gt.0) then

c        check for a bad r.

         rfrobn = f06rjf('frobenius norm','upper','non-unit diagonal',n,
     *            n,r,ldr,w)
         call f06flf(n,r,ldr+1,drmax,drmin)
         cond = adivb(drmax,drmin,overfl)

         if (cond.gt.rcndbd .or. rfrobn.gt.rootn*growth*drmax) then

c           refactorize the hessian and bound the condition estimator.

            if (msgnp.gt.0) then
               write (rec,fmt=99986)
               call x04baf(iprint,rec(1))
            end if
            call e04udr(unitq,n,nfree,nz,ldq,ldr,iw(liperm),iw(lkx),
     *                  w(lgq),r,w(lq),w(lwrk1),w(lres0))
         end if
      end if

c     check the gradients at a feasible x.

      lvrfyc = lverfy
      if (lverfy.ge.10) lvrfyc = -1

      call e04ucy(info,msgnp,nstate,lvlder,nfun,ngrad,ldcj,ldcju,n,
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

      call dcopy (n,w(lgrad),1,w(lgq),1)
      call e04nbw(6,n,nz,nfree,ldq,unitq,iw(lkx),w(lgq),w(lq),w(lwrk1))

c     solve the problem.

      iuser(2) = 1

      if (ncnln.eq.0) then

c        the problem has only linear constraints and bounds.

         call e04ucz(named,names,unitq,inform,iter,n,nclin,ncnln,nctotl,
     *               nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,nfun,ngrad,
     *               istate,iw(lkactv),iw(lkx),objf,fdnorm,xnorm,confun,
     *               objfun,a,w(lax),bl,bu,c,w(lcjac),cjacu,clamda,
     *               w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw,iuser,user)
      else

c        the problem has some nonlinear constraints.

         if (nclin.gt.0) call f06qff('g',nclin,n,a,lda,w(laqp),
     *                               ldaqp)

c        try to add some nonlinear constraint indices to kactiv.

         call e04ucs(cold,n,nclin,ncnln,nctotl,nactiv,nfree,nz,istate,
     *               iw(lkactv),bigbnd,tolact,bl,bu,c)

         call e04ucz(named,names,unitq,inform,iter,n,nclin,ncnln,nctotl,
     *               nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,nfun,ngrad,
     *               istate,iw(lkactv),iw(lkx),objf,fdnorm,xnorm,confun,
     *               objfun,w(laqp),w(lax),bl,bu,c,w(lcjac),cjacu,
     *               clamda,w(lfeatl),w(lgrad),gradu,r,x,iw,w,lenw,
     *               iuser,user)

      end if

c     if required, form the triangular factor of the hessian.

c     first,  form the square matrix  r  such that  h = r'r.
c     compute the  qr  factorization of  r.

      if (lformh.gt.0) then
         lv = lwrk2
         do 60 j = 1, n
            if (j.gt.1) call sload (j-1,zero,w(lv),1)

            lvj = lv + j - 1
            call dcopy (n-j+1,r(j,j),ldr,w(lvj),1)
            call e04nbw(3,n,nz,nfree,ldq,unitq,iw(lkx),w(lv),w(lq),
     *                  w(lwrk1))
            call dcopy (n,w(lv),1,r(j,1),ldr)
   60    continue

         call f01qcf(n,n,r,ldr,w(lwrk1),info)
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

         if (inform.ge.0 .and. inform.lt.7) then
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

      call f06dff(mxparm,ipsvls,1,iprmls,1)
      call dcopy (mxparm,rpsvls,1,rprmls,1)
      call f06dff(mxparm,ipsvnp,1,iprmnp,1)
      call dcopy (mxparm,rpsvnp,1,rprmnp,1)

      if (inform.lt.9) then
         if (ncnln.gt.0) call f06qff('g',ncnln,n,w(lcjac),ldcj,
     *                               cjacu,ldcju)
         call dcopy (n,w(lgrad),1,gradu,1)
      end if

      ifail = p01abf(ifail,inform,srname,0,rec)

99999 format (/' exit nlpsol - user requested termination.')
99998 format (/' exit nlpsol - optimal solution found.')
99997 format (/' exit nlpsol - optimal solution found, but requested a',
     *       'ccuracy not achieved.')
99996 format (/' exit nlpsol - no feasible point for the linear constr',
     *       'aints.')
99995 format (/' exit nlpsol - no feasible point for the nonlinear con',
     *       'straints.')
99994 format (/' exit nlpsol - too many major iterations.             ')
99993 format (/' exit nlpsol - current point cannot be improved upon. ')
99992 format (/' exit nlpsol - large errors found in the derivatives. ')
99991 format (/' exit nlpsol - ',i7,' errors found in the input parame',
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

      subroutine e04mfp(nerror,msglvl,start,liwork,lwork,litotl,lwtotl,
     *                  n,nclin,ncnln,istate,named,names,bigbnd,bl,bu,x,
     *                  lda,ldh,mm,nerr,lprob,ifail)

c     e04mfp   checks the input data for lpsol and e04nff.
c----------------------------------------------------------------------
      implicit none

      double precision  bigbnd
      integer           ifail, lda, ldh, litotl, liwork, lprob, lwork,
     *                  lwtotl, mm, msglvl, n, nclin, ncnln, nerr,
     *                  nerror
      logical           named
      character*4       start
c     .. array arguments ..
      double precision  bl(n+nclin+ncnln), bu(n+nclin+ncnln), x(n)
      integer           istate(n+nclin+ncnln)
      character*8       names(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      double precision  b1, b2
      integer           is, j, k, l
      logical           ok
      character*4       cstart
c     .. local arrays ..
      character*5       id(3)
      character*80      rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
c     .. data statements ..
      data              id(1), id(2), id(3)/'varbl', 'l con', 'n con'/
c----------------------------------------------------------------------
      nerror = 0

      if (n.le.0) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99997) n
            call x04bay(nerr,3,rec)
         end if
      end if

      if (nclin.lt.0 .or. ncnln.lt.0) then
         if (nclin.lt.0) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99996) nclin
               call x04bay(nerr,3,rec)
            end if
         end if

         if (ncnln.lt.0) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99995) ncnln
               call x04bay(nerr,3,rec)
            end if
         end if
      end if

      if (lda.lt.max(1,nclin)) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99994) lda, nclin
            call x04bay(nerr,3,rec)
         end if
      end if

      if ((lprob.eq.1) .or. (lprob.eq.2)) then
c        problem type = fp or lp
         if (ldh.lt.1) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99993) ldh
               call x04bay(nerr,3,rec)
            end if
         end if
      else
c        problem type = qp1, qp2, qp3 or qp4
         if (ldh.lt.mm) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99992) ldh, n, mm
               call x04bay(nerr,4,rec)
            end if
         end if
      end if

c     check if there is enough workspace to solve the problem.

      ok = litotl .le. liwork .and. lwtotl .le. lwork
      if ( .not. ok) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
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

         do 20 j = 1, n + nclin + ncnln
            b1 = bl(j)
            b2 = bu(j)
            ok = b1 .lt. b2 .or. (b1.eq.b2 .and. abs(b1).lt.bigbnd)
            if ( .not. ok) then
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
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99991) names(j), j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99990) names(j), j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               else
                  if (b1.eq.b2) then
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99989) id(l), k, j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99988) id(l), k, j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               end if
            end if
   20    continue
c

c        check  istate.

         if (start.eq.'warm' .or. start.eq.'hot ') then
            if (start.eq.'warm') cstart = 'warm'
            if (start.eq.'hot ') cstart = 'hot '
            do 40 j = 1, n + nclin + ncnln
               is = istate(j)
               ok = is .ge. (-2) .and. is .le. 4
               if ( .not. ok) then
                  nerror = nerror + 1
                  if (ifail.eq.0 .or. ifail.eq.-1) then
                     write (rec,fmt=99987) cstart, j, j, is
                     call x04bay(nerr,3,rec)
                  end if
               end if
   40       continue
         end if
      end if

c     end of  e04mfp.  (cminit)

99999 format (/' workspace provided is     iwork(',i8,'),  work(',i8,
     *       ').',/' to solve problem we need  iwork(',i8,'),  work(',
     *       i8,').')
99998 format (/' ** not enough workspace to solve problem.')
99997 format (/' ** on entry, n.le.0:',/'    n = ',i16)
99996 format (/' ** on entry, nclin.lt.0:',/'    nclin = ',i16)
99995 format (/' ** on entry, ncnln.lt.0:',/'    ncnln = ',i16)
99994 format (/' ** on entry, lda.lt.max(1,nclin):',/'    lda = ',i16,
     *       '   nclin = ',i16)
99993 format (/' ** on entry, ldh.lt.1:',/'    ldh = ',i16)
99992 format (/' ** on entry, either ldh.lt.n or ldh.lt.m',/'    (wher',
     *       'e m is the value of the optional parameter hessian rows):'
     *       ,/'    ldh = ',i16,'   n = ',i16,'   m = ',i16)
99991 format (/' ** on entry, the equal bounds on  ',a8,'  are infinit',
     *       'e (because',/'    bl(',i4,').eq.beta and bu(',i4,').eq.b',
     *       'eta, but abs(beta).ge.bigbnd):',/'    beta =',g16.7,' bi',
     *       'gbnd =',g16.7)
99990 format (/' ** on entry, the bounds on  ',a8,'  are inconsistent:',
     *       /'    bl(',i4,') =',g16.7,'   bu(',i4,') =',g16.7)
99989 format (/' ** on entry, the equal bounds on  ',a5,i3,'  are infi',
     *       'nite (because',/'    bl(',i4,').eq.beta and bu(',i4,').e',
     *       'q.beta, but abs(beta).ge.bigbnd):',/'    beta =',g16.7,
     *       ' bigbnd =',g16.7)
99988 format (/' ** on entry, the bounds on  ',a5,i3,'  are inconsiste',
     *       'nt:',/'    bl(',i4,') =',g16.7,'   bu(',i4,') =',g16.7)
99987 format (/' ** on entry with a ',a4,' start, istate(',i4,') is ou',
     *       't of range:',/'    istate(',i4,') = ',i16)
      end

      subroutine lsmove(hitcon,hitlow,linobj,unitgz,nclin,nrank,nrz,n,
     *                  ldr,jadd,numinf,alfa,ctp,ctx,xnorm,ap,ax,bl,bu,
     *                  gq,hz,p,res,r,x,work)

c     lsmove  changes x to x + alfa*p and updates ctx, ax, res and gq
c     accordingly.
c
c     if a bound was added to the working set,  move x exactly on to it,
c     except when a negative step was taken (e04ucg may have had to move
c     to some other closer constraint.)
c-----------------------------------------------------------------------
      implicit none

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  alfa, ctp, ctx, xnorm
      integer           jadd, ldr, n, nclin, nrank, nrz, numinf
      logical           hitcon, hitlow, linobj, unitgz
c     .. array arguments ..
      double precision  ap(*), ax(*), bl(*), bu(*), gq(*), hz(*), p(n),
     *                  r(ldr,*), res(*), work(*), x(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  bnd
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
      call daxpy (n,alfa,p,1,x,1)
      if (linobj) ctx = ctx + alfa*ctp
c
      if (hitcon .and. jadd.le.n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa.ge.zero) x(jadd) = bnd
      end if
      xnorm = dnrm2(n,x,1)
c
      if (nclin.gt.0) call daxpy (nclin,alfa,ap,1,ax,1)

      if (nrz.le.nrank) then
         if (unitgz) then
            res(nrz) = res(nrz) - alfa*hz(nrz)
         else
            call daxpy (nrz,(-alfa),hz,1,res,1)
         end if

         if (numinf.eq.0) then

c           update the transformed gradient gq so that
c           gq = gq + alfa*r'( hz ).
c                            ( 0  )

            if (unitgz) then
               call daxpy (n-nrz+1,alfa*hz(nrz),r(nrz,nrz),ldr,
     *                     gq(nrz),1)
            else
               call dcopy (nrz,hz,1,work,1)
               call dtrmv ('u','t','n',nrz,r,ldr,work,1)
               if (nrz.lt.n) call dgemv ('t',nrz,n-nrz,one,r(1,nrz+1),
     *                                  ldr,hz,1,zero,work(nrz+1),1)

               call daxpy (n,alfa,work,1,gq,1)

            end if
         end if
      end if

      end

      subroutine chcore(debug,done,first,epsa,epsr,fx,inform,iter,itmax,
     *                  cdest,fdest,sdest,errbnd,f1,f2,h,hopt,hphi)

c     chcore  implements algorithm  fd, the method described in
c     gill, p.e., murray, w., saunders, m.a., and wright, m. h.,
c     computing forward-difference intervals for numerical optimization,
c     siam journal on scientific and statistical computing, vol. 4,
c     pp. 310-321, june 1983.
c
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

c     bndlo, bndup, and rho control the logic of the routine.
c     bndlo and bndup are the lower and upper bounds that define an
c     acceptable value of the bound on the relative condition error in
c     the second derivative estimate.
c
c     the scalar rho is the factor by which the interval is multiplied
c     or divided, and also the multiple of the well-scaled interval
c     that is used as the initial trial interval.
c-----------------------------------------------------------------------
      implicit none

      double precision  bndlo, bndup
      parameter         (bndlo=1.0d-3,bndup=1.0d-1)
      double precision  zero, sixth, fourth
      parameter         (zero=0.0d+0,sixth=1.6d-1,fourth=2.5d-1)
      double precision  half, two
      parameter         (half=5.0d-1,two=2.0d+0)
      double precision  three, four, ten
      parameter         (three=3.0d+0,four=4.0d+0,ten=1.0d+1)
c     .. scalar arguments ..
      double precision  cdest, epsa, epsr, errbnd, f1, f2, fdest, fx, h,
     *                  hopt, hphi, sdest
      integer           inform, iter, itmax
      logical           debug, done, first
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      double precision  afdmin, cdsave, err1, err2, fdcerr, fdest2,
     *                  fdsave, hsave, oldcd, oldh, oldsd, rho, sdcerr,
     *                  sdsave
      logical           ce1big, ce2big, overfl, te2big
c     .. external functions ..
      double precision  adivb
      external          adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
c     .. save statement ..
      save              cdsave, fdsave, hsave, oldh, rho, sdsave,
     *                  ce1big, ce2big, te2big
c----------------------------------------------------------------------
      iter = iter + 1

c     compute the forward-,  backward-,  central-  and second-order
c     difference estimates.

      fdest = adivb(f1-fx,h,overfl)
      fdest2 = adivb(f2-fx,two*h,overfl)

      oldcd = cdest
      cdest = adivb(four*f1-three*fx-f2,two*h,overfl)

      oldsd = sdest
      sdest = adivb(fx-two*f1+f2,h*h,overfl)

c     compute  fdcerr  and  sdcerr,  bounds on the relative condition
c     errors in the first and second derivative estimates.

      afdmin = min(abs(fdest),abs(fdest2))
      fdcerr = adivb(epsa,half*abs(h)*afdmin,overfl)
      sdcerr = adivb(epsa,fourth*abs(sdest)*h*h,overfl)

c     select the correct case.

      if (first) then

c        first time through.
c        check whether sdcerr lies in the acceptable range.

         first = .false.
         done = sdcerr .ge. bndlo .and. sdcerr .le. bndup
         te2big = sdcerr .lt. bndlo
         ce2big = sdcerr .gt. bndup
         ce1big = fdcerr .gt. bndup
c
         if ( .not. ce1big) then
            hsave = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if
c
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

99999 format (/' //chcore//  itn ',i3,' fx     h',11x,1p,2d16.6,/' //c',
     *       'hcore//  f1      fdest',14x,1p,2d16.6,/' //chcore//  f2 ',
     *       '     fdest2',13x,1p,2d16.6,/' //chcore//  cdest   sdest',
     *       14x,1p,2d16.6,/' //chcore//  fdcerr  sdcerr',13x,1p,2d16.6)
99998 format (' //chcore//  ce1big  ce2big  te2big',5x,3l2)
99997 format (' //chcore//  inform  hopt    errbnd',i5,1p,2d16.6)
      end

      logical function e04udx(string)

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

c        name    dimension  type  i/o/s  description
c        e04udx              l      o    set .true. if string appears
c                                        to be numerical data.
c        string              c    i      input data to be tested.

c        (1)  it is assumed that string is a token extracted by
c             e04udv, which will have converted any lower-case
c             characters to upper-case.

c        (2)  e04udv pads string with blanks, so that a genuine
c             number is of the form  '1234        '.
c             hence, the scan of string stops at the first blank.

c        (3)  complex data with parentheses will not look numeric.
c-----------------------------------------------------------------------
      implicit none

      character*(*)           string
c     .. local scalars ..
      integer                 j, length, ndigit, nexp, nminus, nplus,
     *                        npoint
      logical                 number
      character*1             atom
c----------------------------------------------------------------------
      ndigit = 0
      nexp = 0
      nminus = 0
      nplus = 0
      npoint = 0
      number = .true.
      length = len(string)
      j = 0

   20 j = j + 1
      atom = string(j:j)
      if (lge(atom,'0') .and. lle(atom,'9')) then
c        if (atom.ge.'0' .and. atom.le.'9') then
         ndigit = ndigit + 1
      else if (atom.eq.'d' .or. atom.eq.'e'.or.
     *         atom.eq.'D' .or. atom.eq.'E') then
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

      if (number .and. j.lt.length) go to 20

      e04udx = number .and. ndigit .ge. 1 .and. nexp .le. 1 .and.
     *         nminus .le. 2 .and. nplus .le. 2 .and. npoint .le. 1

c     end of  e04udx. (opnumb)

      end

      subroutine e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)
c----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
      implicit none

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  bigalf, bigbnd, palfa1, palfa2, pnorm
      integer           jadd1, jadd2, n, nctotl
      logical           firstv, negstp
c     .. array arguments ..
      double precision  anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer           istate(nctotl)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)
c     .. local scalars ..
      double precision  absatp, atp, atx, res, rownrm
      integer           i, j, js
      logical           lastv
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------
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
            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

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

c     end of  e04uch. (cmalf1)

99999 format (/'    j  js         featol        res             ap    ',
     *       ' jadd1       palfa1     jadd2       palfa2',/)
99998 format (i5,i4,3g15.5,2(i6,g17.7))
      end

      subroutine e04nbv(n,nu,nrank,ldr,lenv,lenw,r,u,v,w,c,s)
c----------------------------------------------------------------------
c     e04nbv  modifies the  nrank*n  upper-triangular matrix  r  so that
c     q*(r + v*w')  is upper triangular,  where  q  is orthogonal,
c     v  and  w  are vectors, and the modified  r  overwrites the old.
c     q  is the product of two sweeps of plane rotations (not stored).
c     if required,  the rotations are applied to the nu columns of
c     the matrix  u.
c
c     the matrix v*w' is an (lenv) by (lenw) matrix.
c     the vector v is overwritten.
c-----------------------------------------------------------------------
      implicit none

      integer           ldr, lenv, lenw, n, nrank, nu
c     .. array arguments ..
      double precision  c(n), r(ldr,*), s(n), u(n,*), v(n), w(n)
c     .. local scalars ..
      integer           j
c----------------------------------------------------------------------
      j = min(lenv,nrank)
      if (nrank.gt.0) then

c        reduce  v to beta*e( j )  using a backward sweep of rotations
c        in planes (j-1, j), (j-2, j), ..., (1, j).

         call f06fqf('fixed','backwards',j-1,v(j),v,1,c,s)

c        apply the sequence of rotations to u.

         if (nu.gt.0) call f06qxf('l','bottom','backwards',j,nu,1,j,
     *                            c,s,u,n)

c        apply the sequence of rotations to r. this generates a spike in
c        the j-th row of r, which is stored in s.

         call f06qwf('l',n,1,j,c,s,r,ldr)

c        form  beta*e(j)*w' + r.  this a spiked matrix, with a row
c        spike in row j.

         call daxpy (min(j-1,lenw),v(j),w,1,s,1)
         call daxpy (lenw-j+1,v(j),w(j),1,r(j,j),ldr)

c        eliminate the spike using a forward sweep of rotations in
c        planes (1, j), (2, j), ..., (j-1, j).

         call f06qsf('l',n,1,j,c,s,r,ldr)

c        apply the rotations to u.

         if (nu.gt.0) call f06qxf('l','bottom','f',j,nu,1,j,c,
     *                            s,u,n)
      end if

c     end of  e04nbv. (cmr1md)

      end

      subroutine e04mfy(nrz,ldr,r,rzz)
c----------------------------------------------------------------------
c     e04mfy  loads the last column of the  nrz x nrz  triangular factor
c     rz  with the multiple  rzz  of the  nrz-th unit vector.
c-----------------------------------------------------------------------
      implicit none
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  rzz
      integer           ldr, nrz
c     .. array arguments ..
      double precision  r(ldr,*)
c----------------------------------------------------------------------
      if (nrz.eq.0) return

      call sload (nrz-1,zero,r(1,nrz),1)
      r(nrz,nrz) = rzz

c     end of  e04mfy.  (lpcolr)

      end

      subroutine e04udt(inform,n,nclin,ncnln,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,anorm,adx,ax,bl,bu,dslk,dx,slk,x)

c     e04udt  finds a step alfa such that the point x + alfa*p reaches
c     one of the slacks or linear constraints.  the step alfa is the
c     maximum step that can be taken without violating one of the slacks
c     or linear constraints that is currently satisfied.
c-----------------------------------------------------------------------
      implicit none

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  alfa, alfmax, alfmin, bigbnd, dxnorm
      integer           inform, n, nclin, ncnln
c     .. array arguments ..
      double precision  adx(*), anorm(*), ax(*), bl(*), bu(*), dslk(*),
     *                  dx(n), slk(*), x(n)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      double precision  adxi, axi, res, rownrm
      integer           i, j
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
c----------------------------------------------------------------------
      alfa = alfmax
      j = 1

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
c
            i = j - n - nclin
            axi = slk(i)
            adxi = dslk(i)
            rownrm = one
         end if
c
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

      end if

c     determine alfa, the bound on the step to be taken.

      alfa = max(alfa,alfmin)

      inform = 0
      if (alfa.ge.alfmax) inform = 1

c     end of  e04udt. (npalf)

99999 format (/' e04udt entered',/'    j            res             ap',
     *       '           alfa ',/)
99998 format (i5,3g15.5)
99997 format (/' //e04udt//  no finite step.',/' //e04udt//           ',
     *       '  alfa',/' //e04udt//  ',g15.4)
      end

      subroutine e04udy(ndict,dictry,alpha,key,entry)
c----------------------------------------------------------------------
c     description and usage:

c     does dictionary lookups.  a pointer is returned if a
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
c     the ordered search stops after the section of the dictionary
c     having the same first letters as the key has been checked, or
c     after a specified number of entries have been examined.  a special
c     dictionary entry, the currency symbol '$', will also terminate the
c     search.  this will speed things up if an appropriate dictionary
c     length parameter cannot be determined.  both types of search are
c     sequential.  see 'notes' below for some suggestions if efficiency
c     is an issue.

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

c         we have assumed that the dictionary is not too big.  if
c         many searches are to be done or if the dictionary has more
c         than a dozen or so entries, it may be advantageous to build
c         an index array of pointers to the beginning of the section
c         of the dictionary containing each letter, then pass in the
c         portion of the dictionary beginning with dictry (index).
c         (this won't generally work for dictionaries with synonyms.)
c         for very large problems, a completely different approach may
c         be advisable, e.g. a binary search for ordered dictionaries.

c         the key need not be left-justified.  any leading (or
c         trailing) characters which are 'non-significant' to e04udw
c         will be ignored.  these include blanks, horizontal tabs,
c         commas, colons, and equal signs.  see e04udw for details.

c         parameter numsig sets a limit on the length of significant
c         dictionary entries.  special applications may require that
c         this be increased.  (it is 16 in the present version.)

c         the handling of ambiguities introduces some ambiguity:

c            alpha = .true.  a potential problem, when one entry
c                            looks like an abbreviation for another
c                            (eg. does 'a' match 'a' or 'ab') was
c                            resolved by dropping out of the search
c                            immediately when an 'exact' match is found.

c            alpha = .false. the programmer must ensure that the above
c                            situation does not arise: each dictionary
c                            entry must be recognizable, at least when
c                            specified to full length.  otherwise, the
c                            result of a search will depend on the
c                            order of entries.
c-----------------------------------------------------------------------
      implicit none

      character         blank, curly
      integer           numsig
      parameter         (blank=' ',curly='$',numsig=16)
c     .. scalar arguments ..
      integer           entry, ndict
      logical           alpha
      character*(*)     key
c     .. array arguments ..
      character*(*)     dictry(ndict)
c     .. local scalars ..
      integer           first, i, ifrst, ifrst1, ilast, ilast1, ilen,
     *                  ilen1, ilst, imark, imark1, last, length, mark
      character*16      flag, target, trgt, trgt1
c----------------------------------------------------------------------

      entry = 0

c     isolate the significant portion of the input key (if any).

      first = 1
      last = min(len(key),numsig)
      call e04udw(key,first,last,mark)

      if (mark.gt.0) then
         target = key(first:mark)

c        look up target in the dictionary.

   20    continue
         length = mark - first + 1

c           select search strategy by cunning choice of termination test
c           flag.  the left curly bracket follows all the alphabetic
c           characters in the ascii collating sequence, but precedes the
c           vertical bar.

         if (alpha) then
            flag = target
         else
            flag = curly
         end if

c           perform search.

         i = 0
   40    continue
         i = i + 1
         if (target(1:length).eq.dictry(i)(1:length)) then
            if (entry.eq.0) then

c                    first 'hit' - must still guard against ambiguities
c                    by searching until we've gone beyond the key
c                    (ordered dictionary) or until the end-of-dictionary
c                    mark is reached (exhaustive search).

               entry = i

c                    special handling if match is exact - terminate
c                    search.  we thus avoid confusion if one dictionary
c                    entry looks like an abbreviation of another.
c                    this fix won't generally work for un-ordered
c                    dictionaries.

               first = 1
               last = numsig
               call e04udw(dictry(entry),first,last,mark)
               if (mark.eq.length) i = ndict
            else
c                    if two hits check if they are attempting to
c                    indicate the same dictionary entry.
c                    extract keyword from first match found

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

c                    extract keyword from next match found

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

c                    if keywords not identical then ambiguity

               if (trgt(1:ilen).ne.trgt1(1:ilen1)) then

c                       oops - two hits.  abnormal termination.

                  entry = -i
                  return
               end if
            end if
         end if

c           check whether we've gone past the appropriate section of the
c           dictionary.  the test on the index provides insurance and an
c           optional means for limiting the extent of the search.

         if (lle(dictry(i)(1:length),flag) .and. i.lt.ndict) go to 40

c           check for a synonym.

         if (entry.gt.0) then

c              look for a second entry 'behind' the first entry.  first
c              and mark were determined above when the hit was detected.

            first = mark + 2
            call e04udw(dictry(entry),first,last,mark)
            if (mark.gt.0) then

c                re-set target and dictionary pointer, then repeat the
c                search for the synonym instead of the original key.

               target = dictry(entry) (first:mark)
               entry = 0
               go to 20

            end if
         end if

      end if

      if (entry.gt.0) key = dictry(entry)

c     end of e04udy.  (cmlook/oplook)

      end

      subroutine e04nbt(mode,nrowt,n,t,y)
c----------------------------------------------------------------------
c     e04nbt  solves equations involving a reverse-triangular matrix  t
c     and a right-hand-side vector  y,  returning the solution in  y.
c-----------------------------------------------------------------------
      implicit none

      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      integer           mode, n, nrowt
c     .. array arguments ..
      double precision  t(nrowt,*), y(n)
c     .. local scalars ..
      double precision  yj
      integer           j, jj, l, n1
c----------------------------------------------------------------------
      n1 = n + 1
      if (mode.eq.1) then

c        mode = 1    solve  t * y(new) = y(old).

         do 20 j = 1, n
            jj = n1 - j
            yj = y(j)/t(j,jj)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.zero) call daxpy (l,(-yj),t(j+1,jj),
     *                                               1,y(j+1),1)
   20    continue
      else

c        mode = 2    solve  t' y(new) = y(old).

         do 40 j = 1, n
            jj = n1 - j
            yj = y(j)/t(jj,j)
            y(j) = yj
            l = jj - 1
            if (l.gt.0 .and. yj.ne.zero) call daxpy (l,(-yj),t(jj,j+1),
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

c     end of  e04nbt. (cmtsol)

      end

      subroutine e04nbu(n,nu,nrank,ldr,i,j,r,u,c,s)
c----------------------------------------------------------------------
c     e04nbu  interchanges the  i-th  and  j-th  (i .lt. j)  columns of
c     an  nrank*n  upper-trapezoidal matrix  r   and restores the
c     resulting matrix to upper-trapezoidal form using two sweeps of
c     plane rotations applied on the left.  r is overwritten.
c
c     if nu .gt. 0,  the rotations are applied to the  nu  columns of
c     the matrix  u.
c-----------------------------------------------------------------------
      implicit none

      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      integer           i, j, ldr, n, nrank, nu
c     .. array arguments ..
      double precision  c(n), r(ldr,*), s(n), u(n,*)
c     .. local scalars ..
      integer           lenj
c----------------------------------------------------------------------
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

         call f06fqf('fixed','backwards',lenj-i-1,r(lenj,j),r(i+1,j),1,
     *               c(i+1),s(i+1))

         if (nu.gt.0) call f06qxf('l','bottom','backwards',n,nu,i+1,
     *                            lenj,c,s,u,n)

c        put zeros into the j-th column of r in positions corresponding
c        to the sub-diagonals of the i-th column.

         s(i) = r(lenj,j)
         call sload (lenj-i,zero,r(i+1,j),1)

c        apply the sequence of rotations to r.  this generates a spike
c        in the lenj-th row of r, which is stored in s.

         call f06qwf('l',n,i+1,lenj,c,s,r,ldr)

c        eliminate the spike using a forward sweep in planes
c        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
c        if necessary, apply the sequence of rotations to u.

         call f06qsf('l',n,i,lenj,c,s,r,ldr)

         if (nu.gt.0) call f06qxf('l','bottom','f',lenj,nu,i,
     *                            lenj,c,s,u,n)
      end if

c     end of  e04nbu

      end

      subroutine e04nct(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,res,r,t,gq,
     *                  zy,c,s)
c----------------------------------------------------------------------
c     e04nct  updates the least-squares factor r and the factorization
c     a(free) (z y) = (0 t) when a regular, temporary or artificial
c     constraint is deleted from the working set.
c-----------------------------------------------------------------------
      implicit none

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      integer           jdel, kdel, lda, ldr, ldt, ldzy, n, nactiv,
     *                  nfree, ngq, nrank, nres, nrz, nz
      logical           unitq
c     .. array arguments ..
      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), zy(ldzy,*)
      integer           kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  cs, sn
      integer           i, ir, itdel, jart, k, ka, ld, npiv, nrz1, nsup,
     *                  nt
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin
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
               if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,nfree,ir,r,
     *                                     res,c,s)
               call dswap(ngq,gq(nfree,1),n,gq(ir,1),n)
            end if
c
            if ( .not. unitq) then

c              copy the incoming column of  a(free)  into the end of t.

               do 20 ka = 1, nactiv
                  i = kactiv(ka)
                  t(ka,nfree) = a(i,jdel)
   20          continue

c              expand q by adding a unit row and column.

               if (nfree.gt.1) then
                  call sload (nfree-1,zero,zy(nfree,1),ldzy)
                  call sload (nfree-1,zero,zy(1,nfree),1)
               end if
               zy(nfree,nfree) = one
            end if
         else

c           case 2.  a general constraint has been deleted.

            itdel = kdel
            nactiv = nactiv - 1

c           delete row  kdel  of t and move up the ones below it.
c           t becomes reverse lower hessenberg.

            do 40 i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld = nfree - i
               call dcopy (i+1,t(i+1,ld),ldt,t(i,ld),ldt)
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
                  call dcopy (nsup-1,t(nactiv-1,nz+1),ldt-1,s(nz+1),1)
                  call f06qzz('remove',nactiv,1,nsup,c(nz+1),s(nz+1),
     *                        t(1,nz+1),ldt)
               end if

               call f06baf(t(nactiv,nz+1),t(nactiv,nz),cs,sn)
               t(nactiv,nz) = zero
               s(nz) = -sn
               c(nz) = cs

               call f06qxf('r','v','backwards',nfree,nfree,
     *                     nz,npiv,c,s,zy,ldzy)
               call f06qxf('left ','v','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gq,n)

               nt = min(nrank,npiv)

               if (nt.lt.npiv .and. nt.gt.0) then

c                 r is upper trapezoidal, pretend r is (nt x n) and
c                 apply the rotations in columns  max(nt,nz)  thru npiv.

                  call f06qxf('r','v','backwards',nt,n,
     *                        max(nt,nz),npiv,c,s,r,ldr)
               end if

c              apply the column transformations to the triangular part
c              of r.  a subdiagonal element is generated that must be
c              eliminated by a row rotation before the next column
c              transformation can be applied.

               if (nz.lt.nt) call f06qtf('r',nt,nz,nt,c,s,r,ldr)


c              apply the row rotations to the remaining rows of r.

               call f06qxf('l','v','backwards',nt,n-nt,nz,nt,
     *                     c,s,r(1,min(nt+1,n)),ldr)

               if (nres.gt.0) call f06qxf('l','v','backwards',
     *                                    nt,nres,nz,nt,c,s,res,n)

            end if
            call f06flf(nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
         end if
      end if

      nrz1 = nrz + 1
c
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
            if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,nrz1,jart,r,
     *                                  res,c,s)
         end if
      end if

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

      subroutine e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,bl,bu,rpq,rpq0,dx,gq,r,t,zy,work)
c----------------------------------------------------------------------
c     e04ucm   defines a point which lies on the initial working set for
c     the qp subproblem.  this routine is similar to e04nch except
c     that advantage is taken of the fact that the initial estimate of
c     the solution of the least-squares subproblem is zero.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  dxnorm, gdx
      integer           ldaqp, ldr, ldt, ldzy, n, nactiv, ncqp, nctotl,
     *                  nfree, nlnx, nz
      logical           unitq
c     .. array arguments ..
      double precision  adx(*), aqp(ldaqp,*), bl(nctotl), bu(nctotl),
     *                  dx(n), gq(n), r(ldr,*), rpq(nlnx), rpq0(nlnx),
     *                  t(ldt,*), work(n), zy(ldzy,*)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  bnd
      integer           i, j, k, nfixed, nr
c     .. external functions ..
      double precision  ddot, dnrm2
      external          ddot, dnrm2
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
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
         if (nactiv.gt.0) call e04nbt(1,ldt,nactiv,t(1,nz+1),work(nz+1))
         call dcopy (nactiv+nfixed,work(nz+1),1,dx(nz+1),1)
         if (nz.gt.0) call sload (nz,zero,dx,1)
c
         gdx = ddot(nactiv+nfixed,gq(nz+1),1,dx(nz+1),1)
c
         if (nz.lt.n) then
            call dgemv ('n',nz,n-nz,-one,r(1,nz+1),ldr,dx(nz+1),1,one,
     *                 rpq,1)
            if (nz.lt.nlnx) then
               nr = ldr
               if (nz+1.eq.n) nr = 1
               call dcopy (nlnx-nz,dx(nz+1),1,rpq(nz+1),1)
               call dscal (nlnx-nz,(-one),rpq(nz+1),1)
               call dtrmv ('u','n','n',nlnx-nz,r(nz+1,nz+1),nr,rpq(nz+1)
     *                    ,1)
               if (nlnx.lt.n) then
                  nr = ldr
                  if (nlnx+1.eq.n) nr = n - nz
                  call dgemv ('n',nlnx-nz,n-nlnx,-one,r(nz+1,nlnx+1),nr,
     *                       dx(nlnx+1),1,one,rpq(nz+1),1)
               end if
            end if
         end if

         call e04nbw(2,n,nz,nfree,ldzy,unitq,kx,dx,zy,work)
      end if

c     compute the 2-norm of  dx.
c     initialize  a*dx.

      dxnorm = dnrm2(n,dx,1)
      if (ncqp.gt.0) call dgemv ('n',ncqp,n,one,aqp,ldaqp,dx,1,zero,adx,
     *                          1)

c     end of  e04ucm. (npsetx)

99999 format (/' //e04ucm// variables after e04ucm ... ')
99998 format (5g12.3)
      end

      subroutine e04nck(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,
     *                  nrz,istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,
     *                  jtiny,jbigst,kbigst,trulam,a,anorms,gq,rlamda,t,
     *                  wtinf)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)
c     .. scalar arguments ..
      double precision  dinky, trulam
      integer           jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, msglvl, n, nactiv, nfree, nrz, numinf,
     *                  nz
      character*2       prbtyp
c     .. array arguments ..
      double precision  a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  anormj, biggst, blam, rlam, scdlam, smllst,
     *                  tinylm
      integer           i, is, j, k, l, nfixed
c     .. local arrays ..
      character*80      rec(80)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
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
c
      end if
c
c     
c     compute jsmlst for regular constraints and temporary bounds.
c     
c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.
c
      if (n.gt.nz) call dcopy (n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call e04nbt(2,ldt,nactiv,t(1,nz+1),rlamda)
c

c     set elements nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.

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

c     find jsmlst and ksmlst.

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

c     if required, print the multipliers.

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

      subroutine e04mfq(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,bl,
     *                  bu,featol,x)
c----------------------------------------------------------------------
c     e04mfq  checks the residuals of the constraints that are believed
c     to be feasible.  the number of constraints violated by more than
c     featol is computed, along with the maximum constraint violation.

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, errmax
      integer           jmax, n, nclin, nviol
c     .. array arguments ..
      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)
c     .. local scalars ..
      double precision  biglow, bigupp, con, feasj, res
      integer           is, j
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
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

      do 40 j = 1, n + nclin
         is = istate(j)

         if (is.ge.0) then
            feasj = featol(j)

            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if

c           check for constraint violations.

            if (bl(j).gt.biglow) then
               res = bl(j) - con
               if (res.gt.feasj) then
                  nviol = nviol + 1
                  go to 20
               end if
            end if

            if (bu(j).lt.bigupp) then
               res = bu(j) - con
               if (res.lt.(-feasj)) then
                  nviol = nviol + 1
                  res = -res
                  go to 20
               end if
            end if

c           this constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.

            res = zero

            if (is.eq.1) then
               res = abs(bl(j)-con)

            else if (is.eq.2) then
               res = abs(bu(j)-con)

            else if (is.eq.3) then
               res = abs(bu(j)-con)
            end if

            if (res.gt.feasj) nviol = nviol + 1

   20       if (res.gt.errmax) then
               jmax = j
               errmax = res
            end if
         end if
   40 continue

c     end of  e04mfq.  (cmfeas)

99999 format (/' //e04mfq//  the maximum violation is ',1p,d14.2,' in ',
     *       'constraint',i5)
      end

      subroutine e04udv(string,numin,numout,list)
c----------------------------------------------------------------------
c     an aid to parsing input data.  the individual 'tokens' in a
c     character string are isolated, converted to uppercase, and stored
c     in an array.  here, a token is a group of significant, contiguous
c     characters.  the following are non-significant, and hence may
c     serve as separators:  blanks, horizontal tabs, commas, colons,
c     and equal signs.  see e04udw for details.  processing continues
c     until the requested number of tokens have been found or the end
c     of the input string is reached.

c     parameters:

c     name    dimension  type  i/o/s  description
c     string              c    i      input string to be analyzed.
c     numin               i    i/o    number of tokens requested (input)
c     numout                          and found (output).
c     (numin and numout were both called number in the original)
c
c     list    numin       c      o    array of tokens, changed to upper
c                                    case.

      character         blank
      parameter         (blank=' ')
c     .. scalar arguments ..
      integer           numin, numout
      character*(*)     string
c     .. array arguments ..
      character*(*)     list(numin)
c     .. local scalars ..
      integer           count, first, i, last, mark
c----------------------------------------------------------------------
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
c           pass token to output string array
c
         list(count) = string(first:mark)
         first = mark + 2
         if (count.lt.numin) go to 20
c
      end if

c     fill the rest of list with blanks and set number for output.

      do 40 i = count + 1, numin
         list(i) = blank
   40 continue

      numout = count

c     end of  e04udv. (optokn)
      end

      subroutine e04uck(first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *                  toltny,alfa,alfbst,fbest,gbest)
c----------------------------------------------------------------------
c     e04uck  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     e04uck  requires both  f(alpha)  and  f'(alpha) to be evaluated at
c     points in the interval.  estimates of the minimizer are computed
c     using safeguarded cubic interpolation.

c     reverse communication is used to allow the calling program to
c     evaluate f and f'.  some of the parameters must be set or tested
c     by the calling program.  the remainder would ordinarily be local
c     variables.

c     input parameters (relevant to the calling program)

c     first         must be .true. on the first entry.
c                   it is subsequently altered by e04uck.

c     debug         specifies whether detailed output is wanted.

c     maxf          is an upper limit on the number of times e04uck is
c                   to be entered consecutively with done = .false.
c                   (following an initial entry with first = .true.).

c     alfa          is the first estimate of a minimizer.  alfa is
c                   subsequently altered by e04uck (see below).

c     alfmax        is the upper limit of the interval to be searched.

c     epsaf         is an estimate of the absolute precision in the
c                   computed value of f(0).

c     ftry, gtry    are the values of f, f'  at the new point
c                   alfa = alfbst + xtry.

c     g0            is the value of f'(0).  g0 must be negative.

c     tolabs,tolrel define a function tol(alfa) = tolrel*alfa + tolabs
c                   such that if f has already been evaluated at alfa,
c                   it will not be evaluated closer than tol(alfa).
c                   these values may be reduced by e04uck.

c     targtg        is the target value of abs(f'(alfa)). the search
c                   is terminated when
c                    abs(f'(alfa)) le targtg and f(alfa) lt 0.

c     toltny        is the smallest value that tolabs is allowed to be
c                   reduced to.

c     output parameters (relevant to the calling program)

c     imprvd        is .true. if the previous alfa was the best point so
c                   far.  any related quantities should be saved by the
c                   calling program (e.g., gradient arrays) before
c                   paying attention to the variable done.

c     done = .false.  means the calling program should evaluate
c                      ftry = f(alfa),  gtry = f'(alfa)
c                   for the new trial alfa, and re-enter e04uck.

c     done = .true.   means that no new alfa was calculated.  the value
c                   of inform gives the result of the search as follows
c
c                   inform = 1 means the search has terminated
c                              successfully with alfbst < alfmax.

c                   inform = 2 means the search has terminated
c                              successfully with alfbst = alfmax.

c                   inform = 3 means that the search failed to find a
c                              point of sufficient decrease in maxf
c                              functions, but a lower point was found.

c                   inform = 4 means alfmax is so small that a search
c                              should not have been attempted.

c                   inform = 5 is never set by e04uck.

c                   inform = 6 means the search has failed to find a
c                              useful step.  the interval of uncertainty
c                              is [0,b] with b < 2*tolabs. a minimizer
c                              lies very close to alfa = 0, or f'(0) is
c                              not sufficiently accurate.

c                   inform = 7 if no better point could be found after
c                              maxf  function calls.

c                   inform = 8 means the input parameters were bad.
c                              alfmax le toltny  or g0 ge zero.
c                              no function evaluations were made.

c     numf          counts the number of times e04uck has been entered
c                   consecutively with done = .false. (i.e., with a new
c                   function value ftry).

c     alfa          is the point at which the next function ftry and
c                   derivative gtry must be computed.

c     alfbst        should be accepted by the calling program as the
c                   approximate minimizer, whenever e04uck returns
c                   inform = 1 or 2 (and possibly 3).

c     fbest, gbest  will be the corresponding values of f, f'.

c     the following parameters retain information between entries

c     braktd        is .false. if f and f' have not been evaluated at
c                   the far end of the interval of uncertainty.  in this
c                   case, the point b will be at alfmax + tol(alfmax).

c     crampd        is .true. if alfmax is very small (le tolabs).  if
c                   the search fails, this indicates that a zero step
c                   should be taken.

c     extrap        is .true. if xw lies outside the interval of
c                   uncertainty.  in this case, extra safeguards are
c                   applied to allow for instability in the polynomial
c                   fit.

c     moved         is .true. if a better point has been found, i.e.,
c                   alfbst gt 0.

c     wset          records whether a second-best point has been
c                   determined it will always be .true. when convergence
c                   is tested.

c     nsamea        is the number of consecutive times that the
c                   left-hand end point of the interval of uncertainty
c                   has remained the same.

c     nsameb        similarly for the right-hand end.

c     a, b, alfbst  define the current interval of uncertainty.
c                   a minimizer lies somewhere in the interval
c                   [alfbst + a, alfbst + b].

c     alfbst        is the best point so far.  it is always at one end
c                   of the interval of uncertainty.  hence we have
c                   either  a lt 0,  b = 0  or  a = 0,  b gt 0.

c     fbest, gbest  are the values of f, f' at the point alfbst.

c     factor        controls the rate at which extrapolated estimates
c                   of alfa may expand into the interval of uncertainty.
c                   factor is not used if a minimizer has been bracketed
c                   (i.e., when the variable braktd is .true.).

c     fw, gw        are the values of f, f' at the point alfbst + xw.
c                   they are not defined until wset is .true..

c     xtry          is the trial point in the shifted interval (a, b).

c     xw            is such that  alfbst + xw  is the second-best point.
c                   it is not defined until  wset  is .true..
c                   in some cases,  xw  will replace a previous  xw
c                   that has a lower function but has just been excluded
c                   from the interval of uncertainty.

c     local variables
c     ===============

c     closef     is .true. if the new function ftry is within epsaf of
c                fbest (up or down).

c     found      is .true. if the sufficient decrease conditions hold at
c                alfbst.

c     quitf      is .true. when  maxf  function calls have been made.

c     quiti      is .true. when the interval of uncertainty is less than
c                2*tol.
c----------------------------------------------------------------------
      double precision  zero, point1, half
      parameter         (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision  one, three, five
      parameter         (one=1.0d+0,three=3.0d+0,five=5.0d+0)
      double precision  ten, eleven
      parameter         (ten=1.0d+1,eleven=1.1d+1)
c     .. scalar arguments ..
      double precision  alfa, alfbst, alfmax, epsaf, fbest, ftry, g0,
     *                  gbest, gtry, targtg, tolabs, tolrel, toltny
      integer           inform, maxf, nout, numf
      logical           debug, done, first, imprvd
c     .. local scalars ..
      double precision  a, absr, artifa, artifb, b, daux, dtry, factor,
     *                  fw, gw, q, r, s, scale, tol, tolmax, truea,
     *                  trueb, xmidpt, xtry, xw
      integer           nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, fitok,
     *                  found, moved, quitf, quiti, setxw, wset
c     .. save statement ..
      save              braktd, crampd, extrap, moved, wset, nsamea,
     *                  nsameb, a, b, factor, xtry, xw, fw, gw, tolmax
c----------------------------------------------------------------------
      badfun = .false.
      quitf = .false.
      quiti = .false.
      imprvd = .false.

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

         if ( .not. braktd) then
            tolmax = tolabs + tolrel*alfmax
            b = alfmax - alfbst + tolmax
         end if

c        see if the new step is better.  if alfa is large enough that
c        ftry can be distinguished numerically from zero,  the function
c        is required to be sufficiently negative.

         closef = abs(ftry-fbest) .le. epsaf
         if (closef) then
            imprvd = abs(gtry) .le. abs(gbest)
         else
            imprvd = ftry .lt. fbest
         end if

         if (imprvd) then

c           we seem to have an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.

            fw = fbest
            fbest = ftry
            gw = gbest
            gbest = gtry
            alfbst = alfa
            moved = .true.

            a = a - xtry
            b = b - xtry
            xw = zero - xtry
            wset = .true.
            extrap = xw .lt. zero .and. gbest .lt. zero .or. xw .gt.
     *               zero .and. gbest .gt. zero

c           decrease the length of the interval of uncertainty.

            if (gtry.le.zero) then
               a = zero
               nsamea = 0
            else
               b = zero
               nsameb = 0
               braktd = .true.
            end if
         else

c           the new function value is not better than the best point so
c           far.  the origin remains unchanged but the new point may
c           qualify as xw.  xtry must be a new bound on the best point.

            if (xtry.le.zero) then
               a = xtry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if

c           if xw has not been set or ftry is better than fw, update the
c           points accordingly.

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

c        check the termination criteria.  wset will always be .true..

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b

         found = abs(gbest) .le. targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol

         if (quiti .and. .not. moved) then

c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.

            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if

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

c              compute  q =  the square root of  r*r - gbest*gw.
c              the method avoids unnecessary underflow and overflow.

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

c                 compute a minimizer of the fitted cubic.

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

c              a minimizer has not been bracketed.  set an artificial
c              upper bound by expanding the interval  xw  by a suitable
c              factor.

               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = five*factor

            else if (extrap) then

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

c              the points are configured for an interpolation.  the
c              default value xtry bisects the interval of uncertainty.
c              the artificial interval is just (a, b).

               xtry = xmidpt

               if (nsamea.ge.3 .or. nsameb.ge.3) then

c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation.

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

c                 accept the polynomial fit.

                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
               end if
            end if
         end if
      end if

      if ( .not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then

c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)

            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
               if (half*(a+b).le.zero) then
                  xtry = -tol
               else
                  xtry = tol
               end if
               alfa = alfbst + xtry
            end if
         else

c           the step is close to, or larger than alfmax, replace it by
c           alfmax to force evaluation of  f  at the boundary.

            braktd = .true.
            xtry = alfmax - alfbst
            alfa = alfmax
         end if
      end if

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
99991 format (' -',/)
      end

      subroutine e04uds(centrl,inform,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,c2,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)

c     e04uds evaluates any missing gradients.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  three, four
      parameter         (three=3.0d+0,four=4.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, cdint, fdint, fdnorm, objf
      integer           inform, ldcj, ldcju, lenw, n, ncnln
      logical           centrl
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  biglow, bigupp, delta, objf1, objf2, stepbl,
     *                  stepbu, xj
      integer           i, j, mode, ncolj, nstate
c     .. common blocks ..
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
      inform = 0

c     use the pre-assigned difference intervals to approximate the
c     derivatives.

c     use either the same interval for each element (lfdset = 1),
c     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).
c
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

            delta = delta*(one+abs(xj))
            fdnorm = max(fdnorm,delta)
            if (half*(stepbl+stepbu).lt.zero) delta = -delta

            x(j) = xj + delta
            if (ncolj.gt.0) then
               call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,nstate,
     *                     iuser,user)
               if (mode.lt.0) go to 100
            end if

            if (gradu(j).eq.rdummy) then
               call objfun(mode,n,x,objf1,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100
            end if

            if (centrl) then

c              central differences.

               x(j) = xj + delta + delta

               if (ncolj.gt.0) then
                  call confun(mode,ncnln,n,ldcju,needc,x,c2,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 100

                  do 40 i = 1, ncnln
                     if (needc(i).eq.1) cjac(i,j) = (four*c1(i)
     *                   -three*c(i)-c2(i))/(delta+delta)
   40             continue
               end if

               if (gradu(j).eq.rdummy) then
                  call objfun(mode,n,x,objf2,gradu,nstate,iuser,user)
                  if (mode.lt.0) go to 100

                  grad(j) = (four*objf1-three*objf-objf2)/(delta+delta)

               end if
            else

c              forward differences.

               if (ncolj.gt.0) then
                  do 60 i = 1, ncnln
                     if (needc(i).eq.1) cjac(i,j) = (c1(i)-c(i))/delta
   60             continue
               end if

               if (gradu(j).eq.rdummy) grad(j) = (objf1-objf)/delta

            end if
         end if
         x(j) = xj

   80 continue
c
      return

  100 inform = mode

c     end of  e04uds. (npfd)

      end

      subroutine e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,nfree,
     *                  nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,a,r,
     *                  t,res,gq,zy,w,c,s,msglvl)
c----------------------------------------------------------------------
c     e04ncv  updates the factorization,  a(free) * (z y) = (0 t),  when
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  condmx
      integer           iadd, ifix, inform, jadd, lda, ldr, ldt, ldzy,
     *                  msglvl, n, nactiv, nfree, ngq, nrank, nres, nz
      logical           unitq
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  cond, condbd, dtnew, tdtmax, tdtmin
      integer           i, nanew, nfmin, npiv, nt
      logical           bound, overfl
c     .. local arrays ..
      character*80      rec(5)
c     .. external functions ..
      double precision  dnrm2, adivb
      external          dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin
c----------------------------------------------------------------------
c     if the condition estimator of the updated factors is greater than
c     condbd,  a warning message is printed.

      condbd = one/epspt9

      overfl = .false.
      bound = jadd .le. n

      if (bound) then

c        a simple bound has entered the working set.  iadd  is not used.

         nanew = nactiv

         if (unitq) then

c           q  is not stored, but kx defines an ordering of the columns
c           of the identity matrix that implicitly define  q.
c           define the sequence of pairwise interchanges p that moves
c           the newly-fixed variable to position nfree.
c           reorder kx accordingly.

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
c           move the  (nfree)-th  row of  q  to position  ifix.

            call dcopy (nfree,zy(ifix,1),ldzy,w,1)
            if (ifix.lt.nfree) then
               call dcopy (nfree,zy(nfree,1),ldzy,zy(ifix,1),ldzy)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else

c        a general constraint has entered the working set.
c        ifix  is not used.

         nanew = nactiv + 1

c        transform the incoming row of  a  by  q'.  use c as workspace.

         call dcopy (n,a(iadd,1),lda,w,1)
         call e04nbw(8,n,nz,nfree,ldzy,unitq,kx,w,zy,c)

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
c           set. update the estimate of the condition number.

            tdtmax = max(dtnew,dtmax)
            tdtmin = min(dtnew,dtmin)
            cond = adivb(tdtmax,tdtmin,overfl)
         end if

         if (cond.gt.condmx .or. overfl) go to 60

         if (unitq) then

c           first general constraint added.  set  q = i.

            call f06qhf('g',nfree,nfree,zero,one,zy,ldzy)
            unitq = .false.
         end if
      end if

      if (bound) then
         npiv = nfree
      else
         npiv = nz
      end if

      nt = min(nrank,npiv)

      if (unitq) then

c        q (i.e., zy) is not stored explicitly.
c        apply the sequence of pairwise interchanges p that moves the
c        newly-fixed variable to position nfree.

         if (ngq.gt.0) call f06qkf('l','t',nfree-1,w,ngq,gq,
     *                             n)

         if (nrank.gt.0) then

c           apply the pairwise interchanges to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1), s(2), ..., s(nt-1).

            call f06qnz('r',n,ifix,nt,s,r,ldr)

            if (nt.lt.npiv) then

c              r is upper trapezoidal.  apply the interchanges in
c              columns  nt  thru  npiv.

               do 40 i = ifix, nt - 1
                  w(i) = i
   40          continue

               call f06qkf('r','n',nfree-1,w,nt,r,ldr)
            end if

c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.

            call f06qrf('l',n,ifix,nt,c,s,r,ldr)
            if (nres.gt.0) call f06qxf('l','v','f',nt,
     *                                 nres,ifix,nt,c,s,res,n)
         end if
      else

c        full matrix q.  define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1,2), (2,3), ...,
c        (npiv-1,npiv).  the rotations must be applied to zy, r, t
c        and gq'.

         call f06fqf('v','f',npiv-1,w(npiv),w,1,c,s)

         if (bound .and. nactiv.gt.0) then

            call dcopy (nactiv,s(nz),1,w(nz),1)

            s(nz) = s(nz)*t(nactiv,nz+1)
            t(nactiv,nz+1) = c(nz)*t(nactiv,nz+1)

            call f06qzz('c',nactiv,1,nactiv,c(nz+1),s(nz+1),
     *                  t(1,nz+1),ldt)
            call dcopy (nactiv,s(nz),1,t(nactiv,nz),ldt-1)

            call dcopy (nactiv,w(nz),1,s(nz),1)
         end if

         if (ngq.gt.0) call f06qxf('l','v','f',npiv,ngq,1,npiv,c,s,gq,n)
         call f06qxf('r','v','f',nfree,nfree,1,npiv,c,s,zy,ldzy)

         if (nrank.gt.0) then

c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nt-1).

            nt = min(nrank,npiv)
            call f06qvf('r',n,1,nt,c,s,r,ldr)

            if (nt.lt.npiv) then

c              r is upper trapezoidal.  pretend r is (nt x n) and
c              apply the rotations in columns  nt  thru  npiv.

               call f06qxf('r','v','f',nt,n,nt,npiv,c,s,r,ldr)

            end if

c           eliminate the subdiagonal elements of r with a left-hand
c           sweep of rotations p2 in planes (1,2), (2,3), ...,(nt-1,nt).
c           apply p2 to res.

            call f06qrf('l',n,1,nt,c,s,r,ldr)
            if (nres.gt.0) call f06qxf('l','v','f',nt,
     *                                 nres,1,nt,c,s,res,n)
         end if

         if (bound) then

c           the last row and column of zy has been transformed to plus
c           or minus the unit vector e(nfree).  reconstitute the
c           columns of gq and r corresponding to the new fixed variable.

            if (w(nfree).lt.zero) then
               nfmin = min(nrank,nfree)
               if (nfmin.gt.0) call dscal (nfmin,-one,r(1,nfree),1)
               if (ngq.gt.0) call dscal (ngq,-one,gq(nfree,1),n)
            end if


c           the diagonals of t have been altered.  recompute the
c           largest and smallest values.

            if (nactiv.gt.0) then
               call f06flf(nactiv,t(nactiv,nz),ldt-1,tdtmax,tdtmin)
               cond = adivb(tdtmax,tdtmin,overfl)
            end if
         else

c           general constraint.  install the new row of t.

            call dcopy (nanew,w(nz),1,t(nanew,nz),ldt)
         end if
      end if

c     prepare to exit.  check the magnitude of the condition estimator.

   60 return

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

      subroutine e04xaw(inform,lvlder,msglvl,ncset,n,ncnln,ldcj,ldcju,
     *                  bigbnd,epsrf,oktol,fdchk,xnorm,confun,needc,bl,
     *                  bu,c,c1,cjac,cjacu,cjdx,dx,err,x,y,iuser,user)
c----------------------------------------------------------------------
c     e04xaw  checks if the gradients of the constraints have been coded
c     correctly.
c
c     on input,  the values of the constraints at the point x are stored
c     in c.  their corresponding gradients are stored in cjacu.  if any
c     jacobian element has not been specified,  it will have a dummy
c     value.  missing values are not checked.
c
c     a cheap test is first undertaken by calculating the directional
c     derivative using two different methods.  if this proves
c     satisfactory and no further information is desired, e04xaw is
c     terminated. otherwise, chcore is called to give optimal
c     step-sizes and a central-difference approximation to each
c     element of the jacobian for which a test is deemed necessary,
c     either by the program or the user.
c
c     lvrfyc has the following meaning...
c
c     -1        do not perform any check.
c     0        do the cheap test only.
c     2 or 3   do both cheap and full test.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, point9
      parameter         (zero=0.0d+0,half=0.5d+0,point9=0.9d+0)
      double precision  one, two, ten
      parameter         (one=1.0d+0,two=2.0d+0,ten=1.0d+1)
      character*4       lbad, lgood
      parameter         (lbad='bad?',lgood='  ok')
c     .. scalar arguments ..
      double precision  bigbnd, epsrf, fdchk, oktol, xnorm
      integer           inform, ldcj, ldcju, lvlder, msglvl, n, ncnln,
     *                  ncset
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  biglow, bigupp, cdest, cij, cjdiff, cjsize,
     *                  colmax, dxj, dxmult, emax, epsaci, errbnd, f1,
     *                  f2, fdest, h, hopt, hphi, sdest, signh, stepbl,
     *                  stepbu, xj
      integer           i, imax, info, irow, iter, itmax, j, j3, j4,
     *                  jcol, mode, ncheck, ncolj, ngood, nstate, nwrong
      logical           const, debug, done, first, headng, needed, ok
      character*4       key
c     .. local arrays ..
      character*18      result(0:4)
      character*120     rec(4)
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy
      common            /fe04uc/inpdbg, npdbg
c     .. data statements ..
      data              result/'                 ', 'constant?      ',
     *                  'linear or odd?   ', 'too nonlinear?',
     *                  'small derivative?'/
c----------------------------------------------------------------------
      inform = 0
      needed = ncnln .gt. 0 .and. lvrfyc .eq. 0 .or. lvrfyc .eq. 2 .or.
     *         lvrfyc .eq. 3
      if ( .not. needed) return

      if (msglvl.gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,4,rec)
      end if

      nstate = 0

      biglow = -bigbnd
      bigupp = bigbnd

c     perform the cheap test.

      h = (one+xnorm)*fdchk

      if (n.le.100) then
         dxmult = 0.9
      else if (n.le.250) then
         dxmult = 0.99
      else
         dxmult = 0.999
      end if

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

         if (half*(stepbl+stepbu).lt.zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
   60 continue

      if (ncheck.eq.0) then
         if (msglvl.gt.0) then
            write (rec,fmt=99995)
            call x04bay(iprint,2,rec)
         end if
      else

c        compute  (jacobian)*dx.

         call dgemv ('n',ncnln,n,one,cjacu,ldcju,dx,1,zero,cjdx,1)

c        make forward-difference approximation along dx.

         call dcopy (n,x,1,y,1)
         call daxpy (n,h,dx,1,y,1)
c
         call f06dbf(ncnln,(1),needc,1)
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
         if (msglvl.gt.0) then
            if (emax.le.oktol) then
               write (rec,fmt=99998)
               call x04bay(iprint,2,rec)
            else
               write (rec,fmt=99997)
               call x04bay(iprint,2,rec)
            end if
            write (rec,fmt=99996) emax, imax
            call x04bay(iprint,3,rec)
         end if
         if (emax.ge.point9) inform = 1
      end if

c     element-wise check.

      if (lvrfyc.ge.2) then
         if (lvlder.eq.3) then

c           recompute the jacobian to find the non-constant elements.

            call f06qhf('g',ncnln,n,rdummy,rdummy,cjacu,ldcju)

            call f06dbf(ncnln,(1),needc,1)
            nstate = 0
            mode = 2

            call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,nstate,
     *                  iuser,user)
            if (mode.lt.0) go to 160

         end if

         call f06dbf(ncnln,(0),needc,1)
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

c        loop over each column.

         do 140 j = j3, j4
            call sload (ncnln,zero,err,1)
            ncolj = 0
            headng = .true.
            xj = x(j)

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

c                 check this jacobian element.

                  ncheck = ncheck + 1
                  ncolj = ncolj + 1
                  needc(i) = 1
c
                  cij = cjac(i,j)
                  cjsize = abs(cij)

c                 find a finite-difference interval by iteration.

                  iter = 0
                  hopt = two*(one+abs(xj))*sqrt(epsrf)
                  h = ten*hopt*signh
                  cdest = zero
                  sdest = zero
                  first = .true.

  100             x(j) = xj + h

                  call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 160
                  f1 = c1(i)

                  x(j) = xj + h + h
                  call confun(mode,ncnln,n,ldcju,needc,x,c1,cjacu,
     *                        nstate,iuser,user)
                  if (mode.lt.0) go to 160
                  f2 = c1(i)

                  call chcore(debug,done,first,epsaci,epsrf,c(i),info,
     *                        iter,itmax,cdest,fdest,sdest,errbnd,f1,f2,
     *                        h,hopt,hphi)

                  if ( .not. done) go to 100

c                 exit for this element.

                  cjdiff = cdest
                  err(i) = abs(cjdiff-cij)/(cjsize+one)

                  ok = err(i) .le. oktol
                  if (ok) then
                     key = lgood
                     ngood = ngood + 1
                  else
                     key = lbad
                     nwrong = nwrong + 1
                  end if
c
                  if (msglvl.gt.0) then
                     const = ok .and. info .eq. 1 .and. abs(cij)
     *                       .lt. epspt8
                     if ( .not. const) then
                        if (headng) then
                           write (rec,fmt=99994)
                           call x04bay(iprint,4,rec)
                           if (ok) then
                              write (rec,fmt=99993) j, xj, hopt, i, cij,
     *                          cjdiff, key, iter
                           else
                              write (rec,fmt=99992) j, xj, hopt, i, cij,
     *                          cjdiff, key, iter, result(info)
                           end if
                           call x04baf(iprint,rec(1))
                           headng = .false.
                        else
                           if (ok) then
                              write (rec,fmt=99991) hopt, i, cij,
     *                          cjdiff, key, iter
                           else
                              write (rec,fmt=99990) hopt, i, cij,
     *                          cjdiff, key, iter, result(info)
                           end if
                           call x04baf(iprint,rec(1))
                        end if
                     end if
                  end if
                  needc(i) = 0
               end if
  120       continue

            if (ncolj.gt.0) then
               imax = idamax(ncnln,err,1)
               emax = abs(err(imax))

               if (emax.ge.colmax) then
                  irow = imax
                  jcol = j
                  colmax = emax
               end if
            end if
            x(j) = xj

  140    continue

         inform = 0
         if (colmax.ge.point9) inform = 1

         if (msglvl.gt.0) then
            if (ncheck.eq.0) then
               write (rec,fmt=99986) ncset
               call x04baf(iprint,rec(1))
            else
               if (nwrong.eq.0) then
                  write (rec,fmt=99989) ngood, ncheck, j3, j4
               else
                  write (rec,fmt=99988) nwrong, ncheck, j3, j4
               end if
               call x04bay(iprint,3,rec)
               write (rec,fmt=99987) colmax, irow, jcol
               call x04bay(iprint,3,rec)
            end if
         end if

      end if

c     copy  ( constants + gradients + dummy values )  back into cjacu.

      call f06qff('g',ncnln,n,cjac,ldcj,cjacu,ldcju)

      return

  160 inform = mode

c     end of  e04xaw. (chcjac)

99999 format (//' verification of the constraint gradients.',/' ',
     *       '--')
99998 format (/' the constraint jacobian seems to be ok.')
99997 format (/' xxx  the constraint jacobian seems to be incorrect.')
99996 format (/' the largest relative error was',1p,d12.2,'  in constr',
     *       'aint',i5,/)
99995 format (/' every column contains a constant or missing element.')
99994 format (//' column    x(j)     dx(j)    row    jacobian value   ',
     *       '   difference approxn  itns',/)
99993 format (i7,1p,2d10.2,i5,1p,2d18.8,2x,a4,i6)
99992 format (i7,1p,2d10.2,i5,1p,2d18.8,2x,a4,i6,2x,a18)
99991 format (17x,1p,d10.2,i5,1p,2d18.8,2x,a4,i6)
99990 format (17x,1p,d10.2,i5,1p,2d18.8,2x,a4,i6,2x,a18)
99989 format (/i7,'  constraint jacobian elements out of the',i6,/9x,
     *       'set in cols',i6,'  through',i6,'  seem to be ok.')
99988 format (/' xxx  there seem to be',i6,'  incorrect jacobian eleme',
     *       'nts out of the',i6,/8x,'set in cols',i6,'  through',i6)
99987 format (/' the largest relative error was',1p,d12.2,'  in row',i5,
     *       ',  column',i5,/)
99986 format (' all',i6,'   assigned jacobian elements are constant.')
      end

      subroutine e04ncq(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,ap,res,
     *                  hz,p,gq,cq,r,zy,work)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      implicit none

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  ctp, pnorm
      integer           lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf
      logical           linobj, singlr, unitgz, unitq
c     .. array arguments ..
      double precision  a(lda,*), ap(*), cq(*), gq(n), hz(*), p(n),
     *                  r(ldr,*), res(*), work(n), zy(ldzy,*)
      integer           kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      double precision  gtp
c     .. external functions ..
      double precision  ddot, dnrm2
      external          ddot, dnrm2
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
c----------------------------------------------------------------------
      if (singlr) then

c        the triangular factor for the current objective function is
c        singular,  i.e., the objective is linear along the last column
c        of z1.  this can only occur when unitgz is true.

         if (nrz.gt.1) then
            call dcopy (nrz-1,r(1,nrz),1,p,1)
            call dtrsv ('u','n','n',nrz-1,r,ldr,p,1)
         end if
         p(nrz) = -one
c
         gtp = ddot(nrz,gq,1,p,1)
         if (gtp.gt.zero) call dscal (nrz,(-one),p,1)
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
               call dcopy (nrz,gq,1,hz,1)
               call dscal (nrz,(-one),hz,1)
               call dtrsv ('u','t','n',nrz,r,ldr,hz,1)
            end if
         else
            call dcopy (nrz,res,1,hz,1)
         end if
c
c        solve  rz1*pz1 = hz1.
c
         call dcopy (nrz,hz,1,p,1)
         call dtrsv ('u','n','n',nrz,r,ldr,p,1)
c
      end if
c
c     compute  p = z1*pz1  and its norm.
c
      if (linobj) ctp = ddot(nrz,cq,1,p,1)
      pnorm = dnrm2(nrz,p,1)
c
      call e04nbw(1,n,nrz,nfree,ldzy,unitq,kx,p,zy,work)

c     compute  ap.

      if (nclin.gt.0) then
         call dgemv ('no transpose',nclin,n,one,a,lda,p,1,zero,ap,1)
      end if

c     end of  e04ncq. (lsgetp)

99999 format (/' //e04ncq//   p ... ')
99998 format (/' //e04ncq//  ap ... ')
99997 format (1p,5d15.5)
      end

      subroutine e04ncp(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,res,featol,gq,cq,r,
     *                  x,wtinf,zy,wrk)
c----------------------------------------------------------------------
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
c                    ( q(free)  0       ),
c                    (   0      i(fixed))
c     where  q(free)  is the orthogonal factor of  a(free)  and  a  is
c     the matrix of constraints in the working set.  the transformed
c     gradients are stored in gq.

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, suminf, tolrnk
      integer           lda, ldr, ldzy, n, nclin, nfree, nrank, nrz,
     *                  numinf, nz
      logical           linobj, singlr, unitgz, unitq
      character*2       prbtyp
c     .. array arguments ..
      double precision  a(lda,*), bl(*), bu(*), cq(*), featol(*), gq(n),
     *                  r(ldr,*), res(*), wrk(n), wtinf(*), x(n),
     *                  zy(ldzy,*)
      integer           istate(*), kx(n)
c     .. local scalars ..
      double precision  biglow, bigupp, ctx, feasj, rownrm, s, weight
      integer           j, k
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           f06klf
      external          ddot, dnrm2, f06klf
c----------------------------------------------------------------------
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
         call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,gq,zy,wrk)
         unitgz = .true.
      else if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         call sload (n,zero,gq,1)
      else
c
c        ready for the optimality phase.
c        set nrz so that rz1 is nonsingular.
c
         if (nrank.eq.0) then
            if (linobj) then
               call dcopy (n,cq,1,gq,1)
            else
               call sload (n,zero,gq,1)
            end if
            nrz = 0
         else

c           compute gq = - r' * (transformed residual)

            call f06fdf(nrank,-one,res,1,gq,1)
            call dtrmv ('u','t','n',nrank,r,ldr,gq,1)
            if (nrank.lt.n) call dgemv ('t',nrank,n-nrank,-one,
     *                                 r(1,nrank+1),ldr,res,1,zero,
     *                                 gq(nrank+1),1)
            if (linobj) call daxpy (n,one,cq,1,gq,1)
            unitgz = .false.
            rownrm = dnrm2(n,r(1,1),ldr)
            if (rownrm.le.tolrnk .or. abs(r(1,1)).le.rownrm*tolrnk) then
               nrz = 0
            else
               nrz = f06klf(min(nrank,nz),r,ldr+1,tolrnk)
            end if
         end if
         singlr = .false.
      end if
c
      return
c
c     end of  e04ncp. (lsgset)
c
      end

      subroutine e04mfn(subr,msg,v,lenv)
c----------------------------------------------------------------------
c     e04mfn  prints the array v in debug format.
c----------------------------------------------------------------------
c     .. scalar arguments ..
      integer           lenv
      character*6       subr
      character*(*)     msg
c     .. array arguments ..
      double precision  v(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      integer           i, ii
c     .. local arrays ..
      character*80      rec(3)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
c----------------------------------------------------------------------

      if (lenv.le.0) then
         write (rec,fmt=99999) subr, msg
         call x04bay(iprint,2,rec)
      else
         write (rec,fmt=99998) subr, msg
         call x04bay(iprint,2,rec)
         do 20 i = 1, lenv, 5
            write (rec,fmt=99997) (v(ii),ii=i,min(i+4,lenv))
            call x04baf(iprint,rec(1))
   20    continue
      end if

c     end of  e04mfn. (cmmsg1)

99999 format (/' //',a6,'//  ',a)
99998 format (/' //',a6,'//  ',a,' ... ')
99997 format (1p,5d15.5)
      end

      subroutine e04ncr(n,nclin,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *                  ax,bl,bu,featol,x,work)

c     e04ncr  computes the following...
c     (1)  the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, cvnorm, errmax
      integer           jmax, n, nclin, nviol
c     .. array arguments ..
      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), work(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  biglow, bigupp, con, feasj, res, tolj
      integer           i, is, j

      double precision  dnrm2
      integer           idamax
      external          dnrm2, idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd

c     compute nviol,  the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint violations and
c     residuals of the constraints in the working set.

      nviol = 0

      do 40 j = 1, n + nclin
         feasj = featol(j)
         is = istate(j)
         res = zero

         if (is.ge.0 .and. is.lt.4) then
            if (j.le.n) then
               con = x(j)
            else
               i = j - n
               con = ax(i)
            end if

            tolj = feasj

c           check for constraint violations.

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

c           this constraint is satisfied,  but count the residual as a
c           violation if the constraint is in the working set.

            if (is.le.0) res = zero
            if (is.eq.1) res = bl(j) - con
            if (is.ge.2) res = bu(j) - con
            if (abs(res).gt.feasj) nviol = nviol + 1
         end if
   20    work(j) = res
   40 continue
c
      jmax = idamax(n+nclin,work,1)
      errmax = abs(work(jmax))

      cvnorm = dnrm2(n+nclin,work,1)

c     end of  e04ncr. (lsfeas)

99999 format (/' //e04ncr//  the maximum violation is ',1p,d14.2,' in ',
     *       'constraint',i5)
      end

      subroutine e04mfl(msglvl,n,nrz,nz,zerolm,notopt,numinf,trusml,
     *                  smllst,jsmlst,tinyst,jtiny,gq)
c----------------------------------------------------------------------
c     e04mfl  updates jsmlst and smllst when there are artificial
c     constraints.
c
c     on input,  jsmlst  is the index of the minimum of the set of
c     adjusted multipliers.
c     on output, a negative jsmlst defines the index in q'g of the
c     artificial constraint to be deleted.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
c     .. scalar arguments ..
      double precision  smllst, tinyst, trusml, zerolm
      integer           jsmlst, jtiny, msglvl, n, notopt, nrz, numinf,
     *                  nz
c     .. array arguments ..
      double precision  gq(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)
c     .. local scalars ..
      double precision  rlam
      integer           j, k, kk, length
c     .. local arrays ..
      character*80      rec(3)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
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
c
      if (msglvl.ge.20) then
         if (isumm.ge.0) then
            write (rec,fmt=99999)
            call x04bay(isumm,2,rec)
            length = nz - nrz
            do 40 k = 1, length, 4
               write (rec,fmt=99998) (gq(kk),kk=k,min(k+3,length))
               call x04baf(isumm,rec(1))
   40       continue
         end if
      end if

c     end of  e04mfl.  (cmmul2)

99999 format (/' multipliers for the artificial constraints        ')
99998 format (4(5x,1p,d11.2))
99997 format (/' //e04mfl//  jsmlst     smllst',/' //e04mfl//  ',i6,1p,
     *       d11.2)
      end

      subroutine e04uct(ktcond,convrg,lsumry,msgnp,msgqp,ldr,ldt,n,
     *                  nclin,ncnln,nctotl,nactiv,linact,nlnact,nz,
     *                  nfree,majit0,majits,minits,istate,alfa,nfun,
     *                  condhz,condh,condt,objalf,objf,gfnorm,gznorm,
     *                  cvnorm,ax,c,r,t,violn,x,work)
c----------------------------------------------------------------------
c     e04uct  prints various levels of output for e04ucz and e04upz.
c
c           msg        cumulative result
c                   --
c
c        le   0        no output.
c
c        eq   1        nothing (but full output later).
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
c     .. scalar arguments ..
      double precision  alfa, condh, condhz, condt, cvnorm, gfnorm,
     *                  gznorm, objalf, objf
      integer           ldr, ldt, linact, majit0, majits, minits, msgnp,
     *                  msgqp, n, nactiv, nclin, ncnln, nctotl, nfree,
     *                  nfun, nlnact, nz
      logical           convrg
      character*5       lsumry
c     .. array arguments ..
      double precision  ax(*), c(*), r(ldr,*), t(ldt,*), violn(*),
     *                  work(n), x(n)
      integer           istate(nctotl)
      logical           ktcond(2)
c     .. scalars in common ..
      double precision  rhodmp, rhomax, rhonrm, scale
      integer           iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  cviols
      integer           i, inct, j, k, mjr, mnr, ndf, neval
      logical           first, newset, nlncon, prthdr
c     .. local arrays ..
      character*132     rec(4)
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
c
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
         nlncon = ncnln .gt. 0
         first = majits .eq. majit0
c

c        if necessary, print a header.
c        print a single line of information.

         if (isumm.ge.0) then

c           terse line for the monitoring file.

            newset = lines1 .ge. 50000
            prthdr = msgqp .gt. 0 .or. first .or. msgnp .ge. 20 .or.
     *               newset

            if (prthdr) then
               if (nlncon) then
                  write (rec,fmt=99998)
                  call x04bay(isumm,3,rec)
               else
                  write (rec,fmt=99996)
                  call x04bay(isumm,3,rec)
               end if
               lines1 = 0
            end if

            if (nlncon) then
               write (rec,fmt=99997) mjr, mnr, alfa, neval, objalf,
     *           cvnorm, gznorm, ndf, n - nfree, linact, nlnact,
     *           scale*rhonrm, gfnorm, condh, condhz, condt, convrg,
     *           ktcond(1), ktcond(2), lsumry
               call x04baf(isumm,rec(1))
            else
               write (rec,fmt=99995) mjr, mnr, alfa, neval, objalf,
     *           gznorm, ndf, n - nfree, linact, gfnorm, condh, condhz,
     *           condt, convrg, ktcond(1), ktcond(2), lsumry
               call x04baf(isumm,rec(1))
            end if
            lines1 = lines1 + 1
         end if

         if (iprint.ge.0 .and. isumm.ne.iprint) then
c           
c           terse line for the print file.
c           
            newset = lines2 .ge. 50000
            prthdr = msgqp .gt. 0 .or. first .or. newset

            if (prthdr) then
               if (nlncon) then
                  write (rec,fmt=99994)
                  call x04bay(iprint,3,rec)
               else
                  write (rec,fmt=99992)
                  call x04bay(iprint,3,rec)
               end if
               lines2 = 0
            end if

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

         if (msgnp.ge.20) then
            if (isumm.ge.0) then
               if (ncnln.eq.0) then
                  write (rec,fmt=99990) objf
                  call x04bay(isumm,2,rec)
               else
                  cviols = dnrm2(ncnln,violn,1)
                  write (rec,fmt=99989) objf, cviols
                  call x04bay(isumm,2,rec)
               end if

c              print the constraint values.

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
               if (ncnln.gt.0) then
                  write (rec,fmt=99985)
                  call x04bay(isumm,2,rec)
                  do 60 i = 1, ncnln, 5
                     write (rec,fmt=99982) (c(k),istate(n+nclin+k),k=i,
     *                 min(i+4,ncnln))
                     call x04baf(isumm,rec(1))
   60             continue
               end if

               if (msgnp.ge.30) then

c                 print the diagonals of  t  and  r.

                  inct = ldt - 1
                  if (nactiv.gt.0) then
                     call dcopy (nactiv,t(nactiv,nz+1),inct,work,1)
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
     *       /' -')
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

      subroutine e04nfm(side,n,k1,k2,s,a,lda)
c----------------------------------------------------------------------
c  e04nfm applies a  sequence  of  pairwise interchanges to either  the
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
c     .. scalar arguments ..
      integer           k1, k2, lda, n
      character*1       side
c     .. array arguments ..
      double precision  a(lda,*), s(*)
c     .. local scalars ..
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

c     end of e04nfm.

      end

      subroutine e04ucj(first,debug,done,imprvd,inform,maxf,numf,nout,
     *                  alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,
     *                  tolrel,toltny,alfa,alfbst,fbest)
c----------------------------------------------------------------------
c     e04ucj  finds a sequence of improving estimates of a minimizer of
c     the univariate function f(alpha) in the interval (0,alfmax].
c     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
c     e04ucj  requires  f(alpha) (but not f'(alpha)) to be evaluated
c     in the interval.  new estimates of a minimizer are computed using
c     safeguarded quadratic interpolation.

c     reverse communication is used to allow the calling program to
c     evaluate f.  some of the parameters must be set or tested by the
c     calling program.  the remainder would ordinarily be local
c     variables.

c     input parameters (relevant to the calling program)

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
c     
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

      double precision  zero, point1, half
      parameter         (zero=0.0d+0,point1=0.1d+0,half=0.5d+0)
      double precision  one, two, five
      parameter         (one=1.0d+0,two=2.0d+0,five=5.0d+0)
      double precision  ten, eleven
      parameter         (ten=1.0d+1,eleven=1.1d+1)
c     .. scalar arguments ..
      double precision  alfa, alfbst, alfmax, alfsml, epsaf, fbest,
     *                  ftry, g0, targtg, tolabs, tolrel, toltny
      integer           inform, maxf, nout, numf
      logical           debug, done, first, imprvd
c     .. local scalars ..
      double precision  a, artifa, artifb, b, daux, dtry, endpnt, fa,
     *                  factor, fv, fw, gv, gw, q, s, tol, tolmax,
     *                  truea, trueb, xmidpt, xtry, xv, xw
      integer           nsamea, nsameb
      logical           badfun, braktd, closef, crampd, extrap, found,
     *                  moved, quitf, quitfz, quiti, quits, setxv, vset,
     *                  wset, xinxw
c     .. save statement ..
      save              braktd, crampd, extrap, moved, vset, wset,
     *                  nsamea, nsameb, a, b, fa, factor, xtry, xw, fw,
     *                  xv, fv, tolmax
c----------------------------------------------------------------------
      imprvd = .false.
      badfun = .false.
      quitf = .false.
      quitfz = .false.
      quits = .false.
      quiti = .false.

      if (first) then

c        first entry.  initialize various quantities, check input data
c        and prepare to evaluate the function at the initial step alfa.

         first = .false.
         numf = 0
         alfbst = zero
         badfun = alfmax .le. toltny .or. g0 .ge. zero
         done = badfun
         moved = .false.

         if ( .not. done) then
            braktd = .false.
            crampd = alfmax .le. tolabs
            extrap = .false.
            vset = .false.
            wset = .false.
            nsamea = 0
            nsameb = 0

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

c        check if xtry is in the interval (xw,0) or (0,xw).

         if (wset) then
            xinxw = zero .lt. xtry .and. xtry .le. xw .or. xw .le.
     *              xtry .and. xtry .lt. zero
         else
            xinxw = .false.
         end if

         imprvd = ftry .lt. fbest
         if (vset) then
            closef = abs(fbest-fv) .le. epsaf
         else
            closef = .false.
         end if

         if (imprvd) then

c           we seem to have an improvement.  the new point becomes the
c           origin and other points are shifted accordingly.

            if (wset) then
               xv = xw - xtry
               fv = fw
               vset = .true.
            end if

            xw = zero - xtry
            fw = fbest
            wset = .true.
            fbest = ftry
            alfbst = alfa
            moved = .true.

            a = a - xtry
            b = b - xtry
            extrap = .not. xinxw

c           decrease the length of (a,b).

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

c           quit if there has been no progress and ftry, fbest, fw
c           and fv are all within epsaf of each other.

            quitfz = .true.
         else

c           the new function value is no better than the current best
c           point.  xtry must an end point of the new (a,b).

            if (xtry.lt.zero) then
               a = xtry
               fa = ftry
               nsamea = 0
            else
               b = xtry
               nsameb = 0
               braktd = .true.
            end if

c           the origin remains unchanged but xtry may qualify as xw.

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

c        check the termination criteria.

         tol = tolabs + tolrel*alfbst
         truea = alfbst + a
         trueb = alfbst + b

         found = moved .and. abs(fa-fbest) .le. -a*targtg
         quitf = numf .ge. maxf
         quiti = b - a .le. tol + tol
         quits = trueb .le. alfsml

         if (quiti .and. .not. moved) then

c           the interval of uncertainty appears to be small enough,
c           but no better point has been found.  check that changing
c           alfa by b-a changes f by less than epsaf.

            tol = tol/ten
            tolabs = tol
            quiti = abs(fw) .le. epsaf .or. tol .le. toltny
         end if

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

c           fit a parabola.
c           see if there are two or three points for the parabolic fit.

            gw = (fw-fbest)/xw
            if (vset .and. moved) then

c              three points available.  use fbest, fw and fv.

               gv = (fv-fbest)/xv
               s = gv - (xv/xw)*gw
               q = two*(gv-gw)

            else

c              only two points available.  use fbest, fw and g0.

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

c              a minimizer has not yet been bracketed.
c              set an artificial upper bound by expanding the interval
c              xw  by a suitable factor.

               xtry = -factor*xw
               artifb = xtry
               if (alfbst+xtry.lt.alfmax) factor = five*factor
            else if (vset .and. moved) then

c              three points exist in the interval of uncertainty.
c              check if the points are configured for an extrapolation
c              or an interpolation.

               if (extrap) then

c                 the points are configured for an extrapolation.

                  if (xw.lt.zero) endpnt = b
                  if (xw.gt.zero) endpnt = a
               else

c                 if the interpolation appears to be overestimating the
c                 distance to a minimizer,  damp the interpolation step.

                  if (nsamea.ge.3 .or. nsameb.ge.3) then
                     factor = factor/five
                     s = factor*s
                  else
                     factor = one
                  end if

c                 the points are configured for an interpolation.  the
c                 artificial interval will be just (a,b).  set endpnt so
c                 that xtry lies in the larger of the intervals (a,b)
c                 and  (0,b).

                  if (xmidpt.gt.zero) then
                     endpnt = b
                  else
                     endpnt = a
                  end if

c                 if a bound has remained the same for three iterations,
c                 set endpnt so that  xtry  is likely to replace the
c                 offending bound.

                  if (nsamea.ge.3) endpnt = a
                  if (nsameb.ge.3) endpnt = b
               end if

c              compute the default value of  xtry.

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

               if (extrap) then
                  if (xtry.le.zero) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else

c              the gradient at the origin is being used for the
c              polynomial fit.  set the default xtry to one tenth xw.

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

c                 accept the polynomial fit.

                  if (abs(s*xw).ge.q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
               end if
            end if
         end if
      end if

      if ( .not. done) then
         alfa = alfbst + xtry
         if (braktd .or. alfa.lt.alfmax-tolmax) then

c           the function must not be evaluated too close to a or b.
c           (it has already been evaluated at both those points.)

            xmidpt = half*(a+b)
            if (xtry.le.a+tol .or. xtry.ge.b-tol) then
               if (xmidpt.le.zero) then
                  xtry = -tol
               else
                  xtry = tol
               end if
            end if

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

c     end of  e04ucj. (srchq)

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
99990 format (' -',/)
      end

      subroutine e04mfm(prbtyp,msglvl,n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,anorms,gq,rlamda,t,wtinf)
c----------------------------------------------------------------------
c     e04mfm  first computes the lagrange multiplier estimates for the
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)
c     .. scalar arguments ..
      double precision  biggst, smllst, tinyst, trubig, trusml, zerolm
      integer           jbigst, jinf, jsmlst, jtiny, kbigst, ksmlst,
     *                  lda, ldt, msglvl, n, nactiv, nfree, notopt,
     *                  numinf, nz
      character*2       prbtyp
c     .. array arguments ..
      double precision  a(lda,*), anorms(*), gq(n), rlamda(n), t(ldt,*),
     *                  wtinf(*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lcdbg
c     .. arrays in common ..
      integer           ilcdbg(ldbg)
c     .. local scalars ..
      double precision  anormj, blam, rlam, scdlam
      integer           i, is, j, k, kk, l, nfixed
c     .. local arrays ..
      character*80      rec(3)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ee04mf/ilcdbg, lcdbg
c----------------------------------------------------------------------
      nfixed = n - nfree

      jtiny = 0
      jsmlst = 0
      ksmlst = 0

      jbigst = 0
      kbigst = 0

c     compute  jsmlst  for regular constraints and temporary bounds.

c     first, compute the lagrange multipliers for the general
c     constraints in the working set, by solving  t'*lamda = y'g.

      if (n.gt.nz) call dcopy (n-nz,gq(nz+1),1,rlamda,1)
      if (nactiv.gt.0) call dtrsv ('u','t','n',nactiv,t(1,nz+1),ldt,
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

         i = j - n
         if (j.le.n) anormj = one
         if (j.gt.n) anormj = anorms(i)

         rlam = rlamda(k)

c        change the sign of the estimate if the constraint is in
c        the working set at its upper bound.

         if (is.eq.2) rlam = -rlam
         if (is.eq.3) rlam = abs(rlam)
         if (is.eq.4) rlam = -abs(rlam)

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
c
         scdlam = rlam/wtinf(j)
         if (scdlam.gt.biggst .and. j.gt.jinf) then
            biggst = scdlam
            trubig = rlamda(k)
            jbigst = j
            kbigst = k
         end if
   60 continue

c     if required, print the multipliers.

      if (msglvl.ge.20) then
         if (isumm.ge.0) then
            if (nfixed.gt.0) then
               write (rec,fmt=99999) prbtyp
               call x04bay(isumm,2,rec)
               do 80 k = 1, nfixed, 4
                  write (rec,fmt=99998) (kx(nfree+kk),rlamda(nactiv+kk),
     *              kk=k,min(k+3,nfixed))
                  call x04baf(isumm,rec(1))
   80          continue
            end if
            if (nactiv.gt.0) then
               write (rec,fmt=99997) prbtyp
               call x04bay(isumm,2,rec)
               do 100 k = 1, nactiv, 4
                  write (rec,fmt=99998) (kactiv(kk),rlamda(nactiv-kk+1),
     *              kk=k,min(k+3,nactiv))
                  call x04baf(isumm,rec(1))
  100          continue
            end if
         end if
      end if

c     end of  e04mfm.  (cmmul1)

99999 format (/' multipliers for the ',a2,' bound  constraints   ')
99998 format (4(i5,1p,d11.2))
99997 format (/' multipliers for the ',a2,' linear constraints   ')
99996 format (/' //e04mfm//  jsmlst     smllst     ksmlst (scaled) ',
     *       /' //e04mfm//  ',i6,1p,d11.2,5x,i6)
99995 format (' //e04mfm//  jbigst     biggst     kbigst (scaled) ',
     *       /' //e04mfm//  ',i6,1p,d11.2,5x,i6)
99994 format (' //e04mfm//   jtiny     tinyst                     ',
     *       /' //e04mfm//  ',i6,1p,d11.2)
      end

      subroutine e04mfx(nout,buffer,key)
c----------------------------------------------------------------------
c     e04mfx   decodes the option contained in  buffer  in order to set
c     a parameter value in the relevant element of  iprmlc  or  rprmlc.

c     input:

c     nout   a unit number for printing error messages.
c            nout  must be a valid unit.

c     output:

c     key    the first keyword contained in buffer.

c     e04mfx  calls e04udx and the subprograms
c             lookup, scannr, tokens (now called e04udy, e04udw, e04udv)
c-----------------------------------------------------------------------
      implicit none

      integer           mxparm
      parameter         (mxparm=30)
      integer           maxkey, maxtie, maxtok, maxtyp
      parameter         (maxkey=18,maxtie=3,maxtok=10,maxtyp=4)
      integer           idummy
      double precision  rdummy
      logical           sorted
      double precision  zero
      parameter         (idummy=-11111,rdummy=-11111.0d+0,sorted=.true.,
     *                  zero=0.0d+0)
c     .. scalar arguments ..
      integer           nout
      character*16      key
      character*(*)     buffer
c     .. scalars in common ..
      double precision  bigbnd, bigdx, bndlow, bndupp, tolact, tolfea,
     *                  tolrnk
      integer           idbglc, iprnt, isumry, itmax1, itmax2, kchk,
     *                  kcycle, lcrash, ldbglc, lprob, maxact, maxnz,
     *                  mm, msglc, mxfree, nn, nnclin, nprob
c     .. arrays in common ..
      double precision  rpadlc(23), rpsvlc(mxparm)
      integer           ipadlc(14), ipsvlc(mxparm)
c     .. local scalars ..
      double precision  rvalue
      integer           i, idbg, lenbuf, loc1, loc2, loc3, msgdbg,
     *                  msglvl, ntoken
      logical           first, more, number
      character*16      key2, key3, value
      character*80      rec
c     .. local arrays ..
      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*16      keys(maxkey), ties(maxtie), token(maxtok),
     *                  type(maxtyp)
c     .. external functions ..
      logical           e04udx
      external          e04udx

      common            /fe04mf/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /ge04mf/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)
c     .. save statement ..
      save              /fe04mf/, /ge04mf/, first
      data              first/.true./
      data              keys(1)/'begin           '/,
     *                  keys(2)/'check           '/,
     *                  keys(3)/'cold            '/,
     *                  keys(4)/'crash           '/,
     *                  keys(5)/'defaults        '/,
     *                  keys(6)/'end             '/,
     *                  keys(7)/'expand          '/,
     *                  keys(8)/'feasibility     '/,
     *                  keys(9)/'infinite        '/,
     *                  keys(10)/'iterations      '/,
     *                  keys(11)/'iters:iterations'/,
     *                  keys(12)/'itns :iterations'/,
     *                  keys(13)/'list            '/,
     *                  keys(14)/'monitoring      '/,
     *                  keys(15)/'nolist          '/,
     *                  keys(16)/'print           '/,
     *                  keys(17)/'problem         '/,
     *                  keys(18)/'warm            '/
      data              ties/'bound           ', 'step            ',
     *                  'type            '/
      data              type/'feasible     :fp', 'fp              ',
     *                  'linear       :lp', 'lp              '/
c----------------------------------------------------------------------
      if (first) then
         first = .false.
         do 20 i = 1, mxparm
            iprmlc(i) = idummy
            rprmlc(i) = rdummy
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

      ntoken = maxtok
      call e04udv(buffer(1:lenbuf),maxtok,ntoken,token)
      key = token(1)
      key2 = token(2)
      key3 = token(3)

c     certain keywords require no action.

      if (key.eq.' ' .or. key.eq.'begin') go to 80
      if (key.eq.'list' .or. key.eq.'nolist') go to 80
      if (key.eq.'end') go to 80

c     most keywords will have an associated integer or real value,
c     so look for it no matter what the keyword.

      i = 1
      number = .false.

   40 if (i.lt.ntoken .and. .not. number) then
         i = i + 1
         value = token(i)
         number = e04udx(value)
         go to 40
      end if

      if (number) then
         read (value,fmt='(bn, e16.0)') rvalue
      else
         rvalue = zero
      end if

c     convert the keywords to their most fundamental form
c     (upper case, no abbreviations).
c     sorted says whether the dictionaries are in alphabetic order.
c     loci   says where the keywords are in the dictionaries.
c     loci = 0 signals that the keyword wasn't there.
c     loci < 0 signals that the keyword is ambiguous.

      call e04udy(maxkey,keys,sorted,key,loc1)
      if (loc1.lt.0) then
         write (rec,fmt=99996) key
         call x04baf(nout,rec)
         return
      end if
      call e04udy(maxtie,ties,sorted,key2,loc2)

c     decide what to do about each keyword.
c     the second keyword (if any) might be needed to break ties.
c     some seemingly redundant testing of more is used
c     to avoid compiler limits on the number of consecutive else ifs.

      more = .true.
      if (more) then
         more = .false.
         if (key.eq.'check       ') then
            kchk = rvalue
         else if (key.eq.'cold        ') then
            lcrash = 0
         else if (key.eq.'crash       ') then
            tolact = rvalue
         else if (key.eq.'defaults    ') then
            do 60 i = 1, mxparm
               iprmlc(i) = idummy
               rprmlc(i) = rdummy
   60       continue
         else if (key.eq.'expand      ') then
            kcycle = rvalue
         else if (key.eq.'feasibility ') then
            tolfea = rvalue
         else
            more = .true.
         end if
      end if

      if (more) then
         more = .false.
         if (key.eq.'infinite    ') then
            if (key2.eq.'bound       ') bigbnd = rvalue*0.99999d+0
            if (key2.eq.'step        ') bigdx = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'iterations  ') then
            itmax2 = rvalue
         else if (key.eq.'monitoring  ') then
            isumry = rvalue
         else if (key.eq.'print       ') then
            msglvl = rvalue
         else if (key.eq.'problem     ') then
            if (key2.eq.'type  ') then

c             recognize     problem type = lp     etc.

               call e04udy(maxtyp,type,sorted,key3,loc3)
               if (key3.eq.'fp') lprob = 1
               if (key3.eq.'lp') lprob = 2
               if (loc3.eq.0) then
                  write (rec,fmt=99997) key3
                  call x04baf(nout,rec)
                  lprob = 10
               end if
            else
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'warm        ') then
            lcrash = 1
         else
            write (rec,fmt=99999) key
            call x04baf(nout,rec)
         end if
      end if

   80 return

c     end of  e04mfx.  (lpkey)

99999 format (' xxx  keyword not recognized:         ',a)
99998 format (' xxx  second keyword not recognized:  ',a)
99997 format (' xxx  third  keyword not recognized:  ',a)
99996 format (' xxx  ambiguous keyword:              ',a)
      end

      subroutine e04xax(inform,msglvl,n,bigbnd,epsrf,oktol,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,y,iuser,user)
c----------------------------------------------------------------------
c     e04xax  checks if the gradients of the objective function have
c     been coded correctly.

c     on input,  the value of the objective function at the point x is
c     stored in objf.  the corresponding gradient is stored in gradu.
c     if any gradient element has not been specified,  it will have a
c     dummy value.  missing values are not checked.

c     a cheap test is first undertaken by calculating the directional
c     derivative using two different methods. if this proves
c     satisfactory and no further information is desired, e04xax is
c     terminated. otherwise, the routine chcore is called to give
c     optimal step-sizes and a forward-difference approximation to
c     each element of the gradient for which a test is deemed
c     necessary, either by the program or the user.

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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  zero, half, point9
      parameter         (zero=0.0d+0,half=0.5d+0,point9=0.9d+0)
      double precision  one, two, ten
      parameter         (one=1.0d+0,two=2.0d+0,ten=1.0d+1)
      character*4       lbad, lgood
      parameter         (lbad='bad?',lgood='  ok')
c     .. scalar arguments ..
      double precision  bigbnd, epsrf, fdchk, objf, oktol, xnorm
      integer           inform, msglvl, n
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  biglow, bigupp, cdest, dxj, dxmult, emax, epsa,
     *                  errbnd, error, f1, f2, fdest, gdiff, gdx, gj,
     *                  gsize, h, hopt, hphi, objf1, sdest, stepbl,
     *                  stepbu, xj
      integer           info, iter, itmax, j, j1, j2, jmax, mode,
     *                  ncheck, ngood, nstate, nwrong
      logical           const, debug, done, first, headng, needed, ok
      character*4       key
c     .. local arrays ..
      character*18      result(0:4)
      character*120     rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy
      common            /fe04uc/inpdbg, npdbg
c     .. data statements ..
      data              result/'                 ', 'constant?      ',
     *                  'linear or odd?   ', 'too nonlinear?',
     *                  'small derivative?'/
c----------------------------------------------------------------------
      inform = 0

      needed = lvrfyc .eq. 0 .or. lvrfyc .eq. 1 .or. lvrfyc .eq. 3
      if ( .not. needed) return

      if (msglvl.gt.0) then
         write (rec,fmt=99999)
         call x04bay(iprint,4,rec)
      end if

      nstate = 0

      biglow = -bigbnd
      bigupp = bigbnd

c     perform the cheap test.

      h = (one+xnorm)*fdchk

      if (n.le.100) then
         dxmult = 0.9
      else if (n.le.250) then
         dxmult = 0.99
      else
         dxmult = 0.999
      end if

      dxj = one/n
      do 20 j = 1, n
         dx(j) = dxj
         dxj = -dxj*dxmult
   20 continue

c     do not perturb x(j) if the  j-th  element is missing.
c     compute the directional derivative.

      ncheck = 0

      do 40 j = 1, n

         if (grad(j).eq.rdummy) then
            dx(j) = zero
         else

            ncheck = ncheck + 1

            xj = x(j)
            stepbl = -one
            stepbu = one

            if (bl(j).gt.biglow) stepbl = max(stepbl,bl(j)-xj)

            if (bu(j).lt.bigupp .and. bu(j).gt.bl(j))
     *          stepbu = min(stepbu,bu(j)-xj)

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
      gdx = ddot(n,gradu,1,dx,1)

c     make forward-difference approximation along  p.

      call dcopy (n,x,1,y,1)
      call daxpy (n,h,dx,1,y,1)

      mode = 0
      call objfun (mode,n,y,objf1,gradu,nstate,iuser,user)
      if (mode.lt.0) go to 100

      gdiff = (objf1-objf)/h
      error = abs(gdiff-gdx)/(abs(gdx)+one)

      ok = error .le. oktol

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

               stepbl = biglow
               stepbu = bigupp
               if (bl(j).gt.biglow) stepbl = bl(j) - xj
               if (bu(j).lt.bigupp) stepbu = bu(j) - xj

               hopt = two*(one+abs(xj))*sqrt(epsrf)
               h = ten*hopt
               if (half*(stepbl+stepbu).lt.zero) h = -h

   60          x(j) = xj + h
               call objfun (mode,n,x,f1,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100

               x(j) = xj + h + h
               call objfun(mode,n,x,f2,gradu,nstate,iuser,user)
               if (mode.lt.0) go to 100

               call chcore (debug,done,first,epsa,epsrf,objf,info,iter,
     *                     itmax,cdest,fdest,sdest,errbnd,f1,f2,h,hopt,
     *                     hphi)

               if ( .not. done) go to 60

c              exit for this variable.

               gdiff = cdest
               x(j) = xj

               error = abs(gdiff-gj)/(gsize+one)
               if (error.ge.emax) then
                  emax = error
                  jmax = j
               end if

               ok = error .le. oktol
               if (ok) then
                  key = lgood
                  ngood = ngood + 1
               else
                  key = lbad
                  nwrong = nwrong + 1
               end if

               if (msglvl.gt.0) then

c                 zero elements are not printed.

                  const = ok .and. info .eq. 1 .and. abs(gj) .lt. epspt8
                  if ( .not. const) then
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
      call dcopy (n,grad,1,gradu,1)
c
      return
c
  100 inform = mode

c     end of  e04xax. (chkgrd)

99999 format (//' verification of the objective gradients.',/' -',
     *       '')
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

      subroutine e04ucu(feasqp,unitq,nqperr,majits,minits,n,nclin,ncnln,
     *                  ldcj,ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,
     *                  numinf,istate,kactiv,kx,dxnorm,gdx,qpcurv,aqp,
     *                  adx,anorm,ax,bl,bu,c,cjac,clamda,cmul,cs,dlam,
     *                  dslk,dx,qpbl,qpbu,qptol,r,rho,slk,violn,x,wtinf,
     *                  w)
c----------------------------------------------------------------------
c     e04ucu   does the following:
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
c     .. scalar arguments ..
      double precision  dxnorm, gdx, qpcurv
      integer           ldaqp, ldcj, ldr, linact, majits, minits, n,
     *                  nactiv, nclin, ncnln, nfree, nlnact, nqperr,
     *                  numinf, nz
      logical           feasqp, unitq
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  amin, biglow, bigupp, blj, buj, con, condmx,
     *                  quotnt, ssq, ssq1, suminf, viol, weight, wscale,
     *                  wtmax, wtmin
      integer           i, idbg, inform, iswap, j, jinf, k, k1,
     *                  k2, kviol, l, lgq, lhpq, lrlam, lrpq, lrpq0, lt,
     *                  lwrk1, lzy, mjrdbg, mnrdbg, msgqp, nartif, ncqp,
     *                  nctotl, ngq, nmajor, nminor, nplin, nrank,
     *                  nrejtd, nrpq, ntry, nviol, nz1
      logical           linobj, overfl
c     .. local arrays ..
      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*8       names(1)
      character*80      rec(3)
c     .. external functions ..
      double precision  ddot, dnrm2, adivb
      external          ddot, dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common            /be04nb/lennam, ldt, ncolt, ldzy
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
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c     .. save statement ..
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/
c----------------------------------------------------------------------

      lrpq = locls(5)
      lrpq0 = locls(6)
      lhpq = locls(8)
      lgq = locls(9)
      lrlam = locls(10)
      lt = locls(11)
      lzy = locls(12)
      lwrk1 = locls(14)

      nrpq = 0
      ngq = 1

      feasqp = .true.
      linobj = .true.

      biglow = -bigbnd
      bigupp = bigbnd
      ssq1 = zero

      nplin = n + nclin
      nctotl = nplin + ncnln
      ncqp = nclin + ncnln
      nrank = n
      nrejtd = 0

      if (msgqp.gt.0) then
         write (rec,fmt=99999) majits
         call x04bay(iprint,3,rec)
      end if

c     generate the upper and lower bounds upon the search direction, the
c     weights on the sum of infeasibilities and the nonlinear constraint
c     violations.

      wscale = -one
      do 40 j = 1, nctotl

         if (j.le.n) then
            con = x(j)
         else if (j.le.nplin) then
            con = ax(j-n)
         else
            con = c(j-nplin)
         end if

         blj = bl(j)
         buj = bu(j)
         if (blj.gt.biglow) blj = blj - con
         if (buj.lt.bigupp) buj = buj - con

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

c           set the vector of nonlinear constraint violations.

   20       violn(i) = viol
         end if

         wtinf(j) = weight
         qpbl(j) = blj
         qpbu(j) = buj

   40 continue

      if (wscale.gt.zero) then
         wscale = one/wscale
         call dscal (nctotl,(wscale),wtinf,1)
      end if
c
      call f06flf(nctotl,wtinf,1,wtmax,wtmin)
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
         call f06qff('g',ncnln,n,cjac,ldcj,aqp(nclin+1,1),ldaqp)
c
         do 80 j = nclin + 1, ncqp
            anorm(j) = dnrm2(n,aqp(j,1),ldaqp)
   80    continue

c        count the number of linear constraints in the working set and
c        move them to the front of kactiv.  compute the norm of the
c        matrix of constraints in the working set.
c        let k1  point to the first nonlinear constraint.  constraints
c        with indices kactiv(k1),..., kactiv(nactiv)  must be
c        refactorized.

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
         if (nactiv.le.1) call f06flf(ncqp,anorm,1,asize,amin)

c        compute the absolute values of the nonlinear constraints in
c        the working set.  use dx as workspace.

         do 120 k = linact + 1, nactiv
            j = n + kactiv(k)
            if (istate(j).eq.1) dx(k) = abs(qpbl(j))
            if (istate(j).ge.2) dx(k) = abs(qpbu(j))
  120    continue

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
         if (k1.le.k2) call e04ncy(unitq,vertex,inform,k1,k2,nactiv,
     *                             nartif,nz,nfree,nrank,nrejtd,nrpq,
     *                             ngq,n,ldzy,ldaqp,ldr,ldt,istate,
     *                             kactiv,kx,condmx,aqp,r,w(lt),w(lrpq),
     *                             w(lgq),w(lzy),w(lwrk1),dx,w(lrlam),
     *                             msgqp)
      end if

c     solve for dx, the vector of minimum two-norm that satisfies the
c     constraints in the working set.

      call e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,ldaqp,
     *            ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,adx,qpbl,qpbu,
     *            w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt),w(lzy),w(lwrk1))

c     solve a quadratic program for the search direction  dx  and
c     multiplier estimates  clamda.

c     if there is no feasible point for the subproblem,  the sum of
c     infeasibilities is minimized subject to the linear constraints
c     (1  thru  jinf)  being satisfied.

      jinf = n + nclin

      ntry = 1

  180 call e04ncz('qp subproblem',qpnamd,names,linobj,unitq,nqperr,
     *            minits,jinf,ncqp,nctotl,nactiv,nfree,nrank,nz,nz1,n,
     *            ldaqp,ldr,istate,kactiv,kx,gdx,ssq,ssq1,suminf,numinf,
     *            dxnorm,qpbl,qpbu,aqp,clamda,adx,qptol,r,dx,w)

      nviol = 0
      if (numinf.gt.0) then

c           count the violated linear constraints.

         do 200 j = 1, nplin
            if (istate(j).lt.0) nviol = nviol + 1
  200    continue

         if (nviol.gt.0) then
            ntry = ntry + 1
            unitq = .true.
            nactiv = 0
            nfree = n
            nz = n
            call f06dbf(nctotl,(0),istate,1)

            call e04ucm(unitq,ncqp,nactiv,nfree,nz,n,nlnx,nctotl,ldzy,
     *                  ldaqp,ldr,ldt,istate,kactiv,kx,dxnorm,gdx,aqp,
     *                  adx,qpbl,qpbu,w(lrpq),w(lrpq0),dx,w(lgq),r,w(lt)
     *                  ,w(lzy),w(lwrk1))
         end if
      end if
      if ( .not. (nviol.eq.0 .or. ntry.gt.2)) go to 180

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
c
      linact = nactiv - nlnact

c     extract various useful quantities from the qp solution.

c     compute  hpq = r'r(pq)  from the transformed gradient of the qp
c     objective function and  r(pq)  from the transformed residual.
c
      call dscal (n,(-one),w(lrpq),1)
      call daxpy (n,(-one),w(lgq),1,w(lhpq),1)
      qpcurv = two*ssq
c
      if (ncnln.gt.0) then
         if (numinf.gt.0) then
            feasqp = .false.
            call sload (nctotl,(zero),clamda,1)

            if (nz.gt.0) then

c              compute a null space element for the search direction
c              as the solution of  z'hz(pz) = -z'g - z'hy(py).

c              overwrite dx with the transformed search direction
c              q'(dx).  the first nz elements of dx are zero.

               call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,dx,w(lzy),w(lwrk1)
     *                     )

c              overwrite the first nz elements of dx with the solution
c              of  (rz)u = -(v + w),  where  (rz)'w = z'g  and  v  is
c              vector of first nz elements of  r(pq).

               call dcopy (nz,w(lgq),1,dx,1)
               call dtrsv ('u','t','n',nz,r,ldr,dx,1)

               call daxpy (nz,(one),w(lrpq),1,dx,1)

               call dtrsv ('u','n','n',nz,r,ldr,dx,1)
               call dscal (nz,(-one),dx,1)

c              recompute rpq, hpq, gdx and qpcurv.

               call dcopy (nlnx,dx,1,w(lrpq),1)
               call dtrmv ('u','n','n',nlnx,r,ldr,w(lrpq),1)
               if (nlnx.lt.n) call dgemv ('n',nlnx,n-nlnx,one,
     *                       r(1,nlnx+1),ldr,dx(nlnx+1),1,one,w(lrpq),1)

               gdx = ddot(n,w(lgq),1,dx,1)
               qpcurv = ddot(n,w(lrpq),1,w(lrpq),1)

               call e04nbw(3,n,nz,nfree,ldzy,unitq,kx,dx,w(lzy),w(lwrk1)
     *                     )

c              recompute adx and the 2-norm of dx.

               dxnorm = dnrm2(n,dx,1)
               if (ncqp.gt.0) call dgemv ('n',ncqp,n,one,aqp,ldaqp,dx,1,
     *                                   zero,adx,1)

            end if

            call dcopy (nlnx,w(lrpq),1,w(lhpq),1)
            call dtrmv ('u','t','n',nlnx,r,ldr,w(lhpq),1)
            if (nlnx.lt.n) call dgemv ('t',nlnx,n-nlnx,one,r(1,nlnx+1),
     *                                ldr,w(lrpq),1,zero,w(lhpq+nlnx),1)
         end if


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
            quotnt = adivb(cmul(i),scale*rho(i),overfl)
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

c           compute the search direction for the slacks and multipliers.

            dslk(i) = adx(nclin+i) + cs(i)

            if (feasqp) then

c              if any constraint is such that  (dlam)*(c - s)  is
c              positive,  the merit function may be reduced immediately
c              by substituting the qp multiplier.

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

c     end of  e04ucu. (npiqp)

99999 format (/1x,79('-'),/' start of major itn',i6)
99998 format (/' //e04ucu // nqperr',/' //e04ucu // ',i6)
99997 format (/' //e04ucu // dx recomputed with null space portion...')
99996 format (/' //e04ucu // violations = ')
99995 format (/' //e04ucu // slacks     = ')
99994 format (1p,5d15.6)
99993 format (5g12.3)
      end

      subroutine e04ncj(prbtyp,isdel,iter,jadd,jdel,msglvl,nactiv,nfree,
     *                  n,nclin,nrank,ldr,ldt,nz,nrz,istate,alfa,condrz,
     *                  condt,gfnorm,gzrnrm,numinf,suminf,ctx,ssq,ax,r,
     *                  t,x,work)
c----------------------------------------------------------------------
c     e04ncj  prints various levels of output for  e04ncz.
c
c           msg        cumulative result
c                   --
c
c        le   0        no output.
c
c        eq   1        nothing (but full output later).
c
c        eq   5        one terse line of output.
c
c        ge  10        same as 5 (but full output later).
c
c        ge  20        constraint status,  x  and  ax.
c
c        ge  30        diagonals of  t  and  r.

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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      integer           mline1, mline2
      parameter         (mline1=50000,mline2=50000)
c     .. scalar arguments ..
      double precision  alfa, condrz, condt, ctx, gfnorm, gzrnrm, ssq,
     *                  suminf
      integer           isdel, iter, jadd, jdel, ldr, ldt, msglvl, n,
     *                  nactiv, nclin, nfree, nrank, nrz, numinf, nz
      character*2       prbtyp
c     .. array arguments ..
      double precision  ax(*), r(ldr,*), t(ldt,*), work(n), x(n)
      integer           istate(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  obj
      integer           i, itn, j, k, kadd, kdel, nart, ndf
      logical           first, linobj, newset, prthdr
      character*2       ladd, ldel
c     .. local arrays ..
      character*2       lstate(0:5)
      character*120     rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
c     .. data statements ..
      data              lstate(0), lstate(1), lstate(2)/'  ', 'l ',
     *                  'u '/
      data              lstate(3), lstate(4), lstate(5)/'e ', 'f ',
     *                  'a '/
c----------------------------------------------------------------------
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

c           terse line for the monitoring file.

            newset = lines1 .ge. mline1
            prthdr = msglvl .ge. 15 .or. first .or. newset
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
         if (iprint.ge.0 .and. isumm.ne.iprint) then
c           
c           terse line for the print file.
c           
            newset = lines2 .ge. mline2
            prthdr = first .or. newset
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

c                 print the diagonals of  t  and  r.

                  if (nactiv.gt.0) then
                     call dcopy (nactiv,t(nactiv,nz+1),ldt-1,work,1)
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
     *                    )
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
99993 format (/' values and status of the ',a2,' constraints',/' --',
     *       '-')
99992 format (/' variables...')
99991 format (/' general linear constraints...')
99990 format (/' diagonals of ',a2,' working set factor t')
99989 format (/' diagonals of ',a2,' triangle r         ')
99988 format (//' -',
     *       '-')
99987 format (1x,5(1p,d15.6,i5))
99986 format (1p,5d15.6)
      end

      subroutine e04mfh(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,a,
     *                  featol,cvec,x,wtinf)
c----------------------------------------------------------------------
c     e04mfh  finds the number and weighted sum of infeasibilities for
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
c     .. scalar arguments ..
      double precision  bigbnd, suminf
      integer           lda, n, nclin, numinf
c     .. array arguments ..
      double precision  a(lda,*), bl(*), bu(*), cvec(n), featol(*),
     *                  wtinf(*), x(n)
      integer           istate(*)
c     .. local scalars ..
      double precision  biglow, bigupp, ctx, feasj, s, weight
      integer           j, k
c     .. external functions ..
      double precision  ddot
      external          ddot
c----------------------------------------------------------------------
      bigupp = bigbnd
      biglow = -bigbnd
c
      numinf = 0
      suminf = zero
      call sload (n,(zero),cvec,1)
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
               call daxpy (n,weight,a(k,1),lda,cvec,1)
            end if
         end if
   40 continue

c     end of  e04mfh.  (cmsinf)

      end

      subroutine e04nfp(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,t,gqm,q,work,c,s)
c----------------------------------------------------------------------
c     e04nfp   updates the matrices  z, y, t, r  and  d  associated with
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
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      integer           it, jdel, kdel, lda, ldq, ldt, n, nactiv, nfree,
     *                  ngq, nrz, nz
      logical           unitq
c     .. array arguments ..
      double precision  a(lda,*), c(n), gqm(n,*), q(ldq,*), s(n),
     *                  t(ldt,*), work(n)
      integer           kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)
c     .. local scalars ..
      double precision  cs, sn
      integer           i, ir, itdel, j, jart, jt, k, l, npiv, nrz1,
     *                  nsup
c     .. external functions ..
      integer           idamax
      external          idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /de04nb/asize, dtmax, dtmin
      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------
      jt = nz + 1

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

            if ( .not. unitq) then

c              copy the incoming column of  a(free)  into the end of  t.

               do 20 k = 1, nactiv
                  i = kactiv(k)
                  t(nactiv-k+1,nfree) = a(i,jdel)
   20          continue

c              expand  q  by adding a unit row and column.

               if (nfree.gt.1) then
                  call sload (nfree-1,zero,q(nfree,1),ldq)
                  call sload (nfree-1,zero,q(1,nfree),1)
               end if
               q(nfree,nfree) = one
            end if
         else

c           case 2.  a general constraint has been deleted.
c           =======

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

            if (nsup.gt.0) then
               npiv = jt + itdel - 1
               if (nsup.gt.1) then
                  call dcopy (nsup-1,t(it+1,jt+1),ldt+1,s(jt+1),1)
                  call f06qrf('r',nactiv,1,nsup,c(jt+1),s(jt+1),
     *                        t(it,jt+1),ldt)
               end if
c
               call f06baf(t(it,jt+1),t(it,jt),cs,sn)
               t(it,jt) = zero
               s(jt) = -sn
               c(jt) = cs
               call f06qxf('r','v','backwards',nfree,nfree,
     *                     nz,npiv,c,s,q,ldq)
               call f06qxf('left ','v','backwards',npiv,ngq,nz,
     *                     npiv,c,s,gqm,n)
            end if
c
            jt = jt + 1
            call f06flf(nactiv,t(it,jt),ldt+1,dtmax,dtmin)
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

c     end of  e04nfp.  (rzdel)

99999 format (/' //e04nfp //  columns nrz and jart swapped.       ',
     *       /' //e04nfp //      nz   nrz   jart                 ',
     *       /' //e04nfp //  ',3i6)
99998 format (/' //e04nfp //  simple bound deleted.               ',/
     *       ' //e04nfp //  nactiv   nrz    nz nfree    ir  jdel unitq',
     *       /' //e04nfp //  ',6i6,l6)
99997 format (/' //e04nfp //  general constraint deleted.         ',
     *       /' //e04nfp //  nactiv   nrz    nz nfree  kdel  jdel',
     *       /' //e04nfp //  ',6i6)
      end

      subroutine e04mfs(firstv,n,nclin,istate,bigalf,bigbnd,pnorm,
     *                  hitlow,move,onbnd,unbndd,alfa,alfap,jhit,anorm,
     *                  ap,ax,bl,bu,featol,featlu,p,x)
c----------------------------------------------------------------------
c     e04mfs  finds a step alfa such that the point x + alfa*p reaches
c     one of the linear constraints (including bounds).
c
c     in this version of e04mfs, when x is infeasible, the number of
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

c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            if alfa reaches a bound as large as bigbnd, it is
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
c     bu     the upper bounds on ditto.
c     x      the values of       ditto.
c     p      the search direction.

c     output parameters

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
c     e04mfs is based on minos 5.2 routine m5chzr, which implements the
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
c     .. scalar arguments ..
      double precision  alfa, alfap, bigalf, bigbnd, pnorm
      integer           jhit, n, nclin
      logical           firstv, hitlow, move, onbnd, unbndd
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  alfai, atp, atpabs, atpmxf, atpmxi, atpscd, atx,
     *                  biglow, bigupp, bound, delta, exact, res,
     *                  stepmn, tolpiv
      integer           i, j, jhitf, jhiti, js
      logical           blockf, blocki
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04mf/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg
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
c
            if (atpscd.le.tolpiv) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

               res = -one
c
            else if (atp.le.zero .and. js.ne.-2) then

c              a'x  is decreasing and the lower bound is not violated.

c              test for smaller alfap.
c              if the upper bound is violated. test for bigger atp.
c
               if (bl(j).gt.biglow) then
                  res = atx - bl(j) + delta
c
                  if (res.lt.alfap*atpabs) alfap = res/atpabs
               end if
c
               if (js.eq.-1) atpmxi = max(atpmxi,atpscd)
c
            else if (atp.gt.zero .and. js.ne.-1) then

c              a'x  is increasing and the upper bound is not violated.

c              test for smaller alfap.
c              if the lower bound is violated. test for bigger atp.
c
               if (bu(j).lt.bigupp) then
                  res = bu(j) - atx + delta
c
                  if (res.lt.alfap*atp) alfap = res/atp
               end if
c
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
c
               if (atpscd.gt.atpmxf) then
c
                  if (bl(j).gt.biglow) then
                     res = atx - bl(j)
c
                     if (res.le.alfap*atpabs) then
                        atpmxf = atpscd
                        jhitf = j
                     end if
                  end if
               end if
c
               if (js.eq.-1) then
c
c                 the upper bound is violated.
c                 test for bigger or smaller alfai,  depending on the
c                 value of firstv.
c
                  if (firstv) then
                     res = atx - bu(j)
c
                     if (res.le.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
c
                  else if (atpscd.ge.atpmxi) then
                     res = atx - bu(j)
c
                     if (res.gt.alfai*atpabs) then
                        alfai = res/atpabs
                        jhiti = j
                     end if
                  end if
               end if
c
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
c
      if (unbndd) go to 60
c
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
c
c     if there is a choice between alfaf and alfai, it is probably best
c     to take alfai.  however, we can't if alfai is bigger than alfap.
c
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

c     end of  e04mfs.  (cmchzr)

99999 format (/' e04mfs entered.  pass 1.',/'    j  js         featol ',
     *       '       res             ap            alfap           atp',
     *       'mxi',/)
99998 format (i5,i4,3g15.5,2g17.7)
99997 format (/'                  pass 2.',/'    j  js         featol ',
     *       '       res             ap     jhitf           atpmxf jhi',
     *       'ti            alfai',/)
99996 format (i5,i4,3g15.5,2(i6,g17.7))
99995 format (/' //e04mfs//  unbounded step.',/' //e04mfs//  jhit     ',
     *       '      alfa',/' //e04mfs//  ',i4,g15.4)
      end

      subroutine e04nbx(msglvl,nfree,nrowa,n,nclin,nctotl,bigbnd,named,
     *                  names,nactiv,istate,kactiv,kx,a,bl,bu,c,clamda,
     *                  rlamda,x)
c----------------------------------------------------------------------
c     e04nbx   creates the expanded lagrange multiplier vector clamda.
c     if msglvl .eq 1 or msglvl .ge. 10,  e04nbx prints  x,  a*x,
c     c(x),  their bounds, the multipliers, and the residuals (distance
c     to the nearer bound).
c
c     e04nbx is called by e04ncz, e04ucz and e04upz just before exiting.

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd
      integer           msglvl, n, nactiv, nclin, nctotl, nfree, nrowa
      logical           named
c     .. array arguments ..
      double precision  a(nrowa,*), bl(nctotl), bu(nctotl), c(*),
     *                  clamda(nctotl), rlamda(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)
c     .. local scalars ..
      double precision  b1, b2, res, res2, v, wlam
      integer           ip, is, j, k, nfixed, nplin, nz
      character*2       ls
      character*5       id3
      character*8       id4
c     .. local arrays ..
      character*2       lstate(7)
      character*5       id(3)
      character*80      rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
c     .. data statements ..
      data              id(1)/'varbl'/
      data              id(2)/'l con'/
      data              id(3)/'n con'/
      data              lstate(1)/'--'/, lstate(2)/'++'/
      data              lstate(3)/'fr'/, lstate(4)/'ll'/
      data              lstate(5)/'ul'/, lstate(6)/'eq'/
      data              lstate(7)/'tf'/
c----------------------------------------------------------------------
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

      if (msglvl.lt.10 .and. msglvl.ne.1) return

      write (rec,fmt=99999)
      call x04bay(iprint,4,rec)
      id3 = id(1)
c
      do 40 j = 1, nctotl
         b1 = bl(j)
         b2 = bu(j)
         wlam = clamda(j)
         is = istate(j)
         ls = lstate(is+3)
         if (j.le.n) then

c           section 1 -- the variables  x.

            k = j
            v = x(j)

         else if (j.le.nplin) then

c           section 2 -- the linear constraints  a*x.

            if (j.eq.n+1) then
               write (rec,fmt=99998)
               call x04bay(iprint,4,rec)
               id3 = id(2)
            end if
c
            k = j - n
            v = ddot(n,a(k,1),nrowa,x,1)
         else
c
c           section 3 -- the nonlinear constraints  c(x).
c           
c
            if (j.eq.nplin+1) then
               write (rec,fmt=99997)
               call x04bay(iprint,4,rec)
               id3 = id(3)
            end if
c
            k = j - nplin
            v = c(k)
         end if
c
c        print a line for the j-th variable or constraint.
c        -
         res = v - b1
         res2 = b2 - v
         if (abs(res).gt.abs(res2)) res = res2
         ip = 1
         if (b1.le.(-bigbnd)) ip = 2
         if (b2.ge.bigbnd) ip = ip + 2
         if (named) then
c
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
c
         else
c
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

      subroutine e04nfr(unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
     *                  nfree,nrz,ngq,n,lda,ldq,ldr,ldt,kx,condmx,drzz,
     *                  a,r,t,gqm,q,w,c,s,msglvl)
c----------------------------------------------------------------------
c     e04nfr  updates the matrices  z, y, t, r  and  d  associated with
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
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  condmx, drzz
      integer           iadd, ifix, inform, it, jadd, lda, ldq, ldr,
     *                  ldt, msglvl, n, nactiv, nfree, ngq, nrz, nz
      logical           rset, unitq
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  cond, condbd, dtnew, tdtmax, tdtmin
      integer           i, j, jt, k, nanew, npiv, nsup
      logical           bound, overfl
c     .. local arrays ..
      character*80      rec(5)
c     .. external functions ..
      double precision  dnrm2, adivb
      external          dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04nb/asize, dtmax, dtmin
      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------
c     if the condition estimator of the updated t is greater than
c     condbd,  a warning message is printed.

      condbd = one/epspt9

      overfl = .false.
      bound = jadd .le. n
      jt = nz + 1

      if (bound) then

c        a simple bound has entered the working set.  iadd is not used.

         nanew = nactiv

         if (unitq) then
c
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
c
            call dcopy (nfree,q(ifix,1),ldq,w,1)
            if (ifix.lt.nfree) then
               call dcopy (nfree,q(nfree,1),ldq,q(ifix,1),ldq)
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else

c        a general constraint has entered the working set.
c        ifix is not used.

         nanew = nactiv + 1

c        transform the incoming row of a by q'.

         call dcopy (n,a(iadd,1),lda,w,1)
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
c
         if (unitq) then

c           first general constraint added.  set  q = i.

            call f06qhf('g',nfree,nfree,zero,one,q,ldq)
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

         if (ngq.gt.0) call f06qkf('l','t',nfree-1,w,ngq,gqm,
     *                             n)

         if (rset) then

c           apply the pairwise interchanges to  rz.
c           the subdiagonal elements generated by this process are
c           stored in  s(ifix), s(2), ..., s(nrz-1).

            nsup = nrz - ifix
            call e04nfm('r',nrz,ifix,nrz,s,r,ldr)
         end if
      else

c        the matrix  q  is stored explicitly.
c        define a sweep of plane rotations p such that
c                           pw = beta*e(npiv).
c        the rotations are applied in the planes (1, 2), (2, 3), ...,
c        (npiv-1, npiv).  the rotations must be applied to q, gqm', r
c        and t.

         call f06fqf('v','f',npiv-1,w(npiv),w,1,c,s)

         if (ngq.gt.0) call f06qxf('l','v','f',npiv,
     *                             ngq,1,npiv,c,s,gqm,n)
         call f06qxf('r','v','f',nfree,nfree,1,npiv,c,s,q,ldq)

         if (rset) then

c           apply the rotations to the triangular part of r.
c           the subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nrz-1).

            nsup = nrz - 1
            call f06qvf('r',nrz,1,nrz,c,s,r,ldr)
         end if
      end if

      if (rset) then

c        eliminate the  nsup  subdiagonal elements of  r  stored in
c        s(nrz-nsup), ..., s(nrz-1)  with a left-hand sweep of rotations
c        in planes (nrz-nsup, nrz-nsup+1), ..., (nrz-1, nrz).

         call f06qrf('left ',nrz,nrz-nsup,nrz,c,s,r,ldr)

         if (nsup.gt.0 .and. drzz.ne.one) then
            drzz = c(nrz-1)**2 + drzz*s(nrz-1)**2
         end if
      end if

      if ( .not. unitq) then
         if (bound) then

c           bound constraint added.   the rotations affect columns
c           nz+1  thru  nfree  of  gqm'  and  t.

c           the last row and column of  q  has been transformed to plus
c           or minus the unit vector  e(nfree).  reconstitute the
c           column of gqm' corresponding to the new fixed variable.
c
            if (w(nfree).lt.zero) then
               if (ngq.gt.0) call dscal (ngq,-one,gqm(nfree,1),n)
            end if
c
            if (nactiv.gt.0) then
               t(it,jt-1) = s(jt-1)*t(it,jt)
               t(it,jt) = c(jt-1)*t(it,jt)
c
               if (nactiv.gt.1) then
                  call f06qvf('r',nactiv,1,nactiv,c(jt),s(jt),
     *                        t(it,jt),ldt)
                  call dcopy (nactiv-1,s(jt),1,t(it+1,jt),ldt+1)
               end if
c
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
            if (cond.ge.condbd) then
               if (msglvl.gt.0) then
                  write (rec,fmt=99997) jadd
                  call x04bay(iprint,5,rec)
               end if
            end if
         else
c           the proposed working set appears to be linearly dependent.

            inform = 1

         end if
      end if

c     end of  e04nfr.  (rzadd)

99999 format (/' //e04nfr //  simple bound added.',/' //e04nfr //  nac',
     *       'tiv   nrz    nz nfree  ifix  jadd unitq         ',/' //e',
     *       '04nfr //  ',6i6,l6)
99998 format (/' //e04nfr //  general constraint added.           ',/
     *       ' //e04nfr //  nactiv   nrz    nz nfree  iadd  jadd unitq',
     *       /' //e04nfr //  ',6i6,l6)
99997 format (/' xxx  serious ill-conditioning in the working set afte',
     *       'r adding constraint ',i5,/' xxx  overflow may occur in s',
     *       'ubsequent iterations.',//)
99996 format (/' //e04nfr //  dependent constraint rejected.')
99995 format (/' //e04nfr //     asize     dtmax     dtmin        ',
     *       /' //e04nfr //',1p,3d10.2)
99994 format (/' //e04nfr //     asize     dtmax     dtmin     dtnew',
     *       /' //e04nfr //',1p,4d10.2)
99993 format (/' //e04nfr //     asize     dtnew',/' //e04nfr //',1p,
     *       2d10.2)
      end

      subroutine e04ucn(feasqp,n,nclin,ncnln,objalf,grdalf,qpcurv,
     *                  istate,cjdx,cmul,cs,dlam,rho,violn,work1,work2)
c----------------------------------------------------------------------
c     e04ucn   computes the value and directional derivative of the
c     augmented lagrangian merit function.  the penalty parameters
c     rho(j) are boosted if the directional derivative of the resulting
c     augmented lagrangian function is not sufficiently negative.  if
c     rho needs to be increased,  the perturbation with minimum two-norm
c     is found that gives a directional derivative equal to  - p'hp.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two
      parameter         (two=2.0d+0)
c     .. scalar arguments ..
      double precision  grdalf, objalf, qpcurv
      integer           n, nclin, ncnln
      logical           feasqp
c     .. array arguments ..
      double precision  cjdx(*), cmul(*), cs(*), dlam(*), rho(*),
     *                  violn(*), work1(*), work2(*)
      integer           istate(*)
c     .. scalars in common ..
      double precision  rhodmp, rhomax, rhonrm, scale
      integer           iprint, isumm, lines1, lines2, nout
      logical           incrun, npdbg
c     .. arrays in common ..

      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  pterm, pterm2, qnorm, rho1, rhoi, rhomin,
     *                  rhonew, rtmin, tscl
      integer           i, nplin
      logical           boost, overfl
c     .. external functions ..
      double precision  ddot, dnrm2, adivb
      external          ddot, dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
c
      if (ncnln.eq.0) return
c
      rtmin = wmach(6)
c
      objalf = objalf - ddot(ncnln,cmul,1,cs,1)
      grdalf = grdalf - ddot(ncnln,dlam,1,cs,1)
c
      call dcopy (ncnln,cs,1,work1,1)
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
         tscl = adivb(grdalf+half*qpcurv,qnorm,overfl)
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
         call f06fcf(ncnln,rho,1,work1,1)
         pterm2 = ddot(ncnln,work1,1,cs,1)
c
         scale = rhomax
         tscl = adivb(grdalf,pterm2,overfl)
         if (tscl.gt.scale .and. tscl.le.rhomax/(one+rhonrm)
     *       .and. .not. overfl) scale = tscl
c
         call dcopy (ncnln,cs,1,work1,1)
      end if
c

c     compute the new value and directional derivative of the
c     merit function.

      call f06fcf(ncnln,rho,1,work1,1)
c
      pterm = ddot(ncnln,work1,1,cs,1)
      objalf = objalf + half*scale*pterm
c
      if (feasqp) pterm2 = pterm
c
      grdalf = grdalf - scale*pterm2

c     end of  e04ucn. (npmrt)

99999 format (/' //e04ucn //        qpcurv        grdalf ',/' //e04ucn',
     *       ' //',1p,2d14.2)
99998 format (/' //e04ucn //         scale        rhonrm        grdalf '
     *       ,/' //e04ucn //',1p,3d14.2)
99997 format (/' //e04ucn //  penalty parameters =       ')
99996 format (1p,5d15.6)
      end

      subroutine e04xay(inform,msglvl,lvlder,n,ncnln,ldcj,ldcju,bigbnd,
     *                  epsrf,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,c2,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,y,
     *                  iuser,user)
c----------------------------------------------------------------------
c     e04xay  computes difference intervals for the missing gradients of
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
      double precision  factor
      parameter         (factor=0.97d+0)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two, four, ten
      parameter         (two=2.0d+0,four=4.0d+0,ten=1.0d+1)
c     .. scalar arguments ..
      double precision  bigbnd, epsrf, fdnorm, objf
      integer           inform, ldcj, ldcju, lvlder, msglvl, n, ncnln
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  biglow, bigupp, cdest, cjdiff, d, dx, epsa,
     *                  errbnd, errmax, errmin, f1, f2, fdest, fx,
     *                  gdiff, h, hcd, hfd, hmax, hmin, hopt, hphi,
     *                  objf2, sdest, signh, stepbl, stepbu, sumeps,
     *                  sumsd, test, xj, yj
      integer           i, info, irow1, irow2, iter, itmax, j, mode,
     *                  nccnst, ncolj, nfcnst, nstate
      logical           debug, done, first, headng, needed
c     .. local arrays ..
      character*80      rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
      inform = 0
      needed = lvlder .eq. 0 .or. lvlder .eq. 2 .or. lvlder .eq. 1 .and.
     *         ncnln .gt. 0
      if ( .not. needed) return

      if (lfdset.eq.0) then
         if (msglvl.gt.0) then
            write (rec,fmt=99999)
            call x04bay(iprint,4,rec)
         end if

         nstate = 0
         itmax = 3
         mode = 0

         nccnst = 0
         nfcnst = 0
         headng = .true.

         fdnorm = zero

c        for each column of the jacobian augmented by the transpose of
c        the objective gradient, rows irow1 thru irow2 are searched for
c        missing elements.

         irow1 = 1
         irow2 = ncnln + 1
         if (lvlder.eq.1) irow2 = ncnln
         if (lvlder.eq.2) irow1 = ncnln + 1

         biglow = -bigbnd
         bigupp = bigbnd

         if (ncnln.gt.0) call f06dbf(ncnln,(0),needc,1)

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

            stepbl = biglow
            stepbu = bigupp
            if (bl(j).gt.biglow) stepbl = bl(j) - xj
            if (bu(j).lt.bigupp) stepbu = bu(j) - xj

            signh = one
            if (half*(stepbl+stepbu).lt.zero) signh = -one

            do 40 i = irow1, irow2

               if (i.le.ncnln) then
                  test = cjacu(i,j)
               else
                  test = gradu(j)
               end if

               if (test.eq.rdummy) then

c                 get the difference interval for this element.

                  ncolj = ncolj + 1

                  if (i.le.ncnln) then
                     needc(i) = 1
                     fx = c(i)
                     epsa = epsrf*(one+abs(c(i)))
                  else
                     fx = objf
                     epsa = epsrf*(one+abs(fx))
                  end if

c                 find a finite-difference interval by iteration.

                  iter = 0
                  hopt = two*(one+abs(xj))*sqrt(epsrf)
                  h = signh*ten*hopt
                  cdest = zero
                  sdest = zero
                  first = .true.

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

                  call chcore(debug,done,first,epsa,epsrf,fx,info,iter,
     *                        itmax,cdest,fdest,sdest,errbnd,f1,f2,h,
     *                        hopt,hphi)

                  if ( .not. done) go to 20

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

                  if (info.eq.0) hcd = max(hcd,hphi)
               end if
   40       continue

            if (ncolj.gt.0) then
               if (hmin.gt.hmax) then
                  hmin = hmax
                  errmin = errmax
               end if

               if (four*sumeps.lt.hmin*hmin*sumsd) then
                  hfd = hmin
                  errmax = errmin
               else if (four*sumeps.gt.hmax*hmax*sumsd) then
                  hfd = hmax
               else
                  hfd = two*sqrt(sumeps/sumsd)
                  errmax = two*sqrt(sumeps*sumsd)
               end if

               if (hcd.eq.zero) hcd = ten*hfd

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

         if (nccnst+nfcnst.gt.0) then

c           check that the constants have been set properly by
c           evaluating the gradients at a strange (but feasible) point.

            d = one/n

            do 80 j = 1, n
               xj = x(j)
               stepbl = -one
               stepbu = one
               if (bl(j).gt.biglow) stepbl = max(stepbl,bl(j)-xj)
               if (bu(j).lt.bigupp .and. bu(j).gt.bl(j))
     *             stepbu = min(stepbu,bu(j)-xj)

               if (half*(stepbl+stepbu).lt.zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if

               d = factor*d
   80       continue

            if (ncnln.gt.0) then
               call f06dbf(ncnln,(1),needc,1)
               call confun(mode,ncnln,n,ldcju,needc,y,c2,cjacu,nstate,
     *                     iuser,user)
               if (mode.lt.0) go to 200
            end if

            call objfun(mode,n,y,objf2,gradu,nstate,iuser,user)
            if (mode.lt.0) go to 200

c           loop over each of the elements of  x.

            do 140 j = 1, n
               yj = y(j)
               dx = half*(x(j)-yj)
               y(j) = yj + dx

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

c              check the objective gradient element.

               if (gradu(j).eq.-rdummy) then

                  call objfun(mode,n,y,f1,gradu,nstate,iuser,user)
                  if (mode.lt.0) go to 200

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
            if (msglvl.gt.0) then
               if (lvlder.lt.2 .and. nccnst.gt.0) then
                  write (rec,fmt=99996) nccnst
                  call x04bay(iprint,2,rec)
               end if
               if (lvlder.ne.1 .and. nfcnst.gt.0) then
                  write (rec,fmt=99995) nfcnst
                  call x04bay(iprint,2,rec)
               end if
            end if
c
            if (ncdiff.eq.0 .and. lvlder.lt.2) then
               if (lvlder.eq.0) lvlder = 2
               if (lvlder.eq.1) lvlder = 3
               if (msglvl.gt.0) then
                  write (rec,fmt=99994) lvlder
                  call x04bay(iprint,4,rec)
               end if
            end if
c
            if (nfdiff.eq.0 .and. lvlder.ne.1) then
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
c
  200 inform = mode

c     end of  e04xay. (chfd)

99999 format (//' computation of the finite-difference intervals',/' -',
     *       '')
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

      subroutine e04udw(string,first,last,mark)
c----------------------------------------------------------------------
c     looks for non-blank fields ('tokens') in a string, where the
c     fields are of arbitrary length, separated by blanks, tabs, commas,
c     colons, or equal signs.  the position of the end of the 1st token
c     is also returned, so this routine may be conveniently used within
c     a loop to process an entire line of text.
c
c     the procedure examines a substring, string (first : last), which
c     may of course be the entire string (in which case just call e04udw
c     with first .le. 1 and last .ge. len (string) ).  the indices
c     returned are relative to string itself, not the substring.

c     parameters:

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


c     (2)  the pseudo-recursive structure was chosen for fun.  it is
c         equivalent to three do loops with embedded go tos in sequence.
c
c     (3)  the variety of separators recognized limits the usefulness of
c         this routine somewhat.  the intent is to facilitate handling
c         such tokens as keywords or numerical values.  in other
c         applications, it may be necessary for all printing characters
c         to be significant.  a simple modification to statement
c         function solid will do the trick.

      character         blank, equal, colon, comma, rparn, lparn
      parameter         (blank=' ',equal='=',colon=':',comma=',',
     *                  rparn=')',lparn='(')
c     .. scalar arguments ..
      integer           first, last, mark
      character*(*)     string
c     .. local scalars ..
      integer           begin, end, length
      character         dummy
c     .. statement functions ..
      logical           solid
c     .. statement function definitions ..
      solid(dummy) = (dummy.ne.blank) .and. (dummy.ne.colon)
     *               .and. (dummy.ne.comma) .and. (dummy.ne.equal)
     *               .and. (dummy.ne.rparn) .and. (dummy.ne.lparn)
c----------------------------------------------------------------------
      mark = 0
      length = len(string)
      begin = max(first,1)
      end = min(length,last)

c     find the first significant character ...

      do 60 first = begin, end, +1
         if (solid(string(first:first))) then

c           ... then the end of the first token ...

            do 40 mark = first, end - 1, +1
               if ( .not. solid(string(mark+1:mark+1))) then
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

      subroutine e04ucw(n,nclin,ncnln,istate,bigbnd,cvnorm,errmax,jmax,
     *                  nviol,ax,bl,bu,c,featol,x,work)
c----------------------------------------------------------------------
c     e04ucw  computes 
c          the number of constraints that are violated by more
c          than  featol  and the 2-norm of the constraint violations.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, cvnorm, errmax
      integer           jmax, n, nclin, ncnln, nviol
c     .. array arguments ..
      double precision  ax(*), bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     *                  c(*), featol(n+nclin+ncnln),
     *                  work(n+nclin+ncnln), x(n)
      integer           istate(n+nclin+ncnln)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  biglow, bigupp, con, feasj, res, tolj
      integer           is, j
c     .. external functions ..
      double precision  dnrm2
      integer           idamax
      external          dnrm2, idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
      biglow = -bigbnd
      bigupp = bigbnd

c     compute nviol, the number of constraints violated by more than
c     featol,  and cvnorm,  the 2-norm of the constraint
c     violations and residuals of the constraints in the qp working set.

      nviol = 0
c
      do 40 j = 1, n + nclin + ncnln
         feasj = featol(j)
         res = zero

         if (j.le.n+nclin) then

c           bound or general linear constraint.

            if (j.le.n) then
               con = x(j)
            else
               con = ax(j-n)
            end if

            tolj = feasj
         else

c           nonlinear constraint.

            con = c(j-n-nclin)
            tolj = zero
         end if

c        check for constraint violations.

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

c     end of  e04ucw. (npfeas)

99999 format (/' //e04ucw//  the maximum violation is ',1p,d14.2,' in ',
     *       'constraint',i5)
      end

      subroutine e04ucr(needfd,inform,n,ncnln,ldcj,ldcju,nfun,ngrad,
     *                  needc,confun,objfun,alfa,alfbnd,alfmax,alfsml,
     *                  dxnorm,epsrf,eta,gdx,grdalf,glf1,glf,objf,
     *                  objalf,qpcurv,xnorm,c,c2,cjac,cjacu,cjdx,cjdx2,
     *                  cmul1,cmul,cs1,cs,dx,dlam,dslk,grad,gradu,qpmul,
     *                  rho,slk1,slk,x1,x,work,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     e04ucr finds the steplength alfa that gives sufficient decrease in
c     the augmented lagrangian merit function.
c
c     on exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
c     with an associated merit function value  objalf  which is lower
c     than that at the base point. if  inform = 4, 5, 6, 7 or 8,  alfa
c     is zero and  objalf will be the merit value at the base point.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
      double precision  two
      parameter         (two=2.0d+0)
      double precision  tolg
      parameter         (tolg=1.0d-1)
      double precision  rmu
      parameter         (rmu=1.0d-4)
c     .. scalar arguments ..
      double precision  alfa, alfbnd, alfmax, alfsml, dxnorm, epsrf,
     *                  eta, gdx, glf, glf1, grdalf, objalf, objf,
     *                  qpcurv, xnorm
      integer           inform, ldcj, ldcju, lenw, n, ncnln, nfun, ngrad
      logical           needfd
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  alfbst, cs1jdx, csjdx, curvc, curvlf, epsaf,
     *                  epsmch, fbest, fterm, ftry, g0, gbest, gtry,
     *                  oldf, oldg, q, rhobfs, s, t, targtg, tgdx, tglf,
     *                  tobj, tobjm, tolabs, tolax, tolrel, tolrx,
     *                  toltny
      integer           j, maxf, mode, nstate, numf
      logical           debug, done, first, imprvd
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /npdebg/inpdbg, npdbg
c----------------------------------------------------------------------
      epsmch = wmach(3)

      if ( .not. needfd .and. ncnln.gt.0) cs1jdx = ddot(ncnln,cs1,1,
     *    cjdx,1)

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

      if (needfd) then
         maxf = 15
      else
         maxf = 10
      end if

      epsaf = epsrf*(one+abs(objalf))
      tolax = epspt8
      tolrx = epspt8

      if (tolrx*xnorm+tolax.lt.dxnorm*alfmax) then
         tolabs = (tolrx*xnorm+tolax)/dxnorm
      else
         tolabs = alfmax
      end if
      tolrel = max(tolrx,epsmch)

      t = zero
      do 20 j = 1, n
         s = abs(dx(j))
         q = abs(x(j))*tolrx + tolax
         if (s.gt.t*q) t = s/q
   20 continue

      if (t*tolabs.gt.one) then
         toltny = one/t
      else
         toltny = tolabs
      end if

      oldf = objalf
      oldg = grdalf
      alfbst = zero
      fbest = zero
      gbest = (one-rmu)*oldg
      targtg = (rmu-eta)*oldg
      g0 = gbest

      if (ncnln.gt.0) call f06dbf(ncnln,(1),needc,1)

      if (needfd) then
         mode = 0
      else
         mode = 2
      end if

      first = .true.

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

   40 if (needfd) then
         call e04ucj(first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,alfsml,epsaf,g0,targtg,ftry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest)
      else
         call e04uck(first,debug,done,imprvd,inform,maxf,numf,iprint,
     *               alfmax,epsaf,g0,targtg,ftry,gtry,tolabs,tolrel,
     *               toltny,alfa,alfbst,fbest,gbest)
      end if

      if (imprvd) then
         objf = tobj
         objalf = tobjm

         if (ncnln.gt.0) call dcopy (ncnln,c2,1,c,1)

         if ( .not. needfd) then
            call dcopy (n,gradu,1,grad,1)
            gdx = tgdx
            glf = tglf

            if (ncnln.gt.0) then
               call dcopy (ncnln,cjdx2,1,cjdx,1)
               call f06qff('g',ncnln,n,cjacu,ldcju,cjac,ldcj)
            end if
         end if
      end if

c     if done = .false.,  the problem functions must be computed for
c     the next entry to e04uck or e04ucj.
c     if done = .true.,   this is the last time through.

      if ( .not. done) then

         call dcopy (n,x1,1,x,1)
         call daxpy (n,alfa,dx,1,x,1)

         if (ncnln.gt.0) then

c           compute the new estimates of the multipliers and slacks.
c           if the step length is greater than one,  the multipliers
c           are fixed as the qp-multipliers.

            if (alfa.le.one) then
               call dcopy (ncnln,cmul1,1,cmul,1)
               call daxpy (ncnln,alfa,dlam,1,cmul,1)
            end if
            call dcopy (ncnln,slk1,1,slk,1)
            call daxpy (ncnln,alfa,dslk,1,slk,1)

c           compute the new constraint vector and jacobian.

            call confun(mode,ncnln,n,ldcju,needc,x,c2,cjacu,nstate,
     *                  iuser,user)
            if (mode.lt.0) go to 60

            call dcopy (ncnln,c2,1,cs,1)
            call daxpy (ncnln,(-one),slk,1,cs,1)

            call dcopy (ncnln,cs,1,work,1)
            call f06fcf(ncnln,rho,1,work,1)

            fterm = ddot(ncnln,cmul,1,cs,1) - half*scale*ddot(ncnln,
     *              work,1,cs,1)
         end if

c        compute the value and gradient of the objective function.

         call objfun(mode,n,x,tobj,gradu,nstate,iuser,user)
         if (mode.lt.0) go to 60

         if (ncnln.gt.0) then
            tobjm = tobj - fterm
         else
            tobjm = tobj
         end if

         ftry = tobjm - oldf - rmu*oldg*alfa

c        compute auxiliary gradient information.

         if ( .not. needfd) then
            gtry = ddot(n,gradu,1,dx,1)
            tgdx = gtry
            tglf = gtry
            if (ncnln.gt.0) then

c              compute the jacobian times the search direction.

               call dgemv ('n',ncnln,n,one,cjacu,ldcju,dx,1,zero,
     *                     cjdx2,1)

               call dcopy (ncnln,cjdx2,1,work,1)
               call daxpy (ncnln,(-one),dslk,1,work,1)

               gtry = gtry - ddot(ncnln,cmul,1,work,1)
               if (alfa.le.one) gtry = gtry - ddot(ncnln,dlam,1,cs,1)

               call f06fcf(ncnln,rho,1,work,1)
               gtry = gtry + scale*ddot(ncnln,work,1,cs,1)

               tglf = tgdx - ddot(ncnln,cjdx2,1,qpmul,1)

c              if alfbnd .le. alfa .lt. alfmax and the norm of the
c              quasi-newton update is bounded, set alfmax to be alfa.
c              this will cause the linesearch to stop if the merit
c              function is decreasing at the boundary.

               if (alfbnd.le.alfa .and. alfa.lt.alfmax) then

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

            gtry = gtry - rmu*oldg

         end if
      end if

      if ( .not. done) go to 40

      nfun = nfun + numf
      if ( .not. needfd) ngrad = ngrad + numf
      alfa = alfbst

      if ( .not. imprvd) then
         call dcopy (n,x1,1,x,1)
         call daxpy (n,alfa,dx,1,x,1)
         if (ncnln.gt.0) then
            if (alfa.le.one) then
               call dcopy (ncnln,cmul1,1,cmul,1)
               call daxpy (ncnln,alfa,dlam,1,cmul,1)
            end if
            call dcopy (ncnln,slk1,1,slk,1)
            call daxpy (ncnln,alfa,dslk,1,slk,1)
            call dcopy (ncnln,c,1,cs,1)
            call daxpy (ncnln,(-one),slk,1,cs,1)
         end if
      end if

      return

   60 inform = mode

c     end of  e04ucr. (npsrch)

99999 format (/' //e04ucr// inform  = ',i4)
99998 format (/' //e04ucr//        alfbnd          alfa        alfmax',
     *       /' //e04ucr//',1p,3d14.2)
99997 format (/' //e04ucr//         csjdx        cs1jdx        curvlf',
     *       /' //e04ucr//',1p,3d14.2)
      end

      subroutine e04ucl(lsumry,unitq,n,ncnln,nfree,nz,ldcj1,ldcj2,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,cjac1,cjac2,cjdx1,
     *                  cjdx2,cs1,cs2,gq1,gq2,hpq,rpq,qpmul,r,omega,zy,
     *                  wrk1,wrk2)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      double precision  tolg
      parameter         (tolg=1.0d-1)
c     .. scalar arguments ..
      double precision  alfa, glf1, glf2, qpcurv
      integer           ldcj1, ldcj2, ldr, ldzy, n, ncnln, nfree, nz
      logical           unitq
      character*5       lsumry
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  beta, curvl, eta, qi, qmax, qnorm, rtgtp, rtyts,
     *                  test, tinycl, trace1, trace2
      integer           i, j
      logical           overfl, ssbfgs
c     .. local arrays ..
      character*80      rec(3)
c     .. external functions ..
      double precision  dnrm2, adivb
      integer           idamax
      external          dnrm2, adivb, idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /de04uc/rhomax, rhonrm, rhodmp, scale, incrun
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
      if (ncnln.gt.0) call sload (ncnln,zero,omega,1)

c     set curvl = (g2 - g1)'dx,  the approximate curvature along dx of
c     the (augmented) lagrangian.  at first, the curvature is not scaled
c     by the steplength alfa.

      curvl = glf2 - glf1
      tinycl = qpcurv*tolg
      ssbfgs = curvl .le. alfa*tinycl

c     test if curvl is sufficiently positive.  if there are no nonlinear
c     constraints,  no update can be performed.

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
            beta = adivb(qmax*test,qnorm*qnorm,overfl)
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

c     compute the difference in the augmented lagrangian gradient.

c     update gq1 to include the augmented lagrangian terms.

      if (ncnln.gt.0) then

         do 80 i = 1, ncnln
            wrk1(i) = -qpmul(i) + omega(i)*cs1(i)
   80    continue
         call dgemv ('t',ncnln,n,one,cjac1,ldcj1,wrk1,1,zero,wrk2,1)

         do 100 i = 1, ncnln
            wrk1(i) = qpmul(i) - omega(i)*cs2(i)
  100    continue
         call dgemv ('t',ncnln,n,one,cjac2,ldcj2,wrk1,1,one,wrk2,1)

         call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,wrk2,zy,wrk1)
         call daxpy (n,one,wrk2,1,gq1,1)
      end if

      if (npdbg .and. inpdbg(1).gt.0) then
         write (rec,fmt=99998) alfa, curvl
         call x04bay(iprint,3,rec)
      end if

      if (curvl.lt.tinycl) curvl = tinycl

      do 120 j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
  120 continue

      rtgtp = sqrt(qpcurv)
      rtyts = sqrt(alfa*curvl)
      eta = one
      if (ssbfgs) eta = rtyts/(rtgtp*alfa)

      trace1 = dnrm2(n,hpq,1)/rtgtp
      trace2 = dnrm2(n,wrk2,1)/(rtyts*eta)
      rfrobn = eta*sqrt(abs((rfrobn-trace1)*(rfrobn+trace1)+trace2**2))

c     update the cholesky factor of  q'hq.

c     normalize the vector  rpq ( = r(pq) ).
c
      call dscal (n,(one/rtgtp),rpq,1)
c
c     do the self-scaled or regular bfgs update.
c     form the vector wrk1 = gamma * (gq2 - gq1) - beta * r'r*pq,
c     where  gamma = 1/sqrt( curv ) = 1/sqrt( (gq2 - gq1)'sq )

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

c     update to  r = r + rpq*wrk1'.
c     rpq is overwritten. arrays gq1 and hpq are used to store the
c     sines and cosines defined by the plane rotations.

      call e04nbv(n,0,n,ldr,n,n,r,hpq,rpq,wrk1,gq1,hpq)

c     end of  e04ucl. (npupdt)

99999 format (/' //e04ucl// ssbfgs    min. curvl         curvl ',/' //',
     *       'e04ucl//   ',l4,1p,2d14.2)
99998 format (/' //e04ucl//          alfa         curvl ',
     *       /' //e04ucl//',1p,2d14.2)
99997 format (/' //e04ucl//   omega(imax)',/' //e04ucl//',1p,d14.2)
99996 format (/' //e04ucl//  penalty parameters = ')
99995 format (1p,5d15.6)
      end

      subroutine e04ucg(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfa,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  anorm,ap,ax,bl,bu,featol,p,x)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      implicit none

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  alfa, atphit, bigalf, bigbnd, palfa, pnorm
      integer           inform, jadd, n, nctotl, numinf
      logical           firstv, hitlow
c     .. array arguments ..
      double precision  anorm(*), ap(*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), p(n), x(n)
      integer           istate(nctotl)
c     .. scalars in common ..
      double precision  epspt3, epspt5, epspt8, epspt9
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)
c     .. local scalars ..
      double precision  absatp, alfa1, alfa2, apmax1, apmax2, atp, atp1,
     *                  atp2, atx, palfa1, palfa2, res, rownrm
      integer           i, j, jadd1, jadd2, js, jsave1, jsave2
      logical           hlow1, hlow2, lastv, negstp, step2
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------
      inform = 0

c     first pass -- find steps to perturbed constraints, so that
c     palfa1 will be slightly larger than the true step, and
c     palfa2 will be slightly smaller than it should be.
c     in degenerate cases, this strategy gives us some freedom in the
c     second pass.  the general idea follows that described by p.m.j.
c     harris, p.21 of mathematical programming 5, 1 (1973), 1--28.

      negstp = .false.
      call e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,palfa1,
     *            palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,featol,p,x)

      jsave1 = jadd1
      jsave2 = jadd2

c     second pass -- recompute step-lengths without perturbation.
c     amongst constraints that are less than the perturbed steps,
c     choose the one (of each type) that makes the largest angle
c     with the search direction.

      alfa1 = bigalf
      alfa2 = zero
      if (firstv) alfa2 = bigalf

      apmax1 = zero
      apmax2 = zero
      atp1 = zero
      atp2 = zero
      hlow1 = .false.
      hlow2 = .false.
      lastv = .not. firstv

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

            if (abs(atp).le.epspt9*rownrm*pnorm) then

c              this constraint appears to be constant along p.  it is
c              not used to compute the step.  give the residual a value
c              that can be spotted in the debug output.

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
         end if
   20 continue

c     determine alfa, the step to be taken.

c     in the infeasible case, check whether to take the step alfa2
c     rather than alfa1...
c
      step2 = numinf .gt. 0 .and. jadd2 .gt. 0
c
c     we do so if alfa2 is less than alfa1 or (if firstv is false)
c     lies in the range  (alfa1, palfa1)  and has a smaller value of
c     atp.

      step2 = step2 .and. (alfa2.lt.alfa1 .or. lastv .and. alfa2.le.
     *        palfa1 .and. apmax2.ge.apmax1)

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

c        if alfa1 is negative, the constraint to be added (jadd)
c        remains unchanged, but alfa may be shortened to the step
c        to the nearest perturbed satisfied constraint along  - p.

         negstp = alfa .lt. zero
         if (negstp) then
            call e04uch(firstv,negstp,bigalf,bigbnd,pnorm,jadd1,jadd2,
     *                  palfa1,palfa2,istate,n,nctotl,anorm,ap,ax,bl,bu,
     *                  featol,p,x)

            alfa = -min(abs(alfa),palfa1)
         end if
      end if

c     test for undefined or infinite step.

      if (jadd.eq.0) then
         alfa = bigalf
         palfa = bigalf
      end if

      if (alfa.ge.bigalf) inform = 3

c     end of  e04ucg. (cmalf)

99999 format (/' e04ucg  entered',/'    j  js         featol        re',
     *       's             ap     jadd1        alfa1     jadd2       ',
     *       ' alfa2 ',/)
99998 format (i5,i4,3g15.5,2(i6,g17.7))
99997 format (/' //e04ucg //  negative step',/' //e04ucg //           ',
     *       'alfa          palfa',/' //e04ucg //',2g15.4)
99996 format (/' //e04ucg //  unbounded step.',/' //e04ucg //  jadd   ',
     *       '        alfa',/' //e04ucg //  ',i4,g15.4)
      end

      subroutine e04ucq(nout,buffer,key)
c----------------------------------------------------------------------
c     e04ucq   decodes the option contained in  buffer  in order to set
c     a parameter value in the relevant element of the parameter arrays.

c     input:
c
c     nout   a unit number for printing error messages.
c            nout  must be a valid unit.
c
c     output:
c
c     key    the first keyword contained in buffer.
c
c     e04ucq  calls e04udx and the subprograms
c                 lookup, scannr, tokens
c     (now called e04udy, e04udw, e04udv)
c----------------------------------------------------------------------
      implicit none

      integer           mxparm
      parameter         (mxparm=30)
      integer           maxkey, maxtie, maxtok
      parameter         (maxkey=41,maxtie=19,maxtok=10)
      integer           idummy
      double precision  rdummy
      logical           sorted
      double precision  zero
      parameter         (idummy=-11111,rdummy=-11111.0d+0,sorted=.true.,
     *                  zero=0.0d+0)
c     .. scalar arguments ..
      integer           nout
      character*16      key
      character*(*)     buffer
c     .. scalars in common ..
      double precision  bigbnd, bigdx, bndlow, bndupp, cdint, ctol,
     *                  dxlim, epsrf, eta, fdint, ftol, hcndbd, tolact,
     *                  tolfea, tolrnk
      integer           idbgls, idbgnp, iprnt, isumry, itmax1, itmax2,
     *                  itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, ksave,
     *                  lcrash, ldbgls, ldbgnp, lformh, lprob, lverfy,
     *                  lvlder, msgls, msgnp, nlnf, nlnj, nlnx, nload,
     *                  nn, nnclin, nncnln, nprob, nsave
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm)
c     .. local scalars ..
      double precision  rvalue
      integer           i, idbg, ivalue, lenbuf, loc1, loc2, mjrdbg,
     *                  mnrdbg, msgqp, nmajor, nminor, ntoken
      logical           first, more, number
      character*16      key2, key3, value
      character*80      rec
c     .. local arrays ..
      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*16      keys(maxkey), ties(maxtie), token(maxtok)
c     .. external functions ..
      logical           e04udx
      external          e04udx
c     .. common blocks ..
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
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c     .. save statement ..
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/, first
c     .. data statements ..
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
c----------------------------------------------------------------------
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

      if (key.eq.' ' .or. key.eq.'begin') go to 80
      if (key.eq.'list' .or. key.eq.'nolist') go to 80
      if (key.eq.'end') go to 80

c     most keywords will have an associated integer or real value,
c     so look for it no matter what the keyword.

      i = 1
      number = .false.

   40 if (i.lt.ntoken .and. .not. number) then
         i = i + 1
         value = token(i)
         number = e04udx(value)
         go to 40
      end if

      if (number) then
         read (value,fmt='(bn, e16.0)') rvalue
      else
         rvalue = zero
      end if

c     convert the keywords to their most fundamental form
c     (upper case, no abbreviations).
c     sorted says whether the dictionaries are in alphabetic order.
c     loci   says where the keywords are in the dictionaries.
c     loci = 0 signals that the keyword wasn't there.
c     loci < 0 signals that the keyword is ambiguous.

      call e04udy(maxkey,keys,sorted,key,loc1)
      if (loc1.lt.0) then
         write (rec,fmt=99996) key
         call x04baf(nout,rec)
         return
      end if
      call e04udy(maxtie,ties,sorted,key2,loc2)

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
            if (ivalue.ge.1 .and. ivalue.le.mxparm) then
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
            if (key2.eq.'iterations  ') nmajor = rvalue
            if (key2.eq.'print       ') msgnp = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'minor       ') then
            if (key2.eq.'iterations  ') nminor = rvalue
            if (key2.eq.'print       ') msgqp = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
         else if (key.eq.'monitoring  ') then
            isumry = rvalue
         else if (key.eq.'nonlinear   ') then
            if (key2.eq.'constraints ') nncnln = rvalue
            if (key2.eq.'feasibility ') ctol = rvalue
            if (key2.eq.'jacobian    ') nlnj = rvalue
            if (key2.eq.'objective   ') nlnf = rvalue
            if (key2.eq.'variables   ') nlnx = rvalue
            if (loc2.eq.0) then
               write (rec,fmt=99998) key2
               call x04baf(nout,rec)
            end if
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
            if (ivalue.ge.1 .and. ivalue.le.mxparm) then
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

   80 return

c     end of  e04ucq. (npkey)

99999 format (' xxx  keyword not recognized:         ',a)
99998 format (' xxx  second keyword not recognized:  ',a)
99997 format (' xxx  the parm subscript is out of range:',i10)
99996 format (' xxx  ambiguous keyword:              ',a)
      end

      subroutine e04nch(linobj,rowerr,unitq,nclin,nactiv,nfree,nrank,nz,
     *                  n,nctotl,ldzy,lda,ldr,ldt,istate,kactiv,kx,jmax,
     *                  errmax,ctx,xnorm,a,ax,bl,bu,cq,res,res0,featol,
     *                  r,t,x,zy,p,work)
c----------------------------------------------------------------------
c     e04nch  computes the point on a working set that is closest to the
c     input vector  x  (in the least-squares sense).  the norm of  x,
c     the transformed residual vector  pr - rq'x,  and the constraint
c     values ax  are also initialized.
c
c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable  rowerr
c     is set.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      integer           ntry
      parameter         (ntry=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  ctx, errmax, xnorm
      integer           jmax, lda, ldr, ldt, ldzy, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nz
      logical           linobj, rowerr, unitq
c     .. array arguments ..
      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl), cq(*),
     *                  featol(nctotl), p(n), r(ldr,*), res(*), res0(*),
     *                  t(ldt,*), work(nctotl), x(n), zy(ldzy,*)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lsdbg
c     .. arrays in common ..
      integer           ilsdbg(ldbg)
c     .. local scalars ..
      double precision  bnd
      integer           i, is, j, k, ktry

      double precision  ddot, dnrm2
      integer           idamax
      external          ddot, dnrm2, idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ce04nc/ilsdbg, lsdbg
c----------------------------------------------------------------------
c     move  x  onto the simple bounds in the working set.

      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = bl(j)
         if (is.ge.2) bnd = bu(j)
         if (is.ne.4) x(j) = bnd
   20 continue

c     move  x  onto the general constraints in the working set.
c     we shall make  ntry  tries at getting acceptable row errors.

      ktry = 1
      jmax = 1
      errmax = zero

   40 if (nactiv.gt.0) then

c        set  work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.

         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(i) = bnd - ddot(n,a(k,1),lda,x,1)
   60    continue

         call e04nbt(1,ldt,nactiv,t(1,nz+1),work)
         call sload (n,zero,p,1)
         call dcopy (nactiv,work,1,p(nz+1),1)

         call e04nbw(2,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
         call daxpy (n,one,p,1,x,1)
      end if

c     compute the 2-norm of  x.
c     initialize  ax  for all the general constraints.

      xnorm = dnrm2(n,x,1)
      if (nclin.gt.0) call dgemv ('n',nclin,n,one,a,lda,x,1,zero,ax,1)

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
         call dcopy (n,x,1,p,1)
         call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,p,zy,work)
      end if
c
      ctx = zero
      if (linobj) ctx = ddot(n,cq,1,p,1)
c
      if (nrank.gt.0) then
         call dtrmv ('u','n','n',nrank,r,ldr,p,1)
         if (nrank.lt.n) call dgemv ('n',nrank,n-nrank,one,r(1,nrank+1),
     *                              ldr,p(nrank+1),1,one,p,1)
         call dcopy (nrank,res0,1,res,1)
         call daxpy (nrank,-one,p,1,res,1)
      end if

c     end of  e04nch. (lssetx)

99999 format (/' //e04nch// variables after refinement ... ')
99998 format (5g12.3)
      end

      subroutine e04ncy(unitq,vertex,inform,k1,k2,nactiv,nartif,nz,
     *                  nfree,nrank,nrejtd,nres,ngq,n,ldzy,lda,ldr,ldt,
     *                  istate,kactiv,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)
c----------------------------------------------------------------------
c     e04ncy  includes general constraints k1 thru k2 as new rows of
c     the tq factorization stored in t, zy.  if nrank is nonzero, the
c     changes in q are reflected in nrank by n triangular factor r such
c     that
c                         c  =  p ( r ) q,
c                                 ( 0 )
c     where  p  is orthogonal.
c----------------------------------------------------------------------
      implicit none

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  condmx
      integer           inform, k1, k2, lda, ldr, ldt, ldzy, msglvl, n,
     *                  nactiv, nartif, nfree, ngq, nrank, nrejtd, nres,
     *                  nz
      logical           unitq, vertex
c     .. array arguments ..
      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer           istate(*), kactiv(n), kx(n)
c     .. scalars in common ..
      double precision  asize, dtmax, dtmin
c     .. arrays in common ..

c     .. local scalars ..
      double precision  cndmax, rnorm, rowmax, rtmax
      integer           i, iadd, iartif, ifix, iswap, jadd, k, l, nzadd
c     .. external functions ..
      double precision  dnrm2
      external          dnrm2
c     .. common blocks ..
      double precision wmach
      common/ cstmch /wmach(10)

      common            /de04nb/asize, dtmax, dtmin
c----------------------------------------------------------------------
      rtmax = wmach(8)
c
c     estimate the condition number of the constraints that are not
c     to be refactorized.
c
      if (nactiv.eq.0) then
         dtmax = zero
         dtmin = one
      else
         call f06flf(nactiv,t(nactiv,nz+1),ldt-1,dtmax,dtmin)
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

      nrejtd = k2 - nactiv

c     end of  e04ncy. (lsadds)

      end

      subroutine e04mfz(prbtyp,msg,cset,named,names,rset,unitq,iter,
     *                  itmax,jinf,nviol,n,nclin,lda,nactiv,nfree,nrz,
     *                  nz,istate,kactiv,kx,e04mfu,obj,numinf,xnorm,a,
     *                  ax,bl,bu,cvec,featol,featlu,x,iw,w)
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
c
c     over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     after kdegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.

c     values of istate(j) for the linear constraints.......

c     istate(j)

c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.

c     constraint j may be violated by as much as featol(j).
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      integer           lenlc
      parameter         (lenlc=20)
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      character*6       empty
      parameter         (empty='      ')
c     .. scalar arguments ..
      double precision  obj, xnorm
      integer           iter, itmax, jinf, lda, n, nactiv, nclin, nfree,
     *                  nrz, numinf, nviol, nz
      logical           cset, named, rset, unitq
      character*2       prbtyp
      character*6       msg
c     .. array arguments ..
      double precision  a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  cvec(*), featlu(n+nclin), featol(n+nclin), w(*),
     *                  x(n)
      integer           istate(n+nclin), iw(*), kactiv(n), kx(n)
      character*8       names(*)
c     .. subroutine arguments ..
      external          e04mfu
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
c     .. local scalars ..
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
c     .. local arrays ..
      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*80      rec(3)
c     .. external functions ..
      double precision  ddot, dnrm2, adivb
      external          ddot, dnrm2, adivb
c     .. common blocks ..
      common            /ae04mf/loclc
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04nb/lennam, ldt, ncolt, ldq
      common            /ce04mf/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04mf/alfa, trulam, isdel, jdel, jadd, header,
     *                  prnt
      common            /de04nb/asize, dtmax, dtmin
      common            /ee04mf/ilcdbg, lcdbg
      common            /fe04mf/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /fe04nb/icmdbg, cmdbg
      common            /ge04mf/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)
c     .. save statement ..
      save              /ae04mf/, /fe04mf/, /ge04mf/, firstv
c----------------------------------------------------------------------

c     specify the machine-dependent parameters.

      epsmch = wmach(3)
      flmax = wmach(7)
      rtmax = wmach(8)
c
      if (cset) then
         ngq = 2
      else
         ngq = 1
      end if
c
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
c
c     we need a temporary array when changing the active set.
c     use the multiplier array.

      ltmp = lrlam

      if (iter.eq.0) then

c        first entry.  initialize.

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

      call e04mfh(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,bu,a,
     *            featol,w(lgq),x,w(lwtinf))

      if (numinf.gt.0) then
         call e04nbw(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),w(lwrk))
      else if (lp) then
         call dcopy (n,w(lcq),1,w(lgq),1)
      end if

      if (numinf.eq.0 .and. lp) then
         obj = ddot(n,cvec,1,x,1)
      else
         obj = suminf
      end if

      msg = empty

   20 if (msg.eq.empty) then

         gznorm = zero
         if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)

         if (nrz.eq.nz) then
            grznrm = gznorm
         else
            grznrm = zero
            if (nrz.gt.0) grznrm = dnrm2(nrz,w(lgq),1)
         end if

         gfnorm = gznorm
         if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),
     *       1)

c        print the details of this iteration.

c        define small quantities that reflect the size of x, r and
c        the constraints in the working set.

         if (prnt) then
            condt = one
            if (nactiv.gt.0) condt = adivb(dtmax,dtmin,overfl)

            call e04mfu(prbtyp,header,rset,msglvl,iter,isdel,jdel,jadd,
     *                  n,nclin,nactiv,nfree,nz,nrz,ldr,ldt,istate,alfa,
     *                  condrz,condt,dzz,gfnorm,gznorm,numinf,suminf,
     *                  notopt,obj,trulam,ax,w(lr),w(lt),x,w(lwrk))
            jdel = 0
            jadd = 0
            alfa = zero
         end if

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
            call sload (nactiv+nfixed,zero,w(lrlam),1)
            go to 20
         end if
c
         if (grznrm.le.dinky) then

c           the point  x  is a constrained stationary point.
c           compute lagrange multipliers.

c           define what we mean by 'tiny' and non-optimal multipliers.
c
            notopt = 0
            jdel = 0
            zerolm = -dinky
            smllst = -dinky
            biggst = dinky + one
            tinyst = dinky
c
            call e04mfm(prbtyp,msglvl,n,lda,ldt,nactiv,nfree,nz,istate,
     *                  kactiv,kx,zerolm,notopt,numinf,trusml,smllst,
     *                  jsmlst,ksmlst,tinyst,jtiny,jinf,trubig,biggst,
     *                  jbigst,kbigst,a,w(lanorm),w(lgq),w(lrlam),w(lt),
     *                  w(lwtinf))
c
            if (nrz.lt.nz) call e04mfl(msglvl,n,nrz,nz,zerolm,notopt,
     *                                 numinf,trusml,smllst,jsmlst,
     *                                 tinyst,jtiny,w(lgq))
c
            if (abs(jsmlst).gt.0) then

c              delete a constraint.

c              e04mfm  or  e04mfl  found a non-optimal multiplier.
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
            end if

c           constraint  jdel  has been deleted.
c           update the  tq  factorization.

            call e04nfp(unitq,it,n,nactiv,nfree,ngq,nz,nrz,lda,ldq,ldt,
     *                  jdel,kdel,kactiv,kx,a,w(lt),w(lgq),w(lq),w(lwrk)
     *                  ,w(ld),w(lrlam))
            if (rset) call e04mfy(nrz,ldr,w(lr),one)
c
            prnt = .false.
         else
c           compute a search direction.

            if (iter.ge.itmax) then
               msg = 'itnlim'
               go to 20
            end if

            prnt = .true.
            iter = iter + 1

            call dcopy (nrz,w(lgq),1,w(ld),1)
            call dscal (nrz,(-one),w(ld),1)

            dnorm = dnrm2(nrz,w(ld),1)

            call e04nbw(1,n,nrz,nfree,ldq,unitq,kx,w(ld),w(lq),w(lwrk))
            call dgemv ('no transpose',nclin,n,one,a,lda,w(ld),1,zero,
     *                 w(lad),1)

c           find the constraint we bump into along d.
c           update  x  and  ax  if the step alfa is nonzero.

c           alfhit is initialized to bigalf. if it remains that value
c           after the call to  e04mfs, it is regarded as infinite.

            bigalf = adivb(bigdx,dnorm,overfl)

            call e04mfs(firstv,n,nclin,istate,bigalf,bigbnd,dnorm,
     *                  hitlow,move,onbnd,unbndd,alfhit,alfap,jadd,
     *                  w(lanorm),w(lad),ax,bl,bu,featol,featlu,w(ld),x)

            if (unbndd) then
               msg = 'unbndd'
               go to 20
            end if

            alfa = alfhit
            call daxpy (n,alfa,w(ld),1,x,1)

            if (nclin.gt.0) call daxpy (nclin,alfa,w(lad),1,ax,1)
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

            call e04nfr(unitq,rset,inform,ifix,iadd,jadd,it,nactiv,nz,
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

            call daxpy (nctotl,tolinc,featlu,1,featol,1)

            if (mod(iter,kchk).eq.0) then

c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  set inform to force iterative
c              refinement and a switch to phase 1.

               call e04mfq(n,nclin,istate,bigbnd,nviol,jmax,errmax,ax,
     *                     bl,bu,featol,x)

               if (nviol.gt.0) then
                  if (msglvl.gt.0) then
                     write (rec,fmt=99999) errmax, jmax
                     call x04baf(iprint,rec(1))
                  end if
               end if
            end if

            if (mod(iter,kdegen).eq.0) then

c              every  kdegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.

               call e04mfr('end of cycle',msglvl,n,nclin,nmoved,iter,
     *                     numinf,istate,bigbnd,ax,bl,bu,featol,featlu,
     *                     x)

               nviol = nviol + nmoved
            end if

            if (nviol.gt.0) then
               msg = 'resetx'
               go to 20
            end if

            if (numinf.ne.0) then
               call e04mfh(n,nclin,lda,istate,bigbnd,numinf,suminf,bl,
     *                     bu,a,featol,w(lgq),x,w(lwtinf))

               if (numinf.gt.0) then
                  call e04nbw(6,n,nz,nfree,ldq,unitq,kx,w(lgq),w(lq),
     *                        w(lwrk))
               else if (lp) then
                  call dcopy (n,w(lcq),1,w(lgq),1)
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

      end if

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
c DEBUG 691                       random assignment (null)
      msglvl = 0

c     end of  e04mfz.  (lpcore)

99999 format (' xxx  iterative refinement.  the max violation is ',1p,
     *       d10.2,' in constraint',i5)
99998 format (/' //e04mfz//       grznrm      dinky',/' //e04mfz//  ',
     *       1p,2d11.2)
      end

      subroutine e04ucp(n,nclin,ncnln,nctotl,litotl,lwtotl)
c----------------------------------------------------------------------
c     e04ucp   allocates the addresses of the work arrays for e04ucz and
c     e04ncz.
c----------------------------------------------------------------------
      implicit none

      integer           lenls
      parameter         (lenls=20)
      integer           lennp
      parameter         (lennp=35)
      integer           ldbg
      parameter         (ldbg=5)
c     .. scalar arguments ..
      integer           litotl, lwtotl, n, nclin, ncnln, nctotl
c     .. scalars in common ..
      integer           ldt, ldzy, lennam, ncolt
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg), locls(lenls), locnp(lennp)
c     .. local scalars ..
      integer           ladx, lanorm, laqp, lbl, lbu, lc1mul, lcjac,
     *                  lcjdx, lcmul, lcs1, lcs2, ldlam, ldslk, ldx,
     *                  lenaqp, lent, lenzy, lfeatl, lgq, lgq1, lgrad,
     *                  lhctrl, lhfrwd, liperm, lkactv, lkx, lneedc,
     *                  lqpadx, lqpdx, lqpgq, lqphz, lqptol, lrho,
     *                  lrlam, lrpq, lrpq0, lslk, lslk1, lt, lwrk1,
     *                  lwrk2, lwrk3, lwtinf, lx1, lzy, miniw, minw
c     .. common blocks ..
      common            /ae04nc/locls
      common            /ae04uc/locnp
      common            /be04nb/lennam, ldt, ncolt, ldzy
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
      miniw = litotl + 1
      minw = lwtotl + 1

c     assign array lengths that depend upon the problem dimensions.

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

c     assign the addresses for the workspace arrays used by  e04ucu.

      lqpadx = minw
      lqpdx = lqpadx + nclin + ncnln
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
      lbl = ladx + nclin + ncnln
      lbu = lbl + nctotl
      ldx = lbu + nctotl
      lgq1 = ldx + n
      lfeatl = lgq1 + n
      lx1 = lfeatl + nctotl
      lwrk2 = lx1 + n
      minw = lwrk2 + nctotl

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

      lcjac = minw
      lgrad = lcjac + ncnln*n
      minw = lgrad + n

      locnp(25) = lhfrwd
      locnp(26) = lhctrl
      locnp(27) = lcjac
      locnp(28) = lgrad
c
      litotl = miniw - 1
      lwtotl = minw - 1

c     end of  e04ucp. (nploc)

      end

      subroutine e04ucz(named,names,unitq,inform,majits,n,nclin,ncnln,
     *                  nctotl,nactiv,nfree,nz,ldcj,ldcju,ldaqp,ldr,
     *                  nfun,ngrad,istate,kactiv,kx,objf,fdnorm,xnorm,
     *                  confun,objfun,aqp,ax,bl,bu,c,cjac,cjacu,clamda,
     *                  featol,grad,gradu,r,x,iw,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     e04ucz  is the core routine for  nlpsol,  a sequential quadratic
c     programming (sqp) method for nonlinearly constrained optimization.
c----------------------------------------------------------------------
      implicit none

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
c     .. scalar arguments ..
      double precision  fdnorm, objf, xnorm
      integer           inform, ldaqp, ldcj, ldcju, ldr, lenw, majits,
     *                  n, nactiv, nclin, ncnln, nctotl, nfree, nfun,
     *                  ngrad, nz
      logical           named, unitq
c     .. array arguments ..
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
c     .. local scalars ..
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
     *                  mnrsum, mode, msgqp, msgsv1, msgsv2, ncqp, nl,
     *                  nlnact, nlserr, nmajor, nminor, nplin, nqperr,
     *                  nqpinf, nstate, numinf, nviol
      logical           centrl, convpt, convrg, done, error, feasqp,
     *                  goodgq, infeas, needfd, newgq, optiml, overfl
      character*5       mjrmsg
c     .. local arrays ..
      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      logical           ktcond(2)
      character*80      rec(5)
c     .. external functions ..
      double precision  ddot, dnrm2, adivb
      external          ddot, dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls
      common            /ae04uc/locnp

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04nb/lennam, ldt, ncolt, ldzy
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
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c     .. save statement ..
      save              /de04nc/, /ee04nc/, /ge04uc/, /he04uc/
c----------------------------------------------------------------------
c     specify machine-dependent parameters.

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

c     initialize

      mjrmsg = '     '
      nqpinf = 0
      mnrsum = 0

      majit0 = majits
      nplin = n + nclin
      ncqp = nclin + ncnln
      nl = min(nplin+1,nctotl)

      ldcj1 = max(ncqp,1)

      needfd = lvlder .eq. 0 .or. lvlder .eq. 2 .or.
     *         (lvlder.eq.1 .and. ncnln.gt.0)

      alfa = zero
      alfdx = zero
      rtftol = sqrt(ftol)
      rootn = sqrt(dble(n))

c     information from the feasibility phase will be used to generate a
c     hot start for the first qp subproblem.

      call dcopy (nctotl,featol,1,w(lqptol),1)

      nstate = 0

      objalf = objf
      if (ncnln.gt.0) then
         objalf = objalf - ddot(ncnln,w(lcmul),1,c,1)
      end if
c
      newgq = .false.

   20 minits = 0
c
   40 centrl = lvldif .eq. 2
c
      if (newgq) then
         if (needfd) then

c           compute any missing gradient elements and the
c           transformed gradient of the objective.

            call e04uds(centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,iw(lneedc),bl,
     *                  bu,c,w(lwrk2),w(lwrk3),cjac,cjacu,grad,gradu,
     *                  w(lhfrwd),w(lhctrl),x,w,lenw,iuser,user)
            inform = mode
            if (mode.lt.0) go to 60
c
         end if
c
         call dcopy (n,grad,1,w(lgq),1)
         call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),w(lwrk1))
         newgq = .false.
      end if

c     (1) solve an inequality quadratic program (iqp) for the
c         search direction and multiplier estimates.
c     (2) for each nonlinear inequality constraint,  compute
c         the slack variable for which the merit function is
c         minimized.
c     (3) compute the search direction for the slack variables
c         and multipliers.

c     the array violn is wrk3.

      call e04ucu(feasqp,unitq,nqperr,majits,mnr,n,nclin,ncnln,ldcj,
     *            ldaqp,ldr,linact,nlnact,nactiv,nfree,nz,numinf,istate,
     *            kactiv,kx,dxnorm,gdx,qpcurv,aqp,w(ladx),w(lanorm),ax,
     *            bl,bu,c,cjac,clamda,w(lcmul),w(lcs1),w(ldlam),w(ldslk)
     *            ,w(ldx),w(lbl),w(lbu),w(lqptol),r,w(lrho),w(lslk),
     *            w(lvioln),x,w(lwtinf),w)

      minits = minits + mnr
      mnrsum = mnrsum + mnr

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

      end if

      if ( .not. goodgq) go to 40

c     (1) compute the number of constraints that are violated by more
c         than featol.
c     (2) compute the 2-norm of the residuals of the constraints in
c         the qp working set.

      call e04ucw(n,nclin,ncnln,istate,bigbnd,cvnorm,errmax,jmax,nviol,
     *            ax,bl,bu,c,featol,x,w(lwrk2))

c     define small quantities that reflect the magnitude of objf and
c     the norm of grad(free).

      objsiz = one + abs(objf)
      xsize = one + xnorm
      gtest = max(objsiz,gfnorm)
      dinky = rtftol*gtest

      if (nactiv.eq.0) then
         condt = zero
      else if (nactiv.eq.1) then
         condt = dtmin
      else
         condt = adivb(dtmax,dtmin,overfl)
      end if

      call f06flf(n,r,ldr+1,drmax,drmin)

      condh = adivb(drmax,drmin,overfl)
      if (condh.lt.rtmax) then
         condh = condh*condh
      else
         condh = flmax
      end if

      if (nz.eq.0) then
         condhz = one
      else if (nz.eq.n) then
         condhz = condh
      else
         call f06flf(nz,r,ldr+1,drzmax,drzmin)
         condhz = adivb(drzmax,drzmin,overfl)
         if (condhz.lt.rtmax) then
            condhz = condhz*condhz
         else
            condhz = flmax
         end if
      end if

c     test for convergence.
c     the point test convpt checks for a k-t point at the initial
c     point or after a large change in x.

      convpt = dxnorm .le. epspt8*gtest .and. nviol .eq. 0 .and.
     *         nqperr .le. 1

      ktcond(1) = gznorm .lt. dinky
      ktcond(2) = nviol .eq. 0
      optiml = ktcond(1) .and. ktcond(2)

      convrg = majits .gt. 0 .and. alfdx .le. rtftol*xsize

      infeas = convrg .and. .not. feasqp .or. nqpinf .gt. 7

      done = convpt .or. (convrg .and. optiml) .or. infeas

      objalf = objf
      grdalf = gdx
      glf1 = gdx
      if (ncnln.gt.0) then
         glf1 = glf1 - ddot(ncnln,w(lcjdx),1,clamda(nl),1)

c        compute the value and directional derivative of the
c        augmented lagrangian merit function.
c        the penalty parameters may be increased or decreased.

         call e04ucn(feasqp,n,nclin,ncnln,objalf,grdalf,qpcurv,istate,
     *               w(lcjdx),w(lcmul),w(lcs1),w(ldlam),w(lrho),
     *               w(lvioln),w(lwrk1),w(lwrk2))
      end if

c     print the details of this iteration.

      call e04uct(ktcond,convrg,mjrmsg,msgnp,msgqp,ldr,ldt,n,nclin,
     *            ncnln,nctotl,nactiv,linact,nlnact,nz,nfree,majit0,
     *            majits,minits,istate,alfa,nfun,condhz,condh,condt,
     *            objalf,objf,gfnorm,gznorm,cvnorm,ax,c,r,w(lt),
     *            w(lvioln),x,w(lwrk1))

      alfa = zero
      error = majits .ge. nmajor

      if ( .not. (done .or. error)) then
         majits = majits + 1

c        make copies of information needed for the bfgs update.

         call dcopy (n,x,1,w(lx1),1)
         call dcopy (n,w(lgq),1,w(lgq1),1)

         if (ncnln.gt.0) then
            call dcopy (ncnln,w(lcjdx),1,w(lcjdx1),1)
            call dcopy (ncnln,w(lcmul),1,w(lc1mul),1)
            call dcopy (ncnln,w(lslk),1,w(lslk1),1)
         end if

c        compute the parameters for the linesearch.

c        alfmin is the smallest allowable step predicted by the qp
c        subproblem.

         alfmin = one
         if ( .not. feasqp) alfmin = zero

c        alfmax is the largest feasible steplength subject to a user-
c        defined limit alflim on the change in x.

         if (ncnln.gt.0 .and. needfd) then
            alfmax = one
         else
            alfmax = adivb(bigdx,dxnorm,overfl)
            call e04udt(info,n,nclin,ncnln,alfa,alfmin,alfmax,bigbnd,
     *                  dxnorm,w(lanorm),w(ladx),ax,bl,bu,w(ldslk),
     *                  w(ldx),w(lslk),x)
            alfmax = alfa
            if (alfmax.lt.one+epspt3 .and. feasqp) alfmax = one
         end if

c        alfbnd is a tentative upper bound on the steplength.  if the
c        merit function is decreasing at alfbnd and certain
c        conditions hold,  alfbnd will be increased in multiples of
c        two (subject to not being greater than alfmax).

         if (ncnln.eq.0) then
            alfbnd = alfmax
         else
            alfbnd = min(one,alfmax)
         end if

c        alfsml trips the computation of central differences.  if a
c        trial steplength falls below alfsml, the linesearch is
c        terminated.

         alfsml = zero
         if (needfd .and. .not. centrl) then
            alfsml = adivb(fdnorm,dxnorm,overfl)
            alfsml = min(alfsml,alfmax)
         end if

c        compute the steplength using safeguarded interpolation.

         alflim = adivb((one+xnorm)*dxlim,dxnorm,overfl)
         alfa = min(alflim,one)

         call e04ucr(needfd,nlserr,n,ncnln,ldcj,ldcju,nfun,ngrad,
     *               iw(lneedc),confun,objfun,alfa,alfbnd,alfmax,alfsml,
     *               dxnorm,epsrf,eta,gdx,grdalf,glf1,glf2,objf,objalf,
     *               qpcurv,xnorm,c,w(lwrk1),cjac,cjacu,w(lcjdx),
     *               w(lwrk3),w(lc1mul),w(lcmul),w(lcs1),w(lcs2),w(ldx),
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

c              compute the missing gradients.

               mode = 1
               ngrad = ngrad + 1

               if (ncnln.gt.0) then
                  call f06dbf(ncnln,(1),iw(lneedc),1)

                  call confun(mode,ncnln,n,ldcju,iw(lneedc),x,w(lwrk1),
     *                        cjacu,nstate,iuser,user)
                  inform = mode
                  if (mode.lt.0) go to 60

                  call f06qff('g',ncnln,n,cjacu,ldcju,cjac,ldcj)
               end if

               call objfun(mode,n,x,obj,gradu,nstate,iuser,user)
               inform = mode
               if (mode.lt.0) go to 60

               call dcopy (n,gradu,1,grad,1)

               call e04uds(centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                     fdint,fdnorm,objf,confun,objfun,iw(lneedc),
     *                     bl,bu,c,w(lwrk2),w(lwrk3),cjac,cjacu,grad,
     *                     gradu,w(lhfrwd),w(lhctrl),x,w,lenw,iuser,
     *                     user)

               inform = mode
               if (mode.lt.0) go to 60

               gdx = ddot(n,grad,1,w(ldx),1)
               glf2 = gdx
               if (ncnln.gt.0) then
                  call dgemv ('n',ncnln,n,one,cjac,ldcj,w(ldx),1,zero,
     *                       w(lcjdx),1)
                  glf2 = glf2 - ddot(ncnln,w(lcjdx),1,clamda(nl),1)
               end if
            end if

            call dcopy (n,grad,1,w(lgq),1)
            call e04nbw(6,n,nz,nfree,ldzy,unitq,kx,w(lgq),w(lzy),
     *                  w(lwrk1))

            xnorm = dnrm2(n,x,1)

            if (ncnln.gt.0 .and. alfa.ge.one) call dcopy (ncnln,
     *          clamda(nl),1,w(lcmul),1)

            if (nclin.gt.0) call daxpy (nclin,alfa,w(ladx),1,ax,1)
            alfdx = alfa*dxnorm

c           update the factors of the approximate hessian of the
c           lagrangian function.

            call e04ucl(mjrmsg,unitq,n,ncnln,nfree,nz,ldcj1,ldcj,ldzy,
     *                  ldr,kx,alfa,glf1,glf2,qpcurv,w(lcjac1),cjac,
     *                  w(lcjdx1),w(lcjdx),w(lcs1),w(lcs2),w(lgq1),
     *                  w(lgq),w(lhpq),w(lrpq),clamda(nl),r,w(lwrk3),
     *                  w(lzy),w(lwrk2),w(lwrk1))

            call f06flf(n,r,ldr+1,drmax,drmin)
            cond = adivb(drmax,drmin,overfl)

            if (cond.gt.rcndbd .or. rfrobn.gt.rootn*growth*drmax) then

c              reset the condition estimator and range-space
c              partition of q'hq.

               mjrmsg(5:5) = 'refactorize hessian'

               call e04udr(unitq,n,nfree,nz,ldzy,ldr,iw(liperm),kx,
     *                     w(lgq),r,w(lzy),w(lwrk1),w(lqrwrk))
            end if
         end if
      end if

      if ( .not. (done .or. error)) go to 20

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

c     set  clamda.  print the full solution.

c DEBUG 691                       random assignment (null)

   60 if (msgnp.gt.0) then
         write (rec,fmt=99998) majits, mnrsum
         call x04bay(iprint,3,rec)
      end if
c
      call e04nbx(msgnp,nfree,ldaqp,n,nclin,nctotl,bigbnd,named,names,
     *            nactiv,istate,kactiv,kx,aqp,bl,bu,c,clamda,w(lrlam),x)
      if (ncnln.gt.0) call dcopy (ncnln,w(lcmul),1,clamda(n+nclin+1),1)

c     end of  e04ucz. (npcore)

99999 format (' mnr itn ',i4,' -- re-solve qp subproblem.')
99998 format (/' exit from np problem after ',i5,' major iterations,',
     *       /'                            ',i5,' minor iterations.')
99997 format (/' //e04ucz//        rfrobn         drmax         drmin',
     *       /' //e04ucz//',1p,3d14.2,/' //e04ucz//          cond     ',
     *       '   rcndbd',/' //e04ucz//',1p,2d14.2)
      end

      subroutine e04mfj(rowerr,unitq,nclin,nactiv,nfree,nz,n,ldq,lda,
     *                  ldt,istate,kactiv,kx,jmax,errmax,xnorm,a,ax,bl,
     *                  bu,featol,t,x,q,p,work)
c-----------------------------------------------------------------------
c     e04mfj  computes the point on a working set that is closest in the
c     least-squares sense to the input vector x.

c     if the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  if  x  is still infeasible,  the logical variable rowerr
c     is set.
c-----------------------------------------------------------------------
      implicit none

      integer           ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
      integer           ntry
      parameter         (ntry=5)
c     .. scalar arguments ..
      double precision  errmax, xnorm
      integer           jmax, lda, ldq, ldt, n, nactiv, nclin, nfree, nz
      logical           rowerr, unitq
c     .. array arguments ..
      double precision  a(lda,*), ax(*), bl(n+nclin), bu(n+nclin),
     *                  featol(n+nclin), p(n), q(ldq,*), t(ldt,*),
     *                  work(n), x(n)
      integer           istate(n+nclin), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)
c     .. local scalars ..
      double precision  bnd
      integer           i, is, j, k, ktry
c     .. external functions ..
      double precision  ddot, dnrm2
      integer           idamax
      external          ddot, dnrm2, idamax
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------

c     move  x  onto the simple bounds in the working set.

      do 20 k = nfree + 1, n
         j = kx(k)
         is = istate(j)
         bnd = bl(j)
         if (is.ge.2) bnd = bu(j)
         if (is.ne.4) x(j) = bnd
   20 continue

c     move  x  onto the general constraints in the working set.
c     ntry  attempts are made to get acceptable row errors.

      ktry = 1
      jmax = 1
      errmax = zero

   40 if (nactiv.gt.0) then

c        set work = residuals for constraints in the working set.
c        solve for p, the smallest correction to x that gives a point
c        on the constraints in the working set.  define  p = y*(py),
c        where  py  solves the triangular system  t*(py) = residuals.

         do 60 i = 1, nactiv
            k = kactiv(i)
            j = n + k
            bnd = bl(j)
            if (istate(j).eq.2) bnd = bu(j)
            work(nactiv-i+1) = bnd - ddot(n,a(k,1),lda,x,1)
   60    continue
c
         call dtrsv ('u','n','n',nactiv,t(1,nz+1),ldt,work,1)
         call sload (n,zero,p,1)
         call dcopy (nactiv,work,1,p(nz+1),1)

         call e04nbw(2,n,nz,nfree,ldq,unitq,kx,p,q,work)
         call daxpy (n,one,p,1,x,1)
      end if

c     initialize  ax  for all the general constraints.

      xnorm = dnrm2(n,x,1)
      if (nclin.gt.0) call dgemv ('n',nclin,n,one,a,lda,x,1,zero,ax,1)

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

      ktry = ktry + 1

      if ( .not. (errmax.le.featol(jmax) .or. ktry.gt.ntry)) go to 40

      rowerr = errmax .gt. featol(jmax)

c     end of  e04mfj.  (cmsetx)

      end

      subroutine e04mfv(cset,n,nclin,litotl,lwtotl)
c----------------------------------------------------------------------
c     e04mfv   allocates the addresses of the work arrays for e04mfz.

c     note that the arrays ( gq, cq ) lie in contiguous areas of
c     workspace.

      integer           lenlc
      parameter         (lenlc=20)
      integer           ldbg
      parameter         (ldbg=5)
c     .. scalar arguments ..
      integer           litotl, lwtotl, n, nclin
      logical           cset
c     .. scalars in common ..
      integer           iprint, isumm, ldq, ldt, lennam, lines1, lines2,
     *                  ncolt, nout
      logical           lcdbg
c     .. arrays in common ..
      integer           ilcdbg(ldbg), loclc(lenlc)
c     .. local scalars ..
      integer           lad, lanorm, lcq, ld, lencq, lenq, lenrt,
     *                  lfeatu, lgq, lkactv, lkx, lq, lr, lrlam, lt,
     *                  lwrk, lwtinf, miniw, minw
c     .. common blocks ..
      common            /ae04mf/loclc
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /be04nb/lennam, ldt, ncolt, ldq
      common            /ee04mf/ilcdbg, lcdbg
c     .. save statement ..
      save              /ae04mf/
c----------------------------------------------------------------------
c     refer to the first free space in the work arrays.

      miniw = litotl + 1
      minw = lwtotl + 1

c     integer workspace.

      lkactv = miniw
      lkx = lkactv + n
      miniw = lkx + n

c     real workspace.
c     assign array lengths that depend upon the problem dimensions.

      lenrt = ldt*ncolt
      if (nclin.eq.0) then
         lenq = 0
      else
         lenq = ldq*ldq
      end if

      if (cset) then
         lencq = n
      else
         lencq = 0
      end if

c     we start with arrays that can be preloaded by smart users.

      lfeatu = minw
      minw = lfeatu + nclin + n

c     next comes stuff used by  e04mfz  and  e04nfz.

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

c     load the addresses in loclc.

      loclc(1) = lkactv
      loclc(2) = lkx

      loclc(3) = lfeatu
      loclc(4) = lanorm
      loclc(5) = lad

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

c     end of  e04mfv.  (lploc)

      end

      subroutine e04mfw(n,nclin,title)
c----------------------------------------------------------------------
c     e04mfw loads the default values of parameters not set by the user.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      integer           mxparm
      parameter         (mxparm=30)
      double precision  zero
      parameter         (zero=0.0d+0)
      double precision  rdummy
      integer           idummy
      parameter         (rdummy=-11111.0d+0,idummy=-11111)
      double precision  gigant
      parameter         (gigant=1.0d+20*0.99999d+0)
      double precision  wrktol
      parameter         (wrktol=1.0d-2)
c     .. scalar arguments ..
      integer           n, nclin
      character*(*)     title
c     .. scalars in common ..
      double precision  bigbnd, bigdx, bndlow, bndupp, epspt3, epspt5,
     *                  epspt8, epspt9, tolact, tolfea, tolinc, tolrnk,
     *                  tolx0
      integer           idbglc, iprint, iprnt, isumm, isumry, itmax1,
     *                  itmax2, itnfix, kchk, kcycle, kdegen, lcrash,
     *                  ldbglc, lines1, lines2, lprob, maxact, maxnz,
     *                  mm, msglc, mxfree, ndegen, nn, nnclin, nout,
     *                  nprob
      logical           cmdbg, lcdbg, newopt
c     .. arrays in common ..
      double precision  rpadlc(23), rpsvlc(mxparm)
      integer           icmdbg(ldbg), ilcdbg(ldbg), ipadlc(14),
     *                  ipsvlc(mxparm), nfix(2)
c     .. local scalars ..
      double precision  epsmch
      integer           i, idbg, j, k, lent, msgdbg, msglvl
      character*16      key
c     .. local arrays ..
      double precision  rprmlc(mxparm)
      integer           iprmlc(mxparm)
      character*4       icrsh(0:2)
      character*7       lptype(1:10)
      character*80      rec(3)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04mf/newopt
      common            /ce04mf/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ee04mf/ilcdbg, lcdbg
      common            /fe04mf/ipsvlc, idbglc, iprnt, isumry, itmax1,
     *                  itmax2, kchk, kcycle, lcrash, lprob, maxact,
     *                  mxfree, maxnz, mm, ldbglc, msglc, nn, nnclin,
     *                  nprob, ipadlc
      common            /fe04nb/icmdbg, cmdbg
      common            /ge04mf/rpsvlc, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadlc
c     .. equivalences ..
      equivalence       (iprmlc(1),idbglc), (rprmlc(1),bigbnd)
      equivalence       (msglc,msglvl), (idbglc,idbg), (ldbglc,msgdbg)
c     .. save statement ..
      save              /be04mf/, /fe04mf/, /ge04mf/
c     .. data statements ..
      data              icrsh(0), icrsh(1), icrsh(2)/'cold', 'warm',
     *                  'hot '/
      data              lptype(1), lptype(2)/'     fp', '     lp'/
      data              lptype(3), lptype(4), lptype(5),
     *                  lptype(6)/'illegal', 'illegal', 'illegal',
     *                  'illegal'/
      data              lptype(7), lptype(8), lptype(9),
     *                  lptype(10)/'       ', '       ', '       ',
     *                  'illegal'/
c----------------------------------------------------------------------
      epsmch = wmach(3)

c     make a dummy call to e04mfx to ensure that the defaults are set.

      call e04mfx (nout,'*',key)
      newopt = .true.

c     save the optional parameters set by the user.  the values in
c     rprmlc and iprmlc may be changed to their default values.

      call f06dff (mxparm,iprmlc,1,ipsvlc,1)
      call dcopy (mxparm,rprmlc,1,rpsvlc,1)

      if (msglvl.eq.idummy) msglvl = 10
      if (iprnt.lt.0) iprnt = nout
      if (isumry.lt.0 .or. msglvl.lt.5) isumry = -1
      iprint = iprnt
      isumm = isumry
      if (kchk.le.0) kchk = 50
      if (kcycle.le.0) kcycle = 10000
      if (kcycle.gt.9999999) kcycle = 9999999
      kdegen = kcycle
      if (lprob.lt.0) lprob = 2
      if (lcrash.lt.0 .or. lcrash.gt.2) lcrash = 0
      if (itmax1.lt.0) itmax1 = max(50,5*(n+nclin))
      if (itmax2.lt.0) itmax2 = max(50,5*(n+nclin))
      if (maxact.lt.0 .or. maxact.gt.n .or. maxact.gt.nclin)
     *    maxact = max(1,min(n,nclin))
      if (maxnz.lt.0 .or. maxnz.gt.n) maxnz = n
      if (mxfree.lt.0 .or. mxfree.gt.n) mxfree = n
      if (mxfree.lt.maxnz) mxfree = maxnz
      if (nclin.lt.n) then
         mxfree = nclin + 1
         maxnz = mxfree
      end if

      if (tolact.lt.zero) tolact = wrktol
      if (tolfea.eq.rdummy .or. (tolfea.ge.zero .and. tolfea.lt.epsmch))
     *    tolfea = epspt5
      if (bigbnd.le.zero) bigbnd = gigant
      if (bigdx.le.zero) bigdx = max(gigant,bigbnd)

      k = 1

      if (msglvl.gt.0) then

c        print the title.

         lent = len(title)
         write (rec,fmt=99993) (title(j:j),j=1,lent)
         call x04bay(iprint,2,rec)

         if (msglvl.ge.5 .and. isumm.ge.0 .and. isumm.ne.iprint) then
            write (rec,fmt=99992) (title(j:j),j=1,lent)
            call x04bay(isumm,2,rec)
         end if
c
         write (rec,fmt=99999)
         call x04bay(iprint,3,rec)
         write (rec,fmt=99998) lptype(lprob)
         call x04bay(iprint,2,rec)
         write (rec,fmt=99997) nclin, tolfea, n, tolact
         call x04bay(iprint,3,rec)
         write (rec,fmt=99996) bigbnd, icrsh(lcrash), bigdx, epsmch
         call x04bay(iprint,3,rec)
         write (rec,fmt=99995) kchk, kdegen
         call x04bay(iprint,2,rec)
         write (rec,fmt=99994) msglvl, itmax2, isumry
         call x04bay(iprint,3,rec)
      end if

c     end of  e04mfw.  (lpdflt)

99999 format (/' parameters',/' -')
99998 format (/' problem type...........',3x,a7)
99997 format (/' linear constraints.....',i10,7x,'feasibility toleranc',
     *       'e..',1p,d10.2,/' variables..............',i10,7x,'crash ',
     *       'tolerance........',1p,d10.2)
99996 format (/' infinite bound size....',1p,d10.2,7x,a4,' start......',
     *       '.......',/' infinite step size.....',1p,d10.2,7x,'eps (m',
     *       'achine precision)',1p,d10.2)
99995 format (/' check frequency........',i10,7x,'expand frequency....',
     *       '...',i10)
99994 format (/' print level............',i10,7x,'iteration limit.....',
     *       '...',i10,/' monitoring file........',i10)
99993 format (/80a1)
99992 format (/11a1,' monitoring information ')
      end

      subroutine e04nfq(unitq,vertex,k1,k2,it,nactiv,nartif,nz,nfree,
     *                  nrejtd,ngq,n,ldq,lda,ldt,istate,kactiv,kx,
     *                  condmx,a,t,gqm,q,w,c,s,msglvl)
c----------------------------------------------------------------------
c     e04nfq  includes general constraints  k1  thru  k2  as new rows of
c     the  tq  factorization:
c              a(free) * q(free)  = (  0 t )
c                        q(free)  = (  z y )
c
c     a) the  nactiv x nactiv  upper-triangular matrix  t  is stored
c        with its (1,1) element in position  (it,jt)  of the array  t.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  condmx
      integer           it, k1, k2, lda, ldq, ldt, msglvl, n, nactiv,
     *                  nartif, nfree, ngq, nrejtd, nz
      logical           unitq, vertex
c     .. array arguments ..
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
c     .. local scalars ..
      double precision  cndmax, cond, delta, drzz, dtnew, rnorm, rowmax,
     *                  rtmax, tdtmax, tdtmin
      integer           i, iadd, iartif, ifix, inform, iswap, j, jadd,
     *                  jt, k, l, nzadd
      logical           overfl, rset
c     .. external functions ..
      double precision  dnrm2, adivb
      external          dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /de04nb/asize, dtmax, dtmin
      common            /fe04nb/icmdbg, cmdbg
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

            call f06qhf('g',nfree,nfree,zero,one,q,ldq)
            unitq = .false.
         end if
      else
         call f06flf(nactiv,t(it,jt),ldt+1,dtmax,dtmin)
      end if
c
      do 20 k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv.lt.nfree) then

            overfl = .false.

c           transform the incoming row of  a  by  q'.

            call dcopy (n,a(iadd,1),lda,w,1)
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
c
c                 note that delta  has to be stored after the reflection
c                 is used.
c
                  delta = w(nz)
                  call f06frf(nz-1,delta,w,1,zero,w(nz))
                  if (w(nz).gt.zero) then
c
                     call dgemv ('n',nfree,nz,one,q,ldq,w,1,zero,s,1)
                     call dger(nfree,nz,(-one),s,1,w,1,q,ldq)
c
                     if (ngq.gt.0) then
                        call dgemv ('t',nz,ngq,one,gqm,n,w,1,zero,s,1)
                        call dger(nz,ngq,(-one),w,1,s,1,gqm,n)
                     end if
                  end if

                  w(nz) = delta
               end if
               it = it - 1
               jt = jt - 1
               nactiv = nactiv + 1
               nz = nz - 1
               call dcopy (nactiv,w(jt),1,t(it,jt),ldt)
               dtmax = tdtmax
               dtmin = tdtmin
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

c        if a vertex is required,  add some temporary bounds.
c        we must accept the resulting condition number of the working
c        set.

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

                  call e04nfr(unitq,rset,inform,ifix,iadd,jadd,it,
     *                        nactiv,nz,nfree,nz,ngq,n,lda,ldq,ldt,ldt,
     *                        kx,cndmax,drzz,a,t,t,gqm,q,w,c,s,msglvl)
               end if
               nfree = nfree - 1
               nz = nz - 1
               nartif = nartif + 1
               istate(jadd) = 4
   80       continue
         end if

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

c     end of  e04nfq.  (rzadds)

99999 format (/' //e04nfq//  constraint added.           ',/' //e04nfq',
     *       '//  nactiv    nz nfree  iadd  jadd',/' //e04nfq//  ',5i6)
99998 format (/' //e04nfq//  dependent constraint rejected.')
99997 format (/' //e04nfq//     asize     dtmax     dtmin     dtnew',
     *       /' //e04nfq//',1p,4d10.2)
99996 format (/' //e04nfq//     asize     dtnew',/' //e04nfq//',1p,
     *       2d10.2)
      end

      subroutine e04ucx(n,nclin,ncnln,title)
c----------------------------------------------------------------------
c     e04ucx  loads the default values of parameters not set in the
c     options file.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
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
c     .. scalar arguments ..
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
      logical           cmdbg, lsdbg, newopt, npdbg
c     .. arrays in common ..
      double precision  rpadls(23), rpadnp(22), rpsvls(mxparm),
     *                  rpsvnp(mxparm)
      integer           icmdbg(ldbg), ilsdbg(ldbg), inpdbg(ldbg),
     *                  ipadls(18), ipadnp(12), ipsvls(mxparm),
     *                  ipsvnp(mxparm), jverfy(4)
c     .. local scalars ..
      double precision  condbd, dctol, epsmch
      integer           idbg, j, k, lent, mjrdbg, mnrdbg, 
     *                  msgqp, nctotl, nmajor, nminor, nplin
      character*16      key
c     .. local arrays ..
      double precision  rprmls(mxparm), rprmnp(mxparm)
      integer           iprmls(mxparm), iprmnp(mxparm)
      character*3       chess(0:1)
      character*4       icrsh(0:2)
      character*80      rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

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
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (iprmnp(1),idbgnp), (rprmnp(1),cdint)
      equivalence       (idbgnp,idbg), (itmxnp,nmajor), (itmax2,nminor)
      equivalence       (ldbgls,mnrdbg), (ldbgnp,mjrdbg), (msgls,msgqp)
c     .. save statement ..
      save              /ee04uc/, /de04nc/, /ee04nc/,
     *                  /ge04uc/, /he04uc/
c     .. data statements ..
      data              icrsh(0), icrsh(1), icrsh(2)/'cold', 'warm',
     *                  'hot '/
      data              chess(0), chess(1)/' no', 'yes'/
c----------------------------------------------------------------------
      epsmch = wmach(3)
      nout = 6

      condbd = max(one/(hundrd*epsmch*dble(n)),tenp6)

      nplin = n + nclin
      nctotl = nplin + ncnln

c     make a dummy call e04ucq to ensure that the defaults are set.

      call e04ucq(nout,'*',key)
      newopt = .true.

c     save the optional parameters set by the user.  the values in
c     iprmls, rprmls, iprmnp and rprmnp may be changed to their
c     default values.
c
      call f06dff(mxparm,iprmls,1,ipsvls,1)
      call dcopy (mxparm,rprmls,1,rpsvls,1)
      call f06dff(mxparm,iprmnp,1,ipsvnp,1)
      call dcopy (mxparm,rprmnp,1,rpsvnp,1)
c
      if (msgnp.eq.idummy) msgnp = 10
      if (msgqp.eq.idummy) msgqp = 0
      if (iprnt.lt.0) iprnt = nout
      if (isumry.lt.0 .or. (msgnp.lt.5 .and. msgqp.lt.5)) isumry = -1
      iprint = iprnt
      isumm = isumry
      if (lcrash.lt.0 .or. lcrash.gt.2) lcrash = 0
      if (lvlder.lt.0 .or. lvlder.gt.3) lvlder = 3
      if (lformh.lt.0 .or. lformh.gt.1) lformh = 0
c
      if (nmajor.lt.0) nmajor = max(50,3*nplin+10*ncnln)
      if (nminor.lt.1) nminor = max(50,3*nctotl)

      nlnf = n
      nlnj = n
      nlnx = n
      if (jvrfy2.le.0 .or. jvrfy2.gt.n) jvrfy2 = n
      if (jvrfy1.le.0 .or. jvrfy1.gt.jvrfy2) jvrfy1 = 1
      if (jvrfy4.le.0 .or. jvrfy4.gt.n) jvrfy4 = n
      if (jvrfy3.le.0 .or. jvrfy3.gt.jvrfy4) jvrfy3 = 1
      if ((lverfy.lt.-1 .or. lverfy.gt.13)
     *    .or. (lverfy.ge.4 .and. lverfy.le.9)) lverfy = 0
c
      if (ksave.le.0) ksave = nmajor + 1
      if (nsave.lt.0) nsave = 0
      if (nsave.eq.0) ksave = nmajor + 1
      if (nload.lt.0) nload = 0
      if (lcrash.le.1) nload = 0
      if (nload.eq.0 .and. lcrash.eq.2) lcrash = 0
c
      if (tolact.lt.zero .or. tolact.ge.one) tolact = wrktol
      if (tolfea.lt.epsmch .or. tolfea.ge.one) tolfea = epspt5
      if (epsrf.lt.epsmch .or. epsrf.ge.one) epsrf = epspt9
      lfdset = 0
      if (fdint.lt.zero) lfdset = 2
      if (fdint.eq.rdummy) lfdset = 0
      if (fdint.ge.epsmch .and. fdint.lt.one) lfdset = 1
      if (lfdset.eq.1 .and. (cdint.lt.epsmch .or. cdint.ge.one))
     *    cdint = epsrf**point3
      if (bigbnd.le.zero) bigbnd = gigant
      if (bigdx.le.zero) bigdx = max(gigant,bigbnd)
      if (dxlim.le.zero) dxlim = two
      if (eta.lt.zero .or. eta.ge.one) eta = point9
      if (ftol.lt.epsrf .or. ftol.ge.one) ftol = epsrf**point8
c
      if (hcndbd.lt.one) hcndbd = condbd
c
      dctol = epspt5
      if (lvlder.lt.2) dctol = epspt3
      if (ctol.lt.epsmch .or. ctol.ge.one) ctol = dctol
c
      itmax1 = max(50,3*(n+nclin+ncnln))
      jverfy(1) = jvrfy1
      jverfy(2) = jvrfy2
      jverfy(3) = jvrfy3
      jverfy(4) = jvrfy4

      k = 1

      if (msgnp.gt.0) then

c        print the title. if no hot start is specified, the parameters
c        are final and can be printed.

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
            write (rec,fmt=99996) ncnln, ctol, nlnf, ftol, nlnj, eta
            call x04bay(iprint,4,rec)
            write (rec,fmt=99995) lvlder, epsrf, lverfy, isumry
            call x04bay(iprint,3,rec)
            if (lverfy.gt.0) then
               write (rec,fmt=99994) jvrfy1, jvrfy2
               call x04bay(iprint,2,rec)
               if (ncnln.gt.0) then
                  write (rec,fmt=99993) jvrfy3, jvrfy4
                  call x04baf(iprint,rec(1))
               end if
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

         end if
      end if

      if (msgnp.ge.5 .or. msgqp.ge.5) then
         if (isumm.ge.0 .and. isumm.ne.iprint) then
            lent = len(title)
            write (rec,fmt=99987) (title(j:j),j=1,lent)
            call x04bay(isumm,2,rec)
         end if
      end if

c     end of  e04ucx. (npdflt)

99999 format (/' parameters',/' -')
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

      subroutine e04ncu(cold,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *                  lda,istate,kactiv,bigbnd,tolact,a,ax,bl,bu,x,wx,
     *                  work)
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      implicit none

      integer ldbg
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
c     .. local scalars ..
      double precision  b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, toobig
      integer           i, imin, is, j, jmin, k, nfixed
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)
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

      call dcopy (n,x,1,wx,1)

      nfixed = 0
      nactiv = 0
      nartif = 0

c     if a cold start is being made, initialize  istate.
c     if  bl(j) = bu(j),  set  istate(j)=3  for all variables and linear
c     constraints.

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

c     initialize nfixed, nfree and kactiv.
c     ensure that the number of bounds and general constraints in the
c     working set does not exceed n.

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

c     if a cold start is required,  attempt to add as many
c     constraints as possible to the working set.

      if (cold) then
c        see if any bounds are violated or nearly satisfied.
c        if so,  add these bounds to the working set and set the
c        variables exactly on their bounds.

         j = n

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

         end if

c        the following loop finds the linear constraint (if any) with
c        smallest residual less than or equal to tolact  and adds it
c        to the working set.  this is repeated until the working set
c        is complete or all the remaining residuals are too large.

c        compute the residuals for all the constraints not in the
c        working set.

         if (nclin.gt.0 .and. nactiv+nfixed.lt.n) then
            do 140 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot(n,a(i,1),lda,wx,1)
  140       continue

            is = 1
            toobig = tolact + tolact

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
c        this is an expensive loop.  later replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).

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

         end if
      end if

      nfree = n - nfixed

c     end of  e04ncu. (lscrsh)

99999 format (/' //e04ncu// cold nclin nctotl',/' //e04ncu// ',l4,i6,i7)
99998 format (/' //e04ncu// variables before crash... ')
99997 format (/' //e04ncu// variables after  crash... ')
99996 format (/' //e04ncu// working set selected ...             ',
     *       /' //e04ncu// nfixed nactiv nartif      ',/' //e04ncu// ',
     *       i6,2i7)
99995 format (5g12.3)
      end

      subroutine e04nbw(mode,n,nz,nfree,nq,unitq,kx,v,zy,wrk)
c----------------------------------------------------------------------
c     e04nbw  transforms the vector  v  in various ways using the
c     matrix  q = ( z  y )  defined by the input parameters.
c
c        mode               result
c        -               
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
c----------------------------------------------------------------------
      implicit none

      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)

      integer           mode, n, nfree, nq, nz
      logical           unitq

      double precision  v(n), wrk(n), zy(nq,*)
      integer           kx(n)

      integer           j, j1, j2, k, l, lenv, nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree
      j1 = 1
      j2 = nfree
      if (mode.eq.1 .or. mode.eq.4) j2 = nz
      if (mode.eq.2 .or. mode.eq.5 .or. mode.eq.7) j1 = nz + 1
      lenv = j2 - j1 + 1
      if (mode.le.3) then

c        mode = 1, 2  or  3.

         if (nfree.gt.0) call sload (nfree,zero,wrk,1)

c        copy  v(fixed)  into the end of  wrk.

         if (mode.ge.2 .and. nfixed.gt.0) call dcopy (nfixed,v(nfree+1),
     *       1,wrk(nfree+1),1)
c
c        set  wrk  =  relevant part of  zy * v.
c
         if (lenv.gt.0) then
            if (unitq) then
               call dcopy (lenv,v(j1),1,wrk(j1),1)
            else
               call dgemv ('n',nfree,j2-j1+1,one,zy(1,j1),nq,v(j1),1,one
     *                    ,wrk,1)
            end if
         end if
c
c        expand  wrk  into  v  as a full n-vector.
c
         call sload (n,zero,v,1)
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
c
      else

c        mode = 4, 5, 6, 7  or  8.

c        put the fixed components of  v  into the end of  wrk.
c
         if (mode.eq.5 .or. mode.eq.6) then
            do 60 l = 1, nfixed
               j = kx(nfree+l)
               wrk(nfree+l) = v(j)
   60       continue
         end if
c
c        put the free  components of  v  into the beginning of  wrk.
c
         if (nfree.gt.0) then
            do 80 k = 1, nfree
               j = kx(k)
               wrk(k) = v(j)
   80       continue
c
c           set  v  =  relevant part of  zy' * wrk.
c
            if (lenv.gt.0) then
               if (unitq) then
                  call dcopy (lenv,wrk(j1),1,v(j1),1)
               else
                  call dgemv ('t',nfree,j2-j1+1,one,zy(1,j1),nq,wrk,1,
     *                       zero,v(j1),1)
               end if
            end if
         end if
c
c        copy the fixed components of  wrk  into the end of  v.
c
         if (nfixed.gt.0 .and. (mode.eq.5 .or. mode.eq.6))
     *       call dcopy (nfixed,wrk(nfree+1),1,v(nfree+1),1)
      end if

c     end of  e04nbw. (cmqmul)

      end

      subroutine e04ncz(prbtyp,named,names,linobj,unitq,inform,iter,
     *                  jinf,nclin,nctotl,nactiv,nfree,nrank,nz,nrz,n,
     *                  lda,ldr,istate,kactiv,kx,ctx,ssq,ssq1,suminf,
     *                  numinf,xnorm,bl,bu,a,clamda,ax,featol,r,x,w)
c----------------------------------------------------------------------
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
c     
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
c----------------------------------------------------------------------
      implicit none

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
c     .. scalar arguments ..
      double precision  ctx, ssq, ssq1, suminf, xnorm
      integer           inform, iter, jinf, lda, ldr, n, nactiv, nclin,
     *                  nctotl, nfree, nrank, nrz, numinf, nz
      logical           linobj, named, unitq
      character*2       prbtyp
c     .. array arguments ..
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
c     .. local scalars ..
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
c     .. local arrays ..
      double precision  rprmls(mxparm)
      integer           iprmls(mxparm)
      character*80      rec(3)
c     .. external functions ..
      double precision  dnrm2, adivb
      external          dnrm2, adivb
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ae04nc/locls

      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04nb/lennam, ldt, ncolt, ldzy
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04nc/ilsdbg, lsdbg
      common            /de04nb/asize, dtmax, dtmin
      common            /de04nc/ipsvls, idbgls, iprnt, isumry, itmax1,
     *                  itmax2, lcrash, ldbgls, lprob, msgls, nn,
     *                  nnclin, nprob, ipadls
      common            /ee04nc/rpsvls, bigbnd, bigdx, bndlow, bndupp,
     *                  tolact, tolfea, tolrnk, rpadls
      common            /fe04nb/icmdbg, cmdbg
c     .. equivalences ..
      equivalence       (iprmls(1),idbgls), (rprmls(1),bigbnd)
      equivalence       (msgls,msglvl), (idbgls,idbg), (ldbgls,msgdbg)
c     .. save statement ..
      save              /de04nc/, /ee04nc/
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

   20 if (needfg) then
         if (nrank.gt.0) then
            resnrm = dnrm2(nrank,w(lres),1)
            ssq = half*(ssq1**2+resnrm**2)
         end if

         if (numinf.ne.0) then

            call e04ncp(prbtyp,linobj,singlr,unitgz,unitq,n,nclin,nfree,
     *                  lda,ldzy,ldr,nrank,nz,nrz,istate,kx,bigbnd,
     *                  tolrnk,numinf,suminf,bl,bu,a,w(lres),featol,
     *                  w(lgq),w(lcq),r,x,w(lwtinf),w(lzy),w(lwrk))
            if (numinf.eq.0 .and. prbtyp.ne.'fp' .and. nphase.eq.1) then
               itmax = iter + itmax2
               nphase = 2
            end if
         end if
      end if

      gznorm = zero
      if (nz.gt.0) gznorm = dnrm2(nz,w(lgq),1)

      if (nrz.eq.nz) then
         grznrm = gznorm
      else
         grznrm = zero
         if (nrz.gt.0) grznrm = dnrm2(nrz,w(lgq),1)
      end if

      gfnorm = gznorm
      if (nfree.gt.0 .and. nactiv.gt.0) gfnorm = dnrm2(nfree,w(lgq),1)

      if (nrz.eq.0) then
         singlr = .false.
      else
         if (numinf.gt.0 .or. nrz.gt.nrank) then
            absrzz = zero
            singlr = .true.
         else
            call f06flf(nrz,r,ldr+1,drzmax,drzmin)
            absrzz = abs(r(nrz,nrz))
            rownrm = dnrm2(n,r(1,1),ldr)
            singlr = absrzz .le. drzmax*tolrnk .or. rownrm .le.
     *               tolrnk .or. abs(r(1,1)) .le. rownrm*tolrnk
         end if

      end if

      condrz = adivb(drzmax,drzmin,overfl)

      condt = one
      if (nactiv.gt.0) condt = adivb(dtmax,dtmin,overfl)

      if (prnt) then
         call e04ncj(prbtyp,isdel,iter,jadd,jdel,msglvl,nactiv,nfree,n,
     *               nclin,nrank,ldr,ldt,nz,nrz,istate,alfa,condrz,
     *               condt,gfnorm,grznrm,numinf,suminf,ctx,ssq,ax,r,
     *               w(lt),x,w(lwrk))

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

      statpt = .not. singlr .and. grznrm .le. dinky .or. irefn .gt.
     *         mrefn

      if ( .not. statpt) then

c        compute a search direction.

         prnt = .true.

         error = iter .ge. itmax
         if ( .not. error) then
c
            irefn = irefn + 1
            iter = iter + 1

            call e04ncq(linobj,singlr,unitgz,unitq,n,nclin,nfree,lda,
     *                  ldzy,ldr,nrank,numinf,nrz,kx,ctp,pnorm,a,w(lap),
     *                  w(lres),w(lhz),w(lpx),w(lgq),w(lcq),r,w(lzy),
     *                  w(lwrk))

            bigalf = adivb(bigdx,pnorm,overfl)
c
            call e04ucg(firstv,hitlow,istate,inform,jadd,n,nctotl,
     *                  numinf,alfhit,palfa,atphit,bigalf,bigbnd,pnorm,
     *                  w(lanorm),w(lap),ax,bl,bu,featol,w(lpx),x)

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
c
            error = unbndd .or. cyclin
            if ( .not. error) then

c              set x = x + alfa*p.  update ax, gq, res and ctx.

               if (alfa.ne.zero) call lsmove(hitcon,hitlow,linobj,
     *                                unitgz,nclin,nrank,nrz,n,ldr,jadd,
     *                                numinf,alfa,ctp,ctx,xnorm,w(lap),
     *                                ax,bl,bu,w(lgq),w(lhz),w(lpx),
     *                                w(lres),r,x,w(lwrk))
c
               if (hitcon) then

c                 add a constraint to the working set.
c                 update the tq factors of the working set.
c                 use p as temporary work space.

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

                  call e04ncv(unitq,inform,ifix,iadd,jadd,nactiv,nz,
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

                  irefn = 0
               end if

c              check the feasibility of constraints with non-
c              negative istate values.  if some violations have
c              occurred.  refine the current x and set inform so
c              that feasibility is checked in e04ncp.

               call e04ncr(n,nclin,istate,bigbnd,cnorm,err1,jmax1,nviol,
     *                     ax,bl,bu,featol,x,w(lwrk))

               if (err1.gt.featol(jmax1)) then
                  call e04nch(linobj,rowerr,unitq,nclin,nactiv,nfree,
     *                        nrank,nz,n,nctotl,ldzy,lda,ldr,ldt,istate,
     *                        kactiv,kx,jmax1,err2,ctx,xnorm,a,ax,bl,bu,
     *                        w(lcq),w(lres),w(lres0),featol,r,w(lt),x,
     *                        w(lzy),w(lpx),w(lwrk))

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

c        until      statpt  .or.  error
      if ( .not. (statpt .or. error)) go to 20


c     try and find the index jdel of a constraint to drop from
c     the working set.

      jdel = 0

      if (numinf.eq.0 .and. prbtyp.eq.'fp') then
         if (n.gt.nz) call sload (n-nz,(zero),w(lrlam),1)
         jtiny = 0
         jsmlst = 0
         jbigst = 0
      else

         call e04nck(prbtyp,msglvl,n,nactiv,nfree,lda,ldt,numinf,nz,nrz,
     *               istate,kactiv,kx,dinky,jsmlst,ksmlst,jinf,jtiny,
     *               jbigst,kbigst,trulam,a,w(lanorm),w(lgq),w(lrlam),
     *               w(lt),w(lwtinf))

      end if

      if ( .not. error) then
         if (jsmlst.gt.0) then

c           e04nck found a regular constraint with multiplier less
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

         if (jdel.ne.0 .and. singlr) then

c           cannot delete a constraint when rz is singular.
c           probably a weak minimum.

            jdel = 0
         else if (jdel.ne.0) then

c           constraint jdel has been deleted.
c           update the matrix factorizations.

            call e04nct(unitq,n,nactiv,nfree,nres,ngq,nz,nrz,lda,ldzy,
     *                  ldr,ldt,nrank,jdel,kdel,kactiv,kx,a,w(lres),r,
     *                  w(lt),w(lgq),w(lzy),w(lwrk),w(lpx))
         end if
      end if

      irefn = 0
      convrg = jdel .eq. 0

      prnt = .false.
      uncon = .false.
      needfg = .false.

      if ( .not. (convrg .or. error)) go to 20

      weak = jtiny .gt. 0 .or. singlr

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

c     set   clamda.  print the full solution.

c DEBUG 691                       random assignment (null)
      msglvl = 0

      if (msglvl.gt.0) then
         write (rec,fmt=99999) prbtyp, iter
         call x04bay(iprint,2,rec)
      end if

      call e04nbx(msglvl,nfree,lda,n,nclin,nctotl,bigbnd,named,names,
     *            nactiv,istate,kactiv,kx,a,bl,bu,x,clamda,w(lrlam),x)

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

      subroutine e04ucy(inform,msgnp,nstate,lvlder,nfun,ngrad,ldcj,
     *                  ldcju,n,ncnln,confun,objfun,needc,bigbnd,epsrf,
     *                  cdint,fdint,fdchk,fdnorm,objf,xnorm,bl,bu,c,c1,
     *                  cjac,cjacu,cjdx,dx,grad,gradu,hforwd,hcntrl,x,
     *                  wrk1,wrk2,w,lenw,iuser,user)
c----------------------------------------------------------------------
c     e04ucy  does the following...
c     (1)  computes the objective and constraint values objf and c.
c     (2)  evaluates the user-provided gradients in cjacu and gradu.
c     (3)  counts the missing gradients.
c     (4)  loads the known gradients into grad and cjac.
c     (5)  checks that the known gradients are programmed correctly.
c     (6)  computes the missing gradient elements.

      double precision  rdummy
      parameter         (rdummy=-11111.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, cdint, epsrf, fdchk, fdint, fdnorm,
     *                  objf, xnorm
      integer           inform, ldcj, ldcju, lenw, lvlder, msgnp, n,
     *                  ncnln, nfun, ngrad, nstate
c     .. array arguments ..
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
c     .. local scalars ..
      integer           i, infog, infoj, j, mode, ncset
      logical           centrl, needfd
c     .. local arrays ..
      character*80      rec(3)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /be04uc/lvldif, ncdiff, nfdiff, lfdset
      common            /ce04nb/epspt3, epspt5, epspt8, epspt9
      common            /ce04uc/lvrfyc, jverfy
c----------------------------------------------------------------------
      infog = 0
      infoj = 0
      nfdiff = 0
      ncdiff = 0
      ncset = n*ncnln
c
      if (ncnln.gt.0) then

c        compute the constraints and jacobian matrix.

c        if some derivatives are missing, load the jacobian with dummy
c        values.  any elements left unaltered after the call to confun
c        must be estimated.  a record of the missing jacobian elements
c        is stored in  cjacu.

         needfd = lvlder .eq. 0 .or. lvlder .eq. 1

         if (needfd) call f06qhf('g',ncnln,n,rdummy,rdummy,cjacu,
     *                           ldcju)

         call f06dbf(ncnln,(1),needc,1)

         mode = 2
         call confun(mode,ncnln,n,ldcju,needc,x,c,cjacu,nstate,iuser,
     *               user)
         if (mode.lt.0) go to 80

         call f06qff('g',ncnln,n,cjacu,ldcju,cjac,ldcj)

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
                  if (msgnp.gt.0) then
                     write (rec,fmt=99999) lvlder
                     call x04bay(iprint,3,rec)
                  end if
               else
                  if (msgnp.gt.0) then
                     write (rec,fmt=99998) ncset, n*ncnln, ncdiff
                     call x04bay(iprint,3,rec)
                  end if
               end if
            end if
         end if
      end if

c     repeat the procedure above for the objective function.

      needfd = lvlder .eq. 0 .or. lvlder .eq. 2

      if (needfd) call sload (n,rdummy,gradu,1)
c                                 output the initial value
      iuser(2) = 1

      mode = 2

      call objfun (mode,n,x,objf,gradu,nstate,iuser,user)

      if (mode.lt.0) go to 80

      call dcopy (n,gradu,1,grad,1)
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

      nfun = nfun + 1
      ngrad = ngrad + 1

c     check whatever gradient elements have been provided.

      if (lvrfyc.ge.0) then
         if (ncset.gt.0) then
            call e04xaw(mode,lvlder,msgnp,ncset,n,ncnln,ldcj,ldcju,
     *                  bigbnd,epsrf,epspt3,fdchk,xnorm,confun,needc,bl,
     *                  bu,c,c1,cjac,cjacu,cjdx,dx,wrk2,x,wrk1,iuser,
     *                  user)
            if (mode.lt.0) go to 80
            infoj = mode
         end if

         if (nfdiff.lt.n) then
            call e04xax(mode,msgnp,n,bigbnd,epsrf,epspt3,fdchk,objf,
     *                  xnorm,objfun,bl,bu,grad,gradu,dx,x,wrk1,iuser,
     *                  user)
            if (mode.lt.0) go to 80
            infog = mode
         end if
      end if

      needfd = ncdiff .gt. 0 .or. nfdiff .gt. 0

      if (needfd) then

c        compute the missing gradient elements.

         call e04xay(mode,msgnp,lvlder,n,ncnln,ldcj,ldcju,bigbnd,epsrf,
     *               fdnorm,objf,confun,objfun,needc,bl,bu,c,c1,cjdx,
     *               cjac,cjacu,grad,gradu,hforwd,hcntrl,x,dx,iuser,
     *               user)

         if (mode.lt.0) go to 80

         if (lfdset.gt.0) then
            centrl = lvldif .eq. 2
            call e04uds(centrl,mode,ldcj,ldcju,n,ncnln,bigbnd,cdint,
     *                  fdint,fdnorm,objf,confun,objfun,needc,bl,bu,c,
     *                  c1,cjdx,cjac,cjacu,grad,gradu,hforwd,hcntrl,x,w,
     *                  lenw,iuser,user)
c
            if (mode.lt.0) go to 80
         end if
      end if
c
      inform = infoj + infog
      return

c     the user requested termination.

   80 inform = mode

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

      subroutine e04mfk(msglvl,nfree,lda,n,nclin,ncnln,nctotl,bigbnd,
     *                  named,names,lennam,nactiv,istate,kactiv,kx,a,bl,
     *                  bu,c,clamda,rlamda,x)
c----------------------------------------------------------------------
c     e04mfk   creates the expanded lagrange multiplier vector clamda.
c     if  msglvl .eq. 1  or  msglvl .ge. 10,  e04mfk  prints  x,  a*x,
c     c(x),  their bounds, the multipliers, and the residuals (distance
c     to the nearer bound).

      integer           lcmdbg
      parameter         (lcmdbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd
      integer           lda, lennam, msglvl, n, nactiv, nclin, ncnln,
     *                  nctotl, nfree
      logical           named
c     .. array arguments ..
      double precision  a(lda,*), bl(nctotl), bu(nctotl), c(*),
     *                  clamda(nctotl), rlamda(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
      character*8       names(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(lcmdbg)
c     .. local scalars ..
      double precision  b1, b2, res, res2, rlam, v, wlam
      integer           ip, is, j, k, nfixed, nplin, nz
      character*2       ls
      character*5       id3
      character*8       id4
c     .. local arrays ..
      character*2       lstate(7)
      character*5       id(3)
      character*80      rec(4)
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04nb/icmdbg, cmdbg
c     .. data statements ..
      data              id(1)/'varbl'/
      data              id(2)/'lncon'/
      data              id(3)/'nlcon'/
      data              lstate(1)/'--'/, lstate(2)/'++'/
      data              lstate(3)/'fr'/, lstate(4)/'ll'/
      data              lstate(5)/'ul'/, lstate(6)/'eq'/
      data              lstate(7)/'tf'/
c----------------------------------------------------------------------
c
      nplin = n + nclin
      nz = nfree - nactiv
c
c     expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.
c
      call sload (nctotl,zero,clamda,1)
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
c
      if (msglvl.ne.1 .and. msglvl.lt.10) return
c
      write (rec,fmt=99999)
      call x04bay(iprint,4,rec)
      id3 = id(1)
c
      do 40 j = 1, nctotl
         b1 = bl(j)
         b2 = bu(j)
         wlam = clamda(j)
         is = istate(j)
         ls = lstate(is+3)
         if (j.le.n) then
c
c           section 1 -- the variables  x.
c           
            k = j
            v = x(j)
c
         else if (j.le.nplin) then
c
c           section 2 -- the linear constraints  a*x.

            if (j.eq.n+1) then
               write (rec,fmt=99998)
               call x04bay(iprint,4,rec)
               id3 = id(2)
            end if
c
            k = j - n
            v = ddot(n,a(k,1),lda,x,1)
         else
c
c           section 3 -- the nonlinear constraints  c(x).
c           
c
            if (j.eq.nplin+1) then
               write (rec,fmt=99997)
               call x04bay(iprint,4,rec)
               id3 = id(3)
            end if
c
            k = j - nplin
            v = c(k)
         end if
c
c        print a line for the j-th variable or constraint.
c        -
         res = v - b1
         res2 = b2 - v
         if (abs(res).gt.abs(res2)) res = res2
         ip = 1
         if (b1.le.(-bigbnd)) ip = 2
         if (b2.ge.bigbnd) ip = ip + 2
         if (named) then
c
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
c
         else
c
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

c     end of  e04mfk.  (cmprnt)

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

      subroutine e04ncx(unitq,inform,nz,nfree,nrank,nres,ngq,n,ldzy,lda,
     *                  ldr,ldt,istate,kx,condmx,a,r,t,res,gq,zy,w,c,s,
     *                  msglvl)
c----------------------------------------------------------------------
c     e04ncx updates the factor r as kx is reordered to reflect the
c     status of the bound constraints given by istate.  kx is reordered
c     so that the fixed variables come last.  one of two alternative
c     are used to reorder kx. one method needs fewer accesses to kx, the
c     other gives a matrix rz with more rows and columns.

      double precision  condmx
      integer           inform, lda, ldr, ldt, ldzy, msglvl, n, nfree,
     *                  ngq, nrank, nres, nz
      logical           unitq
c     .. array arguments ..
      double precision  a(lda,*), c(n), gq(n,*), r(ldr,*), res(n,*),
     *                  s(n), t(ldt,*), w(n), zy(ldzy,*)
      integer           istate(*), kx(n)
c     .. local scalars ..
      integer           iadd, ifix, j, j2, jadd, k, l, lstart, nactv,
     *                  nfixed
c----------------------------------------------------------------------
      nfixed = n - nfree
c
      if (nrank.lt.n .and. nrank.gt.0) then

c        r is specified but singular.  try and keep the dimension of rz
c        as large as possible.

         nactv = 0
         nfree = n
         nz = n

         j = n

   20    if (j.gt.0 .and. n-nfree.lt.nfixed) then
            if (istate(j).gt.0) then
               jadd = j
               do 40 ifix = nfree, 1, -1
                  if (kx(ifix).eq.jadd) go to 60
   40          continue

c              add bound jadd.

   60          call e04ncv(unitq,inform,ifix,iadd,jadd,nactv,nz,nfree,
     *                     nrank,nres,ngq,n,lda,ldzy,ldr,ldt,kx,condmx,
     *                     a,r,t,res,gq,zy,w,c,s,msglvl)

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

c           order kx so that the free variables come first.

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

                  if (nrank.gt.0) call e04nbu(n,nres,nrank,ldr,k,l,r,
     *                                 res,c,s)
               end if
  120       continue

         end if
         nz = nfree
      end if

c     end of  e04ncx. (lsbnds)

      end

      subroutine e04mfr(job,msglvl,n,nclin,nmoved,iter,numinf,istate,
     *                  bigbnd,ax,bl,bu,featol,featlu,x)
c----------------------------------------------------------------------
c     e04mfr does most of the manoeuvres associated with degeneracy.
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
c     if job = 'i', e04mfr initializes the parameters in
c     common block ce04mf:
c
c     tolx0   is the minimum (scaled) feasibility tolerance.
c     tolx1   is the maximum (scaled) feasibility tolerance.
c     tolinc  is the scaled increment to the current featol.
c     kdegen  is the expand frequency (specified by the user).
c             it is the frequency of resetting featol to (scaled) tolx0.
c     ndegen  counts the number of degenerate steps (incremented
c             by e04mfs).
c     itnfix  is the last iteration at which a job = 'e' or 'o' entry
c             caused an x to be put on a constraint.
c     nfix(j) counts the number of times a job = 'o' entry has
c             caused the variables to be placed on the working set,
c             where j=1 if infeasible, j=2 if feasible.
c
c     tolx0*featlu and tolx1*featlu are both close to the feasibility
c     tolerance featlu specified by the user.  (they must both be less
c     than featlu.)

c     if job = 'e',  e04mfr has been called after a cycle of kdegen
c     iterations.  constraints in the working set are examined to see if
c     any are off their bounds by an amount approaching featol.  nmoved
c     returns how many.  if nmoved is positive,  x  is moved onto the
c     constraints in the working set.  it is assumed that the calling
c     routine will then continue iterations.

c     if job = 'o',  e04mfr is being called after a subproblem has been
c     judged optimal, infeasible or unbounded.  constraint violations
c     are examined as above.

      double precision  zero, point6
      parameter         (zero=0.0d+0,point6=0.6d+0)
c     .. scalar arguments ..
      double precision  bigbnd
      integer           iter, msglvl, n, nclin, nmoved, numinf
      character         job
c     .. array arguments ..
      double precision  ax(*), bl(n+nclin), bu(n+nclin),
     *                  featlu(n+nclin), featol(n+nclin), x(n)
      integer           istate(n+nclin)
c     .. scalars in common ..
      double precision  tolinc, tolx0
      integer           iprint, isumm, itnfix, kdegen, lines1, lines2,
     *                  ndegen, nout
c     .. arrays in common ..
      integer           nfix(2)
c     .. local scalars ..
      double precision  d, epsmch, tolx1, tolz
      integer           is, j, maxfix
      character*80      rec
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /ce04mf/tolx0, tolinc, kdegen, ndegen, itnfix,
     *                  nfix
c     .. save statement ..
      save              tolz
c----------------------------------------------------------------------
      nmoved = 0
      if (job.eq.'i' .or. job.eq.'i') then

c        job = 'initialize'.
c        initialize at the start of each linear problem.
c        kdegen  is the expand frequency      and
c        featlu  are the user-supplied feasibility tolerances.
c        they are not changed.

         epsmch = wmach(3)

         ndegen = 0
         itnfix = 0
         nfix(1) = 0
         nfix(2) = 0
         tolx0 = 0.5d+0
         tolx1 = 0.99d+0
         tolz = epsmch**point6

         if (kdegen.lt.9999999) then
            tolinc = (tolx1-tolx0)/kdegen
         else
            tolinc = zero
         end if

         do 20 j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
   20    continue
      else

c        job = 'end of cycle' or 'optimal'.
c        initialize local variables maxfix and tolz.

         maxfix = 2

         if (job.eq.'o' .or. job.eq.'o') then

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

c        reset featol to its minimum value.

         do 40 j = 1, n + nclin
            featol(j) = tolx0*featlu(j)
   40    continue

c        count the number of times a variable is moved a nontrivial
c        distance onto its bound.

         itnfix = iter

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
c
         if (nmoved.gt.0) then
c
c           some variables were moved onto their bounds.
c
            if (msglvl.gt.0) then
               write (rec,fmt=99999) iter, nmoved
               call x04baf(iprint,rec)
            end if
         end if
      end if

c     end of e04mfr.  (cmdgen)

99999 format (' itn',i6,' --',i7,'  variables moved to their bounds.')
      end

      subroutine e04nbz(nerror,msglvl,lcrash,liwork,lwork,litotl,lwtotl,
     *                  n,nclin,ncnln,istate,kx,named,names,bigbnd,bl,
     *                  bu,x,m,lda,ldr,ldcj,ldfj,nerr,ifail)
c----------------------------------------------------------------------
c     e04nbz   checks the input data for nlpsol
c----------------------------------------------------------------------
      double precision  bigbnd
      integer           ifail, lcrash, lda, ldcj, ldfj, ldr, litotl,
     *                  liwork, lwork, lwtotl, m, msglvl, n, nclin,
     *                  ncnln, nerr, nerror
      logical           named
c     .. array arguments ..
      double precision  bl(n+nclin+ncnln), bu(n+nclin+ncnln), x(n)
      integer           istate(n+nclin+ncnln), kx(n)
      character*8       names(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
c     .. local scalars ..
      double precision  b1, b2
      integer           is, j, k, l
      logical           ok
c     .. local arrays ..
      character*5       id(3)
      character*80      rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
c     .. data statements ..
      data              id(1), id(2), id(3)/'varbl', 'l con', 'n con'/
c----------------------------------------------------------------------

      nerror = 0

c     check  m.

      if (m.le.0) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99997) m
            call x04bay(nerr,3,rec)
         end if
      end if

c     check  n.

      if (n.le.0) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99996) n
            call x04bay(nerr,3,rec)
         end if
      end if

c     check  nclin and ncnln.

      if (nclin.lt.0 .or. ncnln.lt.0) then
         if (nclin.lt.0) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99995) nclin
               call x04bay(nerr,3,rec)
            end if
         end if
c
         if (ncnln.lt.0) then
            nerror = nerror + 1
            if (ifail.eq.0 .or. ifail.eq.-1) then
               write (rec,fmt=99994) ncnln
               call x04bay(nerr,3,rec)
            end if
         end if
      end if

c     check  lda.

      if (lda.lt.max(1,nclin)) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99993) lda, nclin
            call x04bay(nerr,3,rec)
         end if
      end if

c     check  ldcj.

      if (ldcj.lt.max(1,ncnln)) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99992) ldcj, ncnln
            call x04bay(nerr,3,rec)
         end if
      end if

c     check  ldfj.

      if (ldfj.lt.m) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99991) ldfj, m
            call x04bay(nerr,3,rec)
         end if
      end if

c     check  ldr.

      if (ldr.lt.n) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
            write (rec,fmt=99990) ldr, n
            call x04bay(nerr,3,rec)
         end if
      end if

c     check if there is enough workspace to solve the problem.

      ok = litotl .le. liwork .and. lwtotl .le. lwork
      if ( .not. ok) then
         nerror = nerror + 1
         if (ifail.eq.0 .or. ifail.eq.-1) then
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

         do 20 j = 1, n + nclin + ncnln
            b1 = bl(j)
            b2 = bu(j)
            ok = b1 .lt. b2 .or. (b1.eq.b2 .and. abs(b1).lt.bigbnd)
            if ( .not. ok) then
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
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99989) names(j), j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99988) names(j), j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               else
                  if (b1.eq.b2) then
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99987) id(l), k, j, j, b1, bigbnd
                        call x04bay(nerr,4,rec)
                     end if
                  else
                     if (ifail.eq.0 .or. ifail.eq.-1) then
                        write (rec,fmt=99986) id(l), k, j, b1, j, b2
                        call x04bay(nerr,3,rec)
                     end if
                  end if
               end if
            end if
   20    continue

c        if warm start, check  istate.

         if (lcrash.eq.1) then
            do 40 j = 1, n + nclin + ncnln
               is = istate(j)
               ok = is .ge. (-2) .and. is .le. 4
               if ( .not. ok) then
                  nerror = nerror + 1
                  if (ifail.eq.0 .or. ifail.eq.-1) then
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

      subroutine e04ucs(cold,n,nclin,ncnln,nctotl,nactiv,nfree,nz,
     *                  istate,kactiv,bigbnd,tolact,bl,bu,c)
c----------------------------------------------------------------------
c     e04ucs  adds indices of nonlinear constraints to the initial
c     working set.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  one
      parameter         (one=1.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, tolact
      integer           n, nactiv, nclin, ncnln, nctotl, nfree, nz
      logical           cold
c     .. array arguments ..
      double precision  bl(nctotl), bu(nctotl), c(*)
      integer           istate(nctotl), kactiv(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  b1, b2, biglow, bigupp, cmin, res, resl, resu,
     *                  toobig
      integer           i, imin, is, j, linact, nfixed, nlnact, nplin
c     .. local arrays ..
      character*80      rec(4)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /fe04uc/inpdbg, npdbg

      double precision wmach
      common/ cstmch /wmach(10)
c----------------------------------------------------------------------
      nfixed = n - nfree
      linact = nactiv
      nplin = n + nclin

c     if a cold start is being made, initialize the status of the qp
c     working set.  first,  if  bl(j) = bu(j),  set istate(j)=3.

      if (cold) then
         do 20 j = nplin + 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
   20    continue
      end if

c     increment nactiv and kactiv.
c     ensure that the number of bounds and general constraints in the
c     qp  working set does not exceed n.

      do 40 j = nplin + 1, nctotl
         if (nfixed+nactiv.eq.n) istate(j) = 0
         if (istate(j).gt.0) then
            nactiv = nactiv + 1
            kactiv(nactiv) = j - n
         end if
   40 continue

      if (cold) then

c        if a cold start is required, an attempt is made to add as many
c        nonlinear constraints as possible to the working set.

c        the following loop finds the most violated constraint.  if
c        there is room in kactiv, it will be added to the working set
c        and the process will be repeated.

         is = 1
         biglow = -bigbnd
         bigupp = bigbnd
         toobig = tolact + tolact

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
      if (npdbg .and. inpdbg(1).gt.0) then
         write (rec,fmt=99999) nfixed, linact, nlnact
         call x04bay(iprint,4,rec)
      end if

c     end of  e04ucs. (npcrsh)

99999 format (/' //e04ucs//  working set selected....',/' //e04ucs// n',
     *       'fixed linact nlnact     ',/' //e04ucs//',3i7)
      end

      subroutine e04mfu(prbtyp,header,rset,msglvl,iter,isdel,jdel,jadd,
     *                  n,nclin,nactiv,nfree,nz,nrz,ldr,ldt,istate,alfa,
     *                  condrz,condt,drzz,gfnorm,grznrm,numinf,suminf,
     *                  notopt,objlp,trusml,ax,r,t,x,work)
c----------------------------------------------------------------------
c     e04mfu  prints various levels of output for e04mfz.

c           msg        cumulative result

c       .le.  0        no output.
c       .eq.  1        nothing now (but full output later).
c       .eq.  5        one terse line of output.
c       .ge. 10        same as 5 (but full output later).
c       .ge. 20        constraint status,  x  and  ax.
c       .ge. 30        diagonals of  t  and  r.

c     debug printing is controlled by the logical variable  lcdbg.
c     lcdbg  is set true after  idbg  iterations.  then, the amount of
c     output is determined by a string of binary digits of the form
c     svt  (stored in the integer array  ilcdbg).

c     s  set 'on'  gives information from the max step routine e04mfs.
c     v  set 'on'  gives various vectors in e04mfz  and its auxiliaries.
c     t  set 'on'  gives a trace of which routine was called and an
c                  indication of the progress of the run.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero
      parameter         (zero=0.0d+0)
      integer           mline1, mline2
      parameter         (mline1=50000,mline2=50000)
c     .. scalar arguments ..
      double precision  alfa, condrz, condt, drzz, gfnorm, grznrm,
     *                  objlp, suminf, trusml
      integer           isdel, iter, jadd, jdel, ldr, ldt, msglvl, n,
     *                  nactiv, nclin, nfree, notopt, nrz, numinf, nz
      logical           header, rset
      character*2       prbtyp
c     .. array arguments ..
      double precision  ax(*), r(ldr,*), t(ldt,*), work(n), x(n)
      integer           istate(*)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           lcdbg
c     .. arrays in common ..
      integer           ilcdbg(ldbg)
c     .. local scalars ..
      double precision  obj
      integer           itn, j, jj, k, kadd, kdel, kk, ndf
      logical           newset, prthdr
      character*2       ladd, ldel
      character*15      lmchar
c     .. local arrays ..
      character*2       lstate(0:5)
      character*120     rec(5)
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2
      common            /ee04mf/ilcdbg, lcdbg
c     .. data statements ..
      data              lstate(0), lstate(1), lstate(2)/'  ', 'l ',
     *                  'u '/
      data              lstate(3), lstate(4), lstate(5)/'e ', 'f ',
     *                  'a '/
c----------------------------------------------------------------------
      if (msglvl.ge.15) then
         if (isumm.ge.0) then
            write (rec,fmt=99999) prbtyp, iter
            call x04bay(isumm,5,rec)
         end if
      end if

      if (msglvl.ge.5) then

c        some printing required.  set up information for the terse line.

         itn = mod(iter,10000)
         ndf = mod(nrz,10000)

         if (jdel.ne.0) then
            if (notopt.gt.0) then
               write (lmchar,fmt='( i5, 1p,d10.2 )') notopt, trusml
            else
               write (lmchar,fmt='( 5x, 1p,d10.2 )') trusml
            end if

            if (jdel.gt.0) then
               kdel = isdel

            else if (jdel.lt.0) then
               jdel = nz - nrz + 1
               kdel = 5
            end if
         else
            jdel = 0
            kdel = 0
            lmchar = '               '
         end if

         if (jadd.gt.0) then
            kadd = istate(jadd)
         else
            kadd = 0
         end if

         ldel = lstate(kdel)
         ladd = lstate(kadd)

         if (numinf.gt.0) then
            obj = suminf
         else
            obj = objlp
         end if

c        if necessary, print a header.
c        print a single line of information.

         if (isumm.ge.0) then

c           terse line for the monitoring file.

            newset = lines1 .ge. mline1
            prthdr = msglvl .ge. 15 .or. header .or. newset

            if (prthdr) then
               write (rec,fmt=99997)
               call x04bay(isumm,3,rec)
               lines1 = 0
            end if

            write (rec,fmt=99995) itn, jdel, ldel, jadd, ladd, alfa,
     *        numinf, obj, n - nfree, nactiv, nz - nrz, ndf, grznrm,
     *        lmchar, condt
            call x04baf(isumm,rec(1))
            lines1 = lines1 + 1
         end if

         if (iprint.ge.0 .and. isumm.ne.iprint) then

c           terse line for the print file.

            newset = lines2 .ge. mline2
            prthdr = header .or. newset

            if (prthdr) then
               write (rec,fmt=99998)
               call x04bay(iprint,3,rec)
               lines2 = 0
            end if

            write (rec,fmt=99996) itn, alfa, numinf, obj, grznrm
            call x04baf(iprint,rec(1))
            lines2 = lines2 + 1
         end if

         if (msglvl.ge.20) then
            if (isumm.ge.0) then
               write (rec,fmt=99994) prbtyp
               call x04bay(isumm,3,rec)
               write (rec,fmt=99993)
               call x04bay(isumm,2,rec)
               do 20 j = 1, n, 5
                  write (rec,fmt=99992) (x(jj),istate(jj),jj=j,
     *              min(j+4,n))
                  call x04baf(isumm,rec(1))
   20          continue
               if (nclin.gt.0) then
                  write (rec,fmt=99991)
                  call x04bay(isumm,2,rec)
                  do 40 k = 1, nclin, 5
                     write (rec,fmt=99992) (ax(kk),istate(n+kk),kk=k,
     *                 min(k+4,nclin))
                     call x04baf(isumm,rec(1))
   40             continue
               end if

               if (msglvl.ge.30) then

c                 print the diagonals of  t  and  r.

                  if (nactiv.gt.0) then
                     call dcopy (nactiv,t(1,nz+1),ldt+1,work,1)
                     write (rec,fmt=99990) prbtyp
                     call x04bay(isumm,2,rec)
                     do 60 j = 1, nactiv, 5
                        write (rec,fmt=99989) (work(jj),jj=j,
     *                    min(j+4,nactiv))
                        call x04baf(isumm,rec(1))
   60                continue
                  end if
                  if (rset .and. nrz.gt.0) then
                     write (rec,fmt=99988) prbtyp
                     call x04bay(isumm,2,rec)
                     do 80 j = 1, nrz, 5
                        write (rec,fmt=99989) (r(jj,jj),jj=j,
     *                    min(j+4,nrz))
                        call x04baf(isumm,rec(1))
   80                continue
                  end if
               end if
               write (rec,fmt=99987)
               call x04bay(isumm,4,rec)
            end if
         end if
      end if

      header = .false.
      jdel = 0
      jadd = 0
      alfa = zero

c     end of e04mfu.  (lpprnt)

99999 format (///' ',a2,' iteration',i5,/' =================')
99998 format (//'  itn     step ninf sinf/objective  norm gz')
99997 format (//'  itn jdel  jadd      step ninf  sinf/objective  bnd ',
     *       ' lin  art   zr  norm gz nopt    min lm   cond t')
99996 format (i5,1p,d9.1,i5,d15.6,d9.1)
99995 format (i5,i5,a1,i5,a1,1p,d9.1,i5,d16.8,4i5,d9.1,a15,d9.1)
99994 format (/' values and status of the ',a2,' constraints',/' --',
     *       '-')
99993 format (/' variables...')
99992 format (1x,5(1p,d15.6,i5))
99991 format (/' general linear constraints...')
99990 format (/' diagonals of ',a2,' working set factor t')
99989 format (1p,5d15.6)
99988 format (/' diagonals of ',a2,' triangle rz        ')
99987 format (///' ','--')
      end

      subroutine e04udr(unitq,n,nfree,nz,nq,nrowr,iperm,kx,gq,r,zy,work,
     *                  qrwork)
c----------------------------------------------------------------------
c     e04udr  bounds the condition estimator of the transformed hessian.
c     on exit, r is of the form
c                  ( drz   0     )
c                  (  0  sigma*i )
c     where d is a diagonal matrix such that drz has a bounded condition
c     number,  i is the identity matrix and sigma  is the geometric mean
c     of the largest and smallest elements of drz. the qr factorization
c     with interchanges is used to give diagonals of drz that are
c     decreasing in modulus.
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, half, one
      parameter         (zero=0.0d+0,half=0.5d+0,one=1.0d+0)
c     .. scalar arguments ..
      integer           n, nfree, nq, nrowr, nz
      logical           unitq
c     .. array arguments ..
      double precision  gq(n), qrwork(2*n), r(nrowr,*), work(n),
     *                  zy(nq,*)
      integer           iperm(n), kx(n)
c     .. scalars in common ..
      double precision  drmax, drmin, rcndbd, rfrobn
      logical           npdbg
c     .. arrays in common ..
      integer           inpdbg(ldbg)
c     .. local scalars ..
      double precision  drgm, drgs, gjmax, scle, sumsq
      integer           info, j, jmax, jsave, nrank
c     .. external functions ..
      double precision  f06bmf
      integer           f06klf
      external          f06bmf, f06klf
c     .. common blocks ..
      common            /ee04nb/rcndbd, rfrobn, drmax, drmin
      common            /fe04uc/inpdbg, npdbg
c----------------------------------------------------------------------
c     bound the condition estimator of q'hq.

      if (nz.gt.1) then

c        refactorize rz.  interchanges are used to give diagonals
c        of decreasing magnitude.

         do 20 j = 1, nz - 1
            call sload (nz-j,zero,r(j+1,j),1)
   20    continue

         call f01qff('column iterchanges',nz,nz,r,nrowr,work,iperm,
     *               qrwork,info)

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

               gjmax = gq(jmax)
               gq(jmax) = gq(j)
               gq(j) = gjmax
            end if
   40    continue
      end if

      drgm = one

      if (nz.gt.0) then
         nrank = f06klf(nz,r,nrowr+1,one/rcndbd)
         drgm = half*sqrt(abs(r(1,1)*r(nrank,nrank)))
         drgs = abs(r(1,1))/rcndbd

         if (nz.gt.nrank) then
            do 60 j = nrank + 1, nz
               call sload (j-1,zero,r(1,j),1)
   60       continue
            call sload (nz-nrank,drgs,r(nrank+1,nrank+1),nrowr+1)
         end if
      end if

c     reset the range-space partition of the hessian.

      if (nz.lt.n) then
         do 80 j = nz + 1, n
            call sload (j,zero,r(1,j),1)
   80    continue
         call sload (n-nz,drgm,r(nz+1,nz+1),nrowr+1)
      end if

c     recompute the frobenius norm of r.

      scle = sqrt(dble(n-nz))*drgm
      sumsq = one
      do 100 j = 1, nz
         call f06fjf(j,r(1,j),1,scle,sumsq)
  100 continue
      rfrobn = f06bmf(scle,sumsq)

c     end of  e04udr. (nprset)

      end

      subroutine e04mft(start,vertex,nclin,nctotl,nactiv,nartif,nfree,n,
     *                  lda,istate,kactiv,kx,bigbnd,tolact,a,ax,bl,bu,
     *                  featol,x,wx,work)
c----------------------------------------------------------------------
c     e04mft  computes the quantities  istate (optionally),  kactiv,
c     nactiv,  nz  and  nfree  associated with the working set at x.

c     the computation depends upon the value of the input parameter
c     start,  as follows...

c     start = 'cold'  an initial working set will be selected. first,
c                     nearly-satisfied or violated bounds are added.
c                     next,  general linear constraints are added that
c                     have small residuals.

c     start = 'warm'  the quantities kactiv, nactiv and nfree are
c                     initialized from istate,  specified by the user.

c     if vertex is true, an artificial vertex is defined by fixing some
c     variables on their bounds.  infeasible variables selected for the
c     artificial vertex are fixed at their nearest bound.  otherwise,
c     the variables are unchanged.

c     values of istate(j)
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c----------------------------------------------------------------------
      implicit none

      integer ldbg
      parameter         (ldbg=5)
      double precision  zero, one
      parameter         (zero=0.0d+0,one=1.0d+0)
c     .. scalar arguments ..
      double precision  bigbnd, tolact
      integer           lda, n, nactiv, nartif, nclin, nctotl, nfree
      logical           vertex
      character*4       start
c     .. array arguments ..
      double precision  a(lda,*), ax(*), bl(nctotl), bu(nctotl),
     *                  featol(nctotl), work(n), wx(n), x(n)
      integer           istate(nctotl), kactiv(n), kx(n)
c     .. scalars in common ..
      integer           iprint, isumm, lines1, lines2, nout
      logical           cmdbg
c     .. arrays in common ..
      integer           icmdbg(ldbg)
c     .. local scalars ..
      double precision  b1, b2, biglow, bigupp, colmin, colsiz, flmax,
     *                  residl, resl, resmin, resu, tol, toobig
      integer           i, imin, is, j, jfix, jfree, jmin, k
c     .. external functions ..
      double precision  ddot
      external          ddot
c     .. common blocks ..
      common            /ae04nb/nout, iprint, isumm, lines1, lines2

      double precision wmach
      common/ cstmch /wmach(10)

      common            /fe04nb/icmdbg, cmdbg
c----------------------------------------------------------------------
      flmax = wmach(7)
      biglow = -bigbnd
      bigupp = bigbnd

c     move the variables inside their bounds.

      do 20 j = 1, n
         b1 = bl(j)
         b2 = bu(j)
         tol = featol(j)
c
         if (b1.gt.biglow) then
            if (x(j).lt.b1-tol) x(j) = b1
         end if
c
         if (b2.lt.bigupp) then
            if (x(j).gt.b2+tol) x(j) = b2
         end if
   20 continue

      call dcopy (n,x,1,wx,1)

      nfree = n
      nactiv = 0
      nartif = 0

      if (start.eq.'cold') then
         do 40 j = 1, nctotl
            istate(j) = 0
            if (bl(j).eq.bu(j)) istate(j) = 3
   40    continue

      else if (start.eq.'warm') then
         do 60 j = 1, nctotl
            if (istate(j).gt.3 .or. istate(j).lt.0) istate(j) = 0
            if (bl(j).ne.bu(j) .and. istate(j).eq.3) istate(j) = 0
   60    continue
      end if

c     define nfree and kactiv.
c     ensure that the number of bounds and general constraints in the
c     working set does not exceed n.

      do 80 j = 1, nctotl
         if (nactiv.eq.nfree) istate(j) = 0

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

c     if a cold start is required,  attempt to add as many
c     constraints as possible to the working set.

      if (start.eq.'cold') then
c
c        see if any bounds are violated or nearly satisfied.
c        if so,  add these bounds to the working set and set the
c        variables exactly on their bounds.

         j = n

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

         end if

c        the following loop finds the linear constraint (if any) with
c        smallest residual less than or equal to tolact  and adds it
c        to the working set.  this is repeated until the working set
c        is complete or all the remaining residuals are too large.

c        first, compute the residuals for all the constraints not in the
c        working set.

         if (nclin.gt.0 .and. nactiv.lt.nfree) then
            do 120 i = 1, nclin
               if (istate(n+i).le.0) ax(i) = ddot(n,a(i,1),lda,wx,1)
  120       continue

            is = 1
            toobig = tolact + tolact

  140       if (is.gt.0 .and. nactiv.lt.nfree) then
               is = 0
               resmin = tolact

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

               if (is.gt.0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j = n + imin
                  istate(j) = is
               end if
               go to 140

            end if
         end if
      end if

      if (vertex .and. nactiv.lt.nfree) then

c        find an initial vertex by temporarily fixing some variables.

c        compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).

         do 200 j = 1, n
            if (istate(j).eq.0) then
               colsiz = zero
               do 180 k = 1, nclin
                  if (istate(n+k).gt.0) colsiz = colsiz + abs(a(k,j))
  180          continue
               work(j) = colsiz
            end if
  200    continue

c        find the  nartif  smallest such columns.
c        this is an expensive loop.  later replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).

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

c           fix x(j) at its current value.

  260       istate(j) = 4
            nartif = nartif + 1
            nfree = nfree - 1
            go to 220

         end if
      end if

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

c     end of  e04mft.  (cmcrsh)

99999 format (/' //e04mft// start nclin nctotl',/' //e04mft// ',a4,i6,
     *       i7)
99998 format (/' //e04mft// working set selected ...  ',/' //e04mft// ',
     *       'nfixed nactiv nartif      ',/' //e04mft// ',i6,2i7)
      end

      subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )
c----------------------------------------------------------------------
c BLAS Level 2 subroutine to the matrix-vector operations

c     x := a*x,   or   x := a'*x,

c  where x is n element vector and a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.
c----------------------------------------------------------------------
      implicit none

      integer            incx, lda, n
      character*1        diag, trans, uplo

      double precision   a( lda, * ), x( * )

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
c----------------------------------------------------------------------
c     test the input parameters.

      info = 0
      if     ( .not.(uplo .eq.'u' ).and.
     *         .not.(uplo .eq.'l' )      )then
         info = 1
      else if( .not.(trans.eq.'n' ).and.
     *         .not.(trans.eq.'t' ).and.
     *         .not.(trans.eq.'c' )      )then
         info = 2
      else if( .not.(diag .eq.'u' ).and.
     $         .not.(diag .eq.'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call f06aaz( 'dtrmv/dtrmv ', info )
         return
      end if

c     quick return if possible.

      if( n.eq.0 ) return

      nounit = (diag.eq.'n' .or. diag.eq.'n')

c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.

      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if

c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.

      if(trans.eq.'n')then

c        form  x := a*x.

         if(uplo.eq.'u')then
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

c        form  x := a'*x.

         if (uplo.eq.'u') then
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
                  if( nounit ) temp = temp*a( j, j )
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

      end

      integer function idamax ( n, x, incx )
c----------------------------------------------------------------------
c BLAS Level 2 function.
c----------------------------------------------------------------------
      implicit none

      integer                  incx, n

      double precision         x( * )

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

      end

      subroutine daxpy ( n, alpha, x, incx, y, incy )
c----------------------------------------------------------------------
c  BLAS Level 2 routine to do the operation y := alpha*x + y
c----------------------------------------------------------------------
      implicit none

      double precision   alpha
      integer            incx, incy, n

      double precision   x( * ), y( * )

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

      end

      subroutine f06aaz ( srname, info )
c----------------------------------------------------------------------
c     .. scalar arguments ..
      integer            info
      character*13       srname

c  f06aaz  is an error handler for the level 2 blas routines.
c
c  it is called by the level 2 blas routines if an input parameter is
c  invalid.

c  srname - character*13.
c           on entry, srname specifies the name of the routine which
c           called f06aaz.

c  info   - integer.
c           on entry, info specifies the position of the invalid
c           parameter in the parameter-list of the calling routine.
      integer            ierr, ifail
      character*4        varbnm
c     .. local arrays ..
      character*80       rec (1)
c     .. external functions ..
      integer            p01acf
      external           p01acf
c----------------------------------------------------------------------
      write (rec (1),99999) srname, info
      if (srname(1:3).eq.'f06') then
         ierr = -1
         varbnm = '    '
      else
         ierr = -info
         varbnm = 'info'
      end if
      ifail = 0
      ifail = p01acf (ifail, ierr, srname(1:6), varbnm, 1, rec)

99999 format ( ' ** on entry to ', a13, ' parameter number ', i2,
     *         ' had an illegal value' )

c     end of f06aaz.

      end

      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,
     *                   beta, y, incy )
c----------------------------------------------------------------------
c BLAS Level 2 routine to do the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c----------------------------------------------------------------------
      implicit none

      double precision   alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans

      double precision   a( lda, * ), x( * ), y( * )

      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   temp, temp1, temp2, temp3, temp4
      integer            i, info, iy, j, jx, kx, ky, lenx, leny, m4, n4
c----------------------------------------------------------------------
c     test the input parameters.

      info = 0
      if     ( .not.(trans.eq.'n' ).and.
     *         .not.(trans.eq.'t' ).and.
     *         .not.(trans.eq.'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call f06aaz( 'dgemv/dgemv ', info )
         return
      end if

      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     *    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) ) return

c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.

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

c     start the operations. in this version the inner loops are all
c     equivalent to axpy operations.

c     first form  y := beta*y.

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
      if( alpha.eq.zero )
     *   return
      jx = kx
      if( trans.eq.'n' )then

c        form  y := alpha*a*x + y.

         if( incy.eq.1 )then

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

      end

      subroutine dtrsv ( uplo, trans, diag, n, a, lda, x, incx )
c----------------------------------------------------------------------
c  dtrsv  solves one of the systems of equations

c     a*x = b,   or   a'*x = b,

c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c----------------------------------------------------------------------
      implicit none

      integer            incx, lda, n
      character*1        diag, trans, uplo

      double precision   a( lda, * ), x( * )

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
c----------------------------------------------------------------------
      info = 0
      if     ( .not.(uplo .eq.'u' ).and.
     *         .not.(uplo .eq.'l' )      )then
         info = 1
      else if( .not.(trans.eq.'n' ).and.
     *         .not.(trans.eq.'t' ).and.
     *         .not.(trans.eq.'c' )      )then
         info = 2
      else if( .not.(diag .eq.'u' ).and.
     *         .not.(diag .eq.'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call f06aaz( 'dtrsv/dtrsv ', info )
         return
      end if

c     quick return if possible.

      if( n.eq.0 ) return

      nounit = (diag.eq.'n')

c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.

      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if

c     start the operations

      if( trans.eq.'n' )then

c        form  x := inv( a )*x.

         if( (uplo.eq.'u' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     *                  x( j ) = x( j )/a( j, j )
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
                     if( nounit )
     *                  x( jx ) = x( jx )/a( j, j )
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
                     if( nounit )
     *                  x( j ) = x( j )/a( j, j )
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
                     if( nounit )
     *                  x( jx ) = x( jx )/a( j, j )
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
                  if( nounit )
     *               temp = temp/a( j, j )
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
                  if( nounit )
     *               temp = temp/a( j, j )
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
                  if( nounit )
     *               temp = temp/a( j, j )
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
                  if( nounit )
     *               temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   - incx
  160          continue
            end if
         end if
      end if

c     end of dtrsv (dtrsv ).

      end

      subroutine f06egf( n, x, incx, y, incy )
c----------------------------------------------------------------------
c  f06egf does the operations temp := x,   x := y,   y := temp.
c----------------------------------------------------------------------
      implicit none
c     .. entry points ..
      entry      dswap ( n, x, incx, y, incy )

      integer            incx, incy, n

      double precision   x( * ), y( * )

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

c     end of f06egf. ( dswap )

      end

      double precision function f06eaf( n, x, incx, y, incy )
c----------------------------------------------------------------------
c  f06eaf returns  x'y
c----------------------------------------------------------------------
      implicit none
c     .. entry points ..
      double precision          ddot
      entry                     ddot  ( n, x, incx, y, incy )

      integer                           incx, incy, n

      double precision                  x( * ), y( * )

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
      f06eaf = sum

c     end of f06eaf. ( ddot )

      end

      subroutine f06pmf( m, n, alpha, x, incx, y, incy, a, lda )
c----------------------------------------------------------------------

c  dger   does the rank 1 operation
c
c     a := alpha*x*y' + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.

      entry      dger  ( m, n, alpha, x, incx, y, incy, a, lda )
c----------------------------------------------------------------------
      implicit none

      double precision   alpha
      integer            incx, incy, lda, m, n

      double precision   a( lda, * ), x( * ), y( * )

      double precision   zero
      parameter        ( zero = 0.0d+0 )

      double precision   temp
      integer            i, info, ix, j, jy, kx
c----------------------------------------------------------------------
c     test the input parameters.

      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call f06aaz( 'f06pmf/dger  ', info )
         return
      end if

      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) ) return

c     start the operations

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

c     end of f06pmf (dger  ).

      end

      subroutine f06eff( n, x, incx, y, incy )
c----------------------------------------------------------------------

c  f06eff does the operation

c     y := x

      entry      dcopy ( n, x, incx, y, incy )
c----------------------------------------------------------------------
      implicit none

      integer            incx, incy, n

      double precision   x( * ), y( * )

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

c     end of f06eff. ( dcopy )

      end

      double precision function f06ejf( n, x, incx )
c----------------------------------------------------------------------
c  f06ejf returns the euclidean norm of a vector via the function
c  name, so that

c     f06ejf := sqrt( x'*x )
c----------------------------------------------------------------------
      implicit none

      double precision          dnrm2
      entry                     dnrm2 ( n, x, incx )
c     .. scalar arguments ..
      integer                           incx, n
c     .. array arguments ..
      double precision                  x( * )

      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      double precision      norm, scale, ssq
c     .. external functions ..
      double precision      f06bmf
      external              f06bmf
c----------------------------------------------------------------------
      if( n.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
         call f06fjf( n, x, incx, scale, ssq )
         norm  = f06bmf( scale, ssq )
      end if

      f06ejf = norm

c     end of f06ejf. ( dnrm2 )

      end

      subroutine dscal( n, alpha, x, incx )
c----------------------------------------------------------------------
c  BLAS Level 2 subroutine to do the operation x := alpha*x
c----------------------------------------------------------------------
      implicit none

c     .. scalar arguments ..
      double precision   alpha
      integer            incx, n
c     .. array arguments ..
      double precision   x( * )
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      integer            ix
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

      end

      subroutine f01qff(pivot,m,n,a,lda,zeta,perm,work,ifail)
c----------------------------------------------------------------------
c  f01qff  finds  a  qr factorization  of  the  real  m by n  matrix  a,
c  incorporating  column interchanges,  so that  a  is reduced to  upper
c  triangular form  by means of  orthogonal transformations  and  column
c  permutations.
c----------------------------------------------------------------------
      implicit none

      double precision  lamda, one, zero
      parameter         (lamda=1.0d-2,one=1.0d+0,zero=0.0d+0)
      character*6       srname
      parameter         (srname='f01qff')
c     .. scalar arguments ..
      integer           ifail, lda, m, n
      character*1       pivot
c     .. array arguments ..
      double precision  a(lda,*), work(*), zeta(*)
      integer           perm(*)
c     .. local scalars ..
      double precision  eps, maxnrm, norm, temp, tol
      integer           ierr, j, jmax, k, la
c     .. local arrays ..
      character*46      rec(1)
c     .. external functions ..
      double precision  dnrm2
      integer           p01abf
      external          dnrm2, p01abf
c----------------------------------------------------------------------
c     check the input parameters.

      ierr = 0
      if ((pivot.ne.'c') 
     *    .and. (pivot.ne.'s')) call p01abw(pivot,'pivot',ifail,ierr,
     *                               srname)
      if (m.lt.0) call p01aby(m,'m',ifail,ierr,srname)
      if (n.lt.0) call p01aby(n,'n',ifail,ierr,srname)
      if (lda.lt.m) call p01aby(lda,'lda',ifail,ierr,srname)
      if (ierr.gt.0) then
         write (rec,fmt=99999) ierr
         ifail = p01abf(ifail,-1,srname,1,rec)
         return
      end if

c     compute eps and the initial column norms.

      if (min(m,n).eq.0) then
         ifail = p01abf(ifail,0,srname,0,rec)
         return
      end if
      eps = 1.11022302462516d-16
      do 20 j = 1, n
         work(j) = dnrm2(m,a(1,j),1)
         work(j+n) = work(j)
   20 continue

c     perform the factorization. tol is the tolerance for f06frf.

      la = lda
      do 120 k = 1, min(m,n)

c        find the pivot column.

         maxnrm = zero
         jmax = k
         if (pivot.eq.'c') then
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

c           use a householder reflection to zero the kth column of a.
c           first set up the reflection.

            call f06frf(m-k,a(k,k),a(k+1,k),1,tol,zeta(k))
            if (k.lt.n) then
               if (zeta(k).gt.zero) then
                  if ((k+1).eq.n) la = m - k + 1

c                 temporarily store beta and put zeta( k ) in a( k, k ).

                  temp = a(k,k)
                  a(k,k) = zeta(k)

c                 perform the operation  a := q( k )*a.

                  call dgemv ('t',m-k+1,n-k,one,a(k,k+1),la,
     *                       a(k,k),1,zero,zeta(k+1),1)

c                 form  b := b - u*work'.

                  call dger(m-k+1,n-k,-one,a(k,k),1,zeta(k+1),1,a(k,k+1)
     *                      ,la)

c                 restore beta.
c
                  a(k,k) = temp
               end if

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

c     set the final  zeta  when  m.le.n.

      if (m.le.n) zeta(m) = zero

      ifail = p01abf(ifail,0,srname,0,rec)

c     end of f01qff. ( sgeqrp )

99999 format ('    the input parameters contained ',i2,' error(s)')
      end

      subroutine f06fdf( n, alpha, x, incx, y, incy )
c----------------------------------------------------------------------
c  f06fdf does the operation y := alpha*x
c----------------------------------------------------------------------
      implicit none

      double precision   alpha
      integer            incx, incy, n
c     .. array arguments ..
      double precision   x( * ), y( * )
      double precision   zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      integer            i, ix, iy
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( ( alpha.eq.zero ).and.( incy.ne.0 ) )then
            call sload ( n, zero, y, abs( incy ) )
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

c     end of f06fdf. ( sscmv )

      end

      subroutine f01qcf(m,n,a,lda,zeta,ifail)
c----------------------------------------------------------------------
c  f01qcf  finds  the  qr factorization  of the real  m by n,  m .ge. n,
c  matrix a,  so that  a is reduced to upper triangular form by means of
c  orthogonal transformations.
c----------------------------------------------------------------------
      implicit none

      double precision  one, zero
      parameter         (one=1.0d+0,zero=0.0d+0)
      character*6       srname
      parameter         (srname='f01qcf')
c     .. scalar arguments ..
      integer           ifail, lda, m, n
c     .. array arguments ..
      double precision  a(lda,*), zeta(*)
c     .. local scalars ..
      double precision  temp
      integer           ierr, k, la
c     .. local arrays ..
      character*46      rec(1)
c     .. external functions ..
      integer           p01abf
      external          p01abf
c----------------------------------------------------------------------
c     check the input parameters.

      ierr = 0
      if (m.lt.n) call p01aby(m,'m',ifail,ierr,srname)
      if (n.lt.0) call p01aby(n,'n',ifail,ierr,srname)
      if (lda.lt.m) call p01aby(lda,'lda',ifail,ierr,srname)
      if (ierr.gt.0) then
         write (rec,fmt=99999) ierr
         ifail = p01abf(ifail,-1,srname,1,rec)
         return
      end if

c     perform the factorization.

      if (n.eq.0) then
         ifail = p01abf(ifail,0,srname,0,rec)
         return
      end if
      la = lda
      do 20 k = 1, min(m-1,n)

c        use a  householder reflection  to  zero the  kth column  of  a.
c        first set up the reflection.

         call f06frf(m-k,a(k,k),a(k+1,k),1,zero,zeta(k))

         if ((zeta(k).gt.zero) .and. (k.lt.n)) then
            if ((k+1).eq.n) la = m - k + 1

c           temporarily  store  beta and  put  zeta( k )  in  a( k, k ).

            temp = a(k,k)
            a(k,k) = zeta(k)

c           perform the operation  a := q( k )*a.

            call dgemv ('t',m-k+1,n-k,one,a(k,k+1),la,a(k,k),1,
     *                 zero,zeta(k+1),1)

c           form  b := b - u*work'.

            call dger(m-k+1,n-k,-one,a(k,k),1,zeta(k+1),1,a(k,k+1),la)

c           restore beta.

            a(k,k) = temp
         end if
   20 continue

c     set the final  zeta  when  m.eq.n.

      if (m.eq.n) zeta(n) = zero

      ifail = p01abf(ifail,0,srname,0,rec)

c     end of f01qcf. ( sgeqr )

99999 format ('    the input parameters contained ',i2,' error(s)')
      end

      double precision function f06rjf(norm,uplo,diag,m,n,a,lda,work)
c----------------------------------------------------------------------
c  dlantr  returns the value of the one norm,  or the frobenius norm, or
c  the  infinity norm,  or the  element of  largest absolute value  of a
c  trapezoidal or triangular matrix a.
c----------------------------------------------------------------------
      implicit none

      double precision                 one, zero
      parameter                        (one=1.0d+0,zero=0.0d+0)

      integer                          lda, m, n
      character                        diag, norm, uplo

      double precision                 a(lda,*), work(*)

      double precision                 scale, sum, value
      integer                          i, j
      logical                          udiag
c----------------------------------------------------------------------
      if (min(m,n).eq.0) then
         value = zero
      else if ( norm.eq.'m' ) then

c        find max(abs(a(i,j))).

         if (diag.eq.'u') then
            value = one
            if (uplo.eq.'u') then
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
            if (uplo.eq.'u') then
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
      else if ((norm.eq.'o') .or. (norm.eq.'1')) then

c        find norm1(a).

         value = zero
         udiag = (diag.eq.'u')
         if (uplo.eq.'u') then
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
      else if (norm.eq.'i') then

c        find normi(a).

         if (uplo.eq.'u') then
            if (diag.eq.'u') then
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
            if ((diag.eq.'u' .or. diag.eq.'u')) then
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

         if (uplo.eq.'u') then
            if (diag.eq.'u') then
               scale = one
               sum = min(m,n)
               do 580 j = 2, n
                  call f06fjf(min(m,j-1),a(1,j),1,scale,sum)
  580          continue
            else
               scale = zero
               sum = one
               do 600 j = 1, n
                  call f06fjf(min(m,j),a(1,j),1,scale,sum)
  600          continue
            end if
         else
            if (diag.eq.'u') then
               scale = one
               sum = min(m,n)
               do 620 j = 1, n
                  call f06fjf(m-j,a(min(m,j+1),j),1,scale,sum)
  620          continue
            else
               scale = zero
               sum = one
               do 640 j = 1, n
                  call f06fjf(m-j+1,a(j,j),1,scale,sum)
  640          continue
            end if
         end if
         value = scale*sqrt(sum)
      end if

      f06rjf = value

c     end of f06rjf (dlantr)

      end

      subroutine f06bcf( t, c, s )
c----------------------------------------------------------------------
c  f06bcf returns values c and s such that
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
c----------------------------------------------------------------------
      implicit none

      double precision   c, s, t

      double precision   one
      parameter        ( one = 1.0d+0 )

      double precision   abst, eps, reps, rrteps, rteps
      logical            first

      save               first, eps, reps, rteps, rrteps

      data               first/ .true. /
c----------------------------------------------------------------------
      if( first )then
         first  = .false.
         eps    =  1.11022302462516d-16
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

c     end of f06bcf. ( scsg )

      end

      subroutine f06fjf( n, x, incx, scale, sumsq )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      double precision   scale, sumsq
      integer            incx, n

      double precision   x( * )

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

c     end of f06fjf. ( sssq )

      end

      subroutine f06flf( n, x, incx, xmax, xmin )
c----------------------------------------------------------------------
c  f06flf returns the values xmax and xmin given by

c     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).

c  if n is less than unity then xmax and xmin are returned as zero.
c----------------------------------------------------------------------
      implicit none

      double precision   xmax, xmin
      integer            incx, n

      double precision   x( * )
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

c     end of f06flf. ( scond )

      end

      subroutine f06dbf( n, const, x, incx )
c----------------------------------------------------------------------
c  f06dbf does the operation x = const*e,   e' = ( 1  1 ... 1 ).
c----------------------------------------------------------------------
      implicit none

      integer            const, incx, n

      integer            x( * )
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

c     end of f06dbf. ( iload )

      end

      integer function f06klf( n, x, incx, tol )
c----------------------------------------------------------------------
c  f06klf finds the first element of the n element vector x for which
c
c     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
c
c  and returns the value ( k - 1 ) in the function name f06klf. if no
c  such k exists then f06klf is returned as n.
c
c  if tol is supplied as less than zero then the value epsmch, where
c  epsmch is the relative machine precision, is used in place of tol.
c----------------------------------------------------------------------
      implicit none

c     .. scalar arguments ..
      double precision         tol
      integer                  incx, n
c     .. array arguments ..
      double precision         x( * )

      double precision         zero
      parameter              ( zero = 0.0d+0 )
c     .. local scalars ..
      double precision         tl, xmax
      integer                  ix, k
c----------------------------------------------------------------------
      k = 0
      if( n.ge.1 )then
         ix = 1
         if( tol.lt.zero )then
            tl = 1.11022302462516d-16
         else
            tl = tol
         end if
         xmax = abs( x( ix ) )

   10    if   ( k.lt.n )then
            if( abs( x( ix ) ).le.tl*xmax )
     $         go to 20
            xmax = max( xmax, abs( x( ix ) ) )
            k    = k  + 1
            ix   = ix + incx
            go to 10
         end if

      end if

   20 f06klf = k

c     end of f06klf. ( isrank )

      end

      subroutine f06qtf( side, n, k1, k2, c, s, a, lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, fill, stemp, temp
      integer            i, i1, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.(k2.le.k1).or.( k2.gt.n )) return

      if( side.eq.'l' )then

         do 20 j = k1 + 1, n

            aij = a( k1, j )
            do 10 i = k1, min( j - 1, k2 - 1 )
               a( i, j ) = s( i )*a( i + 1, j ) + c( i )*aij
               aij = c( i )*a( i + 1, j ) - s( i )*aij
   10       continue
            a( i, j ) = aij
   20    continue

         do 40 j = k1, k2 - 1

            fill = -s( j )*a( j, j )
            a( j, j ) = c( j )*a( j, j )

            call f06baf( a( j + 1, j + 1 ), fill, ctemp, stemp )
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

      else if( side.eq.'r' )then

         do 60 j = k2 - 1, k1, -1

            if( ( c( j ).ne.one ).or.( s( j ).ne.zero ) )then
               ctemp = c( j )
               stemp = s( j )
               do 50 i = 1, j
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
   50          continue

               fill = s( j )*a( j + 1, j + 1 )
               a( j + 1, j + 1 ) = c( j )*a( j + 1, j + 1 )

               call f06baf( a( j, j ), fill, c( j ), s( j ) )
            end if
   60    continue

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

c     end of f06qtf. ( sutsqr )

      end

      subroutine f06fcf( n, d, incd, x, incx )
c----------------------------------------------------------------------
c  f06fcf does the operation x := diag( d )*x
c----------------------------------------------------------------------
      implicit none

      integer            incd, incx, n

      double precision   d( * ), x( * )

      integer            i, id, ix
c----------------------------------------------------------------------
      if( n.gt.0 )then
         if( ( incd.eq.0 ).and.( incx.ne.0 ) )then

            call dscal ( n, d( 1 ), x, abs( incx ) )

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

c     end of f06fcf. ( sdscl )

      end

      subroutine f06fqf( pivot, direct, n, alpha, x, incx, c, s )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      double precision   alpha
      integer            incx, n
      character*1        direct, pivot
      double precision   c( * ), s( * ), x( * )
      integer            i, ix
c----------------------------------------------------------------------
      if( n.gt.0 )then

         if( direct.eq.'b' )then

            ix = 1 + ( n - 1 )*incx

            if( pivot.eq.'v' )then

               do 10, i = n, 2, -1
                  call f06baf( x( ix - incx ), x( ix ), c( i ), s( i ) )
                  ix = ix - incx
   10          continue

               call f06baf( alpha, x( ix ), c( 1 ), s( 1 ) )

            else if( pivot.eq.'f')then

               do 20, i = n, 1, -1

                  call f06baf( alpha, x( ix ), c( i ), s( i ) )
                  s( i )  = -s( i )
                  x( ix ) = -x( ix )
                  ix      =  ix      - incx

   20          continue

            end if

         else if( direct.eq.'f' )then

            ix = 1

            if( pivot.eq.'v' )then

               do 30, i = 1, n - 1
                  call f06baf( x( ix + incx ), x( ix ), c( i ), s( i ) )
                  s( i )  = -s( i )
                  x( ix ) = -x( ix )
                  ix      =  ix      + incx
   30          continue

               call f06baf( alpha, x( ix ), c( n ), s( n ) )
               s( n )  = -s( n )
               x( ix ) = -x( ix )

            else if( pivot.eq.'f' )then

               do 40, i = 1, n
                  call f06baf( alpha, x( ix ), c( i ), s( i ) )
                  ix = ix + incx
   40          continue

            end if

         end if

      end if

c     end of f06fqf. ( ssrotg )

      end

      double precision function adivb ( a, b, fail )
c----------------------------------------------------------------------
c BLAS Level 2 routine (named sdiv) to do scalar division.

c  adivb returns the value div given by
c
c     div = ( a/b                 if a/b does not overflow,
c           (
c           ( 0.0                 if a .eq. 0.0,
c           (
c           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
c
c  where  flmax  is a large value, via the function name. in addition if
c  a/b would overflow then  fail is returned as true, otherwise  fail is
c  returned as false.
c
c  note that when  a and b  are both zero, fail is returned as true, but
c  div  is returned as  0.0. in all other cases of overflow  div is such
c  that  abs( div ) = flmax.
c
c  when  b = 0  then  sign( a/b )  is taken as  sign( a ).
c----------------------------------------------------------------------
      double precision                  a, b
      logical                           fail

      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      double precision      absb, div, flmax, flmin
      logical               first
c     .. save statement ..
      save                  first, flmin, flmax
c     .. data statements ..
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

         if( first )then
            first = .false.
            flmin =  2.22507385850721d-308
            flmax =  1/flmin
         end if

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

      adivb = div

      end

      subroutine f06qsf( side, n, k1, k2, c, s, a, lda )
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, n
      character*1        side

      double precision   a( lda, * ), c( * ), s( * )

      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   aij, ctemp, spike, stemp, temp
      integer            i, j
c----------------------------------------------------------------------
      if((min(n,k1).lt.1).or.(k2.le.k1).or.(k2.gt.n)) return

      if( side.eq.'l' )then

         do 20 j = k1, k2 - 1
            spike = s( j )
            do 10 i = k1, j - 1
               aij = a( i, j )
               a( i, j ) = s( i )*spike + c( i )*aij
               spike = c( i )*spike - s( i )*aij
   10       continue
            call f06baf( a( j, j ), spike, c( j ), s( j ) )
   20    continue

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

         do 70 j = k2, k1 + 1, -1

            call f06baf( a( j, j ), s( j - 1 ), ctemp, stemp )
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

c     end of f06qsf. ( susqr )

      end

      subroutine f06qxf( side, pivot, direct, m, n, k1, k2, c, s, a,
     *                   lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, m, n
      character*1        direct, pivot, side

      double precision   a( lda, * ), c( * ), s( * )
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
            if( direct.eq.'f' )then
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
         else if( pivot.eq.'t')then
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
      else if ( right )then
         if ( pivot.eq.'v' )then
            if (direct.eq.'f' )then
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
            else if( direct.eq.'b')then
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
            else if( ( direct.eq.'b' ).or.( direct.eq.'b' ) )then
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

c     end of f06qxf. ( sgesrc )

      end

      subroutine f06dff( n, x, incx, y, incy )
c----------------------------------------------------------------------
c  f06dff does the operation y := x
c----------------------------------------------------------------------
      implicit none

      integer            incx, incy, n
      integer            x( * ), y( * )
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

c     end of f06dff. ( icopy )

      end

      double precision function f06bmf( scale, ssq )
c----------------------------------------------------------------------
c  f06bmf returns the value norm given by
c
c     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
c            (
c            ( flmax,             scale*sqrt( ssq ) .ge. flmax
c----------------------------------------------------------------------
      implicit none

      double precision                  scale, ssq
      double precision      flmax, flmin, norm, sqt
      logical               first

      save                  first, flmax

      data                  first/ .true. /
c----------------------------------------------------------------------
      if( first )then
         first = .false.
         flmin =  2.22507385850721d-308
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
      f06bmf = norm

c     end of f06bmf. ( snorm )

      end

      subroutine f06qvf( side, n, k1, k2, c, s, a, lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      integer            k1, k2, lda, n
      character*1        side
      double precision   a( lda, * ), c( * ), s( * )
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
      double precision   aij, ctemp, stemp, temp
      integer            i, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return

      if(side.eq.'l')then
         do 20 j = n, k1, -1
            if( j.ge.k2 )then
               aij = a( k2, j )
            else
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

      else if(side.eq.'r')then

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

c     end of f06qvf. ( sutsrh )

      end

      subroutine f06qkf( side, trans, n, perm, k, b, ldb )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      integer            k, ldb, n
      character*1        side, trans
      double precision   perm( * ), b( ldb, * )

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

c     end of f06qkf. ( sgeapr )

      end

      subroutine f06qnz(side,n,k1,k2,s,a,lda)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      double precision  zero
      parameter         (zero=0.0d+0)
      integer           k1, k2, lda, n
      character*1       side
      double precision  a(lda,*), s(*)
      double precision  aij, temp
      integer           i, j
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return
      if ((side.eq.'l') .or. (side.eq.'l')) then

c        apply the permutations to columns n back to k1.

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
      else if ((side.eq.'r') .or. (side.eq.'r')) then

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

c     end of f06qnz. ( sutsrh )

      end

      subroutine f06frf( n, alpha, x, incx, tol, zeta )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      double precision   alpha, tol, zeta
      integer            incx, n
c     .. array arguments ..
      double precision   x( * )
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      double precision   beta, eps, scale, ssq
      logical            first
c     .. save statement ..
      save               eps, first
c     .. data statements ..
      data               first/ .true. /
c----------------------------------------------------------------------
      if( n.lt.1 )then
         zeta = zero
      else if( ( n.eq.1 ).and.( x( 1 ).eq.zero ) )then
         zeta = zero
      else

         if( first )then
            first = .false.
            eps   =  1.11022302462516d-16
         end if

c        treat case where p is a 2 by 2 matrix specially.

         if( n.eq.1 )then

c           deal with cases where  alpha = zero  and
c           abs( x( 1 ) ) .le. max( eps*abs( alpha ), tol )  first.

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
            call f06fjf( n, x, incx, scale, ssq )

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
               call dscal ( n, -1/alpha, x, incx )
            else
               if( scale.lt.abs( alpha ) )then
                  beta = abs( alpha )*sqrt( 1 + ssq*( scale/alpha )**2 )
               else
                  beta = scale       *sqrt( ssq +   ( alpha/scale )**2 )
               end if
               zeta = sqrt( ( beta + abs( alpha ) )/beta )
               if( alpha.gt.zero ) beta = -beta
               call dscal ( n, -1/( zeta*beta ), x, incx )
               alpha = beta
            end if
         end if
      end if

c     end of f06frf. ( sgrfg )

      end

      subroutine f06qzz(hess,n,k1,k2,c,s,a,lda)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      double precision  one, zero
      parameter         (one=1.0d+0,zero=0.0d+0)
      integer           k1, k2, lda, n
      character*1       hess
      double precision  a(lda,*), c(*), s(*)
      double precision  ctemp, stemp, suph, temp
      integer           i, j
c----------------------------------------------------------------------
      if ((min(n,k1).lt.1) .or. (k2.le.k1) .or. (k2.gt.n)) return

      if (hess.eq.'c') then

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

         do 80 j = k2 - 1, k1, -1
            suph = s(j)
            call f06baf(a(n-j,j+1),suph,ctemp,stemp)
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

c     end of f06qzz.

      end

      subroutine f06qwf( side, n, k1, k2, c, s, a, lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none
      integer            k1, k2, lda, n
      character*1        side
      double precision   a( lda, * ), c( * ), s( * )
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
      double precision   aij, ctemp, spike, stemp, temp
      integer            i, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return

      if( side.eq.'l' )then

         do 20 j = n, k2, -1
            temp = a( k2, j )
            do 10 i = k2 - 1, k1, -1
               aij = a( i, j )
               a( i, j ) = s( i )*temp + c( i )*aij
               temp = c( i )*temp - s( i )*aij
   10       continue
            a( k2, j ) = temp
   20    continue

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

c     end of f06qwf. ( sutsrs )

      end

      subroutine f06qhf( matrix, m, n, const, diag, a, lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none

      character*1        matrix
      double precision   const, diag
      integer            lda, m, n
      double precision   a( lda, * )
      integer            i, j
c----------------------------------------------------------------------
      if( ( matrix.eq.'g' ).or.( matrix.eq.'g' ) )then
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
      else if( ( matrix.eq.'u' ).or.( matrix.eq.'u' ) )then
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
      else if( ( matrix.eq.'l' ).or.( matrix.eq.'l' ) )then
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

c     end of f06qhf. ( smload )

      end

      subroutine f06qrf( side, n, k1, k2, c, s, a, lda )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none
      integer            k1, k2, lda, n
      character*1        side
      double precision   a( lda, * ), c( * ), s( * )
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
      double precision   aij, ctemp, stemp, subh, temp
      integer            i, j
c----------------------------------------------------------------------
      if( ( min( n, k1 ).lt.1 ).or.( k2.le.k1 ).or.
     *   ( k2.gt.n ) )return

      if( side.eq.'l' )then

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
               call f06baf( aij, subh, c( j ), s( j ) )
               a( j, j ) = aij
            else
               a( k2, j ) = aij
            end if
   20    continue

      else if( side.eq.'r' )then

         do 40 j = k2 - 1, k1, -1
            subh = s( j )
            call f06baf( a( j + 1, j + 1 ), subh, ctemp, stemp )
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

c     end of f06qrf. ( suhqr )

      end

      subroutine f06qff( matrix, m, n, a, lda, b, ldb )
c----------------------------------------------------------------------

c----------------------------------------------------------------------
      implicit none
      character*1        matrix
      integer            m, n, lda, ldb
      double precision   a( lda, * ), b( ldb, * )
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

c     end of f06qff. ( smcopy )

      end

      subroutine f06baf( a, b, c, s )
c----------------------------------------------------------------------
c  version of srotg blas,  except that c is always
c  returned as non-negative and  b  is overwritten by the tangent of the
c  angle that defines the plane rotation.
c----------------------------------------------------------------------
      implicit none

      double precision   a, b, c, s

      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )

      double precision   t
      logical            fail
c     .. external functions ..
      double precision   adivb
      external           adivb
c----------------------------------------------------------------------
      if( b.eq.zero )then
         c  = one
         s  = zero
      else
         t  = adivb( b, a, fail )
         call f06bcf( t, c, s )
         a  = c*a + s*b
         b  = t
      end if

c     end of f06baf. ( srotgc )

      end

      subroutine sload ( n, const, x, incx )
c----------------------------------------------------------------------
c sload - BLAS level 2 routine does the operation 
c x = const*e,   e' = ( 1  1 ... 1 ).
c----------------------------------------------------------------------
      implicit none
      double precision   const
      integer            incx, n
      double precision   x( * )
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

      end

      subroutine p01aby(n,name,inform,ierr,srname)
c----------------------------------------------------------------------
c     p01aby increases the value of ierr by 1 and, if
c
c        ( mod( inform, 10 ).ne.1 ).or.( mod( inform/10, 10 ).ne.0 )
c
c     writes a message on the current error message channel giving the
c     value of n, a message to say that n is invalid and the strings
c     name and srname.
c
c     name must be the name of the actual argument for n and srname must
c     be the name of the calling routine.
c
c     this routine is intended for use when n is an invalid input
c     parameter to routine srname. for example
c
c        ierr = 0
c        if( n.lt.1 )call p01aby( n, 'n', idiag, ierr, srname )

      integer           ierr, inform, n
      character*(*)     name, srname
c     .. local scalars ..
      integer           nerr
c     .. local arrays ..
      character*65      rec(2)
c----------------------------------------------------------------------
      ierr = ierr + 1
      if ((mod(inform,10).ne.1) .or. (mod(inform/10,10).ne.0)) then
         call x04aaf(0,nerr)
         write (rec,fmt=99999) name, srname, n
         call x04baf(nerr,' ')
         call x04baf(nerr,rec(1))
         call x04baf(nerr,rec(2))
      end if

99999 format (' *****  parameter  ',a,'  is invalid in routine  ',a,
     *  '  ***** ',/8x,'value supplied is ',i6)
      end

      integer function p01abf(ifail,ierror,srname,nrec,rec)
c----------------------------------------------------------------------
c     p01abf is the error-handling routine for the library.
c
c     p01abf either returns the value of ierror through the routine
c     name (soft failure), or terminates execution of the program
c     (hard failure). diagnostic messages may be output.
c
c     if ierror = 0 (successful exit from the calling routine),
c     the value 0 is returned through the routine name, and no
c     message is output
c
c     if ierror is non-zero (abnormal exit from the calling routine),
c     the action taken depends on the value of ifail.
c
c     ifail =  1: soft failure, silent exit (i.e. no messages are
c                 output)
c     ifail = -1: soft failure, noisy exit (i.e. messages are output)
c     ifail =-13: soft failure, noisy exit but standard messages from
c                 p01abf are suppressed
c     ifail =  0: hard failure, noisy exit

c     a = 0: hard failure  a = 1: soft failure
c     b = 0: silent exit   b = 1: noisy exit

      integer                 ierror, ifail, nrec
      character*(*)           srname
c     .. array arguments ..
      character*(*)           rec(*)
c     .. local scalars ..
      integer                 i, nerr
      character*72            mess
c----------------------------------------------------------------------
      if (ierror.ne.0) then
c        abnormal exit from calling routine
         if (ifail.eq.-1 .or. ifail.eq.0 .or. ifail.eq.-13 .or.
     *       (ifail.gt.0 .and. mod(ifail/10,10).ne.0)) then
c           noisy exit
            call x04aaf(0,nerr)


         end if
      end if

      p01abf = ierror

99999 format (' ** abnormal exit from library routine ',a,': ifail',
     *  ' =',i6)
      end

      integer function p01acf(ifail,ierror,srname,varbnm,nrec,rec)
c----------------------------------------------------------------------
c     p01acf either returns the value of ierror through the routine
c     name (soft failure), or terminates execution of the program
c     (hard failure). diagnostic messages may be output.
c
c     if ierror = 0 (successful exit from the calling routine),
c     the value 0 is returned through the routine name, and no
c     message is output
c
c     if ierror is non-zero (abnormal exit from the calling routine),
c     the action taken depends on the value of ifail.
c
c     ifail =  1: soft failure, silent exit (i.e. no messages are
c                 output)
c     ifail = -1: soft failure, noisy exit (i.e. messages are output)
c     ifail =-13: soft failure, noisy exit but standard messages from
c                 p01acf are suppressed
c     ifail =  0: hard failure, noisy exit

c     a = 0: hard failure  a = 1: soft failure
c     b = 0: silent exit   b = 1: noisy exit

      integer                 ierror, ifail, nrec
      character*(*)           srname, varbnm
c     .. array arguments ..
      character*(*)           rec(*)
c     .. local scalars ..
      integer                 i, nerr, varlen
      character*72            mess
c----------------------------------------------------------------------
      if (ierror.ne.0) then
         varlen = 0
         do 20 i = len(varbnm), 1, -1
            if (varbnm(i:i).ne.' ') then
               varlen = i
               go to 40
            end if
   20    continue
   40    continue
c        abnormal exit from calling routine
         if (ifail.eq.-1 .or. ifail.eq.0 .or. ifail.eq.-13 .or.
     *       (ifail.gt.0 .and. mod(ifail/10,10).ne.0)) then
c           noisy exit
            call x04aaf(0,nerr)
            do 60 i = 1, nrec
               call x04baf(nerr,rec(i))
   60       continue
            if (ifail.ne.-13) then
               if (varlen.ne.0) then
                  write (mess,fmt=99999) srname, varbnm(1:varlen),
     *              ierror
               else
                  write (mess,fmt=99998) srname
               end if
               call x04baf(nerr,mess)
               call x04baf(nerr,
     *                        ' ** soft failure - control returned')
            end if
         end if
      end if
      p01acf = ierror

99999 format (' ** abnormal exit from library routine ',a,': ',a,
     *       ' =',i6)
99998 format (' ** abnormal exit from library routine ',a)
      end

      subroutine p01abw(n,name,inform,ierr,srname)

c     p01abw increases the value of ierr by 1 and, if

c        ( mod( inform, 10 ).ne.1 ).or.( mod( inform/10, 10 ).ne.0 )

c     writes a message on the current error message channel giving the
c     value of n, a message to say that n is invalid and the strings
c     name and srname.

c     name must be the name of the actual argument for n and srname must
c     be the name of the calling routine.

      integer           ierr, inform
      character*(*)     n
      character*(*)     name, srname
c     .. local scalars ..
      integer           nerr
c     .. local arrays ..
      character*65      rec(3)
c----------------------------------------------------------------------
      ierr = ierr + 1
      if ((mod(inform,10).ne.1) .or. (mod(inform/10,10).ne.0)) then
         call x04aaf(0,nerr)
         write (rec,fmt=99999) name, srname, n
         call x04baf(nerr,' ')
         call x04baf(nerr,rec(1))
         call x04baf(nerr,rec(2))
         call x04baf(nerr,rec(3))
      end if

c     end of p01abw.

99999 format (' *****  parameter  ',a,'  is invalid in routine  ',a,
     *  '  ***** ',/8x,'value supplied is',/8x,a)
      end

      subroutine x04abf(i,nadv)
c----------------------------------------------------------------------
c      if i = 0, sets nadv to current advisory message unit number
c     (stored in nadv1).
c     if i = 1, changes current advisory message unit number to
c     value specified by nadv.
c
c     .. scalar arguments ..
      integer           i, nadv
c     .. local scalars ..
      integer           nadv1
c     .. save statement ..
      save              nadv1
c     .. data statements ..
      data              nadv1/6/
c----------------------------------------------------------------------
      if (i.eq.0) nadv = nadv1
      if (i.eq.1) nadv1 = nadv

      end

      subroutine x04aaf(i,nerr)
c----------------------------------------------------------------------
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     .. scalar arguments ..
      integer           i, nerr
c     .. local scalars ..
      integer           nerr1
c     .. save statement ..
      save              nerr1
c     .. data statements ..
      data              nerr1/0/
c----------------------------------------------------------------------
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr

      end

      subroutine x04bay(nout,nrec,rec)
c----------------------------------------------------------------------
c     x04bay outputs nrec records on device nout, by calling x04baf.
c     if nrec is 0 then no records are output.
c
c     .. scalar arguments ..
      integer           nout, nrec
c     .. array arguments ..
      character*(*)     rec(*)
c     .. local scalars ..
      integer           i
c----------------------------------------------------------------------
      do 20 i = 1, nrec
         call x04baf(nout,rec(i))
   20 continue

c     end of x04bay.

      end

      subroutine x04baf(nout,rec)
c----------------------------------------------------------------------
c     x04baf writes the contents of rec to the unit defined by nout.
c     trailing blanks are not output, except that if rec is entirely
c     blank, a single blank character is output.
c     if nout.lt.0, i.e. if nout is not a valid fortran unit identifier,
c     then no output occurs.
c----------------------------------------------------------------------
      integer           nout
      character*(*)     rec
c     .. local scalars ..
      integer           i
c----------------------------------------------------------------------
      if (nout.ge.0) then
c        remove trailing blanks
         do 20 i = len(rec), 2, -1
            if (rec(i:i).ne.' ') go to 40
   20    continue
c        write record to external file
   40    write (nout,fmt=99999) rec(1:i)
      end if

99999 format (a)
      end

      subroutine e04mhf(string)
c----------------------------------------------------------------------
c     e04mhf  loads the option supplied in  string  into the relevant
c     element of  iprmlc  or  rprmlc.

      character*(*)     string
c     .. scalars in common ..
      logical           newopt
c     .. local scalars ..
      integer           nout
      logical           first, prnt
      character*16      key
      character*72      buffer
c     .. local arrays ..
      character*80      rec(5)
c     .. common blocks ..
      double precision wmach
      common/ cstmch /wmach(10)

      common            /be04mf/newopt
c     .. save statement ..
      save              /be04mf/, first, nout, prnt
c     .. data statements ..
      data              first/.true./
c----------------------------------------------------------------------
c     if first time in, set  nout.
c     newopt  is true first time into  e04mgf  or  e04mhf
c     and just after a call to the main routine (e.g. e04mff).
c     prnt  is set to  true  whenever  newopt  is true.

      if (first) then
         first = .false.
         newopt = .true.
         nout = 6
      end if
      buffer = string

c     call  e04mfx   to decode the option and set the parameter value.
c     if  newopt  is true, reset  prnt  and test specially for nolist.

      if (newopt) then
         newopt = .false.
         prnt = .true.
         call e04mfx(nout,buffer,key)
c
         if (key.eq.'nolist') then
            prnt = .false.
         else
            write (rec,fmt='(// a / a /)') ' calls to e04mhf',
     *        ' ---------------'
            call x04bay(nout,5,rec)
            write (rec,fmt='( 6x, a )') buffer
            call x04baf(nout,rec(1))
         end if
      else
         if (prnt) then
            write (rec,fmt='( 6x, a )') buffer
            call x04baf(nout,rec(1))
         end if
         call e04mfx(nout,buffer,key)

         if (key.eq.'list') prnt = .true.
         if (key.eq.'nolist') prnt = .false.
      end if

c     end of  e04mhf.  (lpprm)

      end

      subroutine e04uef(string)
c----------------------------------------------------------------------
c     e04uef  loads the option supplied in string into the relevant
c     element of iprmls, rprmls, iprmnp or rprmnp.

      character*(*)     string
c     .. scalars in common ..
      logical           newopt

      double precision wmach
      common/ cstmch /wmach(10)
c     .. local scalars ..
      integer           nout
      logical           first, prnt
      character*16      key
      character*72      buffer
c     .. local arrays ..
      character*80      rec(5)
c     .. common blocks ..
      common            /ee04uc/newopt
c     .. save statement ..
      save prnt
      save              /ee04uc/
c     .. data statements ..
      data              first/.true./
c----------------------------------------------------------------------
c     if first time in, set nout.
c     newopt is true first time into e04udf or e04uef
c     and just after a call to an optimization routine.
c     prnt is set to true whenever newopt is true.

      if (first) then
         first = .false.
         newopt = .true.
      end if

      nout = 6
      buffer = string

c     call e04ucq to decode the option and set the parameter value.
c     if newopt is true, reset prnt and test specially for nolist.

      if (newopt) then
         newopt = .false.
         prnt = .true.
         call e04ucq(nout,buffer,key)

         if (key.eq.'nolist') then
            prnt = .false.
         else
            write (rec,fmt='(// a / a /)') ' calls to e04uef',
     *        ' ---------------'
            call x04bay(nout,5,rec)
            write (rec,fmt='( 6x, a )') buffer
            call x04baf(nout,rec(1))
         end if
      else
         if (prnt) then
            write (rec,fmt='( 6x, a )') buffer
            call x04baf(nout,rec(1))
         end if
         call e04ucq(nout,buffer,key)

         if (key.eq.'list') prnt = .true.
         if (key.eq.'nolist') prnt = .false.
      end if

c     end of  e04uef. (npoptn)

      end
