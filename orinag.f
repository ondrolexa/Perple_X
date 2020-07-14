      SUBROUTINE lpsol (N,NCLIN,A,LDA,BL,BU,CVEC,ISTATE,X,ITER,OBJ,AX,
     *                  CLAMDA,IW,LENIW,W,LENW,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFF  solves the linear programming problem
C
C           minimize               c' x
C              x
C                                 (  x )
C           subject to    bl  .le.(    ).ge.  bu,
C                                 ( Ax )
C
C     where  A  is a constant  nclin by n  matrix.
C     The feasible region is defined by a mixture of linear equality or
C     inequality constraints on  x.
C
C     n  is the number of variables (dimension of x).
C        (n must be positive.)
C
C     nclin  is the number of general linear constraints (rows of  A).
C        (nclin may be zero.)
C
C     The first  n  elements of  bl  and   bu  are lower and upper
C     bounds on the variables.  The next  nclin  elements are
C     lower and upper bounds on the general linear constraints.
C
C     The matrix  A  of coefficients in the general linear constraints
C     is entered as the two-dimensional array  A  (of dimension
C     LDA by n).  If nclin = 0, A is not referenced.
C
C     The vector  x  must contain an initial estimate of the solution,
C     and will contain the computed solution on output.
C
C     Documentation for  E04MFF  is coming real soon now.
C     Wait for the release of  users guide for LPOPT (Version 1.00-6),
C     by P. E. Gill, W. Murray and M. A. Saunders,
C
C     Version 1.0-6    Jun  30, 1991.  (Nag Mk 16 version.)
C
C     Copyright  1989  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04MFF')
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, POINT3
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1)
      DOUBLE PRECISION  POINT8, POINT9, ONE
      PARAMETER         (POINT8=0.8D+0,POINT9=0.9D+0,ONE=1.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ
      INTEGER           IFAIL, ITER, LDA, LENIW, LENW, N, NCLIN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CLAMDA(N+NCLIN), CVEC(*), W(LENW), X(N)
      INTEGER           ISTATE(N+NCLIN), IW(LENIW)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ALFA, ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9,
     *                  TOLACT, TOLFEA, TOLINC, TOLRNK, TOLX0, TRULAM
      INTEGER           IDBGLC, IPRINT, IPRNT, ISDEL, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITNFIX, JADD, JDEL, KCHK,
     *                  KCYCLE, KDEGEN, LCRASH, LDBGLC, LDQ, LDT,
     *                  LENNAM, LINES1, LINES2, LPROB, MAXACT, MAXNZ,
     *                  MM, MSGLC, MXFREE, NCOLT, NDEGEN, NN, NNCLIN,
     *                  NOUT, NPROB
      LOGICAL           HEADER, LCDBG, NEWOPT, PRNT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           ILCDBG(LDBG), IPADLC(14), IPSVLC(MXPARM),
     *                  LOCLC(LENLC), NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, CONDMX, EPSMCH, ERRMAX, FEAMAX, FEAMIN,
     *                  RTEPS, XNORM
      INTEGER           IANRMJ, IDBG, INFORM, IT, ITMAX, J, JINF, JMAX,
     *                  LANORM, LCQ, LD, LDH, LDR, LFEATU, LGQ, LITOTL,
     *                  LKACTV, LKX, LLPTYP, LQ, LR, LRLAM, LT, LWRK,
     *                  LWTINF, LWTOTL, MINACT, MINFXD, MSGDBG, MSGLVL,
     *                  NACT1, NACTIV, NARTIF, NCNLN, NCTOTL, NERR,
     *                  NERROR, NFREE, NGQ, NMOVED, NREJTD, NRZ, NUMINF,
     *                  NVIOL, NZ
      LOGICAL           COLD, CSET, DONE, FOUND, HALTED, HOT, NAMED,
     *                  ROWERR, RSET, UNITQ, VERTEX, WARM
      CHARACTER*2       PRBTYP
      CHARACTER*4       START
      CHARACTER*6       MSG
      CHARACTER*11      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      ERRREC(2), REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           P01ABF
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04MFJ, E04MFK, E04MFP, E04MFR,
     *                  E04MFT, E04MFU, E04MFV, E04MFW, E04MFZ, E04NBW,
     *                  E04NFQ, F06DFF, F06FBF, F06FLF, X02ZAZ, X04BAF,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04MF/NEWOPT
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /CE04MF/TOLX0, TOLINC, KDEGEN, NDEGEN, ITNFIX,
     *                  NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04MF/ALFA, TRULAM, ISDEL, JDEL, JADD, HEADER,
     *                  PRNT
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /EE04MF/ILCDBG, LCDBG
      COMMON            /FE04MF/IPSVLC, IDBGLC, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, KCHK, KCYCLE, LCRASH, LPROB, MAXACT,
     *                  MXFREE, MAXNZ, MM, LDBGLC, MSGLC, NN, NNCLIN,
     *                  NPROB, IPADLC
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLC(1),IDBGLC), (RPRMLC(1),BIGBND)
      EQUIVALENCE       (MSGLC,MSGLVL), (IDBGLC,IDBG), (LDBGLC,MSGDBG)
C     .. Save statement ..
      SAVE              /AX02ZA/, /AE04MF/, /BE04MF/, /FE04MF/, /GE04MF/
C     .. Data statements ..
      DATA              TITLE/' *** E04MFF'/
C     .. Executable Statements ..
C
C     Set the machine-dependent constants.
C
      CALL X02ZAZ
C
      EPSMCH = WMACH(3)
      RTEPS = WMACH(4)
      NOUT = WMACH(11)
      NERR = WMACH(12)
C
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
C
      NAMED = .FALSE.
C
      INFORM = 0
      ITER = 0
C
      HEADER = .TRUE.
      PRNT = .TRUE.
C
      CONDMX = MAX(ONE/EPSPT5,HUNDRD)
C
C     Set the default values of the parameters.
C
      CALL E04MFW(N,NCLIN,TITLE)
C
      LLPTYP = LPROB
      NCTOTL = N + NCLIN
C
C     Set all parameters determined by the problem type.
C
      IF (LLPTYP.EQ.1) THEN
         PRBTYP = 'FP'
         CSET = .FALSE.
      ELSE IF (LLPTYP.EQ.2) THEN
         PRBTYP = 'LP'
         CSET = .TRUE.
      ELSE
         PRBTYP = 'illegal'
         MSG = 'noprob'
         GO TO 60
      END IF
C
C     Assign the dimensions of arrays in the parameter list of E04MFZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C     If a linear program is being solved and the matrix of general
C     constraints has fewer rows than columns, i.e.,  nclin .lt. n,
C     a non-zero value is known for MINFXD.  Note that in this case,
C     VERTEX must be set  .true..
C
      VERTEX = NCLIN .LT. N
C
      MINFXD = N - MXFREE
      MINACT = MXFREE - MAXNZ
C
      LDT = MAX(MAXNZ,MAXACT)
      NCOLT = MXFREE
      IF (NCLIN.EQ.0) THEN
         LDQ = 1
      ELSE
         LDQ = MAX(1,MXFREE)
      END IF
      LDR = LDT
C
      NCNLN = 0
      LENNAM = 1
      LDH = 1
      MM = 0
C
C     ==================================================================
C     Cold start:  Only  x  is provided.
C     Warm start:  Initial working set is specified in  ISTATE.
C     Hot  start:  The work arrays  IW  and  W  are assumed to have been
C                  initialized during a previous run.
C                  The first three elements of  IW  contain details
C                  on the dimension of the initial working set.
C     ==================================================================
      IF (LCRASH.EQ.0) THEN
         START = 'COLD'
      ELSE IF (LCRASH.EQ.1) THEN
         START = 'WARM'
      ELSE IF (LCRASH.EQ.2) THEN
         START = 'HOT '
      END IF
C
      COLD = LCRASH .EQ. 0
      WARM = LCRASH .EQ. 1
      HOT = LCRASH .EQ. 2
C
C     Allocate remaining work arrays.
C
      LITOTL = 3
      LWTOTL = 0
      CALL E04MFV(CSET,N,NCLIN,LITOTL,LWTOTL)
C
C     Check input parameters and storage limits.
C
      CALL E04MFP(NERROR,MSGLVL,START,LENIW,LENW,LITOTL,LWTOTL,N,NCLIN,
     *            NCNLN,ISTATE,NAMED,NAMES,BIGBND,BL,BU,X,LDA,LDH,MM,
     *            NERR,LLPTYP,IFAIL)
C
      IF (NERROR.GT.0) THEN
         MSG = 'errors'
         GO TO 60
      END IF
C
      LKACTV = LOCLC(1)
      LKX = LOCLC(2)
C
      LFEATU = LOCLC(3)
      LANORM = LOCLC(4)
      LD = LOCLC(7)
      LGQ = LOCLC(8)
      LCQ = LOCLC(9)
      LRLAM = LOCLC(10)
      LR = LOCLC(11)
      LT = LOCLC(12)
      LQ = LOCLC(13)
      LWTINF = LOCLC(14)
      LWRK = LOCLC(15)
C
C     ------------------------------------------------------------------
C     Define the initial feasibility tolerances in CLAMDA.
C     ------------------------------------------------------------------
      IF (TOLFEA.GT.ZERO) CALL F06FBF(N+NCLIN,(TOLFEA),W(LFEATU),1)
C
      CALL E04MFR('Initialize anti-cycling variables',MSGLVL,N,NCLIN,
     *            NMOVED,ITER,NUMINF,ISTATE,BIGBND,AX,BL,BU,CLAMDA,
     *            W(LFEATU),X)
C
      IF (COLD .OR. WARM) THEN
C        ---------------------------------------------------------------
C        Cold or warm start.  Just about everything must be initialized.
C        The only exception is ISTATE during a warm start.
C        ---------------------------------------------------------------
         IANRMJ = LANORM
         DO 20 J = 1, NCLIN
            W(IANRMJ) = DNRM2(N,A(J,1),LDA)
            IANRMJ = IANRMJ + 1
   20    CONTINUE
         IF (NCLIN.GT.0) CALL F06FLF(NCLIN,W(LANORM),1,ASIZE,AMIN)
C
         CALL F06FLF(NCTOTL,W(LFEATU),1,FEAMAX,FEAMIN)
         CALL DCOPY(NCTOTL,W(LFEATU),1,W(LWTINF),1)
         CALL DSCAL(NCTOTL,(ONE/FEAMIN),W(LWTINF),1)
C
C        ---------------------------------------------------------------
C        Define the initial working set.
C               NFREE ,  NACTIV,  KACTIV, KX,
C               ISTATE (if START  = 'COLD')
C               NARTIF (if VERTEX = 'TRUE')
C        ---------------------------------------------------------------
         CALL E04MFT(START,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,
     *               LDA,ISTATE,IW(LKACTV),IW(LKX),BIGBND,TOLACT,A,AX,
     *               BL,BU,CLAMDA,X,W(LGQ),W(LWRK))
C
C        ---------------------------------------------------------------
C        Compute the TQ factorization of the working set matrix.
C        ---------------------------------------------------------------
         UNITQ = .TRUE.
         NZ = NFREE
C
         IF (NACTIV.GT.0) THEN
            IT = NACTIV + 1
            NACT1 = NACTIV
            NACTIV = 0
            NGQ = 0
C
            CALL E04NFQ(UNITQ,VERTEX,1,NACT1,IT,NACTIV,NARTIF,NZ,NFREE,
     *                  NREJTD,NGQ,N,LDQ,LDA,LDT,ISTATE,IW(LKACTV),
     *                  IW(LKX),CONDMX,A,W(LT),W(LGQ),W(LQ),W(LWRK),
     *                  W(LD),W(LRLAM),MSGLVL)
         END IF
      ELSE IF (HOT) THEN
C        ---------------------------------------------------------------
C        Arrays  IW  and  W  have been defined in a previous run.
C        The first three elements of  IW  are  UNITQ,  NFREE and NACTIV.
C        ---------------------------------------------------------------
         UNITQ = IW(1) .EQ. 1
         NFREE = IW(2)
         NACTIV = IW(3)
C
         NZ = NFREE - NACTIV
      END IF
C
      IF (CSET) THEN
C
C        Install the transformed linear term in CQ.
C
         CALL DCOPY(N,CVEC,1,W(LCQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LCQ),W(LQ),W(LWRK)
     *               )
      END IF
C
      RSET = .FALSE.
      ITMAX = ITMAX2
      JINF = 0
C
C     +    Take your pick when minimizing the sum of infeasibilities:
C     +    NRZ    =  NZ  implies steepest-descent in the two-norm.
C     +    NRZ    =  0   implies steepest-descent in the infinity norm.
      NRZ = 0
C
C     ==================================================================
C     repeat               (until working set residuals are acceptable)
C     ---------------------------------------------------------------
C     Move x onto the constraints in the working set.
C     ---------------------------------------------------------------
   40 CALL E04MFJ(ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NZ,N,LDQ,LDA,LDT,
     *            ISTATE,IW(LKACTV),IW(LKX),JMAX,ERRMAX,XNORM,A,AX,BL,
     *            BU,W(LFEATU),W(LT),X,W(LQ),W(LD),W(LWRK))
C
      IF (ROWERR) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99981)
            CALL X04BAF(IPRINT,REC(1))
         END IF
         MSG = 'infeas'
         NUMINF = 1
         OBJ = ERRMAX
         GO TO 60
      END IF
C
      CALL E04MFZ(PRBTYP,MSG,CSET,NAMED,NAMES,RSET,UNITQ,ITER,ITMAX,
     *            JINF,NVIOL,N,NCLIN,LDA,NACTIV,NFREE,NRZ,NZ,ISTATE,
     *            IW(LKACTV),IW(LKX),E04MFU,OBJ,NUMINF,XNORM,A,AX,BL,BU,
     *            CVEC,CLAMDA,W(LFEATU),X,IW,W)
C
      FOUND = MSG .EQ. 'feasbl' .OR. MSG .EQ. 'optiml' .OR. MSG .EQ.
     *        'weak  ' .OR. MSG .EQ. 'unbndd' .OR. MSG .EQ. 'infeas'
      HALTED = MSG .EQ. 'itnlim'
C
      IF (FOUND) THEN
         CALL E04MFR('Optimal',MSGLVL,N,NCLIN,NMOVED,ITER,NUMINF,ISTATE,
     *               BIGBND,AX,BL,BU,CLAMDA,W(LFEATU),X)
      END IF
C
      DONE = FOUND .AND. NVIOL .EQ. 0 .AND. NMOVED .EQ. 0
C
C     until      done  .or.  halted
      IF ( .NOT. (DONE .OR. HALTED)) GO TO 40
C     ===========================================================
C     Set   CLAMDA.  Print the full solution.
C     Clean up.  Save values for a subsequent hot start.
C     ------------------------------------------------------------------
      CALL E04MFK(MSGLVL,NFREE,LDA,N,NCLIN,NCNLN,NCTOTL,BIGBND,NAMED,
     *            NAMES,LENNAM,NACTIV,ISTATE,IW(LKACTV),IW(LKX),A,BL,BU,
     *            X,CLAMDA,W(LRLAM),X)
C
      IW(1) = 0
      IF (UNITQ) IW(1) = 1
      IW(2) = NFREE
      IW(3) = NACTIV
C
C     ==================================================================
C     Print messages if required.
C     Recover the optional parameters set by the user.
C     ==================================================================
   60 IF (MSG.EQ.'optiml') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99997) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'feasbl') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'weak  ') THEN
         INFORM = 1
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99996) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99995)
     *       PRBTYP
C
      ELSE IF (MSG.EQ.'unbndd') THEN
         INFORM = 2
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99994) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99993)
     *       PRBTYP
C
      ELSE IF (MSG.EQ.'infeas') THEN
         INFORM = 3
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99992)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99991)
C
      ELSE IF (MSG.EQ.'itnlim') THEN
         INFORM = 4
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99990)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99989)
C
      ELSE IF (MSG.EQ.'errors') THEN
         INFORM = 6
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99988) NERROR
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99987)
     *       NERROR
C
      ELSE IF (MSG.EQ.'noprob') THEN
         INFORM = 7
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99986)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99985)
      END IF
C
      IF (MSGLVL.GT.0) THEN
C
         IF (INFORM.LT.6) THEN
            IF (NUMINF.EQ.0) THEN
               IF (PRBTYP.NE.'FP') THEN
                  WRITE (REC,FMT=99984) PRBTYP, OBJ
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE IF (INFORM.EQ.3) THEN
               WRITE (REC,FMT=99983) OBJ
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               WRITE (REC,FMT=99982) OBJ
               CALL X04BAY(IPRINT,2,REC)
            END IF
         END IF
      END IF
C
      IF (INFORM.LT.6) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
            CALL X04BAY(IPRINT,2,REC)
         END IF
      END IF
C
      CALL F06DFF(MXPARM,IPSVLC,1,IPRMLC,1)
      CALL DCOPY(MXPARM,RPSVLC,1,RPRMLC,1)
C
      IF ((INFORM.GE.1 .AND. INFORM.LE.7)
     *    .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) CALL X04BAY(NERR,2,
     *    ERRREC)
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
      RETURN
C
C
C
C     End of  E04MFF.  (LPOPT)
C
99999 FORMAT (/' Exit from ',A2,' problem after ',I5,' iterations.')
99998 FORMAT (/' Exit E04MFF - Feasible point found.     ')
99997 FORMAT (/' Exit E04MFF - Optimal ',A2,' solution.')
99996 FORMAT (/' Exit E04MFF - Weak ',A2,' solution.')
99995 FORMAT (/' ** Weak ',A2,' solution.')
99994 FORMAT (/' Exit E04MFF - ',A2,' solution is unbounded.')
99993 FORMAT (/' ** ',A2,' solution is unbounded.')
99992 FORMAT (/' Exit E04MFF - No feasible point for the linear constr',
     *       'aints.')
99991 FORMAT (/' ** No feasible point for the linear constraints.')
99990 FORMAT (/' Exit E04MFF - Too many iterations.')
99989 FORMAT (/' ** Too many iterations.')
99988 FORMAT (/' Exit E04MFF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99987 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
99986 FORMAT (/' Exit E04MFF - Problem type not recognized.  Problem a',
     *       'bandoned.')
99985 FORMAT (/' ** Problem type not recognized.  Problem abandoned.')
99984 FORMAT (/' Final ',A2,' objective value =',G16.7)
99983 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99982 FORMAT (/' Final sum of infeasibilities =',G16.7)
99981 FORMAT (' XXX  Cannot satisfy the working set constraints to the',
     *       ' accuracy requested.')
      END

      SUBROUTINE nlpsol (N,NCLIN,NCNLN,LDA,LDCJU,LDR,A,BL,BU,CONFUN,
     *                  OBJFUN,ITER,ISTATE,C,CJACU,CLAMDA,OBJF,GRADU,R,
     *                  X,IW,LENIW,W,LENW,IUSER,USER,IFAIL,jprint)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-951 (NOV 1991).
C     MARK 16 REVISED. IER-1077 (JUL 1993).
C
C     ==================================================================
C     E04UCF   solves the nonlinear program
C
C            minimize                   F(x)
C
C                                    (    x  )
C            subject to    bl  .le.  (  A*x  )  .le.  bu
C                                    (  c(x) )
C
C     where  F(x)  is a smooth scalar function,  A  is a constant matrix
C     and  c(x)  is a vector of smooth nonlinear functions.  The
C     feasible region is defined by a mixture of linear and nonlinear
C     equality or inequality constraints on  x.
C
C     The dimensions of the problem are...
C
C     N        the number of variables (dimension of  x),
C
C     NCLIN    the number of linear constraints (rows of the matrix  A),
C
C     NCNLN    the number of nonlinear constraints (dimension of  c(x)),
C
C
C     E04UCF   uses a sequential quadratic programming algorithm, with a
C     positive-definite quasi-Newton approximation to the transformed
C     Hessian  Q'HQ  of the Lagrangian function (which will be stored in
C     the array  R).
C
C
C     Complete documentation for  E04UCF  is contained in Report
C     SOL 86-2, Users guide for E04UCF (Version 4.0), by P.E. Gill,
C     W. Murray, M.A. Saunders and M.H. Wright, Department of Operations
C     Research,  Stanford University, Stanford, California 94305.
C
C     Systems Optimization Laboratory, Stanford University.
C     Version 1.1,  April     12, 1983. The less said about this one...
C     Version 2.0,  April     30, 1984.
C     Version 3.0,  March     20, 1985. First Fortran 77 version
C     Version 3.2,  August    20, 1985.
C     Version 4.0,  April     16, 1986. First version with differences
C     Version 4.01, June      30, 1986. Level 2 BLAS + F77 linesearch
C     Version 4.02, August     5, 1986. Reset SSBFGS. One call to E04XAY
C     Version 4.03, June      14, 1987. Step limit
C     Version 4.04, June      28, 1989. Vectorizable BLAS
C     Version 4.05, November  28, 1989. Load and save files added
C                                       (but not to NAG version!)
C     Version 4.06, November   5, 1991. E04UCJ and E04UCK updated
C                   October   17, 1992. Summary/print file option.
C
C     Copyright  1983  Stanford University.
C
C     This material may be reproduced by or for the U.S. Government
C     pursuant to the copyright license under DAR Clause 7-104.9(a)
C     (1979 Mar).
C
C     This material is based upon work partially supported by the
C     National Science Foundation under Grants MCS-7926009 and
C     ECS-8312142; the Department of Energy Contract AM03-76SF00326,
C     PA No. DE-AT03-76ER72018; the Army Research Office Contract
C     DAA29-84-K-0156; and the Office of Naval Research Grant
C     N00014-75-C-0267.
C     ==================================================================
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04UCF')
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, POINT3, POINT8
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1,POINT8=0.8D+0)
      DOUBLE PRECISION  POINT9, ONE
      PARAMETER         (POINT9=0.9D+0,ONE=1.0D+0)
      DOUBLE PRECISION  GROWTH
      PARAMETER         (GROWTH=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJF
      INTEGER           IFAIL, ITER, LDA, LDCJU, LDR, LENIW, LENW, N,
     *                  NCLIN, NCNLN, jprint
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN),
     *                  C(*), CJACU(LDCJU,*), CLAMDA(N+NCLIN+NCNLN),
     *                  GRADU(N), R(LDR,*), USER(*), W(LENW), X(N)
      INTEGER           ISTATE(N+NCLIN+NCNLN), IUSER(*), IW(LENIW)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DRMAX, DRMIN, DTMAX, DTMIN, DXLIM, EPSPT3,
     *                  EPSPT5, EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL,
     *                  HCNDBD, RCNDBD, RFROBN, RHODMP, RHOMAX, RHONRM,
     *                  SCALE, TOLACT, TOLFEA, TOLRNK
      INTEGER           IDBGLS, IDBGNP, IPRINT, IPRNT, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, KSAVE, LCRASH, LDBGLS, LDBGNP, LDQ, LDT,
     *                  LENNAM, LFDSET, LFORMH, LINES1, LINES2, LPROB,
     *                  LVERFY, LVLDER, LVLDIF, LVRFYC, MSGLS, MSGNP,
     *                  NACTIV, NCDIFF, NCOLT, NFDIFF, NFREE, NLNF,
     *                  NLNJ, NLNX, NLOAD, NN, NNCLIN, NNCNLN, NOUT,
     *                  NPROB, NSAVE, NZ
      LOGICAL           CMDBG, INCRUN, LSDBG, NPDBG, UNITQ
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), ILSDBG(LDBG), INPDBG(LDBG),
     *                  IPADLS(18), IPADNP(12), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM), JVERFY(4), LOCLS(LENLS),
     *                  LOCNP(LENNP)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, COND, CONDMX, CTX, EPSMCH, ERRMAX, FDCHK,
     *                  FDNORM, FEAMAX, FEAMIN, OBJ, ROOTN, RTEPS, SSQ1,
     *                  SUMINF, XNORM
      INTEGER           I, IANRMJ, IDBG, IDBGSV, IKX, INFO, INFORM,
     *                  ITMXSV, ITNS, J, JINF, JMAX, LANORM, LAQP, LAX,
     *                  LCJAC, LCJDX, LCLAM, LCMUL, LDAQP, LDCJ, LDFJU,
     *                  LDX, LFEATL, LGQ, LGRAD, LHCTRL, LHFRWD, LIPERM,
     *                  LITOTL, LKACTV, LKX, LNEEDC, LQ, LRES, LRES0,
     *                  LRHO, LRLAM, LT, LV, LVJ, LWRK1, LWRK2, LWRK3,
     *                  LWTINF, LWTOTL, M, MAXACT, MAXNZ, MINACT,
     *                  MINFXD, MJRDBG, MNRDBG, MSGQP, MXFREE, NACT1,
     *                  NARTIF, NCTOTL, NERR, NERROR, NFUN, NGQ, NGRAD,
     *                  NLPERR, NMAJOR, NMINOR, NPLIN, NRANK, NREJTD,
     *                  NRES, NSTATE, NUMINF, NZ1
      LOGICAL           COLD, LINOBJ, NAMED, NEEDFD, OVERFL, ROWERR,
     *                  VERTEX
      CHARACTER*11      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF, F06RJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, F06BLF, F06RJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBW, E04NBZ, E04NCH, E04NCU,
     *                  E04NCX, E04NCY, E04NCZ, E04UCP, E04UCS, E04UCX,
     *                  E04UCY, E04UCZ, E04UDR, F01QCF, F06DFF, F06FBF,
     *                  F06FLF, F06QFF, F06QHF, X02ZAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IDBGLS, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LPROB, MSGLS, NN,
     *                  NNCLIN, NPROB, IPADLS
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /FE04NC/NACTIV, NFREE, NZ, UNITQ
      COMMON            /FE04UC/INPDBG, NPDBG
      COMMON            /GE04UC/IPSVNP, IDBGNP, ITMXNP, JVRFY1, JVRFY2,
     *                  JVRFY3, JVRFY4, LDBGNP, LFORMH, LVLDER, LVERFY,
     *                  MSGNP, NLNF, NLNJ, NLNX, NNCNLN, NSAVE, NLOAD,
     *                  KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IDBGLS), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),IDBGNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IDBGNP,IDBG), (ITMXNP,NMAJOR), (ITMAX2,NMINOR)
      EQUIVALENCE       (LDBGLS,MNRDBG), (LDBGNP,MJRDBG), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/,
     *                  /AX02ZA/, /FE04NC/
C     .. Data statements ..
      DATA              TITLE/' *** E04UCF'/
C     .. Executable Statements ..
C
C     Set the machine-dependent constants.
C
      CALL X02ZAZ
C
      EPSMCH = WMACH(3)
      RTEPS = WMACH(4)
      NOUT = WMACH(11)
      NERR = WMACH(12)
C
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
C
      RHOMAX = ONE/EPSMCH
      ROOTN = SQRT(DBLE(N))
C
C     Default names will be provided for variables during printing.
C
      NAMED = .FALSE.
      INFORM = 0
C
C     Set the default values for the parameters.
C
      CALL E04UCX(N,NCLIN,NCNLN,TITLE)
C
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR.
     *         (LVLDER.EQ.1 .AND. NCNLN.GT.0)
      COLD = LCRASH .EQ. 0
      LVLDIF = 0
      IF (NEEDFD) LVLDIF = 1
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
C
C     Assign the dimensions of arrays in the parameter list of E04UCZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C
      MINACT = 0
      MINFXD = 0
C
      MXFREE = N - MINFXD
      MAXACT = MAX(1,MIN(N,NCLIN))
      MAXNZ = N - (MINFXD+MINACT)
C
      IF (NCLIN+NCNLN.EQ.0) THEN
         LDQ = 1
         LDT = 1
         NCOLT = 1
      ELSE
         LDQ = MAX(1,MXFREE)
         LDT = MAX(MAXNZ,MAXACT)
         NCOLT = MXFREE
      END IF
C
      LENNAM = 1
      M = 1
      LDFJU = 2
C
      LDAQP = MAX(NCLIN+NCNLN,1)
      IF (NCNLN.EQ.0 .AND. NCLIN.GT.0) LDAQP = LDA
C
C     E04UCP  defines the arrays that contain the locations of various
C     work arrays within  W  and  IW.
C
      LITOTL = 0
      LWTOTL = 0
      CALL E04UCP(N,NCLIN,NCNLN,NCTOTL,LITOTL,LWTOTL)
C
C     Allocate certain addresses that are not allocated in E04UCP.
C
      LAX = LWTOTL + 1
      LWTOTL = LAX + NCLIN - 1
      LAX = MIN(LAX,LWTOTL)
C
C     Check input parameters and storage limits.
C
      CALL E04NBZ(NERROR,MSGNP,LCRASH,LENIW,LENW,LITOTL,LWTOTL,N,NCLIN,
     *            NCNLN,ISTATE,IW,NAMED,NAMES,BIGBND,BL,BU,X,M,LDA,LDR,
     *            LDCJU,LDFJU,NERR,IFAIL)
C
      IF (NERROR.GT.0) THEN
         INFORM = 9
         GO TO 80
      END IF
C
      LKACTV = LOCLS(1)
      LANORM = LOCLS(2)
      LCJDX = LOCLS(3)
      LRES = LOCLS(5)
      LRES0 = LOCLS(6)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LQ = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1 = LOCLS(14)
C
      LKX = LOCNP(1)
      LIPERM = LOCNP(2)
      LAQP = LOCNP(3)
      LDX = LOCNP(7)
      LFEATL = LOCNP(10)
      LWRK2 = LOCNP(12)
C
      LCMUL = LOCNP(16)
      LRHO = LOCNP(20)
      LWRK3 = LOCNP(21)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
      LCJAC = LOCNP(27)
      LGRAD = LOCNP(28)
C
      LDCJ = MAX(NCNLN,1)
C
      TOLRNK = ZERO
      RCNDBD = SQRT(HCNDBD)
C
C     ==================================================================
C     If a unit number for a load file has been set, read initial values
C     from an old run.  These values override existing settings.
C     ==================================================================
C     IF (NLOAD.GT.0) THEN
C        CALL NPGETR(NERROR,UNITQ,N,NCLIN,NCNLN,LDR,LDQ,NFREE0,ITER,
C     *               ISTATE,IW(LKX),W(LHFRWD),W(LHCTRL),W(LCMUL),R,
C     *               W(LRHO),X,W(LQ))
C
C        IF (NERROR.GT.0) THEN
C           INFORM = 9
C           GO TO 80
C        END IF
C     END IF
C
C     ==================================================================
C     Load the arrays of feasibility tolerances.
C     ==================================================================
      IF (TOLFEA.GT.ZERO) CALL F06FBF(NPLIN,TOLFEA,W(LFEATL),1)
C
      IF (NCNLN.GT.0 .AND. CTOL.GT.ZERO) CALL F06FBF(NCNLN,CTOL,
     *    W(LFEATL+NPLIN),1)
C
      IF (LFDSET.EQ.0) THEN
         FDCHK = SQRT(EPSRF)
      ELSE IF (LFDSET.EQ.1) THEN
         FDCHK = FDINT
      ELSE
         FDCHK = W(LHFRWD)
      END IF
C
      NFUN = 0
      NGRAD = 0
      NSTATE = 1
C
C     ------------------------------------------------------------------
C     If required,  compute the problem functions.
C     If the constraints are nonlinear,  the first call of confun
C     sets up any constant elements in the Jacobian matrix.  A copy of
C     the Jacobian (with constant elements set) is placed in  CJACU.
C     ------------------------------------------------------------------
      IF (LVERFY.GE.10) THEN
         XNORM = DNRM2(N,X,1)
         LVRFYC = LVERFY - 10
C
         CALL E04UCY(INFO,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,LDCJU,N,
     *               NCNLN,CONFUN,OBJFUN,IW(LNEEDC),BIGBND,EPSRF,CDINT,
     *               FDINT,FDCHK,FDNORM,OBJF,XNORM,BL,BU,C,W(LWRK3),
     *               W(LCJAC),CJACU,W(LCJDX),W(LDX),W(LGRAD),GRADU,
     *               W(LHFRWD),W(LHCTRL),X,W(LWRK1),W(LWRK2),W,LENW,
     *               IUSER,USER)
C
         IF (INFO.NE.0) THEN
            IF (INFO.GT.0) INFORM = 7
            IF (INFO.LT.0) INFORM = INFO
            GO TO 80
         END IF
         NSTATE = 0
      END IF
C
      CALL F06DFF(LDBG,ILSDBG,1,ICMDBG,1)
C
      IF (NCLIN.GT.0) THEN
         IANRMJ = LANORM
         DO 20 J = 1, NCLIN
            W(IANRMJ) = DNRM2(N,A(J,1),LDA)
            IANRMJ = IANRMJ + 1
   20    CONTINUE
         CALL F06FLF(NCLIN,W(LANORM),1,ASIZE,AMIN)
      END IF
C
      CALL F06FLF(NPLIN,W(LFEATL),1,FEAMAX,FEAMIN)
      CALL DCOPY(NPLIN,W(LFEATL),1,W(LWTINF),1)
      CALL DSCAL(NPLIN,(ONE/FEAMIN),W(LWTINF),1)
C
C     ==================================================================
C     The input values of x and (optionally)  ISTATE are used by
C     E04NCU  to define an initial working set.
C     ==================================================================
      VERTEX = .FALSE.
      CALL E04NCU(COLD,VERTEX,NCLIN,NPLIN,NACTIV,NARTIF,NFREE,N,LDA,
     *            ISTATE,IW(LKACTV),BIGBND,TOLACT,A,W(LAX),BL,BU,X,
     *            W(LWRK1),W(LWRK2))
C
      NRES = 0
      NGQ = 0
      CONDMX = ONE/EPSPT5
C
      IF (LCRASH.LE.1) THEN
C        ===============================================================
C        Cold or warm start. The upper-triangular matrix R is the factor
C        of an approximate Lagrangian Hessian.
C        ===============================================================
         UNITQ = .TRUE.
         ITER = 0
C
         IKX = LKX
         DO 40 I = 1, N
            IW(IKX) = I
            IKX = IKX + 1
   40    CONTINUE
C
         IF (COLD) THEN
            CALL F06QHF('Upper-triangular',N,N,ZERO,ONE,R,LDR)
            RFROBN = ROOTN
C
            NRANK = 0
            IF (NCNLN.GT.0) CALL F06FBF(NCNLN,(ZERO),W(LCMUL),1)
         ELSE
C
C           R will be updated while finding a feasible x.
C
            NRANK = NLNX
            CALL F06FBF(NLNX,(ZERO),W(LRES0),1)
            IF (NCNLN.GT.0) CALL DCOPY(NCNLN,CLAMDA(NPLIN+1),1,W(LCMUL),
     *                                 1)
C
         END IF
C
         INCRUN = .TRUE.
         RHONRM = ZERO
         RHODMP = ONE
         SCALE = ONE
         CALL F06FBF(NCNLN,(ZERO),W(LRHO),1)
C
C        ---------------------------------------------------------------
C        Re-order KX so that the free variables come first.
C        If a warm start is required, NRANK will be nonzero and the
C        factor R will be updated.
C        ---------------------------------------------------------------
         CALL E04NCX(UNITQ,INFORM,NZ,NFREE,NRANK,NRES,NGQ,N,LDQ,LDA,LDR,
     *               LDT,ISTATE,IW(LKX),CONDMX,A,R,W(LT),W(LRES0),W(LGQ)
     *               ,W(LQ),W(LWRK1),W(LWRK2),W(LRLAM),MSGNP)
C
      ELSE
C        ===============================================================
C        Hot start.
C        Stop if the computed and input values of NFREE don't match.
C        ===============================================================
C        IF (NFREE0.NE.NFREE) THEN
C           NERROR = 1
C           INFORM = 9
C           GO TO 80
C        END IF
      END IF
C
C     ------------------------------------------------------------------
C     Factorize the linear constraints in the initial working set.
C     ------------------------------------------------------------------
      IF (NACTIV.GT.0) THEN
         NACT1 = NACTIV
         NACTIV = 0
C
         CALL E04NCY(UNITQ,VERTEX,INFORM,1,NACT1,NACTIV,NARTIF,NZ,NFREE,
     *               NRANK,NREJTD,NRES,NGQ,N,LDQ,LDA,LDR,LDT,ISTATE,
     *               IW(LKACTV),IW(LKX),CONDMX,A,R,W(LT),W(LRES0),W(LGQ)
     *               ,W(LQ),W(LWRK1),W(LWRK2),W(LRLAM),MSGNP)
      END IF
C
      IF (LCRASH.LE.1) THEN
C        ===============================================================
C        Cold or warm start.  Move  x  on to the linear constraints and
C        find a feasible point.
C        ===============================================================
         SSQ1 = ZERO
         LINOBJ = .FALSE.
         CALL E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NRANK,NZ,N,
     *               NPLIN,LDQ,LDA,LDR,LDT,ISTATE,IW(LKACTV),IW(LKX),
     *               JMAX,ERRMAX,CTX,XNORM,A,W(LAX),BL,BU,W(LGQ),W(LRES)
     *               ,W(LRES0),W(LFEATL),R,W(LT),X,W(LQ),W(LWRK1),
     *               W(LWRK2))
C
C        ---------------------------------------------------------------
C        Call  E04NCZ  to find a feasible  x.
C        ---------------------------------------------------------------
C        Use  WORK2  as the multiplier vector.
C
         JINF = 0
         LCLAM = LWRK2
C
         IDBGSV = IDBG
         IF (IDBG.GT.0) THEN
            IDBG = NMINOR + 1
         END IF
C
         ITMXSV = ITMAX1
         ITMAX1 = NMINOR
C
         CALL E04NCZ('FP problem',NAMED,NAMES,LINOBJ,UNITQ,NLPERR,ITNS,
     *               JINF,NCLIN,NPLIN,NACTIV,NFREE,NRANK,NZ,NZ1,N,LDA,
     *               LDR,ISTATE,IW(LKACTV),IW(LKX),CTX,OBJ,SSQ1,SUMINF,
     *               NUMINF,XNORM,BL,BU,A,W(LCLAM),W(LAX),W(LFEATL),R,X,
     *               W)
C
         ITMAX1 = ITMXSV
C
         IF (NLPERR.GT.0) THEN
            INFORM = 2
            GO TO 80
         ELSE IF (MSGQP.GT.0) THEN
            WRITE (REC,FMT=99987)
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
         IDBG = IDBGSV
         CALL F06DFF(LDBG,INPDBG,1,ICMDBG,1)
      ELSE
C        ---------------------------------------------------------------
C        Hot start.
C        The point  x  is preassigned.  Compute the 2-norm of  x.
C        Initialize  Ax  for the linear constraints.
C        ---------------------------------------------------------------
C        NRANK = NLNX
C        XNORM = DNRM2(N,X,1)
C        IF (NCLIN.GT.0) CALL DGEMV('N',NCLIN,N,ONE,A,LDA,X,1,ZERO,
C        *                            W(LAX),1)
      END IF
C
      IF (LCRASH.GT.0) THEN
C
C        Check for a bad R.
C
         RFROBN = F06RJF('Frobenius norm','Upper','Non-unit diagonal',N,
     *            N,R,LDR,W)
         CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
         COND = F06BLF(DRMAX,DRMIN,OVERFL)
C
         IF (COND.GT.RCNDBD .OR. RFROBN.GT.ROOTN*GROWTH*DRMAX) THEN
C           ------------------------------------------------------------
C           Refactorize the Hessian and bound the condition estimator.
C           ------------------------------------------------------------
            IF (MSGNP.GT.0) THEN
               WRITE (REC,FMT=99986)
               CALL X04BAF(IPRINT,REC(1))
            END IF
            CALL E04UDR(UNITQ,N,NFREE,NZ,LDQ,LDR,IW(LIPERM),IW(LKX),
     *                  W(LGQ),R,W(LQ),W(LWRK1),W(LRES0))
         END IF
      END IF
C
C     ==================================================================
C     Now we can check the gradients at a feasible x.
C     ==================================================================
      LVRFYC = LVERFY
      IF (LVERFY.GE.10) LVRFYC = -1
C
      CALL E04UCY(INFO,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,LDCJU,N,
     *            NCNLN,CONFUN,OBJFUN,IW(LNEEDC),BIGBND,EPSRF,CDINT,
     *            FDINT,FDCHK,FDNORM,OBJF,XNORM,BL,BU,C,W(LWRK3),
     *            W(LCJAC),CJACU,W(LCJDX),W(LDX),W(LGRAD),GRADU,
     *            W(LHFRWD),W(LHCTRL),X,W(LWRK1),W(LWRK2),W,LENW,IUSER,
     *            USER)
C
      IF (INFO.NE.0) THEN
         IF (INFO.GT.0) INFORM = 7
         IF (INFO.LT.0) INFORM = INFO
         GO TO 80
      END IF
C
      CALL DCOPY(N,W(LGRAD),1,W(LGQ),1)
      CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LGQ),W(LQ),W(LWRK1))

      iuser(2) = 1
C     ==================================================================
C     Solve the problem.
C     ==================================================================
      IF (NCNLN.EQ.0) THEN
C        ---------------------------------------------------------------
C        The problem has only linear constraints and bounds.
C        ---------------------------------------------------------------
         CALL E04UCZ(NAMED,NAMES,UNITQ,INFORM,ITER,N,NCLIN,NCNLN,NCTOTL,
     *               NACTIV,NFREE,NZ,LDCJ,LDCJU,LDAQP,LDR,NFUN,NGRAD,
     *               ISTATE,IW(LKACTV),IW(LKX),OBJF,FDNORM,XNORM,CONFUN,
     *               OBJFUN,A,W(LAX),BL,BU,C,W(LCJAC),CJACU,CLAMDA,
     *               W(LFEATL),W(LGRAD),GRADU,R,X,IW,W,LENW,IUSER,USER)
      ELSE
C        ---------------------------------------------------------------
C        The problem has some nonlinear constraints.
C        ---------------------------------------------------------------
         IF (NCLIN.GT.0) CALL F06QFF('General',NCLIN,N,A,LDA,W(LAQP),
     *                               LDAQP)
C
C        Try and add some nonlinear constraint indices to KACTIV.
C
         CALL E04UCS(COLD,N,NCLIN,NCNLN,NCTOTL,NACTIV,NFREE,NZ,ISTATE,
     *               IW(LKACTV),BIGBND,TOLACT,BL,BU,C)
C
         CALL E04UCZ(NAMED,NAMES,UNITQ,INFORM,ITER,N,NCLIN,NCNLN,NCTOTL,
     *               NACTIV,NFREE,NZ,LDCJ,LDCJU,LDAQP,LDR,NFUN,NGRAD,
     *               ISTATE,IW(LKACTV),IW(LKX),OBJF,FDNORM,XNORM,CONFUN,
     *               OBJFUN,W(LAQP),W(LAX),BL,BU,C,W(LCJAC),CJACU,
     *               CLAMDA,W(LFEATL),W(LGRAD),GRADU,R,X,IW,W,LENW,
     *               IUSER,USER)
C
      END IF
C
C     ==================================================================
C     If a unit number for a save file has been set, save the details of
C     this run.
C     ==================================================================
C     IF (NSAVE.GT.0 .AND. KSAVE.GT.NMAJOR) THEN
C        CALL NPSAVR(UNITQ,N,NCLIN,NCNLN,LDR,LDQ,NFREE,NSAVE,ITER,
C     *               ISTATE,IW(LKX),W(LHFRWD),W(LHCTRL),W(LCMUL),R,
C     *               W(LRHO),X,W(LQ))
C     END IF
C
C     ------------------------------------------------------------------
C     If required, form the triangular factor of the Hessian.
C     ------------------------------------------------------------------
C     First,  form the square matrix  R  such that  H = R'R.
C     Compute the  QR  factorization of  R.
C
      IF (LFORMH.GT.0) THEN
         LV = LWRK2
         DO 60 J = 1, N
            IF (J.GT.1) CALL F06FBF(J-1,ZERO,W(LV),1)
C
            LVJ = LV + J - 1
            CALL DCOPY(N-J+1,R(J,J),LDR,W(LVJ),1)
            CALL E04NBW(3,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LV),W(LQ),
     *                  W(LWRK1))
            CALL DCOPY(N,W(LV),1,R(J,1),LDR)
   60    CONTINUE
C
         CALL F01QCF(N,N,R,LDR,W(LWRK1),INFO)
      END IF
C
C     ==================================================================
C     Print messages if required.
C     ==================================================================
   80 IF (MSGNP.GT.0) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99999)
         IF (INFORM.EQ.0) WRITE (REC,FMT=99998)
         IF (INFORM.EQ.1) WRITE (REC,FMT=99997)
         IF (INFORM.EQ.2) WRITE (REC,FMT=99996)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99995)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99994)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99993)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99992)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99991) NERROR
         CALL X04BAY(IPRINT,2,REC)
C
         IF (INFORM.GE.0 .AND. INFORM.LT.7) THEN
            IF (NLPERR.EQ.0) THEN
               WRITE (REC,FMT=99990) OBJF
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               IF (NLPERR.EQ.3) THEN
                  WRITE (REC,FMT=99989) SUMINF
                  CALL X04BAY(IPRINT,2,REC)
               ELSE
                  WRITE (REC,FMT=99988) SUMINF
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            END IF
         END IF
      END IF
C
C     Recover the optional parameters set by the user.
C
      CALL F06DFF(MXPARM,IPSVLS,1,IPRMLS,1)
      CALL DCOPY(MXPARM,RPSVLS,1,RPRMLS,1)
      CALL F06DFF(MXPARM,IPSVNP,1,IPRMNP,1)
      CALL DCOPY(MXPARM,RPSVNP,1,RPRMNP,1)
C
      IF (INFORM.LT.9) THEN
         IF (NCNLN.GT.0) CALL F06QFF('General',NCNLN,N,W(LCJAC),LDCJ,
     *                               CJACU,LDCJU)
         CALL DCOPY(N,W(LGRAD),1,GRADU,1)
      END IF
C
      IF (INFORM.NE.0 .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99985)
         IF (INFORM.EQ.1) WRITE (REC,FMT=99984)
         IF (INFORM.EQ.2) WRITE (REC,FMT=99983)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99982)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99981)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99980)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99979)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99978) NERROR
         CALL X04BAY(NERR,2,REC)
      END IF
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
      RETURN
C
C
C
C     End of  E04UCF. (NPSOL)
C
99999 FORMAT (/' Exit E04UCF - User requested termination.')
99998 FORMAT (/' Exit E04UCF - Optimal solution found.')
99997 FORMAT (/' Exit E04UCF - Optimal solution found, but requested a',
     *       'ccuracy not achieved.')
99996 FORMAT (/' Exit E04UCF - No feasible point for the linear constr',
     *       'aints.')
99995 FORMAT (/' Exit E04UCF - No feasible point for the nonlinear con',
     *       'straints.')
99994 FORMAT (/' Exit E04UCF - Too many major iterations.             ')
99993 FORMAT (/' Exit E04UCF - Current point cannot be improved upon. ')
99992 FORMAT (/' Exit E04UCF - Large errors found in the derivatives. ')
99991 FORMAT (/' Exit E04UCF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99990 FORMAT (/' Final objective value =',G16.7)
99989 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99988 FORMAT (/' Final sum of infeasibilities =',G16.7)
99987 FORMAT (/' The linear constraints are feasible.')
99986 FORMAT (' XXX  Bad initial Hessian,   R  refactorized.')
99985 FORMAT (/' ** User requested termination.')
99984 FORMAT (/' ** Optimal solution found, but requested accuracy not',
     *       ' achieved.')
99983 FORMAT (/' ** No feasible point for the linear constraints.')
99982 FORMAT (/' ** No feasible point for the nonlinear constraints.')
99981 FORMAT (/' ** Too many major iterations.             ')
99980 FORMAT (/' ** Current point cannot be improved upon. ')
99979 FORMAT (/' ** Large errors found in the derivatives. ')
99978 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
      END

      SUBROUTINE A00AAF
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Writes information about the particular implementation of the
C     NAG Library in use.
C
C     The output channel is given by a call to X04ABF.
C
C     .. Local Scalars ..
      INTEGER          I, NADV
C     .. Local Arrays ..
      CHARACTER*80     MSG(20)
C     .. External Subroutines ..
      EXTERNAL         A00AAZ, X04ABF, X04BAF
C     .. Executable Statements ..
      CALL A00AAZ(MSG)
      CALL X04ABF(0,NADV)
      DO 20 I = 1, 7
         CALL X04BAF(NADV,MSG(I))
   20 CONTINUE
      CALL X04BAF(NADV,MSG(20))
      RETURN
      END

      SUBROUTINE A00AAZ(MSG)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Returns information about the particular implementation of the
C     NAG Fortran Library in use.
C
C     **************************************************************
C
C     Implementors must insert the correct details for each
C     distinct implementation.
C
C     **************************************************************
C
C     .. Array Arguments ..
      CHARACTER*80      MSG(20)
C     .. Executable Statements ..
      MSG(1) = ' *** Start of NAG Library implementation details ***'
      MSG(2) = ' '
      MSG(3) = ' Implementation title: Sun Solaris'
      MSG(4) = '            Precision: Fortran Double Precision'
      MSG(5) = '         Product Code: FLSOL16D'
      MSG(6) = '                 Mark: 16A'
      MSG(7) = ' '
      MSG(8) = ' Created using:'
      MSG(9) = '     hardware -   Sun SPARCstation 1'
      MSG(10) = '     op. sys. -   Solaris 2.1'
      MSG(11) = '     compiler -   Sun f77 2.0'
      MSG(12) = '       using optimisation = -O'
      MSG(13) = '       BLAS information: NAG BLAS used'
      MSG(14) = ' '
      MSG(15) = ' Applicable to:'
      MSG(16) = '     hardware -   Sun Sparcstation'
      MSG(17) = '     op. sys. -   Solaris 2.x'
      MSG(18) = '     compiler -   Sun f77 2.0'
      MSG(19) = ' '
      MSG(20) = ' *** End of NAG Library implementation details ***'
      RETURN
      END

      SUBROUTINE E04MFP(NERROR,MSGLVL,START,LIWORK,LWORK,LITOTL,LWTOTL,
     *                  N,NCLIN,NCNLN,ISTATE,NAMED,NAMES,BIGBND,BL,BU,X,
     *                  LDA,LDH,MM,NERR,LPROB,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFP   checks the input data for E04MFF and E04NFF.
C
C     First version written by PEG,  15-Nov-90.
C     This version of E04MFP dated   21-Mar-93.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           IFAIL, LDA, LDH, LITOTL, LIWORK, LPROB, LWORK,
     *                  LWTOTL, MM, MSGLVL, N, NCLIN, NCNLN, NERR,
     *                  NERROR
      LOGICAL           NAMED
      CHARACTER*4       START
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN), X(N)
      INTEGER           ISTATE(N+NCLIN+NCNLN)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2
      INTEGER           IS, J, K, L
      LOGICAL           OK
      CHARACTER*4       CSTART
C     .. Local Arrays ..
      CHARACTER*5       ID(3)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Data statements ..
      DATA              ID(1), ID(2), ID(3)/'Varbl', 'L Con', 'N Con'/
C     .. Executable Statements ..
C
      NERROR = 0
C
C     ------------------------------------------------------------------
C     Check  N.
C     ------------------------------------------------------------------
      IF (N.LE.0) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99997) N
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  NCLIN and NCNLN.
C     ------------------------------------------------------------------
      IF (NCLIN.LT.0 .OR. NCNLN.LT.0) THEN
         IF (NCLIN.LT.0) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99996) NCLIN
               CALL X04BAY(NERR,3,REC)
            END IF
         END IF
C
         IF (NCNLN.LT.0) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99995) NCNLN
               CALL X04BAY(NERR,3,REC)
            END IF
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDA.
C     ------------------------------------------------------------------
      IF (LDA.LT.MAX(1,NCLIN)) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99994) LDA, NCLIN
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDH.
C     ------------------------------------------------------------------
      IF ((LPROB.EQ.1) .OR. (LPROB.EQ.2)) THEN
C        Problem Type = FP or LP
         IF (LDH.LT.1) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99993) LDH
               CALL X04BAY(NERR,3,REC)
            END IF
         END IF
      ELSE
C        Problem Type = QP1, QP2, QP3 or QP4
         IF (LDH.LT.MM) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99992) LDH, N, MM
               CALL X04BAY(NERR,4,REC)
            END IF
         END IF
      END IF
C
C     ---------------------------------------------------------------
C     Check if there is enough workspace to solve the problem.
C     ---------------------------------------------------------------
      OK = LITOTL .LE. LIWORK .AND. LWTOTL .LE. LWORK
      IF ( .NOT. OK) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99999) LIWORK, LWORK, LITOTL, LWTOTL
            CALL X04BAY(NERR,3,REC)
            WRITE (REC,FMT=99998)
            CALL X04BAY(NERR,2,REC)
         END IF
      ELSE IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999) LIWORK, LWORK, LITOTL, LWTOTL
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (NERROR.EQ.0) THEN
C
C        ---------------------------------------------------------------
C        Check the bounds on all variables and constraints.
C        ---------------------------------------------------------------
         DO 20 J = 1, N + NCLIN + NCNLN
            B1 = BL(J)
            B2 = BU(J)
            OK = B1 .LT. B2 .OR. (B1.EQ.B2 .AND. ABS(B1).LT.BIGBND)
            IF ( .NOT. OK) THEN
               NERROR = NERROR + 1
               IF (J.GT.N+NCLIN) THEN
                  K = J - N - NCLIN
                  L = 3
               ELSE IF (J.GT.N) THEN
                  K = J - N
                  L = 2
               ELSE
                  K = J
                  L = 1
               END IF
               IF (NAMED) THEN
                  IF (B1.EQ.B2) THEN
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99991) NAMES(J), J, J, B1, BIGBND
                        CALL X04BAY(NERR,4,REC)
                     END IF
                  ELSE
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99990) NAMES(J), J, B1, J, B2
                        CALL X04BAY(NERR,3,REC)
                     END IF
                  END IF
               ELSE
                  IF (B1.EQ.B2) THEN
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99989) ID(L), K, J, J, B1, BIGBND
                        CALL X04BAY(NERR,4,REC)
                     END IF
                  ELSE
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99988) ID(L), K, J, B1, J, B2
                        CALL X04BAY(NERR,3,REC)
                     END IF
                  END IF
               END IF
            END IF
   20    CONTINUE
C
C        ---------------------------------------------------------------
C        Check  ISTATE.
C        ---------------------------------------------------------------
         IF (START.EQ.'WARM' .OR. START.EQ.'HOT ') THEN
            IF (START.EQ.'WARM') CSTART = 'Warm'
            IF (START.EQ.'HOT ') CSTART = 'Hot '
            DO 40 J = 1, N + NCLIN + NCNLN
               IS = ISTATE(J)
               OK = IS .GE. (-2) .AND. IS .LE. 4
               IF ( .NOT. OK) THEN
                  NERROR = NERROR + 1
                  IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                     WRITE (REC,FMT=99987) CSTART, J, J, IS
                     CALL X04BAY(NERR,3,REC)
                  END IF
               END IF
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of  E04MFP.  (CMINIT)
C
99999 FORMAT (/' Workspace provided is     IWORK(',I8,'),  WORK(',I8,
     *       ').',/' To solve problem we need  IWORK(',I8,'),  WORK(',
     *       I8,').')
99998 FORMAT (/' ** Not enough workspace to solve problem.')
99997 FORMAT (/' ** On entry, N.le.0:',/'    N = ',I16)
99996 FORMAT (/' ** On entry, NCLIN.lt.0:',/'    NCLIN = ',I16)
99995 FORMAT (/' ** On entry, NCNLN.lt.0:',/'    NCNLN = ',I16)
99994 FORMAT (/' ** On entry, LDA.lt.max(1,NCLIN):',/'    LDA = ',I16,
     *       '   NCLIN = ',I16)
99993 FORMAT (/' ** On entry, LDH.lt.1:',/'    LDH = ',I16)
99992 FORMAT (/' ** On entry, either LDH.lt.N or LDH.lt.m',/'    (wher',
     *       'e m is the value of the optional parameter Hessian Rows):'
     *       ,/'    LDH = ',I16,'   N = ',I16,'   m = ',I16)
99991 FORMAT (/' ** On entry, the equal bounds on  ',A8,'  are infinit',
     *       'e (because',/'    BL(',I4,').eq.beta and BU(',I4,').eq.b',
     *       'eta, but abs(beta).ge.bigbnd):',/'    beta =',G16.7,' bi',
     *       'gbnd =',G16.7)
99990 FORMAT (/' ** On entry, the bounds on  ',A8,'  are inconsistent:',
     *       /'    BL(',I4,') =',G16.7,'   BU(',I4,') =',G16.7)
99989 FORMAT (/' ** On entry, the equal bounds on  ',A5,I3,'  are infi',
     *       'nite (because',/'    BL(',I4,').eq.beta and BU(',I4,').e',
     *       'q.beta, but abs(beta).ge.bigbnd):',/'    beta =',G16.7,
     *       ' bigbnd =',G16.7)
99988 FORMAT (/' ** On entry, the bounds on  ',A5,I3,'  are inconsiste',
     *       'nt:',/'    BL(',I4,') =',G16.7,'   BU(',I4,') =',G16.7)
99987 FORMAT (/' ** On entry with a ',A4,' Start, ISTATE(',I4,') is ou',
     *       't of range:',/'    ISTATE(',I4,') = ',I16)
      END
      SUBROUTINE E04NCL(HITCON,HITLOW,LINOBJ,UNITGZ,NCLIN,NRANK,NRZ,N,
     *                  LDR,JADD,NUMINF,ALFA,CTP,CTX,XNORM,AP,AX,BL,BU,
     *                  GQ,HZ,P,RES,R,X,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1064 (JUL 1993).
C
C     ******************************************************************
C     E04NCL  changes X to X + ALFA*P and updates CTX, AX, RES and GQ
C     accordingly.
C
C     If a bound was added to the working set,  move X exactly on to it,
C     except when a negative step was taken (E04UCG may have had to move
C     to some other closer constraint.)
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 27-December-1985.
C     This version dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CTP, CTX, XNORM
      INTEGER           JADD, LDR, N, NCLIN, NRANK, NRZ, NUMINF
      LOGICAL           HITCON, HITLOW, LINOBJ, UNITGZ
C     .. Array Arguments ..
      DOUBLE PRECISION  AP(*), AX(*), BL(*), BU(*), GQ(*), HZ(*), P(N),
     *                  R(LDR,*), RES(*), WORK(*), X(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRMV
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Executable Statements ..
      CALL DAXPY(N,ALFA,P,1,X,1)
      IF (LINOBJ) CTX = CTX + ALFA*CTP
C
      IF (HITCON .AND. JADD.LE.N) THEN
         BND = BU(JADD)
         IF (HITLOW) BND = BL(JADD)
         IF (ALFA.GE.ZERO) X(JADD) = BND
      END IF
      XNORM = DNRM2(N,X,1)
C
      IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,AP,1,AX,1)
C
      IF (NRZ.LE.NRANK) THEN
         IF (UNITGZ) THEN
            RES(NRZ) = RES(NRZ) - ALFA*HZ(NRZ)
         ELSE
            CALL DAXPY(NRZ,(-ALFA),HZ,1,RES,1)
         END IF
C
         IF (NUMINF.EQ.0) THEN
C
C           Update the transformed gradient GQ so that
C           GQ = GQ + ALFA*R'( HZ ).
C                            ( 0  )
C
            IF (UNITGZ) THEN
               CALL DAXPY(N-NRZ+1,ALFA*HZ(NRZ),R(NRZ,NRZ),LDR,GQ(NRZ),1)
            ELSE
               CALL DCOPY(NRZ,HZ,1,WORK,1)
               CALL DTRMV('U','T','N',NRZ,R,LDR,WORK,1)
               IF (NRZ.LT.N) CALL DGEMV('T',NRZ,N-NRZ,ONE,R(1,NRZ+1),
     *                                  LDR,HZ,1,ZERO,WORK(NRZ+1),1)
C
               CALL DAXPY(N,ALFA,WORK,1,GQ,1)
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of  E04NCL. (LSMOVE)
C
      END
      SUBROUTINE E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSR,FX,INFORM,ITER,ITMAX,
     *                  CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,HPHI)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1108 (JUL 1993).
C
C     ******************************************************************
C     E04XAZ  implements algorithm  FD, the method described in
C     Gill, P.E., Murray, W., Saunders, M.A., and Wright, M. H.,
C     Computing Forward-Difference Intervals for Numerical Optimization,
C     Siam Journal on Scientific and Statistical Computing, vol. 4,
C     pp. 310-321, June 1983.
C
C     The procedure is based on finding an interval (HPHI) that
C     produces an acceptable estimate of the second derivative, and
C     then using that estimate to compute an interval that should
C     produce a reasonable forward-difference approximation.
C
C     One-sided difference estimates are used to ensure feasibility with
C     respect to an upper or lower bound on X. If X is close to an upper
C     bound, the trial intervals will be negative. The final interval is
C     always positive.
C
C     E04XAZ has been designed to use a reverse communication
C     control structure, i.e., all evaluations of the function occur
C     outside this routine. The calling routine repeatedly calls  E04XAZ
C     after computing the indicated function values.
C
C     E04XAZ  is similar to subroutine FDCORE described in Report
C     SOL 83-6, Documentation of FDCORE and FDCALC, by P.E. Gill,
C     W. Murray,  M.A. Saunders, and M.H. Wright, Department of
C     Operations Research,  Stanford University, Stanford, California
C     94305, June 1983.
C
C     Systems Optimization Laboratory, Stanford University.
C     Based on Fortran 66 Version 2.1 of  FDCORE  written June 1983.
C     Fortran 77 Version written 25-May-1985.
C     This version of  E04XAZ  dated  27-Oct-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  BNDLO, BNDUP
      PARAMETER         (BNDLO=1.0D-3,BNDUP=1.0D-1)
      DOUBLE PRECISION  ZERO, SIXTH, FOURTH
      PARAMETER         (ZERO=0.0D+0,SIXTH=1.6D-1,FOURTH=2.5D-1)
      DOUBLE PRECISION  HALF, TWO
      PARAMETER         (HALF=5.0D-1,TWO=2.0D+0)
      DOUBLE PRECISION  THREE, FOUR, TEN
      PARAMETER         (THREE=3.0D+0,FOUR=4.0D+0,TEN=1.0D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CDEST, EPSA, EPSR, ERRBND, F1, F2, FDEST, FX, H,
     *                  HOPT, HPHI, SDEST
      INTEGER           INFORM, ITER, ITMAX
      LOGICAL           DEBUG, DONE, FIRST
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  AFDMIN, CDSAVE, ERR1, ERR2, FDCERR, FDEST2,
     *                  FDSAVE, HSAVE, OLDCD, OLDH, OLDSD, RHO, SDCERR,
     *                  SDSAVE
      LOGICAL           CE1BIG, CE2BIG, OVERFL, TE2BIG
C     .. Local Arrays ..
      CHARACTER*80      REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF
      EXTERNAL          F06BLF
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Save statement ..
      SAVE              CDSAVE, FDSAVE, HSAVE, OLDH, RHO, SDSAVE,
     *                  CE1BIG, CE2BIG, TE2BIG
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Explanation of local variables...
C
C     BNDLO, BNDUP, and RHO control the logic of the routine.
C     BNDLO and BNDUP are the lower and upper bounds that define an
C     acceptable value of the bound on the relative condition error in
C     the second derivative estimate.
C
C     The scalar RHO is the factor by which the interval is multiplied
C     or divided, and also the multiple of the well-scaled interval
C     that is used as the initial trial interval.
C
C     All these values are discussed in the documentation.
C     ------------------------------------------------------------------
C
      ITER = ITER + 1
C
C     Compute the forward-,  backward-,  central-  and second-order
C     difference estimates.
C
      FDEST = F06BLF(F1-FX,H,OVERFL)
      FDEST2 = F06BLF(F2-FX,TWO*H,OVERFL)
C
      OLDCD = CDEST
      CDEST = F06BLF(FOUR*F1-THREE*FX-F2,TWO*H,OVERFL)
C
      OLDSD = SDEST
      SDEST = F06BLF(FX-TWO*F1+F2,H*H,OVERFL)
C
C     Compute  FDCERR  and  SDCERR,  bounds on the relative condition
C     errors in the first and second derivative estimates.
C
      AFDMIN = MIN(ABS(FDEST),ABS(FDEST2))
      FDCERR = F06BLF(EPSA,HALF*ABS(H)*AFDMIN,OVERFL)
      SDCERR = F06BLF(EPSA,FOURTH*ABS(SDEST)*H*H,OVERFL)
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99999) ITER, FX, H, F1, FDEST, F2, FDEST2,
     *     CDEST, SDEST, FDCERR, SDCERR
         CALL X04BAY(IPRINT,6,REC)
      END IF
C
C     ==================================================================
C     Select the correct case.
C     ==================================================================
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First time through.
C        Check whether SDCERR lies in the acceptable range.
C        ------------------------------------------------------------
         FIRST = .FALSE.
         DONE = SDCERR .GE. BNDLO .AND. SDCERR .LE. BNDUP
         TE2BIG = SDCERR .LT. BNDLO
         CE2BIG = SDCERR .GT. BNDUP
         CE1BIG = FDCERR .GT. BNDUP
C
         IF ( .NOT. CE1BIG) THEN
            HSAVE = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF
C
         RHO = EPSR**(-SIXTH)/FOUR
         IF (TE2BIG) THEN
C
C           The truncation error may be too big  (same as saying
C           SDCERR is too small).  Decrease the trial interval.
C
            RHO = TEN*RHO
            OLDH = H
            H = H/RHO
         ELSE IF (CE2BIG) THEN
C
C           SDCERR is too large.  Increase the trial interval.
C
            OLDH = H
            H = H*RHO
         END IF
      ELSE IF (CE2BIG) THEN
C        ---------------------------------------------------------------
C        During the last iteration,  the trial interval was
C        increased in order to decrease SDCERR.
C        ---------------------------------------------------------------
         IF (CE1BIG .AND. FDCERR.LE.BNDUP) THEN
            CE1BIG = .FALSE.
            HSAVE = H
            FDSAVE = FDEST
            CDSAVE = CDEST
            SDSAVE = SDEST
         END IF
C
C        If SDCERR is small enough, accept H.  Otherwise,
C        increase H again.
C
         DONE = SDCERR .LE. BNDUP
         IF ( .NOT. DONE) THEN
            OLDH = H
            H = H*RHO
         END IF
      ELSE IF (TE2BIG) THEN
C        ---------------------------------------------------------------
C        During the last iteration,  the interval was decreased in order
C        to reduce the truncation error.
C        ---------------------------------------------------------------
         DONE = SDCERR .GT. BNDUP
         IF (DONE) THEN
C
C           SDCERR has jumped from being too small to being too
C           large.  Accept the previous value of H.
C
            H = OLDH
            SDEST = OLDSD
            CDEST = OLDCD
         ELSE
C
C           Test whether FDCERR is sufficiently small.
C
            IF (FDCERR.LE.BNDUP) THEN
               CE1BIG = .FALSE.
               HSAVE = H
               FDSAVE = FDEST
               CDSAVE = CDEST
               SDSAVE = SDEST
            END IF
C
C           Check whether SDCERR is in range.
C
            DONE = SDCERR .GE. BNDLO
C
            IF ( .NOT. DONE) THEN
C
C              SDCERR is still too small, decrease H again.
C
               OLDH = H
               H = H/RHO
            END IF
         END IF
C
      END IF
C
C     ==================================================================
C     We have either finished or have a new estimate of H.
C     ==================================================================
      IF (DONE) THEN
C
C        Sufficiently good second-derivative estimate found.
C        Compute the optimal interval.
C
         HPHI = ABS(H)
         HOPT = TWO*SQRT(EPSA)/SQRT(ABS(SDEST))
C
C        ERR1 is the error bound on the forward-difference estimate
C        with the final value of H.  ERR2 is the difference of FDEST
C        and the central-difference estimate with HPHI.
C
         ERR1 = HOPT*ABS(SDEST)
         ERR2 = ABS(FDEST-CDEST)
         ERRBND = MAX(ERR1,ERR2)
C
C        Set INFORM = 4  if the forward- and central-difference
C        estimates are not close.
C
         INFORM = 0
         IF (ERRBND.GT.HALF*ABS(FDEST)) INFORM = 4
      ELSE
C        ---------------------------------------------------------------
C        Check whether the maximum number of iterations has been
C        exceeded.  If not, exit.
C        ---------------------------------------------------------------
         DONE = ITER .GE. ITMAX
         IF (DONE) THEN
            IF (CE1BIG) THEN
C
C              FDCERR was never small.  Probably a constant function.
C
               INFORM = 1
               HPHI = HOPT
               FDEST = ZERO
               CDEST = ZERO
               SDEST = ZERO
               ERRBND = ZERO
            ELSE IF (CE2BIG) THEN
C
C              FDCERR was small,  but SDCERR was never small.
C              Probably a linear or odd function.
C
               INFORM = 2
               HPHI = ABS(HSAVE)
               HOPT = HPHI
               FDEST = FDSAVE
               CDEST = CDSAVE
               SDEST = ZERO
               ERRBND = TWO*EPSA/HOPT
            ELSE
C
C              The only remaining case occurs when the second
C              derivative is changing too rapidly for an adequate
C              interval to be found (SDCERR remained small even
C              though H was decreased ITMAX times).
C
               INFORM = 3
               HPHI = ABS(HSAVE)
               HOPT = HPHI
               FDEST = FDSAVE
               CDEST = CDSAVE
               SDEST = SDSAVE
               ERRBND = HOPT*ABS(SDEST)/TWO + TWO*EPSA/HOPT
            END IF
         END IF
      END IF
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99998) CE1BIG, CE2BIG, TE2BIG
         CALL X04BAF(IPRINT,REC(1))
         IF (DONE) THEN
            WRITE (REC,FMT=99997) INFORM, HOPT, ERRBND
            CALL X04BAF(IPRINT,REC(1))
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04XAZ. (CHCORE)
C
99999 FORMAT (/' //E04XAZ//  ITN ',I3,' FX     H',11X,1P,2D16.6,/' //E',
     *       '04XAZ//  F1      FDEST',14X,1P,2D16.6,/' //E04XAZ//  F2 ',
     *       '     FDEST2',13X,1P,2D16.6,/' //E04XAZ//  CDEST   SDEST',
     *       14X,1P,2D16.6,/' //E04XAZ//  FDCERR  SDCERR',13X,1P,2D16.6)
99998 FORMAT (' //E04XAZ//  CE1BIG  CE2BIG  TE2BIG',5X,3L2)
99997 FORMAT (' //E04XAZ//  INFORM  HOPT    ERRBND',I5,1P,2D16.6)
      END
      LOGICAL FUNCTION E04UDX(STRING)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C***********************************************************************
C     Description and usage:
C
C        A simple(-minded) test for numeric data is implemented by
C        searching an input string for legitimate characters:
C                digits 0 to 9, D, E, -, + and .
C        Insurance is provided by requiring that a numeric string
C        have at least one digit, at most one D, E or .
C        and at most two -s or +s.  Note that a few ambiguities remain:
C
C           (a)  A string might have the form of numeric data but be
C                intended as text.  No general test can hope to detect
C                such cases.
C
C           (b)  There is no check for correctness of the data format.
C                For example a meaningless string such as 'E1.+2-'
C                will be accepted as numeric.
C
C        Despite these weaknesses, the method should work in the
C        majority of cases.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        E04UDX              L      O    Set .TRUE. if STRING appears
C                                        to be numerical data.
C        STRING              C    I      Input data to be tested.
C
C
C     Environment:  ANSI FORTRAN 77.
C
C
C     Notes:
C
C        (1)  It is assumed that STRING is a token extracted by
C             E04UDV, which will have converted any lower-case
C             characters to upper-case.
C
C        (2)  E04UDV pads STRING with blanks, so that a genuine
C             number is of the form  '1234        '.
C             Hence, the scan of STRING stops at the first blank.
C
C        (3)  COMPLEX data with parentheses will not look numeric.
C
C
C     Systems Optimization Laboratory, Stanford University.
C     12 Nov  1985    Initial design and coding, starting from the
C                     routine ALPHA from Informatics General, Inc.
C***********************************************************************
C
C     .. Scalar Arguments ..
      CHARACTER*(*)           STRING
C     .. Local Scalars ..
      INTEGER                 J, LENGTH, NDIGIT, NEXP, NMINUS, NPLUS,
     *                        NPOINT
      LOGICAL                 NUMBER
      CHARACTER*1             ATOM
C     .. Intrinsic Functions ..
      INTRINSIC               LEN, LGE, LLE
C     .. Executable Statements ..
      NDIGIT = 0
      NEXP = 0
      NMINUS = 0
      NPLUS = 0
      NPOINT = 0
      NUMBER = .TRUE.
      LENGTH = LEN(STRING)
      J = 0
C
   20 J = J + 1
      ATOM = STRING(J:J)
      IF (LGE(ATOM,'0') .AND. LLE(ATOM,'9')) THEN
C        IF (ATOM.GE.'0' .AND. ATOM.LE.'9') THEN
         NDIGIT = NDIGIT + 1
      ELSE IF (ATOM.EQ.'D' .OR. ATOM.EQ.'E' .or.
     *         ATOM.EQ.'d' .OR. ATOM.EQ.'e') THEN
         NEXP = NEXP + 1
      ELSE IF (ATOM.EQ.'-') THEN
         NMINUS = NMINUS + 1
      ELSE IF (ATOM.EQ.'+') THEN
         NPLUS = NPLUS + 1
      ELSE IF (ATOM.EQ.'.') THEN
         NPOINT = NPOINT + 1
      ELSE IF (ATOM.EQ.' ') THEN
         J = LENGTH
      ELSE
         NUMBER = .FALSE.
      END IF
C
      IF (NUMBER .AND. J.LT.LENGTH) GO TO 20
C
      E04UDX = NUMBER .AND. NDIGIT .GE. 1 .AND. NEXP .LE. 1 .AND.
     *         NMINUS .LE. 2 .AND. NPLUS .LE. 2 .AND. NPOINT .LE. 1
C
      RETURN
C
C     End of  E04UDX. (OPNUMB)
      END
      SUBROUTINE E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,
     *                  PALFA1,PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,
     *                  FEATOL,P,X)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1079 (JUL 1993).
C
C     ******************************************************************
C     E04UCH  finds steps PALFA1, PALFA2 such that
C        X + PALFA1*P  reaches a linear constraint that is currently not
C                      in the working set but is satisfied.
C        X + PALFA2*P  reaches a linear constraint that is currently not
C                      in the working set but is violated.
C     The constraints are perturbed by an amount FEATOL, so that PALFA1
C     is slightly larger than it should be,  and PALFA2 is slightly
C     smaller than it should be.  This gives some leeway later when the
C     exact steps are computed by E04UCG.
C
C     Constraints in the working set are ignored  (ISTATE(j) .GE. 1).
C
C     If NEGSTP is true, the search direction will be taken to be  - P.
C
C
C     Values of ISTATE(j)....
C
C        - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     The values  -2  and  -1  do not occur once a feasible point has
C     been found.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written  May 1980.
C     This version of E04UCH dated 26-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGALF, BIGBND, PALFA1, PALFA2, PNORM
      INTEGER           JADD1, JADD2, N, NCTOTL
      LOGICAL           FIRSTV, NEGSTP
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), P(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSATP, ATP, ATX, RES, ROWNRM
      INTEGER           I, J, JS
      LOGICAL           LASTV
C     .. Local Arrays ..
      CHARACTER*120     REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
      END IF
      LASTV = .NOT. FIRSTV
      JADD1 = 0
      JADD2 = 0
      PALFA1 = BIGALF
C
      PALFA2 = ZERO
      IF (FIRSTV) PALFA2 = BIGALF
C
      DO 20 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS.LE.0) THEN
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ROWNRM = ONE
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ROWNRM = ONE + ANORM(I)
            END IF
            IF (NEGSTP) ATP = -ATP
C
            IF (ABS(ATP).LE.EPSPT9*ROWNRM*PNORM) THEN
C
C              This constraint appears to be constant along P.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C
               RES = -ONE
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing and the lower bound is not violated.
C              ---------------------------------------------------------
C              First test for smaller PALFA1.
C
               ABSATP = -ATP
               IF (BL(J).GT.(-BIGBND)) THEN
                  RES = ATX - BL(J) + FEATOL(J)
                  IF (BIGALF*ABSATP.GT.ABS(RES)) THEN
                     IF (PALFA1*ABSATP.GT.RES) THEN
                        PALFA1 = RES/ABSATP
                        JADD1 = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 The upper bound is violated.  Test for either larger
C                 or smaller PALFA2, depending on the value of FIRSTV.
C
                  RES = ATX - BU(J) - FEATOL(J)
                  IF (BIGALF*ABSATP.GT.ABS(RES)) THEN
                     IF (FIRSTV .AND. PALFA2*ABSATP.GT.RES .OR.
     *                   LASTV .AND. PALFA2*ABSATP.LT.RES) THEN
                        PALFA2 = RES/ABSATP
                        JADD2 = J
                     END IF
                  END IF
               END IF
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              Test for smaller PALFA1.
C
               IF (BU(J).LT.BIGBND) THEN
                  RES = BU(J) - ATX + FEATOL(J)
                  IF (BIGALF*ATP.GT.ABS(RES)) THEN
                     IF (PALFA1*ATP.GT.RES) THEN
                        PALFA1 = RES/ATP
                        JADD1 = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 The lower bound is violated.  Test for a new PALFA2.
C
                  RES = BL(J) - ATX - FEATOL(J)
                  IF (BIGALF*ATP.GT.ABS(RES)) THEN
                     IF (FIRSTV .AND. PALFA2*ATP.GT.RES .OR. LASTV .AND.
     *                   PALFA2*ATP.LT.RES) THEN
                        PALFA2 = RES/ATP
                        JADD2 = J
                     END IF
                  END IF
               END IF
            END IF
C
            IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
               WRITE (REC,FMT=99998) J, JS, FEATOL(J), RES, ATP, JADD1,
     *           PALFA1, JADD2, PALFA2
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
   20 CONTINUE
C
      RETURN
C
C
C     End of  E04UCH. (CMALF1)
C
99999 FORMAT (/'    J  JS         FEATOL        RES             AP    ',
     *       ' JADD1       PALFA1     JADD2       PALFA2',/)
99998 FORMAT (I5,I4,3G15.5,2(I6,G17.7))
      END
      SUBROUTINE E04UDU(STRING)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     PURPOSE:  This subroutine changes all lower case letters in the
C               character string to upper case.
C
C     METHOD:   Each character in STRING is treated in turn.  The
C               intrinsic function INDEX effectively allows a table
C               lookup, with the local strings LOW and UPP acting as
C               two tables. This method avoids the use of CHAR and
C               ICHAR, which appear be different on ASCII and EBCDIC
C               machines.
C
C     ARGUMENTS
C     ARG       DIM     TYPE I/O/S DESCRIPTION
C     STRING       *       C   I/O   Character string possibly
C                                   containing some lower-case
C                                   letters  on input; strictly
C                                   upper-case letters on output
C                                   with no change to any
C                                   non-alphabetic characters.
C
C     EXTERNAL REFERENCES:
C     LEN    - Returns the declared length of a CHARACTER variable.
C     INDEX  - Returns the position of second string within first.
C
C     ENVIRONMENT:  ANSI FORTRAN 77
C
C     DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C     06/28/83   CLH    Initial design.
C     01/03/84   RAK    Eliminated NCHAR input.
C     06/14/84   RAK    Used integer PARAMETERs in comparison.
C     04/21/85   RAK    Eliminated DO/END DO in favor of standard code.
C     09/10/85   MAS    Eliminated CHAR,ICHAR in favor of LOW, UPP,
C                       INDEX.
C
C     AUTHOR: Charles Hooper, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      CHARACTER*(*)     STRING
C     .. Local Scalars ..
      INTEGER           I, J
      CHARACTER*1       C
      CHARACTER*26      LOW, UPP
C     .. Intrinsic Functions ..
      INTRINSIC         INDEX, LEN, LGE, LLE
C     .. Data statements ..
      DATA              LOW/'abcdefghijklmnopqrstuvwxyz'/,
     *                  UPP/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C     .. Executable Statements ..
C
      DO 20 J = 1, LEN(STRING)
         C = STRING(J:J)
         IF (LGE(C,'a') .AND. LLE(C,'z')) THEN
C           IF (C.GE.'a' .AND. C.LE.'z') THEN
            I = INDEX(LOW,C)
            IF (I.GT.0) STRING(J:J) = UPP(I:I)
         END IF
   20 CONTINUE
      RETURN
C
C     End of  E04UDU. (OPUPPR)
C
      END

      SUBROUTINE E04NBV(N,NU,NRANK,LDR,LENV,LENW,R,U,V,W,C,S)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-589 (MAR 1988).
C     MARK 16 REVISED. IER-1057 (JUL 1993).
C
C     ==================================================================
C     E04NBV  modifies the  nrank*n  upper-triangular matrix  R  so that
C     Q*(R + v*w')  is upper triangular,  where  Q  is orthogonal,
C     v  and  w  are vectors, and the modified  R  overwrites the old.
C     Q  is the product of two sweeps of plane rotations (not stored).
C     If required,  the rotations are applied to the NU columns of
C     the matrix  U.
C
C     The matrix v*w' is an (LENV) by (LENW) matrix.
C     The vector v is overwritten.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   October  1984.
C     Level-2 matrix routines added 22-Apr-1988.
C     This version of  E04NBV  dated 22-Apr-1988.
C     ==================================================================
C     .. Scalar Arguments ..
      INTEGER           LDR, LENV, LENW, N, NRANK, NU
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), R(LDR,*), S(N), U(N,*), V(N), W(N)
C     .. Local Scalars ..
      INTEGER           J
C     .. External Subroutines ..
      EXTERNAL          DAXPY, F06FQF, F06QSF, F06QWF, F06QXF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      J = MIN(LENV,NRANK)
      IF (NRANK.GT.0) THEN
C        ---------------------------------------------------------------
C        Reduce  v to beta*e( j )  using a backward sweep of rotations
C        in planes (j-1, j), (j-2, j), ..., (1, j).
C        ---------------------------------------------------------------
         CALL F06FQF('Fixed','Backwards',J-1,V(J),V,1,C,S)
C
C        ---------------------------------------------------------------
C        Apply the sequence of rotations to U.
C        ---------------------------------------------------------------
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Backwards',J,NU,1,J,
     *                            C,S,U,N)
C
C        ---------------------------------------------------------------
C        Apply the sequence of rotations to R. This generates a spike in
C        the j-th row of R, which is stored in s.
C        ---------------------------------------------------------------
         CALL F06QWF('Left',N,1,J,C,S,R,LDR)
C
C        ---------------------------------------------------------------
C        Form  beta*e(j)*w' + R.  This a spiked matrix, with a row
C        spike in row j.
C        ---------------------------------------------------------------
         CALL DAXPY(MIN(J-1,LENW),V(J),W,1,S,1)
         CALL DAXPY(LENW-J+1,V(J),W(J),1,R(J,J),LDR)
C
C        ---------------------------------------------------------------
C        Eliminate the spike using a forward sweep of rotations in
C        planes (1, j), (2, j), ..., (j-1, j).
C        ---------------------------------------------------------------
         CALL F06QSF('Left',N,1,J,C,S,R,LDR)
C
C        ---------------------------------------------------------------
C        Apply the rotations to U.
C        ---------------------------------------------------------------
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Forwards',J,NU,1,J,C,
     *                            S,U,N)
      END IF
C
C     End of  E04NBV. (CMR1MD)
C
      END
      SUBROUTINE E04MFY(NRZ,LDR,R,RZZ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFY  loads the last column of the  NRZ x NRZ  triangular factor
C     Rz  with the multiple  RZZ  of the  NRZ-th unit vector.
C
C     Original version written by PEG,  23-Jul-87.
C     This version of  E04MFY  dated 17-Jul-90.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RZZ
      INTEGER           LDR, NRZ
C     .. Array Arguments ..
      DOUBLE PRECISION  R(LDR,*)
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Executable Statements ..
C
      IF (NRZ.EQ.0) RETURN
C
      CALL F06FBF(NRZ-1,ZERO,R(1,NRZ),1)
      R(NRZ,NRZ) = RZZ
C
      RETURN
C
C     End of  E04MFY.  (LPCOLR)
C
      END
      SUBROUTINE E04UDT(INFORM,N,NCLIN,NCNLN,ALFA,ALFMIN,ALFMAX,BIGBND,
     *                  DXNORM,ANORM,ADX,AX,BL,BU,DSLK,DX,SLK,X)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1096 (JUL 1993).
C
C     ******************************************************************
C     E04UDT  finds a step ALFA such that the point X + ALFA*P reaches
C     one of the slacks or linear constraints.  The step ALFA is the
C     maximum step that can be taken without violating one of the slacks
C     or linear constraints that is currently satisfied.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written  June 1986.
C     This version of E04UDT dated  13-Jun-1987.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFMAX, ALFMIN, BIGBND, DXNORM
      INTEGER           INFORM, N, NCLIN, NCNLN
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), ANORM(*), AX(*), BL(*), BU(*), DSLK(*),
     *                  DX(N), SLK(*), X(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ADXI, AXI, RES, ROWNRM
      INTEGER           I, J
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,4,REC)
      END IF
C
      ALFA = ALFMAX
      J = 1
C
C     +    WHILE (J .LE. N+NCLIN+NCNLN .AND. ALFA .GT. ALFMIN) DO
   20 IF (J.LE.N+NCLIN+NCNLN .AND. ALFA.GT.ALFMIN) THEN
C
         IF (J.LE.N) THEN
            AXI = X(J)
            ADXI = DX(J)
            ROWNRM = ONE
         ELSE IF (J.LE.N+NCLIN) THEN
C
            I = J - N
            AXI = AX(I)
            ADXI = ADX(I)
            ROWNRM = ANORM(I) + ONE
         ELSE
C
            I = J - N - NCLIN
            AXI = SLK(I)
            ADXI = DSLK(I)
            ROWNRM = ONE
         END IF
C
         RES = -ONE
         IF (ADXI.LE.-EPSPT9*ROWNRM*DXNORM) THEN
C
C           Constraint decreasing.
C
            ADXI = -ADXI
            IF (BL(J).GT.-BIGBND) RES = AXI - BL(J)
         ELSE IF (ADXI.GT.EPSPT9*ROWNRM*DXNORM) THEN
C
C           Constraint increasing.
C
            IF (BU(J).LT.BIGBND) RES = BU(J) - AXI
         END IF
C
         IF (RES.GT.ZERO .AND. ALFA*ADXI.GT.RES) ALFA = RES/ADXI
C
         IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
            WRITE (REC,FMT=99998) J, RES, ADXI, ALFA
            CALL X04BAF(IPRINT,REC(1))
         END IF
C
         J = J + 1
         GO TO 20
C        +    END WHILE
      END IF
C
C     ==================================================================
C     Determine ALFA, the bound on the step to be taken.
C     ==================================================================
      ALFA = MAX(ALFA,ALFMIN)
C
      INFORM = 0
      IF (ALFA.GE.ALFMAX) INFORM = 1
C
      IF (CMDBG .AND. ICMDBG(1).GT.0 .AND. INFORM.GT.0) THEN
         WRITE (REC,FMT=99997) ALFA
         CALL X04BAY(IPRINT,4,REC)
      END IF
C
      RETURN
C
C
C     End of  E04UDT. (NPALF)
C
99999 FORMAT (/' E04UDT entered',/'    J            RES             AP',
     *       '           ALFA ',/)
99998 FORMAT (I5,3G15.5)
99997 FORMAT (/' //E04UDT//  No finite step.',/' //E04UDT//           ',
     *       '  ALFA',/' //E04UDT//  ',G15.4)
      END
      SUBROUTINE E04UDY(NDICT,DICTRY,ALPHA,KEY,ENTRY)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-724 (DEC 1989).
C
C     Description and usage:
C
C       Performs dictionary lookups.  A pointer is returned if a
C     match is found between the input key and the corresponding
C     initial characters of one of the elements of the dictionary.
C     If a 'synonym' has been provided for an entry, the search is
C     continued until a match to a primary dictionary entry is found.
C     Cases of no match, or multiple matches, are also provided for.
C
C     Dictionary entries must be left-justified, and may be alphabetized
C     for faster searches.  Secondary entries, if any, are composed of
C     two words separated by one or more characters such as blank, tab,
C     comma, colon, or equal sign which are treated as non-significant
C     by E04UDW.  The first entry of each such pair serves as a synonym
C     for the second, more fundamental keyword.
C
C       The ordered search stops after the section of the dictionary
C     having the same first letters as the key has been checked, or
C     after a specified number of entries have been examined.  A special
C     dictionary entry, the currency symbol '$', will also terminate the
C     search.  This will speed things up if an appropriate dictionary
C     length parameter cannot be determined.  Both types of search are
C     sequential.  See 'Notes' below for some suggestions if efficiency
C     is an issue.
C
C
C     Parameters:
C
C     Name    Dimension  Type  I/O/S  Description
C     NDICT               I    I     Number of dictionary entries to be
C                                    examined.
C     DICTRY  NDICT       C    I     Array of dictionary entries,
C                                    left-justified in their fields.
C                                    May be alphabetized for efficiency,
C                                    in which case ALPHA should be
C                                    .TRUE.  Entries with synonyms are
C                                    of the form
C                                    'ENTRY : SYNONYM', where 'SYNONYM'
C                                    is a more fundamental entry in the
C                                    same dictionary.  NOTE: Don't build
C                                    'circular' dictionaries.
C     ALPHA               L    I     Indicates whether the dictionary
C                                    is in alphabetical order, in which
C                                    case the search can be terminated
C                                    sooner.
C     KEY                 C    I/O   String to be compared against the
C                                    dictionary.  Abbreviations are OK
C                                    if they correspond to a unique
C                                    entry in the dictionary.  KEY is
C                                    replaced on termination by its most
C                                    fundamental equivalent dictionary
C                                    entry (uppercase, left-justified)
C                                    if a match was found.
C     ENTRY               I      O   Dictionary pointer.  If .GT. 0, it
C                                    indicates which entry matched KEY.
C                                    In case of trouble, a negative
C                                    value means that a UNIQUE match
C                                    was not found - the absolute value
C                                    of ENTRY points to the second
C                                    dictionary entry that matched KEY.
C                                    Zero means that NO match could be
C                                    found.  ENTRY always refers to the
C                                    last search performed -
C                                    in searching a chain of synonyms,
C                                    a non-positive value will be
C                                    returned if there is any break,
C                                    even if the original input key
C                                    was found.
C
C
C     External references:
C
C     Name    Description
C     E04UDW  Finds first and last significant characters.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C     Notes:
C
C     (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C     (2)  We have assumed that the dictionary is not too big.  If
C         many searches are to be done or if the dictionary has more
C         than a dozen or so entries, it may be advantageous to build
C         an index array of pointers to the beginning of the section
C         of the dictionary containing each letter, then pass in the
C         portion of the dictionary beginning with DICTRY (INDEX).
C         (This won't generally work for dictionaries with synonyms.)
C         For very large problems, a completely different approach may
C         be advisable, e.g. a binary search for ordered dictionaries.
C
C     (3)  E04UDY is case sensitive.  In most applications it will be
C         necessary to use an uppercase dictionary, and to convert the
C         input key to uppercase before calling E04UDY.  Companion
C         routines E04UDV and PAIRS, available from the author, already
C         take care of this.
C
C     (4)  The key need not be left-justified.  Any leading (or
C         trailing) characters which are 'non-significant' to E04UDW
C         will be ignored.  These include blanks, horizontal tabs,
C         commas, colons, and equal signs.  See E04UDW for details.
C
C     (5)  The ASCII collating sequence for character data is assumed.
C         (N.B. This means the numerals precede the alphabet, unlike
C         common practice.)  This should not cause trouble on EBCDIC
C         machines if DICTRY just contains alphabetic keywords.
C         Otherwise it may be necessary to use the FORTRAN lexical
C         library routines to force use of the ASCII sequence.
C
C     (6)  Parameter NUMSIG sets a limit on the length of significant
C         dictionary entries.  Special applications may require that
C         this be increased.  (It is 16 in the present version.)
C
C     (7)  No protection against 'circular' dictionaries is provided:
C         don't claim that A is B, and that B is A.  All synonym chains
C         must terminate.  Other potential errors not checked for
C         include duplicate or mis-ordered entries.
C
C     (8)  The handling of ambiguities introduces some ambiguity:
C
C            ALPHA = .TRUE.  A potential problem, when one entry
C                            looks like an abbreviation for another
C                            (eg. does 'A' match 'A' or 'AB') was
C                            resolved by dropping out of the search
C                            immediately when an 'exact' match is found.
C
C            ALPHA = .FALSE. The programmer must ensure that the above
C                            situation does not arise: each dictionary
C                            entry must be recognizable, at least when
C                            specified to full length.  Otherwise, the
C                            result of a search will depend on the
C                            order of entries.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C     24 Feb. 1984  RAK/DAS  Initial design and coding.
C     25 Feb. 1984    RAK    Combined the two searches by suitable
C                            choice of terminator FLAG.
C     28 Feb. 1984    RAK    Optional synonyms in dictionary, no
C                            longer update KEY.
C     29 Mar. 1984    RAK    Put back replacement of KEY by its
C                            corresponding entry.
C     21 June 1984    RAK    Corrected bug in error handling for cases
C                            where no match was found.
C     23 Apr. 1985    RAK    Introduced test for exact matches, which
C                            permits use of dictionary entries which
C                            would appear to be ambiguous (for ordered
C                            case).  Return -I to point to the entry
C                            which appeared ambiguous (had been -1).
C                            Repaired loop termination - had to use
C                            equal length strings or risk quitting too
C                            soon when one entry is an abbreviation
C                            for another.  Eliminated HIT, reduced
C                            NUMSIG to 16.
C     28 May 1986     MPH    Changed to test if an apparently
C                            ambiguous keyword is in fact unique.
C                            If ALPHA is .FALSE. FLAG is now '$'.
C                            Local CHARACTER scalars declared as
C                            CHARACTER*16.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER         BLANK, CURLY
      INTEGER           NUMSIG
      PARAMETER         (BLANK=' ',CURLY='$',NUMSIG=16)
C     .. Scalar Arguments ..
      INTEGER           ENTRY, NDICT
      LOGICAL           ALPHA
      CHARACTER*(*)     KEY
C     .. Array Arguments ..
      CHARACTER*(*)     DICTRY(NDICT)
C     .. Local Scalars ..
      INTEGER           FIRST, I, IFRST, IFRST1, ILAST, ILAST1, ILEN,
     *                  ILEN1, ILST, IMARK, IMARK1, LAST, LENGTH, MARK
      CHARACTER*16      FLAG, TARGET, TRGT, TRGT1
C     .. External Subroutines ..
      EXTERNAL          E04UDW
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, LEN, LLE
C     .. Executable Statements ..
C
      ENTRY = 0
C
C     Isolate the significant portion of the input key (if any).
C
      FIRST = 1
      LAST = MIN(LEN(KEY),NUMSIG)
      CALL E04UDW(KEY,FIRST,LAST,MARK)
C
      IF (MARK.GT.0) THEN
         TARGET = KEY(FIRST:MARK)
C
C        Look up TARGET in the dictionary.
C
   20    CONTINUE
         LENGTH = MARK - FIRST + 1
C
C           Select search strategy by cunning choice of termination test
C           flag.  The left curly bracket follows all the alphabetic
C           characters in the ASCII collating sequence, but precedes the
C           vertical bar.
C
         IF (ALPHA) THEN
            FLAG = TARGET
         ELSE
            FLAG = CURLY
         END IF
C
C
C           Perform search.
C           ---------------
C
         I = 0
   40    CONTINUE
         I = I + 1
         IF (TARGET(1:LENGTH).EQ.DICTRY(I)(1:LENGTH)) THEN
            IF (ENTRY.EQ.0) THEN
C
C                    First 'hit' - must still guard against ambiguities
C                    by searching until we've gone beyond the key
C                    (ordered dictionary) or until the end-of-dictionary
C                    mark is reached (exhaustive search).
C
               ENTRY = I
C
C                    Special handling if match is exact - terminate
C                    search.  We thus avoid confusion if one dictionary
C                    entry looks like an abbreviation of another.
C                    This fix won't generally work for un-ordered
C                    dictionaries.
C
               FIRST = 1
               LAST = NUMSIG
               CALL E04UDW(DICTRY(ENTRY),FIRST,LAST,MARK)
               IF (MARK.EQ.LENGTH) I = NDICT
            ELSE
C                    If two hits check if they are attempting to
C                    indicate the same dictionary entry.
C
C                    Extract keyword from first match found
C
               ILST = NUMSIG
               IFRST = MARK + 2
               CALL E04UDW(DICTRY(ENTRY),IFRST,ILST,IMARK)
               IF (IMARK.GT.0) THEN
                  TRGT = DICTRY(ENTRY) (IFRST:IMARK)
                  ILEN = IMARK - IFRST + 1
               ELSE
                  TRGT = DICTRY(ENTRY) (FIRST:MARK)
                  ILEN = MARK - FIRST + 1
               END IF
C
C                    Extract keyword from next match found
C
               IFRST = 1
               ILAST = NUMSIG
               CALL E04UDW(DICTRY(I),IFRST,ILAST,IMARK)
               ILAST1 = NUMSIG
               IFRST1 = IMARK + 2
               CALL E04UDW(DICTRY(I),IFRST1,ILAST1,IMARK1)
               IF (IMARK1.GT.0) THEN
                  TRGT1 = DICTRY(I) (IFRST1:IMARK1)
                  ILEN1 = IMARK1 - IFRST1 + 1
               ELSE
                  TRGT1 = DICTRY(I) (IFRST:IMARK)
                  ILEN1 = IMARK - IFRST + 1
               END IF
C
C                    If keywords not identical then ambiguity
C
               IF (TRGT(1:ILEN).NE.TRGT1(1:ILEN1)) THEN
C
C                       Oops - two hits.  Abnormal termination.
C                       ---------------------------------------
C
                  ENTRY = -I
                  RETURN
               END IF
            END IF
         END IF
C
C           Check whether we've gone past the appropriate section of the
C           dictionary.  The test on the index provides insurance and an
C           optional means for limiting the extent of the search.
C
         IF (LLE(DICTRY(I)(1:LENGTH),FLAG) .AND. I.LT.NDICT) GO TO 40
C
C
C           Check for a synonym.
C           --------------------
C
         IF (ENTRY.GT.0) THEN
C
C              Look for a second entry 'behind' the first entry.  FIRST
C              and MARK were determined above when the hit was detected.
C
            FIRST = MARK + 2
            CALL E04UDW(DICTRY(ENTRY),FIRST,LAST,MARK)
            IF (MARK.GT.0) THEN
C
C                 Re-set target and dictionary pointer, then repeat the
C                 search for the synonym instead of the original key.
C
               TARGET = DICTRY(ENTRY) (FIRST:MARK)
               ENTRY = 0
               GO TO 20
C
            END IF
         END IF
C
      END IF
      IF (ENTRY.GT.0) KEY = DICTRY(ENTRY)
C
C
C     Normal termination.
C     -------------------
C
      RETURN
C
C     End of E04UDY.  (CMLOOK/OPLOOK)
      END
      SUBROUTINE E04NBT(MODE,NROWT,N,T,Y)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-587 (MAR 1988).
C
C     ******************************************************************
C     E04NBT  solves equations involving a reverse-triangular matrix  T
C     and a right-hand-side vector  y,  returning the solution in  y.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written February-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NROWT
C     .. Array Arguments ..
      DOUBLE PRECISION  T(NROWT,*), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  YJ
      INTEGER           J, JJ, L, N1
C     .. External Subroutines ..
      EXTERNAL          DAXPY
C     .. Executable Statements ..
C
      N1 = N + 1
      IF (MODE.EQ.1) THEN
C
C        Mode = 1  ---  Solve  T * y(new) = y(old).
C
         DO 20 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(J,JJ)
            Y(J) = YJ
            L = JJ - 1
            IF (L.GT.0 .AND. YJ.NE.ZERO) CALL DAXPY(L,(-YJ),T(J+1,JJ),1,
     *          Y(J+1),1)
   20    CONTINUE
      ELSE
C
C        Mode = 2  ---  Solve  T' y(new) = y(old).
C
         DO 40 J = 1, N
            JJ = N1 - J
            YJ = Y(J)/T(JJ,J)
            Y(J) = YJ
            L = JJ - 1
            IF (L.GT.0 .AND. YJ.NE.ZERO) CALL DAXPY(L,(-YJ),T(JJ,J+1),
     *          NROWT,Y(J+1),1)
   40    CONTINUE
      END IF
C
C     Reverse the solution vector.
C
      IF (N.GT.1) THEN
         L = N/2
         DO 60 J = 1, L
            JJ = N1 - J
            YJ = Y(J)
            Y(J) = Y(JJ)
            Y(JJ) = YJ
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of  E04NBT. (CMTSOL)
C
      END
      SUBROUTINE E04NBU(N,NU,NRANK,LDR,I,J,R,U,C,S)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     This version dated 8-June-1988. (F06 routines included.)
C
C
C***********************************************************************
C     E04NBU  interchanges the  I-th  and  J-th  (I .LT. J)  columns of
C     an  NRANK*N  upper-trapezoidal matrix  R   and restores the
C     resulting matrix to upper-trapezoidal form using two sweeps of
C     plane rotations applied on the left.  R is overwritten.
C
C     If NU .GT. 0,  the rotations are applied to the  nu  columns of
C     the matrix  U.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level-2 matrix routines added 13-May-1988.
C     This version of  E04NBU  dated  30-May-1988.
C***********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           I, J, LDR, N, NRANK, NU
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), R(LDR,*), S(N), U(N,*)
C     .. Local Scalars ..
      INTEGER           LENJ
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F06FQF, F06QSF, F06QWF, F06QXF, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Swap the elements of the i-th and j-th columns of R on, or above,
C     the main diagonal.
C
      CALL DSWAP(MIN(I,NRANK),R(1,I),1,R(1,J),1)
      LENJ = MIN(J,NRANK)
C
      IF (LENJ.GT.I) THEN
C        ---------------------------------------------------------------
C        Reduce elements  r(i+1,j), ..., r(lenj,j)  to  beta*e(lenj)
C        using a backward sweep in planes
C        (lenj-1,lenj), (lenj-2,lenj), ..., (i+1,lenj).
C        If required, apply the sequence of rotations to U.
C        ---------------------------------------------------------------
         CALL F06FQF('Fixed','Backwards',LENJ-I-1,R(LENJ,J),R(I+1,J),1,
     *               C(I+1),S(I+1))
C
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Backwards',N,NU,I+1,
     *                            LENJ,C,S,U,N)
C
C        Put zeros into the j-th column of R in positions corresponding
C        to the sub-diagonals of the i-th column.
C
         S(I) = R(LENJ,J)
         CALL F06FBF(LENJ-I,ZERO,R(I+1,J),1)
C
C        Apply the sequence of rotations to R.  This generates a spike
C        in the lenj-th row of R, which is stored in S.
C
         CALL F06QWF('Left',N,I+1,LENJ,C,S,R,LDR)
C
C        Eliminate the spike using a forward sweep in planes
C        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
C        If necessary, apply the sequence of rotations to U.
C
         CALL F06QSF('Left',N,I,LENJ,C,S,R,LDR)
C
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Forwards',LENJ,NU,I,
     *                            LENJ,C,S,U,N)
      END IF
C
      RETURN
C
C     End of  E04NBU
C
      END
      SUBROUTINE E04NCT(UNITQ,N,NACTIV,NFREE,NRES,NGQ,NZ,NRZ,LDA,LDZY,
     *                  LDR,LDT,NRANK,JDEL,KDEL,KACTIV,KX,A,RES,R,T,GQ,
     *                  ZY,C,S)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-950 (NOV 1991).
C     MARK 16 REVISED. IER-1071 (JUL 1993).
C     ******************************************************************
C     E04NCT  updates the least-squares factor R and the factorization
C     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
C     constraint is deleted from the working set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level-2 matrix routines added 25-Apr-1988.
C     This version of E04NCT dated 10-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           JDEL, KDEL, LDA, LDR, LDT, LDZY, N, NACTIV,
     *                  NFREE, NGQ, NRANK, NRES, NRZ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), ZY(LDZY,*)
      INTEGER           KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, SN
      INTEGER           I, IR, ITDEL, JART, K, KA, LD, NPIV, NRZ1, NSUP,
     *                  NT
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, E04NBU, F06BAF, F06FBF, F06FLF,
     *                  F06QTF, F06QXF, F06QZZ, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
      IF (JDEL.GT.0) THEN
C        ---------------------------------------------------------------
C        Regular constraint or temporary bound deleted.
C        ---------------------------------------------------------------
C
         IF (JDEL.LE.N) THEN
C
C           Case 1.  A simple bound has been deleted.
C           =======  Columns NFREE+1 and IR of R must be swapped.
C
            IR = NZ + KDEL
            IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
               WRITE (REC,FMT=99998) NACTIV, NZ, NFREE, IR, JDEL, UNITQ
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
            ITDEL = 1
            NFREE = NFREE + 1
C
            IF (NFREE.LT.IR) THEN
               KX(IR) = KX(NFREE)
               KX(NFREE) = JDEL
               IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,NFREE,IR,R,
     *                                     RES,C,S)
               CALL DSWAP(NGQ,GQ(NFREE,1),N,GQ(IR,1),N)
            END IF
C
            IF ( .NOT. UNITQ) THEN
C
C              Copy the incoming column of  A(free)  into the end of T.
C
               DO 20 KA = 1, NACTIV
                  I = KACTIV(KA)
                  T(KA,NFREE) = A(I,JDEL)
   20          CONTINUE
C
C              Expand Q by adding a unit row and column.
C
               IF (NFREE.GT.1) THEN
                  CALL F06FBF(NFREE-1,ZERO,ZY(NFREE,1),LDZY)
                  CALL F06FBF(NFREE-1,ZERO,ZY(1,NFREE),1)
               END IF
               ZY(NFREE,NFREE) = ONE
            END IF
         ELSE
C
C           Case 2.  A general constraint has been deleted.
C           =======
C
            IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
               WRITE (REC,FMT=99997) NACTIV, NZ, NFREE, KDEL, JDEL,
     *           UNITQ
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
            ITDEL = KDEL
            NACTIV = NACTIV - 1
C
C           Delete row  kdel  of T and move up the ones below it.
C           T becomes reverse lower Hessenberg.
C
            DO 40 I = KDEL, NACTIV
               KACTIV(I) = KACTIV(I+1)
               LD = NFREE - I
               CALL DCOPY(I+1,T(I+1,LD),LDT,T(I,LD),LDT)
   40       CONTINUE
         END IF
C
         NZ = NZ + 1
C
         IF (NACTIV.EQ.0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
C           ------------------------------------------------------------
C           Restore the NACTIV x (NACTIV+1) reverse-Hessenberg matrix  T
C           to reverse-triangular form.  The last NACTIV super-diagonal
C           elements are removed using a backward sweep of plane
C           rotations.  The rotation for the singleton in the first
C           column is generated separately.
C           ------------------------------------------------------------
            NSUP = NACTIV - ITDEL + 1
C
            IF (NSUP.GT.0) THEN
               NPIV = NFREE - ITDEL + 1
               IF (NSUP.GT.1) THEN
                  CALL DCOPY(NSUP-1,T(NACTIV-1,NZ+1),LDT-1,S(NZ+1),1)
                  CALL F06QZZ('Remove',NACTIV,1,NSUP,C(NZ+1),S(NZ+1),
     *                        T(1,NZ+1),LDT)
               END IF
C
               CALL F06BAF(T(NACTIV,NZ+1),T(NACTIV,NZ),CS,SN)
               T(NACTIV,NZ) = ZERO
               S(NZ) = -SN
               C(NZ) = CS
C
               CALL F06QXF('Right','Variable','Backwards',NFREE,NFREE,
     *                     NZ,NPIV,C,S,ZY,LDZY)
               CALL F06QXF('Left ','Variable','Backwards',NPIV,NGQ,NZ,
     *                     NPIV,C,S,GQ,N)
C
               NT = MIN(NRANK,NPIV)
C
               IF (NT.LT.NPIV .AND. NT.GT.0) THEN
C
C                 R is upper trapezoidal, pretend R is (NT x n) and
C                 apply the rotations in columns  max(NT,NZ)  thru NPIV.
C
                  CALL F06QXF('Right','Variable','Backwards',NT,N,
     *                        MAX(NT,NZ),NPIV,C,S,R,LDR)
               END IF
C
C              Apply the column transformations to the triangular part
C              of R.  A subdiagonal element is generated that must be
C              eliminated by a row rotation before the next column
C              transformation can be applied.
C
               IF (NZ.LT.NT) THEN
                  CALL F06QTF('Right',NT,NZ,NT,C,S,R,LDR)
               END IF
C
C              Apply the row rotations to the remaining rows of R.
C
               CALL F06QXF('Left','Variable','Backwards',NT,N-NT,NZ,NT,
     *                     C,S,R(1,MIN(NT+1,N)),LDR)
C
               IF (NRES.GT.0) CALL F06QXF('Left','Variable','Backwards',
     *                                    NT,NRES,NZ,NT,C,S,RES,N)
C
            END IF
            CALL F06FLF(NACTIV,T(NACTIV,NZ+1),LDT-1,DTMAX,DTMIN)
         END IF
      END IF
C
      NRZ1 = NRZ + 1
C
      IF (NZ.GT.NRZ) THEN
         IF (JDEL.GT.0) THEN
            JART = NRZ1 - 1 + IDAMAX(NZ-NRZ1+1,GQ(NRZ1,1),1)
         ELSE
            JART = -JDEL
         END IF
C
         IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) NZ, NRZ1, JART
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         IF (JART.GT.NRZ1) THEN
C
C           Swap columns NRZ1 and JART of R.
C
            IF (UNITQ) THEN
               K = KX(NRZ1)
               KX(NRZ1) = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP(NFREE,ZY(1,NRZ1),1,ZY(1,JART),1)
            END IF
C
            CALL DSWAP(NGQ,GQ(NRZ1,1),N,GQ(JART,1),N)
            IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,NRZ1,JART,R,
     *                                  RES,C,S)
         END IF
      END IF
C
      NRZ = NRZ1
C
      RETURN
C
C
C     End of  E04NCT (LSDEL).
C
99999 FORMAT (/' //E04NCT //  Artificial constraint deleted.      ',
     *       /' //E04NCT //      NZ   NRZ   JART                 ',
     *       /' //E04NCT //  ',3I6)
99998 FORMAT (/' //E04NCT //  Simple bound deleted.               ',
     *       /' //E04NCT //  NACTIV    NZ NFREE    IR  JDEL UNITQ',
     *       /' //E04NCT //  ',5I6,L6)
99997 FORMAT (/' //E04NCT //  General constraint deleted.         ',
     *       /' //E04NCT //  NACTIV    NZ NFREE  KDEL  JDEL UNITQ',
     *       /' //E04NCT //  ',5I6,L6)
      END
      SUBROUTINE E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,
     *                  LDAQP,LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,
     *                  ADX,BL,BU,RPQ,RPQ0,DX,GQ,R,T,ZY,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1082 (JUL 1993).
C
C     ******************************************************************
C     E04UCM   defines a point which lies on the initial working set for
C     the QP subproblem.  This routine is similar to E04NCH except
C     that advantage is taken of the fact that the initial estimate of
C     the solution of the least-squares subproblem is zero.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level 2 BLAS added 12-June-1986.
C     This version of E04UCM dated 11-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DXNORM, GDX
      INTEGER           LDAQP, LDR, LDT, LDZY, N, NACTIV, NCQP, NCTOTL,
     *                  NFREE, NLNX, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), AQP(LDAQP,*), BL(NCTOTL), BU(NCTOTL),
     *                  DX(N), GQ(N), R(LDR,*), RPQ(NLNX), RPQ0(NLNX),
     *                  T(LDT,*), WORK(N), ZY(LDZY,*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, J, K, NFIXED, NR
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRMV, E04NBT, E04NBW,
     *                  F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
C
      GDX = ZERO
      CALL F06FBF(N,ZERO,DX,1)
      CALL F06FBF(NLNX,ZERO,RPQ,1)
      CALL F06FBF(NLNX,ZERO,RPQ0,1)
C
      IF (NACTIV+NFIXED.GT.0) THEN
C
C        Set  work = residuals for constraints in the working set.
C        Solve for  dx,  the smallest correction to  x  that gives a
C        point on the constraints in the working set.
C        Set the fixed variables on their bounds,  solve the triangular
C        system  T*(dxy) = residuals,  and define  dx = Y*(dxy).
C        Use  (dxy)  to update  d(=Pr)  as  d = d - R'(  0  ).
C                                                     ( dxy )
C
         DO 20 I = 1, NFIXED
            J = KX(NFREE+I)
            IF (ISTATE(J).LE.3) THEN
               BND = BL(J)
               IF (ISTATE(J).EQ.2) BND = BU(J)
               DX(J) = BND
               WORK(NFREE+I) = BND
            ELSE
               WORK(NFREE+I) = ZERO
            END IF
   20    CONTINUE
C
         DO 40 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(NZ+I) = BND - DDOT(N,AQP(K,1),LDAQP,DX,1)
   40    CONTINUE
C
         IF (NACTIV.GT.0) CALL E04NBT(1,LDT,NACTIV,T(1,NZ+1),WORK(NZ+1))
         CALL DCOPY(NACTIV+NFIXED,WORK(NZ+1),1,DX(NZ+1),1)
         IF (NZ.GT.0) CALL F06FBF(NZ,ZERO,DX,1)
C
         GDX = DDOT(NACTIV+NFIXED,GQ(NZ+1),1,DX(NZ+1),1)
C
         IF (NZ.LT.N) THEN
            CALL DGEMV('N',NZ,N-NZ,-ONE,R(1,NZ+1),LDR,DX(NZ+1),1,ONE,
     *                 RPQ,1)
            IF (NZ.LT.NLNX) THEN
               NR = LDR
               IF (NZ+1.EQ.N) NR = 1
               CALL DCOPY(NLNX-NZ,DX(NZ+1),1,RPQ(NZ+1),1)
               CALL DSCAL(NLNX-NZ,(-ONE),RPQ(NZ+1),1)
               CALL DTRMV('U','N','N',NLNX-NZ,R(NZ+1,NZ+1),NR,RPQ(NZ+1),
     *                    1)
               IF (NLNX.LT.N) THEN
                  NR = LDR
                  IF (NLNX+1.EQ.N) NR = N - NZ
                  CALL DGEMV('N',NLNX-NZ,N-NLNX,-ONE,R(NZ+1,NLNX+1),NR,
     *                       DX(NLNX+1),1,ONE,RPQ(NZ+1),1)
               END IF
            END IF
         END IF
C
         CALL E04NBW(2,N,NZ,NFREE,LDZY,UNITQ,KX,DX,ZY,WORK)
      END IF
C
C     ------------------------------------------------------------------
C     Compute the 2-norm of  DX.
C     Initialize  A*DX.
C     ------------------------------------------------------------------
      DXNORM = DNRM2(N,DX,1)
      IF (NCQP.GT.0) CALL DGEMV('N',NCQP,N,ONE,AQP,LDAQP,DX,1,ZERO,ADX,
     *                          1)
C
      IF (NPDBG .AND. INPDBG(2).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,2,REC)
         DO 60 I = 1, N, 5
            WRITE (REC,FMT=99998) (DX(J),J=I,MIN(I+4,N))
            CALL X04BAF(IPRINT,REC(1))
   60    CONTINUE
      END IF
C
      RETURN
C
C
C     End of  E04UCM. (NPSETX)
C
99999 FORMAT (/' //E04UCM// Variables after E04UCM ... ')
99998 FORMAT (5G12.3)
      END
      SUBROUTINE E04NCK(PRBTYP,MSGLVL,N,NACTIV,NFREE,LDA,LDT,NUMINF,NZ,
     *                  NRZ,ISTATE,KACTIV,KX,DINKY,JSMLST,KSMLST,JINF,
     *                  JTINY,JBIGST,KBIGST,TRULAM,A,ANORMS,GQ,RLAMDA,T,
     *                  WTINF)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1063 (JUL 1993).
C
C     ******************************************************************
C     E04NCK  first computes the Lagrange multiplier estimates for the
C     given working set.  It then determines the values and indices of
C     certain significant multipliers.  In this process, the multipliers
C     for inequalities at their upper bounds are adjusted so that a
C     negative multiplier for an inequality constraint indicates non-
C     optimality.  All adjusted multipliers are scaled by the 2-norm
C     of the associated constraint row.  In the following, the term
C     minimum refers to the ordering of numbers on the real line,  and
C     not to their magnitude.
C
C     JSMLST  is the index of the minimum of the set of adjusted
C             multipliers with values less than  - DINKY.  A negative
C             JSMLST defines the index in Q'g of the artificial
C             constraint to be deleted.
C     KSMLST  marks the position of general constraint JSMLST in KACTIV.
C
C     JBIGST  is the index of the largest of the set of adjusted
C             multipliers with values greater than (1 + DINKY).
C     KBIGST  marks its position in KACTIV.
C
C     On exit,  elements 1 thru NACTIV of RLAMDA contain the unadjusted
C     multipliers for the general constraints.  Elements NACTIV onwards
C     of RLAMDA contain the unadjusted multipliers for the bounds.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCK dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DINKY, TRULAM
      INTEGER           JBIGST, JINF, JSMLST, JTINY, KBIGST, KSMLST,
     *                  LDA, LDT, MSGLVL, N, NACTIV, NFREE, NRZ, NUMINF,
     *                  NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ANORMS(*), GQ(N), RLAMDA(N), T(LDT,*),
     *                  WTINF(*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORMJ, BIGGST, BLAM, RLAM, SCDLAM, SMLLST,
     *                  TINYLM
      INTEGER           I, IS, J, K, L, NFIXED
C     .. Local Arrays ..
      CHARACTER*80      REC(80)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04NBT, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Executable Statements ..
      NFIXED = N - NFREE
C
      JSMLST = 0
      KSMLST = 0
      SMLLST = -DINKY
C
      TINYLM = DINKY
      JTINY = 0
C
      JBIGST = 0
      KBIGST = 0
      BIGGST = ONE + DINKY
C
      IF (NRZ.LT.NZ) THEN
C        ---------------------------------------------------------------
C        Compute JSMLST for the artificial constraints.
C        ---------------------------------------------------------------
         DO 20 J = NRZ + 1, NZ
            RLAM = -ABS(GQ(J))
            IF (RLAM.LT.SMLLST) THEN
               SMLLST = RLAM
               JSMLST = -J
            ELSE IF (RLAM.LT.TINYLM) THEN
               TINYLM = RLAM
               JTINY = J
            END IF
   20    CONTINUE
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99999)
               CALL X04BAY(ISUMM,2,REC)
               DO 40 J = NRZ + 1, NZ, 4
                  WRITE (REC,FMT=99993) (GQ(K),K=J,MIN(J+3,NZ))
                  CALL X04BAF(ISUMM,REC(1))
   40          CONTINUE
            END IF
         END IF
C
      END IF
C
C     ---------------------------------------------------------------
C     Compute JSMLST for regular constraints and temporary bounds.
C     ---------------------------------------------------------------
C     First, compute the Lagrange multipliers for the general
C     constraints in the working set, by solving  T'*lamda = Y'g.
C
      IF (N.GT.NZ) CALL DCOPY(N-NZ,GQ(NZ+1),1,RLAMDA,1)
      IF (NACTIV.GT.0) CALL E04NBT(2,LDT,NACTIV,T(1,NZ+1),RLAMDA)
C
C     --------------------------------------------------------------
C     Now set elements NACTIV, NACTIV+1,... of  RLAMDA  equal to
C     the multipliers for the bound constraints.
C     --------------------------------------------------------------
      DO 80 L = 1, NFIXED
         J = KX(NFREE+L)
         BLAM = RLAMDA(NACTIV+L)
         DO 60 K = 1, NACTIV
            I = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(K)
   60    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
   80 CONTINUE
C
C     --------------------------------------------------------------
C     Find JSMLST and KSMLST.
C     --------------------------------------------------------------
      DO 100 K = 1, N - NZ
         IF (K.GT.NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(K) + N
         END IF
C
         IS = ISTATE(J)
C
         I = J - N
         IF (J.LE.N) ANORMJ = ONE
         IF (J.GT.N) ANORMJ = ANORMS(I)
C
         RLAM = RLAMDA(K)
C
C        Change the sign of the estimate if the constraint is in
C        the working set at its upper bound.
C
         IF (IS.EQ.2) RLAM = -RLAM
         IF (IS.EQ.3) RLAM = ABS(RLAM)
         IF (IS.EQ.4) RLAM = -ABS(RLAM)
C
         IF (IS.NE.3) THEN
            SCDLAM = RLAM*ANORMJ
            IF (SCDLAM.LT.SMLLST) THEN
               SMLLST = SCDLAM
               JSMLST = J
               KSMLST = K
            ELSE IF (SCDLAM.LT.TINYLM) THEN
               TINYLM = SCDLAM
               JTINY = J
            END IF
         END IF
C
         IF (NUMINF.GT.0 .AND. J.GT.JINF) THEN
            SCDLAM = RLAM/WTINF(J)
            IF (SCDLAM.GT.BIGGST) THEN
               BIGGST = SCDLAM
               TRULAM = RLAMDA(K)
               JBIGST = J
               KBIGST = K
            END IF
         END IF
  100 CONTINUE
C
C     --------------------------------------------------------------
C     If required, print the multipliers.
C     --------------------------------------------------------------
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            IF (NFIXED.GT.0) THEN
               WRITE (REC,FMT=99998) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 120 J = 1, NFIXED, 4
                  WRITE (REC,FMT=99992) (KX(NFREE+K),RLAMDA(NACTIV+K),
     *              K=J,MIN(J+3,NFIXED))
                  CALL X04BAF(ISUMM,REC(1))
  120          CONTINUE
            END IF
            IF (NACTIV.GT.0) THEN
               WRITE (REC,FMT=99997) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 140 J = 1, NACTIV, 4
                  WRITE (REC,FMT=99992) (KACTIV(K),RLAMDA(K),K=J,
     *              MIN(J+3,NACTIV))
                  CALL X04BAF(ISUMM,REC(1))
  140          CONTINUE
            END IF
         END IF
      END IF
C
      IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
         WRITE (REC,FMT=99996) JSMLST, SMLLST, KSMLST
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99995) JBIGST, BIGGST, KBIGST
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99994) JTINY, TINYLM
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      RETURN
C
C
C     End of  E04NCK. (LSMULS)
C
99999 FORMAT (/' Multipliers for the artificial constraints        ')
99998 FORMAT (/' Multipliers for the ',A2,' bound  constraints   ')
99997 FORMAT (/' Multipliers for the ',A2,' linear constraints   ')
99996 FORMAT (/' //E04NCK//  JSMLST     SMLLST     KSMLST (Scaled) ',
     *       /' //E04NCK//  ',I6,1P,D11.2,5X,I6)
99995 FORMAT (' //E04NCK//  JBIGST     BIGGST     KBIGST (Scaled) ',
     *       /' //E04NCK//  ',I6,1P,D11.2,5X,I6)
99994 FORMAT (' //E04NCK//   JTINY     TINYLM                     ',
     *       /' //E04NCK//  ',I6,1P,D11.2)
99993 FORMAT (4(5X,1P,D11.2))
99992 FORMAT (4(I5,1P,D11.2))
      END
      SUBROUTINE E04MFQ(N,NCLIN,ISTATE,BIGBND,NVIOL,JMAX,ERRMAX,AX,BL,
     *                  BU,FEATOL,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFQ  checks the residuals of the constraints that are believed
C     to be feasible.  The number of constraints violated by more than
C     featol is computed, along with the maximum constraint violation.
C
C     Original version written by PEG,   April    1984.
C     This version of  E04MFQ  dated  30-Jun-1988.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, ERRMAX
      INTEGER           JMAX, N, NCLIN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES
      INTEGER           IS, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Compute the number of constraints (NVIOL) violated by more than
C     FEATOL and  the maximum constraint violation (ERRMAX).
C     (The residual of a constraint in the working set is treated as if
C     it were an equality constraint fixed at that bound.)
C     ==================================================================
      NVIOL = 0
      JMAX = 0
      ERRMAX = ZERO
C
      DO 40 J = 1, N + NCLIN
         IS = ISTATE(J)
C
         IF (IS.GE.0) THEN
            FEASJ = FEATOL(J)
C
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               CON = AX(J-N)
            END IF
C
C           Check for constraint violations.
C
            IF (BL(J).GT.BIGLOW) THEN
               RES = BL(J) - CON
               IF (RES.GT.FEASJ) THEN
                  NVIOL = NVIOL + 1
                  GO TO 20
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               RES = BU(J) - CON
               IF (RES.LT.(-FEASJ)) THEN
                  NVIOL = NVIOL + 1
                  RES = -RES
                  GO TO 20
               END IF
            END IF
C
C           this constraint is satisfied,  but count a large residual
C           as a violation if the constraint is in the working set.
C
            RES = ZERO
C
            IF (IS.EQ.1) THEN
               RES = ABS(BL(J)-CON)
C
            ELSE IF (IS.EQ.2) THEN
               RES = ABS(BU(J)-CON)
C
            ELSE IF (IS.EQ.3) THEN
               RES = ABS(BU(J)-CON)
            END IF
C
            IF (RES.GT.FEASJ) NVIOL = NVIOL + 1
C
   20       IF (RES.GT.ERRMAX) THEN
               JMAX = J
               ERRMAX = RES
            END IF
         END IF
   40 CONTINUE
C
      IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) ERRMAX, JMAX
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      RETURN
C
C     End of  E04MFQ.  (CMFEAS)
C
99999 FORMAT (/' //E04MFQ//  the maximum violation is ',1P,D14.2,' in ',
     *       'constraint',I5)
      END
      SUBROUTINE E04UDV(STRING,NUMIN,NUMOUT,LIST)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Description and usage:
C
C       An aid to parsing input data.  The individual 'tokens' in a
C     character string are isolated, converted to uppercase, and stored
C     in an array.  Here, a token is a group of significant, contiguous
C     characters.  The following are NON-significant, and hence may
C     serve as separators:  blanks, horizontal tabs, commas, colons,
C     and equal signs.  See E04UDW for details.  Processing continues
C     until the requested number of tokens have been found or the end
C     of the input string is reached.
C
C
C     Parameters:
C
C     Name    Dimension  Type  I/O/S  Description
C     STRING              C    I      Input string to be analyzed.
C     NUMIN               I    I/O    Number of tokens requested (input)
C     NUMOUT                          and found (output).
C     (NUMIN and NUMOUT were both called NUMBER in the original)
C
C     LIST    NUMIN       C      O    Array of tokens, changed to upper
C                                    case.
C
C
C     External references:
C
C     Name    Description
C     E04UDW  Finds positions of first and last significant characters.
C     E04UDU  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               Appears to satisfy the ANSI Fortran 77 standard.
C
C
C     Notes:
C
C     (1)  IMPLICIT NONE is non-standard.  (Has been commented out.)
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C     16 Jan. 1984    RAK    Initial design and coding.
C     16 Mar. 1984    RAK    Revised header to reflect full list of
C                            separators, repaired faulty WHILE clause
C                            in '10' loop.
C     18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                            at a time, leaving STRING unchanged.
C
C-----------------------------------------------------------------------
C
C     .. Parameters ..
      CHARACTER         BLANK
      PARAMETER         (BLANK=' ')
C     .. Scalar Arguments ..
      INTEGER           NUMIN, NUMOUT
      CHARACTER*(*)     STRING
C     .. Array Arguments ..
      CHARACTER*(*)     LIST(NUMIN)
C     .. Local Scalars ..
      INTEGER           COUNT, FIRST, I, LAST, MARK
C     .. External Subroutines ..
      EXTERNAL          E04UDU, E04UDW
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      FIRST = 1
      LAST = LEN(STRING)
C
      COUNT = 0
   20 CONTINUE
C
C        Get delimiting indices of next token, if any.
C
      CALL E04UDW(STRING,FIRST,LAST,MARK)
      IF (LAST.GT.0) THEN
         COUNT = COUNT + 1
C
C           Pass token to output string array, then change case.
C
         LIST(COUNT) = STRING(FIRST:MARK)
c        CALL E04UDU(LIST(COUNT))
         FIRST = MARK + 2
         IF (COUNT.LT.NUMIN) GO TO 20
C
      END IF
C
C
C     Fill the rest of LIST with blanks and set NUMBER for output.
C
      DO 40 I = COUNT + 1, NUMIN
         LIST(I) = BLANK
   40 CONTINUE
C
      NUMOUT = COUNT
C
C
C     Termination.
C     ------------
C
      RETURN
C
C     End of  E04UDV. (OPTOKN)
      END
      SUBROUTINE E04UCK(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,NOUT,
     *                  ALFMAX,EPSAF,G0,TARGTG,FTRY,GTRY,TOLABS,TOLREL,
     *                  TOLTNY,ALFA,ALFBST,FBEST,GBEST)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1080 (JUL 1993).
C
C     ==================================================================
C     E04UCK  finds a sequence of improving estimates of a minimizer of
C     the univariate function f(alpha) in the interval (0,ALFMAX].
C     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
C     E04UCK  requires both  f(alpha)  and  f'(alpha) to be evaluated at
C     points in the interval.  Estimates of the minimizer are computed
C     using safeguarded cubic interpolation.
C
C     Reverse communication is used to allow the calling program to
C     evaluate f and f'.  Some of the parameters must be set or tested
C     by the calling program.  The remainder would ordinarily be local
C     variables.
C
C     Input parameters (relevant to the calling program)
C     --------------------------------------------------
C
C     FIRST         must be .TRUE. on the first entry.
C                   It is subsequently altered by E04UCK.
C
C     DEBUG         specifies whether detailed output is wanted.
C
C     MAXF          is an upper limit on the number of times E04UCK is
C                   to be entered consecutively with DONE = .FALSE.
C                   (following an initial entry with FIRST = .TRUE.).
C
C     ALFA          is the first estimate of a minimizer.  ALFA is
C                   subsequently altered by E04UCK (see below).
C
C     ALFMAX        is the upper limit of the interval to be searched.
C
C     EPSAF         is an estimate of the absolute precision in the
C                   computed value of f(0).
C
C     FTRY, GTRY    are the values of f, f'  at the new point
C                   ALFA = ALFBST + XTRY.
C
C     G0            is the value of f'(0).  G0 must be negative.
C
C     TOLABS,TOLREL define a function TOL(ALFA) = TOLREL*ALFA + TOLABS
C                   such that if f has already been evaluated at ALFA,
C                   it will not be evaluated closer than TOL(ALFA).
C                   These values may be reduced by E04UCK.
C
C     TARGTG        is the target value of abs(f'(ALFA)). The search
C                   is terminated when
C                    abs(f'(ALFA)) le TARGTG and f(ALFA) lt 0.
C
C     TOLTNY        is the smallest value that TOLABS is allowed to be
C                   reduced to.
C
C     Output parameters (relevant to the calling program)
C     ---------------------------------------------------
C
C     IMPRVD        is .TRUE. if the previous ALFA was the best point so
C                   far.  Any related quantities should be saved by the
C                   calling program (e.g., gradient arrays) before
C                   paying attention to the variable DONE.
C
C     DONE = .FALSE.  means the calling program should evaluate
C                      FTRY = f(ALFA),  GTRY = f'(ALFA)
C                   for the new trial ALFA, and re-enter E04UCK.
C
C     DONE = .TRUE.   means that no new ALFA was calculated.  The value
C                   of INFORM gives the result of the search as follows
C
C                   INFORM = 1 means the search has terminated
C                              successfully with ALFBST < ALFMAX.
C
C                   INFORM = 2 means the search has terminated
C                              successfully with ALFBST = ALFMAX.
C
C                   INFORM = 3 means that the search failed to find a
C                              point of sufficient decrease in MAXF
C                              functions, but a lower point was found.
C
C                   INFORM = 4 means ALFMAX is so small that a search
C                              should not have been attempted.
C
C                   INFORM = 5 is never set by E04UCK.
C
C                   INFORM = 6 means the search has failed to find a
C                              useful step.  The interval of uncertainty
C                              is [0,B] with B < 2*TOLABS. A minimizer
C                              lies very close to ALFA = 0, or f'(0) is
C                              not sufficiently accurate.
C
C                   INFORM = 7 if no better point could be found after
C                              MAXF  function calls.
C
C                   INFORM = 8 means the input parameters were bad.
C                              ALFMAX le TOLTNY  or G0 ge zero.
C                              No function evaluations were made.
C
C     NUMF          counts the number of times E04UCK has been entered
C                   consecutively with DONE = .FALSE. (i.e., with a new
C                   function value FTRY).
C
C     ALFA          is the point at which the next function FTRY and
C                   derivative GTRY must be computed.
C
C     ALFBST        should be accepted by the calling program as the
C                   approximate minimizer, whenever E04UCK returns
C                   INFORM = 1 or 2 (and possibly 3).
C
C     FBEST, GBEST  will be the corresponding values of f, f'.
C
C
C     The following parameters retain information between entries
C     -----------------------------------------------------------
C
C     BRAKTD        is .FALSE. if f and f' have not been evaluated at
C                   the far end of the interval of uncertainty.  In this
C                   case, the point B will be at ALFMAX + TOL(ALFMAX).
C
C     CRAMPD        is .TRUE. if ALFMAX is very small (le TOLABS).  If
C                   the search fails, this indicates that a zero step
C                   should be taken.
C
C     EXTRAP        is .TRUE. if XW lies outside the interval of
C                   uncertainty.  In this case, extra safeguards are
C                   applied to allow for instability in the polynomial
C                   fit.
C
C     MOVED         is .TRUE. if a better point has been found, i.e.,
C                   ALFBST gt 0.
C
C     WSET          records whether a second-best point has been
C                   determined it will always be .TRUE. when convergence
C                   is tested.
C
C     NSAMEA        is the number of consecutive times that the
C                   left-hand end point of the interval of uncertainty
C                   has remained the same.
C
C     NSAMEB        similarly for the right-hand end.
C
C     A, B, ALFBST  define the current interval of uncertainty.
C                   A minimizer lies somewhere in the interval
C                   [ALFBST + A, ALFBST + B].
C
C     ALFBST        is the best point so far.  It is always at one end
C                   of the interval of uncertainty.  hence we have
C                   either  A lt 0,  B = 0  or  A = 0,  B gt 0.
C
C     FBEST, GBEST  are the values of f, f' at the point ALFBST.
C
C     FACTOR        controls the rate at which extrapolated estimates
C                   of ALFA may expand into the interval of uncertainty.
C                   FACTOR is not used if a minimizer has been bracketed
C                   (i.e., when the variable BRAKTD is .TRUE.).
C
C     FW, GW        are the values of f, f' at the point ALFBST + XW.
C                   they are not defined until WSET is .TRUE..
C
C     XTRY          is the trial point in the shifted interval (A, B).
C
C     XW            is such that  ALFBST + XW  is the second-best point.
C                   it is not defined until  WSET  is .TRUE..
C                   in some cases,  XW  will replace a previous  XW
C                   that has a lower function but has just been excluded
C                   from the interval of uncertainty.
C
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version February 1982.  Rev. May 1983.
C     Original f77 version 22-August-1985.
C     This version of E04UCK dated  14-Sep-92.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT1, HALF
      PARAMETER         (ZERO=0.0D+0,POINT1=0.1D+0,HALF=0.5D+0)
      DOUBLE PRECISION  ONE, THREE, FIVE
      PARAMETER         (ONE=1.0D+0,THREE=3.0D+0,FIVE=5.0D+0)
      DOUBLE PRECISION  TEN, ELEVEN
      PARAMETER         (TEN=1.0D+1,ELEVEN=1.1D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBST, ALFMAX, EPSAF, FBEST, FTRY, G0,
     *                  GBEST, GTRY, TARGTG, TOLABS, TOLREL, TOLTNY
      INTEGER           INFORM, MAXF, NOUT, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ABSR, ARTIFA, ARTIFB, B, DAUX, DTRY, FACTOR,
     *                  FW, GW, Q, R, S, SCALE, TOL, TOLMAX, TRUEA,
     *                  TRUEB, XMIDPT, XTRY, XW
      INTEGER           NSAMEA, NSAMEB
      LOGICAL           BADFUN, BRAKTD, CLOSEF, CRAMPD, EXTRAP, FITOK,
     *                  FOUND, MOVED, QUITF, QUITI, SETXW, WSET
C     .. Local Arrays ..
      CHARACTER*120     REC(7)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Save statement ..
      SAVE              BRAKTD, CRAMPD, EXTRAP, MOVED, WSET, NSAMEA,
     *                  NSAMEB, A, B, FACTOR, XTRY, XW, FW, GW, TOLMAX
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Local variables
C     ===============
C
C     CLOSEF     is .TRUE. if the new function FTRY is within EPSAF of
C                FBEST (up or down).
C
C     FOUND      is .TRUE. if the sufficient decrease conditions hold at
C                ALFBST.
C
C     QUITF      is .TRUE. when  MAXF  function calls have been made.
C
C     QUITI      is .TRUE. when the interval of uncertainty is less than
C                2*TOL.
C     ------------------------------------------------------------------
C
      BADFUN = .FALSE.
      QUITF = .FALSE.
      QUITI = .FALSE.
      IMPRVD = .FALSE.
C
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First entry.  Initialize various quantities, check input data
C        and prepare to evaluate the function at the initial ALFA.
C        ---------------------------------------------------------------
         FIRST = .FALSE.
         NUMF = 0
         ALFBST = ZERO
         BADFUN = ALFMAX .LE. TOLTNY .OR. G0 .GE. ZERO
         DONE = BADFUN
         MOVED = .FALSE.
C
         IF ( .NOT. DONE) THEN
            BRAKTD = .FALSE.
            CRAMPD = ALFMAX .LE. TOLABS
            EXTRAP = .FALSE.
            WSET = .FALSE.
            NSAMEA = 0
            NSAMEB = 0
C
            TOLMAX = TOLABS + TOLREL*ALFMAX
            A = ZERO
            B = ALFMAX + TOLMAX
            FACTOR = FIVE
            TOL = TOLABS
            XTRY = ALFA
            IF (DEBUG) THEN
               WRITE (REC,FMT=99999) G0, TOLABS, ALFMAX, TARGTG, TOLREL,
     *           EPSAF, CRAMPD
               CALL X04BAY(NOUT,4,REC)
            END IF
         END IF
      ELSE
C        ---------------------------------------------------------------
C        Subsequent entries. The function has just been evaluated at
C        ALFA = ALFBST + XTRY,  giving FTRY and GTRY.
C        ---------------------------------------------------------------
         IF (DEBUG) THEN
            WRITE (REC,FMT=99998) ALFA, FTRY, GTRY
            CALL X04BAY(NOUT,2,REC)
         END IF
C
         NUMF = NUMF + 1
         NSAMEA = NSAMEA + 1
         NSAMEB = NSAMEB + 1
C
         IF ( .NOT. BRAKTD) THEN
            TOLMAX = TOLABS + TOLREL*ALFMAX
            B = ALFMAX - ALFBST + TOLMAX
         END IF
C
C        See if the new step is better.  If ALFA is large enough that
C        FTRY can be distinguished numerically from zero,  the function
C        is required to be sufficiently negative.
C
         CLOSEF = ABS(FTRY-FBEST) .LE. EPSAF
         IF (CLOSEF) THEN
            IMPRVD = ABS(GTRY) .LE. ABS(GBEST)
         ELSE
            IMPRVD = FTRY .LT. FBEST
         END IF
C
         IF (IMPRVD) THEN
C
C           We seem to have an improvement.  The new point becomes the
C           origin and other points are shifted accordingly.
C
            FW = FBEST
            FBEST = FTRY
            GW = GBEST
            GBEST = GTRY
            ALFBST = ALFA
            MOVED = .TRUE.
C
            A = A - XTRY
            B = B - XTRY
            XW = ZERO - XTRY
            WSET = .TRUE.
            EXTRAP = XW .LT. ZERO .AND. GBEST .LT. ZERO .OR. XW .GT.
     *               ZERO .AND. GBEST .GT. ZERO
C
C           Decrease the length of the interval of uncertainty.
C
            IF (GTRY.LE.ZERO) THEN
               A = ZERO
               NSAMEA = 0
            ELSE
               B = ZERO
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
         ELSE
C
C           The new function value is not better than the best point so
C           far.  The origin remains unchanged but the new point may
C           qualify as XW.  XTRY must be a new bound on the best point.
C
            IF (XTRY.LE.ZERO) THEN
               A = XTRY
               NSAMEA = 0
            ELSE
               B = XTRY
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
C
C           If XW has not been set or FTRY is better than FW, update the
C           points accordingly.
C
            IF (WSET) THEN
               SETXW = FTRY .LT. FW .OR. .NOT. EXTRAP
            ELSE
               SETXW = .TRUE.
            END IF
C
            IF (SETXW) THEN
               XW = XTRY
               FW = FTRY
               GW = GTRY
               WSET = .TRUE.
               EXTRAP = .FALSE.
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Check the termination criteria.  WSET will always be .TRUE..
C        ---------------------------------------------------------------
         TOL = TOLABS + TOLREL*ALFBST
         TRUEA = ALFBST + A
         TRUEB = ALFBST + B
C
         FOUND = ABS(GBEST) .LE. TARGTG
         QUITF = NUMF .GE. MAXF
         QUITI = B - A .LE. TOL + TOL
C
         IF (QUITI .AND. .NOT. MOVED) THEN
C
C           The interval of uncertainty appears to be small enough,
C           but no better point has been found.  Check that changing
C           ALFA by B-A changes f by less than EPSAF.
C
            TOL = TOL/TEN
            TOLABS = TOL
            QUITI = ABS(FW) .LE. EPSAF .OR. TOL .LE. TOLTNY
         END IF
C
         DONE = QUITF .OR. QUITI .OR. FOUND
C
         IF (DEBUG) THEN
            WRITE (REC,FMT=99997) TRUEA, TRUEB, B - A, TOL, NSAMEA,
     *        NSAMEB, NUMF, BRAKTD, EXTRAP, CLOSEF, IMPRVD, FOUND,
     *        QUITI, ALFBST, FBEST, GBEST, ALFBST + XW, FW, GW
            CALL X04BAY(NOUT,7,REC)
         END IF
C
C        ---------------------------------------------------------------
C        Proceed with the computation of a trial steplength.
C        The choices are...
C        1. Parabolic fit using derivatives only, if the f values are
C           close.
C        2. Cubic fit for a minimizer, using both f and f'.
C        3. Damped cubic or parabolic fit if the regular fit appears to
C           be consistently overestimating the distance to a minimizer.
C        4. Bisection, geometric bisection, or a step of  TOL  if
C           choices 2 or 3 are unsatisfactory.
C        ---------------------------------------------------------------
         IF ( .NOT. DONE) THEN
            XMIDPT = HALF*(A+B)
            S = ZERO
            Q = ZERO
C
            IF (CLOSEF) THEN
C              ---------------------------------------------------------
C              Fit a parabola to the two best gradient values.
C              ---------------------------------------------------------
               S = GBEST
               Q = GBEST - GW
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99995)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            ELSE
C              ---------------------------------------------------------
C              Fit cubic through  FBEST  and  FW.
C              ---------------------------------------------------------
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99996)
                  CALL X04BAF(NOUT,REC(1))
               END IF
               FITOK = .TRUE.
               R = THREE*(FBEST-FW)/XW + GBEST + GW
               ABSR = ABS(R)
               S = SQRT(ABS(GBEST))*SQRT(ABS(GW))
C
C              Compute  Q =  the square root of  R*R - GBEST*GW.
C              The method avoids unnecessary underflow and overflow.
C
               IF ((GW.LT.ZERO .AND. GBEST.GT.ZERO)
     *             .OR. (GW.GT.ZERO .AND. GBEST.LT.ZERO)) THEN
                  SCALE = ABSR + S
                  IF (SCALE.EQ.ZERO) THEN
                     Q = ZERO
                  ELSE
                     Q = SCALE*SQRT((ABSR/SCALE)**2+(S/SCALE)**2)
                  END IF
               ELSE IF (ABSR.GE.S) THEN
                  Q = SQRT(ABSR+S)*SQRT(ABSR-S)
               ELSE
                  FITOK = .FALSE.
               END IF
C
               IF (FITOK) THEN
C
C                 Compute a minimizer of the fitted cubic.
C
                  IF (XW.LT.ZERO) Q = -Q
                  S = GBEST - R - Q
                  Q = GBEST - GW - Q - Q
               END IF
            END IF
C           ------------------------------------------------------------
C           Construct an artificial interval  (ARTIFA, ARTIFB)  in which
C           the new estimate of a minimizer must lie.  Set a default
C           value of XTRY that will be used if the polynomial fit fails.
C           ------------------------------------------------------------
            ARTIFA = A
            ARTIFB = B
            IF ( .NOT. BRAKTD) THEN
C
C              A minimizer has not been bracketed.  Set an artificial
C              upper bound by expanding the interval  XW  by a suitable
C              FACTOR.
C
               XTRY = -FACTOR*XW
               ARTIFB = XTRY
               IF (ALFBST+XTRY.LT.ALFMAX) FACTOR = FIVE*FACTOR
C
            ELSE IF (EXTRAP) THEN
C
C              The points are configured for an extrapolation.
C              Set a default value of  XTRY  in the interval  (A, B)
C              that will be used if the polynomial fit is rejected.  In
C              the following,  DTRY  and  DAUX  denote the lengths of
C              the intervals  (A, B)  and  (0, XW)  (or  (XW, 0),  if
C              appropriate).  The value of  XTRY is the point at which
C              the exponents of  DTRY  and  DAUX  are approximately
C              bisected.
C
               DAUX = ABS(XW)
               DTRY = B - A
               IF (DAUX.GE.DTRY) THEN
                  XTRY = FIVE*DTRY*(POINT1+DTRY/DAUX)/ELEVEN
               ELSE
                  XTRY = HALF*SQRT(DAUX)*SQRT(DTRY)
               END IF
               IF (XW.GT.ZERO) XTRY = -XTRY
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99993) XTRY, DAUX, DTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
C
C              Reset the artificial bounds.  If the point computed by
C              extrapolation is rejected,  XTRY will remain at the
C              relevant artificial bound.
C
               IF (XTRY.LE.ZERO) ARTIFA = XTRY
               IF (XTRY.GT.ZERO) ARTIFB = XTRY
            ELSE
C
C              The points are configured for an interpolation.  The
C              default value XTRY bisects the interval of uncertainty.
C              the artificial interval is just (A, B).
C
               XTRY = XMIDPT
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99994) XTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
               IF (NSAMEA.GE.3 .OR. NSAMEB.GE.3) THEN
C
C                 If the interpolation appears to be overestimating the
C                 distance to a minimizer,  damp the interpolation.
C
                  FACTOR = FACTOR/FIVE
                  S = FACTOR*S
               ELSE
                  FACTOR = ONE
               END IF
            END IF
C           ------------------------------------------------------------
C           The polynomial fits give  (S/Q)*XW  as the new step.
C           Reject this step if it lies outside  (ARTIFA, ARTIFB).
C           ------------------------------------------------------------
            IF (Q.NE.ZERO) THEN
               IF (Q.LT.ZERO) S = -S
               IF (Q.LT.ZERO) Q = -Q
               IF (S*XW.GE.Q*ARTIFA .AND. S*XW.LE.Q*ARTIFB) THEN
C
C                 Accept the polynomial fit.
C
                  IF (ABS(S*XW).GE.Q*TOL) THEN
                     XTRY = (S/Q)*XW
                  ELSE
                     XTRY = ZERO
                  END IF
                  IF (DEBUG) THEN
                     WRITE (REC,FMT=99992) XTRY
                     CALL X04BAF(NOUT,REC(1))
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     ==================================================================
C
      IF ( .NOT. DONE) THEN
         ALFA = ALFBST + XTRY
         IF (BRAKTD .OR. ALFA.LT.ALFMAX-TOLMAX) THEN
C
C           The function must not be evaluated too close to A or B.
C           (It has already been evaluated at both those points.)
C
            IF (XTRY.LE.A+TOL .OR. XTRY.GE.B-TOL) THEN
               IF (HALF*(A+B).LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
               ALFA = ALFBST + XTRY
            END IF
         ELSE
C
C           The step is close to, or larger than ALFMAX, replace it by
C           ALFMAX to force evaluation of  f  at the boundary.
C
            BRAKTD = .TRUE.
            XTRY = ALFMAX - ALFBST
            ALFA = ALFMAX
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Exit.
C     ------------------------------------------------------------------
      IF (DONE) THEN
         IF (BADFUN) THEN
            INFORM = 8
         ELSE IF (FOUND) THEN
            IF (ALFBST.LT.ALFMAX) THEN
               INFORM = 1
            ELSE
               INFORM = 2
            END IF
         ELSE IF (MOVED) THEN
            INFORM = 3
         ELSE IF (QUITF) THEN
            INFORM = 7
         ELSE IF (CRAMPD) THEN
            INFORM = 4
         ELSE
            INFORM = 6
         END IF
      END IF
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAY(NOUT,2,REC)
      END IF
      RETURN
C
C
C     End of E04UCK. (SRCHC)
C
99999 FORMAT (/'     G0  TOLABS  ALFMAX        ',1P,2D22.14,D16.8,/' T',
     *       'ARGTG  TOLREL   EPSAF        ',1P,2D22.14,D16.8,/' CRAMP',
     *       'D                        ',L3)
99998 FORMAT (/' ALFA    FTRY    GTRY          ',1P,2D22.14,D16.8)
99997 FORMAT (/' A       B       B - A   TOL   ',1P,2D22.14,2D16.8,
     *       /' NSAMEA  NSAMEB  NUMF          ',3I3,/' BRAKTD  EXTRAP ',
     *       ' CLOSEF  IMPRVD',4L3,/' FOUND   QUITI                 ',
     *       2L3,/' ALFBST  FBEST   GBEST         ',1P,3D22.14,/' ALFA',
     *       'W   FW      GW            ',1P,3D22.14)
99996 FORMAT (' Cubic.   ')
99995 FORMAT (' Parabola.')
99994 FORMAT (' Bisection.              XMIDPT',1P,D22.14)
99993 FORMAT (' Geo. bisection. XTRY,DAUX,DTRY',1P,3D22.14)
99992 FORMAT (' Polynomial fit accepted.  XTRY',1P,D22.14)
99991 FORMAT (' ----------------------------------------------------',/)
      END

      SUBROUTINE E04UDS(CENTRL,INFORM,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,
     *                  C1,C2,CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,W,
     *                  LENW,IUSER,USER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1095 (JUL 1993).
C
C     ******************************************************************
C     E04UDS evaluates any missing gradients.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 3-July-1986.
C     This version of E04UDS dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  THREE, FOUR
      PARAMETER         (THREE=3.0D+0,FOUR=4.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CDINT, FDINT, FDNORM, OBJF
      INTEGER           INFORM, LDCJ, LDCJU, LENW, N, NCNLN
      LOGICAL           CENTRL
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), GRAD(N), GRADU(N), HCNTRL(N),
     *                  HFORWD(N), USER(*), W(LENW), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           LFDSET, LVLDIF, NCDIFF, NFDIFF
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, DELTA, OBJF1, OBJF2, STEPBL,
     *                  STEPBU, XJ
      INTEGER           I, J, MODE, NCOLJ, NSTATE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      INFORM = 0
C
C     ==================================================================
C     Use the pre-assigned difference intervals to approximate the
C     derivatives.
C     ==================================================================
C     Use either the same interval for each element (LFDSET = 1),
C     or the intervals already in HFORWD or HCNTRL (LFDSET = 0 or 2).
C
      NSTATE = 0
      MODE = 0
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
      FDNORM = ZERO
C
      DO 80 J = 1, N
         XJ = X(J)
         NCOLJ = 0
         IF (NCDIFF.GT.0) THEN
            DO 20 I = 1, NCNLN
               IF (CJACU(I,J).EQ.RDUMMY) THEN
                  NEEDC(I) = 1
                  NCOLJ = NCOLJ + 1
               ELSE
                  NEEDC(I) = 0
               END IF
   20       CONTINUE
         END IF
C
         IF (NCOLJ.GT.0 .OR. GRADU(J).EQ.RDUMMY) THEN
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            IF (CENTRL) THEN
               IF (LFDSET.EQ.1) THEN
                  DELTA = CDINT
               ELSE
                  DELTA = HCNTRL(J)
               END IF
            ELSE
               IF (LFDSET.EQ.1) THEN
                  DELTA = FDINT
               ELSE
                  DELTA = HFORWD(J)
               END IF
            END IF
C
            DELTA = DELTA*(ONE+ABS(XJ))
            FDNORM = MAX(FDNORM,DELTA)
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) DELTA = -DELTA
C
            X(J) = XJ + DELTA
            IF (NCOLJ.GT.0) THEN
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 100
            END IF
C
            IF (GRADU(J).EQ.RDUMMY) THEN
               CALL OBJFUN(MODE,N,X,OBJF1,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
            END IF
C
            IF (CENTRL) THEN
C              ---------------------------------------------------------
C              Central differences.
C              ---------------------------------------------------------
               X(J) = XJ + DELTA + DELTA
C
               IF (NCOLJ.GT.0) THEN
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
C
                  DO 40 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (FOUR*C1(I)
     *                   -THREE*C(I)-C2(I))/(DELTA+DELTA)
   40             CONTINUE
               END IF
C
               IF (GRADU(J).EQ.RDUMMY) THEN
                  CALL OBJFUN(MODE,N,X,OBJF2,GRADU,NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 100
C
                  GRAD(J) = (FOUR*OBJF1-THREE*OBJF-OBJF2)/(DELTA+DELTA)
C
               END IF
            ELSE
C              ---------------------------------------------------------
C              Forward Differences.
C              ---------------------------------------------------------
               IF (NCOLJ.GT.0) THEN
                  DO 60 I = 1, NCNLN
                     IF (NEEDC(I).EQ.1) CJAC(I,J) = (C1(I)-C(I))/DELTA
   60             CONTINUE
               END IF
C
               IF (GRADU(J).EQ.RDUMMY) GRAD(J) = (OBJF1-OBJF)/DELTA
C
            END IF
         END IF
         X(J) = XJ
C
   80 CONTINUE
C
      RETURN
C
  100 INFORM = MODE
      RETURN
C
C     End of  E04UDS. (NPFD)
C
      END
      SUBROUTINE E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,NFREE,
     *                  NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,A,R,
     *                  T,RES,GQ,ZY,W,C,S,MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1072 (JUL 1993).
C
C     ******************************************************************
C     E04NCV  updates the factorization,  A(free) * (Z Y) = (0 T),  when
C     a constraint is added to the working set.  If  NRANK .gt. 0, the
C     factorization  ( R ) = PCQ  is also updated,  where  C  is the
C                    ( 0 )
C     least squares matrix,  R  is upper-triangular,  and  P  is an
C     orthogonal matrix.  The matrices  C  and  P  are not stored.
C
C     There are three separate cases to consider (although each case
C     shares code with another)...
C
C     (1) A free variable becomes fixed on one of its bounds when there
C         are already some general constraints in the working set.
C
C     (2) A free variable becomes fixed on one of its bounds when there
C         are only bound constraints in the working set.
C
C     (3) A general constraint (corresponding to row  IADD  of  A) is
C         added to the working set.
C
C     In cases (1) and (2), we assume that  KX(IFIX) = JADD.
C     In all cases,  JADD  is the index of the constraint being added.
C
C     If there are no general constraints in the working set,  the
C     matrix  Q = (Z Y)  is the identity and will not be touched.
C
C     If  NRES .GT. 0,  the row transformations are applied to the rows
C     of the  (N by NRES)  matrix  RES.
C     If  NGQ .GT. 0,  the column transformations are applied to the
C     columns of the  (NGQ by N)  matrix  GQ'.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October--1984.
C     Level-2 matrix routines  added 25-Apr-1988.
C     This version of  E04NCV  dated 28-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           IADD, IFIX, INFORM, JADD, LDA, LDR, LDT, LDZY,
     *                  MSGLVL, N, NACTIV, NFREE, NGQ, NRANK, NRES, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8,
     *                  EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, CONDBD, DTNEW, TDTMAX, TDTMIN
      INTEGER           I, NANEW, NFMIN, NPIV, NT
      LOGICAL           BOUND, OVERFL
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBW, F06FLF, F06FQF, F06QHF,
     *                  F06QKF, F06QNZ, F06QRF, F06QVF, F06QXF, F06QZZ,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Executable Statements ..
C
C     If the condition estimator of the updated factors is greater than
C     CONDBD,  a warning message is printed.
C
      CONDBD = ONE/EPSPT9
C
      OVERFL = .FALSE.
      BOUND = JADD .LE. N
C
      IF (BOUND) THEN
C        ===============================================================
C        A simple bound has entered the working set.  IADD  is not used.
C        ===============================================================
         IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) NACTIV, NZ, NFREE, IFIX, JADD, UNITQ
            CALL X04BAY(IPRINT,4,REC)
         END IF
         NANEW = NACTIV
C
         IF (UNITQ) THEN
C
C           Q  is not stored, but KX defines an ordering of the columns
C           of the identity matrix that implicitly define  Q.
C           Define the sequence of pairwise interchanges P that moves
C           the newly-fixed variable to position NFREE.
C           Reorder KX accordingly.
C
            DO 20 I = 1, NFREE - 1
               IF (I.GE.IFIX) THEN
                  W(I) = I + 1
                  KX(I) = KX(I+1)
               ELSE
                  W(I) = I
               END IF
   20       CONTINUE
C
         ELSE
C           ------------------------------------------------------------
C           Q  is stored explicitly.
C           ------------------------------------------------------------
C           Set  W = the  (IFIX)-th  row of  Q.
C           Move the  (NFREE)-th  row of  Q  to position  IFIX.
C
            CALL DCOPY(NFREE,ZY(IFIX,1),LDZY,W,1)
            IF (IFIX.LT.NFREE) THEN
               CALL DCOPY(NFREE,ZY(NFREE,1),LDZY,ZY(IFIX,1),LDZY)
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
C        ===============================================================
C        A general constraint has entered the working set.
C        IFIX  is not used.
C        ===============================================================
         IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99998) NACTIV, NZ, NFREE, IADD, JADD, UNITQ
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         NANEW = NACTIV + 1
C
C        Transform the incoming row of  A  by  Q'.  Use C as workspace.
C
         CALL DCOPY(N,A(IADD,1),LDA,W,1)
         CALL E04NBW(8,N,NZ,NFREE,LDZY,UNITQ,KX,W,ZY,C)
C
C        Check that the incoming row is not dependent upon those
C        already in the working set.
C
         DTNEW = DNRM2(NZ,W,1)
         IF (NACTIV.EQ.0) THEN
C
C           This is the only general constraint in the working set.
C
            COND = F06BLF(ASIZE,DTNEW,OVERFL)
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE
C
C           There are already some general constraints in the working
C           set. Update the estimate of the condition number.
C
            TDTMAX = MAX(DTNEW,DTMAX)
            TDTMIN = MIN(DTNEW,DTMIN)
            COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
         END IF
C
         IF (COND.GT.CONDMX .OR. OVERFL) GO TO 60
C
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,ZY,LDZY)
            UNITQ = .FALSE.
         END IF
      END IF
C
      IF (BOUND) THEN
         NPIV = NFREE
      ELSE
         NPIV = NZ
      END IF
C
      NT = MIN(NRANK,NPIV)
C
      IF (UNITQ) THEN
C        ---------------------------------------------------------------
C        Q (i.e., ZY) is not stored explicitly.
C        Apply the sequence of pairwise interchanges P that moves the
C        newly-fixed variable to position NFREE.
C        ---------------------------------------------------------------
         IF (NGQ.GT.0) CALL F06QKF('Left','Transpose',NFREE-1,W,NGQ,GQ,
     *                             N)
C
         IF (NRANK.GT.0) THEN
C
C           Apply the pairwise interchanges to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  s(1), s(2), ..., s(nt-1).
C
            CALL F06QNZ('Right',N,IFIX,NT,S,R,LDR)
C
            IF (NT.LT.NPIV) THEN
C
C              R is upper trapezoidal.  Apply the interchanges in
C              columns  nt  thru  npiv.
C
               DO 40 I = IFIX, NT - 1
                  W(I) = I
   40          CONTINUE
C
               CALL F06QKF('Right','Normal',NFREE-1,W,NT,R,LDR)
            END IF
C
C           Eliminate the subdiagonal elements of R with a left-hand
C           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
C           Apply P2 to RES.
C
            CALL F06QRF('Left ',N,IFIX,NT,C,S,R,LDR)
            IF (NRES.GT.0) CALL F06QXF('Left','Variable','Forwards',NT,
     *                                 NRES,IFIX,NT,C,S,RES,N)
         END IF
      ELSE
C        ---------------------------------------------------------------
C        Full matrix Q.  Define a sweep of plane rotations P such that
C                           Pw = beta*e(npiv).
C        The rotations are applied in the planes (1,2), (2,3), ...,
C        (npiv-1,npiv).  The rotations must be applied to ZY, R, T
C        and GQ'.
C        ---------------------------------------------------------------
         CALL F06FQF('Varble','Forwrds',NPIV-1,W(NPIV),W,1,C,S)
C
         IF (BOUND .AND. NACTIV.GT.0) THEN
C
            CALL DCOPY(NACTIV,S(NZ),1,W(NZ),1)
C
            S(NZ) = S(NZ)*T(NACTIV,NZ+1)
            T(NACTIV,NZ+1) = C(NZ)*T(NACTIV,NZ+1)
C
            CALL F06QZZ('Create',NACTIV,1,NACTIV,C(NZ+1),S(NZ+1),
     *                  T(1,NZ+1),LDT)
            CALL DCOPY(NACTIV,S(NZ),1,T(NACTIV,NZ),LDT-1)
C
            CALL DCOPY(NACTIV,W(NZ),1,S(NZ),1)
         END IF
C
         IF (NGQ.GT.0) CALL F06QXF('Left ','Variable','Forwards',NPIV,
     *                             NGQ,1,NPIV,C,S,GQ,N)
         CALL F06QXF('Right','Variable','Forwards',NFREE,NFREE,1,NPIV,C,
     *               S,ZY,LDZY)
C
         IF (NRANK.GT.0) THEN
C
C           Apply the rotations to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  s(1),  s(2), ..., s(nt-1).
C
            NT = MIN(NRANK,NPIV)
            CALL F06QVF('Right',N,1,NT,C,S,R,LDR)
C
            IF (NT.LT.NPIV) THEN
C
C              R is upper trapezoidal.  Pretend R is (nt x n) and
C              apply the rotations in columns  nt  thru  npiv.
C
               CALL F06QXF('Right','Variable','Forwards',NT,N,NT,NPIV,C,
     *                     S,R,LDR)
            END IF
C
C           Eliminate the subdiagonal elements of R with a left-hand
C           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
C           Apply P2 to RES.
C
            CALL F06QRF('Left ',N,1,NT,C,S,R,LDR)
            IF (NRES.GT.0) CALL F06QXF('Left','Variable','Forwards',NT,
     *                                 NRES,1,NT,C,S,RES,N)
         END IF
C
         IF (BOUND) THEN
C
C           The last row and column of ZY has been transformed to plus
C           or minus the unit vector E(NFREE).  We can reconstitute the
C           columns of GQ and R corresponding to the new fixed variable.
C
            IF (W(NFREE).LT.ZERO) THEN
               NFMIN = MIN(NRANK,NFREE)
               IF (NFMIN.GT.0) CALL DSCAL(NFMIN,-ONE,R(1,NFREE),1)
               IF (NGQ.GT.0) CALL DSCAL(NGQ,-ONE,GQ(NFREE,1),N)
            END IF
C
C           ------------------------------------------------------------
C           The diagonals of T have been altered.  Recompute the
C           largest and smallest values.
C           ------------------------------------------------------------
            IF (NACTIV.GT.0) THEN
               CALL F06FLF(NACTIV,T(NACTIV,NZ),LDT-1,TDTMAX,TDTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
         ELSE
C           ------------------------------------------------------------
C           General constraint.  Install the new row of T.
C           ------------------------------------------------------------
            CALL DCOPY(NANEW,W(NZ),1,T(NANEW,NZ),LDT)
         END IF
      END IF
C
C     ==================================================================
C     Prepare to exit.  Check the magnitude of the condition estimator.
C     ==================================================================
   60 IF (NANEW.GT.0) THEN
         IF (COND.LT.CONDMX .AND. .NOT. OVERFL) THEN
C
C           The factorization has been successfully updated.
C
            INFORM = 0
            DTMAX = TDTMAX
            DTMIN = TDTMIN
            IF (COND.GE.CONDBD) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99997) JADD
                  CALL X04BAY(IPRINT,5,REC)
               END IF
            END IF
         ELSE
C
C           The proposed working set appears to be linearly dependent.
C
            INFORM = 1
            IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
               WRITE (REC,FMT=99996)
               CALL X04BAY(IPRINT,2,REC)
               IF (BOUND) THEN
                  WRITE (REC,FMT=99995) ASIZE, DTMAX, DTMIN
                  CALL X04BAY(IPRINT,3,REC)
               ELSE
                  IF (NACTIV.GT.0) THEN
                     WRITE (REC,FMT=99994) ASIZE, DTMAX, DTMIN, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  ELSE
                     WRITE (REC,FMT=99993) ASIZE, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCV. (LSADD)
C
99999 FORMAT (/' //E04NCV //  Simple bound added.',/' //E04NCV //  NAC',
     *       'TIV    NZ NFREE  IFIX  JADD UNITQ',/' //E04NCV //  ',5I6,
     *       L6)
99998 FORMAT (/' //E04NCV //  General constraint added.           ',
     *       /' //E04NCV //  NACTIV    NZ NFREE  IADD  JADD UNITQ',
     *       /' //E04NCV //  ',5I6,L6)
99997 FORMAT (/' XXX  Serious ill-conditioning in the working set afte',
     *       'r adding constraint ',I5,/' XXX  Overflow may occur in s',
     *       'ubsequent iterations.',//)
99996 FORMAT (/' //E04NCV //  Dependent constraint rejected.')
99995 FORMAT (/' //E04NCV //     ASIZE     DTMAX     DTMIN        ',
     *       /' //E04NCV //',1P,3D10.2)
99994 FORMAT (/' //E04NCV //     ASIZE     DTMAX     DTMIN     DTNEW',
     *       /' //E04NCV //',1P,4D10.2)
99993 FORMAT (/' //E04NCV //     ASIZE     DTNEW',/' //E04NCV //',1P,
     *       2D10.2)
      END
      SUBROUTINE E04XAW(INFORM,LVLDER,MSGLVL,NCSET,N,NCNLN,LDCJ,LDCJU,
     *                  BIGBND,EPSRF,OKTOL,FDCHK,XNORM,CONFUN,NEEDC,BL,
     *                  BU,C,C1,CJAC,CJACU,CJDX,DX,ERR,X,Y,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16A REVISED. IER-999 (JUN 1993).
C
C     ******************************************************************
C     E04XAW  checks if the gradients of the constraints have been coded
C     correctly.
C
C     On input,  the values of the constraints at the point X are stored
C     in C.  Their corresponding gradients are stored in CJACU.  If any
C     Jacobian element has not been specified,  it will have a dummy
C     value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods.  If this proves
C     satisfactory and no further information is desired, E04XAW is
C     terminated. Otherwise, E04XAZ is called to give optimal
C     step-sizes and a central-difference approximation to each
C     element of the Jacobian for which a test is deemed necessary,
C     either by the program or the user.
C
C     LVRFYC has the following meaning...
C
C     -1        do not perform any check.
C     0        do the cheap test only.
C     2 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     Level-2 matrix routines added 18-May-1988.
C     This version of E04XAW dated 14-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  ZERO, HALF, POINT9
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,POINT9=0.9D+0)
      DOUBLE PRECISION  ONE, TWO, TEN
      PARAMETER         (ONE=1.0D+0,TWO=2.0D+0,TEN=1.0D+1)
      CHARACTER*4       LBAD, LGOOD
      PARAMETER         (LBAD='BAD?',LGOOD='  OK')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, EPSRF, FDCHK, OKTOL, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LVLDER, MSGLVL, N, NCNLN,
     *                  NCSET
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), CJDX(*), DX(N), ERR(*), USER(*),
     *                  X(N), Y(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, LVRFYC, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG), JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, CIJ, CJDIFF, CJSIZE,
     *                  COLMAX, DXJ, DXMULT, EMAX, EPSACI, ERRBND, F1,
     *                  F2, FDEST, H, HOPT, HPHI, SDEST, SIGNH, STEPBL,
     *                  STEPBU, XJ
      INTEGER           I, IMAX, INFO, IROW, ITER, ITMAX, J, J3, J4,
     *                  JCOL, MODE, NCHECK, NCOLJ, NGOOD, NSTATE, NWRONG
      LOGICAL           CONST, DEBUG, DONE, FIRST, HEADNG, NEEDED, OK
      CHARACTER*4       KEY
C     .. Local Arrays ..
      CHARACTER*18      RESULT(0:4)
      CHARACTER*120     REC(4)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04XAZ, F06DBF, F06FBF,
     *                  F06QFF, F06QHF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Data statements ..
      DATA              RESULT/'                 ', 'Constant?      ',
     *                  'Linear or odd?   ', 'Too nonlinear?',
     *                  'Small derivative?'/
C     .. Executable Statements ..
C
      INFORM = 0
      NEEDED = NCNLN .GT. 0 .AND. LVRFYC .EQ. 0 .OR. LVRFYC .EQ. 2 .OR.
     *         LVRFYC .EQ. 3
      IF ( .NOT. NEEDED) RETURN
C
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,4,REC)
      END IF
      DEBUG = NPDBG .AND. INPDBG(5) .GT. 0
      NSTATE = 0
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Perform the cheap test.
C     ==================================================================
      H = (ONE+XNORM)*FDCHK
C
      IF (N.LE.100) THEN
         DXMULT = 0.9
      ELSE IF (N.LE.250) THEN
         DXMULT = 0.99
      ELSE
         DXMULT = 0.999
      END IF
C
      DXJ = ONE/N
      DO 20 J = 1, N
         DX(J) = DXJ
         DXJ = -DXJ*DXMULT
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Do not perturb  X(J)  if the  J-th  column contains any
C     unknown elements.  Compute the directional derivative for each
C     constraint gradient.
C     ------------------------------------------------------------------
      NCHECK = 0
      DO 60 J = 1, N
         DO 40 I = 1, NCNLN
            IF (CJAC(I,J).EQ.RDUMMY) THEN
               DX(J) = ZERO
               GO TO 60
            END IF
   40    CONTINUE
         NCHECK = NCHECK + 1
C
         XJ = X(J)
         STEPBL = -ONE
         STEPBU = ONE
         IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
         IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J)) STEPBU = MIN(STEPBU,
     *       BU(J)-XJ)
C
         IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
            DX(J) = DX(J)*STEPBL
         ELSE
            DX(J) = DX(J)*STEPBU
         END IF
   60 CONTINUE
C
      IF (NCHECK.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99995)
            CALL X04BAY(IPRINT,2,REC)
         END IF
      ELSE
C
C        Compute  (Jacobian)*DX.
C
         CALL DGEMV('Normal',NCNLN,N,ONE,CJACU,LDCJU,DX,1,ZERO,CJDX,1)
C
C        ---------------------------------------------------------------
C        Make forward-difference approximation along DX.
C        ---------------------------------------------------------------
         CALL DCOPY(N,X,1,Y,1)
         CALL DAXPY(N,H,DX,1,Y,1)
C
         CALL F06DBF(NCNLN,(1),NEEDC,1)
C
         MODE = 0
         CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C1,CJACU,NSTATE,IUSER,
     *               USER)
         IF (MODE.LT.0) GO TO 160
C
C        Set  ERR = (C1 - C)/H  - Jacobian*DX.  This should be small.
C
         DO 80 I = 1, NCNLN
            ERR(I) = (C1(I)-C(I))/H - CJDX(I)
   80    CONTINUE
         IMAX = IDAMAX(NCNLN,ERR,1)
         EMAX = ABS(ERR(IMAX))/(ABS(CJDX(IMAX))+ONE)
C
         IF (MSGLVL.GT.0) THEN
            IF (EMAX.LE.OKTOL) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,2,REC)
            END IF
            WRITE (REC,FMT=99996) EMAX, IMAX
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (EMAX.GE.POINT9) INFORM = 1
      END IF
C
C     ==================================================================
C     Element-wise check.
C     ==================================================================
      IF (LVRFYC.GE.2) THEN
         IF (LVLDER.EQ.3) THEN
C
C           Recompute the Jacobian to find the non-constant elements.
C
            CALL F06QHF('General',NCNLN,N,RDUMMY,RDUMMY,CJACU,LDCJU)
C
            CALL F06DBF(NCNLN,(1),NEEDC,1)
            NSTATE = 0
            MODE = 2
C
            CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,NSTATE,
     *                  IUSER,USER)
            IF (MODE.LT.0) GO TO 160
C
         END IF
C
         CALL F06DBF(NCNLN,(0),NEEDC,1)
C
         ITMAX = 3
         NCHECK = 0
         NWRONG = 0
         NGOOD = 0
         COLMAX = -ONE
         JCOL = 0
         IROW = 0
         MODE = 0
         J3 = JVERFY(3)
         J4 = JVERFY(4)
C
C        ---------------------------------------------------------------
C        Loop over each column.
C        ---------------------------------------------------------------
         DO 140 J = J3, J4
C
            CALL F06FBF(NCNLN,ZERO,ERR,1)
            NCOLJ = 0
            HEADNG = .TRUE.
            XJ = X(J)
C
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            SIGNH = ONE
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) SIGNH = -ONE
C
            DO 120 I = 1, NCNLN
               EPSACI = EPSRF*(ONE+ABS(C(I)))
C
               IF (CJACU(I,J).NE.RDUMMY) THEN
C                 ------------------------------------------------------
C                 Check this Jacobian element.
C                 ------------------------------------------------------
                  NCHECK = NCHECK + 1
                  NCOLJ = NCOLJ + 1
                  NEEDC(I) = 1
C
                  CIJ = CJAC(I,J)
                  CJSIZE = ABS(CIJ)
C                 ------------------------------------------------------
C                 Find a finite-difference interval by iteration.
C                 ------------------------------------------------------
                  ITER = 0
                  HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
                  H = TEN*HOPT*SIGNH
                  CDEST = ZERO
                  SDEST = ZERO
                  FIRST = .TRUE.
C
C                 +                REPEAT
  100             X(J) = XJ + H
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 160
                  F1 = C1(I)
C
                  X(J) = XJ + H + H
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                        NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 160
                  F2 = C1(I)
C
                  CALL E04XAZ(DEBUG,DONE,FIRST,EPSACI,EPSRF,C(I),INFO,
     *                        ITER,ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,
     *                        H,HOPT,HPHI)
C
C                 +                UNTIL     DONE
                  IF ( .NOT. DONE) GO TO 100
C
C                 ------------------------------------------------------
C                 Exit for this element.
C                 ------------------------------------------------------
                  CJDIFF = CDEST
                  ERR(I) = ABS(CJDIFF-CIJ)/(CJSIZE+ONE)
C
                  OK = ERR(I) .LE. OKTOL
                  IF (OK) THEN
                     KEY = LGOOD
                     NGOOD = NGOOD + 1
                  ELSE
                     KEY = LBAD
                     NWRONG = NWRONG + 1
                  END IF
C
                  IF (MSGLVL.GT.0) THEN
                     CONST = OK .AND. INFO .EQ. 1 .AND. ABS(CIJ)
     *                       .LT. EPSPT8
                     IF ( .NOT. CONST) THEN
                        IF (HEADNG) THEN
                           WRITE (REC,FMT=99994)
                           CALL X04BAY(IPRINT,4,REC)
                           IF (OK) THEN
                              WRITE (REC,FMT=99993) J, XJ, HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER
                           ELSE
                              WRITE (REC,FMT=99992) J, XJ, HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER, RESULT(INFO)
                           END IF
                           CALL X04BAF(IPRINT,REC(1))
                           HEADNG = .FALSE.
                        ELSE
                           IF (OK) THEN
                              WRITE (REC,FMT=99991) HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER
                           ELSE
                              WRITE (REC,FMT=99990) HOPT, I, CIJ,
     *                          CJDIFF, KEY, ITER, RESULT(INFO)
                           END IF
                           CALL X04BAF(IPRINT,REC(1))
                        END IF
                     END IF
                  END IF
                  NEEDC(I) = 0
               END IF
  120       CONTINUE
C
C           ------------------------------------------------------------
C           Finished with this column.
C           ------------------------------------------------------------
            IF (NCOLJ.GT.0) THEN
               IMAX = IDAMAX(NCNLN,ERR,1)
               EMAX = ABS(ERR(IMAX))
C
               IF (EMAX.GE.COLMAX) THEN
                  IROW = IMAX
                  JCOL = J
                  COLMAX = EMAX
               END IF
            END IF
            X(J) = XJ
C
  140    CONTINUE
C
         INFORM = 0
         IF (COLMAX.GE.POINT9) INFORM = 1
C
         IF (MSGLVL.GT.0) THEN
            IF (NCHECK.EQ.0) THEN
               WRITE (REC,FMT=99986) NCSET
               CALL X04BAF(IPRINT,REC(1))
            ELSE
               IF (NWRONG.EQ.0) THEN
                  WRITE (REC,FMT=99989) NGOOD, NCHECK, J3, J4
               ELSE
                  WRITE (REC,FMT=99988) NWRONG, NCHECK, J3, J4
               END IF
               CALL X04BAY(IPRINT,3,REC)
               WRITE (REC,FMT=99987) COLMAX, IROW, JCOL
               CALL X04BAY(IPRINT,3,REC)
            END IF
         END IF
C
      END IF
C
C     Copy  ( constants + gradients + dummy values )  back into CJACU.
C
      CALL F06QFF('General',NCNLN,N,CJAC,LDCJ,CJACU,LDCJU)
C
      RETURN
C
  160 INFORM = MODE
      RETURN
C
C
C     End of  E04XAW. (CHCJAC)
C
99999 FORMAT (//' Verification of the constraint gradients.',/' ------',
     *       '-----------------------------------')
99998 FORMAT (/' The constraint Jacobian seems to be ok.')
99997 FORMAT (/' XXX  The constraint Jacobian seems to be incorrect.')
99996 FORMAT (/' The largest relative error was',1P,D12.2,'  in constr',
     *       'aint',I5,/)
99995 FORMAT (/' Every column contains a constant or missing element.')
99994 FORMAT (//' Column    X(J)     DX(J)    Row    Jacobian Value   ',
     *       '   Difference Approxn  Itns',/)
99993 FORMAT (I7,1P,2D10.2,I5,1P,2D18.8,2X,A4,I6)
99992 FORMAT (I7,1P,2D10.2,I5,1P,2D18.8,2X,A4,I6,2X,A18)
99991 FORMAT (17X,1P,D10.2,I5,1P,2D18.8,2X,A4,I6)
99990 FORMAT (17X,1P,D10.2,I5,1P,2D18.8,2X,A4,I6,2X,A18)
99989 FORMAT (/I7,'  constraint Jacobian elements out of the',I6,/9X,
     *       'set in cols',I6,'  through',I6,'  seem to be ok.')
99988 FORMAT (/' XXX  There seem to be',I6,'  incorrect Jacobian eleme',
     *       'nts out of the',I6,/8X,'set in cols',I6,'  through',I6)
99987 FORMAT (/' The largest relative error was',1P,D12.2,'  in row',I5,
     *       ',  column',I5,/)
99986 FORMAT (' All',I6,'   assigned Jacobian elements are constant.')
      END
      SUBROUTINE E04NCQ(LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,LDA,
     *                  LDZY,LDR,NRANK,NUMINF,NRZ,KX,CTP,PNORM,A,AP,RES,
     *                  HZ,P,GQ,CQ,R,ZY,WORK)
C     MARK 16 REVISED. IER-1068 (JUL 1993).
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C     E04NCQ  computes the following quantities for  E04NCZ.
C     (1) The vector  (hz1) = (Rz1)(pz1).
C         If X is not yet feasible,  the product is computed directly.
C         If  Rz1 is singular,  hz1  is zero.  Otherwise  hz1  satisfies
C         the equations
C                        Rz1'hz1 = -gz1,
C         where  g  is the total gradient.  If there is no linear term
C         in the objective,  hz1  is set to  dz1  directly.
C     (2) The search direction P (and its 2-norm).  The vector P is
C         defined as  Z*(pz1), where  (pz1)  depends upon whether or
C         not X is feasible and the nonsingularity of  (Rz1).
C         If  NUMINF .GT. 0,  (pz1)  is the steepest-descent direction.
C         Otherwise,  x  is the solution of the  NRZ*NRZ  triangular
C         system   (Rz1)*(pz1) = (hz1).
C     (3) The vector Ap,  where A is the matrix of linear constraints.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level 2 BLAS added 11-June-1986.
C     This version of E04NCQ dated 28-July-1987.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTP, PNORM
      INTEGER           LDA, LDR, LDZY, N, NCLIN, NFREE, NRANK, NRZ,
     *                  NUMINF
      LOGICAL           LINOBJ, SINGLR, UNITGZ, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AP(*), CQ(*), GQ(N), HZ(*), P(N),
     *                  R(LDR,*), RES(*), WORK(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  GTP
      INTEGER           I, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRSV, E04NBW, F06FBF,
     *                  X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Executable Statements ..
      IF (SINGLR) THEN
C        ---------------------------------------------------------------
C        The triangular factor for the current objective function is
C        singular,  i.e., the objective is linear along the last column
C        of Z1.  This can only occur when UNITGZ is TRUE.
C        ---------------------------------------------------------------
         IF (NRZ.GT.1) THEN
            CALL DCOPY(NRZ-1,R(1,NRZ),1,P,1)
            CALL DTRSV('U','N','N',NRZ-1,R,LDR,P,1)
         END IF
         P(NRZ) = -ONE
C
         GTP = DDOT(NRZ,GQ,1,P,1)
         IF (GTP.GT.ZERO) CALL DSCAL(NRZ,(-ONE),P,1)
C
         IF (NRZ.LE.NRANK) THEN
            IF (NUMINF.EQ.0) THEN
               IF (UNITGZ) THEN
                  HZ(NRZ) = R(NRZ,NRZ)*P(NRZ)
               ELSE
                  CALL F06FBF(NRZ,(ZERO),HZ,1)
               END IF
            ELSE
               HZ(1) = R(1,1)*P(1)
            END IF
         END IF
      ELSE
C        ---------------------------------------------------------------
C        The objective is quadratic in the space spanned by Z1.
C        ---------------------------------------------------------------
         IF (LINOBJ) THEN
            IF (UNITGZ) THEN
               IF (NRZ.GT.1) CALL F06FBF(NRZ-1,(ZERO),HZ,1)
               HZ(NRZ) = -GQ(NRZ)/R(NRZ,NRZ)
            ELSE
               CALL DCOPY(NRZ,GQ,1,HZ,1)
               CALL DSCAL(NRZ,(-ONE),HZ,1)
               CALL DTRSV('U','T','N',NRZ,R,LDR,HZ,1)
            END IF
         ELSE
            CALL DCOPY(NRZ,RES,1,HZ,1)
         END IF
C
C        Solve  Rz1*pz1 = hz1.
C
         CALL DCOPY(NRZ,HZ,1,P,1)
         CALL DTRSV('U','N','N',NRZ,R,LDR,P,1)
C
      END IF
C
C     Compute  p = Z1*pz1  and its norm.
C
      IF (LINOBJ) CTP = DDOT(NRZ,CQ,1,P,1)
      PNORM = DNRM2(NRZ,P,1)
C
      CALL E04NBW(1,N,NRZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
C
      IF (LSDBG .AND. ILSDBG(2).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,2,REC)
         DO 20 I = 1, N, 5
            WRITE (REC,FMT=99997) (P(J),J=I,MIN(I+4,N))
            CALL X04BAF(IPRINT,REC(1))
   20    CONTINUE
      END IF
C
C     Compute  Ap.
C
      IF (NCLIN.GT.0) THEN
         CALL DGEMV('No transpose',NCLIN,N,ONE,A,LDA,P,1,ZERO,AP,1)
         IF (LSDBG .AND. ILSDBG(2).GT.0) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
            DO 40 I = 1, N, 5
               WRITE (REC,FMT=99997) (AP(J),J=I,MIN(I+4,N))
               CALL X04BAF(IPRINT,REC(1))
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCQ. (LSGETP)
C
99999 FORMAT (/' //E04NCQ//   P ... ')
99998 FORMAT (/' //E04NCQ//  AP ... ')
99997 FORMAT (1P,5D15.5)
      END
      SUBROUTINE E04NCP(PRBTYP,LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDZY,LDR,NRANK,NZ,NRZ,ISTATE,KX,BIGBND,
     *                  TOLRNK,NUMINF,SUMINF,BL,BU,A,RES,FEATOL,GQ,CQ,R,
     *                  X,WTINF,ZY,WRK)
C     MARK 16 REVISED. IER-1067 (JUL 1993).
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     ******************************************************************
C     E04NCP  finds the number and weighted sum of infeasibilities for
C     the bounds and linear constraints.   An appropriate transformed
C     gradient vector is returned in  GQ.
C
C     Positive values of  ISTATE(j)  will not be altered.  These mean
C     the following...
C
C               1             2           3
C           a'x = bl      a'x = bu     bl = bu
C
C     Other values of  ISTATE(j)  will be reset as follows...
C           a'x lt bl     a'x gt bu     a'x free
C              - 2           - 1           0
C
C     If  x  is feasible,  E04NCP computes the vector Q(free)'g(free),
C     where  g  is the gradient of the the sum of squares plus the
C     linear term.  The matrix Q is of the form
C                    ( Q(free)  0       ),
C                    (   0      I(fixed))
C     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
C     the matrix of constraints in the working set.  The transformed
C     gradients are stored in GQ.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCP dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, SUMINF, TOLRNK
      INTEGER           LDA, LDR, LDZY, N, NCLIN, NFREE, NRANK, NRZ,
     *                  NUMINF, NZ
      LOGICAL           LINOBJ, SINGLR, UNITGZ, UNITQ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(*), BU(*), CQ(*), FEATOL(*), GQ(N),
     *                  R(LDR,*), RES(*), WRK(N), WTINF(*), X(N),
     *                  ZY(LDZY,*)
      INTEGER           ISTATE(*), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CTX, FEASJ, ROWNRM, S, WEIGHT
      INTEGER           J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           F06KLF
      EXTERNAL          DDOT, DNRM2, F06KLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRMV, E04NBW, F06FBF,
     *                  F06FDF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      BIGUPP = BIGBND
      BIGLOW = -BIGBND
C
      NUMINF = 0
      SUMINF = ZERO
      CALL F06FBF(N,ZERO,GQ,1)
C
      DO 40 J = 1, N + NCLIN
         IF (ISTATE(J).LE.0) THEN
            FEASJ = FEATOL(J)
            IF (J.LE.N) THEN
               CTX = X(J)
            ELSE
               K = J - N
               CTX = DDOT(N,A(K,1),LDA,X,1)
            END IF
            ISTATE(J) = 0
C
C           See if the lower bound is violated.
C
            IF (BL(J).GT.BIGLOW) THEN
               S = BL(J) - CTX
               IF (S.GT.FEASJ) THEN
                  ISTATE(J) = -2
                  WEIGHT = -WTINF(J)
                  GO TO 20
               END IF
            END IF
C
C           See if the upper bound is violated.
C
            IF (BU(J).GE.BIGUPP) GO TO 40
            S = CTX - BU(J)
            IF (S.LE.FEASJ) GO TO 40
            ISTATE(J) = -1
            WEIGHT = WTINF(J)
C
C           Add the infeasibility.
C
   20       NUMINF = NUMINF + 1
            SUMINF = SUMINF + ABS(WEIGHT)*S
            IF (J.LE.N) THEN
               GQ(J) = WEIGHT
            ELSE
               CALL DAXPY(N,WEIGHT,A(K,1),LDA,GQ,1)
            END IF
         END IF
   40 CONTINUE
C
C     ------------------------------------------------------------------
C     Install  GQ,  the transformed gradient.
C     ------------------------------------------------------------------
      SINGLR = .FALSE.
      UNITGZ = .TRUE.
C
      IF (NUMINF.GT.0) THEN
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,GQ,ZY,WRK)
         UNITGZ = .TRUE.
      ELSE IF (NUMINF.EQ.0 .AND. PRBTYP.EQ.'FP') THEN
         CALL F06FBF(N,ZERO,GQ,1)
      ELSE
C
C        Ready for the optimality phase.
C        Set NRZ so that Rz1 is nonsingular.
C
         IF (NRANK.EQ.0) THEN
            IF (LINOBJ) THEN
               CALL DCOPY(N,CQ,1,GQ,1)
            ELSE
               CALL F06FBF(N,ZERO,GQ,1)
            END IF
            NRZ = 0
         ELSE
C
C           Compute GQ = - R' * (transformed residual)
C
            CALL F06FDF(NRANK,-ONE,RES,1,GQ,1)
            CALL DTRMV('U','T','N',NRANK,R,LDR,GQ,1)
            IF (NRANK.LT.N) CALL DGEMV('T',NRANK,N-NRANK,-ONE,
     *                                 R(1,NRANK+1),LDR,RES,1,ZERO,
     *                                 GQ(NRANK+1),1)
            IF (LINOBJ) CALL DAXPY(N,ONE,CQ,1,GQ,1)
            UNITGZ = .FALSE.
            ROWNRM = DNRM2(N,R(1,1),LDR)
            IF (ROWNRM.LE.TOLRNK .OR. ABS(R(1,1)).LE.ROWNRM*TOLRNK) THEN
               NRZ = 0
            ELSE
               NRZ = F06KLF(MIN(NRANK,NZ),R,LDR+1,TOLRNK)
            END IF
         END IF
         SINGLR = .FALSE.
      END IF
C
      RETURN
C
C     End of  E04NCP. (LSGSET)
C
      END

      SUBROUTINE E04MFN(SUBR,MSG,V,LENV)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFN  prints the array V in debug format.
C
C     Original version dated 17-Jul-1987.
C     This version of  E04MFN  dated  31-Jan-1988.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LENV
      CHARACTER*6       SUBR
      CHARACTER*(*)     MSG
C     .. Array Arguments ..
      DOUBLE PRECISION  V(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      INTEGER           I, II
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      IF (LENV.LE.0) THEN
         WRITE (REC,FMT=99999) SUBR, MSG
         CALL X04BAY(IPRINT,2,REC)
      ELSE
         WRITE (REC,FMT=99998) SUBR, MSG
         CALL X04BAY(IPRINT,2,REC)
         DO 20 I = 1, LENV, 5
            WRITE (REC,FMT=99997) (V(II),II=I,MIN(I+4,LENV))
            CALL X04BAF(IPRINT,REC(1))
   20    CONTINUE
      END IF
C
      RETURN
C
C     End of  E04MFN. (CMMSG1)
C
99999 FORMAT (/' //',A6,'//  ',A)
99998 FORMAT (/' //',A6,'//  ',A,' ... ')
99997 FORMAT (1P,5D15.5)
      END
      SUBROUTINE E04NCR(N,NCLIN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,NVIOL,
     *                  AX,BL,BU,FEATOL,X,WORK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1069 (JUL 1993).
C
C     ******************************************************************
C     E04NCR  computes the following...
C     (1)  The number of constraints that are violated by more
C          than  FEATOL  and the 2-norm of the constraint violations.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      April    1984.
C     This version of  E04NCR  dated  17-October-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CVNORM, ERRMAX
      INTEGER           JMAX, N, NCLIN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), WORK(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES, TOLJ
      INTEGER           I, IS, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Executable Statements ..
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Compute NVIOL,  the number of constraints violated by more than
C     FEATOL,  and CVNORM,  the 2-norm of the constraint violations and
C     residuals of the constraints in the working set.
C     ==================================================================
      NVIOL = 0
C
      DO 40 J = 1, N + NCLIN
         FEASJ = FEATOL(J)
         IS = ISTATE(J)
         RES = ZERO
C
         IF (IS.GE.0 .AND. IS.LT.4) THEN
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               I = J - N
               CON = AX(I)
            END IF
C
            TOLJ = FEASJ
C
C           Check for constraint violations.
C
            IF (BL(J).GT.BIGLOW) THEN
               RES = BL(J) - CON
               IF (RES.GT.FEASJ) NVIOL = NVIOL + 1
               IF (RES.GT.TOLJ) GO TO 20
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               RES = BU(J) - CON
               IF (RES.LT.(-FEASJ)) NVIOL = NVIOL + 1
               IF (RES.LT.(-TOLJ)) GO TO 20
            END IF
C
C           This constraint is satisfied,  but count the residual as a
C           violation if the constraint is in the working set.
C
            IF (IS.LE.0) RES = ZERO
            IF (IS.EQ.1) RES = BL(J) - CON
            IF (IS.GE.2) RES = BU(J) - CON
            IF (ABS(RES).GT.FEASJ) NVIOL = NVIOL + 1
         END IF
   20    WORK(J) = RES
   40 CONTINUE
C
      JMAX = IDAMAX(N+NCLIN,WORK,1)
      ERRMAX = ABS(WORK(JMAX))
C
      IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) ERRMAX, JMAX
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      CVNORM = DNRM2(N+NCLIN,WORK,1)
C
      RETURN
C
C
C     End of  E04NCR. (LSFEAS)
C
99999 FORMAT (/' //E04NCR//  The maximum violation is ',1P,D14.2,' in ',
     *       'constraint',I5)
      END
      SUBROUTINE E04MFL(MSGLVL,N,NRZ,NZ,ZEROLM,NOTOPT,NUMINF,TRUSML,
     *                  SMLLST,JSMLST,TINYST,JTINY,GQ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFL  updates JSMLST and SMLLST when there are artificial
C     constraints.
C
C     On input,  JSMLST  is the index of the minimum of the set of
C     adjusted multipliers.
C     On output, a negative JSMLST defines the index in Q'g of the
C     artificial constraint to be deleted.
C
C     Original version written 17-Jan-1988.
C     This version of E04MFL dated  23-Jul-1991.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SMLLST, TINYST, TRUSML, ZEROLM
      INTEGER           JSMLST, JTINY, MSGLVL, N, NOTOPT, NRZ, NUMINF,
     *                  NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  GQ(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  RLAM
      INTEGER           J, K, KK, LENGTH
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      DO 20 J = NRZ + 1, NZ
         RLAM = -ABS(GQ(J))
C
         IF (RLAM.LT.ZEROLM) THEN
            IF (NUMINF.EQ.0) NOTOPT = NOTOPT + 1
C
            IF (RLAM.LT.SMLLST) THEN
               TRUSML = GQ(J)
               SMLLST = RLAM
               JSMLST = -J
            END IF
C
         ELSE IF (RLAM.LT.TINYST) THEN
            TINYST = RLAM
            JTINY = -J
         END IF
   20 CONTINUE
C
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAY(ISUMM,2,REC)
            LENGTH = NZ - NRZ
            DO 40 K = 1, LENGTH, 4
               WRITE (REC,FMT=99998) (GQ(KK),KK=K,MIN(K+3,LENGTH))
               CALL X04BAF(ISUMM,REC(1))
   40       CONTINUE
         END IF
      END IF
C
      IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
         WRITE (REC,FMT=99997) JSMLST, SMLLST
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      RETURN
C
C     End of  E04MFL.  (CMMUL2)
C
99999 FORMAT (/' Multipliers for the artificial constraints        ')
99998 FORMAT (4(5X,1P,D11.2))
99997 FORMAT (/' //E04MFL//  JSMLST     SMLLST',/' //E04MFL//  ',I6,1P,
     *       D11.2)
      END
      SUBROUTINE E04UCT(KTCOND,CONVRG,LSUMRY,MSGNP,MSGQP,LDR,LDT,N,
     *                  NCLIN,NCNLN,NCTOTL,NACTIV,LINACT,NLNACT,NZ,
     *                  NFREE,MAJIT0,MAJITS,MINITS,ISTATE,ALFA,NFUN,
     *                  CONDHZ,CONDH,CONDT,OBJALF,OBJF,GFNORM,GZNORM,
     *                  CVNORM,AX,C,R,T,VIOLN,X,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-723 (DEC 1989).
C     MARK 16 REVISED. IER-1088 (JUL 1993).
C
C     ==================================================================
C     E04UCT  prints various levels of output for E04UCZ and E04UPZ.
C
C           Msg        Cumulative result
C           ---        -----------------
C
C        le   0        no output.
C
C        eq   1        nothing now (but full output later).
C
C        eq   5        one terse line of output.
C
C        ge  10        same as 5 (but full output later).
C
C        ge  20        objective function,  x,  Ax  and  c.
C
C        ge  30        diagonals of  T  and  R.
C
C     Debug print is performed depending on the logical variable NPDBG.
C     NPDBG is set true when IDBG major iterations have been performed.
C     At this point,  printing is done according to a string of binary
C     digits of the form CLSVT (stored in the integer array INPDBG).
C
C     C  set 'on' gives detailed information from the checking routines.
C     L  set 'on' gives information from the linesearch.
C     S  set 'on' gives information from the maximum step routine E04UDT
C     V  set 'on' gives various vectors in  E04UCZ  and its auxiliaries.
C     T  set 'on' gives a trace of which routine was called and an
C                 indication of the progress of the run.
C     For example, `Major debug level 11000' gives much output from
C                  the checking routines and the linesearch.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written November-1982.
C     This version of  E04UCT  dated  21-Oct-92.
C     ==================================================================
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDH, CONDHZ, CONDT, CVNORM, GFNORM,
     *                  GZNORM, OBJALF, OBJF
      INTEGER           LDR, LDT, LINACT, MAJIT0, MAJITS, MINITS, MSGNP,
     *                  MSGQP, N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE,
     *                  NFUN, NLNACT, NZ
      LOGICAL           CONVRG
      CHARACTER*5       LSUMRY
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), C(*), R(LDR,*), T(LDT,*), VIOLN(*),
     *                  WORK(N), X(N)
      INTEGER           ISTATE(NCTOTL)
      LOGICAL           KTCOND(2)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RHODMP, RHOMAX, RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           INCRUN, NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  CVIOLS
      INTEGER           I, INCT, J, K, MJR, MNR, NDF, NEVAL
      LOGICAL           FIRST, NEWSET, NLNCON, PRTHDR
C     .. Local Arrays ..
      CHARACTER*132     REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      IF (MSGNP.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) MAJITS
            CALL X04BAY(ISUMM,4,REC)
         END IF
      END IF
C
      IF (MSGNP.GE.5) THEN
C
         MJR = MOD(MAJITS,1000)
         MNR = MOD(MINITS,1000)
         NEVAL = MOD(NFUN,1000)
         NDF = MOD(NZ,1000)
         NLNCON = NCNLN .GT. 0
         FIRST = MAJITS .EQ. MAJIT0
C
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. 50000
            PRTHDR = MSGQP .GT. 0 .OR. FIRST .OR. MSGNP .GE. 20 .OR.
     *               NEWSET
C
            IF (PRTHDR) THEN
               IF (NLNCON) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAY(ISUMM,3,REC)
               ELSE
                  WRITE (REC,FMT=99996)
                  CALL X04BAY(ISUMM,3,REC)
               END IF
               LINES1 = 0
            END IF
C
            IF (NLNCON) THEN
               WRITE (REC,FMT=99997) MJR, MNR, ALFA, NEVAL, OBJALF,
     *           CVNORM, GZNORM, NDF, N - NFREE, LINACT, NLNACT,
     *           SCALE*RHONRM, GFNORM, CONDH, CONDHZ, CONDT, CONVRG,
     *           KTCOND(1), KTCOND(2), LSUMRY
               CALL X04BAF(ISUMM,REC(1))
            ELSE
               WRITE (REC,FMT=99995) MJR, MNR, ALFA, NEVAL, OBJALF,
     *           GZNORM, NDF, N - NFREE, LINACT, GFNORM, CONDH, CONDHZ,
     *           CONDT, CONVRG, KTCOND(1), KTCOND(2), LSUMRY
               CALL X04BAF(ISUMM,REC(1))
            END IF
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. 50000
            PRTHDR = MSGQP .GT. 0 .OR. FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               IF (NLNCON) THEN
                  WRITE (REC,FMT=99994)
                  CALL X04BAY(IPRINT,3,REC)
               ELSE
                  WRITE (REC,FMT=99992)
                  CALL X04BAY(IPRINT,3,REC)
               END IF
               LINES2 = 0
            END IF
C
            IF (NLNCON) THEN
               WRITE (REC,FMT=99993) MJR, MNR, ALFA, OBJALF, CVNORM,
     *           GZNORM, CONDHZ, LSUMRY
               CALL X04BAF(IPRINT,REC(1))
            ELSE
               WRITE (REC,FMT=99991) MJR, MNR, ALFA, OBJALF, GZNORM,
     *           CONDHZ, LSUMRY
               CALL X04BAF(IPRINT,REC(1))
            END IF
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGNP.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               IF (NCNLN.EQ.0) THEN
                  WRITE (REC,FMT=99990) OBJF
                  CALL X04BAY(ISUMM,2,REC)
               ELSE
                  CVIOLS = DNRM2(NCNLN,VIOLN,1)
                  WRITE (REC,FMT=99989) OBJF, CVIOLS
                  CALL X04BAY(ISUMM,2,REC)
               END IF
C
C              ---------------------------------------------------------
C              Print the constraint values.
C              ---------------------------------------------------------
               WRITE (REC,FMT=99988)
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99987)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 I = 1, N, 5
                  WRITE (REC,FMT=99982) (X(J),ISTATE(J),J=I,MIN(I+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99986)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 I = 1, NCLIN, 5
                     WRITE (REC,FMT=99982) (AX(K),ISTATE(N+K),K=I,
     *                 MIN(I+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
               IF (NCNLN.GT.0) THEN
                  WRITE (REC,FMT=99985)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 60 I = 1, NCNLN, 5
                     WRITE (REC,FMT=99982) (C(K),ISTATE(N+NCLIN+K),K=I,
     *                 MIN(I+4,NCNLN))
                     CALL X04BAF(ISUMM,REC(1))
   60             CONTINUE
               END IF
C
               IF (MSGNP.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  INCT = LDT - 1
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(NACTIV,NZ+1),INCT,WORK,1)
                     WRITE (REC,FMT=99984)
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 I = 1, NACTIV, 5
                        WRITE (REC,FMT=99981) (WORK(J),J=I,
     *                    MIN(I+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
                  WRITE (REC,FMT=99983)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 100 I = 1, N, 5
                     WRITE (REC,FMT=99981) (R(J,J),J=I,MIN(I+4,N))
                     CALL X04BAF(ISUMM,REC(1))
  100             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      IF (MSGNP.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99980)
            CALL X04BAY(ISUMM,3,REC)
         END IF
      END IF
C
      LSUMRY(1:2) = '  '
      LSUMRY(4:5) = '  '
C
      RETURN
C
C
C     End of E04UCT. (NPPRT)
C
99999 FORMAT (//' Major iteration',I5,/' ====================')
99998 FORMAT (//'  Maj  Mnr    Step Nfun  Merit function  Violtn Norm ',
     *       'Gz   Nz  Bnd  Lin  Nln Penalty Norm Gf  Cond H Cond Hz  ',
     *       'Cond T Conv')
99997 FORMAT (2I5,1P,D8.1,I5,D16.8,2D8.1,4I5,5D8.1,1X,L1,1X,2L1,A5)
99996 FORMAT (//'  Maj  Mnr    Step Nfun       Objective Norm Gz   Nz ',
     *       ' Bnd  Lin Norm Gf  Cond H Cond Hz  Cond T Conv')
99995 FORMAT (2I5,1P,D8.1,I5,D16.8,D8.1,3I5,4D8.1,1X,L1,1X,2L1,A5)
99994 FORMAT (//'  Maj  Mnr    Step Merit function  Violtn Norm Gz Con',
     *       'd Hz')
99993 FORMAT (2I5,1P,D8.1,D15.6,3D8.1,2X,A5)
99992 FORMAT (//'  Maj  Mnr    Step      Objective Norm Gz Cond Hz')
99991 FORMAT (2I5,1P,D8.1,D15.6,2D8.1,2X,A5)
99990 FORMAT (/' Nonlinear objective value = ',1P,D15.6)
99989 FORMAT (/' Nonlinear objective value = ',1P,D15.6,'   Norm of th',
     *       'e nonlinear constraint violations = ',D15.6)
99988 FORMAT (/' Values of the constraints and their predicted status',
     *       /' ----------------------------------------------------')
99987 FORMAT (/' Variables                  ')
99986 FORMAT (/' General linear constraints ')
99985 FORMAT (/' Nonlinear constraints      ')
99984 FORMAT (/' Diagonals of  T  =         ')
99983 FORMAT (/' Diagonals of  R  =         ')
99982 FORMAT (1X,5(1P,D15.6,I4))
99981 FORMAT (1P,5D15.6)
99980 FORMAT (' ======================================================',
     *       '===================================',//)
      END
      SUBROUTINE E04NFM(SIDE,N,K1,K2,S,A,LDA)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C  E04NFM applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C  Based on F06QNF.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K1, K2, LDA, N
      CHARACTER*1       SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AIJ, TEMP
      INTEGER           I, J
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF ((MIN(N,K1).LT.1) .OR. (K2.LE.K1) .OR. (K2.GT.N)) RETURN
      IF ((SIDE.EQ.'L') .OR. (SIDE.EQ.'l')) THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 40 J = N, K1, -1
            IF (J.GE.K2) THEN
               AIJ = A(K2,J)
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = ZERO
               S(J) = A(J,J)
            END IF
            DO 20 I = MIN(K2,J) - 1, K1, -1
               TEMP = A(I,J)
               A(I+1,J) = TEMP
               AIJ = AIJ
   20       CONTINUE
            A(K1,J) = AIJ
   40    CONTINUE
      ELSE IF ((SIDE.EQ.'R') .OR. (SIDE.EQ.'r')) THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 80 J = K1, K2 - 1
            DO 60 I = 1, J
               TEMP = A(I,J+1)
               A(I,J+1) = A(I,J)
               A(I,J) = TEMP
   60       CONTINUE
            S(J) = A(J+1,J+1)
            A(J+1,J+1) = ZERO
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of E04NFM.
C
      END
      SUBROUTINE E04UCJ(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,NOUT,
     *                  ALFMAX,ALFSML,EPSAF,G0,TARGTG,FTRY,TOLABS,
     *                  TOLREL,TOLTNY,ALFA,ALFBST,FBEST)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-997 (JUN 1993).
C
C     ==================================================================
C     E04UCJ  finds a sequence of improving estimates of a minimizer of
C     the univariate function f(alpha) in the interval (0,ALFMAX].
C     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
C     E04UCJ  requires  f(alpha) (but not f'(alpha)) to be evaluated
C     in the interval.  New estimates of a minimizer are computed using
C     safeguarded quadratic interpolation.
C
C     Reverse communication is used to allow the calling program to
C     evaluate f.  Some of the parameters must be set or tested by the
C     calling program.  The remainder would ordinarily be local
C     variables.
C
C     Input parameters (relevant to the calling program)
C     --------------------------------------------------
C
C     FIRST         must be .TRUE. on the first entry.
C                   It is subsequently altered by E04UCJ.
C
C     DEBUG         specifies whether detailed output is wanted.
C
C     MAXF          is an upper limit on the number of times E04UCJ is
C                   to be entered consecutively with DONE = .FALSE.
C                   (following an initial entry with FIRST = .TRUE.).
C
C     ALFA          is the first estimate of a minimizer.  ALFA is
C                   subsequently altered by E04UCJ (see below).
C
C     ALFMAX        is the upper limit of the interval to be searched.
C
C     ALFSML        is intended to prevent inefficiency when a minimizer
C                   is very small, for cases where the calling program
C                   would prefer to redefine f'(ALFA).  ALFSML is
C                   allowed to be zero.  Early termination will occur if
C                   E04UCJ determines that a minimizer lies somewhere in
C                   the interval [0, ALFSML) (but not if ALFMAX is
C                   smaller that ALFSML).
C
C     EPSAF         is an estimate of the absolute precision in the
C                   computed value of f(0).
C
C     FTRY          the value of f at the new point
C                   ALFA = ALFBST + XTRY.
C
C     G0            is the value of f'(0).  G0 must be negative.
C
C     TOLABS,TOLREL define a function TOL(ALFA) = TOLREL*ALFA + TOLABS
C                   such that if f has already been evaluated at ALFA,
C                   it will not be evaluated closer than TOL(ALFA).
C                   These values may be reduced by E04UCK.
C
C     TARGTG        is the target value of abs(f'(ALFA)). The search
C                   is terminated when
C                    abs(f'(ALFA)) le TARGTG and f(ALFA) lt 0.
C
C     TOLTNY        is the smallest value that TOLABS is allowed to be
C                   reduced to.
C
C     Output parameters (relevant to the calling program)
C     ---------------------------------------------------
C
C     IMPRVD        is .TRUE. if the previous ALFA was the best point so
C                   far.  Any related quantities should be saved by the
C                   calling program (e.g., arrays) before paying
C                   attention to the variable DONE.
C
C     DONE = .FALSE.  means the calling program should evaluate FTRY
C                   for the new trial step ALFA, and reenter E04UCJ.
C
C     DONE = .TRUE.   means that no new ALFA was calculated.  The value
C                   of INFORM gives the result of the search as follows
C
C                   INFORM = 1 means the search has terminated
C                              successfully with ALFBST < ALFMAX.
C
C                   INFORM = 2 means the search has terminated
C                              successfully with ALFBST = ALFMAX.
C
C                   INFORM = 3 means that the search failed to find a
C                              point of sufficient decrease in MAXF
C                              functions, but a lower point was found.
C
C                   INFORM = 4 means ALFMAX is so small that a search
C                              should not have been attempted.
C
C                   INFORM = 5 means that the search was terminated
C                              because of ALFSML (see above).
C
C                   INFORM = 6 means the search has failed to find a
C                              useful step.  The interval of uncertainty
C                              is [0,B] with B < 2*TOLABS. A minimizer
C                              lies very close to ALFA = 0, or f'(0) is
C                              not sufficiently accurate.
C
C                   INFORM = 7 if no better point could be found after
C                              MAXF  function calls.
C
C                   INFORM = 8 means the input parameters were bad.
C                              ALFMAX le TOLTNY  or  G0 ge zero.
C                              No function evaluations were made.
C
C     NUMF          counts the number of times E04UCJ has been entered
C                   consecutively with DONE = .FALSE. (i.e., with a new
C                   function value FTRY).
C
C     ALFA          is the point at which the next function FTRY must
C                   be computed.
C
C     ALFBST        should be accepted by the calling program as the
C                   approximate minimizer, whenever E04UCJ returns
C                   INFORM = 1, 2 or 3.
C
C     FBEST         will be the corresponding value of f.
C
C     The following parameters retain information between entries
C     -----------------------------------------------------------
C
C     BRAKTD        is .FALSE. if f has not been evaluated at the far
C                   end of the interval of uncertainty.  In this case,
C                   the point B will be at ALFMAX + TOL(ALFMAX).
C
C     CRAMPD        is .TRUE. if ALFMAX is very small (le TOLABS).  If
C                   the search fails, this indicates that a zero step
C                   should be taken.
C
C     EXTRAP        is .TRUE. if ALFBST has MOVED at least once and XV
C                   lies outside the interval of uncertainty.  In this
C                   case, extra safeguards are applied to allow for
C                   instability in the polynomial fit.
C
C     MOVED         is .TRUE. if a better point has been found, i.e.,
C                   ALFBST gt 0.
C
C     VSET          records whether a third-best point has been defined.
C
C     WSET          records whether a second-best point has been
C                   defined.  It will always be .TRUE. by the time the
C                   convergence test is applied.
C
C     NSAMEA        is the number of consecutive times that the
C                   left-hand end point of the interval of uncertainty
C                   has remained the same.
C
C     NSAMEB        similarly for the right-hand end.
C
C     A, B, ALFBST  define the current interval of uncertainty.
C                   A minimizer lies somewhere in the  interval
C                   [ALFBST + A, ALFBST + B].
C
C     ALFBST        is the best point so far.  It lies strictly within
C                   [ATRUE,BTRUE]  (except when ALFBST has not been
C                   MOVED, in which case it lies at the left-hand end
C                   point).  Hence we have A .le. 0 and B .gt. 0.
C
C     FBEST         is the value of f at the point ALFBST.
C
C     FA            is the value of f at the point ALFBST + A.
C
C     FACTOR        controls the rate at which extrapolated estimates of
C                   ALFA  may expand into the interval of uncertainty.
C                   FACTOR is not used if a minimizer has been bracketed
C                   (i.e., when the variable BRAKTD is .TRUE.).
C
C     FV, FW        are the values of f at the points ALFBST + XV  and
C                   ALFBST + XW.  They are not defined until  VSET  or
C                   WSET  are .TRUE..
C
C     XTRY          is the trial point within the shifted interval
C                   (A, B).  The new trial function value must be
C                   computed at the point ALFA = ALFBST + XTRY.
C
C     XV            is such that ALFBST + XV is the third-best point.
C                   It is not defined until VSET is .TRUE..
C
C     XW            is such that ALFBST + XW is the second-best point.
C                   It is not defined until WSET is .TRUE..  In some
C                   cases,  XW will replace a previous XW that has a
C                   lower function but has just been excluded from
C                   (A,B).
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version February 1982.  Rev. May 1983.
C     Original F77 version 22-August-1985.
C     This version of E04UCJ dated  24-Oct-91.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT1, HALF
      PARAMETER         (ZERO=0.0D+0,POINT1=0.1D+0,HALF=0.5D+0)
      DOUBLE PRECISION  ONE, TWO, FIVE
      PARAMETER         (ONE=1.0D+0,TWO=2.0D+0,FIVE=5.0D+0)
      DOUBLE PRECISION  TEN, ELEVEN
      PARAMETER         (TEN=1.0D+1,ELEVEN=1.1D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBST, ALFMAX, ALFSML, EPSAF, FBEST,
     *                  FTRY, G0, TARGTG, TOLABS, TOLREL, TOLTNY
      INTEGER           INFORM, MAXF, NOUT, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ARTIFA, ARTIFB, B, DAUX, DTRY, ENDPNT, FA,
     *                  FACTOR, FV, FW, GV, GW, Q, S, TOL, TOLMAX,
     *                  TRUEA, TRUEB, XMIDPT, XTRY, XV, XW
      INTEGER           NSAMEA, NSAMEB
      LOGICAL           BADFUN, BRAKTD, CLOSEF, CRAMPD, EXTRAP, FOUND,
     *                  MOVED, QUITF, QUITFZ, QUITI, QUITS, SETXV, VSET,
     *                  WSET, XINXW
C     .. Local Arrays ..
      CHARACTER*120     REC(7)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Save statement ..
      SAVE              BRAKTD, CRAMPD, EXTRAP, MOVED, VSET, WSET,
     *                  NSAMEA, NSAMEB, A, B, FA, FACTOR, XTRY, XW, FW,
     *                  XV, FV, TOLMAX
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Local variables
C     ===============
C
C     CLOSEF     is .TRUE. if the worst function FV is within EPSAF of
C                FBEST (up or down).
C
C     FOUND      is .TRUE. if the sufficient decrease conditions holds
C                at ALFBST.
C
C     QUITF      is .TRUE. when  MAXF  function calls have been made.
C
C     QUITFZ     is .TRUE. when the three best function values are
C                within EPSAF of each other, and the new point satisfies
C                FBEST le FTRY le FBEST+EPSAF.
C
C     QUITI      is .TRUE. when the interval of uncertainty is less than
C                2*TOL.
C
C     QUITS      is .TRUE. as soon as ALFA is too small to be useful;
C                i.e., BTRUE le ALFSML.
C
C     XINXW      is .TRUE. if XTRY is in (XW,0) or (0,XW).
C     ------------------------------------------------------------------
C
      IMPRVD = .FALSE.
      BADFUN = .FALSE.
      QUITF = .FALSE.
      QUITFZ = .FALSE.
      QUITS = .FALSE.
      QUITI = .FALSE.
C
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First entry.  Initialize various quantities, check input data
C        and prepare to evaluate the function at the initial step ALFA.
C        ---------------------------------------------------------------
         FIRST = .FALSE.
         NUMF = 0
         ALFBST = ZERO
         BADFUN = ALFMAX .LE. TOLTNY .OR. G0 .GE. ZERO
         DONE = BADFUN
         MOVED = .FALSE.
C
         IF ( .NOT. DONE) THEN
            BRAKTD = .FALSE.
            CRAMPD = ALFMAX .LE. TOLABS
            EXTRAP = .FALSE.
            VSET = .FALSE.
            WSET = .FALSE.
            NSAMEA = 0
            NSAMEB = 0
C
            TOLMAX = TOLREL*ALFMAX + TOLABS
            A = ZERO
            B = ALFMAX + TOLMAX
            FA = ZERO
            FACTOR = FIVE
            TOL = TOLABS
            XTRY = ALFA
            IF (DEBUG) THEN
               WRITE (REC,FMT=99999) G0, TOLABS, ALFMAX, TARGTG, TOLREL,
     *           EPSAF, CRAMPD
               CALL X04BAY(NOUT,4,REC)
            END IF
         END IF
      ELSE
C        ---------------------------------------------------------------
C        Subsequent entries.  The function has just been evaluated at
C        ALFA = ALFBST + XTRY,  giving FTRY.
C        ---------------------------------------------------------------
         IF (DEBUG) THEN
            WRITE (REC,FMT=99998) ALFA, FTRY
            CALL X04BAY(NOUT,2,REC)
         END IF
C
         NUMF = NUMF + 1
         NSAMEA = NSAMEA + 1
         NSAMEB = NSAMEB + 1
C
         IF ( .NOT. BRAKTD) THEN
            TOLMAX = TOLABS + TOLREL*ALFMAX
            B = ALFMAX - ALFBST + TOLMAX
         END IF
C
C        Check if XTRY is in the interval (XW,0) or (0,XW).
C
         IF (WSET) THEN
            XINXW = ZERO .LT. XTRY .AND. XTRY .LE. XW .OR. XW .LE.
     *              XTRY .AND. XTRY .LT. ZERO
         ELSE
            XINXW = .FALSE.
         END IF
C
         IMPRVD = FTRY .LT. FBEST
         IF (VSET) THEN
            CLOSEF = ABS(FBEST-FV) .LE. EPSAF
         ELSE
            CLOSEF = .FALSE.
         END IF
C
         IF (IMPRVD) THEN
C
C           We seem to have an improvement.  The new point becomes the
C           origin and other points are shifted accordingly.
C
            IF (WSET) THEN
               XV = XW - XTRY
               FV = FW
               VSET = .TRUE.
            END IF
C
            XW = ZERO - XTRY
            FW = FBEST
            WSET = .TRUE.
            FBEST = FTRY
            ALFBST = ALFA
            MOVED = .TRUE.
C
            A = A - XTRY
            B = B - XTRY
            EXTRAP = .NOT. XINXW
C
C           Decrease the length of (A,B).
C
            IF (XTRY.GE.ZERO) THEN
               A = XW
               FA = FW
               NSAMEA = 0
            ELSE
               B = XW
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
         ELSE IF (CLOSEF .AND. FTRY-FBEST.LT.EPSAF) THEN
C
C           Quit if there has been no progress and FTRY, FBEST, FW
C           and FV are all within EPSAF of each other.
C
            QUITFZ = .TRUE.
         ELSE
C
C           The new function value is no better than the current best
C           point.  XTRY must an end point of the new (A,B).
C
            IF (XTRY.LT.ZERO) THEN
               A = XTRY
               FA = FTRY
               NSAMEA = 0
            ELSE
               B = XTRY
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
C
C           The origin remains unchanged but XTRY may qualify as XW.
C
            IF (WSET) THEN
               IF (FTRY.LT.FW) THEN
                  XV = XW
                  FV = FW
                  VSET = .TRUE.
C
                  XW = XTRY
                  FW = FTRY
                  IF (MOVED) EXTRAP = XINXW
               ELSE IF (MOVED) THEN
                  IF (VSET) THEN
                     SETXV = FTRY .LT. FV .OR. .NOT. EXTRAP
                  ELSE
                     SETXV = .TRUE.
                  END IF
C
                  IF (SETXV) THEN
                     IF (VSET .AND. XINXW) THEN
                        XW = XV
                        FW = FV
                     END IF
                     XV = XTRY
                     FV = FTRY
                     VSET = .TRUE.
                  END IF
               ELSE
                  XW = XTRY
                  FW = FTRY
               END IF
            ELSE
               XW = XTRY
               FW = FTRY
               WSET = .TRUE.
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Check the termination criteria.
C        ---------------------------------------------------------------
         TOL = TOLABS + TOLREL*ALFBST
         TRUEA = ALFBST + A
         TRUEB = ALFBST + B
C
         FOUND = MOVED .AND. ABS(FA-FBEST) .LE. -A*TARGTG
         QUITF = NUMF .GE. MAXF
         QUITI = B - A .LE. TOL + TOL
         QUITS = TRUEB .LE. ALFSML
C
         IF (QUITI .AND. .NOT. MOVED) THEN
C
C           The interval of uncertainty appears to be small enough,
C           but no better point has been found.  Check that changing
C           ALFA by B-A changes f by less than EPSAF.
C
            TOL = TOL/TEN
            TOLABS = TOL
            QUITI = ABS(FW) .LE. EPSAF .OR. TOL .LE. TOLTNY
         END IF
C
         DONE = QUITF .OR. QUITFZ .OR. QUITS .OR. QUITI .OR. FOUND
C
         IF (DEBUG) THEN
            WRITE (REC,FMT=99997) TRUEA, TRUEB, B - A, TOL, NSAMEA,
     *        NSAMEB, NUMF, BRAKTD, EXTRAP, CLOSEF, IMPRVD, FOUND,
     *        QUITI, QUITFZ, QUITS, ALFBST, FBEST, ALFBST + XW, FW
            CALL X04BAY(NOUT,7,REC)
            IF (VSET) THEN
               WRITE (REC,FMT=99996) ALFBST + XV, FV
               CALL X04BAY(NOUT,2,REC)
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Proceed with the computation of an estimate of a minimizer.
C        The choices are...
C        1. Parabolic fit using function values only.
C        2. Damped parabolic fit if the regular fit appears to be
C           consistently overestimating the distance to a minimizer.
C        3. Bisection, geometric bisection, or a step of TOL if the
C           parabolic fit is unsatisfactory.
C        ---------------------------------------------------------------
         IF ( .NOT. DONE) THEN
            XMIDPT = HALF*(A+B)
            S = ZERO
            Q = ZERO
C
C           ============================================================
C           Fit a parabola.
C           ============================================================
C           See if there are two or three points for the parabolic fit.
C
            GW = (FW-FBEST)/XW
            IF (VSET .AND. MOVED) THEN
C
C              Three points available.  Use FBEST, FW and FV.
C
               GV = (FV-FBEST)/XV
               S = GV - (XV/XW)*GW
               Q = TWO*(GV-GW)
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99994)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            ELSE
C
C              Only two points available.  Use FBEST, FW and G0.
C
               IF (MOVED) THEN
                  S = G0 - TWO*GW
               ELSE
                  S = G0
               END IF
               Q = TWO*(G0-GW)
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99995)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            END IF
C
C           ------------------------------------------------------------
C           Construct an artificial interval (ARTIFA, ARTIFB) in which
C           the new estimate of the steplength must lie.  Set a default
C           value of  XTRY  that will be used if the polynomial fit is
C           rejected. In the following, the interval (A,B) is considered
C           the sum of two intervals of lengths  DTRY  and  DAUX, with
C           common end point the best point (zero).  DTRY is the length
C           of the interval into which the default XTRY will be placed
C           and ENDPNT denotes its non-zero end point.  The magnitude of
C           XTRY is computed so that the exponents of DTRY and DAUX are
C           approximately bisected.
C           ------------------------------------------------------------
            ARTIFA = A
            ARTIFB = B
            IF ( .NOT. BRAKTD) THEN
C
C              A minimizer has not yet been bracketed.
C              Set an artificial upper bound by expanding the interval
C              XW  by a suitable FACTOR.
C
               XTRY = -FACTOR*XW
               ARTIFB = XTRY
               IF (ALFBST+XTRY.LT.ALFMAX) FACTOR = FIVE*FACTOR
            ELSE IF (VSET .AND. MOVED) THEN
C
C              Three points exist in the interval of uncertainty.
C              Check if the points are configured for an extrapolation
C              or an interpolation.
C
               IF (EXTRAP) THEN
C
C                 The points are configured for an extrapolation.
C
                  IF (XW.LT.ZERO) ENDPNT = B
                  IF (XW.GT.ZERO) ENDPNT = A
               ELSE
C
C                 If the interpolation appears to be overestimating the
C                 distance to a minimizer,  damp the interpolation step.
C
                  IF (NSAMEA.GE.3 .OR. NSAMEB.GE.3) THEN
                     FACTOR = FACTOR/FIVE
                     S = FACTOR*S
                  ELSE
                     FACTOR = ONE
                  END IF
C
C                 The points are configured for an interpolation.  The
C                 artificial interval will be just (A,B).  Set ENDPNT so
C                 that XTRY lies in the larger of the intervals (A,B)
C                 and  (0,B).
C
                  IF (XMIDPT.GT.ZERO) THEN
                     ENDPNT = B
                  ELSE
                     ENDPNT = A
                  END IF
C
C                 If a bound has remained the same for three iterations,
C                 set ENDPNT so that  XTRY  is likely to replace the
C                 offending bound.
C
                  IF (NSAMEA.GE.3) ENDPNT = A
                  IF (NSAMEB.GE.3) ENDPNT = B
               END IF
C
C              Compute the default value of  XTRY.
C
               DTRY = ABS(ENDPNT)
               DAUX = B - A - DTRY
               IF (DAUX.GE.DTRY) THEN
                  XTRY = FIVE*DTRY*(POINT1+DTRY/DAUX)/ELEVEN
               ELSE
                  XTRY = HALF*SQRT(DAUX)*SQRT(DTRY)
               END IF
               IF (ENDPNT.LT.ZERO) XTRY = -XTRY
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99992) XTRY, DAUX, DTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
C
C              If the points are configured for an extrapolation set the
C              artificial bounds so that the artificial interval lies
C              within (A,B).  If the polynomial fit is rejected,  XTRY
C              will remain at the relevant artificial bound.
C
               IF (EXTRAP) THEN
                  IF (XTRY.LE.ZERO) THEN
                     ARTIFA = XTRY
                  ELSE
                     ARTIFB = XTRY
                  END IF
               END IF
            ELSE
C
C              The gradient at the origin is being used for the
C              polynomial fit.  Set the default XTRY to one tenth XW.
C
               IF (EXTRAP) THEN
                  XTRY = -XW
               ELSE
                  XTRY = XW/TEN
               END IF
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99993) XTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
            END IF
C
C           ------------------------------------------------------------
C           The polynomial fits give (S/Q)*XW as the new step.  Reject
C           this step if it lies outside (ARTIFA, ARTIFB).
C           ------------------------------------------------------------
            IF (Q.NE.ZERO) THEN
               IF (Q.LT.ZERO) S = -S
               IF (Q.LT.ZERO) Q = -Q
               IF (S*XW.GE.Q*ARTIFA .AND. S*XW.LE.Q*ARTIFB) THEN
C
C                 Accept the polynomial fit.
C
                  IF (ABS(S*XW).GE.Q*TOL) THEN
                     XTRY = (S/Q)*XW
                  ELSE
                     XTRY = ZERO
                  END IF
                  IF (DEBUG) THEN
                     WRITE (REC,FMT=99991) XTRY
                     CALL X04BAF(NOUT,REC(1))
                  END IF
               END IF
            END IF
         END IF
      END IF
C     ==================================================================
C
      IF ( .NOT. DONE) THEN
         ALFA = ALFBST + XTRY
         IF (BRAKTD .OR. ALFA.LT.ALFMAX-TOLMAX) THEN
C
C           The function must not be evaluated too close to A or B.
C           (It has already been evaluated at both those points.)
C
            XMIDPT = HALF*(A+B)
            IF (XTRY.LE.A+TOL .OR. XTRY.GE.B-TOL) THEN
               IF (XMIDPT.LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
            END IF
C
            IF (ABS(XTRY).LT.TOL) THEN
               IF (XMIDPT.LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
            END IF
            ALFA = ALFBST + XTRY
         ELSE
C
C           The step is close to or larger than ALFMAX, replace it by
C           ALFMAX to force evaluation of the function at the boundary.
C
            BRAKTD = .TRUE.
            XTRY = ALFMAX - ALFBST
            ALFA = ALFMAX
         END IF
      END IF
C     ------------------------------------------------------------------
C     Exit.
C     ------------------------------------------------------------------
      IF (DONE) THEN
         IF (BADFUN) THEN
            INFORM = 8
         ELSE IF (QUITS) THEN
            INFORM = 5
         ELSE IF (FOUND) THEN
            IF (ALFBST.LT.ALFMAX) THEN
               INFORM = 1
            ELSE
               INFORM = 2
            END IF
         ELSE IF (MOVED) THEN
            INFORM = 3
         ELSE IF (QUITF) THEN
            INFORM = 7
         ELSE IF (CRAMPD) THEN
            INFORM = 4
         ELSE
            INFORM = 6
         END IF
      END IF
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99990)
         CALL X04BAY(NOUT,2,REC)
      END IF
      RETURN
C
C
C     End of  E04UCJ. (SRCHQ)
C
99999 FORMAT (/'     G0  TOLABS  ALFMAX        ',1P,2D22.14,D16.8,/' T',
     *       'ARGTG  TOLREL   EPSAF        ',1P,2D22.14,D16.8,/' CRAMP',
     *       'D                        ',L3)
99998 FORMAT (/' ALFA    FTRY                  ',1P,2D22.14)
99997 FORMAT (/' A       B       B - A   TOL   ',1P,2D22.14,2D16.8,
     *       /' NSAMEA  NSAMEB  NUMF          ',3I3,/' BRAKTD  EXTRAP ',
     *       ' CLOSEF  IMPRVD',4L3,/' FOUND   QUITI   QUITFZ  QUITS ',
     *       4L3,/' ALFBST  FBEST                 ',1P,2D22.14,/' ALFA',
     *       'W   FW                    ',1P,2D22.14)
99996 FORMAT (' ALFAV   FV                    ',1P,2D22.14,/)
99995 FORMAT (' Parabolic fit,    two points. ')
99994 FORMAT (' Parabolic fit,  three points. ')
99993 FORMAT (' Exponent reduced.  Trial point',1P,D22.14)
99992 FORMAT (' Geo. bisection. XTRY,DAUX,DTRY',1P,3D22.14)
99991 FORMAT (' Polynomial fit accepted.  XTRY',1P,D22.14)
99990 FORMAT (' ----------------------------------------------------',/)
      END

      SUBROUTINE E04MFM(PRBTYP,MSGLVL,N,LDA,LDT,NACTIV,NFREE,NZ,ISTATE,
     *                  KACTIV,KX,ZEROLM,NOTOPT,NUMINF,TRUSML,SMLLST,
     *                  JSMLST,KSMLST,TINYST,JTINY,JINF,TRUBIG,BIGGST,
     *                  JBIGST,KBIGST,A,ANORMS,GQ,RLAMDA,T,WTINF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFM  first computes the Lagrange multiplier estimates for the
C     given working set.  It then determines the values and indices of
C     certain significant multipliers.  In this process, the multipliers
C     for inequalities at their upper bounds are adjusted so that a
C     negative multiplier for an inequality constraint indicates non-
C     optimality.  All adjusted multipliers are scaled by the 2-norm
C     of the associated constraint row.  In the following, the term
C     minimum refers to the ordering of numbers on the real line,  and
C     not to their magnitude.
C
C     JSMLST          is the index of the constraint whose multiplier is
C                     the minimum of the set of adjusted multipliers
C                     with values less than  small.
C     RLAMDA(KSMLST)  is the associated multiplier.
C
C     JBIGST          is the index of the constraint whose multiplier is
C                     the largest of the set of adjusted multipliers
C                     with values greater than (1 + small).
C     RLAMDA(KBIGST)  is the associated multiplier.
C
C     On exit,  elements  1  thru  NACTIV  of  RLAMDA  contain the
C     unadjusted multipliers for the general constraints.  Elements
C     NACTIV  onwards of  RLAMDA  contain the unadjusted multipliers
C     for the bounds.
C
C     Original version written 31-October-1984.
C     Based on a version of  LSMULS  dated 30-June-1986.
C     This version of  E04MFM  dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGGST, SMLLST, TINYST, TRUBIG, TRUSML, ZEROLM
      INTEGER           JBIGST, JINF, JSMLST, JTINY, KBIGST, KSMLST,
     *                  LDA, LDT, MSGLVL, N, NACTIV, NFREE, NOTOPT,
     *                  NUMINF, NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ANORMS(*), GQ(N), RLAMDA(N), T(LDT,*),
     *                  WTINF(*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LCDBG
C     .. Arrays in Common ..
      INTEGER           ILCDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORMJ, BLAM, RLAM, SCDLAM
      INTEGER           I, IS, J, K, KK, L, NFIXED
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTRSV, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /EE04MF/ILCDBG, LCDBG
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
C
      JTINY = 0
      JSMLST = 0
      KSMLST = 0
C
      JBIGST = 0
      KBIGST = 0
C
C     ------------------------------------------------------------------
C     Compute  JSMLST  for regular constraints and temporary bounds.
C     ------------------------------------------------------------------
C     First, compute the Lagrange multipliers for the general
C     constraints in the working set, by solving  T'*lamda = Y'g.
C
      IF (N.GT.NZ) CALL DCOPY(N-NZ,GQ(NZ+1),1,RLAMDA,1)
      IF (NACTIV.GT.0) CALL DTRSV('U','T','N',NACTIV,T(1,NZ+1),LDT,
     *                            RLAMDA,1)
C
C     -----------------------------------------------------------------
C     Now set elements  NACTIV, NACTIV+1,... of  RLAMDA  equal to
C     the multipliers for the bound constraints.
C     -----------------------------------------------------------------
      DO 40 L = 1, NFIXED
         J = KX(NFREE+L)
         BLAM = RLAMDA(NACTIV+L)
         DO 20 K = 1, NACTIV
            I = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(NACTIV-K+1)
   20    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
   40 CONTINUE
C
C     -----------------------------------------------------------------
C     Find  JSMLST  and  KSMLST.
C     -----------------------------------------------------------------
      DO 60 K = 1, N - NZ
         IF (K.GT.NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(NACTIV-K+1) + N
         END IF
C
         IS = ISTATE(J)
C
         I = J - N
         IF (J.LE.N) ANORMJ = ONE
         IF (J.GT.N) ANORMJ = ANORMS(I)
C
         RLAM = RLAMDA(K)
C
C        Change the sign of the estimate if the constraint is in
C        the working set at its upper bound.
C
         IF (IS.EQ.2) RLAM = -RLAM
         IF (IS.EQ.3) RLAM = ABS(RLAM)
         IF (IS.EQ.4) RLAM = -ABS(RLAM)
C
         IF (IS.NE.3) THEN
            SCDLAM = RLAM*ANORMJ
C
            IF (SCDLAM.LT.ZEROLM) THEN
               IF (NUMINF.EQ.0) NOTOPT = NOTOPT + 1
C
               IF (SCDLAM.LT.SMLLST) THEN
                  SMLLST = SCDLAM
                  TRUSML = RLAMDA(K)
                  JSMLST = J
                  KSMLST = K
               END IF
            ELSE IF (SCDLAM.LT.TINYST) THEN
               TINYST = SCDLAM
               JTINY = J
            END IF
         END IF
C
         SCDLAM = RLAM/WTINF(J)
         IF (SCDLAM.GT.BIGGST .AND. J.GT.JINF) THEN
            BIGGST = SCDLAM
            TRUBIG = RLAMDA(K)
            JBIGST = J
            KBIGST = K
         END IF
   60 CONTINUE
C
C     -----------------------------------------------------------------
C     If required, print the multipliers.
C     -----------------------------------------------------------------
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            IF (NFIXED.GT.0) THEN
               WRITE (REC,FMT=99999) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 80 K = 1, NFIXED, 4
                  WRITE (REC,FMT=99998) (KX(NFREE+KK),RLAMDA(NACTIV+KK),
     *              KK=K,MIN(K+3,NFIXED))
                  CALL X04BAF(ISUMM,REC(1))
   80          CONTINUE
            END IF
            IF (NACTIV.GT.0) THEN
               WRITE (REC,FMT=99997) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 100 K = 1, NACTIV, 4
                  WRITE (REC,FMT=99998) (KACTIV(KK),RLAMDA(NACTIV-KK+1),
     *              KK=K,MIN(K+3,NACTIV))
                  CALL X04BAF(ISUMM,REC(1))
  100          CONTINUE
            END IF
         END IF
      END IF
C
      IF (LCDBG .AND. ILCDBG(1).GT.0) THEN
         WRITE (REC,FMT=99996) JSMLST, SMLLST, KSMLST
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99995) JBIGST, BIGGST, KBIGST
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99994) JTINY, TINYST
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      RETURN
C
C     End of  E04MFM.  (CMMUL1)
C
99999 FORMAT (/' Multipliers for the ',A2,' bound  constraints   ')
99998 FORMAT (4(I5,1P,D11.2))
99997 FORMAT (/' Multipliers for the ',A2,' linear constraints   ')
99996 FORMAT (/' //E04MFM//  JSMLST     SMLLST     KSMLST (scaled) ',
     *       /' //E04MFM//  ',I6,1P,D11.2,5X,I6)
99995 FORMAT (' //E04MFM//  JBIGST     BIGGST     KBIGST (scaled) ',
     *       /' //E04MFM//  ',I6,1P,D11.2,5X,I6)
99994 FORMAT (' //E04MFM//   JTINY     TINYST                     ',
     *       /' //E04MFM//  ',I6,1P,D11.2)
      END
      subroutine e04mfx(nout,buffer,key)
c     mark 16 release. nag copyright 1992.
c
c     ******************************************************************
c     e04mfx   decodes the option contained in  buffer  in order to set
c     a parameter value in the relevant element of  iprmlc  or  rprmlc.
c
c
c     input:
c
c     nout   a unit number for printing error messages.
c            nout  must be a valid unit.
c
c     output:
c
c     key    the first keyword contained in buffer.
c
c
c     e04mfx  calls e04udx and the subprograms
c                 lookup, scannr, tokens, upcase
c     (now called e04udy, e04udw, e04udv, e04udu)
c     supplied by informatics general, inc., palo alto, california.
c
c     this version of e04mfx dated 20-apr-93.
c     ******************************************************************
c
c     .. parameters ..
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
c     .. external subroutines ..
      external          e04udv, e04udy, x04baf
c     .. intrinsic functions ..
      intrinsic         index, len
c     .. common blocks ..
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
c     .. data statements ..
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
c     .. executable statements ..
c     ------------------------------------------------------------------
c
      if (first) then
         first = .false.
         do 20 i = 1, mxparm
            iprmlc(i) = idummy
            rprmlc(i) = rdummy
   20    continue
      end if
c
c     eliminate comments and empty lines.
c     a '*' appearing anywhere in buffer terminates the string.
c
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
c
c     ------------------------------------------------------------------
c     extract up to maxtok tokens from the record.
c     ntoken returns how many were actually found.
c     key, key2, key3 are the first tokens if any, otherwise blank.
c     ------------------------------------------------------------------
      ntoken = maxtok
      call e04udv(buffer(1:lenbuf),maxtok,ntoken,token)
      key = token(1)
      key2 = token(2)
      key3 = token(3)
c
c     certain keywords require no action.
c
      if (key.eq.' ' .or. key.eq.'begin') go to 80
      if (key.eq.'list' .or. key.eq.'nolist') go to 80
      if (key.eq.'end') go to 80
c
c     most keywords will have an associated integer or real value,
c     so look for it no matter what the keyword.
c
      i = 1
      number = .false.
c
   40 if (i.lt.ntoken .and. .not. number) then
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
c     ------------------------------------------------------------------
c     decide what to do about each keyword.
c     the second keyword (if any) might be needed to break ties.
c     some seemingly redundant testing of more is used
c     to avoid compiler limits on the number of consecutive else ifs.
c     ------------------------------------------------------------------
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
c
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
c
c              recognize     problem type = lp     etc.
c
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
c
   80 return
c
c
C     End of  E04MFX.  (LPKEY)
C
99999 FORMAT (' XXX  Keyword not recognized:         ',A)
99998 FORMAT (' XXX  Second keyword not recognized:  ',A)
99997 FORMAT (' XXX  Third  keyword not recognized:  ',A)
99996 FORMAT (' XXX  Ambiguous keyword:              ',A)
      END
      SUBROUTINE E04XAX(INFORM,MSGLVL,N,BIGBND,EPSRF,OKTOL,FDCHK,OBJF,
     *                  XNORM,OBJFUN,BL,BU,GRAD,GRADU,DX,X,Y,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-882 (NOV 1990).
C     MARK 16A REVISED. IER-1000 (JUN 1993).
C
C     ******************************************************************
C     E04XAX  checks if the gradients of the objective function have
C     been coded correctly.
C
C     On input,  the value of the objective function at the point X is
C     stored in OBJF.  The corresponding gradient is stored in GRADU.
C     If any gradient element has not been specified,  it will have a
C     dummy value.  Missing values are not checked.
C
C     A cheap test is first undertaken by calculating the directional
C     derivative using two different methods. If this proves
C     satisfactory and no further information is desired, E04XAX is
C     terminated. Otherwise, the routine E04XAZ is called to give
C     optimal step-sizes and a forward-difference approximation to
C     each element of the gradient for which a test is deemed
C     necessary, either by the program or the user.
C
C     Other inputs:
C
C        X         The n-dimensional point at which the
C                  gradient is to be verified.
C        EPSRF     The positive bound on the relative error
C                  associated with computing the function at
C                  the point x.
C        OKTOL     The desired relative accuracy which the
C                  elements of the gradient should satisfy.
C
C     LVRFYC has the following meaning...
C
C     -1        do not perform any check.
C     0        do the cheap test only.
C     1 or 3   do both cheap and full test.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  19-May-1985.
C     This version of E04XAX  dated  19-Oct-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  ZERO, HALF, POINT9
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,POINT9=0.9D+0)
      DOUBLE PRECISION  ONE, TWO, TEN
      PARAMETER         (ONE=1.0D+0,TWO=2.0D+0,TEN=1.0D+1)
      CHARACTER*4       LBAD, LGOOD
      PARAMETER         (LBAD='BAD?',LGOOD='  OK')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, EPSRF, FDCHK, OBJF, OKTOL, XNORM
      INTEGER           INFORM, MSGLVL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), DX(N), GRAD(N), GRADU(N), USER(*),
     *                  X(N), Y(N)
      INTEGER           IUSER(*)
C     .. Subroutine Arguments ..
      EXTERNAL          OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, LVRFYC, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG), JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, DXJ, DXMULT, EMAX, EPSA,
     *                  ERRBND, ERROR, F1, F2, FDEST, GDIFF, GDX, GJ,
     *                  GSIZE, H, HOPT, HPHI, OBJF1, SDEST, STEPBL,
     *                  STEPBU, XJ
      INTEGER           INFO, ITER, ITMAX, J, J1, J2, JMAX, MODE,
     *                  NCHECK, NGOOD, NSTATE, NWRONG
      LOGICAL           CONST, DEBUG, DONE, FIRST, HEADNG, NEEDED, OK
      CHARACTER*4       KEY
C     .. Local Arrays ..
      CHARACTER*18      RESULT(0:4)
      CHARACTER*120     REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, E04XAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Data statements ..
      DATA              RESULT/'                 ', 'Constant?      ',
     *                  'Linear or odd?   ', 'Too nonlinear?',
     *                  'Small derivative?'/
C     .. Executable Statements ..
C
      INFORM = 0
      NEEDED = LVRFYC .EQ. 0 .OR. LVRFYC .EQ. 1 .OR. LVRFYC .EQ. 3
      IF ( .NOT. NEEDED) RETURN
C
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,4,REC)
      END IF
      DEBUG = NPDBG .AND. INPDBG(5) .GT. 0
      NSTATE = 0
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Perform the cheap test.
C     ==================================================================
      H = (ONE+XNORM)*FDCHK
C
      IF (N.LE.100) THEN
         DXMULT = 0.9
      ELSE IF (N.LE.250) THEN
         DXMULT = 0.99
      ELSE
         DXMULT = 0.999
      END IF
C
      DXJ = ONE/N
      DO 20 J = 1, N
         DX(J) = DXJ
         DXJ = -DXJ*DXMULT
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Do not perturb X(J) if the  J-th  element is missing.
C     Compute the directional derivative.
C     ------------------------------------------------------------------
      NCHECK = 0
      DO 40 J = 1, N
         IF (GRAD(J).EQ.RDUMMY) THEN
            DX(J) = ZERO
         ELSE
            NCHECK = NCHECK + 1
C
            XJ = X(J)
            STEPBL = -ONE
            STEPBU = ONE
            IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
            IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J))
     *          STEPBU = MIN(STEPBU,BU(J)-XJ)
C
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
               DX(J) = DX(J)*STEPBL
            ELSE
               DX(J) = DX(J)*STEPBU
            END IF
         END IF
   40 CONTINUE
C
      IF (NCHECK.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99989)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         RETURN
      END IF
      GDX = DDOT(N,GRADU,1,DX,1)
C
C     ------------------------------------------------------------------
C     Make forward-difference approximation along  p.
C     ------------------------------------------------------------------
      CALL DCOPY(N,X,1,Y,1)
      CALL DAXPY(N,H,DX,1,Y,1)
C
      MODE = 0
      CALL OBJFUN(MODE,N,Y,OBJF1,GRADU,NSTATE,IUSER,USER)
      IF (MODE.LT.0) GO TO 100
C
      GDIFF = (OBJF1-OBJF)/H
      ERROR = ABS(GDIFF-GDX)/(ABS(GDX)+ONE)
C
      OK = ERROR .LE. OKTOL
C
      IF (MSGLVL.GT.0) THEN
         IF (OK) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
         ELSE
            WRITE (REC,FMT=99997)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         WRITE (REC,FMT=99996) GDX, GDIFF
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (ERROR.GE.POINT9) INFORM = 1
C
C     ==================================================================
C     Element-wise check.
C     ==================================================================
      IF (LVRFYC.EQ.1 .OR. LVRFYC.EQ.3) THEN
         HEADNG = .TRUE.
         ITMAX = 3
         NCHECK = 0
         NWRONG = 0
         NGOOD = 0
         JMAX = 0
         EMAX = ZERO
         J1 = JVERFY(1)
         J2 = JVERFY(2)
C
C        ---------------------------------------------------------------
C        Loop over each of the elements of  x.
C        ---------------------------------------------------------------
         DO 80 J = J1, J2
C
            IF (GRAD(J).NE.RDUMMY) THEN
C              ---------------------------------------------------------
C              Check this gradient element.
C              ---------------------------------------------------------
               NCHECK = NCHECK + 1
               GJ = GRAD(J)
               GSIZE = ABS(GJ)
               XJ = X(J)
C              ---------------------------------------------------------
C              Find a finite-difference interval by iteration.
C              ---------------------------------------------------------
               ITER = 0
               EPSA = EPSRF*(ONE+ABS(OBJF))
               CDEST = ZERO
               SDEST = ZERO
               FIRST = .TRUE.
C
               STEPBL = BIGLOW
               STEPBU = BIGUPP
               IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
               IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
               HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
               H = TEN*HOPT
               IF (HALF*(STEPBL+STEPBU).LT.ZERO) H = -H
C
C              +             REPEAT
   60          X(J) = XJ + H
               CALL OBJFUN(MODE,N,X,F1,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               X(J) = XJ + H + H
               CALL OBJFUN(MODE,N,X,F2,GRADU,NSTATE,IUSER,USER)
               IF (MODE.LT.0) GO TO 100
C
               CALL E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSRF,OBJF,INFO,ITER,
     *                     ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,H,HOPT,
     *                     HPHI)
C
C              +             UNTIL     DONE
               IF ( .NOT. DONE) GO TO 60
C
C              ---------------------------------------------------------
C              Exit for this variable.
C              ---------------------------------------------------------
               GDIFF = CDEST
               X(J) = XJ
C
               ERROR = ABS(GDIFF-GJ)/(GSIZE+ONE)
               IF (ERROR.GE.EMAX) THEN
                  EMAX = ERROR
                  JMAX = J
               END IF
C
               OK = ERROR .LE. OKTOL
               IF (OK) THEN
                  KEY = LGOOD
                  NGOOD = NGOOD + 1
               ELSE
                  KEY = LBAD
                  NWRONG = NWRONG + 1
               END IF
C
               IF (MSGLVL.GT.0) THEN
C
C                 Zero elements are not printed.
C
                  CONST = OK .AND. INFO .EQ. 1 .AND. ABS(GJ) .LT. EPSPT8
                  IF ( .NOT. CONST) THEN
                     IF (HEADNG) THEN
                        WRITE (REC,FMT=99995)
                        CALL X04BAY(IPRINT,4,REC)
                        HEADNG = .FALSE.
                     END IF
                     IF (OK) THEN
                        WRITE (REC,FMT=99994) J, XJ, HOPT, GJ, GDIFF,
     *                    KEY, ITER
                     ELSE
                        WRITE (REC,FMT=99993) J, XJ, HOPT, GJ, GDIFF,
     *                    KEY, ITER, RESULT(INFO)
                     END IF
                     CALL X04BAF(IPRINT,REC(1))
                  END IF
               END IF
            END IF
   80    CONTINUE
C
C        ===============================================================
C        Done.
C        ===============================================================
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            IF (NWRONG.EQ.0) THEN
               WRITE (REC,FMT=99992) NGOOD, NCHECK, J1, J2
               CALL X04BAY(IPRINT,3,REC)
            ELSE
               WRITE (REC,FMT=99991) NWRONG, NCHECK, J1, J2
               CALL X04BAY(IPRINT,3,REC)
            END IF
            WRITE (REC,FMT=99990) EMAX, JMAX
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (ERROR.GE.POINT9) INFORM = 1
      END IF
C
      CALL DCOPY(N,GRAD,1,GRADU,1)
C
      RETURN
C
  100 INFORM = MODE
      RETURN
C
C
C     End of  E04XAX. (CHKGRD)
C
99999 FORMAT (//' Verification of the objective gradients.',/' -------',
     *       '---------------------------------')
99998 FORMAT (/' The objective gradients seem to be ok.')
99997 FORMAT (/' XXX  The objective gradients seem to be incorrect.')
99996 FORMAT (/' Directional derivative of the objective',1P,D18.8,
     *       /' Difference approximation               ',1P,D18.8)
99995 FORMAT (//4X,'J',4X,'X(J)',5X,'DX(J)',11X,'G(J)',11X,'Difference',
     *       ' approxn  Itns',/)
99994 FORMAT (I5,1P,2D10.2,1P,2D18.8,2X,A4,I6)
99993 FORMAT (I5,1P,2D10.2,1P,2D18.8,2X,A4,I6,2X,A18)
99992 FORMAT (/I7,'  Objective gradients out of the',I6,/9X,'set in co',
     *       'ls',I6,'  through',I6,'  seem to be ok.')
99991 FORMAT (/' XXX  There seem to be',I6,'  incorrect objective grad',
     *       'ients out of the',I6,/8X,'set in cols',I6,'  through',I6)
99990 FORMAT (/' The largest relative error was',1P,D12.2,'   in eleme',
     *       'nt',I6,/)
99989 FORMAT (/' No gradient elements assigned.')
      END
      SUBROUTINE E04UCU(FEASQP,UNITQ,NQPERR,MAJITS,MINITS,N,NCLIN,NCNLN,
     *                  LDCJ,LDAQP,LDR,LINACT,NLNACT,NACTIV,NFREE,NZ,
     *                  NUMINF,ISTATE,KACTIV,KX,DXNORM,GDX,QPCURV,AQP,
     *                  ADX,ANORM,AX,BL,BU,C,CJAC,CLAMDA,CMUL,CS,DLAM,
     *                  DSLK,DX,QPBL,QPBU,QPTOL,R,RHO,SLK,VIOLN,X,WTINF,
     *                  W)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-881 (NOV 1990).
C     MARK 16 REVISED. IER-1089 (JUL 1993).
C
C     ******************************************************************
C     E04UCU   does the following:
C
C     (1)  Generate the upper and lower bounds for the QP  subproblem.
C
C     (2)  Compute the  TQ  factors of the rows of  AQP  specified by
C          the array  ISTATE.  The part of the factorization defined by
C          the first contiguous group of linear constraints does not
C          need to be recomputed.  The remaining rows (which could be
C          comprised of both linear and nonlinear constraints) are
C          included as new rows of the  TQ  factorization stored in
C          T and ZY.  Note that if there are no nonlinear constraints,
C          no factorization is required.
C
C     (3)  Solve the  QP  subproblem.
C                 minimize     1/2 (W p - d)'(Wp - d) + g'p
C
C                 subject to   qpbl .le. (  p ) .le. qpbu,
C                                        ( Ap )
C
C          where  W  is a matrix (not stored) such that  W'W = H  and
C          WQ = R,  d  is the zero vector,  and  g  is the gradient.
C          If the subproblem is infeasible, compute the point which
C          minimizes the sum of infeasibilities.
C
C     (4)   Find the value of each slack variable for which the merit
C          function is minimized.
C
C     (5)   Compute  DSLK,  DLAM  and  DX,  the search directions for
C          the slack variables, the multipliers and the variables.
C
C     Systems Optimization Laboratory, Stanford University.
C     Fortran 66 version written 10-January-1983.
C     Level-2 matrix routines added 18-May-1988.
C     This version of E04UCU dated 23-Dec-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      LOGICAL           QPNAMD, VERTEX
      PARAMETER         (QPNAMD=.FALSE.,VERTEX=.FALSE.)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DXNORM, GDX, QPCURV
      INTEGER           LDAQP, LDCJ, LDR, LINACT, MAJITS, MINITS, N,
     *                  NACTIV, NCLIN, NCNLN, NFREE, NLNACT, NQPERR,
     *                  NUMINF, NZ
      LOGICAL           FEASQP, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), ANORM(*), AQP(LDAQP,*), AX(*), BL(*),
     *                  BU(*), C(*), CJAC(LDCJ,*), CLAMDA(*), CMUL(*),
     *                  CS(*), DLAM(*), DSLK(*), DX(N), QPBL(*),
     *                  QPBU(*), QPTOL(*), R(LDR,*), RHO(*), SLK(*),
     *                  VIOLN(*), W(*), WTINF(*), X(N)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DRMAX, DRMIN, DTMAX, DTMIN, DXLIM, EPSPT3,
     *                  EPSPT5, EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL,
     *                  HCNDBD, RCNDBD, RFROBN, RHODMP, RHOMAX, RHONRM,
     *                  SCALE, TOLACT, TOLFEA, TOLRNK
      INTEGER           IDBGLS, IDBGNP, IPRINT, IPRNT, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, KSAVE, LCRASH, LDBGLS, LDBGNP, LDT,
     *                  LDZY, LENNAM, LFORMH, LINES1, LINES2, LPROB,
     *                  LVERFY, LVLDER, MSGLS, MSGNP, NCOLT, NLNF, NLNJ,
     *                  NLNX, NLOAD, NN, NNCLIN, NNCNLN, NOUT, NPROB,
     *                  NSAVE
      LOGICAL           CMDBG, INCRUN, LSDBG, NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM)
      INTEGER           ICMDBG(LDBG), ILSDBG(LDBG), INPDBG(LDBG),
     *                  IPADLS(18), IPADNP(12), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM), LOCLS(LENLS)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, BIGLOW, BIGUPP, BLJ, BUJ, CON, CONDMX,
     *                  QUOTNT, SSQ, SSQ1, SUMINF, VIOL, WEIGHT, WSCALE,
     *                  WTMAX, WTMIN
      INTEGER           I, IDBG, IDBGSV, INFORM, ISWAP, J, JINF, K, K1,
     *                  K2, KVIOL, L, LGQ, LHPQ, LRLAM, LRPQ, LRPQ0, LT,
     *                  LWRK1, LZY, MJRDBG, MNRDBG, MSGQP, NARTIF, NCQP,
     *                  NCTOTL, NGQ, NMAJOR, NMINOR, NPLIN, NRANK,
     *                  NREJTD, NRPQ, NTRY, NVIOL, NZ1
      LOGICAL           LINOBJ, OVERFL
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, DTRMV, DTRSV,
     *                  E04NBW, E04NCY, E04NCZ, E04UCM, F06DBF, F06DFF,
     *                  F06FBF, F06FLF, F06QFF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IDBGLS, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LPROB, MSGLS, NN,
     *                  NNCLIN, NPROB, IPADLS
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /FE04UC/INPDBG, NPDBG
      COMMON            /GE04UC/IPSVNP, IDBGNP, ITMXNP, JVRFY1, JVRFY2,
     *                  JVRFY3, JVRFY4, LDBGNP, LFORMH, LVLDER, LVERFY,
     *                  MSGNP, NLNF, NLNJ, NLNX, NNCNLN, NSAVE, NLOAD,
     *                  KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IDBGLS), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),IDBGNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IDBGNP,IDBG), (ITMXNP,NMAJOR), (ITMAX2,NMINOR)
      EQUIVALENCE       (LDBGLS,MNRDBG), (LDBGNP,MJRDBG), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/
C     .. Executable Statements ..
C
      IDBGSV = IDBG
      IF (NPDBG) THEN
         IDBG = 0
      ELSE
         IDBG = NMINOR + 1
      END IF
      LSDBG = NPDBG
      CMDBG = NPDBG
      CALL F06DFF(LDBG,ILSDBG,1,ICMDBG,1)
C
      LRPQ = LOCLS(5)
      LRPQ0 = LOCLS(6)
      LHPQ = LOCLS(8)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWRK1 = LOCLS(14)
C
      NRPQ = 0
      NGQ = 1
C
      FEASQP = .TRUE.
      LINOBJ = .TRUE.
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
      SSQ1 = ZERO
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
      NCQP = NCLIN + NCNLN
      NRANK = N
      NREJTD = 0
C
      IF (MSGQP.GT.0) THEN
         WRITE (REC,FMT=99999) MAJITS
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
C     ==================================================================
C     Generate the upper and lower bounds upon the search direction, the
C     weights on the sum of infeasibilities and the nonlinear constraint
C     violations.
C     ==================================================================
      WSCALE = -ONE
      DO 40 J = 1, NCTOTL
C
         IF (J.LE.N) THEN
            CON = X(J)
         ELSE IF (J.LE.NPLIN) THEN
            CON = AX(J-N)
         ELSE
            CON = C(J-NPLIN)
         END IF
C
         BLJ = BL(J)
         BUJ = BU(J)
         IF (BLJ.GT.BIGLOW) BLJ = BLJ - CON
         IF (BUJ.LT.BIGUPP) BUJ = BUJ - CON
C
         WEIGHT = ONE
         IF (J.LE.NPLIN) THEN
            IF (ABS(BLJ).LE.QPTOL(J)) BLJ = ZERO
            IF (ABS(BUJ).LE.QPTOL(J)) BUJ = ZERO
         ELSE
            I = J - NPLIN
            VIOL = ZERO
            IF (BL(J).GT.BIGLOW) THEN
               IF (BLJ.GT.ZERO) THEN
                  VIOL = BLJ
                  IF (RHO(I).GT.ZERO) THEN
                     WEIGHT = VIOL*RHO(I)
                  ELSE
                     WEIGHT = VIOL
                  END IF
                  WSCALE = MAX(WSCALE,WEIGHT)
                  GO TO 20
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               IF (BUJ.LT.ZERO) THEN
                  VIOL = BUJ
                  IF (RHO(I).GT.ZERO) THEN
                     WEIGHT = -VIOL*RHO(I)
                  ELSE
                     WEIGHT = -VIOL
                  END IF
                  WSCALE = MAX(WSCALE,WEIGHT)
               END IF
            END IF
C
C           Set the vector of nonlinear constraint violations.
C
   20       VIOLN(I) = VIOL
         END IF
C
         WTINF(J) = WEIGHT
         QPBL(J) = BLJ
         QPBU(J) = BUJ
C
   40 CONTINUE
C
      IF (WSCALE.GT.ZERO) THEN
         WSCALE = ONE/WSCALE
         CALL DSCAL(NCTOTL,(WSCALE),WTINF,1)
      END IF
C
      CALL F06FLF(NCTOTL,WTINF,1,WTMAX,WTMIN)
      WTMIN = EPSPT9*WTMAX
      DO 60 J = 1, NCTOTL
         WTINF(J) = MAX(WTINF(J),WTMIN)
   60 CONTINUE
C
C     Set the maximum allowable condition estimator of the constraints
C     in the working set.  Note that a relatively well-conditioned
C     working set is used to start the QP iterations.
C
      CONDMX = MAX(ONE/EPSPT3,HUNDRD)
C
      IF (NCNLN.GT.0) THEN
C        ===============================================================
C        Refactorize part of the  QP  constraint matrix.
C        ===============================================================
C        Load the new Jacobian into the  QP  matrix  A.  Compute the
C        2-norms of the rows of the Jacobian.
C
         CALL F06QFF('General',NCNLN,N,CJAC,LDCJ,AQP(NCLIN+1,1),LDAQP)
C
         DO 80 J = NCLIN + 1, NCQP
            ANORM(J) = DNRM2(N,AQP(J,1),LDAQP)
   80    CONTINUE
C
C        Count the number of linear constraints in the working set and
C        move them to the front of KACTIV.  Compute the norm of the
C        matrix of constraints in the working set.
C        Let K1  point to the first nonlinear constraint.  Constraints
C        with indices KACTIV(K1),..., KACTIV(NACTIV)  must be
C        refactorized.
C
         ASIZE = ZERO
         LINACT = 0
         K1 = NACTIV + 1
         DO 100 K = 1, NACTIV
            I = KACTIV(K)
            ASIZE = MAX(ASIZE,ANORM(I))
C
            IF (I.LE.NCLIN) THEN
               LINACT = LINACT + 1
               IF (LINACT.NE.K) THEN
                  ISWAP = KACTIV(LINACT)
                  KACTIV(LINACT) = I
                  KACTIV(K) = ISWAP
               END IF
            ELSE
C
C              Record the old position of the 1st. nonlinear constraint.
C
               IF (K1.GT.NACTIV) K1 = K
            END IF
  100    CONTINUE
C
         IF (NACTIV.LE.1) CALL F06FLF(NCQP,ANORM,1,ASIZE,AMIN)
C
C        Compute the absolute values of the nonlinear constraints in
C        the working set.  Use DX as workspace.
C
         DO 120 K = LINACT + 1, NACTIV
            J = N + KACTIV(K)
            IF (ISTATE(J).EQ.1) DX(K) = ABS(QPBL(J))
            IF (ISTATE(J).GE.2) DX(K) = ABS(QPBU(J))
  120    CONTINUE
C
C        Sort the elements of KACTIV corresponding to nonlinear
C        constraints in descending order of violation (i.e.,
C        the first element of KACTIV for a nonlinear constraint
C        is associated with the most violated constraint.)
C        In this way, the rows of the Jacobian corresponding
C        to the more violated constraints tend to be included
C        in the  TQ  factorization.
C
C        The sorting procedure is taken from the simple insertion
C        sort in D. Knuth, ACP Volume 3, Sorting and Searching,
C        Page 81.  It should be replaced by a faster sort if the
C        number of active nonlinear constraints becomes large.
C
         DO 160 K = LINACT + 2, NACTIV
            L = K
            VIOL = DX(L)
            KVIOL = KACTIV(L)
C           WHILE (L .GT. LINACT+1  .AND.  DX(L-1) .LT. VIOL) DO
  140       IF (L.GT.LINACT+1) THEN
               IF (DX(L-1).LT.VIOL) THEN
                  DX(L) = DX(L-1)
                  KACTIV(L) = KACTIV(L-1)
                  L = L - 1
                  GO TO 140
               END IF
C              END WHILE
            END IF
            DX(L) = VIOL
            KACTIV(L) = KVIOL
  160    CONTINUE
C
         K2 = NACTIV
         NACTIV = K1 - 1
         NZ = NFREE - NACTIV
C
C        Update the factors  R,  T  and  Q  to include constraints
C        K1  through  K2.
C
         IF (K1.LE.K2) CALL E04NCY(UNITQ,VERTEX,INFORM,K1,K2,NACTIV,
     *                             NARTIF,NZ,NFREE,NRANK,NREJTD,NRPQ,
     *                             NGQ,N,LDZY,LDAQP,LDR,LDT,ISTATE,
     *                             KACTIV,KX,CONDMX,AQP,R,W(LT),W(LRPQ),
     *                             W(LGQ),W(LZY),W(LWRK1),DX,W(LRLAM),
     *                             MSGQP)
      END IF
C
C     ==================================================================
C     Solve for DX, the vector of minimum two-norm that satisfies the
C     constraints in the working set.
C     ==================================================================
      CALL E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,LDAQP,
     *            LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,ADX,QPBL,QPBU,
     *            W(LRPQ),W(LRPQ0),DX,W(LGQ),R,W(LT),W(LZY),W(LWRK1))
C
C     ==================================================================
C     Solve a quadratic program for the search direction  DX  and
C     multiplier estimates  CLAMDA.
C     ==================================================================
C     If there is no feasible point for the subproblem,  the sum of
C     infeasibilities is minimized subject to the linear constraints
C     (1  thru  JINF)  being satisfied.
C
      JINF = N + NCLIN
C
      NTRY = 1
C     +    REPEAT
  180 CALL E04NCZ('QP subproblem',QPNAMD,NAMES,LINOBJ,UNITQ,NQPERR,
     *            MINITS,JINF,NCQP,NCTOTL,NACTIV,NFREE,NRANK,NZ,NZ1,N,
     *            LDAQP,LDR,ISTATE,KACTIV,KX,GDX,SSQ,SSQ1,SUMINF,NUMINF,
     *            DXNORM,QPBL,QPBU,AQP,CLAMDA,ADX,QPTOL,R,DX,W)
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99998) NQPERR
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      NVIOL = 0
      IF (NUMINF.GT.0) THEN
C
C           Count the violated linear constraints.
C
         DO 200 J = 1, NPLIN
            IF (ISTATE(J).LT.0) NVIOL = NVIOL + 1
  200    CONTINUE
C
         IF (NVIOL.GT.0) THEN
            NTRY = NTRY + 1
            UNITQ = .TRUE.
            NACTIV = 0
            NFREE = N
            NZ = N
            CALL F06DBF(NCTOTL,(0),ISTATE,1)
C
            CALL E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,
     *                  LDAQP,LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,
     *                  ADX,QPBL,QPBU,W(LRPQ),W(LRPQ0),DX,W(LGQ),R,W(LT)
     *                  ,W(LZY),W(LWRK1))
         END IF
      END IF
      IF ( .NOT. (NVIOL.EQ.0 .OR. NTRY.GT.2)) GO TO 180
C     +    UNTIL (    NVIOL .EQ. 0  .OR.  NTRY .GT. 2)
C
C     ==================================================================
C     Count the number of nonlinear constraint gradients in the  QP
C     working set.  Make sure that all small  QP  multipliers associated
C     with nonlinear inequality constraints have the correct sign.
C     ==================================================================
      NLNACT = 0
      IF (NACTIV.GT.0 .AND. NCNLN.GT.0) THEN
         DO 220 K = 1, NACTIV
            L = KACTIV(K)
            IF (L.GT.NCLIN) THEN
               NLNACT = NLNACT + 1
               J = N + L
               IF (ISTATE(J).EQ.1) CLAMDA(J) = MAX(ZERO,CLAMDA(J))
               IF (ISTATE(J).EQ.2) CLAMDA(J) = MIN(ZERO,CLAMDA(J))
            END IF
  220    CONTINUE
      END IF
C
      LINACT = NACTIV - NLNACT
C
C     ------------------------------------------------------------------
C     Extract various useful quantities from the QP solution.
C     ------------------------------------------------------------------
C     Compute  HPQ = R'R(pq)  from the transformed gradient of the QP
C     objective function and  R(pq)  from the transformed residual.
C
      CALL DSCAL(N,(-ONE),W(LRPQ),1)
      CALL DAXPY(N,(-ONE),W(LGQ),1,W(LHPQ),1)
      QPCURV = TWO*SSQ
C
      IF (NCNLN.GT.0) THEN
         IF (NUMINF.GT.0) THEN
            FEASQP = .FALSE.
            CALL F06FBF(NCTOTL,(ZERO),CLAMDA,1)
C
            IF (NZ.GT.0) THEN
C              ---------------------------------------------------------
C              Compute a null space element for the search direction
C              as the solution of  Z'HZ(pz) = -Z'g - Z'HY(py).
C              ---------------------------------------------------------
C              Overwrite DX with the transformed search direction
C              Q'(dx).  The first NZ elements of DX are zero.
C
               CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,DX,W(LZY),W(LWRK1)
     *                     )
C
C              Overwrite the first NZ elements of DX with the solution
C              of  (Rz)u = -(v + w),  where  (Rz)'w = Z'g  and  v  is
C              vector of first NZ elements of  R(pq).
C
               CALL DCOPY(NZ,W(LGQ),1,DX,1)
               CALL DTRSV('U','T','N',NZ,R,LDR,DX,1)
C
               CALL DAXPY(NZ,(ONE),W(LRPQ),1,DX,1)
C
               CALL DTRSV('U','N','N',NZ,R,LDR,DX,1)
               CALL DSCAL(NZ,(-ONE),DX,1)
C
C              Recompute RPQ, HPQ, GDX and QPCURV.
C
               CALL DCOPY(NLNX,DX,1,W(LRPQ),1)
               CALL DTRMV('U','N','N',NLNX,R,LDR,W(LRPQ),1)
               IF (NLNX.LT.N) CALL DGEMV('N',NLNX,N-NLNX,ONE,R(1,NLNX+1)
     *                                   ,LDR,DX(NLNX+1),1,ONE,W(LRPQ),
     *                                   1)
C
               GDX = DDOT(N,W(LGQ),1,DX,1)
               QPCURV = DDOT(N,W(LRPQ),1,W(LRPQ),1)
C
               CALL E04NBW(3,N,NZ,NFREE,LDZY,UNITQ,KX,DX,W(LZY),W(LWRK1)
     *                     )
C
C              ---------------------------------------------------------
C              Recompute ADX and the 2-norm of DX.
C              ---------------------------------------------------------
               DXNORM = DNRM2(N,DX,1)
               IF (NCQP.GT.0) CALL DGEMV('N',NCQP,N,ONE,AQP,LDAQP,DX,1,
     *                                   ZERO,ADX,1)
C
               IF (NPDBG .AND. INPDBG(2).GT.0) THEN
                  WRITE (REC,FMT=99997)
                  CALL X04BAY(IPRINT,2,REC)
                  DO 240 I = 1, N, 5
                     WRITE (REC,FMT=99993) (DX(J),J=I,MIN(I+4,N))
                     CALL X04BAF(IPRINT,REC(1))
  240             CONTINUE
               END IF
            END IF
C
            CALL DCOPY(NLNX,W(LRPQ),1,W(LHPQ),1)
            CALL DTRMV('U','T','N',NLNX,R,LDR,W(LHPQ),1)
            IF (NLNX.LT.N) CALL DGEMV('T',NLNX,N-NLNX,ONE,R(1,NLNX+1),
     *                                LDR,W(LRPQ),1,ZERO,W(LHPQ+NLNX),1)
         END IF
C
C        ===============================================================
C        For given values of the objective function and constraints,
C        attempt to minimize the merit function with respect to each
C        slack variable.
C        ===============================================================
         DO 280 I = 1, NCNLN
            J = NPLIN + I
            CON = C(I)
C
            IF ( .NOT. FEASQP .AND. VIOLN(I).NE.ZERO .AND. RHO(I)
     *          .LE.ZERO) RHO(I) = ONE
C
            QUOTNT = F06BLF(CMUL(I),SCALE*RHO(I),OVERFL)
C
C           Define the slack variable to be  CON - MULT / RHO.
C           Force each slack to lie within its upper and lower bounds.
C
            IF (BL(J).GT.BIGLOW) THEN
               IF (QPBL(J).GE.-QUOTNT) THEN
                  SLK(I) = BL(J)
                  GO TO 260
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               IF (QPBU(J).LE.-QUOTNT) THEN
                  SLK(I) = BU(J)
                  GO TO 260
               END IF
            END IF
C
            SLK(I) = CON - QUOTNT
C
C           The slack has been set within its bounds.
C
  260       CS(I) = CON - SLK(I)
C
C           ------------------------------------------------------------
C           Compute the search direction for the slacks and multipliers.
C           ------------------------------------------------------------
            DSLK(I) = ADX(NCLIN+I) + CS(I)
C
            IF (FEASQP) THEN
C
C              If any constraint is such that  (DLAM)*(C - S)  is
C              positive,  the merit function may be reduced immediately
C              by substituting the QP multiplier.
C
               DLAM(I) = CLAMDA(J) - CMUL(I)
               IF (DLAM(I)*CS(I).GE.ZERO) THEN
                  CMUL(I) = CLAMDA(J)
                  DLAM(I) = ZERO
               END IF
            ELSE
C
C              The  QP  subproblem was infeasible.
C
               DLAM(I) = ZERO
C
               IF (ISTATE(J).LT.0 .OR. VIOLN(I).NE.ZERO) DSLK(I) = ZERO
C
            END IF
  280    CONTINUE
C
         IF ( .NOT. FEASQP) RHONRM = DNRM2(NCNLN,RHO,1)
C
         IF (NPDBG .AND. INPDBG(2).GT.0) THEN
            WRITE (REC,FMT=99996)
            CALL X04BAY(IPRINT,2,REC)
            DO 300 J = 1, NCNLN, 5
               WRITE (REC,FMT=99994) (VIOLN(I),I=J,MIN(J+4,NCNLN))
               CALL X04BAF(IPRINT,REC(1))
  300       CONTINUE
            WRITE (REC,FMT=99995)
            CALL X04BAY(IPRINT,2,REC)
            DO 320 J = 1, NCNLN, 5
               WRITE (REC,FMT=99994) (SLK(I),I=J,MIN(J+4,NCNLN))
               CALL X04BAF(IPRINT,REC(1))
  320       CONTINUE
         END IF
      END IF
C
      CALL F06DFF(LDBG,INPDBG,1,ICMDBG,1)
      IDBG = IDBGSV
C
      RETURN
C
C
C     End of  E04UCU. (NPIQP)
C
99999 FORMAT (/1X,79('-'),/' Start of major itn',I6)
99998 FORMAT (/' //E04UCU // NQPERR',/' //E04UCU // ',I6)
99997 FORMAT (/' //E04UCU // DX recomputed with null space portion...')
99996 FORMAT (/' //E04UCU // Violations = ')
99995 FORMAT (/' //E04UCU // Slacks     = ')
99994 FORMAT (1P,5D15.6)
99993 FORMAT (5G12.3)
      END
      SUBROUTINE E04NCJ(PRBTYP,ISDEL,ITER,JADD,JDEL,MSGLVL,NACTIV,NFREE,
     *                  N,NCLIN,NRANK,LDR,LDT,NZ,NRZ,ISTATE,ALFA,CONDRZ,
     *                  CONDT,GFNORM,GZRNRM,NUMINF,SUMINF,CTX,SSQ,AX,R,
     *                  T,X,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1062 (JUL 1993).
C
C     ******************************************************************
C     E04NCJ  prints various levels of output for  E04NCZ.
C
C           Msg        Cumulative result
C           ---        -----------------
C
C        le   0        no output.
C
C        eq   1        nothing now (but full output later).
C
C        eq   5        one terse line of output.
C
C        ge  10        same as 5 (but full output later).
C
C        ge  20        constraint status,  x  and  Ax.
C
C        ge  30        diagonals of  T  and  R.
C
C
C     Debug printing is performed depending on the logical variable
C     LSDBG.  LSDBG  is set true when  IDBG  major iterations have
C     been performed. At this point,  printing is done according to
C     a string of binary digits of the form  SVT  (stored in the
C     integer array  ILSDBG).
C
C     S  set 'on'  gives information from the maximum step routine
C                  E04UCG.
C     V  set 'on'  gives various vectors in  E04NCZ  and its
C                  auxiliaries.
C     T  set 'on'  gives a trace of which routine was called and an
C                  indication of the progress of the run.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCJ dated 21-Oct-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MLINE1, MLINE2
      PARAMETER         (MLINE1=50000,MLINE2=50000)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDRZ, CONDT, CTX, GFNORM, GZRNRM, SSQ,
     *                  SUMINF
      INTEGER           ISDEL, ITER, JADD, JDEL, LDR, LDT, MSGLVL, N,
     *                  NACTIV, NCLIN, NFREE, NRANK, NRZ, NUMINF, NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), R(LDR,*), T(LDT,*), WORK(N), X(N)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  OBJ
      INTEGER           I, ITN, J, K, KADD, KDEL, NART, NDF
      LOGICAL           FIRST, LINOBJ, NEWSET, PRTHDR
      CHARACTER*2       LADD, LDEL
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(0:5)
      CHARACTER*120     REC(4)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Data statements ..
      DATA              LSTATE(0), LSTATE(1), LSTATE(2)/'  ', 'L ',
     *                  'U '/
      DATA              LSTATE(3), LSTATE(4), LSTATE(5)/'E ', 'F ',
     *                  'A '/
C     .. Executable Statements ..
      IF (MSGLVL.GE.15) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
            CALL X04BAY(ISUMM,4,REC)
         END IF
      END IF
C
      IF (MSGLVL.GE.5) THEN
C
         FIRST = ITER .EQ. 0
         LINOBJ = NRANK .EQ. 0
C
         ITN = MOD(ITER,10000)
         NDF = MOD(NRZ,10000)
C
         NART = NZ - NRZ
C
         IF (JDEL.GT.0) THEN
            KDEL = ISDEL
         ELSE IF (JDEL.LT.0) THEN
            JDEL = NART + 1
            KDEL = 5
         ELSE
            KDEL = 0
         END IF
C
         IF (JADD.GT.0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF
C
         LDEL = LSTATE(KDEL)
         LADD = LSTATE(KADD)
C
         IF (NUMINF.GT.0) THEN
            OBJ = SUMINF
         ELSE
            OBJ = SSQ + CTX
         END IF
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. MLINE1
            PRTHDR = MSGLVL .GE. 15 .OR. FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               IF (LINOBJ) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAY(ISUMM,3,REC)
               ELSE
                  WRITE (REC,FMT=99996)
                  CALL X04BAY(ISUMM,3,REC)
               END IF
               LINES1 = 0
            END IF
C
            IF (LINOBJ) THEN
               WRITE (REC,FMT=99994) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *           NUMINF, OBJ, N - NFREE, NACTIV, NART, NDF, GZRNRM,
     *           GFNORM, CONDT
               CALL X04BAF(ISUMM,REC(1))
            ELSE
               WRITE (REC,FMT=99994) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *           NUMINF, OBJ, N - NFREE, NACTIV, NART, NDF, GZRNRM,
     *           GFNORM, CONDT, CONDRZ
               CALL X04BAF(ISUMM,REC(1))
            END IF
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. MLINE2
            PRTHDR = FIRST .OR. NEWSET
C
            IF (PRTHDR) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,3,REC)
               LINES2 = 0
            END IF
C
            WRITE (REC,FMT=99995) ITN, ALFA, NUMINF, OBJ, GZRNRM
            CALL X04BAF(IPRINT,REC(1))
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99993) PRBTYP
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99992)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 I = 1, N, 5
                  WRITE (REC,FMT=99987) (X(J),ISTATE(J),J=I,MIN(I+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99991)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 I = 1, NCLIN, 5
                     WRITE (REC,FMT=99987) (AX(K),ISTATE(N+K),K=I,
     *                 MIN(I+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
C
               IF (MSGLVL.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(NACTIV,NZ+1),LDT-1,WORK,1)
                     WRITE (REC,FMT=99990) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 60 I = 1, NACTIV, 5
                        WRITE (REC,FMT=99986) (WORK(J),J=I,
     *                    MIN(I+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   60                CONTINUE
                  END IF
                  IF (NRANK.GT.0) THEN
                     WRITE (REC,FMT=99989) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 I = 1, NRANK, 5
                        WRITE (REC,FMT=99986) (R(J,J),J=I,MIN(I+4,NRANK)
     *                    )
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
               END IF
               WRITE (REC,FMT=99988)
               CALL X04BAY(ISUMM,3,REC)
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCJ. (LSPRT)
C
99999 FORMAT (//' ',A2,' iteration',I5,/' =================')
99998 FORMAT (//' Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd  ',
     *       'Lin  Art   Zr  Norm Gz  Norm Gf   Cond T')
99997 FORMAT (//' Itn     Step Ninf Sinf/Objective  Norm Gz')
99996 FORMAT (//' Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd  ',
     *       'Lin  Art   Zr  Norm Gz  Norm Gf   Cond T  Cond Rz')
99995 FORMAT (I4,1P,D9.1,I5,D15.6,D9.1)
99994 FORMAT (I4,I5,A1,I5,A1,1P,D9.1,I5,D16.8,4I5,4D9.1)
99993 FORMAT (/' Values and status of the ',A2,' constraints',/' -----',
     *       '----------------------------------')
99992 FORMAT (/' Variables...')
99991 FORMAT (/' General linear constraints...')
99990 FORMAT (/' Diagonals of ',A2,' working set factor T')
99989 FORMAT (/' Diagonals of ',A2,' triangle R         ')
99988 FORMAT (//' ----------------------------------------------------',
     *       '-------------------------------------------')
99987 FORMAT (1X,5(1P,D15.6,I5))
99986 FORMAT (1P,5D15.6)
      END
      SUBROUTINE E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,BU,A,
     *                  FEATOL,CVEC,X,WTINF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFH  finds the number and weighted sum of infeasibilities for
C     the bounds and linear constraints.   An appropriate gradient
C     is returned in cvec.
C
C     Positive values of  istate(j)  will not be altered.  These mean
C     the following...
C
C               1             2           3
C           a'x = bl      a'x = bu     bl = bu
C
C     Other values of  istate(j)  will be reset as follows...
C           a'x lt bl     a'x gt bu     a'x free
C              - 2           - 1           0
C
C     Original version written 31-October-1984.
C     This version of E04MFH dated  1-January-1987.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, SUMINF
      INTEGER           LDA, N, NCLIN, NUMINF
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(*), BU(*), CVEC(N), FEATOL(*),
     *                  WTINF(*), X(N)
      INTEGER           ISTATE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CTX, FEASJ, S, WEIGHT
      INTEGER           J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      BIGUPP = BIGBND
      BIGLOW = -BIGBND
C
      NUMINF = 0
      SUMINF = ZERO
      CALL F06FBF(N,(ZERO),CVEC,1)
C
      DO 40 J = 1, N + NCLIN
         IF (ISTATE(J).LE.0) THEN
            FEASJ = FEATOL(J)
            IF (J.LE.N) THEN
               CTX = X(J)
            ELSE
               K = J - N
               CTX = DDOT(N,A(K,1),LDA,X,1)
            END IF
            ISTATE(J) = 0
C
C           See if the lower bound is violated.
C
            IF (BL(J).GT.BIGLOW) THEN
               S = BL(J) - CTX
               IF (S.GT.FEASJ) THEN
                  ISTATE(J) = -2
                  WEIGHT = -WTINF(J)
                  GO TO 20
               END IF
            END IF
C
C           See if the upper bound is violated.
C
            IF (BU(J).GE.BIGUPP) GO TO 40
            S = CTX - BU(J)
            IF (S.LE.FEASJ) GO TO 40
            ISTATE(J) = -1
            WEIGHT = WTINF(J)
C
C           Add the infeasibility.
C
   20       NUMINF = NUMINF + 1
            SUMINF = SUMINF + ABS(WEIGHT)*S
            IF (J.LE.N) THEN
               CVEC(J) = WEIGHT
            ELSE
               CALL DAXPY(N,WEIGHT,A(K,1),LDA,CVEC,1)
            END IF
         END IF
   40 CONTINUE
      RETURN
C
C     End of  E04MFH.  (CMSINF)
C
      END
      SUBROUTINE E04NFP(UNITQ,IT,N,NACTIV,NFREE,NGQ,NZ,NRZ,LDA,LDQ,LDT,
     *                  JDEL,KDEL,KACTIV,KX,A,T,GQM,Q,WORK,C,S)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04NFP   updates the matrices  Z, Y, T, R  and  D  associated with
C     factorizations
C
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C
C     when a regular, temporary or artificial constraint is deleted
C     from the working set.
C
C     The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C     with its (1,1) element in position  (IT,JT)  of the array  T.
C
C     Original version written by PEG,  31-October-1984.
C     This version of  E04NFP  dated 14-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           IT, JDEL, KDEL, LDA, LDQ, LDT, N, NACTIV, NFREE,
     *                  NGQ, NRZ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), S(N),
     *                  T(LDT,*), WORK(N)
      INTEGER           KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, SN
      INTEGER           I, IR, ITDEL, J, JART, JT, K, L, NPIV, NRZ1,
     *                  NSUP
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSWAP, F06BAF, F06FBF, F06FLF, F06QRF,
     *                  F06QXF, X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      JT = NZ + 1
C
      IF (JDEL.GT.0) THEN
C
C        Regular constraint or temporary bound deleted.
C
         IF (JDEL.LE.N) THEN
C
C           Case 1.  A simple bound has been deleted.
C           =======  Columns  NFREE+1  and  IR  of GQM' must be swapped.
C
            IR = NZ + KDEL
            IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
               WRITE (REC,FMT=99998) NACTIV, NRZ, NZ, NFREE, IR, JDEL,
     *           UNITQ
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
            ITDEL = NACTIV + 1
            NFREE = NFREE + 1
            IF (NFREE.LT.IR) THEN
               KX(IR) = KX(NFREE)
               KX(NFREE) = JDEL
               CALL DSWAP(NGQ,GQM(NFREE,1),N,GQM(IR,1),N)
            END IF
C
            IF ( .NOT. UNITQ) THEN
C
C              Copy the incoming column of  A(free)  into the end of  T.
C
               DO 20 K = 1, NACTIV
                  I = KACTIV(K)
                  T(NACTIV-K+1,NFREE) = A(I,JDEL)
   20          CONTINUE
C
C              Expand  Q  by adding a unit row and column.
C
               IF (NFREE.GT.1) THEN
                  CALL F06FBF(NFREE-1,ZERO,Q(NFREE,1),LDQ)
                  CALL F06FBF(NFREE-1,ZERO,Q(1,NFREE),1)
               END IF
               Q(NFREE,NFREE) = ONE
            END IF
         ELSE
C
C           Case 2.  A general constraint has been deleted.
C           =======
C
            IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
               WRITE (REC,FMT=99997) NACTIV, NRZ, NZ, NFREE, KDEL, JDEL
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
C           Delete row  ITDEL  of  T  and move up the ones below it.
C           T  becomes lower Hessenberg.
C
            ITDEL = KDEL
            DO 60 K = ITDEL, NACTIV
               J = JT + K - 1
               DO 40 L = ITDEL, K - 1
                  I = IT + L - 1
                  T(I,J) = T(I+1,J)
   40          CONTINUE
   60       CONTINUE
C
            DO 80 I = NACTIV - ITDEL + 1, NACTIV - 1
               KACTIV(I) = KACTIV(I+1)
   80       CONTINUE
            NACTIV = NACTIV - 1
         END IF
C
         NZ = NZ + 1
C
         IF (NACTIV.EQ.0) THEN
            DTMAX = ONE
            DTMIN = ONE
         ELSE
C           ------------------------------------------------------------
C           Restore the NACTIV x (NACTIV+1) upper-Hessenberg matrix  T
C           to upper-triangular form.  The  NSUP  super-diagonal
C           elements are removed by a backward sweep of rotations.
C           The rotation for the  (1,1)-th  element of  T  is generated
C           separately.
C           ------------------------------------------------------------
            NSUP = ITDEL - 1
C
            IF (NSUP.GT.0) THEN
               NPIV = JT + ITDEL - 1
               IF (NSUP.GT.1) THEN
                  CALL DCOPY(NSUP-1,T(IT+1,JT+1),LDT+1,S(JT+1),1)
                  CALL F06QRF('Right',NACTIV,1,NSUP,C(JT+1),S(JT+1),
     *                        T(IT,JT+1),LDT)
               END IF
C
               CALL F06BAF(T(IT,JT+1),T(IT,JT),CS,SN)
               T(IT,JT) = ZERO
               S(JT) = -SN
               C(JT) = CS
               CALL F06QXF('Right','Variable','Backwards',NFREE,NFREE,
     *                     NZ,NPIV,C,S,Q,LDQ)
               CALL F06QXF('Left ','Variable','Backwards',NPIV,NGQ,NZ,
     *                     NPIV,C,S,GQM,N)
            END IF
C
            JT = JT + 1
            CALL F06FLF(NACTIV,T(IT,JT),LDT+1,DTMAX,DTMIN)
         END IF
      END IF
C
      NRZ1 = NRZ + 1
C
      IF (NZ.GT.NRZ) THEN
         IF (JDEL.GT.0) THEN
            JART = NRZ1 - 1 + IDAMAX(NZ-NRZ1+1,GQM(NRZ1,1),1)
         ELSE
            JART = -JDEL
         END IF
C
         IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) NZ, NRZ1, JART
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         IF (JART.GT.NRZ1) THEN
C
C           Swap columns  NRZ1  and  JART  of  Q  and  GQM.
C
            IF (UNITQ) THEN
               K = KX(NRZ1)
               KX(NRZ1) = KX(JART)
               KX(JART) = K
            ELSE
               CALL DSWAP(NFREE,Q(1,NRZ1),1,Q(1,JART),1)
            END IF
C
            CALL DSWAP(NGQ,GQM(NRZ1,1),N,GQM(JART,1),N)
         END IF
      END IF
C
      NRZ = NRZ1
C
      RETURN
C
C     End of  E04NFP.  (RZDEL)
C
99999 FORMAT (/' //E04NFP //  Columns NRZ and JART swapped.       ',
     *       /' //E04NFP //      NZ   NRZ   JART                 ',
     *       /' //E04NFP //  ',3I6)
99998 FORMAT (/' //E04NFP //  Simple bound deleted.               ',/
     *       ' //E04NFP //  NACTIV   NRZ    NZ NFREE    IR  JDEL UNITQ',
     *       /' //E04NFP //  ',6I6,L6)
99997 FORMAT (/' //E04NFP //  General constraint deleted.         ',
     *       /' //E04NFP //  NACTIV   NRZ    NZ NFREE  KDEL  JDEL',
     *       /' //E04NFP //  ',6I6)
      END

      SUBROUTINE E04MFS(FIRSTV,N,NCLIN,ISTATE,BIGALF,BIGBND,PNORM,
     *                  HITLOW,MOVE,ONBND,UNBNDD,ALFA,ALFAP,JHIT,ANORM,
     *                  AP,AX,BL,BU,FEATOL,FEATLU,P,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFS  finds a step ALFA such that the point x + ALFA*p reaches
C     one of the linear constraints (including bounds).
C
C     In this version of E04MFS, when x is infeasible, the number of
C     infeasibilities will never increase.  If the number stays the
C     same, the sum of infeasibilities will decrease.  If the number
C     decreases by one or more,  the sum of infeasibilities will usually
C     decrease also, but occasionally it will increase after the step
C     ALFA  is taken.  (Convergence is still assured because the number
C     has decreased.)
C
C     Three possible steps are computed as follows:
C
C     alfaf = the maximum step that can be taken without violating
C              one of the constraints that are currently satisfied.
C
C     ALFAI = reaches a linear constraint that is currently violated.
C              Usually this will be the furthest such constraint along
C              p, subject to the angle between the constraint normal and
C              p being reasonably close to the maximum value among
C              infeasible constraints,  but if FIRSTV = .true. it will
C              be the first one along p.  The latter case applies only
C              when the problem has been determined to be infeasible,
C              and the sum of infeasibilities are being minimized.
C              (ALFAI is not defined when x is feasible.)
C
C     ALFAI is needed occasionally when infeasible, to prevent
C     going unnecessarily far when alfaf is quite large.  It will
C     always come into effect when x is about to become feasible.
C     (The sum of infeasibilities will decrease initially as ALFA
C     increases from zero, but may start increasing for larger steps.
C     Choosing a large ALFAI allows several elements of  x  to
C     become feasible at the same time.
C
C     In the end, we take  ALFA = alfaf  if x is feasible, or if
C     ALFAI > ALFAP (where  ALFAP  is the perturbed step from pass 1).
C     Otherwise,  we take  ALFA = ALFAI.
C
C     Input parameters
C     ----------------
C     BIGALF defines what should be treated as an unbounded step.
C     BIGBND provides insurance for detecting unboundedness.
C            If ALFA reaches a bound as large as BIGBND, it is
C            classed as an unbounded step.
C     FEATOL is the array of current feasibility tolerances used by
C            E04MFH.  Typically in the range 0.5*TOLX to 0.99*TOLX,
C            where TOLX is the FEATOL specified by the user.
C     TOLINC (in common) is used to determine STEPMN (see below),
C            the minimum positive step.
C     ISTATE is set as follows:
C            ISTATE(j) = -2  if a'x .lt. BL - FEATOL
C                      = -1  if a'x .gt. BU + FEATOL
C                      =  0  if a'x is not in the working set
C                      =  1  if a'x is in the working set at BL
C                      =  2  if a'x is in the working set at BU
C                      =  3  if a'x is in the working set (an equality)
C                      =  4  if x(j) is temporarily fixed.
C            values -2 and -1 do not occur once feasible.
C     BL     the lower bounds on the variables.
C     BU     the upper bounds on ditto.
C     x      the values of       ditto.
C     p      the search direction.
C
C
C     Output Parameters
C     -----------------
C     HITLOW  = true  if a lower bound restricted ALFA.
C             = false otherwise.
C     MOVE    = true  if  EXACT ge STEPMN  (defined at end of code).
C     ONBND   = true  if  ALFA = EXACT.  This means that the step  ALFA
C                     moves x  exactly onto one of its constraints,
C                     namely  bound.
C             = false if the exact step would be too small
C                     ( EXACT .lt. STEPMN ).
C               (with these definitions,  MOVE = ONBND).
C     UNBNDD  = true  if ALFA = BIGALF.  JHIT may possibly be zero.
C               The parameters HITLOW, MOVE, ONBND, BOUND and EXACT
C               should not be used.
C     JHIT    = the index (if any) such that constraint JHIT reaches
C               a bound.
C     BOUND   = the bound value BL(JHIT) or BU(JHIT) corresponding
C               to HITLOW.
C     EXACT   = the step that would take constraint JHIT exactly onto
C               BOUND.
C     ALFA    = an allowable, positive step.
C               if UNBNDD is true,  ALFA = STEPMX.
C               otherwise,          ALFA = max( STEPMN, EXACT ).
C
C
C     E04MFS is based on MINOS 5.2 routine M5CHZR, which implements the
C     expand procedure to deal with degeneracy. The step alfaf is
C     chosen as in the two-pass approach of Paula Harris (1973), except
C     that this version insists on returning a positive step, ALFA.
C     Two features make this possible:
C
C        1. FEATOL increases slightly each iteration.
C
C        2. The blocking constraint, when added to the working set,
C           retains the value Ax(JHIT) + ALFA * Ap(JHIT),
C           even if this is not exactly on the blocking bound.
C
C     For infeasible variables moving towards their bound, we require
C     the rate of change of the chosen constraint to be at least GAMMA
C     times as large as the biggest available.  This still gives us
C     freedom in pass 2.
C     GAMMA = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
C     GAMMA = 0.001 seems to be safe.
C
C
C     Original version written by PEG,  19-April 1988.
C     This version of  E04MFS  dated   6-Jul-1988.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  GAMMA
      PARAMETER         (GAMMA=1.0D-3)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFAP, BIGALF, BIGBND, PNORM
      INTEGER           JHIT, N, NCLIN
      LOGICAL           FIRSTV, HITLOW, MOVE, ONBND, UNBNDD
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(N+NCLIN),
     *                  BU(N+NCLIN), FEATLU(N+NCLIN), FEATOL(N+NCLIN),
     *                  P(N), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9, TOLINC, TOLX0
      INTEGER           IPRINT, ISUMM, ITNFIX, KDEGEN, LINES1, LINES2,
     *                  NDEGEN, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG), NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFAI, ATP, ATPABS, ATPMXF, ATPMXI, ATPSCD, ATX,
     *                  BIGLOW, BIGUPP, BOUND, DELTA, EXACT, RES,
     *                  STEPMN, TOLPIV
      INTEGER           I, J, JHITF, JHITI, JS
      LOGICAL           BLOCKF, BLOCKI
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04MF/TOLX0, TOLINC, KDEGEN, NDEGEN, ITNFIX,
     *                  NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
C     TOLPIV is a tolerance to exclude negligible elements of a'p.
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
      TOLPIV = EPSPT9*PNORM
C
C     ------------------------------------------------------------------
C     First pass -- find steps to perturbed constraints, so that
C     ALFAP will be slightly larger than the true step.
C     In degenerate cases, this strategy gives us some freedom in the
C     second pass.  The general idea follows that described by P.M.J.
C     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
C     ------------------------------------------------------------------
      IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      ATPMXI = ZERO
      ALFAP = BIGALF
C
      DO 20 J = 1, N + NCLIN
         JS = ISTATE(J)
C
         IF (JS.LE.0) THEN
            DELTA = FEATOL(J)
C
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS/(ONE+ANORM(I))
            END IF
C
            IF (ATPSCD.LE.TOLPIV) THEN
C              ---------------------------------------------------------
C              This constraint appears to be constant along p.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C              ---------------------------------------------------------
               RES = -ONE
C
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing and the lower bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C              if the upper bound is violated. test for bigger ATP.
C
               IF (BL(J).GT.BIGLOW) THEN
                  RES = ATX - BL(J) + DELTA
C
                  IF (RES.LT.ALFAP*ATPABS) ALFAP = RES/ATPABS
               END IF
C
               IF (JS.EQ.-1) ATPMXI = MAX(ATPMXI,ATPSCD)
C
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C              if the lower bound is violated. test for bigger ATP.
C
               IF (BU(J).LT.BIGUPP) THEN
                  RES = BU(J) - ATX + DELTA
C
                  IF (RES.LT.ALFAP*ATP) ALFAP = RES/ATP
               END IF
C
               IF (JS.EQ.-2) ATPMXI = MAX(ATPMXI,ATPSCD)
            END IF
C
            IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
               WRITE (REC,FMT=99998) J, JS, DELTA, RES, ATP, ALFAP,
     *           ATPMXI
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Second pass.
C     For feasible variables, recompute steps without perturbation.
C     amongst constraints that are closer than ALFAP, choose the one
C     That makes the largest angle with the search direction.
C     For infeasible variables, find the largest step subject to a'p
C     being no smaller than GAMMA * max(a'p).
C     ------------------------------------------------------------------
      IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
         WRITE (REC,FMT=99997)
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (FIRSTV) THEN
         ALFAI = BIGALF
      ELSE
         ALFAI = ZERO
      END IF
C
      ATPMXF = ZERO
      ATPMXI = GAMMA*ATPMXI
      JHITF = 0
      JHITI = 0
C
      DO 40 J = 1, N + NCLIN
         JS = ISTATE(J)
C
         IF (JS.LE.0) THEN
C
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS/(ONE+ANORM(I))
            END IF
C
            IF (ATPSCD.LE.TOLPIV) THEN
C              ---------------------------------------------------------
C              this constraint appears to be constant along p.  it is
C              not used to compute the step.  give the residual a value
C              that can be spotted in the debug output.
C              ---------------------------------------------------------
               RES = -ONE
C
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing.
C              ---------------------------------------------------------
C              test for bigger a'p if the lower bound is satisfied.
C              test for smaller alfaf.
C
               IF (ATPSCD.GT.ATPMXF) THEN
C
                  IF (BL(J).GT.BIGLOW) THEN
                     RES = ATX - BL(J)
C
                     IF (RES.LE.ALFAP*ATPABS) THEN
                        ATPMXF = ATPSCD
                        JHITF = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 the upper bound is violated.
C                 test for bigger or smaller ALFAI,  depending on the
C                 value of FIRSTV.
C
                  IF (FIRSTV) THEN
                     RES = ATX - BU(J)
C
                     IF (RES.LE.ALFAI*ATPABS) THEN
                        ALFAI = RES/ATPABS
                        JHITI = J
                     END IF
C
                  ELSE IF (ATPSCD.GE.ATPMXI) THEN
                     RES = ATX - BU(J)
C
                     IF (RES.GT.ALFAI*ATPABS) THEN
                        ALFAI = RES/ATPABS
                        JHITI = J
                     END IF
                  END IF
               END IF
C
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C
               IF (ATPSCD.GT.ATPMXF) THEN
C
                  IF (BU(J).LT.BIGUPP) THEN
                     RES = BU(J) - ATX
C
                     IF (RES.LE.ALFAP*ATP) THEN
                        ATPMXF = ATPSCD
                        JHITF = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 the lower bound is violated.
C                 test for bigger or smaller ALFAI,  depending on the
C                 value of FIRSTV.
C
                  IF (FIRSTV) THEN
                     RES = BL(J) - ATX
C
                     IF (RES.LE.ALFAI*ATP) THEN
                        ALFAI = RES/ATP
                        JHITI = J
                     END IF
                  ELSE IF (ATPSCD.GE.ATPMXI) THEN
                     RES = BL(J) - ATX
C
                     IF (RES.GT.ALFAI*ATP) THEN
                        ALFAI = RES/ATP
                        JHITI = J
                     END IF
                  END IF
               END IF
            END IF
C
            IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
               WRITE (REC,FMT=99996) J, JS, DELTA, RES, ATP, JHITF,
     *           ATPMXF, JHITI, ALFAI
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
   40 CONTINUE
C
C     ------------------------------------------------------------------
C     See if a feasible and/or infeasible constraint blocks.
C     ------------------------------------------------------------------
      BLOCKF = JHITF .GT. 0
      BLOCKI = JHITI .GT. 0
      UNBNDD = .NOT. (BLOCKF .OR. BLOCKI)
C
      IF (UNBNDD) GO TO 60
C
      IF (BLOCKF) THEN
C        ---------------------------------------------------------------
C        A constraint is hit which is currently feasible.
C        The corresponding step alfaf is not used, so no need to get it,
C        but we know that alfaf .le. ALFAP, the step from pass 1.
C        ---------------------------------------------------------------
         JHIT = JHITF
         IF (JHIT.LE.N) THEN
            ATP = P(JHIT)
         ELSE
            ATP = AP(JHIT-N)
         END IF
         HITLOW = ATP .LT. ZERO
      END IF
C
C     If there is a choice between alfaf and ALFAI, it is probably best
C     to take ALFAI.  However, we can't if ALFAI is bigger than ALFAP.
C
      IF (BLOCKI .AND. ALFAI.LE.ALFAP) THEN
C        ---------------------------------------------------------------
C        An infeasible variable reaches its violated bound.
C        ---------------------------------------------------------------
         JHIT = JHITI
         IF (JHIT.LE.N) THEN
            ATP = P(JHIT)
         ELSE
            ATP = AP(JHIT-N)
         END IF
         HITLOW = ATP .GT. ZERO
      END IF
C
      IF (JHIT.LE.N) THEN
         ATX = X(JHIT)
      ELSE
         ATX = AX(JHIT-N)
      END IF
C
C     ------------------------------------------------------------------
C     Try to step exactly onto bound, but make sure the exact step
C     is sufficiently positive.  (Exact will be alfaf or ALFAI.)
C     Since FEATOL increases by  TOLINC  each iteration, we know that
C     a step as large as  STEPMN  (below) will not cause any feasible
C     variables to become infeasible (where feasibility is measured
C     by the current FEATOL).
C     ------------------------------------------------------------------
      IF (HITLOW) THEN
         BOUND = BL(JHIT)
      ELSE
         BOUND = BU(JHIT)
      END IF
C
      UNBNDD = ABS(BOUND) .GE. BIGBND
      IF (UNBNDD) GO TO 60
C
      STEPMN = TOLINC*FEATLU(JHIT)/ABS(ATP)
      EXACT = (BOUND-ATX)/ATP
      ALFA = MAX(STEPMN,EXACT)
      ONBND = ALFA .EQ. EXACT
      MOVE = EXACT .GE. STEPMN
      IF ( .NOT. MOVE) NDEGEN = NDEGEN + 1
C
      RETURN
C     ------------------------------------------------------------------
C     Unbounded.
C     ------------------------------------------------------------------
   60 ALFA = BIGALF
      MOVE = .TRUE.
      ONBND = .FALSE.
C
      IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
         WRITE (REC,FMT=99995) JHIT, ALFA
         CALL X04BAY(IPRINT,4,REC)
      END IF
C
      RETURN
C
C     End of  E04MFS.  (CMCHZR)
C
99999 FORMAT (/' E04MFS entered.  Pass 1.',/'    J  JS         FEATOL ',
     *       '       RES             AP            ALFAP           ATP',
     *       'MXI',/)
99998 FORMAT (I5,I4,3G15.5,2G17.7)
99997 FORMAT (/'                  Pass 2.',/'    J  JS         FEATOL ',
     *       '       RES             AP     JHITF           ATPMXF JHI',
     *       'TI            ALFAI',/)
99996 FORMAT (I5,I4,3G15.5,2(I6,G17.7))
99995 FORMAT (/' //E04MFS//  Unbounded step.',/' //E04MFS//  JHIT     ',
     *       '      ALFA',/' //E04MFS//  ',I4,G15.4)
      END
      SUBROUTINE E04NBX(MSGLVL,NFREE,NROWA,N,NCLIN,NCTOTL,BIGBND,NAMED,
     *                  NAMES,NACTIV,ISTATE,KACTIV,KX,A,BL,BU,C,CLAMDA,
     *                  RLAMDA,X)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-591 (MAR 1988).
C     MARK 16 REVISED. IER-1058 (JUL 1993).
C
C     ******************************************************************
C     E04NBX   creates the expanded Lagrange multiplier vector CLAMDA.
C     If MSGLVL .EQ 1 or MSGLVL .GE. 10,  E04NBX prints  x,  A*x,
C     c(x),  their bounds, the multipliers, and the residuals (distance
C     to the nearer bound).
C
C     E04NBX is called by E04NCZ, E04UCZ and E04UPZ just before exiting.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written  October 1984.
C     This version of  E04NBX  dated  30-Mar-1993.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           MSGLVL, N, NACTIV, NCLIN, NCTOTL, NFREE, NROWA
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NROWA,*), BL(NCTOTL), BU(NCTOTL), C(*),
     *                  CLAMDA(NCTOTL), RLAMDA(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, V, WLAM
      INTEGER           IP, IS, J, K, NFIXED, NPLIN, NZ
      CHARACTER*2       LS
      CHARACTER*5       ID3
      CHARACTER*8       ID4
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(7)
      CHARACTER*5       ID(3)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Data statements ..
      DATA              ID(1)/'Varbl'/
      DATA              ID(2)/'L Con'/
      DATA              ID(3)/'N Con'/
      DATA              LSTATE(1)/'--'/, LSTATE(2)/'++'/
      DATA              LSTATE(3)/'FR'/, LSTATE(4)/'LL'/
      DATA              LSTATE(5)/'UL'/, LSTATE(6)/'EQ'/
      DATA              LSTATE(7)/'TF'/
C     .. Executable Statements ..
C
C
      NPLIN = N + NCLIN
      NZ = NFREE - NACTIV
C
C     Expand multipliers for bounds, linear and nonlinear constraints
C     into the  CLAMDA  array.
C
      CALL F06FBF(NCTOTL,ZERO,CLAMDA,1)
      NFIXED = N - NFREE
      DO 20 K = 1, NACTIV + NFIXED
         IF (K.LE.NACTIV) J = KACTIV(K) + N
         IF (K.GT.NACTIV) J = KX(NZ+K)
         CLAMDA(J) = RLAMDA(K)
   20 CONTINUE
C
      IF (MSGLVL.LT.10 .AND. MSGLVL.NE.1) RETURN
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(IPRINT,4,REC)
      ID3 = ID(1)
C
      DO 40 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
         WLAM = CLAMDA(J)
         IS = ISTATE(J)
         LS = LSTATE(IS+3)
         IF (J.LE.N) THEN
C
C           Section 1 -- the variables  x.
C           ------------------------------
            K = J
            V = X(J)
C
         ELSE IF (J.LE.NPLIN) THEN
C
C           Section 2 -- the linear constraints  A*x.
C           -----------------------------------------
            IF (J.EQ.N+1) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(2)
            END IF
C
            K = J - N
            V = DDOT(N,A(K,1),NROWA,X,1)
         ELSE
C
C           Section 3 -- the nonlinear constraints  c(x).
C           ---------------------------------------------
C
            IF (J.EQ.NPLIN+1) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(3)
            END IF
C
            K = J - NPLIN
            V = C(K)
         END IF
C
C        Print a line for the j-th variable or constraint.
C        -------------------------------------------------
         RES = V - B1
         RES2 = B2 - V
         IF (ABS(RES).GT.ABS(RES2)) RES = RES2
         IP = 1
         IF (B1.LE.(-BIGBND)) IP = 2
         IF (B2.GE.BIGBND) IP = IP + 2
         IF (NAMED) THEN
C
            ID4 = NAMES(J)
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99996) ID4, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99995) ID4, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99994) ID4, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99993) ID4, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
C
         ELSE
C
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99992) ID3, K, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99991) ID3, K, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99990) ID3, K, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99989) ID3, K, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
         END IF
   40 CONTINUE
      RETURN
C
C
C     End of  E04NBX. (CMPRT)
C
99999 FORMAT (//1X,'Varbl',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99998 FORMAT (//1X,'L Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99997 FORMAT (//1X,'N Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99996 FORMAT (1X,A4,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99995 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,2G12.4)
99994 FORMAT (1X,A4,4X,A2,1X,1P,2G14.6,5X,'None',5X,1P,2G12.4)
99993 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,2G12.4)
99992 FORMAT (1X,A1,I3,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99991 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,
     *       2G12.4)
99990 FORMAT (1X,A1,I3,4X,A2,1X,1P,2G14.6,5X,'None',5X,2G12.4)
99989 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,
     *       2G12.4)
      END
      SUBROUTINE E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,NACTIV,NZ,
     *                  NFREE,NRZ,NGQ,N,LDA,LDQ,LDR,LDT,KX,CONDMX,DRZZ,
     *                  A,R,T,GQM,Q,W,C,S,MSGLVL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04NFR  updates the matrices  Z, Y, T, R  and  D  associated with
C     factorizations
C
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C                      R' *D * R  =   Hz
C
C     a) The matrices  R  and  T  are upper triangular.
C     b) The arrays  T  and  R  may be the same array.
C     c) The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C        with its (1,1) element in position  (IT,JT) of the
C        array  T.   The integer  JT  is always  NZ+1.  During regular
C        changes to the working set,  IT = 1;  when several constraints
C        are added simultaneously,  IT  points to the first row of the
C        existing  T.
C     d) The matrix  R  is stored in the first  NZ x NZ  rows
C        and columns of the  NFREE x NFREE  leading principal minor of
C        the array  R.
C     e) If  RSET  is  false,   R  is not touched.
C
C     There are three separate cases to consider (although each case
C     shares code with another)...
C
C     (1) A free variable becomes fixed on one of its bounds when there
C         are already some general constraints in the working set.
C
C     (2) A free variable becomes fixed on one of its bounds when there
C         are only bound constraints in the working set.
C
C     (3) A general constraint (corresponding to row  IADD  of  A) is
C         added to the working set.
C
C     In cases (1) and (2), we assume that  KX(IFIX) = JADD.
C     In all cases,  JADD  is the index of the constraint being added.
C
C     If there are no general constraints in the working set,  the
C     matrix  Q = (Z Y)  is the identity and will not be touched.
C
C     If  NGQ .gt. 0,  the column transformations are applied to the
C     columns of the  (NGQ x N)  matrix  GQM'.
C
C     Original version written by PEG,  31-October-1984.
C     This version of  E04NFR  dated  7-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX, DRZZ
      INTEGER           IADD, IFIX, INFORM, IT, JADD, LDA, LDQ, LDR,
     *                  LDT, MSGLVL, N, NACTIV, NFREE, NGQ, NRZ, NZ
      LOGICAL           RSET, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), R(LDR,*),
     *                  S(N), T(LDT,*), W(N)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8,
     *                  EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  COND, CONDBD, DTNEW, TDTMAX, TDTMIN
      INTEGER           I, J, JT, K, NANEW, NPIV, NSUP
      LOGICAL           BOUND, OVERFL
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04NBW, E04NFM, F06FLF, F06FQF,
     *                  F06QHF, F06QKF, F06QRF, F06QVF, F06QXF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
C     If the condition estimator of the updated T is greater than
C     CONDBD,  a warning message is printed.
C
      CONDBD = ONE/EPSPT9
C
      OVERFL = .FALSE.
      BOUND = JADD .LE. N
      JT = NZ + 1
C
      IF (BOUND) THEN
C        ===============================================================
C        A simple bound has entered the working set.  IADD is not used.
C        ===============================================================
         IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) NACTIV, NRZ, NZ, NFREE, IFIX, JADD,
     *        UNITQ
            CALL X04BAY(IPRINT,4,REC)
         END IF
         NANEW = NACTIV
C
         IF (UNITQ) THEN
C
C           Q is not stored, but  KX  defines an ordering of the columns
C           of the identity matrix that implicitly define Q.
C           Define the sequence of pairwise interchanges P that moves
C           the newly-fixed variable to position  NFREE.
C           Reorder  KX  accordingly.
C
            DO 20 I = 1, NFREE - 1
               IF (I.GE.IFIX) THEN
                  W(I) = I + 1
                  KX(I) = KX(I+1)
               ELSE
                  W(I) = I
               END IF
   20       CONTINUE
         ELSE
C           ------------------------------------------------------------
C           Q  is stored explicitly.
C           ------------------------------------------------------------
C           Set  W = the  (IFIX)-th  row of  Q.
C           Move the  (NFREE)-th  row of  Q  to position IFIX.
C
            CALL DCOPY(NFREE,Q(IFIX,1),LDQ,W,1)
            IF (IFIX.LT.NFREE) THEN
               CALL DCOPY(NFREE,Q(NFREE,1),LDQ,Q(IFIX,1),LDQ)
               KX(IFIX) = KX(NFREE)
            END IF
         END IF
         KX(NFREE) = JADD
      ELSE
C        ===============================================================
C        A general constraint has entered the working set.
C        IFIX is not used.
C        ===============================================================
         IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
            WRITE (REC,FMT=99998) NACTIV, NRZ, NZ, NFREE, IADD, JADD,
     *        UNITQ
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         NANEW = NACTIV + 1
C
C        Transform the incoming row of A by Q'.
C
         CALL DCOPY(N,A(IADD,1),LDA,W,1)
         CALL E04NBW(8,N,NZ,NFREE,LDQ,UNITQ,KX,W,Q,C)
C
C        Check that the incoming row is not dependent upon those
C        already in the working set.
C
         DTNEW = DNRM2(NZ,W,1)
         IF (NACTIV.EQ.0) THEN
C
C           This is the only general constraint in the working set.
C
            COND = F06BLF(ASIZE,DTNEW,OVERFL)
            TDTMAX = DTNEW
            TDTMIN = DTNEW
         ELSE
C
C           There are already some general constraints in the working
C           set.  Update the estimate of the condition number.
C
            TDTMAX = MAX(DTNEW,DTMAX)
            TDTMIN = MIN(DTNEW,DTMIN)
            COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
         END IF
C
         IF (COND.GT.CONDMX .OR. OVERFL) GO TO 80
C
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,Q,LDQ)
            UNITQ = .FALSE.
            IT = 0
         END IF
      END IF
C
      IF (BOUND) THEN
         NPIV = NFREE
      ELSE
         NPIV = NZ
      END IF
C
      IF (UNITQ) THEN
C        ---------------------------------------------------------------
C        The orthogonal matrix  Q  (i.e.,  Q) is not stored explicitly.
C        Apply  P, the sequence of pairwise interchanges that moves the
C        newly-fixed variable to position  NFREE.
C        ---------------------------------------------------------------
         IF (NGQ.GT.0) CALL F06QKF('Left','Transpose',NFREE-1,W,NGQ,GQM,
     *                             N)
C
         IF (RSET) THEN
C
C           Apply the pairwise interchanges to  Rz.
C           The subdiagonal elements generated by this process are
C           stored in  S(IFIX), S(2), ..., S(NRZ-1).
C
            NSUP = NRZ - IFIX
            CALL E04NFM('Right',NRZ,IFIX,NRZ,S,R,LDR)
         END IF
      ELSE
C        ---------------------------------------------------------------
C        The matrix  Q  is stored explicitly.
C        Define a sweep of plane rotations P such that
C                           PW = beta*e(NPIV).
C        The rotations are applied in the planes (1, 2), (2, 3), ...,
C        (NPIV-1, NPIV).  The rotations must be applied to Q, GQM', R
C        and T.
C        ---------------------------------------------------------------
         CALL F06FQF('Varble','Forwrds',NPIV-1,W(NPIV),W,1,C,S)
C
         IF (NGQ.GT.0) CALL F06QXF('Left ','Variable','Forwards',NPIV,
     *                             NGQ,1,NPIV,C,S,GQM,N)
         CALL F06QXF('Right','Variable','Forwards',NFREE,NFREE,1,NPIV,C,
     *               S,Q,LDQ)
C
         IF (RSET) THEN
C
C           Apply the rotations to the triangular part of R.
C           The subdiagonal elements generated by this process are
C           stored in  S(1),  S(2), ..., S(NRZ-1).
C
            NSUP = NRZ - 1
            CALL F06QVF('Right',NRZ,1,NRZ,C,S,R,LDR)
         END IF
      END IF
C
      IF (RSET) THEN
C        ---------------------------------------------------------------
C        Eliminate the  NSUP  subdiagonal elements of  R  stored in
C        S(NRZ-NSUP), ..., S(NRZ-1)  with a left-hand sweep of rotations
C        in planes (NRZ-NSUP, NRZ-NSUP+1), ..., (NRZ-1, NRZ).
C        ---------------------------------------------------------------
         CALL F06QRF('Left ',NRZ,NRZ-NSUP,NRZ,C,S,R,LDR)
C
         IF (NSUP.GT.0 .AND. DRZZ.NE.ONE) THEN
            DRZZ = C(NRZ-1)**2 + DRZZ*S(NRZ-1)**2
         END IF
      END IF
C
      IF ( .NOT. UNITQ) THEN
         IF (BOUND) THEN
C           ------------------------------------------------------------
C           Bound constraint added.   The rotations affect columns
C           NZ+1  thru  NFREE  of  GQM'  and  T.
C           ------------------------------------------------------------
C           The last row and column of  Q  has been transformed to plus
C           or minus the unit vector  e(NFREE).  We can reconstitute the
C           column of GQM' corresponding to the new fixed variable.
C
            IF (W(NFREE).LT.ZERO) THEN
               IF (NGQ.GT.0) CALL DSCAL(NGQ,-ONE,GQM(NFREE,1),N)
            END IF
C
            IF (NACTIV.GT.0) THEN
               T(IT,JT-1) = S(JT-1)*T(IT,JT)
               T(IT,JT) = C(JT-1)*T(IT,JT)
C
               IF (NACTIV.GT.1) THEN
                  CALL F06QVF('Right',NACTIV,1,NACTIV,C(JT),S(JT),
     *                        T(IT,JT),LDT)
                  CALL DCOPY(NACTIV-1,S(JT),1,T(IT+1,JT),LDT+1)
               END IF
C
               JT = JT - 1
               CALL F06FLF(NACTIV,T(IT,JT),LDT+1,TDTMAX,TDTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
         ELSE
C           ------------------------------------------------------------
C           General constraint added.  Install  W  at the front of  T.
C           If there is no room,  shift all the rows down one position.
C           ------------------------------------------------------------
            IT = IT - 1
            IF (IT.LE.0) THEN
               IT = 1
               DO 60 K = 1, NACTIV
                  J = JT + K - 1
                  DO 40 I = K, 1, -1
                     T(I+1,J) = T(I,J)
   40             CONTINUE
   60          CONTINUE
            END IF
            JT = JT - 1
            CALL DCOPY(NANEW,W(JT),1,T(IT,JT),LDT)
         END IF
      END IF
C
C     ==================================================================
C     Prepare to exit.  Check the magnitude of the condition estimator.
C     ==================================================================
   80 IF (NANEW.GT.0) THEN
         IF (COND.LT.CONDMX .AND. .NOT. OVERFL) THEN
C
C           The factorization has been successfully updated.
C
            INFORM = 0
            DTMAX = TDTMAX
            DTMIN = TDTMIN
            IF (COND.GE.CONDBD) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99997) JADD
                  CALL X04BAY(IPRINT,5,REC)
               END IF
            END IF
         ELSE
C
C           The proposed working set appears to be linearly dependent.
C
            INFORM = 1
            IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
               WRITE (REC,FMT=99996)
               CALL X04BAY(IPRINT,2,REC)
               IF (BOUND) THEN
                  WRITE (REC,FMT=99995) ASIZE, DTMAX, DTMIN
                  CALL X04BAY(IPRINT,3,REC)
               ELSE
                  IF (NACTIV.GT.0) THEN
                     WRITE (REC,FMT=99994) ASIZE, DTMAX, DTMIN, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  ELSE
                     WRITE (REC,FMT=99993) ASIZE, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of  E04NFR.  (RZADD)
C
99999 FORMAT (/' //E04NFR //  Simple bound added.',/' //E04NFR //  NAC',
     *       'TIV   NRZ    NZ NFREE  IFIX  JADD UNITQ         ',/' //E',
     *       '04NFR //  ',6I6,L6)
99998 FORMAT (/' //E04NFR //  General constraint added.           ',/
     *       ' //E04NFR //  NACTIV   NRZ    NZ NFREE  IADD  JADD UNITQ',
     *       /' //E04NFR //  ',6I6,L6)
99997 FORMAT (/' XXX  Serious ill-conditioning in the working set afte',
     *       'r adding constraint ',I5,/' XXX  Overflow may occur in s',
     *       'ubsequent iterations.',//)
99996 FORMAT (/' //E04NFR //  Dependent constraint rejected.')
99995 FORMAT (/' //E04NFR //     ASIZE     DTMAX     DTMIN        ',
     *       /' //E04NFR //',1P,3D10.2)
99994 FORMAT (/' //E04NFR //     ASIZE     DTMAX     DTMIN     DTNEW',
     *       /' //E04NFR //',1P,4D10.2)
99993 FORMAT (/' //E04NFR //     ASIZE     DTNEW',/' //E04NFR //',1P,
     *       2D10.2)
      END
      SUBROUTINE E04UCN(FEASQP,N,NCLIN,NCNLN,OBJALF,GRDALF,QPCURV,
     *                  ISTATE,CJDX,CMUL,CS,DLAM,RHO,VIOLN,WORK1,WORK2)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1083 (JUL 1993).
C
C     ******************************************************************
C     E04UCN   computes the value and directional derivative of the
C     augmented Lagrangian merit function.  The penalty parameters
C     RHO(j) are boosted if the directional derivative of the resulting
C     augmented Lagrangian function is not sufficiently negative.  If
C     RHO needs to be increased,  the perturbation with minimum two-norm
C     is found that gives a directional derivative equal to  - p'Hp.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written  27-May-1985.
C     This version of  E04UCN  dated 14-November-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GRDALF, OBJALF, QPCURV
      INTEGER           N, NCLIN, NCNLN
      LOGICAL           FEASQP
C     .. Array Arguments ..
      DOUBLE PRECISION  CJDX(*), CMUL(*), CS(*), DLAM(*), RHO(*),
     *                  VIOLN(*), WORK1(*), WORK2(*)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RHODMP, RHOMAX, RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           INCRUN, NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  PTERM, PTERM2, QNORM, RHO1, RHOI, RHOMIN,
     *                  RHONEW, RTMIN, TSCL
      INTEGER           I, L, NPLIN
      LOGICAL           BOOST, OVERFL
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FCF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      IF (NCNLN.EQ.0) RETURN
C
      RTMIN = WMACH(6)
C
      OBJALF = OBJALF - DDOT(NCNLN,CMUL,1,CS,1)
      GRDALF = GRDALF - DDOT(NCNLN,DLAM,1,CS,1)
C
      CALL DCOPY(NCNLN,CS,1,WORK1,1)
C
      IF ( .NOT. FEASQP) THEN
         NPLIN = N + NCLIN
C
         DO 20 I = 1, NCNLN
            IF (ISTATE(NPLIN+I).LT.0 .OR. VIOLN(I).NE.ZERO) WORK1(I)
     *          = -CJDX(I)
   20    CONTINUE
      END IF
C
      GRDALF = GRDALF + DDOT(NCNLN,WORK1,1,CMUL,1)
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) QPCURV, GRDALF
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (FEASQP) THEN
C
C        Find the quantities that define  rhomin, the vector of minimum
C        two-norm such that the directional derivative is one half of
C        approximate curvature   - (dx)'H(dx).
C
         DO 40 I = 1, NCNLN
            IF (ABS(CS(I)).LE.RTMIN) THEN
               WORK2(I) = ZERO
            ELSE
               WORK2(I) = CS(I)**2
            END IF
   40    CONTINUE
C
         QNORM = DNRM2(NCNLN,WORK2,1)
         TSCL = F06BLF(GRDALF+HALF*QPCURV,QNORM,OVERFL)
         IF (ABS(TSCL).LE.RHOMAX .AND. .NOT. OVERFL) THEN
C           ------------------------------------------------------------
C           Bounded  RHOMIN  found.  The final value of  RHO(J)  will
C           never be less than  RHOMIN(j).  If the  QP  was feasible,  a
C           trial value  RHONEW  is computed that is equal to the
C           geometric mean of the previous  RHO  and a damped value of
C           RHOMIN.  The new  RHO  is defined as  RHONEW  if it is less
C           than half the previous  RHO  and greater than  RHOMIN.
C           ------------------------------------------------------------
            SCALE = ONE
            DO 60 I = 1, NCNLN
               RHOMIN = MAX((WORK2(I)/QNORM)*TSCL,ZERO)
               RHOI = RHO(I)
C
               RHONEW = SQRT(RHOI*(RHODMP+RHOMIN))
               IF (RHONEW.LT.HALF*RHOI) RHOI = RHONEW
               IF (RHOI.LT.RHOMIN) RHOI = RHOMIN
               RHO(I) = RHOI
   60       CONTINUE
C
            RHO1 = RHONRM
            RHONRM = DNRM2(NCNLN,RHO,1)
C
C           ------------------------------------------------------------
C           If  INCRUN = .TRUE.,  there has been a run of iterations in
C           which the norm of  RHO  has not decreased.  Conversely,
C           INCRUN = false  implies that there has been a run of
C           iterations in which the norm of RHO has not increased.  If
C           INCRUN changes during this iteration the damping parameter
C           RHODMP is increased by a factor of two.  This ensures that
C           RHO(j) will oscillate only a finite number of times.
C           ------------------------------------------------------------
            BOOST = .FALSE.
            IF (INCRUN .AND. RHONRM.LT.RHO1) BOOST = .TRUE.
            IF ( .NOT. INCRUN .AND. RHONRM.GT.RHO1) BOOST = .TRUE.
            IF (BOOST) THEN
               RHODMP = TWO*RHODMP
               INCRUN = .NOT. INCRUN
            END IF
         END IF
C
         IF (NPDBG .AND. INPDBG(2).GT.0) THEN
            WRITE (REC,FMT=99997)
            CALL X04BAY(IPRINT,2,REC)
            DO 80 I = 1, NCNLN, 5
               WRITE (REC,FMT=99996) (RHO(L),L=I,MIN(I+4,NCNLN))
               CALL X04BAF(IPRINT,REC(1))
   80       CONTINUE
         END IF
C
      ELSE
C
C        The  QP  was infeasible.  Do not alter the penalty parameters,
C        but compute the scale factor so that the constraint violations
C        are reduced.
C
         CALL F06FCF(NCNLN,RHO,1,WORK1,1)
         PTERM2 = DDOT(NCNLN,WORK1,1,CS,1)
C
         SCALE = RHOMAX
         TSCL = F06BLF(GRDALF,PTERM2,OVERFL)
         IF (TSCL.GT.SCALE .AND. TSCL.LE.RHOMAX/(ONE+RHONRM)
     *       .AND. .NOT. OVERFL) SCALE = TSCL
C
         CALL DCOPY(NCNLN,CS,1,WORK1,1)
      END IF
C
C     ------------------------------------------------------------------
C     Compute the new value and directional derivative of the
C     merit function.
C     ------------------------------------------------------------------
      CALL F06FCF(NCNLN,RHO,1,WORK1,1)
C
      PTERM = DDOT(NCNLN,WORK1,1,CS,1)
      OBJALF = OBJALF + HALF*SCALE*PTERM
C
      IF (FEASQP) PTERM2 = PTERM
C
      GRDALF = GRDALF - SCALE*PTERM2
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99998) SCALE, RHONRM, GRDALF
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      RETURN
C
C
C     End of  E04UCN. (NPMRT)
C
99999 FORMAT (/' //E04UCN //        QPCURV        GRDALF ',/' //E04UCN',
     *       ' //',1P,2D14.2)
99998 FORMAT (/' //E04UCN //         SCALE        RHONRM        GRDALF '
     *       ,/' //E04UCN //',1P,3D14.2)
99997 FORMAT (/' //E04UCN //  Penalty parameters =       ')
99996 FORMAT (1P,5D15.6)
      END
      SUBROUTINE E04XAY(INFORM,MSGLVL,LVLDER,N,NCNLN,LDCJ,LDCJU,BIGBND,
     *                  EPSRF,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,
     *                  C1,C2,CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,Y,
     *                  IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1107 (JUL 1993).
C
C     ******************************************************************
C     E04XAY  computes difference intervals for the missing gradients of
C     F(x) and c(x). Intervals are computed using a procedure that
C     usually requires about two function evaluations if the function
C     is well scaled.  Central-difference gradients are obtained as a
C     by-product of the algorithm.
C
C     On entry...
C     OBJF and C contain the problem functions at the point X.
C     An element of CJAC or GRAD not equal to RDUMMY signifies a known
C     gradient value.  Such values are not estimated by differencing.
C     CJACU and GRADU have dummy elements in the same positions as
C     CJAC and GRADU.
C
C     On exit...
C     CJAC and GRAD contain central-difference derivative estimates.
C     Elements of CJACU and GRADU are unaltered except for those
C     corresponding to constant derivatives, which are given the same
C     values as CJAC or GRAD.
C
C     Systems Optimization Laboratory, Department of Operations Research
C     Stanford University, Stanford, California 94305
C     Original version written 28-July-1985.
C     This version of E04XAY   dated 13-Sep-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
      DOUBLE PRECISION  FACTOR
      PARAMETER         (FACTOR=0.97D+0)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TWO, FOUR, TEN
      PARAMETER         (TWO=2.0D+0,FOUR=4.0D+0,TEN=1.0D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, EPSRF, FDNORM, OBJF
      INTEGER           INFORM, LDCJ, LDCJU, LVLDER, MSGLVL, N, NCNLN
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), C2(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), GRAD(N), GRADU(N), HCNTRL(*),
     *                  HFORWD(*), USER(*), X(N), Y(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LFDSET, LINES1, LINES2, LVLDIF,
     *                  NCDIFF, NFDIFF, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CDEST, CJDIFF, D, DX, EPSA,
     *                  ERRBND, ERRMAX, ERRMIN, F1, F2, FDEST, FX,
     *                  GDIFF, H, HCD, HFD, HMAX, HMIN, HOPT, HPHI,
     *                  OBJF2, SDEST, SIGNH, STEPBL, STEPBU, SUMEPS,
     *                  SUMSD, TEST, XJ, YJ
      INTEGER           I, INFO, IROW1, IROW2, ITER, ITMAX, J, MODE,
     *                  NCCNST, NCOLJ, NFCNST, NSTATE
      LOGICAL           DEBUG, DONE, FIRST, HEADNG, NEEDED
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          E04XAZ, F06DBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
      INFORM = 0
      NEEDED = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR. LVLDER .EQ. 1 .AND.
     *         NCNLN .GT. 0
      IF ( .NOT. NEEDED) RETURN
C
      DEBUG = NPDBG .AND. INPDBG(5) .GT. 0
      IF (LFDSET.EQ.0) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAY(IPRINT,4,REC)
         END IF
C
         NSTATE = 0
         ITMAX = 3
         MODE = 0
C
         NCCNST = 0
         NFCNST = 0
         HEADNG = .TRUE.
C
         FDNORM = ZERO
C
C        ===============================================================
C        For each column of the Jacobian augmented by the transpose of
C        the objective gradient, rows IROW1 thru IROW2 are searched for
C        missing elements.
C        ===============================================================
         IROW1 = 1
         IROW2 = NCNLN + 1
         IF (LVLDER.EQ.1) IROW2 = NCNLN
         IF (LVLDER.EQ.2) IROW1 = NCNLN + 1
C
         BIGLOW = -BIGBND
         BIGUPP = BIGBND
C
         IF (NCNLN.GT.0) CALL F06DBF(NCNLN,(0),NEEDC,1)
C
         DO 60 J = 1, N
            XJ = X(J)
            NCOLJ = 0
            SUMSD = ZERO
            SUMEPS = ZERO
            HFD = ZERO
            HCD = ZERO
            HMAX = ZERO
            HMIN = ONE/EPSPT3
            ERRMAX = ZERO
            ERRMIN = ZERO
C
            STEPBL = BIGLOW
            STEPBU = BIGUPP
            IF (BL(J).GT.BIGLOW) STEPBL = BL(J) - XJ
            IF (BU(J).LT.BIGUPP) STEPBU = BU(J) - XJ
C
            SIGNH = ONE
            IF (HALF*(STEPBL+STEPBU).LT.ZERO) SIGNH = -ONE
C
            DO 40 I = IROW1, IROW2
C
               IF (I.LE.NCNLN) THEN
                  TEST = CJACU(I,J)
               ELSE
                  TEST = GRADU(J)
               END IF
C
               IF (TEST.EQ.RDUMMY) THEN
C                 ======================================================
C                 Get the difference interval for this element.
C                 ======================================================
                  NCOLJ = NCOLJ + 1
C
                  IF (I.LE.NCNLN) THEN
                     NEEDC(I) = 1
                     FX = C(I)
                     EPSA = EPSRF*(ONE+ABS(C(I)))
                  ELSE
                     FX = OBJF
                     EPSA = EPSRF*(ONE+ABS(FX))
                  END IF
C
C                 ------------------------------------------------------
C                 Find a finite-difference interval by iteration.
C                 ------------------------------------------------------
                  ITER = 0
                  HOPT = TWO*(ONE+ABS(XJ))*SQRT(EPSRF)
                  H = SIGNH*TEN*HOPT
                  CDEST = ZERO
                  SDEST = ZERO
                  FIRST = .TRUE.
C
C                 +                REPEAT
   20             X(J) = XJ + H
                  IF (I.LE.NCNLN) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 200
                     F1 = C1(I)
                  ELSE
                     CALL OBJFUN(MODE,N,X,F1,GRADU,NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 200
                  END IF
C
                  X(J) = XJ + H + H
                  IF (I.LE.NCNLN) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 200
                     F2 = C1(I)
                  ELSE
                     CALL OBJFUN(MODE,N,X,F2,GRADU,NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 200
                  END IF
C
                  CALL E04XAZ(DEBUG,DONE,FIRST,EPSA,EPSRF,FX,INFO,ITER,
     *                        ITMAX,CDEST,FDEST,SDEST,ERRBND,F1,F2,H,
     *                        HOPT,HPHI)
C
C                 +                UNTIL     DONE
                  IF ( .NOT. DONE) GO TO 20
C
                  IF (I.LE.NCNLN) THEN
                     CJAC(I,J) = CDEST
                     IF (INFO.EQ.1 .OR. INFO.EQ.2) THEN
                        NCCNST = NCCNST + 1
                        NCDIFF = NCDIFF - 1
                        CJACU(I,J) = -RDUMMY
                     END IF
                  ELSE
                     GRAD(J) = CDEST
                     IF (INFO.EQ.1 .OR. INFO.EQ.2) THEN
                        NFCNST = NFCNST + 1
                        NFDIFF = NFDIFF - 1
                        GRADU(J) = -RDUMMY
                     END IF
                  END IF
C
                  SUMSD = SUMSD + ABS(SDEST)
                  SUMEPS = SUMEPS + EPSA
                  IF (HOPT.GT.HMAX) THEN
                     HMAX = HOPT
                     ERRMAX = ERRBND
                  END IF
                  IF (HOPT.LT.HMIN) THEN
                     HMIN = HOPT
                     ERRMIN = ERRBND
                  END IF
C
                  IF (INFO.EQ.0) HCD = MAX(HCD,HPHI)
               END IF
   40       CONTINUE
C
            IF (NCOLJ.GT.0) THEN
               IF (HMIN.GT.HMAX) THEN
                  HMIN = HMAX
                  ERRMIN = ERRMAX
               END IF
C
               IF (FOUR*SUMEPS.LT.HMIN*HMIN*SUMSD) THEN
                  HFD = HMIN
                  ERRMAX = ERRMIN
               ELSE IF (FOUR*SUMEPS.GT.HMAX*HMAX*SUMSD) THEN
                  HFD = HMAX
               ELSE
                  HFD = TWO*SQRT(SUMEPS/SUMSD)
                  ERRMAX = TWO*SQRT(SUMEPS*SUMSD)
               END IF
C
               IF (HCD.EQ.ZERO) HCD = TEN*HFD
C
               IF (MSGLVL.GT.0) THEN
                  IF (HEADNG) THEN
                     WRITE (REC,FMT=99998)
                     CALL X04BAY(IPRINT,4,REC)
                  END IF
                  WRITE (REC,FMT=99997) J, XJ, HFD, HCD, ERRMAX
                  CALL X04BAF(IPRINT,REC(1))
                  HEADNG = .FALSE.
               END IF
               FDNORM = MAX(FDNORM,HFD)
               HFORWD(J) = HFD/(ONE+ABS(XJ))
               HCNTRL(J) = HCD/(ONE+ABS(XJ))
            END IF
            X(J) = XJ
   60    CONTINUE
C
         IF (NCCNST+NFCNST.GT.0) THEN
C
C           Check that the constants have been set properly by
C           evaluating the gradients at a strange (but feasible) point.
C
            D = ONE/N
C
            DO 80 J = 1, N
               XJ = X(J)
               STEPBL = -ONE
               STEPBU = ONE
               IF (BL(J).GT.BIGLOW) STEPBL = MAX(STEPBL,BL(J)-XJ)
               IF (BU(J).LT.BIGUPP .AND. BU(J).GT.BL(J))
     *             STEPBU = MIN(STEPBU,BU(J)-XJ)
C
               IF (HALF*(STEPBL+STEPBU).LT.ZERO) THEN
                  Y(J) = XJ + D*STEPBL
               ELSE
                  Y(J) = XJ + D*STEPBU
               END IF
C
               D = FACTOR*D
   80       CONTINUE
C
            IF (NCNLN.GT.0) THEN
               CALL F06DBF(NCNLN,(1),NEEDC,1)
               CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C2,CJACU,NSTATE,
     *                     IUSER,USER)
               IF (MODE.LT.0) GO TO 200
            END IF
C
            CALL OBJFUN(MODE,N,Y,OBJF2,GRADU,NSTATE,IUSER,USER)
            IF (MODE.LT.0) GO TO 200
C
C           ------------------------------------------------------------
C           Loop over each of the elements of  x.
C           ------------------------------------------------------------
            DO 140 J = 1, N
               YJ = Y(J)
               DX = HALF*(X(J)-YJ)
               Y(J) = YJ + DX
C
               IF (NCNLN.GT.0) THEN
                  NCOLJ = 0
                  DO 100 I = 1, NCNLN
                     IF (CJACU(I,J).EQ.-RDUMMY) THEN
                        NEEDC(I) = 1
                        NCOLJ = NCOLJ + 1
                     ELSE
                        NEEDC(I) = 0
                     END IF
  100             CONTINUE
C
                  IF (NCOLJ.GT.0) THEN
                     CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,Y,C1,CJACU,
     *                           NSTATE,IUSER,USER)
                     IF (MODE.LT.0) GO TO 200
C
                     DO 120 I = 1, NCNLN
                        IF (NEEDC(I).EQ.1) THEN
                           CJDIFF = (C1(I)-C2(I))/DX
                           IF (CJDIFF.EQ.CJAC(I,J)) THEN
                              CJACU(I,J) = CJDIFF
                           ELSE
                              CJACU(I,J) = RDUMMY
                              NCCNST = NCCNST - 1
                              NCDIFF = NCDIFF + 1
                           END IF
                        END IF
  120                CONTINUE
                  END IF
               END IF
C
C              Now check the objective gradient element.
C
               IF (GRADU(J).EQ.-RDUMMY) THEN
C
                  CALL OBJFUN(MODE,N,Y,F1,GRADU,NSTATE,IUSER,USER)
                  IF (MODE.LT.0) GO TO 200
C
                  GDIFF = (F1-OBJF2)/DX
                  IF (GDIFF.EQ.GRAD(J)) THEN
                     GRADU(J) = GDIFF
                  ELSE
                     GRADU(J) = RDUMMY
                     NFDIFF = NFDIFF + 1
                     NFCNST = NFCNST - 1
                  END IF
               END IF
C
               Y(J) = YJ
  140       CONTINUE
C
            IF (MSGLVL.GT.0) THEN
               IF (LVLDER.LT.2 .AND. NCCNST.GT.0) THEN
                  WRITE (REC,FMT=99996) NCCNST
                  CALL X04BAY(IPRINT,2,REC)
               END IF
               IF (LVLDER.NE.1 .AND. NFCNST.GT.0) THEN
                  WRITE (REC,FMT=99995) NFCNST
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            END IF
C
            IF (NCDIFF.EQ.0 .AND. LVLDER.LT.2) THEN
               IF (LVLDER.EQ.0) LVLDER = 2
               IF (LVLDER.EQ.1) LVLDER = 3
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99994) LVLDER
                  CALL X04BAY(IPRINT,4,REC)
               END IF
            END IF
C
            IF (NFDIFF.EQ.0 .AND. LVLDER.NE.1) THEN
               IF (LVLDER.EQ.0) LVLDER = 1
               IF (LVLDER.EQ.2) LVLDER = 3
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99993) LVLDER
                  CALL X04BAY(IPRINT,4,REC)
               END IF
            END IF
         END IF
      ELSE IF (LFDSET.EQ.2) THEN
C
C        The user has supplied HFORWD and HCNTRL.
C        Check for wild values.
C
         DO 160 J = 1, N
            IF (HFORWD(J).LE.ZERO) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99992) J, HFORWD(J), EPSPT5
                  CALL X04BAF(IPRINT,REC(1))
               END IF
               HFORWD(J) = EPSPT5
            END IF
  160    CONTINUE
         DO 180 J = 1, N
            IF (HCNTRL(J).LE.ZERO) THEN
               IF (MSGLVL.GT.0) THEN
                  WRITE (REC,FMT=99991) J, HCNTRL(J), EPSPT3
                  CALL X04BAF(IPRINT,REC(1))
               END IF
               HCNTRL(J) = EPSPT3
            END IF
  180    CONTINUE
      END IF
C
      RETURN
C
  200 INFORM = MODE
      RETURN
C
C
C     End of  E04XAY. (CHFD)
C
99999 FORMAT (//' Computation of the finite-difference intervals',/' -',
     *       '---------------------------------------------')
99998 FORMAT (//'    J      X(J)   Forward DX(J)   Central DX(J)      ',
     *       'Error est.',/)
99997 FORMAT (I5,1P,D10.2,1P,D16.6,1P,2D16.6)
99996 FORMAT (/I5,'  constant constraint gradient elements assigned.')
99995 FORMAT (/I5,'  constant  objective gradient elements assigned.')
99994 FORMAT (//' All missing Jacobian elements are constants.',/' Der',
     *       'ivative level increased to ',I4)
99993 FORMAT (//' All missing objective gradients are constants.',/' D',
     *       'erivative level increased to ',I4)
99992 FORMAT (' XXX  ',I4,'-th difference interval ',1P,D10.2,' replac',
     *       'ed by ',1P,D10.2)
99991 FORMAT (' XXX  ',I4,'-th central-difference interval ',1P,D10.2,
     *       ' replaced by ',1P,D10.2)
      END
      SUBROUTINE E04UDW(STRING,FIRST,LAST,MARK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Description and usage:
C
C       Looks for non-blank fields ('tokens') in a string, where the
C     fields are of arbitrary length, separated by blanks, tabs, commas,
C     colons, or equal signs.  The position of the end of the 1st token
C     is also returned, so this routine may be conveniently used within
C     a loop to process an entire line of text.
C
C       The procedure examines a substring, STRING (FIRST : LAST), which
C     may of course be the entire string (in which case just call E04UDW
C     with FIRST .LE. 1 and LAST .GE. LEN (STRING) ).  The indices
C     returned are relative to STRING itself, not the substring.
C
C
C     Parameters:
C
C     Name    Dimension  Type  I/O/S  Description
C     STRING              C    I      Text string containing data to be
C                                    scanned.
C     FIRST               I    I/O    Index of beginning of substring.
C                                    If .LE. 1, the search begins with
C                                    1.
C                                    Output is index of beginning of
C                                    first non-blank field, or 0 if no
C                                    token was found.
C     LAST                I    I/O    Index of end of substring.
C                                    If .GE. LEN (STRING), the search
C                                    begins with LEN (STRING).  Output
C                                    is index of end of last non-blank
C                                    field, or 0 if no token was found.
C     MARK                I      O    Points to end of first non-blank
C                                    field in the specified substring.
C                                    Set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C               ANSI Fortran 77, except for the tab character HT.
C
C     Notes:
C
C     (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C         in a non-standard way:  the CHAR function is not permitted
C         in a PARAMETER declaration (OK on VAX, though).  For Absoft
C         FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C         may be best to declare HT as a variable and assign
C         HT = CHAR(9) on ASCII machines, or CHAR(5) for EBCDIC.
C
C     (2)  The pseudo-recursive structure was chosen for fun.  It is
C         equivalent to three DO loops with embedded GO TOs in sequence.
C
C     (3)  The variety of separators recognized limits the usefulness of
C         this routine somewhat.  The intent is to facilitate handling
C         such tokens as keywords or numerical values.  In other
C         applications, it may be necessary for ALL printing characters
C         to be significant.  A simple modification to statement
C         function SOLID will do the trick.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C     29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                           based on SCAN_STRING by Ralph Carmichael.
C     25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C     16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                           (previous re-use of STRING was ambiguous).
C
C-----------------------------------------------------------------------
C
C     .. Parameters ..
      CHARACTER         BLANK, EQUAL, COLON, COMMA, RPARN, LPARN
      PARAMETER         (BLANK=' ',EQUAL='=',COLON=':',COMMA=',',
     *                  RPARN=')',LPARN='(')
C     .. Scalar Arguments ..
      INTEGER           FIRST, LAST, MARK
      CHARACTER*(*)     STRING
C     .. Local Scalars ..
      INTEGER           BEGIN, END, LENGTH
      CHARACTER         DUMMY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, LEN
C     .. Statement Functions ..
      LOGICAL           SOLID
C     .. Statement Function definitions ..
      SOLID(DUMMY) = (DUMMY.NE.BLANK) .AND. (DUMMY.NE.COLON)
     *               .AND. (DUMMY.NE.COMMA) .AND. (DUMMY.NE.EQUAL)
     *               .AND. (DUMMY.NE.RPARN) .AND. (DUMMY.NE.LPARN)
C     .. Executable Statements ..
      MARK = 0
      LENGTH = LEN(STRING)
      BEGIN = MAX(FIRST,1)
      END = MIN(LENGTH,LAST)
C
C     Find the first significant character ...
C
      DO 60 FIRST = BEGIN, END, +1
         IF (SOLID(STRING(FIRST:FIRST))) THEN
C
C           ... then the end of the first token ...
C
            DO 40 MARK = FIRST, END - 1, +1
               IF ( .NOT. SOLID(STRING(MARK+1:MARK+1))) THEN
C
C                 ... and finally the last significant character.
C
                  DO 20 LAST = END, MARK, -1
                     IF (SOLID(STRING(LAST:LAST))) THEN
                        RETURN
                     END IF
   20             CONTINUE
C
C                 Everything past the first token was a separator.
C
                  LAST = LAST + 1
                  RETURN
               END IF
   40       CONTINUE
C
C           There was nothing past the first token.
C
            LAST = MARK
            RETURN
         END IF
   60 CONTINUE
C
C     Whoops - the entire substring STRING (BEGIN : END) was composed of
C     separators .
C
      FIRST = 0
      MARK = 0
      LAST = 0
      RETURN
C
C     End of  E04UDW. (OPSCAN)
      END
      SUBROUTINE E04UCW(N,NCLIN,NCNLN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,
     *                  NVIOL,AX,BL,BU,C,FEATOL,X,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1090 (JUL 1993).
C
C     ******************************************************************
C     E04UCW  computes the following...
C     (1)  The number of constraints that are violated by more
C          than  FEATOL  and the 2-norm of the constraint violations.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      April    1984.
C     This version of  E04UCW  dated  16-October-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CVNORM, ERRMAX
      INTEGER           JMAX, N, NCLIN, NCNLN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN),
     *                  C(*), FEATOL(N+NCLIN+NCNLN),
     *                  WORK(N+NCLIN+NCNLN), X(N)
      INTEGER           ISTATE(N+NCLIN+NCNLN)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES, TOLJ
      INTEGER           IS, J
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Compute NVIOL, the number of constraints violated by more than
C     FEATOL,  and CVNORM,  the 2-norm of the constraint
C     violations and residuals of the constraints in the QP working set.
C     ==================================================================
      NVIOL = 0
C
      DO 40 J = 1, N + NCLIN + NCNLN
         FEASJ = FEATOL(J)
         RES = ZERO
C
         IF (J.LE.N+NCLIN) THEN
C
C           Bound or general linear constraint.
C
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               CON = AX(J-N)
            END IF
C
            TOLJ = FEASJ
         ELSE
C
C           Nonlinear constraint.
C
            CON = C(J-N-NCLIN)
            TOLJ = ZERO
         END IF
C
C        Check for constraint violations.
C
         IF (BL(J).GT.BIGLOW) THEN
            RES = BL(J) - CON
            IF (RES.GT.FEASJ) NVIOL = NVIOL + 1
            IF (RES.GT.TOLJ) GO TO 20
         END IF
C
         IF (BU(J).LT.BIGUPP) THEN
            RES = BU(J) - CON
            IF (RES.LT.(-FEASJ)) NVIOL = NVIOL + 1
            IF (RES.LT.(-TOLJ)) GO TO 20
         END IF
C
C        This constraint is satisfied,  but count the residual as a
C        violation if the constraint is in the working set.
C
         IS = ISTATE(J)
C
         IF (IS.EQ.0) THEN
            RES = ZERO
         ELSE IF (IS.EQ.1 .OR. IS.LE.-2) THEN
            RES = BL(J) - CON
         ELSE IF (IS.GE.2 .OR. IS.EQ.-1) THEN
            RES = BU(J) - CON
         END IF
C
         IF (ABS(RES).GT.FEASJ) NVIOL = NVIOL + 1
C
C        Set the array of violations.
C
   20    WORK(J) = RES
   40 CONTINUE
C
      JMAX = IDAMAX(N+NCLIN+NCNLN,WORK,1)
      ERRMAX = ABS(WORK(JMAX))
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) ERRMAX, JMAX
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      CVNORM = DNRM2(N+NCLIN+NCNLN,WORK,1)
C
      RETURN
C
C
C     End of  E04UCW. (NPFEAS)
C
99999 FORMAT (/' //E04UCW//  The maximum violation is ',1P,D14.2,' in ',
     *       'constraint',I5)
      END
      SUBROUTINE E04UCR(NEEDFD,INFORM,N,NCNLN,LDCJ,LDCJU,NFUN,NGRAD,
     *                  NEEDC,CONFUN,OBJFUN,ALFA,ALFBND,ALFMAX,ALFSML,
     *                  DXNORM,EPSRF,ETA,GDX,GRDALF,GLF1,GLF,OBJF,
     *                  OBJALF,QPCURV,XNORM,C,C2,CJAC,CJACU,CJDX,CJDX2,
     *                  CMUL1,CMUL,CS1,CS,DX,DLAM,DSLK,GRAD,GRADU,QPMUL,
     *                  RHO,SLK1,SLK,X1,X,WORK,W,LENW,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1086 (JUL 1993).
C
C     ==================================================================
C     E04UCR finds the steplength ALFA that gives sufficient decrease in
C     the augmented Lagrangian merit function.
C
C     On exit, if INFORM = 1, 2 or 3,  ALFA will be a nonzero steplength
C     with an associated merit function value  OBJALF  which is lower
C     than that at the base point. If  INFORM = 4, 5, 6, 7 or 8,  ALFA
C     is zero and  OBJALF will be the merit value at the base point.
C
C     Original version written  27-May-1985.
C     Level 2 BLAS added 12-June-1986.
C     This version of E04UCR dated  14-Sep-92.
C     ==================================================================
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D+0)
      DOUBLE PRECISION  TOLG
      PARAMETER         (TOLG=1.0D-1)
      DOUBLE PRECISION  RMU
      PARAMETER         (RMU=1.0D-4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBND, ALFMAX, ALFSML, DXNORM, EPSRF,
     *                  ETA, GDX, GLF, GLF1, GRDALF, OBJALF, OBJF,
     *                  QPCURV, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LENW, N, NCNLN, NFUN, NGRAD
      LOGICAL           NEEDFD
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), C2(*), CJAC(LDCJ,*), CJACU(LDCJU,*),
     *                  CJDX(*), CJDX2(*), CMUL(*), CMUL1(*), CS(*),
     *                  CS1(*), DLAM(*), DSLK(*), DX(N), GRAD(N),
     *                  GRADU(N), QPMUL(*), RHO(*), SLK(*), SLK1(*),
     *                  USER(*), W(LENW), WORK(*), X(N), X1(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9, RHODMP, RHOMAX,
     *                  RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LFDSET, LINES1, LINES2, LVLDIF,
     *                  NCDIFF, NFDIFF, NOUT
      LOGICAL           INCRUN, NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFBST, CS1JDX, CSJDX, CURVC, CURVLF, EPSAF,
     *                  EPSMCH, FBEST, FTERM, FTRY, G0, GBEST, GTRY,
     *                  OLDF, OLDG, Q, RHOBFS, S, T, TARGTG, TGDX, TGLF,
     *                  TOBJ, TOBJM, TOLABS, TOLAX, TOLREL, TOLRX,
     *                  TOLTNY
      INTEGER           J, MAXF, MODE, NSTATE, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04UCJ, E04UCK, F06DBF,
     *                  F06FCF, F06QFF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /NPDEBG/INPDBG, NPDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
      IF ( .NOT. NEEDFD .AND. NCNLN.GT.0) CS1JDX = DDOT(NCNLN,CS1,1,
     *    CJDX,1)
C
C     ------------------------------------------------------------------
C     Set the input parameters and tolerances for E04UCK and E04UCJ.
C
C     TOLRX   is the tolerance on relative changes in DX resulting from
C             changes in ALFA.
C
C     TOLAX   is the tolerance on absolute changes in DX resulting from
C             changes in ALFA.
C
C     TOLABS  is the tolerance on absolute changes in ALFA.
C
C     TOLREL  is the tolerance on relative changes in ALFA.
C
C     TOLTNY  is the magnitude of the smallest allowable value of ALFA.
C             if  M(TOLABS) - M(0) .gt. EPSAF,  the linesearch tries
C             steps in the range  TOLTNY .le. ALFA .le. TOLABS.
C     ------------------------------------------------------------------
      NSTATE = 0
      DEBUG = NPDBG .AND. INPDBG(4) .GT. 0
C
      IF (NEEDFD) THEN
         MAXF = 15
      ELSE
         MAXF = 10
      END IF
C
      EPSAF = EPSRF*(ONE+ABS(OBJALF))
      TOLAX = EPSPT8
      TOLRX = EPSPT8
C
      IF (TOLRX*XNORM+TOLAX.LT.DXNORM*ALFMAX) THEN
         TOLABS = (TOLRX*XNORM+TOLAX)/DXNORM
      ELSE
         TOLABS = ALFMAX
      END IF
      TOLREL = MAX(TOLRX,EPSMCH)
C
      T = ZERO
      DO 20 J = 1, N
         S = ABS(DX(J))
         Q = ABS(X(J))*TOLRX + TOLAX
         IF (S.GT.T*Q) T = S/Q
   20 CONTINUE
C
      IF (T*TOLABS.GT.ONE) THEN
         TOLTNY = ONE/T
      ELSE
         TOLTNY = TOLABS
      END IF
C
      OLDF = OBJALF
      OLDG = GRDALF
      ALFBST = ZERO
      FBEST = ZERO
      GBEST = (ONE-RMU)*OLDG
      TARGTG = (RMU-ETA)*OLDG
      G0 = GBEST
C
      IF (NCNLN.GT.0) CALL F06DBF(NCNLN,(1),NEEDC,1)
C
      IF (NEEDFD) THEN
         MODE = 0
      ELSE
         MODE = 2
      END IF
C
      FIRST = .TRUE.
C
C     ------------------------------------------------------------------
C     Commence main loop, entering E04UCK or E04UCJ two or more times.
C     FIRST = .TRUE. for the first entry, .FALSE. for subsequent entries
C     DONE  = .TRUE. indicates termination, in which case the value of
C     INFORM gives the result of the search.
C     INFORM = 1 if the search is successful and ALFA < ALFMAX.
C            = 2 if the search is successful and ALFA = ALFMAX.
C            = 3 if a better point was found but too many functions
C                were needed (not sufficient decrease).
C            = 4 if ALFMAX < TOLABS (too small to do a search).
C            = 5 if ALFA < ALFSML (E04UCJ only -- maybe want to switch
C                to central differences to get a better direction).
C            = 6 if the search found that there is no useful step.
C                The interval of uncertainty is less than 2*TOLABS.
C                The minimizer is very close to ALFA = zero
C                or the gradients are not sufficiently accurate.
C            = 7 if there were too many function calls.
C            = 8 if the input parameters were bad
C                (ALFMAX le TOLTNY  or  OLDG ge 0).
C     ------------------------------------------------------------------
C     +    repeat
   40 IF (NEEDFD) THEN
         CALL E04UCJ(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,ALFSML,EPSAF,G0,TARGTG,FTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST)
      ELSE
         CALL E04UCK(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,EPSAF,G0,TARGTG,FTRY,GTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST,GBEST)
      END IF
C
      IF (IMPRVD) THEN
         OBJF = TOBJ
         OBJALF = TOBJM
C
         IF (NCNLN.GT.0) CALL DCOPY(NCNLN,C2,1,C,1)
C
         IF ( .NOT. NEEDFD) THEN
            CALL DCOPY(N,GRADU,1,GRAD,1)
            GDX = TGDX
            GLF = TGLF
C
            IF (NCNLN.GT.0) THEN
               CALL DCOPY(NCNLN,CJDX2,1,CJDX,1)
               CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
            END IF
         END IF
      END IF
C
C     ---------------------------------------------------------------
C     If DONE = .FALSE.,  the problem functions must be computed for
C     the next entry to E04UCK or E04UCJ.
C     If DONE = .TRUE.,   this is the last time through.
C     ---------------------------------------------------------------
      IF ( .NOT. DONE) THEN
C
         CALL DCOPY(N,X1,1,X,1)
         CALL DAXPY(N,ALFA,DX,1,X,1)
C
         IF (NCNLN.GT.0) THEN
C
C           Compute the new estimates of the multipliers and slacks.
C           If the step length is greater than one,  the multipliers
C           are fixed as the QP-multipliers.
C
            IF (ALFA.LE.ONE) THEN
               CALL DCOPY(NCNLN,CMUL1,1,CMUL,1)
               CALL DAXPY(NCNLN,ALFA,DLAM,1,CMUL,1)
            END IF
            CALL DCOPY(NCNLN,SLK1,1,SLK,1)
            CALL DAXPY(NCNLN,ALFA,DSLK,1,SLK,1)
C
C           ---------------------------------------------------------
C           Compute the new constraint vector and Jacobian.
C           ---------------------------------------------------------
            CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,NSTATE,
     *                  IUSER,USER)
            IF (MODE.LT.0) GO TO 60
C
            CALL DCOPY(NCNLN,C2,1,CS,1)
            CALL DAXPY(NCNLN,(-ONE),SLK,1,CS,1)
C
            CALL DCOPY(NCNLN,CS,1,WORK,1)
            CALL F06FCF(NCNLN,RHO,1,WORK,1)
C
            FTERM = DDOT(NCNLN,CMUL,1,CS,1) - HALF*SCALE*DDOT(NCNLN,
     *              WORK,1,CS,1)
         END IF
C
C        ------------------------------------------------------------
C        Compute the value and gradient of the objective function.
C        ------------------------------------------------------------
         CALL OBJFUN(MODE,N,X,TOBJ,GRADU,NSTATE,IUSER,USER)
         IF (MODE.LT.0) GO TO 60
C
         IF (NCNLN.GT.0) THEN
            TOBJM = TOBJ - FTERM
         ELSE
            TOBJM = TOBJ
         END IF
C
         FTRY = TOBJM - OLDF - RMU*OLDG*ALFA
C
C        ------------------------------------------------------------
C        Compute auxiliary gradient information.
C        ------------------------------------------------------------
         IF ( .NOT. NEEDFD) THEN
            GTRY = DDOT(N,GRADU,1,DX,1)
            TGDX = GTRY
            TGLF = GTRY
            IF (NCNLN.GT.0) THEN
C
C              Compute the Jacobian times the search direction.
C
               CALL DGEMV('N',NCNLN,N,ONE,CJACU,LDCJU,DX,1,ZERO,CJDX2,1)
C
               CALL DCOPY(NCNLN,CJDX2,1,WORK,1)
               CALL DAXPY(NCNLN,(-ONE),DSLK,1,WORK,1)
C
               GTRY = GTRY - DDOT(NCNLN,CMUL,1,WORK,1)
               IF (ALFA.LE.ONE) GTRY = GTRY - DDOT(NCNLN,DLAM,1,CS,1)
C
               CALL F06FCF(NCNLN,RHO,1,WORK,1)
               GTRY = GTRY + SCALE*DDOT(NCNLN,WORK,1,CS,1)
C
               TGLF = TGDX - DDOT(NCNLN,CJDX2,1,QPMUL,1)
C
C              ------------------------------------------------------
C              If ALFBND .le. ALFA .lt. ALFMAX and the norm of the
C              quasi-Newton update is bounded, set ALFMAX to be ALFA.
C              This will cause the linesearch to stop if the merit
C              function is decreasing at the boundary.
C              ------------------------------------------------------
               IF (ALFBND.LE.ALFA .AND. ALFA.LT.ALFMAX) THEN
C
                  CSJDX = DDOT(NCNLN,CS,1,CJDX2,1)
C
                  IF (NPDBG .AND. INPDBG(1).GT.0) THEN
                     WRITE (REC,FMT=99997) CSJDX, CS1JDX, CURVLF
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
C
                  CURVLF = TGLF - GLF1
                  CURVC = ABS(CSJDX-CS1JDX)
                  RHOBFS = MAX(QPCURV*TOLG-CURVLF,ZERO)
                  IF (RHOBFS.LE.CURVC*RHOMAX) THEN
                     ALFMAX = ALFA
                  ELSE
                     ALFBND = MIN(TWO*ALFA,ALFMAX)
                  END IF
                  IF (NPDBG .AND. INPDBG(1).GT.0) THEN
                     WRITE (REC,FMT=99998) ALFBND, ALFA, ALFMAX
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
            END IF
C
            GTRY = GTRY - RMU*OLDG
C
         END IF
      END IF
C     +    until (      done)
      IF ( .NOT. DONE) GO TO 40
C
      NFUN = NFUN + NUMF
      IF ( .NOT. NEEDFD) NGRAD = NGRAD + NUMF
      ALFA = ALFBST
C
      IF ( .NOT. IMPRVD) THEN
         CALL DCOPY(N,X1,1,X,1)
         CALL DAXPY(N,ALFA,DX,1,X,1)
         IF (NCNLN.GT.0) THEN
            IF (ALFA.LE.ONE) THEN
               CALL DCOPY(NCNLN,CMUL1,1,CMUL,1)
               CALL DAXPY(NCNLN,ALFA,DLAM,1,CMUL,1)
            END IF
            CALL DCOPY(NCNLN,SLK1,1,SLK,1)
            CALL DAXPY(NCNLN,ALFA,DSLK,1,SLK,1)
            CALL DCOPY(NCNLN,C,1,CS,1)
            CALL DAXPY(NCNLN,(-ONE),SLK,1,CS,1)
         END IF
      END IF
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) INFORM
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      RETURN
C
C     The user wants to stop.  Who am I to object?
C
   60 INFORM = MODE
      RETURN
C
C
C     End of  E04UCR. (NPSRCH)
C
99999 FORMAT (/' //E04UCR// INFORM  = ',I4)
99998 FORMAT (/' //E04UCR//        ALFBND          ALFA        ALFMAX',
     *       /' //E04UCR//',1P,3D14.2)
99997 FORMAT (/' //E04UCR//         CSJDX        CS1JDX        CURVLF',
     *       /' //E04UCR//',1P,3D14.2)
      END
      SUBROUTINE E04UCL(LSUMRY,UNITQ,N,NCNLN,NFREE,NZ,LDCJ1,LDCJ2,LDZY,
     *                  LDR,KX,ALFA,GLF1,GLF2,QPCURV,CJAC1,CJAC2,CJDX1,
     *                  CJDX2,CS1,CS2,GQ1,GQ2,HPQ,RPQ,QPMUL,R,OMEGA,ZY,
     *                  WRK1,WRK2)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-880 (NOV 1990).
C     MARK 16 REVISED. IER-1081 (JUL 1993).
C
C     ******************************************************************
C     E04UCL  computes the BFGS update for the approximate Hessian of
C     the Lagrangian.  If the approximate curvature of the Lagrangian
C     function is negative,  a nonnegative penalty vector OMEGA(i) of
C     minimum two norm is computed such that the approximate curvature
C     of the augmented Lagrangian will be positive. If no finite penalty
C     vector exists,  the BFGS update is performed with the approximate
C     curvature modified to be a small positive value.
C
C     On entry,  GQ1 and GQ2 contain the transformed objective gradients
C     at X1 and X2,  HPQ contains  R'R(pq), the transformed Hessian
C     times the transformed search direction.  The vectors GQ1 and HPQ
C     are not saved.  If the regular BFGS quasi-Newton update could not
C     be performed, the first character of LSUMRY is loaded with 'M'.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written April 1984.
C     Level 2 BLAS added 12-June-1986.
C     This version of NPUPTD dated  10-May-1988.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  TOLG
      PARAMETER         (TOLG=1.0D-1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, GLF1, GLF2, QPCURV
      INTEGER           LDCJ1, LDCJ2, LDR, LDZY, N, NCNLN, NFREE, NZ
      LOGICAL           UNITQ
      CHARACTER*5       LSUMRY
C     .. Array Arguments ..
      DOUBLE PRECISION  CJAC1(LDCJ1,*), CJAC2(LDCJ2,*), CJDX1(*),
     *                  CJDX2(*), CS1(*), CS2(*), GQ1(N), GQ2(N),
     *                  HPQ(N), OMEGA(*), QPMUL(*), R(LDR,*), RPQ(N),
     *                  WRK1(N+NCNLN), WRK2(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DRMAX, DRMIN, RCNDBD, RFROBN, RHODMP, RHOMAX,
     *                  RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           INCRUN, NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA, CURVL, ETA, QI, QMAX, QNORM, RTGTP, RTYTS,
     *                  TEST, TINYCL, TRACE1, TRACE2
      INTEGER           I, IMAX, J
      LOGICAL           OVERFL, SSBFGS
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      INTEGER           IDAMAX
      EXTERNAL          DNRM2, F06BLF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DGEMV, DSCAL, E04NBV, E04NBW, F06FBF,
     *                  X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      IF (NCNLN.GT.0) CALL F06FBF(NCNLN,ZERO,OMEGA,1)
C
C     ------------------------------------------------------------------
C     Set CURVL = (G2 - G1)'DX,  the approximate curvature along DX of
C     the (augmented) Lagrangian.  At first, the curvature is not scaled
C     by the steplength ALFA.
C     ------------------------------------------------------------------
      CURVL = GLF2 - GLF1
      TINYCL = QPCURV*TOLG
      SSBFGS = CURVL .LE. ALFA*TINYCL
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) SSBFGS, TINYCL, CURVL
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
C     ------------------------------------------------------------------
C     Test if CURVL is sufficiently positive.  If there are no nonlinear
C     constraints,  no update can be performed.
C     ------------------------------------------------------------------
      IF (CURVL.LT.TINYCL) THEN
         LSUMRY(1:1) = 'Modified BFGS'
         IF (NCNLN.GT.0) THEN
            QMAX = ZERO
            DO 20 I = 1, NCNLN
               QI = CJDX2(I)*CS2(I) - CJDX1(I)*CS1(I)
               QMAX = MAX(QMAX,QI)
               IF (QI.LE.ZERO) WRK1(I) = ZERO
               IF (QI.GT.ZERO) WRK1(I) = QI
   20       CONTINUE
C
            QNORM = DNRM2(NCNLN,WRK1,1)
C
            TEST = MAX(TINYCL-CURVL,ZERO)
            BETA = F06BLF(QMAX*TEST,QNORM*QNORM,OVERFL)
            IF (BETA.LT.RHOMAX .AND. .NOT. OVERFL) THEN
               LSUMRY(1:1) = ' '
               BETA = TEST/(QNORM*QNORM)
               DO 40 I = 1, NCNLN
                  QI = WRK1(I)
                  OMEGA(I) = BETA*QI
                  CURVL = CURVL + BETA*QI*QI
   40          CONTINUE
C
               IF (NPDBG) THEN
                  IMAX = IDAMAX(NCNLN,OMEGA,1)
                  IF (INPDBG(1).GT.0) THEN
                     WRITE (REC,FMT=99997) OMEGA(IMAX)
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
C
                  IF (INPDBG(2).GT.0) THEN
                     WRITE (REC,FMT=99996)
                     CALL X04BAY(IPRINT,2,REC)
                     DO 60 J = 1, NCNLN, 5
                        WRITE (REC,FMT=99995) (OMEGA(I),I=J,
     *                    MIN(J+4,NCNLN))
                        CALL X04BAF(IPRINT,REC(1))
   60                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Compute the difference in the augmented Lagrangian gradient.
C     ------------------------------------------------------------------
C     Update GQ1 to include the augmented Lagrangian terms.
C
      IF (NCNLN.GT.0) THEN
C
         DO 80 I = 1, NCNLN
            WRK1(I) = -QPMUL(I) + OMEGA(I)*CS1(I)
   80    CONTINUE
         CALL DGEMV('T',NCNLN,N,ONE,CJAC1,LDCJ1,WRK1,1,ZERO,WRK2,1)
C
         DO 100 I = 1, NCNLN
            WRK1(I) = QPMUL(I) - OMEGA(I)*CS2(I)
  100    CONTINUE
         CALL DGEMV('T',NCNLN,N,ONE,CJAC2,LDCJ2,WRK1,1,ONE,WRK2,1)
C
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,WRK2,ZY,WRK1)
         CALL DAXPY(N,ONE,WRK2,1,GQ1,1)
      END IF
C
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99998) ALFA, CURVL
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (CURVL.LT.TINYCL) CURVL = TINYCL
C
      DO 120 J = 1, N
         WRK2(J) = GQ2(J) - GQ1(J)
  120 CONTINUE
C
      RTGTP = SQRT(QPCURV)
      RTYTS = SQRT(ALFA*CURVL)
      ETA = ONE
      IF (SSBFGS) ETA = RTYTS/(RTGTP*ALFA)
C
      TRACE1 = DNRM2(N,HPQ,1)/RTGTP
      TRACE2 = DNRM2(N,WRK2,1)/(RTYTS*ETA)
      RFROBN = ETA*SQRT(ABS((RFROBN-TRACE1)*(RFROBN+TRACE1)+TRACE2**2))
C
C     ==================================================================
C     Update the Cholesky factor of  Q'HQ.
C     ==================================================================
C     Normalize the vector  RPQ ( = R(pq) ).
C
      CALL DSCAL(N,(ONE/RTGTP),RPQ,1)
C
C     Do the self-scaled or regular BFGS update.
C     Form the vector WRK1 = gamma * (GQ2 - GQ1) - beta * R'R*PQ,
C     where  gamma = 1/SQRT( CURV ) = 1/SQRT( (GQ2 - GQ1)'SQ )
C
      CALL DSCAL(N,(ONE/RTGTP),HPQ,1)
C
      IF (SSBFGS) THEN
         DO 140 J = 1, N
            CALL DSCAL(J,ETA,R(1,J),1)
            WRK1(J) = WRK2(J)/RTYTS - ETA*HPQ(J)
  140    CONTINUE
      ELSE
         DO 160 J = 1, N
            WRK1(J) = WRK2(J)/RTYTS - HPQ(J)
  160    CONTINUE
      END IF
C
C     Perform the update to  R = R + RPQ*WRK1'.
C     RPQ is overwritten. Arrays GQ1 and HPQ are used to store the
C     sines and cosines defined by the plane rotations.
C
      CALL E04NBV(N,0,N,LDR,N,N,R,HPQ,RPQ,WRK1,GQ1,HPQ)
C
      RETURN
C
C
C     End of  E04UCL. (NPUPDT)
C
99999 FORMAT (/' //E04UCL// SSBFGS    min. CURVL         CURVL ',/' //',
     *       'E04UCL//   ',L4,1P,2D14.2)
99998 FORMAT (/' //E04UCL//          ALFA         CURVL ',
     *       /' //E04UCL//',1P,2D14.2)
99997 FORMAT (/' //E04UCL//   OMEGA(IMAX)',/' //E04UCL//',1P,D14.2)
99996 FORMAT (/' //E04UCL//  Penalty parameters = ')
99995 FORMAT (1P,5D15.6)
      END
      SUBROUTINE E04UCG(FIRSTV,HITLOW,ISTATE,INFORM,JADD,N,NCTOTL,
     *                  NUMINF,ALFA,PALFA,ATPHIT,BIGALF,BIGBND,PNORM,
     *                  ANORM,AP,AX,BL,BU,FEATOL,P,X)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1078 (JUL 1993).
C
C     ******************************************************************
C     E04UCG finds a step ALFA such that the point x + ALFA*P reaches
C     one of the linear constraints (including bounds).  Two possible
C     steps are defined as follows...
C
C     ALFA1   is the maximum step that can be taken without violating
C             one of the linear constraints that is currently satisfied.
C     ALFA2   reaches a linear constraint that is currently violated.
C             Usually this will be the furthest such constraint along P,
C             but if FIRSTV = .TRUE. it will be the first one along P.
C             This is used only when the problem has been determined to
C             be infeasible, and the sum of infeasibilities are being
C             minimized.  (ALFA2  is not defined if NUMINF = 0.)
C
C     ALFA will usually be the minimum of ALFA1 and ALFA2.
C     ALFA could be negative (since we allow inactive constraints
C     to be violated by as much as FEATOL).  In such cases, a
C     third possible step is computed, to find the nearest satisfied
C     constraint (perturbed by FEATOL) along the direction  - P.
C     ALFA  will be reset to this step if it is shorter.  This is the
C     only case for which the final step  ALFA  does not move X exactly
C     onto a constraint (the one denoted by JADD).
C
C     Constraints in the working set are ignored  (ISTATE(j) ge 1).
C
C     JADD    denotes which linear constraint is reached.
C
C     HITLOW  indicates whether it is the lower or upper bound that
C             has restricted ALFA.
C
C     Values of ISTATE(j)....
C
C     - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     The values -2 and -1 do not occur once a feasible point has been
C     found.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written  May 1980.
C     This version of  E04UCG  dated  10-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ATPHIT, BIGALF, BIGBND, PALFA, PNORM
      INTEGER           INFORM, JADD, N, NCTOTL, NUMINF
      LOGICAL           FIRSTV, HITLOW
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), P(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSATP, ALFA1, ALFA2, APMAX1, APMAX2, ATP, ATP1,
     *                  ATP2, ATX, PALFA1, PALFA2, RES, ROWNRM
      INTEGER           I, J, JADD1, JADD2, JS, JSAVE1, JSAVE2
      LOGICAL           HLOW1, HLOW2, LASTV, NEGSTP, STEP2
C     .. Local Arrays ..
      CHARACTER*120     REC(4)
C     .. External Subroutines ..
      EXTERNAL          E04UCH, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
      INFORM = 0
C
C     ------------------------------------------------------------------
C     First pass -- find steps to perturbed constraints, so that
C     PALFA1 will be slightly larger than the true step, and
C     PALFA2 will be slightly smaller than it should be.
C     In degenerate cases, this strategy gives us some freedom in the
C     second pass.  The general idea follows that described by P.M.J.
C     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
C     ------------------------------------------------------------------
C
      NEGSTP = .FALSE.
      CALL E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,PALFA1,
     *            PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,FEATOL,P,X)
C
      JSAVE1 = JADD1
      JSAVE2 = JADD2
C
C     ------------------------------------------------------------------
C     Second pass -- recompute step-lengths without perturbation.
C     Amongst constraints that are less than the perturbed steps,
C     choose the one (of each type) that makes the largest angle
C     with the search direction.
C     ------------------------------------------------------------------
      IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,4,REC)
      END IF
      ALFA1 = BIGALF
      ALFA2 = ZERO
      IF (FIRSTV) ALFA2 = BIGALF
C
      APMAX1 = ZERO
      APMAX2 = ZERO
      ATP1 = ZERO
      ATP2 = ZERO
      HLOW1 = .FALSE.
      HLOW2 = .FALSE.
      LASTV = .NOT. FIRSTV
C
      DO 20 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS.LE.0) THEN
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ROWNRM = ONE
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ROWNRM = ANORM(I) + ONE
            END IF
C
            IF (ABS(ATP).LE.EPSPT9*ROWNRM*PNORM) THEN
C
C              This constraint appears to be constant along P.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C
               RES = -ONE
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing.
C              ---------------------------------------------------------
C              The lower bound is satisfied.  Test for smaller ALFA1.
C
               ABSATP = -ATP
               IF (BL(J).GT.(-BIGBND)) THEN
                  RES = ATX - BL(J)
                  IF (PALFA1*ABSATP.GE.RES .OR. J.EQ.JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM.LT.ABSATP) THEN
                        APMAX1 = ABSATP/(ROWNRM*PNORM)
                        ALFA1 = RES/ABSATP
                        JADD1 = J
                        ATP1 = ATP
                        HLOW1 = .TRUE.
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 The upper bound is violated.  Test for either a bigger
C                 or smaller ALFA2,  depending on the value of FIRSTV.
C
                  RES = ATX - BU(J)
                  IF ((FIRSTV .AND. PALFA2*ABSATP.GE.RES .OR.
     *                LASTV .AND. PALFA2*ABSATP.LE.RES)
     *                .OR. J.EQ.JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM.LT.ABSATP) THEN
                        APMAX2 = ABSATP/(ROWNRM*PNORM)
                        IF (ABSATP.GE.ONE) THEN
                           ALFA2 = RES/ABSATP
                        ELSE IF (RES.LT.BIGALF*ABSATP) THEN
                           ALFA2 = RES/ABSATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2 = J
                        ATP2 = ATP
                        HLOW2 = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              Test for smaller ALFA1.
C
               IF (BU(J).LT.BIGBND) THEN
                  RES = BU(J) - ATX
                  IF (PALFA1*ATP.GE.RES .OR. J.EQ.JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM.LT.ATP) THEN
                        APMAX1 = ATP/(ROWNRM*PNORM)
                        ALFA1 = RES/ATP
                        JADD1 = J
                        ATP1 = ATP
                        HLOW1 = .FALSE.
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 The lower bound is violated.  Test for a new ALFA2.
C
                  RES = BL(J) - ATX
                  IF ((FIRSTV .AND. PALFA2*ATP.GE.RES .OR. LASTV .AND.
     *                PALFA2*ATP.LE.RES) .OR. J.EQ.JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM.LT.ATP) THEN
                        APMAX2 = ATP/(ROWNRM*PNORM)
                        IF (ATP.GE.ONE) THEN
                           ALFA2 = RES/ATP
                        ELSE IF (RES.LT.BIGALF*ATP) THEN
                           ALFA2 = RES/ATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2 = J
                        ATP2 = ATP
                        HLOW2 = .TRUE.
                     END IF
                  END IF
               END IF
            END IF
C
            IF (CMDBG .AND. ICMDBG(3).GT.0) THEN
               WRITE (REC,FMT=99998) J, JS, FEATOL(J), RES, ATP, JADD1,
     *           ALFA1, JADD2, ALFA2
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
   20 CONTINUE
C
C     ==================================================================
C     Determine ALFA, the step to be taken.
C     ==================================================================
C     In the infeasible case, check whether to take the step ALFA2
C     rather than ALFA1...
C
      STEP2 = NUMINF .GT. 0 .AND. JADD2 .GT. 0
C
C     We do so if ALFA2 is less than ALFA1 or (if FIRSTV is false)
C     lies in the range  (ALFA1, PALFA1)  and has a smaller value of
C     ATP.
C
      STEP2 = STEP2 .AND. (ALFA2.LT.ALFA1 .OR. LASTV .AND. ALFA2.LE.
     *        PALFA1 .AND. APMAX2.GE.APMAX1)
C
      IF (STEP2) THEN
         ALFA = ALFA2
         PALFA = PALFA2
         JADD = JADD2
         ATPHIT = ATP2
         HITLOW = HLOW2
      ELSE
         ALFA = ALFA1
         PALFA = PALFA1
         JADD = JADD1
         ATPHIT = ATP1
         HITLOW = HLOW1
C
C        If ALFA1 is negative, the constraint to be added (JADD)
C        remains unchanged, but ALFA may be shortened to the step
C        to the nearest perturbed satisfied constraint along  - P.
C
         NEGSTP = ALFA .LT. ZERO
         IF (NEGSTP) THEN
            CALL E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,
     *                  PALFA1,PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,
     *                  FEATOL,P,X)
C
            IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
               WRITE (REC,FMT=99997) ALFA, PALFA1
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
            ALFA = -MIN(ABS(ALFA),PALFA1)
         END IF
      END IF
C
C     Test for undefined or infinite step.
C
      IF (JADD.EQ.0) THEN
         ALFA = BIGALF
         PALFA = BIGALF
      END IF
C
      IF (ALFA.GE.BIGALF) INFORM = 3
      IF (CMDBG .AND. ICMDBG(1).GT.0 .AND. INFORM.GT.0) THEN
         WRITE (REC,FMT=99996) JADD, ALFA
         CALL X04BAY(IPRINT,4,REC)
      END IF
      RETURN
C
C
C     End of  E04UCG. (CMALF)
C
99999 FORMAT (/' E04UCG  entered',/'    J  JS         FEATOL        RE',
     *       'S             AP     JADD1        ALFA1     JADD2       ',
     *       ' ALFA2 ',/)
99998 FORMAT (I5,I4,3G15.5,2(I6,G17.7))
99997 FORMAT (/' //E04UCG //  Negative step',/' //E04UCG //           ',
     *       'ALFA          PALFA',/' //E04UCG //',2G15.4)
99996 FORMAT (/' //E04UCG //  Unbounded step.',/' //E04UCG //  JADD   ',
     *       '        ALFA',/' //E04UCG //  ',I4,G15.4)
      END
      subroutine e04ucq(nout,buffer,key)
c     mark 13 re-issue.  nag copyright 1988.
c     mark 16 revised. ier-1085 (jul 1993).
c
c     ******************************************************************
c     e04ucq   decodes the option contained in  buffer  in order to set
c     a parameter value in the relevant element of the parameter arrays.
c
c
c     input:
c
c     nout   a unit number for printing error messages.
c            nout  must be a valid unit.
c
c     output:
c
c     key    the first keyword contained in buffer.
c
c
c     e04ucq  calls e04udx and the subprograms
c                 lookup, scannr, tokens, upcase
c     (now called e04udy, e04udw, e04udv, e04udu)
c     supplied by informatics general, inc., palo alto, california.
c
c     systems optimization laboratory, stanford university.
c     this version of e04ucq  dated 19-oct-92.
c     ******************************************************************
c
c     .. parameters ..
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
c     .. external subroutines ..
      external          e04udv, e04udy, x04baf
c     .. intrinsic functions ..
      intrinsic         index, len
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
c     .. executable statements ..
c
      if (first) then
         first = .false.
         do 20 i = 1, mxparm
            rprmls(i) = rdummy
            iprmls(i) = idummy
            rprmnp(i) = rdummy
            iprmnp(i) = idummy
   20    continue
      end if
c
c     eliminate comments and empty lines.
c     a '*' appearing anywhere in buffer terminates the string.
c
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
c
c     ------------------------------------------------------------------
c     extract up to maxtok tokens from the record.
c     ntoken returns how many were actually found.
c     key, key2, key3 are the first tokens if any, otherwise blank.
c     ------------------------------------------------------------------
      call e04udv(buffer(1:lenbuf),maxtok,ntoken,token)
      key = token(1)
      key2 = token(2)
      key3 = token(3)
c
c     certain keywords require no action.
c
      if (key.eq.' ' .or. key.eq.'begin') go to 80
      if (key.eq.'list' .or. key.eq.'nolist') go to 80
      if (key.eq.'end') go to 80
c
c     most keywords will have an associated integer or real value,
c     so look for it no matter what the keyword.
c
      i = 1
      number = .false.
c
   40 if (i.lt.ntoken .and. .not. number) then
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
c     ------------------------------------------------------------------
c     decide what to do about each keyword.
c     the second keyword (if any) might be needed to break ties.
c     some seemingly redundant testing of more is used
c     to avoid compiler limits on the number of consecutive else ifs.
c     ------------------------------------------------------------------
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
      SUBROUTINE E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NRANK,NZ,
     *                  N,NCTOTL,LDZY,LDA,LDR,LDT,ISTATE,KACTIV,KX,JMAX,
     *                  ERRMAX,CTX,XNORM,A,AX,BL,BU,CQ,RES,RES0,FEATOL,
     *                  R,T,X,ZY,P,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1061 (JUL 1993).
C
C     ******************************************************************
C     E04NCH  computes the point on a working set that is closest to the
C     input vector  x  (in the least-squares sense).  The norm of  x,
C     the transformed residual vector  Pr - RQ'x,  and the constraint
C     values Ax  are also initialized.
C
C     If the computed point gives a row error of more than the
C     feasibility tolerance, an extra step of iterative refinement is
C     used.  If  x  is still infeasible,  the logical variable  ROWERR
C     is set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version dated 29-December-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           NTRY
      PARAMETER         (NTRY=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTX, ERRMAX, XNORM
      INTEGER           JMAX, LDA, LDR, LDT, LDZY, N, NACTIV, NCLIN,
     *                  NCTOTL, NFREE, NRANK, NZ
      LOGICAL           LINOBJ, ROWERR, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL), CQ(*),
     *                  FEATOL(NCTOTL), P(N), R(LDR,*), RES(*), RES0(*),
     *                  T(LDT,*), WORK(NCTOTL), X(N), ZY(LDZY,*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, IS, J, K, KTRY
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DDOT, DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRMV, E04NBT, E04NBW,
     *                  F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Move  x  onto the simple bounds in the working set.
C     ------------------------------------------------------------------
      DO 20 K = NFREE + 1, N
         J = KX(K)
         IS = ISTATE(J)
         BND = BL(J)
         IF (IS.GE.2) BND = BU(J)
         IF (IS.NE.4) X(J) = BND
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Move  x  onto the general constraints in the working set.
C     We shall make  NTRY  tries at getting acceptable row errors.
C     ------------------------------------------------------------------
      KTRY = 1
      JMAX = 1
      ERRMAX = ZERO
C
C     REPEAT
   40 IF (NACTIV.GT.0) THEN
C
C        Set  work = residuals for constraints in the working set.
C        Solve for p, the smallest correction to x that gives a point
C        on the constraints in the working set.  Define  p = Y*(py),
C        where  py  solves the triangular system  T*(py) = residuals.
C
         DO 60 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(I) = BND - DDOT(N,A(K,1),LDA,X,1)
   60    CONTINUE
C
         CALL E04NBT(1,LDT,NACTIV,T(1,NZ+1),WORK)
         CALL F06FBF(N,ZERO,P,1)
         CALL DCOPY(NACTIV,WORK,1,P(NZ+1),1)
C
         CALL E04NBW(2,N,NZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
         CALL DAXPY(N,ONE,P,1,X,1)
      END IF
C
C     ---------------------------------------------------------------
C     Compute the 2-norm of  x.
C     Initialize  Ax  for all the general constraints.
C     ---------------------------------------------------------------
      XNORM = DNRM2(N,X,1)
      IF (NCLIN.GT.0) CALL DGEMV('N',NCLIN,N,ONE,A,LDA,X,1,ZERO,AX,1)
C
C     ---------------------------------------------------------------
C     Check the row residuals.
C     ---------------------------------------------------------------
      IF (NACTIV.GT.0) THEN
         DO 80 K = 1, NACTIV
            I = KACTIV(K)
            J = N + I
            IS = ISTATE(J)
            IF (IS.EQ.1) WORK(K) = BL(J) - AX(I)
            IF (IS.GE.2) WORK(K) = BU(J) - AX(I)
   80    CONTINUE
C
         JMAX = IDAMAX(NACTIV,WORK,1)
         ERRMAX = ABS(WORK(JMAX))
      END IF
C
      KTRY = KTRY + 1
C     UNTIL    (ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY
      IF ( .NOT. (ERRMAX.LE.FEATOL(JMAX) .OR. KTRY.GT.NTRY)) GO TO 40
C
      ROWERR = ERRMAX .GT. FEATOL(JMAX)
C
C     ==================================================================
C     Compute the linear objective value  c'x  and the transformed
C     residual  Pr  -  RQ'x = RES0  -  RQ'x.
C     ==================================================================
      IF (NRANK.GT.0 .OR. LINOBJ) THEN
         CALL DCOPY(N,X,1,P,1)
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
      END IF
C
      CTX = ZERO
      IF (LINOBJ) CTX = DDOT(N,CQ,1,P,1)
C
      IF (NRANK.GT.0) THEN
         CALL DTRMV('U','N','N',NRANK,R,LDR,P,1)
         IF (NRANK.LT.N) CALL DGEMV('N',NRANK,N-NRANK,ONE,R(1,NRANK+1),
     *                              LDR,P(NRANK+1),1,ONE,P,1)
         CALL DCOPY(NRANK,RES0,1,RES,1)
         CALL DAXPY(NRANK,-ONE,P,1,RES,1)
      END IF
C
      IF (LSDBG .AND. ILSDBG(2).GT.0) THEN
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,2,REC)
         DO 100 I = 1, N, 5
            WRITE (REC,FMT=99998) (X(J),J=I,MIN(I+4,N))
            CALL X04BAF(IPRINT,REC(1))
  100    CONTINUE
      END IF
C
      RETURN
C
C
C     End of  E04NCH. (LSSETX)
C
99999 FORMAT (/' //E04NCH// Variables after refinement ... ')
99998 FORMAT (5G12.3)
      END
      SUBROUTINE E04NCY(UNITQ,VERTEX,INFORM,K1,K2,NACTIV,NARTIF,NZ,
     *                  NFREE,NRANK,NREJTD,NRES,NGQ,N,LDZY,LDA,LDR,LDT,
     *                  ISTATE,KACTIV,KX,CONDMX,A,R,T,RES,GQ,ZY,W,C,S,
     *                  MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1075 (JUL 1993).
C
C     ******************************************************************
C     E04NCY  includes general constraints K1 thru K2 as new rows of
C     the TQ factorization stored in T, ZY.  If NRANK is nonzero, the
C     changes in Q are reflected in NRANK by N triangular factor R such
C     that
C                         C  =  P ( R ) Q,
C                                 ( 0 )
C     where  P  is orthogonal.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  October-31-1984.
C     Level-2 matrix routines added 18-May-1988.
C     This version dated 27-April-1993.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           INFORM, K1, K2, LDA, LDR, LDT, LDZY, MSGLVL, N,
     *                  NACTIV, NARTIF, NFREE, NGQ, NRANK, NREJTD, NRES,
     *                  NZ
      LOGICAL           UNITQ, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNDMAX, RNORM, ROWMAX, RTMAX
      INTEGER           I, IADD, IARTIF, IFIX, ISWAP, JADD, K, L, NZADD
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          E04NCV, F06FLF
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
      RTMAX = WMACH(8)
C
C     Estimate the condition number of the constraints that are not
C     to be refactorized.
C
      IF (NACTIV.EQ.0) THEN
         DTMAX = ZERO
         DTMIN = ONE
      ELSE
         CALL F06FLF(NACTIV,T(NACTIV,NZ+1),LDT-1,DTMAX,DTMIN)
      END IF
C
      DO 20 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV.LT.NFREE) THEN
C
            CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,NFREE,
     *                  NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,A,R,
     *                  T,RES,GQ,ZY,W,C,S,MSGLVL)
C
            IF (INFORM.EQ.0) THEN
               NACTIV = NACTIV + 1
               NZ = NZ - 1
            ELSE
               ISTATE(JADD) = 0
               KACTIV(K) = -KACTIV(K)
            END IF
         END IF
   20 CONTINUE
C
      IF (NACTIV.LT.K2) THEN
C
C        Some of the constraints were classed as dependent and not
C        included in the factorization.  Re-order the part of  KACTIV
C        that holds the indices of the general constraints in the
C        working set.  Move accepted indices to the front and shift
C        rejected indices (with negative values) to the end.
C
         L = K1 - 1
         DO 40 K = K1, K2
            I = KACTIV(K)
            IF (I.GE.0) THEN
               L = L + 1
               IF (L.NE.K) THEN
                  ISWAP = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
   40    CONTINUE
C
C        If a vertex is required, add some temporary bounds.
C        We must accept the resulting condition number of the working
C        set.
C
         IF (VERTEX) THEN
            CNDMAX = RTMAX
            NZADD = NZ
            DO 80 IARTIF = 1, NZADD
               IF (UNITQ) THEN
                  IFIX = NFREE
                  JADD = KX(IFIX)
               ELSE
                  ROWMAX = ZERO
                  DO 60 I = 1, NFREE
                     RNORM = DNRM2(NZ,ZY(I,1),LDZY)
                     IF (ROWMAX.LT.RNORM) THEN
                        ROWMAX = RNORM
                        IFIX = I
                     END IF
   60             CONTINUE
                  JADD = KX(IFIX)
C
                  CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,
     *                        NFREE,NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,
     *                        KX,CNDMAX,A,R,T,RES,GQ,ZY,W,C,S,MSGLVL)
C
               END IF
               NFREE = NFREE - 1
               NZ = NZ - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
   80       CONTINUE
         END IF
      END IF
C
      NREJTD = K2 - NACTIV
C
      RETURN
C
C     End of  E04NCY. (LSADDS)
C
      END
      SUBROUTINE E04MFZ(PRBTYP,MSG,CSET,NAMED,NAMES,RSET,UNITQ,ITER,
     *                  ITMAX,JINF,NVIOL,N,NCLIN,LDA,NACTIV,NFREE,NRZ,
     *                  NZ,ISTATE,KACTIV,KX,E04MFU,OBJ,NUMINF,XNORM,A,
     *                  AX,BL,BU,CVEC,FEATOL,FEATLU,X,IW,W)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFZ  is a subroutine for linear programming.
C     On entry, it is assumed that an initial working set of
C     linear constraints and bounds is available.  The arrays  ISTATE,
C     KACTIV  and  KX  will have been set accordingly
C     and the arrays  T  and  Q  will contain the TQ factorization of
C     the matrix whose rows are the gradients of the active linear
C     constraints with the columns corresponding to the active bounds
C     removed.  The TQ factorization of the resulting (NACTIV by NFREE)
C     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
C     is upper-triangular.
C
C     Over a cycle of iterations, the feasibility tolerance FEATOL
C     increases slightly (from tolx0 to tolx1 in steps of tolinc).
C     this ensures that all steps taken will be positive.
C
C     After KDEGEN consecutive iterations, variables within FEATOL of
C     their bounds are set exactly on their bounds and iterative
C     refinement is used to satisfy the constraints in the working set.
C     FEATOL is then reduced to tolx0 for the next cycle of iterations.
C
C
C     Values of ISTATE(j) for the linear constraints.......
C
C     ISTATE(j)
C     ---------
C          0    constraint j is not in the working set.
C          1    constraint j is in the working set at its lower bound.
C          2    constraint j is in the working set at its upper bound.
C          3    constraint j is in the working set as an equality.
C
C     Constraint j may be violated by as much as FEATOL(j).
C
C     This version of  E04MFZ  dated  13-Dec-91.
C
C     Copyright  1988  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       EMPTY
      PARAMETER         (EMPTY='      ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ, XNORM
      INTEGER           ITER, ITMAX, JINF, LDA, N, NACTIV, NCLIN, NFREE,
     *                  NRZ, NUMINF, NVIOL, NZ
      LOGICAL           CSET, NAMED, RSET, UNITQ
      CHARACTER*2       PRBTYP
      CHARACTER*6       MSG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CVEC(*), FEATLU(N+NCLIN), FEATOL(N+NCLIN), W(*),
     *                  X(N)
      INTEGER           ISTATE(N+NCLIN), IW(*), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Subroutine Arguments ..
      EXTERNAL          E04MFU
C     .. Scalars in Common ..
      DOUBLE PRECISION  ALFA, ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9,
     *                  TOLACT, TOLFEA, TOLINC, TOLRNK, TOLX0, TRULAM
      INTEGER           IDBGLC, IPRINT, IPRNT, ISDEL, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITNFIX, JADD, JDEL, KCHK,
     *                  KCYCLE, KDEGEN, LCRASH, LDBGLC, LDQ, LDT,
     *                  LENNAM, LINES1, LINES2, LPROB, MAXACT, MAXNZ,
     *                  MM, MSGLC, MXFREE, NCOLT, NDEGEN, NN, NNCLIN,
     *                  NOUT, NPROB
      LOGICAL           CMDBG, HEADER, LCDBG, PRNT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), ILCDBG(LDBG), IPADLC(14),
     *                  IPSVLC(MXPARM), LOCLC(LENLC), NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFAP, ALFHIT, BIGALF, BIGGST, CONDMX, CONDRZ,
     *                  CONDT, DINKY, DNORM, DZZ, EPSMCH, ERRMAX, FLMAX,
     *                  GFNORM, GRZNRM, GZNORM, OBJSIZ, RTMAX, SMLLST,
     *                  SUMINF, TINYST, TRUBIG, TRUSML, WSSIZE, ZEROLM
      INTEGER           IADD, IDBG, IFIX, INFORM, IS, IT, J, JBIGST,
     *                  JMAX, JSMLST, JTINY, KBIGST, KDEL, KSMLST, LAD,
     *                  LANORM, LCQ, LD, LDR, LGQ, LQ, LR, LRLAM, LT,
     *                  LTMP, LWRK, LWTINF, MSGDBG, MSGLVL, MSGSVD,
     *                  NCTOTL, NFIXED, NGQ, NMOVED, NOTOPT, NTFIXD
      LOGICAL           FIRSTV, FP, HITLOW, LP, MOVE, ONBND, OVERFL,
     *                  UNBNDD
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, E04MFH, E04MFL,
     *                  E04MFM, E04MFN, E04MFQ, E04MFR, E04MFS, E04MFY,
     *                  E04NBW, E04NFP, E04NFR, F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MOD
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /CE04MF/TOLX0, TOLINC, KDEGEN, NDEGEN, ITNFIX,
     *                  NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04MF/ALFA, TRULAM, ISDEL, JDEL, JADD, HEADER,
     *                  PRNT
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /EE04MF/ILCDBG, LCDBG
      COMMON            /FE04MF/IPSVLC, IDBGLC, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, KCHK, KCYCLE, LCRASH, LPROB, MAXACT,
     *                  MXFREE, MAXNZ, MM, LDBGLC, MSGLC, NN, NNCLIN,
     *                  NPROB, IPADLC
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLC(1),IDBGLC), (RPRMLC(1),BIGBND)
      EQUIVALENCE       (MSGLC,MSGLVL), (IDBGLC,IDBG), (LDBGLC,MSGDBG)
C     .. Save statement ..
      SAVE              /AX02ZA/, /AE04MF/, /FE04MF/, /GE04MF/, FIRSTV
C     .. Executable Statements ..
C
C     Specify the machine-dependent parameters.
C
      EPSMCH = WMACH(3)
      FLMAX = WMACH(7)
      RTMAX = WMACH(8)
C
      IF (CSET) THEN
         NGQ = 2
      ELSE
         NGQ = 1
      END IF
C
      LP = PRBTYP .EQ. 'lp' .OR. PRBTYP .EQ. 'LP'
      FP = .NOT. LP
C
      LDR = LDT
      IT = 1
C
      LANORM = LOCLC(4)
      LAD = LOCLC(5)
C
      LD = LOCLC(7)
      LGQ = LOCLC(8)
      LCQ = LOCLC(9)
      LRLAM = LOCLC(10)
C
      LR = LOCLC(11)
      LT = LOCLC(12)
      LQ = LOCLC(13)
      LWTINF = LOCLC(14)
      LWRK = LOCLC(15)
C
C     We need a temporary array when changing the active set.
C     Use the multiplier array.
C
      LTMP = LRLAM
C
      IF (ITER.EQ.0) THEN
C        -------------------------
C        First entry.  Initialize.
C        -------------------------
         JADD = 0
         JDEL = 0
         ISDEL = 0
         FIRSTV = .FALSE.
C
         ALFA = ZERO
         DZZ = ONE
      END IF
C
      NCTOTL = N + NCLIN
      NVIOL = 0
C
      CONDMX = FLMAX
C
C     If debug output is required, print nothing until iteration IDBG.
C
      MSGSVD = MSGLVL
      IF (IDBG.GT.ITER .AND. IDBG.LE.ITMAX) THEN
         MSGLVL = 0
      END IF
C
      CALL E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,BU,A,
     *            FEATOL,W(LGQ),X,W(LWTINF))
C
      IF (NUMINF.GT.0) THEN
         CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,W(LGQ),W(LQ),W(LWRK))
      ELSE IF (LP) THEN
         CALL DCOPY(N,W(LCQ),1,W(LGQ),1)
      END IF
C
      IF (NUMINF.EQ.0 .AND. LP) THEN
         OBJ = DDOT(N,CVEC,1,X,1)
      ELSE
         OBJ = SUMINF
      END IF
C
      MSG = EMPTY
C
C*    ======================Start of main loop==========================
C     +    do while (msg .eq. empty)
   20 IF (MSG.EQ.EMPTY) THEN
C
         GZNORM = ZERO
         IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
C
         IF (NRZ.EQ.NZ) THEN
            GRZNRM = GZNORM
         ELSE
            GRZNRM = ZERO
            IF (NRZ.GT.0) GRZNRM = DNRM2(NRZ,W(LGQ),1)
         END IF
C
         GFNORM = GZNORM
         IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),
     *       1)
C
C        ---------------------------------------------------------------
C        Print the details of this iteration.
C        ---------------------------------------------------------------
C        Define small quantities that reflect the size of x, R and
C        the constraints in the working set.
C
         IF (PRNT) THEN
            CONDT = ONE
            IF (NACTIV.GT.0) CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
C
            CALL E04MFU(PRBTYP,HEADER,RSET,MSGLVL,ITER,ISDEL,JDEL,JADD,
     *                  N,NCLIN,NACTIV,NFREE,NZ,NRZ,LDR,LDT,ISTATE,ALFA,
     *                  CONDRZ,CONDT,DZZ,GFNORM,GZNORM,NUMINF,SUMINF,
     *                  NOTOPT,OBJ,TRULAM,AX,W(LR),W(LT),X,W(LWRK))
            JDEL = 0
            JADD = 0
            ALFA = ZERO
         END IF
C
         IF (NUMINF.GT.0) THEN
            DINKY = EPSPT8*ABS(SUMINF)
         ELSE
            OBJSIZ = ONE + ABS(OBJ)
            WSSIZE = ZERO
            IF (NACTIV.GT.0) WSSIZE = DTMAX
            DINKY = EPSPT8*MAX(WSSIZE,OBJSIZ,GFNORM)
         END IF
C
         IF (LCDBG .AND. ILCDBG(1).GT.0) THEN
            WRITE (REC,FMT=99998) GRZNRM, DINKY
            CALL X04BAY(IPRINT,3,REC)
         END IF
C
C        If the reduced gradient Z'g is small enough,
C        Lagrange multipliers will be computed.
C
         IF (NUMINF.EQ.0 .AND. FP) THEN
            MSG = 'feasbl'
            NFIXED = N - NFREE
            CALL F06FBF(NACTIV+NFIXED,ZERO,W(LRLAM),1)
            GO TO 20
         END IF
C
         IF (GRZNRM.LE.DINKY) THEN
C           =========================================================
C           The point  x  is a constrained stationary point.
C           Compute Lagrange multipliers.
C           =========================================================
C           Define what we mean by 'tiny' and non-optimal multipliers.
C
            NOTOPT = 0
            JDEL = 0
            ZEROLM = -DINKY
            SMLLST = -DINKY
            BIGGST = DINKY + ONE
            TINYST = DINKY
C
            CALL E04MFM(PRBTYP,MSGLVL,N,LDA,LDT,NACTIV,NFREE,NZ,ISTATE,
     *                  KACTIV,KX,ZEROLM,NOTOPT,NUMINF,TRUSML,SMLLST,
     *                  JSMLST,KSMLST,TINYST,JTINY,JINF,TRUBIG,BIGGST,
     *                  JBIGST,KBIGST,A,W(LANORM),W(LGQ),W(LRLAM),W(LT),
     *                  W(LWTINF))
C
            IF (NRZ.LT.NZ) CALL E04MFL(MSGLVL,N,NRZ,NZ,ZEROLM,NOTOPT,
     *                                 NUMINF,TRUSML,SMLLST,JSMLST,
     *                                 TINYST,JTINY,W(LGQ))
C
            IF (ABS(JSMLST).GT.0) THEN
C              ------------------------------------------------------
C              Delete a constraint.
C              ------------------------------------------------------
C              E04MFM  or  E04MFL  found a non-optimal multiplier.
C
               TRULAM = TRUSML
               JDEL = JSMLST
C
               IF (JSMLST.GT.0) THEN
C
C                 Regular constraint.
C
                  KDEL = KSMLST
                  ISDEL = ISTATE(JDEL)
                  ISTATE(JDEL) = 0
               END IF
            ELSE
               IF (NUMINF.GT.0 .AND. JBIGST.GT.0) THEN
C
C                 No feasible point exists for the constraints but
C                 the sum of the constraint violations can be reduced
C                 by moving off constraints with multipliers greater
C                 than 1.
C
                  JDEL = JBIGST
                  KDEL = KBIGST
                  ISDEL = ISTATE(JDEL)
                  IF (TRUBIG.LE.ZERO) IS = -1
                  IF (TRUBIG.GT.ZERO) IS = -2
                  ISTATE(JDEL) = IS
                  TRULAM = TRUBIG
                  FIRSTV = .TRUE.
                  NUMINF = NUMINF + 1
               END IF
            END IF
C
            IF (JDEL.EQ.0) THEN
               IF (NUMINF.GT.0) THEN
                  MSG = 'infeas'
               ELSE
                  MSG = 'optiml'
               END IF
               GO TO 20
            END IF
C
C           Constraint  JDEL  has been deleted.
C           Update the  TQ  factorization.
C
            CALL E04NFP(UNITQ,IT,N,NACTIV,NFREE,NGQ,NZ,NRZ,LDA,LDQ,LDT,
     *                  JDEL,KDEL,KACTIV,KX,A,W(LT),W(LGQ),W(LQ),W(LWRK)
     *                  ,W(LD),W(LRLAM))
            IF (RSET) CALL E04MFY(NRZ,LDR,W(LR),ONE)
C
            PRNT = .FALSE.
         ELSE
C           ============================================================
C           Compute a search direction.
C           ============================================================
            IF (ITER.GE.ITMAX) THEN
               MSG = 'itnlim'
               GO TO 20
            END IF
C
            PRNT = .TRUE.
            ITER = ITER + 1
C
            IF (ITER.EQ.IDBG) THEN
               LCDBG = .TRUE.
               CMDBG = LCDBG
               MSGLVL = MSGSVD
            END IF
C
            CALL DCOPY(NRZ,W(LGQ),1,W(LD),1)
            CALL DSCAL(NRZ,(-ONE),W(LD),1)
C
            DNORM = DNRM2(NRZ,W(LD),1)
C
            CALL E04NBW(1,N,NRZ,NFREE,LDQ,UNITQ,KX,W(LD),W(LQ),W(LWRK))
            CALL DGEMV('No transpose',NCLIN,N,ONE,A,LDA,W(LD),1,ZERO,
     *                 W(LAD),1)
C
            IF (LCDBG) THEN
               IF (ILCDBG(2).GT.0) THEN
                  CALL E04MFN('E04MFZ','     d',W(LD),N)
                  CALL E04MFN('E04MFZ','    Ad',W(LAD),NCLIN)
               END IF
            END IF
C
C           ---------------------------------------------------------
C           Find the constraint we bump into along d.
C           Update  x  and  Ax  if the step alfa is nonzero.
C           ---------------------------------------------------------
C           ALFHIT is initialized to BIGALF. If it remains that value
C           after the call to  E04MFS, it is regarded as infinite.
C
            BIGALF = F06BLF(BIGDX,DNORM,OVERFL)
C
            CALL E04MFS(FIRSTV,N,NCLIN,ISTATE,BIGALF,BIGBND,DNORM,
     *                  HITLOW,MOVE,ONBND,UNBNDD,ALFHIT,ALFAP,JADD,
     *                  W(LANORM),W(LAD),AX,BL,BU,FEATOL,FEATLU,W(LD),X)
C
            IF (UNBNDD) THEN
               MSG = 'unbndd'
               GO TO 20
            END IF
C
            ALFA = ALFHIT
            CALL DAXPY(N,ALFA,W(LD),1,X,1)
C
            IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,W(LAD),1,AX,1)
            XNORM = DNRM2(N,X,1)
C
C           ---------------------------------------------------------
C           Add a constraint to the working set.
C           Update the  TQ  factors of the working set.
C           Use  d  as temporary work space.
C           ---------------------------------------------------------
            IF (BL(JADD).EQ.BU(JADD)) THEN
               ISTATE(JADD) = 3
            ELSE IF (HITLOW) THEN
               ISTATE(JADD) = 1
            ELSE
               ISTATE(JADD) = 2
            END IF
C
            IF (JADD.GT.N) THEN
               IADD = JADD - N
            ELSE
               IF (ALFA.GE.ZERO) THEN
                  IF (HITLOW) THEN
                     X(JADD) = BL(JADD)
                  ELSE
                     X(JADD) = BU(JADD)
                  END IF
               END IF
               DO 40 IFIX = 1, NFREE
                  IF (KX(IFIX).EQ.JADD) GO TO 60
   40          CONTINUE
   60       END IF
C
            CALL E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,NACTIV,NZ,
     *                  NFREE,NRZ,NGQ,N,LDA,LDQ,LDR,LDT,KX,CONDMX,DZZ,A,
     *                  W(LR),W(LT),W(LGQ),W(LQ),W(LWRK),W(LRLAM),W(LD),
     *                  MSGLVL)
C
            NZ = NZ - 1
            NRZ = NRZ - 1
C
            IF (JADD.LE.N) THEN
C
C              A simple bound has been added.
C
               NFREE = NFREE - 1
            ELSE
C
C              A general constraint has been added.
C
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = IADD
            END IF
C
C           Increment FEATOL.
C
            CALL DAXPY(NCTOTL,TOLINC,FEATLU,1,FEATOL,1)
C
            IF (MOD(ITER,KCHK).EQ.0) THEN
C              ------------------------------------------------------
C              Check the feasibility of constraints with non-
C              negative ISTATE values.  If some violations have
C              occurred.  Set INFORM to force iterative
C              refinement and a switch to phase 1.
C              ------------------------------------------------------
               CALL E04MFQ(N,NCLIN,ISTATE,BIGBND,NVIOL,JMAX,ERRMAX,AX,
     *                     BL,BU,FEATOL,X)
C
               IF (NVIOL.GT.0) THEN
                  IF (MSGLVL.GT.0) THEN
                     WRITE (REC,FMT=99999) ERRMAX, JMAX
                     CALL X04BAF(IPRINT,REC(1))
                  END IF
               END IF
            END IF
C
            IF (MOD(ITER,KDEGEN).EQ.0) THEN
C
C              Every  KDEGEN  iterations, reset  FEATOL  and
C              move  x  on to the working set if it is close.
C
               CALL E04MFR('End of cycle',MSGLVL,N,NCLIN,NMOVED,ITER,
     *                     NUMINF,ISTATE,BIGBND,AX,BL,BU,FEATOL,FEATLU,
     *                     X)
C
               NVIOL = NVIOL + NMOVED
            END IF
C
            IF (NVIOL.GT.0) THEN
               MSG = 'resetx'
               GO TO 20
            END IF
C
            IF (NUMINF.NE.0) THEN
               CALL E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,
     *                     BU,A,FEATOL,W(LGQ),X,W(LWTINF))
C
               IF (NUMINF.GT.0) THEN
                  CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,W(LGQ),W(LQ),
     *                        W(LWRK))
               ELSE IF (LP) THEN
                  CALL DCOPY(N,W(LCQ),1,W(LGQ),1)
               END IF
            END IF
C
            IF (NUMINF.EQ.0 .AND. LP) THEN
               OBJ = DDOT(N,CVEC,1,X,1)
            ELSE
               OBJ = SUMINF
            END IF
         END IF
         GO TO 20
C        +    end while
      END IF
C     ======================end of main loop============================
C
      IF (MSG.EQ.'optiml') THEN
         IF (LP) THEN
            IF (NRZ.LT.NZ) THEN
               MSG = 'weak  '
            ELSE
               NTFIXD = 0
               DO 80 J = 1, N
                  IF (ISTATE(J).EQ.4) NTFIXD = NTFIXD + 1
   80          CONTINUE
               IF (NTFIXD.GT.0) MSG = 'weak  '
            END IF
            IF (ABS(JTINY).GT.0) MSG = 'weak  '
         END IF
      ELSE IF (MSG.EQ.'unbndd' .AND. NUMINF.GT.0) THEN
         MSG = 'infeas'
      END IF
C
      MSGLVL = MSGSVD
C
      RETURN
C
C
C     End of  E04MFZ.  (LPCORE)
C
99999 FORMAT (' XXX  Iterative refinement.  The max violation is ',1P,
     *       D10.2,' in constraint',I5)
99998 FORMAT (/' //E04MFZ//       GRZNRM      DINKY',/' //E04MFZ//  ',
     *       1P,2D11.2)
      END
      SUBROUTINE E04UCP(N,NCLIN,NCNLN,NCTOTL,LITOTL,LWTOTL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-721 (DEC 1989).
C     MARK 16 REVISED. IER-1084 (JUL 1993).
C
C     ******************************************************************
C     E04UCP   allocates the addresses of the work arrays for E04UCZ and
C     E04NCZ.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   14-February-1985.
C     This version of  E04UCP  dated 12-July-1986.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LWTOTL, N, NCLIN, NCNLN, NCTOTL
C     .. Scalars in Common ..
      INTEGER           LDT, LDZY, LENNAM, NCOLT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG), LOCLS(LENLS), LOCNP(LENNP)
C     .. Local Scalars ..
      INTEGER           LADX, LANORM, LAQP, LBL, LBU, LC1MUL, LCJAC,
     *                  LCJDX, LCMUL, LCS1, LCS2, LDLAM, LDSLK, LDX,
     *                  LENAQP, LENT, LENZY, LFEATL, LGQ, LGQ1, LGRAD,
     *                  LHCTRL, LHFRWD, LIPERM, LKACTV, LKX, LNEEDC,
     *                  LQPADX, LQPDX, LQPGQ, LQPHZ, LQPTOL, LRHO,
     *                  LRLAM, LRPQ, LRPQ0, LSLK, LSLK1, LT, LWRK1,
     *                  LWRK2, LWRK3, LWTINF, LX1, LZY, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
C     Assign array lengths that depend upon the problem dimensions.
C
      IF (NCLIN+NCNLN.EQ.0) THEN
         LENT = 0
         LENZY = 0
      ELSE
         LENT = LDT*NCOLT
         LENZY = LDZY*LDZY
      END IF
C
      IF (NCNLN.EQ.0) THEN
         LENAQP = 0
      ELSE
         LENAQP = (NCLIN+NCNLN)*N
      END IF
C
      LKACTV = MINIW
      LKX = LKACTV + N
      LNEEDC = LKX + N
      LIPERM = LNEEDC + NCNLN
      MINIW = LIPERM + NCTOTL
C
      LHFRWD = MINW
      LHCTRL = LHFRWD + N
      LANORM = LHCTRL + N
      LQPGQ = LANORM + NCLIN + NCNLN
      LGQ = LQPGQ + N
      LRLAM = LGQ + N
      LT = LRLAM + N
      LZY = LT + LENT
      MINW = LZY + LENZY
C
      LOCLS(1) = LKACTV
      LOCLS(2) = LANORM
      LOCLS(8) = LQPGQ
      LOCLS(9) = LGQ
      LOCLS(10) = LRLAM
      LOCLS(11) = LT
      LOCLS(12) = LZY
C
C     Assign the addresses for the workspace arrays used by  E04UCU.
C
      LQPADX = MINW
      LQPDX = LQPADX + NCLIN + NCNLN
      LRPQ = LQPDX + N
      LRPQ0 = LRPQ + N
      LQPHZ = LRPQ0 + N
      LWTINF = LQPHZ + N
      LWRK1 = LWTINF + NCTOTL
      LQPTOL = LWRK1 + NCTOTL
      MINW = LQPTOL + NCTOTL
C
      LOCLS(3) = LQPADX
      LOCLS(4) = LQPDX
      LOCLS(5) = LRPQ
      LOCLS(6) = LRPQ0
      LOCLS(7) = LQPHZ
      LOCLS(13) = LWTINF
      LOCLS(14) = LWRK1
      LOCLS(15) = LQPTOL
C
C     Assign the addresses for arrays used in E04UCZ.
C
      LAQP = MINW
      LADX = LAQP + LENAQP
      LBL = LADX + NCLIN + NCNLN
      LBU = LBL + NCTOTL
      LDX = LBU + NCTOTL
      LGQ1 = LDX + N
      LFEATL = LGQ1 + N
      LX1 = LFEATL + NCTOTL
      LWRK2 = LX1 + N
      MINW = LWRK2 + NCTOTL
C
      LOCNP(1) = LKX
      LOCNP(2) = LIPERM
      LOCNP(3) = LAQP
      LOCNP(4) = LADX
      LOCNP(5) = LBL
      LOCNP(6) = LBU
      LOCNP(7) = LDX
      LOCNP(8) = LGQ1
      LOCNP(10) = LFEATL
      LOCNP(11) = LX1
      LOCNP(12) = LWRK2
C
      LCS1 = MINW
      LCS2 = LCS1 + NCNLN
      LC1MUL = LCS2 + NCNLN
      LCMUL = LC1MUL + NCNLN
      LCJDX = LCMUL + NCNLN
      LDLAM = LCJDX + NCNLN
      LDSLK = LDLAM + NCNLN
      LRHO = LDSLK + NCNLN
      LWRK3 = LRHO + NCNLN
      LSLK1 = LWRK3 + NCNLN
      LSLK = LSLK1 + NCNLN
      MINW = LSLK + NCNLN
C
      LOCNP(13) = LCS1
      LOCNP(14) = LCS2
      LOCNP(15) = LC1MUL
      LOCNP(16) = LCMUL
      LOCNP(17) = LCJDX
      LOCNP(18) = LDLAM
      LOCNP(19) = LDSLK
      LOCNP(20) = LRHO
      LOCNP(21) = LWRK3
      LOCNP(22) = LSLK1
      LOCNP(23) = LSLK
      LOCNP(24) = LNEEDC
C
      LCJAC = MINW
      LGRAD = LCJAC + NCNLN*N
      MINW = LGRAD + N
C
      LOCNP(25) = LHFRWD
      LOCNP(26) = LHCTRL
      LOCNP(27) = LCJAC
      LOCNP(28) = LGRAD
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04UCP. (NPLOC)
C
      END
      SUBROUTINE E04UCZ(NAMED,NAMES,UNITQ,INFORM,MAJITS,N,NCLIN,NCNLN,
     *                  NCTOTL,NACTIV,NFREE,NZ,LDCJ,LDCJU,LDAQP,LDR,
     *                  NFUN,NGRAD,ISTATE,KACTIV,KX,OBJF,FDNORM,XNORM,
     *                  CONFUN,OBJFUN,AQP,AX,BL,BU,C,CJAC,CJACU,CLAMDA,
     *                  FEATOL,GRAD,GRADU,R,X,IW,W,LENW,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1093 (JUL 1993).
C
C     ******************************************************************
C     E04UCZ  is the core routine for  E04UCF,  a sequential quadratic
C     programming (SQP) method for nonlinearly constrained optimization.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      February-1982.
C     This version of E04UCZ dated 23-Dec-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  GROWTH
      PARAMETER         (GROWTH=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FDNORM, OBJF, XNORM
      INTEGER           INFORM, LDAQP, LDCJ, LDCJU, LDR, LENW, MAJITS,
     *                  N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE, NFUN,
     *                  NGRAD, NZ
      LOGICAL           NAMED, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  AQP(LDAQP,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  C(*), CJAC(LDCJ,*), CJACU(LDCJU,*),
     *                  CLAMDA(NCTOTL), FEATOL(NCTOTL), GRAD(N),
     *                  GRADU(N), R(LDR,*), USER(*), W(LENW), X(N)
      INTEGER           ISTATE(*), IUSER(*), IW(*), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DRMAX, DRMIN, DTMAX, DTMIN, DXLIM, EPSPT3,
     *                  EPSPT5, EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL,
     *                  HCNDBD, RCNDBD, RFROBN, RHODMP, RHOMAX, RHONRM,
     *                  SCALE, TOLACT, TOLFEA, TOLRNK
      INTEGER           IDBGLS, IDBGNP, IPRINT, IPRNT, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, KSAVE, LCRASH, LDBGLS, LDBGNP, LDT,
     *                  LDZY, LENNAM, LFDSET, LFORMH, LINES1, LINES2,
     *                  LPROB, LVERFY, LVLDER, LVLDIF, LVRFYC, MSGLS,
     *                  MSGNP, NCDIFF, NCOLT, NFDIFF, NLNF, NLNJ, NLNX,
     *                  NLOAD, NN, NNCLIN, NNCNLN, NOUT, NPROB, NSAVE
      LOGICAL           CMDBG, INCRUN, NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), INPDBG(LDBG), IPADLS(18),
     *                  IPADNP(12), IPSVLS(MXPARM), IPSVNP(MXPARM),
     *                  JVERFY(4), LOCLS(LENLS), LOCNP(LENNP)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFA, ALFBND, ALFDX, ALFLIM, ALFMAX, ALFMIN,
     *                  ALFSML, CNORM, COND, CONDH, CONDHZ, CONDT,
     *                  CVNORM, DINKY, DRZMAX, DRZMIN, DXNORM, ERRMAX,
     *                  FLMAX, GDX, GFNORM, GLF1, GLF2, GLNORM, GLTEST,
     *                  GRDALF, GTEST, GZNORM, OBJ, OBJALF, OBJSIZ,
     *                  QPCURV, ROOTN, RTFTOL, RTMAX, XSIZE
      INTEGER           IDBG, INFO, JMAX, LADX, LANORM, LAQP, LBL, LBU,
     *                  LC1MUL, LCJAC1, LCJDX, LCJDX1, LCMUL, LCS1,
     *                  LCS2, LDCJ1, LDLAM, LDSLK, LDX, LGQ, LGQ1,
     *                  LHCTRL, LHFRWD, LHPQ, LINACT, LIPERM, LNEEDC,
     *                  LQPTOL, LQRWRK, LRHO, LRLAM, LRPQ, LSLK, LSLK1,
     *                  LT, LVIOLN, LWRK1, LWRK2, LWRK3, LWTINF, LX1,
     *                  LZY, MAJIT0, MINITS, MJRDBG, MNR, MNRDBG,
     *                  MNRSUM, MODE, MSGQP, MSGSV1, MSGSV2, NCQP, NL,
     *                  NLNACT, NLSERR, NMAJOR, NMINOR, NPLIN, NQPERR,
     *                  NQPINF, NSTATE, NUMINF, NVIOL
      LOGICAL           CENTRL, CONVPT, CONVRG, DONE, ERROR, FEASQP,
     *                  GOODGQ, INFEAS, NEEDFD, NEWGQ, OPTIML, OVERFL
      CHARACTER*5       MJRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      LOGICAL           KTCOND(2)
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04NBW, E04NBX, E04UCL,
     *                  E04UCN, E04UCR, E04UCT, E04UCU, E04UCW, E04UDR,
     *                  E04UDS, E04UDT, F06DBF, F06FLF, F06QFF, X04BAF,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IDBGLS, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LPROB, MSGLS, NN,
     *                  NNCLIN, NPROB, IPADLS
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /FE04UC/INPDBG, NPDBG
      COMMON            /GE04UC/IPSVNP, IDBGNP, ITMXNP, JVRFY1, JVRFY2,
     *                  JVRFY3, JVRFY4, LDBGNP, LFORMH, LVLDER, LVERFY,
     *                  MSGNP, NLNF, NLNJ, NLNX, NNCNLN, NSAVE, NLOAD,
     *                  KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IDBGLS), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),IDBGNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IDBGNP,IDBG), (ITMXNP,NMAJOR), (ITMAX2,NMINOR)
      EQUIVALENCE       (LDBGLS,MNRDBG), (LDBGNP,MJRDBG), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /AX02ZA/, /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/
C     .. Executable Statements ..
C
C     Specify machine-dependent parameters.
C
      FLMAX = WMACH(7)
      RTMAX = WMACH(8)
C
      LANORM = LOCLS(2)
      LRPQ = LOCLS(5)
      LQRWRK = LOCLS(6)
      LHPQ = LOCLS(8)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1 = LOCLS(14)
      LQPTOL = LOCLS(15)
C
      LIPERM = LOCNP(2)
      LAQP = LOCNP(3)
      LADX = LOCNP(4)
      LBL = LOCNP(5)
      LBU = LOCNP(6)
      LDX = LOCNP(7)
      LGQ1 = LOCNP(8)
      LX1 = LOCNP(11)
      LWRK2 = LOCNP(12)
      LCS1 = LOCNP(13)
      LCS2 = LOCNP(14)
      LC1MUL = LOCNP(15)
      LCMUL = LOCNP(16)
      LCJDX1 = LOCNP(17)
      LDLAM = LOCNP(18)
      LDSLK = LOCNP(19)
      LRHO = LOCNP(20)
      LWRK3 = LOCNP(21)
      LSLK1 = LOCNP(22)
      LSLK = LOCNP(23)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
C
      LCJAC1 = LAQP + NCLIN
      LCJDX = LADX + NCLIN
      LVIOLN = LWRK3
C
C     Initialize
C
      MJRMSG = '     '
      NQPINF = 0
      MNRSUM = 0
C
      MAJIT0 = MAJITS
      NPLIN = N + NCLIN
      NCQP = NCLIN + NCNLN
      NL = MIN(NPLIN+1,NCTOTL)
C
      LDCJ1 = MAX(NCQP,1)
C
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR.
     *         (LVLDER.EQ.1 .AND. NCNLN.GT.0)
C
      ALFA = ZERO
      ALFDX = ZERO
      RTFTOL = SQRT(FTOL)
      ROOTN = SQRT(DBLE(N))
C
C     If debug printing is required,  turn off any extensive printing
C     until iteration  IDBG.
C
      MSGSV1 = MSGNP
      MSGSV2 = MSGQP
      IF (IDBG.LE.NMAJOR .AND. IDBG.GT.0) THEN
         MSGNP = 0
         IF (MSGSV1.GE.5) MSGNP = 5
         MSGQP = 0
         IF (MSGSV2.GE.5) MSGQP = 5
      END IF
C
C     ------------------------------------------------------------------
C     Information from the feasibility phase will be used to generate a
C     hot start for the first QP subproblem.
C     ------------------------------------------------------------------
      CALL DCOPY(NCTOTL,FEATOL,1,W(LQPTOL),1)
C
      NSTATE = 0
C
      OBJALF = OBJF
      IF (NCNLN.GT.0) THEN
         OBJALF = OBJALF - DDOT(NCNLN,W(LCMUL),1,C,1)
      END IF
C
      NEWGQ = .FALSE.
C
C*    ==================================================================
C+    repeat                             (until converged or error exit)
C
C     ===============================================================
C     See if we want to save the details of this iteration.
C     ===============================================================
C     20 IF (MOD(MAJITS,KSAVE).EQ.0 .AND. MAJITS.NE.MAJIT0) THEN
C        CALL NPSAVR(UNITQ,N,NCLIN,NCNLN,LDR,LDQ,NFREE,NSAVE,MAJITS,
C    *               ISTATE,KX,W(LHFRWD),W(LHCTRL),W(LCMUL),R,W(LRHO),X,
C    *               X,W(LQ))
C     END IF
C
C   *    ===============================================================
C   +    repeat                         (Until a good gradient is found)
C
   20 MINITS = 0
C
   40 CENTRL = LVLDIF .EQ. 2
C
      IF (NEWGQ) THEN
         IF (NEEDFD) THEN
C           ------------------------------------------------------
C           Compute any missing gradient elements and the
C           transformed gradient of the objective.
C           ------------------------------------------------------
            CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,IW(LNEEDC),BL,
     *                  BU,C,W(LWRK2),W(LWRK3),CJAC,CJACU,GRAD,GRADU,
     *                  W(LHFRWD),W(LHCTRL),X,W,LENW,IUSER,USER)
            INFORM = MODE
            IF (MODE.LT.0) GO TO 60
C
         END IF
C
         CALL DCOPY(N,GRAD,1,W(LGQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,W(LGQ),W(LZY),W(LWRK1))
         NEWGQ = .FALSE.
      END IF
C
C     ============================================================
C     (1) Solve an inequality quadratic program (IQP) for the
C         search direction and multiplier estimates.
C     (2) For each nonlinear inequality constraint,  compute
C         the slack variable for which the merit function is
C         minimized.
C     (3) Compute the search direction for the slack variables
C         and multipliers.
C
C     Note that the array VIOLN is WRK3.
C     ============================================================
      CALL E04UCU(FEASQP,UNITQ,NQPERR,MAJITS,MNR,N,NCLIN,NCNLN,LDCJ,
     *            LDAQP,LDR,LINACT,NLNACT,NACTIV,NFREE,NZ,NUMINF,ISTATE,
     *            KACTIV,KX,DXNORM,GDX,QPCURV,AQP,W(LADX),W(LANORM),AX,
     *            BL,BU,C,CJAC,CLAMDA,W(LCMUL),W(LCS1),W(LDLAM),W(LDSLK)
     *            ,W(LDX),W(LBL),W(LBU),W(LQPTOL),R,W(LRHO),W(LSLK),
     *            W(LVIOLN),X,W(LWTINF),W)
C
      MINITS = MINITS + MNR
      MNRSUM = MNRSUM + MNR
C
      IF (FEASQP) THEN
         NQPINF = 0
      ELSE
         NQPINF = NQPINF + 1
         MJRMSG(2:2) = 'Infeasible subproblem'
      END IF
C
C     ============================================================
C     Compute quantities needed for the convergence test.
C     ============================================================
C     Compute the norms of the projected gradient and the
C     gradient with respect to the free variables.
C
      GZNORM = ZERO
      IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
      GFNORM = GZNORM
      IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),1)
C
C     If the forward-difference estimate of the transformed
C     gradient of the Lagrangian function is small,  switch to
C     central differences, recompute the derivatives and re-solve
C     the QP.
C
      GOODGQ = .TRUE.
      IF (NEEDFD .AND. .NOT. CENTRL) THEN
         GLNORM = DNRM2(N,W(LHPQ),1)
         IF (NCNLN.EQ.0) THEN
            CNORM = ZERO
         ELSE
            CNORM = DNRM2(NCNLN,C,1)
         END IF
C
         GLTEST = (ONE+ABS(OBJF)+ABS(CNORM))*EPSRF/FDNORM
         IF (GLNORM.LE.GLTEST) THEN
            GOODGQ = .FALSE.
            MJRMSG(3:3) = 'Central differences'
            LVLDIF = 2
            NEWGQ = .TRUE.
            IF (MSGNP.GE.5) THEN
               IF (MINITS.GT.0) THEN
                  WRITE (REC,FMT=99999) MINITS
                  CALL X04BAF(IPRINT,REC(1))
               END IF
            END IF
         END IF
C
      END IF
C
C     +       UNTIL     (GOODGQ)
      IF ( .NOT. GOODGQ) GO TO 40
C
C     ===============================================================
C     (1) Compute the number of constraints that are violated by more
C         than FEATOL.
C     (2) Compute the 2-norm of the residuals of the constraints in
C         the QP working set.
C     ===============================================================
      CALL E04UCW(N,NCLIN,NCNLN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,NVIOL,
     *            AX,BL,BU,C,FEATOL,X,W(LWRK2))
C
C     Define small quantities that reflect the magnitude of OBJF and
C     the norm of GRAD(free).
C
      OBJSIZ = ONE + ABS(OBJF)
      XSIZE = ONE + XNORM
      GTEST = MAX(OBJSIZ,GFNORM)
      DINKY = RTFTOL*GTEST
C
      IF (NACTIV.EQ.0) THEN
         CONDT = ZERO
      ELSE IF (NACTIV.EQ.1) THEN
         CONDT = DTMIN
      ELSE
         CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
      END IF
C
      CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
C
      CONDH = F06BLF(DRMAX,DRMIN,OVERFL)
      IF (CONDH.LT.RTMAX) THEN
         CONDH = CONDH*CONDH
      ELSE
         CONDH = FLMAX
      END IF
C
      IF (NZ.EQ.0) THEN
         CONDHZ = ONE
      ELSE IF (NZ.EQ.N) THEN
         CONDHZ = CONDH
      ELSE
         CALL F06FLF(NZ,R,LDR+1,DRZMAX,DRZMIN)
         CONDHZ = F06BLF(DRZMAX,DRZMIN,OVERFL)
         IF (CONDHZ.LT.RTMAX) THEN
            CONDHZ = CONDHZ*CONDHZ
         ELSE
            CONDHZ = FLMAX
         END IF
      END IF
C
C     ---------------------------------------------------------------
C     Test for convergence.
C     The point test CONVPT checks for a K-T point at the initial
C     point or after a large change in X.
C     ---------------------------------------------------------------
      CONVPT = DXNORM .LE. EPSPT8*GTEST .AND. NVIOL .EQ. 0 .AND.
     *         NQPERR .LE. 1
C
      KTCOND(1) = GZNORM .LT. DINKY
      KTCOND(2) = NVIOL .EQ. 0
      OPTIML = KTCOND(1) .AND. KTCOND(2)
C
      CONVRG = MAJITS .GT. 0 .AND. ALFDX .LE. RTFTOL*XSIZE
C
      INFEAS = CONVRG .AND. .NOT. FEASQP .OR. NQPINF .GT. 7
C
      DONE = CONVPT .OR. (CONVRG .AND. OPTIML) .OR. INFEAS
C
      OBJALF = OBJF
      GRDALF = GDX
      GLF1 = GDX
      IF (NCNLN.GT.0) THEN
         GLF1 = GLF1 - DDOT(NCNLN,W(LCJDX),1,CLAMDA(NL),1)
C
C        Compute the value and directional derivative of the
C        augmented Lagrangian merit function.
C        The penalty parameters may be increased or decreased.
C
         CALL E04UCN(FEASQP,N,NCLIN,NCNLN,OBJALF,GRDALF,QPCURV,ISTATE,
     *               W(LCJDX),W(LCMUL),W(LCS1),W(LDLAM),W(LRHO),
     *               W(LVIOLN),W(LWRK1),W(LWRK2))
      END IF
C
C     ===============================================================
C     Print the details of this iteration.
C     ===============================================================
      CALL E04UCT(KTCOND,CONVRG,MJRMSG,MSGNP,MSGQP,LDR,LDT,N,NCLIN,
     *            NCNLN,NCTOTL,NACTIV,LINACT,NLNACT,NZ,NFREE,MAJIT0,
     *            MAJITS,MINITS,ISTATE,ALFA,NFUN,CONDHZ,CONDH,CONDT,
     *            OBJALF,OBJF,GFNORM,GZNORM,CVNORM,AX,C,R,W(LT),
     *            W(LVIOLN),X,W(LWRK1))
C
      ALFA = ZERO
      ERROR = MAJITS .GE. NMAJOR
C
      IF ( .NOT. (DONE .OR. ERROR)) THEN
         MAJITS = MAJITS + 1
C
         IF (MAJITS.EQ.IDBG) THEN
            NPDBG = .TRUE.
            CMDBG = NPDBG
            MSGNP = MSGSV1
            MSGQP = MSGSV2
         END IF
C
C        Make copies of information needed for the BFGS update.
C
         CALL DCOPY(N,X,1,W(LX1),1)
         CALL DCOPY(N,W(LGQ),1,W(LGQ1),1)
C
         IF (NCNLN.GT.0) THEN
            CALL DCOPY(NCNLN,W(LCJDX),1,W(LCJDX1),1)
            CALL DCOPY(NCNLN,W(LCMUL),1,W(LC1MUL),1)
            CALL DCOPY(NCNLN,W(LSLK),1,W(LSLK1),1)
         END IF
C
C        ============================================================
C        Compute the parameters for the linesearch.
C        ============================================================
C        ALFMIN is the smallest allowable step predicted by the QP
C        subproblem.
C
         ALFMIN = ONE
         IF ( .NOT. FEASQP) ALFMIN = ZERO
C
C        ------------------------------------------------------------
C        ALFMAX is the largest feasible steplength subject to a user-
C        defined limit ALFLIM on the change in X.
C        ------------------------------------------------------------
         IF (NCNLN.GT.0 .AND. NEEDFD) THEN
            ALFMAX = ONE
         ELSE
            ALFMAX = F06BLF(BIGDX,DXNORM,OVERFL)
            CALL E04UDT(INFO,N,NCLIN,NCNLN,ALFA,ALFMIN,ALFMAX,BIGBND,
     *                  DXNORM,W(LANORM),W(LADX),AX,BL,BU,W(LDSLK),
     *                  W(LDX),W(LSLK),X)
            ALFMAX = ALFA
            IF (ALFMAX.LT.ONE+EPSPT3 .AND. FEASQP) ALFMAX = ONE
         END IF
C
C        ------------------------------------------------------------
C        ALFBND is a tentative upper bound on the steplength.  If the
C        merit function is decreasing at ALFBND and certain
C        conditions hold,  ALFBND will be increased in multiples of
C        two (subject to not being greater than ALFMAX).
C        ------------------------------------------------------------
         IF (NCNLN.EQ.0) THEN
            ALFBND = ALFMAX
         ELSE
            ALFBND = MIN(ONE,ALFMAX)
         END IF
C
C        ------------------------------------------------------------
C        ALFSML trips the computation of central differences.  If a
C        trial steplength falls below ALFSML, the linesearch is
C        terminated.
C        ------------------------------------------------------------
         ALFSML = ZERO
         IF (NEEDFD .AND. .NOT. CENTRL) THEN
            ALFSML = F06BLF(FDNORM,DXNORM,OVERFL)
            ALFSML = MIN(ALFSML,ALFMAX)
         END IF
C
C        ============================================================
C        Compute the steplength using safeguarded interpolation.
C        ============================================================
         ALFLIM = F06BLF((ONE+XNORM)*DXLIM,DXNORM,OVERFL)
         ALFA = MIN(ALFLIM,ONE)
C
         CALL E04UCR(NEEDFD,NLSERR,N,NCNLN,LDCJ,LDCJU,NFUN,NGRAD,
     *               IW(LNEEDC),CONFUN,OBJFUN,ALFA,ALFBND,ALFMAX,ALFSML,
     *               DXNORM,EPSRF,ETA,GDX,GRDALF,GLF1,GLF2,OBJF,OBJALF,
     *               QPCURV,XNORM,C,W(LWRK1),CJAC,CJACU,W(LCJDX),
     *               W(LWRK3),W(LC1MUL),W(LCMUL),W(LCS1),W(LCS2),W(LDX),
     *               W(LDLAM),W(LDSLK),GRAD,GRADU,CLAMDA(NL),W(LRHO),
     *               W(LSLK1),W(LSLK),W(LX1),X,W(LWRK2),W,LENW,IUSER,
     *               USER)
C
C           ------------------------------------------------------------
C           E04UCR  sets NLSERR to the following values...
C
C           < 0  if the user wants to stop.
C             1  if the search is successful and ALFA < ALFMAX.
C             2  if the search is successful and ALFA = ALFMAX.
C             3  if a better point was found but too many functions
C                were needed (not sufficient decrease).
C
C           Values of NLSERR occurring with a nonzero value of ALFA.
C             4  if ALFMAX < TOLABS (too small to do a search).
C             5  if ALFA  < ALFSML (E04UCJ only -- maybe want to switch
C                to central differences to get a better direction).
C             6  if the search found that there is no useful step.
C                The interval of uncertainty is less than 2*TOLABS.
C                The minimizer is very close to ALFA = zero
C                or the gradients are not sufficiently accurate.
C             7  if there were too many function calls.
C             8  if the input parameters were bad
C                (ALFMAX le TOLTNY  or  uphill).
C           ------------------------------------------------------------
         IF (NLSERR.LT.0) THEN
            INFORM = NLSERR
            GO TO 60
         END IF
C
         IF (ALFA.GT.ALFLIM) MJRMSG(4:4) = 'L'
C
         ERROR = NLSERR .GE. 4
         IF (ERROR) THEN
C           ---------------------------------------------------------
C           The linesearch failed to find a better point.
C           If exact gradients or central differences are being used,
C           or the KT conditions are satisfied, stop.  Otherwise,
C           switch to central differences and solve the QP again.
C           ---------------------------------------------------------
            IF (NEEDFD .AND. .NOT. CENTRL) THEN
               IF ( .NOT. OPTIML) THEN
                  ERROR = .FALSE.
                  MJRMSG(3:3) = 'Central differences'
                  LVLDIF = 2
                  NEWGQ = .TRUE.
                  IF (MSGNP.GE.5) THEN
                     IF (MINITS.GT.0) THEN
                        WRITE (REC,FMT=99999) MINITS
                        CALL X04BAF(IPRINT,REC(1))
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            IF (NEEDFD) THEN
C              ======================================================
C              Compute the missing gradients.
C              ======================================================
               MODE = 1
               NGRAD = NGRAD + 1
C
               IF (NCNLN.GT.0) THEN
                  CALL F06DBF(NCNLN,(1),IW(LNEEDC),1)
C
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,IW(LNEEDC),X,W(LWRK1),
     *                        CJACU,NSTATE,IUSER,USER)
                  INFORM = MODE
                  IF (MODE.LT.0) GO TO 60
C
                  CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
               END IF
C
               CALL OBJFUN(MODE,N,X,OBJ,GRADU,NSTATE,IUSER,USER)
               INFORM = MODE
               IF (MODE.LT.0) GO TO 60
C
               CALL DCOPY(N,GRADU,1,GRAD,1)
C
               CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                     FDINT,FDNORM,OBJF,CONFUN,OBJFUN,IW(LNEEDC),
     *                     BL,BU,C,W(LWRK2),W(LWRK3),CJAC,CJACU,GRAD,
     *                     GRADU,W(LHFRWD),W(LHCTRL),X,W,LENW,IUSER,
     *                     USER)
C
               INFORM = MODE
               IF (MODE.LT.0) GO TO 60
C
               GDX = DDOT(N,GRAD,1,W(LDX),1)
               GLF2 = GDX
               IF (NCNLN.GT.0) THEN
                  CALL DGEMV('N',NCNLN,N,ONE,CJAC,LDCJ,W(LDX),1,ZERO,
     *                       W(LCJDX),1)
                  GLF2 = GLF2 - DDOT(NCNLN,W(LCJDX),1,CLAMDA(NL),1)
               END IF
            END IF
C
            CALL DCOPY(N,GRAD,1,W(LGQ),1)
            CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,W(LGQ),W(LZY),
     *                  W(LWRK1))
C
            XNORM = DNRM2(N,X,1)
C
            IF (NCNLN.GT.0 .AND. ALFA.GE.ONE) CALL DCOPY(NCNLN,
     *          CLAMDA(NL),1,W(LCMUL),1)
C
            IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,W(LADX),1,AX,1)
            ALFDX = ALFA*DXNORM
C
C           =========================================================
C           Update the factors of the approximate Hessian of the
C           Lagrangian function.
C           =========================================================
            CALL E04UCL(MJRMSG,UNITQ,N,NCNLN,NFREE,NZ,LDCJ1,LDCJ,LDZY,
     *                  LDR,KX,ALFA,GLF1,GLF2,QPCURV,W(LCJAC1),CJAC,
     *                  W(LCJDX1),W(LCJDX),W(LCS1),W(LCS2),W(LGQ1),
     *                  W(LGQ),W(LHPQ),W(LRPQ),CLAMDA(NL),R,W(LWRK3),
     *                  W(LZY),W(LWRK2),W(LWRK1))
C
            CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
            COND = F06BLF(DRMAX,DRMIN,OVERFL)
C
            IF (COND.GT.RCNDBD .OR. RFROBN.GT.ROOTN*GROWTH*DRMAX) THEN
C              ------------------------------------------------------
C              Reset the condition estimator and range-space
C              partition of Q'HQ.
C              ------------------------------------------------------
               IF (NPDBG .AND. INPDBG(1).GT.0) THEN
                  WRITE (REC,FMT=99997) RFROBN, DRMAX, DRMIN, COND,
     *              RCNDBD
                  CALL X04BAY(IPRINT,5,REC)
               END IF
C
               MJRMSG(5:5) = 'Refactorize Hessian'
C
               CALL E04UDR(UNITQ,N,NFREE,NZ,LDZY,LDR,IW(LIPERM),KX,
     *                     W(LGQ),R,W(LZY),W(LWRK1),W(LQRWRK))
            END IF
         END IF
      END IF
C
C     +    UNTIL     (DONE  .OR.  ERROR)
      IF ( .NOT. (DONE .OR. ERROR)) GO TO 20
C
C     ======================end of main loop============================
C
      IF (DONE) THEN
         IF (CONVRG .AND. OPTIML) THEN
            INFORM = 0
         ELSE IF (CONVPT) THEN
            INFORM = 1
         ELSE IF (INFEAS) THEN
            INFORM = 3
         END IF
      ELSE IF (ERROR) THEN
         IF (MAJITS.GE.NMAJOR) THEN
            INFORM = 4
         ELSE IF (OPTIML) THEN
            INFORM = 1
         ELSE
            INFORM = 6
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Set  CLAMDA.  Print the full solution.
C     ------------------------------------------------------------------
   60 MSGNP = MSGSV1
      MSGQP = MSGSV2
      IF (MSGNP.GT.0) THEN
         WRITE (REC,FMT=99998) MAJITS, MNRSUM
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      CALL E04NBX(MSGNP,NFREE,LDAQP,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,
     *            NACTIV,ISTATE,KACTIV,KX,AQP,BL,BU,C,CLAMDA,W(LRLAM),X)
      IF (NCNLN.GT.0) CALL DCOPY(NCNLN,W(LCMUL),1,CLAMDA(N+NCLIN+1),1)
C
      RETURN
C
C
C
C     End of  E04UCZ. (NPCORE)
C
99999 FORMAT (' Mnr itn ',I4,' -- Re-solve QP subproblem.')
99998 FORMAT (/' Exit from NP problem after ',I5,' major iterations,',
     *       /'                            ',I5,' minor iterations.')
99997 FORMAT (/' //E04UCZ//        RFROBN         DRMAX         DRMIN',
     *       /' //E04UCZ//',1P,3D14.2,/' //E04UCZ//          COND     ',
     *       '   RCNDBD',/' //E04UCZ//',1P,2D14.2)
      END
      SUBROUTINE E04MFJ(ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NZ,N,LDQ,LDA,
     *                  LDT,ISTATE,KACTIV,KX,JMAX,ERRMAX,XNORM,A,AX,BL,
     *                  BU,FEATOL,T,X,Q,P,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFJ  computes the point on a working set that is closest in the
C     least-squares sense to the input vector X.
C
C     If the computed point gives a row error of more than the
C     feasibility tolerance, an extra step of iterative refinement is
C     used.  If  X  is still infeasible,  the logical variable ROWERR
C     is set.
C
C     Original version derived from LSSETX January-1987.
C     This version of  E04MFJ  dated   5-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      INTEGER           NTRY
      PARAMETER         (NTRY=5)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERRMAX, XNORM
      INTEGER           JMAX, LDA, LDQ, LDT, N, NACTIV, NCLIN, NFREE, NZ
      LOGICAL           ROWERR, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), P(N), Q(LDQ,*), T(LDT,*),
     *                  WORK(N), X(N)
      INTEGER           ISTATE(N+NCLIN), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, IS, J, K, KTRY
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DDOT, DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRSV, E04MFN, E04NBW,
     *                  F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Move  x  onto the simple bounds in the working set.
C     ------------------------------------------------------------------
      DO 20 K = NFREE + 1, N
         J = KX(K)
         IS = ISTATE(J)
         BND = BL(J)
         IF (IS.GE.2) BND = BU(J)
         IF (IS.NE.4) X(J) = BND
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Move  x  onto the general constraints in the working set.
C     ntry  attempts are made to get acceptable row errors.
C     ------------------------------------------------------------------
      KTRY = 1
      JMAX = 1
      ERRMAX = ZERO
C
C     repeat
   40 IF (NACTIV.GT.0) THEN
C
C        Set work = residuals for constraints in the working set.
C        Solve for P, the smallest correction to x that gives a point
C        on the constraints in the working set.  Define  P = Y*(py),
C        where  py  solves the triangular system  T*(py) = residuals.
C
         DO 60 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(NACTIV-I+1) = BND - DDOT(N,A(K,1),LDA,X,1)
   60    CONTINUE
C
         CALL DTRSV('U','N','N',NACTIV,T(1,NZ+1),LDT,WORK,1)
         CALL F06FBF(N,ZERO,P,1)
         CALL DCOPY(NACTIV,WORK,1,P(NZ+1),1)
C
         CALL E04NBW(2,N,NZ,NFREE,LDQ,UNITQ,KX,P,Q,WORK)
         CALL DAXPY(N,ONE,P,1,X,1)
      END IF
C
C     ---------------------------------------------------------------
C     Compute the 2-norm of  x.
C     Initialize  Ax  for all the general constraints.
C     ---------------------------------------------------------------
      XNORM = DNRM2(N,X,1)
      IF (NCLIN.GT.0) CALL DGEMV('N',NCLIN,N,ONE,A,LDA,X,1,ZERO,AX,1)
C
C     ---------------------------------------------------------------
C     Check the row residuals.
C     ---------------------------------------------------------------
      IF (NACTIV.GT.0) THEN
         DO 80 K = 1, NACTIV
            I = KACTIV(K)
            J = N + I
            IS = ISTATE(J)
            IF (IS.EQ.1) WORK(K) = BL(J) - AX(I)
            IF (IS.GE.2) WORK(K) = BU(J) - AX(I)
   80    CONTINUE
C
         JMAX = IDAMAX(NACTIV,WORK,1)
         ERRMAX = ABS(WORK(JMAX))
      END IF
C
      KTRY = KTRY + 1
C     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      IF ( .NOT. (ERRMAX.LE.FEATOL(JMAX) .OR. KTRY.GT.NTRY)) GO TO 40
C
      ROWERR = ERRMAX .GT. FEATOL(JMAX)
C
      IF (CMDBG .AND. ICMDBG(2).GT.0) CALL E04MFN('E04MFJ','     x',X,N)
C
      RETURN
C
C     End of  E04MFJ.  (CMSETX)
C
      END
      SUBROUTINE E04MFV(CSET,N,NCLIN,LITOTL,LWTOTL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFV   allocates the addresses of the work arrays for E04MFZ.
C
C     Note that the arrays ( GQ, CQ ) lie in contiguous areas of
C     workspace.
C
C     Original version written  2-January-1987.
C     This version of  E04MFV  dated  18-Nov-1990.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
C     .. Scalar Arguments ..
      INTEGER           LITOTL, LWTOTL, N, NCLIN
      LOGICAL           CSET
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LDQ, LDT, LENNAM, LINES1, LINES2,
     *                  NCOLT, NOUT
      LOGICAL           LCDBG
C     .. Arrays in Common ..
      INTEGER           ILCDBG(LDBG), LOCLC(LENLC)
C     .. Local Scalars ..
      INTEGER           LAD, LANORM, LCQ, LD, LENCQ, LENQ, LENRT,
     *                  LFEATU, LGQ, LKACTV, LKX, LQ, LR, LRLAM, LT,
     *                  LWRK, LWTINF, MINIW, MINW
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /EE04MF/ILCDBG, LCDBG
C     .. Save statement ..
      SAVE              /AE04MF/
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Refer to the first free space in the work arrays.
C     ------------------------------------------------------------------
      MINIW = LITOTL + 1
      MINW = LWTOTL + 1
C
C     ------------------------------------------------------------------
C     Integer workspace.
C     ------------------------------------------------------------------
      LKACTV = MINIW
      LKX = LKACTV + N
      MINIW = LKX + N
C
C     ------------------------------------------------------------------
C     Real workspace.
C     Assign array lengths that depend upon the problem dimensions.
C     ------------------------------------------------------------------
      LENRT = LDT*NCOLT
      IF (NCLIN.EQ.0) THEN
         LENQ = 0
      ELSE
         LENQ = LDQ*LDQ
      END IF
C
      IF (CSET) THEN
         LENCQ = N
      ELSE
         LENCQ = 0
      END IF
C
C     ------------------------------------------------------------------
C     We start with arrays that can be preloaded by smart users.
C     ------------------------------------------------------------------
      LFEATU = MINW
      MINW = LFEATU + NCLIN + N
C
C     Next comes stuff used by  E04MFZ  and  E04NFZ.
C
      LANORM = MINW
      LAD = LANORM + NCLIN
      LD = LAD + NCLIN
      LGQ = LD + N
      LCQ = LGQ + N
      LRLAM = LCQ + LENCQ
      LR = LRLAM + N
      LT = LR
      LQ = LT + LENRT
      LWTINF = LQ + LENQ
      LWRK = LWTINF + N + NCLIN
      MINW = LWRK + N + NCLIN
C
C     Load the addresses in LOCLC.
C
      LOCLC(1) = LKACTV
      LOCLC(2) = LKX
C
      LOCLC(3) = LFEATU
      LOCLC(4) = LANORM
      LOCLC(5) = LAD
C
      LOCLC(7) = LD
      LOCLC(8) = LGQ
      LOCLC(9) = LCQ
      LOCLC(10) = LRLAM
      LOCLC(11) = LR
      LOCLC(12) = LT
      LOCLC(13) = LQ
      LOCLC(14) = LWTINF
      LOCLC(15) = LWRK
C
      LITOTL = MINIW - 1
      LWTOTL = MINW - 1
C
      RETURN
C
C     End of  E04MFV.  (LPLOC)
C
      END
      SUBROUTINE E04MFW(N,NCLIN,TITLE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFW loads the default values of parameters not set by the user.
C
C     Original Fortran 77 version written 30-December-1986.
C     This version of  E04MFW  dated  11-Nov-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D+0)
      DOUBLE PRECISION  WRKTOL
      PARAMETER         (WRKTOL=1.0D-2)
C     .. Scalar Arguments ..
      INTEGER           N, NCLIN
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, EPSPT3, EPSPT5,
     *                  EPSPT8, EPSPT9, TOLACT, TOLFEA, TOLINC, TOLRNK,
     *                  TOLX0
      INTEGER           IDBGLC, IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1,
     *                  ITMAX2, ITNFIX, KCHK, KCYCLE, KDEGEN, LCRASH,
     *                  LDBGLC, LINES1, LINES2, LPROB, MAXACT, MAXNZ,
     *                  MM, MSGLC, MXFREE, NDEGEN, NN, NNCLIN, NOUT,
     *                  NPROB
      LOGICAL           CMDBG, LCDBG, NEWOPT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), ILCDBG(LDBG), IPADLC(14),
     *                  IPSVLC(MXPARM), NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH
      INTEGER           I, IDBG, J, K, LENT, MSG, MSGDBG, MSGLVL
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*4       ICRSH(0:2)
      CHARACTER*7       LPTYPE(1:10)
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04MFX, F06DFF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         LEN, MAX, MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04MF/NEWOPT
      COMMON            /CE04MF/TOLX0, TOLINC, KDEGEN, NDEGEN, ITNFIX,
     *                  NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /EE04MF/ILCDBG, LCDBG
      COMMON            /FE04MF/IPSVLC, IDBGLC, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, KCHK, KCYCLE, LCRASH, LPROB, MAXACT,
     *                  MXFREE, MAXNZ, MM, LDBGLC, MSGLC, NN, NNCLIN,
     *                  NPROB, IPADLC
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLC(1),IDBGLC), (RPRMLC(1),BIGBND)
      EQUIVALENCE       (MSGLC,MSGLVL), (IDBGLC,IDBG), (LDBGLC,MSGDBG)
C     .. Save statement ..
      SAVE              /AX02ZA/, /BE04MF/, /FE04MF/, /GE04MF/
C     .. Data statements ..
      DATA              ICRSH(0), ICRSH(1), ICRSH(2)/'COLD', 'WARM',
     *                  'HOT '/
      DATA              LPTYPE(1), LPTYPE(2)/'     FP', '     LP'/
      DATA              LPTYPE(3), LPTYPE(4), LPTYPE(5),
     *                  LPTYPE(6)/'ILLEGAL', 'ILLEGAL', 'ILLEGAL',
     *                  'ILLEGAL'/
      DATA              LPTYPE(7), LPTYPE(8), LPTYPE(9),
     *                  LPTYPE(10)/'       ', '       ', '       ',
     *                  'ILLEGAL'/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
C     Make a dummy call to E04MFX to ensure that the defaults are set.
C
      CALL E04MFX(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     RPRMLC and IPRMLC may be changed to their default values.
C
      CALL F06DFF(MXPARM,IPRMLC,1,IPSVLC,1)
      CALL DCOPY(MXPARM,RPRMLC,1,RPSVLC,1)
C
      IF (MSGLVL.EQ.IDUMMY) MSGLVL = 10
      IF (IPRNT.LT.0) IPRNT = NOUT
      IF (ISUMRY.LT.0 .OR. MSGLVL.LT.5) ISUMRY = -1
      IPRINT = IPRNT
      ISUMM = ISUMRY
      IF (KCHK.LE.0) KCHK = 50
      IF (KCYCLE.LE.0) KCYCLE = 10000
      IF (KCYCLE.GT.9999999) KCYCLE = 9999999
      KDEGEN = KCYCLE
      IF (LPROB.LT.0) LPROB = 2
      IF (LCRASH.LT.0 .OR. LCRASH.GT.2) LCRASH = 0
      IF (ITMAX1.LT.0) ITMAX1 = MAX(50,5*(N+NCLIN))
      IF (ITMAX2.LT.0) ITMAX2 = MAX(50,5*(N+NCLIN))
      IF (MAXACT.LT.0 .OR. MAXACT.GT.N .OR. MAXACT.GT.NCLIN)
     *    MAXACT = MAX(1,MIN(N,NCLIN))
      IF (MAXNZ.LT.0 .OR. MAXNZ.GT.N) MAXNZ = N
      IF (MXFREE.LT.0 .OR. MXFREE.GT.N) MXFREE = N
      IF (MXFREE.LT.MAXNZ) MXFREE = MAXNZ
      IF (NCLIN.LT.N) THEN
         MXFREE = NCLIN + 1
         MAXNZ = MXFREE
      END IF
C
      IF (IDBG.LT.0 .OR. IDBG.GT.ITMAX1+ITMAX2) IDBG = 0
      IF (MSGDBG.LT.0) MSGDBG = 0
      IF (MSGDBG.EQ.0) IDBG = ITMAX1 + ITMAX2 + 1
      IF (TOLACT.LT.ZERO) TOLACT = WRKTOL
      IF (TOLFEA.EQ.RDUMMY .OR. (TOLFEA.GE.ZERO .AND. TOLFEA.LT.EPSMCH))
     *    TOLFEA = EPSPT5
      IF (BIGBND.LE.ZERO) BIGBND = GIGANT
      IF (BIGDX.LE.ZERO) BIGDX = MAX(GIGANT,BIGBND)
C
      LCDBG = IDBG .EQ. 0
      CMDBG = LCDBG
      K = 1
      MSG = MSGDBG
      DO 20 I = 1, LDBG
         ILCDBG(I) = MOD(MSG/K,10)
         ICMDBG(I) = ILCDBG(I)
         K = K*10
   20 CONTINUE
C
      IF (MSGLVL.GT.0) THEN
C
C        Print the title.
C
         LENT = LEN(TITLE)
         WRITE (REC,FMT=99993) (TITLE(J:J),J=1,LENT)
         CALL X04BAY(IPRINT,2,REC)
         CALL A00AAF
C
         IF (MSGLVL.GE.5 .AND. ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
            WRITE (REC,FMT=99992) (TITLE(J:J),J=1,LENT)
            CALL X04BAY(ISUMM,2,REC)
         END IF
C
         WRITE (REC,FMT=99999)
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99998) LPTYPE(LPROB)
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99997) NCLIN, TOLFEA, N, TOLACT
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99996) BIGBND, ICRSH(LCRASH), BIGDX, EPSMCH
         CALL X04BAY(IPRINT,3,REC)
         WRITE (REC,FMT=99995) KCHK, KDEGEN
         CALL X04BAY(IPRINT,2,REC)
         WRITE (REC,FMT=99994) MSGLVL, ITMAX2, ISUMRY
         CALL X04BAY(IPRINT,3,REC)
      END IF
      RETURN
C
C
C     End of  E04MFW.  (LPDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Problem type...........',3X,A7)
99997 FORMAT (/' Linear constraints.....',I10,7X,'Feasibility toleranc',
     *       'e..',1P,D10.2,/' Variables..............',I10,7X,'Crash ',
     *       'tolerance........',1P,D10.2)
99996 FORMAT (/' Infinite bound size....',1P,D10.2,7X,A4,' start......',
     *       '.......',/' Infinite step size.....',1P,D10.2,7X,'EPS (m',
     *       'achine precision)',1P,D10.2)
99995 FORMAT (/' Check frequency........',I10,7X,'Expand frequency....',
     *       '...',I10)
99994 FORMAT (/' Print level............',I10,7X,'Iteration limit.....',
     *       '...',I10,/' Monitoring file........',I10)
99993 FORMAT (/80A1)
99992 FORMAT (/11A1,' monitoring information ')
      END
      SUBROUTINE E04NFQ(UNITQ,VERTEX,K1,K2,IT,NACTIV,NARTIF,NZ,NFREE,
     *                  NREJTD,NGQ,N,LDQ,LDA,LDT,ISTATE,KACTIV,KX,
     *                  CONDMX,A,T,GQM,Q,W,C,S,MSGLVL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04NFQ  includes general constraints  K1  thru  K2  as new rows of
C     the  TQ  factorization:
C              A(free) * Q(free)  = (  0 T )
C                        Q(free)  = (  Z Y )
C
C     a) The  NACTIV x NACTIV  upper-triangular matrix  T  is stored
C        with its (1,1) element in position  (IT,JT)  of the array  T.
C
C     Original version written by PEG,  October-31-1984.
C     This version of  E04NFQ  dated  7-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           IT, K1, K2, LDA, LDQ, LDT, MSGLVL, N, NACTIV,
     *                  NARTIF, NFREE, NGQ, NREJTD, NZ
      LOGICAL           UNITQ, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQM(N,*), Q(LDQ,*), S(N),
     *                  T(LDT,*), W(N)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8,
     *                  EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNDMAX, COND, DELTA, DRZZ, DTNEW, RNORM, ROWMAX,
     *                  RTMAX, TDTMAX, TDTMIN
      INTEGER           I, IADD, IARTIF, IFIX, INFORM, ISWAP, J, JADD,
     *                  JT, K, L, NZADD
      LOGICAL           OVERFL, RSET
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DGER, E04NBW, E04NFR, F06FLF,
     *                  F06FRF, F06QHF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      RTMAX = WMACH(8)
C
      JT = NZ + 1
C
C     Estimate the condition number of the constraints already
C     factorized.
C
      IF (NACTIV.EQ.0) THEN
         DTMAX = ZERO
         DTMIN = ONE
         IF (UNITQ) THEN
C
C           First general constraint added.  Set  Q = I.
C
            CALL F06QHF('General',NFREE,NFREE,ZERO,ONE,Q,LDQ)
            UNITQ = .FALSE.
         END IF
      ELSE
         CALL F06FLF(NACTIV,T(IT,JT),LDT+1,DTMAX,DTMIN)
      END IF
C
      DO 20 K = K1, K2
         IADD = KACTIV(K)
         JADD = N + IADD
         IF (NACTIV.LT.NFREE) THEN
C
            IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
               WRITE (REC,FMT=99999) NACTIV, NZ, NFREE, IADD, JADD
               CALL X04BAY(IPRINT,4,REC)
            END IF
C
            OVERFL = .FALSE.
C
C           Transform the incoming row of  A  by  Q'.
C
            CALL DCOPY(N,A(IADD,1),LDA,W,1)
            CALL E04NBW(8,N,NZ,NFREE,LDQ,UNITQ,KX,W,Q,S)
C
C           Check that the incoming row is not dependent upon those
C           already in the working set.
C
            DTNEW = DNRM2(NZ,W,1)
            IF (NACTIV.EQ.0) THEN
C
C              This is the first general constraint in the working set.
C
               COND = F06BLF(ASIZE,DTNEW,OVERFL)
               TDTMAX = DTNEW
               TDTMIN = DTNEW
            ELSE
C
C              There are already some general constraints in the working
C              set. Update the estimate of the condition number.
C
               TDTMAX = MAX(DTNEW,DTMAX)
               TDTMIN = MIN(DTNEW,DTMIN)
               COND = F06BLF(TDTMAX,TDTMIN,OVERFL)
            END IF
C
            IF (COND.GE.CONDMX .OR. OVERFL) THEN
C              ---------------------------------------------------------
C              This constraint appears to be dependent on those already
C              in the working set.  Skip it.
C              ---------------------------------------------------------
               IF (CMDBG .AND. ICMDBG(1).GT.0) THEN
                  WRITE (REC,FMT=99998)
                  CALL X04BAY(IPRINT,2,REC)
                  IF (NACTIV.GT.0) THEN
                     WRITE (REC,FMT=99997) ASIZE, DTMAX, DTMIN, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  ELSE IF (NACTIV.EQ.0) THEN
                     WRITE (REC,FMT=99996) ASIZE, DTNEW
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
C
               ISTATE(JADD) = 0
               KACTIV(K) = -KACTIV(K)
            ELSE
               IF (NZ.GT.1) THEN
C                 ------------------------------------------------------
C                 Use a single column transformation to reduce the first
C                 NZ-1  elements of  W  to zero.
C                 ------------------------------------------------------
C                 Apply the Householder reflection  I  -  W W'.
C                 The reflection is applied to  Z  and GQM so that
C                    Y  =    Z  * W,   Z    =  Z    -  Y W'  and
C                    Y  =  GQM' * W,   GQM  =  GQM  -  W Y',
C                 where  W = WRK1 (from Householder),
C                 and    Y = WRK2 (workspace).
C
C                 Note that DELTA  has to be stored after the reflection
C                 is used.
C
                  DELTA = W(NZ)
                  CALL F06FRF(NZ-1,DELTA,W,1,ZERO,W(NZ))
                  IF (W(NZ).GT.ZERO) THEN
C
                     CALL DGEMV('N',NFREE,NZ,ONE,Q,LDQ,W,1,ZERO,S,1)
                     CALL DGER(NFREE,NZ,(-ONE),S,1,W,1,Q,LDQ)
C
                     IF (NGQ.GT.0) THEN
                        CALL DGEMV('T',NZ,NGQ,ONE,GQM,N,W,1,ZERO,S,1)
                        CALL DGER(NZ,NGQ,(-ONE),W,1,S,1,GQM,N)
                     END IF
                  END IF
C
                  W(NZ) = DELTA
               END IF
               IT = IT - 1
               JT = JT - 1
               NACTIV = NACTIV + 1
               NZ = NZ - 1
               CALL DCOPY(NACTIV,W(JT),1,T(IT,JT),LDT)
               DTMAX = TDTMAX
               DTMIN = TDTMIN
            END IF
         END IF
   20 CONTINUE
C
      IF (NACTIV.LT.K2) THEN
C
C        Some of the constraints were classed as dependent and not
C        included in the factorization.  Re-order the part of  KACTIV
C        that holds the indices of the general constraints in the
C        working set.  Move accepted indices to the front and shift
C        rejected indices (with negative values) to the end.
C
         L = K1 - 1
         DO 40 K = K1, K2
            I = KACTIV(K)
            IF (I.GE.0) THEN
               L = L + 1
               IF (L.NE.K) THEN
                  ISWAP = KACTIV(L)
                  KACTIV(L) = I
                  KACTIV(K) = ISWAP
               END IF
            END IF
   40    CONTINUE
C
C        If a vertex is required,  add some temporary bounds.
C        We must accept the resulting condition number of the working
C        set.
C
         IF (VERTEX) THEN
            RSET = .FALSE.
            CNDMAX = RTMAX
            DRZZ = ONE
            NZADD = NZ
            DO 80 IARTIF = 1, NZADD
               IF (UNITQ) THEN
                  IFIX = NFREE
                  JADD = KX(IFIX)
               ELSE
                  ROWMAX = ZERO
                  DO 60 I = 1, NFREE
                     RNORM = DNRM2(NZ,Q(I,1),LDQ)
                     IF (ROWMAX.LT.RNORM) THEN
                        ROWMAX = RNORM
                        IFIX = I
                     END IF
   60             CONTINUE
                  JADD = KX(IFIX)
C
                  CALL E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,
     *                        NACTIV,NZ,NFREE,NZ,NGQ,N,LDA,LDQ,LDT,LDT,
     *                        KX,CNDMAX,DRZZ,A,T,T,GQM,Q,W,C,S,MSGLVL)
               END IF
               NFREE = NFREE - 1
               NZ = NZ - 1
               NARTIF = NARTIF + 1
               ISTATE(JADD) = 4
   80       CONTINUE
         END IF
C
         IF (IT.GT.1) THEN
C           ------------------------------------------------------------
C           If some dependent constraints were rejected,  move  T  to
C           the top of the array  T.
C           ------------------------------------------------------------
            DO 120 K = 1, NACTIV
               J = NZ + K
               DO 100 I = 1, K
                  T(I,J) = T(IT+I-1,J)
  100          CONTINUE
  120       CONTINUE
         END IF
      END IF
C
      NREJTD = K2 - NACTIV
C
      RETURN
C
C     End of  E04NFQ.  (RZADDS)
C
99999 FORMAT (/' //E04NFQ//  Constraint added.           ',/' //E04NFQ',
     *       '//  NACTIV    NZ NFREE  IADD  JADD',/' //E04NFQ//  ',5I6)
99998 FORMAT (/' //E04NFQ//  Dependent constraint rejected.')
99997 FORMAT (/' //E04NFQ//     ASIZE     DTMAX     DTMIN     DTNEW',
     *       /' //E04NFQ//',1P,4D10.2)
99996 FORMAT (/' //E04NFQ//     ASIZE     DTNEW',/' //E04NFQ//',1P,
     *       2D10.2)
      END
      SUBROUTINE E04UCX(N,NCLIN,NCNLN,TITLE)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1091 (JUL 1993).
C
C     ******************************************************************
C     E04UCX  loads the default values of parameters not set in the
C     options file.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 77 version written 10-September-1985.
C     This version of E04UCX dated 1-Feb-93.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  POINT3, POINT8
      PARAMETER         (POINT3=3.3D-1,POINT8=0.8D+0)
      DOUBLE PRECISION  POINT9, TWO
      PARAMETER         (POINT9=0.9D+0,TWO=2.0D+0)
      DOUBLE PRECISION  TENP6, HUNDRD
      PARAMETER         (TENP6=1.0D+6,HUNDRD=10.0D+1)
      DOUBLE PRECISION  RDUMMY
      INTEGER           IDUMMY
      PARAMETER         (RDUMMY=-11111.D+0,IDUMMY=-11111)
      DOUBLE PRECISION  GIGANT
      PARAMETER         (GIGANT=1.0D+20*0.99999D+0)
      DOUBLE PRECISION  WRKTOL
      PARAMETER         (WRKTOL=1.0D-2)
C     .. Scalar Arguments ..
      INTEGER           N, NCLIN, NCNLN
      CHARACTER*(*)     TITLE
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT, CTOL,
     *                  DXLIM, EPSPT3, EPSPT5, EPSPT8, EPSPT9, EPSRF,
     *                  ETA, FDINT, FTOL, HCNDBD, TOLACT, TOLFEA, TOLRNK
      INTEGER           IDBGLS, IDBGNP, IPRINT, IPRNT, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, KSAVE, LCRASH, LDBGLS, LDBGNP, LFDSET,
     *                  LFORMH, LINES1, LINES2, LPROB, LVERFY, LVLDER,
     *                  LVLDIF, LVRFYC, MSGLS, MSGNP, NCDIFF, NFDIFF,
     *                  NLNF, NLNJ, NLNX, NLOAD, NN, NNCLIN, NNCNLN,
     *                  NOUT, NPROB, NSAVE
      LOGICAL           CMDBG, LSDBG, NEWOPT, NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), ILSDBG(LDBG), INPDBG(LDBG),
     *                  IPADLS(18), IPADNP(12), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM), JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONDBD, DCTOL, EPSMCH
      INTEGER           I, IDBG, J, K, LENT, MJRDBG, MNRDBG, MSG1, MSG2,
     *                  MSGQP, NCTOTL, NMAJOR, NMINOR, NPLIN
      CHARACTER*16      KEY
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      CHARACTER*3       CHESS(0:1)
      CHARACTER*4       ICRSH(0:2)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          A00AAF, DCOPY, E04UCQ, F06DFF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LEN, MAX, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /DE04NC/IPSVLS, IDBGLS, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LPROB, MSGLS, NN,
     *                  NNCLIN, NPROB, IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /EE04UC/NEWOPT
      COMMON            /FE04NB/ICMDBG, CMDBG
      COMMON            /FE04UC/INPDBG, NPDBG
      COMMON            /GE04UC/IPSVNP, IDBGNP, ITMXNP, JVRFY1, JVRFY2,
     *                  JVRFY3, JVRFY4, LDBGNP, LFORMH, LVLDER, LVERFY,
     *                  MSGNP, NLNF, NLNJ, NLNX, NNCNLN, NSAVE, NLOAD,
     *                  KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IDBGLS), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),IDBGNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IDBGNP,IDBG), (ITMXNP,NMAJOR), (ITMAX2,NMINOR)
      EQUIVALENCE       (LDBGLS,MNRDBG), (LDBGNP,MJRDBG), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /AX02ZA/, /EE04UC/, /DE04NC/, /EE04NC/,
     *                  /GE04UC/, /HE04UC/
C     .. Data statements ..
      DATA              ICRSH(0), ICRSH(1), ICRSH(2)/'COLD', 'WARM',
     *                  'HOT '/
      DATA              CHESS(0), CHESS(1)/' NO', 'YES'/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
      NOUT = WMACH(11)
C
      CONDBD = MAX(ONE/(HUNDRD*EPSMCH*DBLE(N)),TENP6)
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
C
C     Make a dummy call E04UCQ to ensure that the defaults are set.
C
      CALL E04UCQ(NOUT,'*',KEY)
      NEWOPT = .TRUE.
C
C     Save the optional parameters set by the user.  The values in
C     IPRMLS, RPRMLS, IPRMNP and RPRMNP may be changed to their
C     default values.
C
      CALL F06DFF(MXPARM,IPRMLS,1,IPSVLS,1)
      CALL DCOPY(MXPARM,RPRMLS,1,RPSVLS,1)
      CALL F06DFF(MXPARM,IPRMNP,1,IPSVNP,1)
      CALL DCOPY(MXPARM,RPRMNP,1,RPSVNP,1)
C
      IF (MSGNP.EQ.IDUMMY) MSGNP = 10
      IF (MSGQP.EQ.IDUMMY) MSGQP = 0
      IF (IPRNT.LT.0) IPRNT = NOUT
      IF (ISUMRY.LT.0 .OR. (MSGNP.LT.5 .AND. MSGQP.LT.5)) ISUMRY = -1
      IPRINT = IPRNT
      ISUMM = ISUMRY
      IF (LCRASH.LT.0 .OR. LCRASH.GT.2) LCRASH = 0
      IF (LVLDER.LT.0 .OR. LVLDER.GT.3) LVLDER = 3
      IF (LFORMH.LT.0 .OR. LFORMH.GT.1) LFORMH = 0
C
      IF (NMAJOR.LT.0) NMAJOR = MAX(50,3*NPLIN+10*NCNLN)
      IF (NMINOR.LT.1) NMINOR = MAX(50,3*NCTOTL)
      IF (MJRDBG.LT.0) MJRDBG = 0
      IF (MNRDBG.LT.0) MNRDBG = 0
      IF (IDBG.LT.0 .OR. IDBG.GT.NMAJOR) IDBG = 0
      IF (MJRDBG.EQ.0 .AND. MNRDBG.EQ.0) IDBG = NMAJOR + 1
      NLNF = N
      NLNJ = N
      NLNX = N
      IF (JVRFY2.LE.0 .OR. JVRFY2.GT.N) JVRFY2 = N
      IF (JVRFY1.LE.0 .OR. JVRFY1.GT.JVRFY2) JVRFY1 = 1
      IF (JVRFY4.LE.0 .OR. JVRFY4.GT.N) JVRFY4 = N
      IF (JVRFY3.LE.0 .OR. JVRFY3.GT.JVRFY4) JVRFY3 = 1
      IF ((LVERFY.LT.-1 .OR. LVERFY.GT.13)
     *    .OR. (LVERFY.GE.4 .AND. LVERFY.LE.9)) LVERFY = 0
C
      IF (KSAVE.LE.0) KSAVE = NMAJOR + 1
      IF (NSAVE.LT.0) NSAVE = 0
      IF (NSAVE.EQ.0) KSAVE = NMAJOR + 1
      IF (NLOAD.LT.0) NLOAD = 0
      IF (LCRASH.LE.1) NLOAD = 0
      IF (NLOAD.EQ.0 .AND. LCRASH.EQ.2) LCRASH = 0
C
      IF (TOLACT.LT.ZERO .OR. TOLACT.GE.ONE) TOLACT = WRKTOL
      IF (TOLFEA.LT.EPSMCH .OR. TOLFEA.GE.ONE) TOLFEA = EPSPT5
      IF (EPSRF.LT.EPSMCH .OR. EPSRF.GE.ONE) EPSRF = EPSPT9
      LFDSET = 0
      IF (FDINT.LT.ZERO) LFDSET = 2
      IF (FDINT.EQ.RDUMMY) LFDSET = 0
      IF (FDINT.GE.EPSMCH .AND. FDINT.LT.ONE) LFDSET = 1
      IF (LFDSET.EQ.1 .AND. (CDINT.LT.EPSMCH .OR. CDINT.GE.ONE))
     *    CDINT = EPSRF**POINT3
      IF (BIGBND.LE.ZERO) BIGBND = GIGANT
      IF (BIGDX.LE.ZERO) BIGDX = MAX(GIGANT,BIGBND)
      IF (DXLIM.LE.ZERO) DXLIM = TWO
      IF (ETA.LT.ZERO .OR. ETA.GE.ONE) ETA = POINT9
      IF (FTOL.LT.EPSRF .OR. FTOL.GE.ONE) FTOL = EPSRF**POINT8
C
      IF (HCNDBD.LT.ONE) HCNDBD = CONDBD
C
      DCTOL = EPSPT5
      IF (LVLDER.LT.2) DCTOL = EPSPT3
      IF (CTOL.LT.EPSMCH .OR. CTOL.GE.ONE) CTOL = DCTOL
C
      ITMAX1 = MAX(50,3*(N+NCLIN+NCNLN))
      JVERFY(1) = JVRFY1
      JVERFY(2) = JVRFY2
      JVERFY(3) = JVRFY3
      JVERFY(4) = JVRFY4
C
      NPDBG = IDBG .EQ. 0
      CMDBG = NPDBG
      LSDBG = NPDBG
C
      K = 1
      MSG1 = MJRDBG
      MSG2 = MNRDBG
      DO 20 I = 1, LDBG
         INPDBG(I) = MOD(MSG1/K,10)
         ICMDBG(I) = INPDBG(I)
         ILSDBG(I) = MOD(MSG2/K,10)
         K = K*10
   20 CONTINUE
C
      IF (MSGNP.GT.0) THEN
C
C        Print the title. If no hot start is specified, the parameters
C        are final and can be printed.
C
         LENT = LEN(TITLE)
         WRITE (REC,FMT=99988) (TITLE(J:J),J=1,LENT)
         CALL X04BAY(IPRINT,2,REC)
         CALL A00AAF
C
         IF (LCRASH.LE.1) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAY(IPRINT,3,REC)
            WRITE (REC,FMT=99998) NCLIN, TOLFEA, N, TOLACT
            CALL X04BAY(IPRINT,3,REC)
            WRITE (REC,FMT=99997) BIGBND, ICRSH(LCRASH), BIGDX, EPSMCH,
     *        DXLIM, CHESS(LFORMH)
            CALL X04BAY(IPRINT,4,REC)
            WRITE (REC,FMT=99996) NCNLN, CTOL, NLNF, FTOL, NLNJ, ETA
            CALL X04BAY(IPRINT,4,REC)
            WRITE (REC,FMT=99995) LVLDER, EPSRF, LVERFY, ISUMRY
            CALL X04BAY(IPRINT,3,REC)
            IF (LVERFY.GT.0) THEN
               WRITE (REC,FMT=99994) JVRFY1, JVRFY2
               CALL X04BAY(IPRINT,2,REC)
               IF (NCNLN.GT.0) THEN
                  WRITE (REC,FMT=99993) JVRFY3, JVRFY4
                  CALL X04BAF(IPRINT,REC(1))
               END IF
            END IF
            WRITE (REC,FMT=99992) NMAJOR, MSGNP, NMINOR, MSGQP
            CALL X04BAY(IPRINT,3,REC)
C
            IF (LVLDER.LT.3) THEN
               IF (LFDSET.EQ.0) THEN
                  WRITE (REC,FMT=99991)
                  CALL X04BAY(IPRINT,2,REC)
               ELSE IF (LFDSET.EQ.1) THEN
                  WRITE (REC,FMT=99990) FDINT, CDINT
                  CALL X04BAY(IPRINT,2,REC)
               ELSE IF (LFDSET.EQ.2) THEN
                  WRITE (REC,FMT=99989)
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            END IF
C
         END IF
      END IF
C
      IF (MSGNP.GE.5 .OR. MSGQP.GE.5) THEN
         IF (ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
            LENT = LEN(TITLE)
            WRITE (REC,FMT=99987) (TITLE(J:J),J=1,LENT)
            CALL X04BAY(ISUMM,2,REC)
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04UCX. (NPDFLT)
C
99999 FORMAT (/' Parameters',/' ----------')
99998 FORMAT (/' Linear constraints.....',I10,7X,'Linear feasibility..',
     *       '...',1P,D10.2,/' Variables..............',I10,7X,'Crash ',
     *       'tolerance........',1P,D10.2)
99997 FORMAT (/' Infinite bound size....',1P,D10.2,7X,A4,' start......',
     *       '.......',/' Infinite step size.....',1P,D10.2,7X,'EPS (m',
     *       'achine precision)',1P,D10.2,/' Step limit.............',
     *       1P,D10.2,7X,'Hessian................',7X,A3)
99996 FORMAT (/' Nonlinear constraints..',I10,7X,'Nonlinear feasibilit',
     *       'y..',1P,D10.2,/' Nonlinear objectiv vars',I10,7X,'Optima',
     *       'lity tolerance...',1P,D10.2,/' Nonlinear Jacobian vars',
     *       I10,7X,'Linesearch tolerance...',1P,D10.2)
99995 FORMAT (/' Derivative level.......',I10,7X,'Function precision..',
     *       '...',1P,D10.2,/' Verify level...........',I10,7X,'Monito',
     *       'ring file........',I10)
99994 FORMAT (/' Start obj chk at varble',I10,7X,'Stop obj chk at varb',
     *       'le.',I10)
99993 FORMAT (' Start con chk at varble',I10,7X,'Stop con chk at varbl',
     *       'e.',I10)
99992 FORMAT (/' Major iterations limit.',I10,7X,'Major print level...',
     *       '...',I10,/' Minor iterations limit.',I10,7X,'Minor print',
     *       ' level......',I10)
99991 FORMAT (/' Difference intervals to be computed.')
99990 FORMAT (/' Difference interval....',1P,D10.2,7X,'Central diffce ',
     *       'interval',1P,D10.2)
99989 FORMAT (/' User-supplied difference intervals.')
99988 FORMAT (/80A1)
99987 FORMAT (/11A1,' monitoring information ')
      END
      SUBROUTINE E04NCU(COLD,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,
     *                  LDA,ISTATE,KACTIV,BIGBND,TOLACT,A,AX,BL,BU,X,WX,
     *                  WORK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-720 (DEC 1989).
C     MARK 14C REVISED. IER-890 (NOV 1990).
C     MARK 16A REVISED. IER-996 (JUN 1993).
C     ******************************************************************
C     E04NCU  computes the quantities  ISTATE (optionally), KACTIV,
C     NACTIV, nz and NFREE  associated with the working set at X.
C     The computation depends upon the value of the input parameter
C     COLD,  as follows...
C
C     COLD = TRUE.  An initial working set will be selected. First,
C                   nearly-satisfied or violated bounds are added.
C                   Next,  general linear constraints are added that
C                   have small residuals.
C
C     COLD = FALSE. The quantities KACTIV, NACTIV, nz and NFREE are
C                   computed from ISTATE,  specified by the user.
C
C     Values of ISTATE(j)....
C
C        - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCU dated 14-May-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLACT
      INTEGER           LDA, N, NACTIV, NARTIF, NCLIN, NCTOTL, NFREE
      LOGICAL           COLD, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  WORK(N), WX(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LSDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           ILSDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, BIGLOW, BIGUPP, COLMIN, COLSIZ, FLMAX,
     *                  RESIDL, RESL, RESMIN, RESU, TOOBIG
      INTEGER           I, IMIN, IS, J, JMIN, K, NFIXED
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04NC/ILSDBG, LSDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      FLMAX = WMACH(7)
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ------------------------------------------------------------------
C     Move the variables inside their bounds.
C     ------------------------------------------------------------------
C
      DO 20 J = 1, N
         B1 = BL(J)
         B2 = BU(J)
C
         IF (B1.GT.BIGLOW) THEN
            IF (X(J).LT.B1) X(J) = B1
         END IF
C
         IF (B2.LT.BIGUPP) THEN
            IF (X(J).GT.B2) X(J) = B2
         END IF
   20 CONTINUE
C
      CALL DCOPY(N,X,1,WX,1)
C
      IF (LSDBG) THEN
         IF (ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) COLD, NCLIN, NCTOTL
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (ILSDBG(2).GT.0) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
            DO 40 I = 1, N, 5
               WRITE (REC,FMT=99995) (WX(J),J=I,MIN(I+4,N))
               CALL X04BAF(IPRINT,REC(1))
   40       CONTINUE
         END IF
      END IF
C
      NFIXED = 0
      NACTIV = 0
      NARTIF = 0
C
C     If a cold start is being made, initialize  ISTATE.
C     If  BL(j) = BU(j),  set  ISTATE(j)=3  for all variables and linear
C     constraints.
C
      IF (COLD) THEN
         DO 60 J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J).EQ.BU(J)) ISTATE(J) = 3
   60    CONTINUE
      ELSE
         DO 80 J = 1, NCTOTL
            IF (ISTATE(J).GT.3 .OR. ISTATE(J).LT.0) ISTATE(J) = 0
            IF (BL(J).NE.BU(J) .AND. ISTATE(J).EQ.3) ISTATE(J) = 0
   80    CONTINUE
      END IF
C
C     Initialize NFIXED, NFREE and KACTIV.
C     Ensure that the number of bounds and general constraints in the
C     working set does not exceed N.
C
      DO 100 J = 1, NCTOTL
         IF (NFIXED+NACTIV.EQ.N) ISTATE(J) = 0
         IF (ISTATE(J).GT.0) THEN
            IF (J.LE.N) THEN
               NFIXED = NFIXED + 1
               IF (ISTATE(J).EQ.1) WX(J) = BL(J)
               IF (ISTATE(J).GE.2) WX(J) = BU(J)
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
  100 CONTINUE
C
C     ------------------------------------------------------------------
C     If a cold start is required,  attempt to add as many
C     constraints as possible to the working set.
C     ------------------------------------------------------------------
      IF (COLD) THEN
C
C        See if any bounds are violated or nearly satisfied.
C        If so,  add these bounds to the working set and set the
C        variables exactly on their bounds.
C
         J = N
C        +       WHILE (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) DO
  120    IF (J.GE.1 .AND. NFIXED+NACTIV.LT.N) THEN
            IF (ISTATE(J).EQ.0) THEN
               B1 = BL(J)
               B2 = BU(J)
               IS = 0
               IF (B1.GT.BIGLOW) THEN
                  IF (WX(J)-B1.LE.(ONE+ABS(B1))*TOLACT) IS = 1
               END IF
               IF (B2.LT.BIGUPP) THEN
                  IF (B2-WX(J).LE.(ONE+ABS(B2))*TOLACT) IS = 2
               END IF
               IF (IS.GT.0) THEN
                  ISTATE(J) = IS
                  IF (IS.EQ.1) WX(J) = B1
                  IF (IS.EQ.2) WX(J) = B2
                  NFIXED = NFIXED + 1
               END IF
            END IF
            J = J - 1
            GO TO 120
C           +       END WHILE
         END IF
C
C        ---------------------------------------------------------------
C        The following loop finds the linear constraint (if any) with
C        smallest residual less than or equal to TOLACT  and adds it
C        to the working set.  This is repeated until the working set
C        is complete or all the remaining residuals are too large.
C        ---------------------------------------------------------------
C        First, compute the residuals for all the constraints not in the
C        working set.
C
         IF (NCLIN.GT.0 .AND. NACTIV+NFIXED.LT.N) THEN
            DO 140 I = 1, NCLIN
               IF (ISTATE(N+I).LE.0) AX(I) = DDOT(N,A(I,1),LDA,WX,1)
  140       CONTINUE
C
            IS = 1
            TOOBIG = TOLACT + TOLACT
C
C           + WHILE (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) DO
  160       IF (IS.GT.0 .AND. NFIXED+NACTIV.LT.N) THEN
               IS = 0
               RESMIN = TOLACT
C
               DO 180 I = 1, NCLIN
                  J = N + I
                  IF (ISTATE(J).EQ.0) THEN
                     B1 = BL(J)
                     B2 = BU(J)
                     RESL = TOOBIG
                     RESU = TOOBIG
                     IF (B1.GT.BIGLOW) RESL = ABS(AX(I)-B1)/(ONE+ABS(B1)
     *                                        )
                     IF (B2.LT.BIGUPP) RESU = ABS(AX(I)-B2)/(ONE+ABS(B2)
     *                                        )
                     RESIDL = MIN(RESL,RESU)
                     IF (RESIDL.LT.RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN = I
                        IS = 1
                        IF (RESL.GT.RESU) IS = 2
                     END IF
                  END IF
  180          CONTINUE
C
               IF (IS.GT.0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 160
C              +          END WHILE
            END IF
         END IF
      END IF
C
      IF (VERTEX .AND. NACTIV+NFIXED.LT.N) THEN
C        ---------------------------------------------------------------
C        Find an initial vertex by temporarily fixing some variables.
C        ---------------------------------------------------------------
C        Compute lengths of columns of selected linear constraints
C        (just the ones corresponding to variables eligible to be
C        temporarily fixed).
C
         DO 220 J = 1, N
            IF (ISTATE(J).EQ.0) THEN
               COLSIZ = ZERO
               DO 200 K = 1, NCLIN
                  IF (ISTATE(N+K).GT.0) COLSIZ = COLSIZ + ABS(A(K,J))
  200          CONTINUE
               WORK(J) = COLSIZ
            END IF
  220    CONTINUE
C
C        Find the  NARTIF  smallest such columns.
C        This is an expensive loop.  Later we can replace it by a
C        4-pass process (say), accepting the first col that is within
C        t  of  COLMIN, where  t = 0.0, 0.001, 0.01, 0.1 (say).
C        (This comment written in 1980).
C
C        +       WHILE (NFIXED + NACTIV .LT. N) DO
  240    IF (NFIXED+NACTIV.LT.N) THEN
            COLMIN = FLMAX
            DO 260 J = 1, N
               IF (ISTATE(J).EQ.0) THEN
                  IF (NCLIN.EQ.0) GO TO 280
                  COLSIZ = WORK(J)
                  IF (COLMIN.GT.COLSIZ) THEN
                     COLMIN = COLSIZ
                     JMIN = J
                  END IF
               END IF
  260       CONTINUE
            J = JMIN
  280       ISTATE(J) = 4
            NARTIF = NARTIF + 1
            NFIXED = NFIXED + 1
            GO TO 240
C           +       END WHILE
         END IF
      END IF
C
      NFREE = N - NFIXED
C
      IF (LSDBG) THEN
         IF (ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99996) NFIXED, NACTIV, NARTIF
            CALL X04BAY(IPRINT,4,REC)
         END IF
         IF (ILSDBG(2).GT.0) THEN
            WRITE (REC,FMT=99997)
            CALL X04BAY(IPRINT,2,REC)
            DO 300 I = 1, N, 5
               WRITE (REC,FMT=99995) (WX(J),J=I,MIN(I+4,N))
               CALL X04BAF(IPRINT,REC(1))
  300       CONTINUE
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCU. (LSCRSH)
C
99999 FORMAT (/' //E04NCU// COLD NCLIN NCTOTL',/' //E04NCU// ',L4,I6,I7)
99998 FORMAT (/' //E04NCU// Variables before crash... ')
99997 FORMAT (/' //E04NCU// Variables after  crash... ')
99996 FORMAT (/' //E04NCU// Working set selected ...             ',
     *       /' //E04NCU// NFIXED NACTIV NARTIF      ',/' //E04NCU// ',
     *       I6,2I7)
99995 FORMAT (5G12.3)
      END
      SUBROUTINE E04NBW(MODE,N,NZ,NFREE,NQ,UNITQ,KX,V,ZY,WRK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-590 (MAR 1988).
C
C     ******************************************************************
C     E04NBW  transforms the vector  v  in various ways using the
C     matrix  Q = ( Z  Y )  defined by the input parameters.
C
C        MODE               result
C        ----               ------
C
C          1                v = Z v
C          2                v = Y v
C          3                v = Q v
C
C     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
C     on output, v  is a full n-vector.
C
C
C          4                v = Z'v
C          5                v = Y'v
C          6                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is ordered as  ( v(free)  v(fixed) ).
C
C          7                v = Y'v
C          8                v = Q'v
C
C     On input,  v  is a full n-vector.
C     On output, v  is as in modes 5 and 6 except that v(fixed) is not
C     set.
C
C     Modes  1, 4, 7 and 8  do not involve  v(fixed).
C     Original F66 version  April 1983.
C     Fortran 77 version written  9-February-1985.
C     Level 2 BLAS added 10-June-1986.
C     This version of E04NBW dated 10-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NFREE, NQ, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  V(N), WRK(N), ZY(NQ,*)
      INTEGER           KX(N)
C     .. Local Scalars ..
      INTEGER           J, J1, J2, K, L, LENV, NFIXED
C     .. External Subroutines ..
      EXTERNAL          F06FBF, DCOPY, DGEMV
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      J1 = 1
      J2 = NFREE
      IF (MODE.EQ.1 .OR. MODE.EQ.4) J2 = NZ
      IF (MODE.EQ.2 .OR. MODE.EQ.5 .OR. MODE.EQ.7) J1 = NZ + 1
      LENV = J2 - J1 + 1
      IF (MODE.LE.3) THEN
C        ===============================================================
C        Mode = 1, 2  or  3.
C        ===============================================================
C
         IF (NFREE.GT.0) CALL F06FBF(NFREE,ZERO,WRK,1)
C
C        Copy  v(fixed)  into the end of  wrk.
C
         IF (MODE.GE.2 .AND. NFIXED.GT.0) CALL DCOPY(NFIXED,V(NFREE+1),
     *       1,WRK(NFREE+1),1)
C
C        Set  WRK  =  relevant part of  ZY * V.
C
         IF (LENV.GT.0) THEN
            IF (UNITQ) THEN
               CALL DCOPY(LENV,V(J1),1,WRK(J1),1)
            ELSE
               CALL DGEMV('N',NFREE,J2-J1+1,ONE,ZY(1,J1),NQ,V(J1),1,ONE,
     *                    WRK,1)
            END IF
         END IF
C
C        Expand  WRK  into  V  as a full n-vector.
C
         CALL F06FBF(N,ZERO,V,1)
         DO 20 K = 1, NFREE
            J = KX(K)
            V(J) = WRK(K)
   20    CONTINUE
C
C        Copy  WRK(fixed)  into the appropriate parts of  V.
C
         IF (MODE.GT.1) THEN
            DO 40 L = 1, NFIXED
               J = KX(NFREE+L)
               V(J) = WRK(NFREE+L)
   40       CONTINUE
         END IF
C
      ELSE
C        ===============================================================
C        Mode = 4, 5, 6, 7  or  8.
C        ===============================================================
C        Put the fixed components of  V  into the end of  WRK.
C
         IF (MODE.EQ.5 .OR. MODE.EQ.6) THEN
            DO 60 L = 1, NFIXED
               J = KX(NFREE+L)
               WRK(NFREE+L) = V(J)
   60       CONTINUE
         END IF
C
C        Put the free  components of  V  into the beginning of  WRK.
C
         IF (NFREE.GT.0) THEN
            DO 80 K = 1, NFREE
               J = KX(K)
               WRK(K) = V(J)
   80       CONTINUE
C
C           Set  V  =  relevant part of  ZY' * WRK.
C
            IF (LENV.GT.0) THEN
               IF (UNITQ) THEN
                  CALL DCOPY(LENV,WRK(J1),1,V(J1),1)
               ELSE
                  CALL DGEMV('T',NFREE,J2-J1+1,ONE,ZY(1,J1),NQ,WRK,1,
     *                       ZERO,V(J1),1)
               END IF
            END IF
         END IF
C
C        Copy the fixed components of  WRK  into the end of  V.
C
         IF (NFIXED.GT.0 .AND. (MODE.EQ.5 .OR. MODE.EQ.6))
     *       CALL DCOPY(NFIXED,WRK(NFREE+1),1,V(NFREE+1),1)
      END IF
C
      RETURN
C
C     End of  E04NBW. (CMQMUL)
C
      END
      SUBROUTINE E04NCZ(PRBTYP,NAMED,NAMES,LINOBJ,UNITQ,INFORM,ITER,
     *                  JINF,NCLIN,NCTOTL,NACTIV,NFREE,NRANK,NZ,NRZ,N,
     *                  LDA,LDR,ISTATE,KACTIV,KX,CTX,SSQ,SSQ1,SUMINF,
     *                  NUMINF,XNORM,BL,BU,A,CLAMDA,AX,FEATOL,R,X,W)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1076 (JUL 1993).
C
C     ******************************************************************
C     E04NCZ  is a subroutine for linearly constrained linear-least
C     squares.  On entry, it is assumed that an initial working set of
C     linear constraints and bounds is available.
C     The arrays ISTATE, KACTIV and KX will have been set accordingly
C     and the arrays T and ZY will contain the TQ factorization of
C     the matrix whose rows are the gradients of the active linear
C     constraints with the columns corresponding to the active bounds
C     removed.  the TQ factorization of the resulting (NACTIV by NFREE)
C     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
C     is reverse-triangular.
C
C     Values of ISTATE(J) for the linear constraints.......
C
C     ISTATE(J)
C     ---------
C          0    constraint J is not in the working set.
C          1    constraint J is in the working set at its lower bound.
C          2    constraint J is in the working set at its upper bound.
C          3    constraint J is in the working set as an equality.
C
C     Constraint J may be violated by as much as FEATOL(J).
C
C     Systems Optimization Laboratory, Stanford University.
C     This version of  E04NCZ  dated 14-Sep-1992.
C
C     Copyright  1984  Stanford University.
C
C     This material may be reproduced by or for the U.S. Government
C     pursuant to the copyright license under DAR clause 7-104.9(a)
C     (1979 Mar).
C
C     This material is based upon work partially supported by the
C     National Science Foundation under grants MCS-7926009 and
C     ECS-8012974; the Department of Energy Contract AM03-76SF00326, PA
C     No. DE-AT03-76ER72018; and the Army Research Office Contract DAA29
C     -79-C-0110.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      INTEGER           MSTALL, MREFN
      PARAMETER         (MSTALL=50,MREFN=1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTX, SSQ, SSQ1, SUMINF, XNORM
      INTEGER           INFORM, ITER, JINF, LDA, LDR, N, NACTIV, NCLIN,
     *                  NCTOTL, NFREE, NRANK, NRZ, NUMINF, NZ
      LOGICAL           LINOBJ, NAMED, UNITQ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  CLAMDA(NCTOTL), FEATOL(NCTOTL), R(LDR,*), W(*),
     *                  X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, DTMAX,
     *                  DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9, TOLACT,
     *                  TOLFEA, TOLRNK
      INTEGER           IDBGLS, IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LDT, LDZY, LENNAM,
     *                  LINES1, LINES2, LPROB, MSGLS, NCOLT, NN, NNCLIN,
     *                  NOUT, NPROB
      LOGICAL           CMDBG, LSDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPSVLS(MXPARM), WMACH(15)
      INTEGER           ICMDBG(LDBG), ILSDBG(LDBG), IPADLS(18),
     *                  IPSVLS(MXPARM), LOCLS(LENLS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSRZZ, ALFA, ALFHIT, ATPHIT, BIGALF, CNORM,
     *                  CONDMX, CONDRZ, CONDT, CTP, DINKY, DRZMAX,
     *                  DRZMIN, ERR1, ERR2, FLMAX, GFNORM, GRZNRM,
     *                  GZNORM, OBJSIZ, PALFA, PNORM, RESNRM, ROWNRM,
     *                  TRULAM, WSSIZE
      INTEGER           IADD, IDBG, IFIX, IREFN, IS, ISDEL, ITMAX, JADD,
     *                  JBIGST, JDEL, JMAX1, JSMLST, JTINY, KBIGST,
     *                  KDEL, KSMLST, LANORM, LAP, LCQ, LGQ, LHZ, LPX,
     *                  LRES, LRES0, LRLAM, LT, LWRK, LWTINF, LZY,
     *                  MSGDBG, MSGLVL, MSGSVD, NGQ, NPHASE, NRES,
     *                  NSTALL, NVIOL
      LOGICAL           CONVRG, CYCLIN, ERROR, FIRSTV, HITCON, HITLOW,
     *                  NEEDFG, OVERFL, PRNT, ROWERR, SINGLR, STALL,
     *                  STATPT, UNBNDD, UNCON, UNITGZ, WEAK
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM)
      INTEGER           IPRMLS(MXPARM)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF
      EXTERNAL          DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          E04NBX, E04NCH, E04NCJ, E04NCK, E04NCL, E04NCP,
     *                  E04NCQ, E04NCR, E04NCT, E04NCV, E04UCG, F06FBF,
     *                  F06FLF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04NC/ILSDBG, LSDBG
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IDBGLS, IPRNT, ISUMRY, ITMAX1,
     *                  ITMAX2, LCRASH, LDBGLS, LPROB, MSGLS, NN,
     *                  NNCLIN, NPROB, IPADLS
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IDBGLS), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (MSGLS,MSGLVL), (IDBGLS,IDBG), (LDBGLS,MSGDBG)
C     .. Save statement ..
      SAVE              /AX02ZA/, /DE04NC/, /EE04NC/
C     .. Executable Statements ..
      FLMAX = WMACH(7)
C
      LANORM = LOCLS(2)
      LAP = LOCLS(3)
      LPX = LOCLS(4)
      LRES = LOCLS(5)
      LRES0 = LOCLS(6)
      LHZ = LOCLS(7)
      LGQ = LOCLS(8)
      LCQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK = LOCLS(14)
C
C     Set up the adresses of the contiguous arrays  ( RES0, RES )
C     and  ( GQ, CQ ).
C
      NRES = 0
      IF (NRANK.GT.0) NRES = 2
      NGQ = 1
      IF (LINOBJ) NGQ = 2
C
C     Initialize.
C
      IREFN = 0
      ITER = 0
      JADD = 0
      JDEL = 0
      NPHASE = 1
      NSTALL = 0
      NUMINF = -1
      NRZ = 0
C
      IF (PRBTYP.EQ.'FP') THEN
         ITMAX = ITMAX2
      ELSE
         ITMAX = ITMAX1
      END IF
C
      ALFA = ZERO
      CONDMX = FLMAX
      DRZMAX = ONE
      DRZMIN = ONE
      SSQ = ZERO
C
      CYCLIN = .FALSE.
      ERROR = .FALSE.
      FIRSTV = .FALSE.
      PRNT = .TRUE.
      NEEDFG = .TRUE.
      STALL = .TRUE.
      UNCON = .FALSE.
      UNITGZ = .TRUE.
      UNBNDD = .FALSE.
C
C     If debug output is required,  print nothing until iteration IDBG.
C
      MSGSVD = MSGLVL
      IF (IDBG.GT.0 .AND. IDBG.LE.ITMAX) THEN
         MSGLVL = 0
      END IF
C
C     =================== start of the main loop =======================
C
C      cyclin = false
C      unbndd = false
C      error  = false
C      k      = 0
C
C      repeat
C            repeat
C                  compute Z'g,  print details of this iteration
C                  stat pt = (Z'g .eq. 0)
C                  if (not stat pt) then
C                     error =  k .ge. itmax
C                     if (not error) then
C                        compute p, alfa
C                        error = unbndd  or  cyclin
C                        if (not error) then
C                           k = k + 1
C                           x = x + alfa p
C                           if (feasible) update Z'g
C                           if necessary, add a constraint
C                        end if
C                     end if
C                  end if
C            until  stat pt  or  error
C
C            compute lam1, lam2, smllst
C            optmul =  smllst .gt. 0
C            if ( not (optmul .or. error) ) then
C                  delete an artificial or regular constraint
C            end if
C      until optmul  or  error
C
C     ==================================================================
C
C     REPEAT
C        REPEAT
   20 IF (NEEDFG) THEN
         IF (NRANK.GT.0) THEN
            RESNRM = DNRM2(NRANK,W(LRES),1)
            SSQ = HALF*(SSQ1**2+RESNRM**2)
         END IF
C
         IF (NUMINF.NE.0) THEN
C
C           Compute the transformed gradient of either the sum of
C           infeasibilities or the objective.  Initialize
C           SINGLR and UNITGZ.
C
            CALL E04NCP(PRBTYP,LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDZY,LDR,NRANK,NZ,NRZ,ISTATE,KX,BIGBND,
     *                  TOLRNK,NUMINF,SUMINF,BL,BU,A,W(LRES),FEATOL,
     *                  W(LGQ),W(LCQ),R,X,W(LWTINF),W(LZY),W(LWRK))
            IF (NUMINF.EQ.0 .AND. PRBTYP.NE.'FP' .AND. NPHASE.EQ.1) THEN
               ITMAX = ITER + ITMAX2
               NPHASE = 2
            END IF
         END IF
      END IF
C
      GZNORM = ZERO
      IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
C
      IF (NRZ.EQ.NZ) THEN
         GRZNRM = GZNORM
      ELSE
         GRZNRM = ZERO
         IF (NRZ.GT.0) GRZNRM = DNRM2(NRZ,W(LGQ),1)
      END IF
C
      GFNORM = GZNORM
      IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),1)
C
C        ------------------------------------------------------------
C        Print the details of this iteration.
C        ------------------------------------------------------------
C        Define small quantities that reflect the magnitude of  x,
C        R  and  the matrix of constraints in the working set.
C        Use the largest and smallest diagonals of  R  to estimate
C        the condition number of  Rz1.
C
      IF (NRZ.EQ.0) THEN
         SINGLR = .FALSE.
      ELSE
         IF (NUMINF.GT.0 .OR. NRZ.GT.NRANK) THEN
            ABSRZZ = ZERO
            SINGLR = .TRUE.
         ELSE
            CALL F06FLF(NRZ,R,LDR+1,DRZMAX,DRZMIN)
            ABSRZZ = ABS(R(NRZ,NRZ))
            ROWNRM = DNRM2(N,R(1,1),LDR)
            SINGLR = ABSRZZ .LE. DRZMAX*TOLRNK .OR. ROWNRM .LE.
     *               TOLRNK .OR. ABS(R(1,1)) .LE. ROWNRM*TOLRNK
         END IF
C
         IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
            WRITE (REC,FMT=99995) SINGLR, ABSRZZ, DRZMAX, DRZMIN
            CALL X04BAY(IPRINT,3,REC)
         END IF
C
      END IF
C
      CONDRZ = F06BLF(DRZMAX,DRZMIN,OVERFL)
C
      CONDT = ONE
      IF (NACTIV.GT.0) CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
C
      IF (PRNT) THEN
         CALL E04NCJ(PRBTYP,ISDEL,ITER,JADD,JDEL,MSGLVL,NACTIV,NFREE,N,
     *               NCLIN,NRANK,LDR,LDT,NZ,NRZ,ISTATE,ALFA,CONDRZ,
     *               CONDT,GFNORM,GRZNRM,NUMINF,SUMINF,CTX,SSQ,AX,R,
     *               W(LT),X,W(LWRK))
C
         JDEL = 0
         JADD = 0
         ALFA = ZERO
      END IF
C
      IF (NUMINF.GT.0) THEN
         DINKY = ZERO
      ELSE
         OBJSIZ = ONE + ABS(SSQ+CTX)
         WSSIZE = ZERO
         IF (NACTIV.GT.0) WSSIZE = DTMAX
         DINKY = EPSPT8*MAX(WSSIZE,OBJSIZ,GFNORM)
         IF (UNCON) THEN
            UNITGZ = GRZNRM .LE. DINKY
         END IF
      END IF
C
      IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
         WRITE (REC,FMT=99996) UNITGZ, IREFN, GRZNRM, DINKY
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
C     If the projected gradient  Z'g  is small and Rz is of full
C     rank, X is a minimum on the working set.  An additional
C     refinement step is allowed to take care of an inaccurate
C     value of DINKY.
C
      STATPT = .NOT. SINGLR .AND. GRZNRM .LE. DINKY .OR. IREFN .GT.
     *         MREFN
C
      IF ( .NOT. STATPT) THEN
C        ---------------------------------------------------------
C        Compute a search direction.
C        ---------------------------------------------------------
         PRNT = .TRUE.
C
         ERROR = ITER .GE. ITMAX
         IF ( .NOT. ERROR) THEN
C
            IREFN = IREFN + 1
            ITER = ITER + 1
C
            IF (ITER.EQ.IDBG) THEN
               LSDBG = .TRUE.
               CMDBG = LSDBG
               MSGLVL = MSGSVD
            END IF
C
            CALL E04NCQ(LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,LDA,
     *                  LDZY,LDR,NRANK,NUMINF,NRZ,KX,CTP,PNORM,A,W(LAP),
     *                  W(LRES),W(LHZ),W(LPX),W(LGQ),W(LCQ),R,W(LZY),
     *                  W(LWRK))
C
C           ------------------------------------------------------
C           Find the constraint we bump into along P.
C           Update X and AX if the step ALFA is nonzero.
C           ------------------------------------------------------
C           ALFHIT is initialized to BIGALF.  If it remains
C           that way after the call to E04UCG, it will be
C           regarded as infinite.
C
            BIGALF = F06BLF(BIGDX,PNORM,OVERFL)
C
            CALL E04UCG(FIRSTV,HITLOW,ISTATE,INFORM,JADD,N,NCTOTL,
     *                  NUMINF,ALFHIT,PALFA,ATPHIT,BIGALF,BIGBND,PNORM,
     *                  W(LANORM),W(LAP),AX,BL,BU,FEATOL,W(LPX),X)
C
C           If  Rz1  is nonsingular,  ALFA = 1.0  will be the
C           step to the least-squares minimizer on the
C           current subspace. If the unit step does not violate
C           the nearest constraint by more than FEATOL,  the
C           constraint is not added to the working set.
C
            HITCON = SINGLR .OR. PALFA .LE. ONE
            UNCON = .NOT. HITCON
C
            IF (HITCON) THEN
               ALFA = ALFHIT
            ELSE
               JADD = 0
               ALFA = ONE
            END IF
C
C           Check for an unbounded solution or negligible step.
C
            UNBNDD = ALFA .GE. BIGALF
            STALL = ABS(ALFA*PNORM) .LE. EPSPT9*XNORM
            IF (STALL) THEN
               NSTALL = NSTALL + 1
               CYCLIN = NSTALL .GT. MSTALL
            ELSE
               NSTALL = 0
            END IF
C
            ERROR = UNBNDD .OR. CYCLIN
            IF ( .NOT. ERROR) THEN
C              ---------------------------------------------------
C              Set X = X + ALFA*P.  Update AX, GQ, RES and CTX.
C              ---------------------------------------------------
               IF (ALFA.NE.ZERO) CALL E04NCL(HITCON,HITLOW,LINOBJ,
     *                                UNITGZ,NCLIN,NRANK,NRZ,N,LDR,JADD,
     *                                NUMINF,ALFA,CTP,CTX,XNORM,W(LAP),
     *                                AX,BL,BU,W(LGQ),W(LHZ),W(LPX),
     *                                W(LRES),R,X,W(LWRK))
C
               IF (HITCON) THEN
C                 ------------------------------------------------
C                 Add a constraint to the working set.
C                 Update the TQ factors of the working set.
C                 Use P as temporary work space.
C                 ------------------------------------------------
C                 Update  ISTATE.
C
                  IF (BL(JADD).EQ.BU(JADD)) THEN
                     ISTATE(JADD) = 3
                  ELSE IF (HITLOW) THEN
                     ISTATE(JADD) = 1
                  ELSE
                     ISTATE(JADD) = 2
                  END IF
                  IADD = JADD - N
                  IF (JADD.LE.N) THEN
C
                     DO 40 IFIX = 1, NFREE
                        IF (KX(IFIX).EQ.JADD) GO TO 60
   40                CONTINUE
                  END IF
   60             CONTINUE
C
                  CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTIV,NZ,
     *                        NFREE,NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,
     *                        KX,CONDMX,A,R,W(LT),W(LRES),W(LGQ),W(LZY),
     *                        W(LWRK),W(LRLAM),W(LPX),MSGLVL)
C
                  NRZ = NRZ - 1
                  NZ = NZ - 1
C
                  IF (JADD.LE.N) THEN
C
C                    A simple bound has been added.
C
                     NFREE = NFREE - 1
                  ELSE
C
C                    A general constraint has been added.
C
                     NACTIV = NACTIV + 1
                     KACTIV(NACTIV) = IADD
                  END IF
C
                  IREFN = 0
               END IF
C
C              ---------------------------------------------------
C              Check the feasibility of constraints with non-
C              negative ISTATE values.  If some violations have
C              occurred.  Refine the current X and set INFORM so
C              that feasibility is checked in E04NCP.
C              ---------------------------------------------------
               CALL E04NCR(N,NCLIN,ISTATE,BIGBND,CNORM,ERR1,JMAX1,NVIOL,
     *                     AX,BL,BU,FEATOL,X,W(LWRK))
C
               IF (ERR1.GT.FEATOL(JMAX1)) THEN
                  CALL E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,
     *                        NRANK,NZ,N,NCTOTL,LDZY,LDA,LDR,LDT,ISTATE,
     *                        KACTIV,KX,JMAX1,ERR2,CTX,XNORM,A,AX,BL,BU,
     *                        W(LCQ),W(LRES),W(LRES0),FEATOL,R,W(LT),X,
     *                        W(LZY),W(LPX),W(LWRK))
C
                  IF (LSDBG .AND. ILSDBG(1).GT.0) THEN
                     WRITE (REC,FMT=99998) ERR1, ERR2
                     CALL X04BAY(IPRINT,2,REC)
                  END IF
                  IF (ROWERR) THEN
                     IF (MSGLVL.GT.0) THEN
                        WRITE (REC,FMT=99997)
                        CALL X04BAF(IPRINT,REC(1))
                     END IF
                     NUMINF = 1
                     ERROR = .TRUE.
                  ELSE
                     NUMINF = -1
                     UNCON = .FALSE.
                     IREFN = 0
                  END IF
               END IF
               NEEDFG = ALFA .NE. ZERO
            END IF
         END IF
      END IF
C
C        UNTIL      STATPT  .OR.  ERROR
      IF ( .NOT. (STATPT .OR. ERROR)) GO TO 20
C
C     ===============================================================
C     Try and find the index JDEL of a constraint to drop from
C     the working set.
C     ===============================================================
      JDEL = 0
C
      IF (NUMINF.EQ.0 .AND. PRBTYP.EQ.'FP') THEN
         IF (N.GT.NZ) CALL F06FBF(N-NZ,(ZERO),W(LRLAM),1)
         JTINY = 0
         JSMLST = 0
         JBIGST = 0
      ELSE
C
         CALL E04NCK(PRBTYP,MSGLVL,N,NACTIV,NFREE,LDA,LDT,NUMINF,NZ,NRZ,
     *               ISTATE,KACTIV,KX,DINKY,JSMLST,KSMLST,JINF,JTINY,
     *               JBIGST,KBIGST,TRULAM,A,W(LANORM),W(LGQ),W(LRLAM),
     *               W(LT),W(LWTINF))
C
      END IF
C
      IF ( .NOT. ERROR) THEN
         IF (JSMLST.GT.0) THEN
C
C           E04NCK found a regular constraint with multiplier less
C           than (-DINKY).
C
            JDEL = JSMLST
            KDEL = KSMLST
            ISDEL = ISTATE(JDEL)
            ISTATE(JDEL) = 0
C
         ELSE IF (JSMLST.LT.0) THEN
C
            JDEL = JSMLST
C
         ELSE IF (NUMINF.GT.0 .AND. JBIGST.GT.0) THEN
C
C           No feasible point exists for the constraints but the
C           sum of the constraint violations may be reduced by
C           moving off constraints with multipliers greater than 1.
C
            JDEL = JBIGST
            KDEL = KBIGST
            ISDEL = ISTATE(JDEL)
            IF (TRULAM.LE.ZERO) IS = -1
            IF (TRULAM.GT.ZERO) IS = -2
            ISTATE(JDEL) = IS
            FIRSTV = .TRUE.
            NUMINF = NUMINF + 1
         END IF
C
         IF (JDEL.NE.0 .AND. SINGLR) THEN
C
C           Cannot delete a constraint when Rz is singular.
C           Probably a weak minimum.
C
            JDEL = 0
         ELSE IF (JDEL.NE.0) THEN
C
C           Constraint JDEL has been deleted.
C           Update the matrix factorizations.
C
            CALL E04NCT(UNITQ,N,NACTIV,NFREE,NRES,NGQ,NZ,NRZ,LDA,LDZY,
     *                  LDR,LDT,NRANK,JDEL,KDEL,KACTIV,KX,A,W(LRES),R,
     *                  W(LT),W(LGQ),W(LZY),W(LWRK),W(LPX))
         END IF
      END IF
C
      IREFN = 0
      CONVRG = JDEL .EQ. 0
C
      PRNT = .FALSE.
      UNCON = .FALSE.
      NEEDFG = .FALSE.
C
C     until       convrg  .or.  error
      IF ( .NOT. (CONVRG .OR. ERROR)) GO TO 20
C
C     .....................End of main loop........................
C
      WEAK = JTINY .GT. 0 .OR. SINGLR
C
      IF (ERROR) THEN
         IF (UNBNDD) THEN
            INFORM = 2
            IF (NUMINF.GT.0) INFORM = 3
         ELSE IF (ITER.GE.ITMAX) THEN
            INFORM = 4
         ELSE IF (CYCLIN) THEN
            INFORM = 5
         END IF
      ELSE IF (CONVRG) THEN
         INFORM = 0
         IF (NUMINF.GT.0) THEN
            INFORM = 3
         ELSE IF (PRBTYP.NE.'FP' .AND. WEAK) THEN
            INFORM = 1
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Set   CLAMDA.  Print the full solution.
C     ------------------------------------------------------------------
      MSGLVL = MSGSVD
      IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999) PRBTYP, ITER
         CALL X04BAY(IPRINT,2,REC)
      END IF
C
      CALL E04NBX(MSGLVL,NFREE,LDA,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,
     *            NACTIV,ISTATE,KACTIV,KX,A,BL,BU,X,CLAMDA,W(LRLAM),X)
C
      RETURN
C
C
C     End of  E04NCZ. (LSCORE)
C
99999 FORMAT (/' Exit from ',A2,' problem after ',I5,' iterations.')
99998 FORMAT (' XXX  Iterative refinement.  Maximum errors before and ',
     *       'after refinement are',/6X,1P,2D14.2)
99997 FORMAT (' XXX  Warning.  Cannot satisfy the constraints to the a',
     *       'ccuracy requested.')
99996 FORMAT (/' //E04NCZ//  UNITGZ IREFN     GRZNRM      DINKY   ',
     *       /' //E04NCZ//  ',L6,I6,1P,2D11.2)
99995 FORMAT (/' //E04NCZ//  SINGLR   ABS(RZZ1)      DRZMAX      DRZMIN'
     *       ,/' //E04NCZ//  ',L6,1P,3D12.4)
      END
      SUBROUTINE E04UCY(INFORM,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,
     *                  LDCJU,N,NCNLN,CONFUN,OBJFUN,NEEDC,BIGBND,EPSRF,
     *                  CDINT,FDINT,FDCHK,FDNORM,OBJF,XNORM,BL,BU,C,C1,
     *                  CJAC,CJACU,CJDX,DX,GRAD,GRADU,HFORWD,HCNTRL,X,
     *                  WRK1,WRK2,W,LENW,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1092 (JUL 1993).
C
C     ******************************************************************
C     E04UCY  performs the following...
C     (1)  Computes the objective and constraint values OBJF and C.
C     (2)  Evaluates the user-provided gradients in CJACU and GRADU.
C     (3)  Counts the missing gradients.
C     (4)  Loads the known gradients into GRAD and CJAC.
C     (5)  Checks that the known gradients are programmed correctly.
C     (6)  Computes the missing gradient elements.
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version written 4-September-1985.
C     This version of E04UCY dated 26-Nov-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  RDUMMY
      PARAMETER         (RDUMMY=-11111.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CDINT, EPSRF, FDCHK, FDINT, FDNORM,
     *                  OBJF, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, LENW, LVLDER, MSGNP, N,
     *                  NCNLN, NFUN, NGRAD, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), C(*), C1(*), CJAC(LDCJ,*),
     *                  CJACU(LDCJU,*), CJDX(*), DX(N), GRAD(N),
     *                  GRADU(N), HCNTRL(*), HFORWD(*), USER(*),
     *                  W(LENW), WRK1(N+NCNLN), WRK2(N+NCNLN), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LFDSET, LINES1, LINES2, LVLDIF,
     *                  LVRFYC, NCDIFF, NFDIFF, NOUT
C     .. Arrays in Common ..
      INTEGER           JVERFY(4)
C     .. Local Scalars ..
      INTEGER           I, INFOG, INFOJ, J, MODE, NCSET
      LOGICAL           CENTRL, NEEDFD
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04UDS, E04XAW, E04XAX, E04XAY, F06DBF,
     *                  F06FBF, F06QFF, F06QHF, X04BAY
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
C     .. Executable Statements ..
C
      INFOG = 0
      INFOJ = 0
      NFDIFF = 0
      NCDIFF = 0
      NCSET = N*NCNLN
C
      IF (NCNLN.GT.0) THEN
C        ===============================================================
C        Compute the constraints and Jacobian matrix.
C        ===============================================================
C        If some derivatives are missing, load the Jacobian with dummy
C        values.  Any elements left unaltered after the call to CONFUN
C        must be estimated.  A record of the missing Jacobian elements
C        is stored in  CJACU.
C
         NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 1
C
         IF (NEEDFD) CALL F06QHF('General',NCNLN,N,RDUMMY,RDUMMY,CJACU,
     *                           LDCJU)
C
         CALL F06DBF(NCNLN,(1),NEEDC,1)
C
         MODE = 2
         CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C,CJACU,NSTATE,IUSER,
     *               USER)
         IF (MODE.LT.0) GO TO 80
C
         CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
C
         IF (NEEDFD) THEN
C
C           Count the number of missing Jacobian elements.
C
            DO 40 J = 1, N
               DO 20 I = 1, NCNLN
                  IF (CJACU(I,J).EQ.RDUMMY) NCDIFF = NCDIFF + 1
   20          CONTINUE
   40       CONTINUE
C
            NCSET = NCSET - NCDIFF
            IF (NSTATE.EQ.1) THEN
               IF (NCDIFF.EQ.0) THEN
                  IF (LVLDER.EQ.0) LVLDER = 2
                  IF (LVLDER.EQ.1) LVLDER = 3
                  IF (MSGNP.GT.0) THEN
                     WRITE (REC,FMT=99999) LVLDER
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               ELSE
                  IF (MSGNP.GT.0) THEN
                     WRITE (REC,FMT=99998) NCSET, N*NCNLN, NCDIFF
                     CALL X04BAY(IPRINT,3,REC)
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     ==================================================================
C     Repeat the procedure above for the objective function.
C     ==================================================================
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2
C
      IF (NEEDFD) CALL F06FBF(N,RDUMMY,GRADU,1)
c                                 output the initial value
      iuser(2) = 1

      MODE = 2
      CALL OBJFUN(MODE,N,X,OBJF,GRADU,NSTATE,IUSER,USER)
      IF (MODE.LT.0) GO TO 80
C
      CALL DCOPY(N,GRADU,1,GRAD,1)
c                                 shut off output for derivative
c                                 evaluation
      iuser(2) = 0

      IF (NEEDFD) THEN
C
C        Count the number of missing gradient elements.
C
         DO 60 J = 1, N
            IF (GRADU(J).EQ.RDUMMY) NFDIFF = NFDIFF + 1
   60    CONTINUE
C
         IF (NSTATE.EQ.1) THEN
            IF (NFDIFF.EQ.0) THEN
               IF (LVLDER.EQ.0) LVLDER = 1
               IF (LVLDER.EQ.2) LVLDER = 3
               IF (MSGNP.GT.0) THEN
                  WRITE (REC,FMT=99997) LVLDER
                  CALL X04BAY(IPRINT,3,REC)
               END IF
            ELSE
               IF (MSGNP.GT.0) THEN
                  WRITE (REC,FMT=99996) N - NFDIFF, N, NFDIFF
                  CALL X04BAY(IPRINT,3,REC)
               END IF
            END IF
         END IF
      END IF
C
      NFUN = NFUN + 1
      NGRAD = NGRAD + 1
C
C     ==================================================================
C     Check whatever gradient elements have been provided.
C     ==================================================================
      IF (LVRFYC.GE.0) THEN
         IF (NCSET.GT.0) THEN
            CALL E04XAW(MODE,LVLDER,MSGNP,NCSET,N,NCNLN,LDCJ,LDCJU,
     *                  BIGBND,EPSRF,EPSPT3,FDCHK,XNORM,CONFUN,NEEDC,BL,
     *                  BU,C,C1,CJAC,CJACU,CJDX,DX,WRK2,X,WRK1,IUSER,
     *                  USER)
            IF (MODE.LT.0) GO TO 80
            INFOJ = MODE
         END IF
C
         IF (NFDIFF.LT.N) THEN
            CALL E04XAX(MODE,MSGNP,N,BIGBND,EPSRF,EPSPT3,FDCHK,OBJF,
     *                  XNORM,OBJFUN,BL,BU,GRAD,GRADU,DX,X,WRK1,IUSER,
     *                  USER)
            IF (MODE.LT.0) GO TO 80
            INFOG = MODE
         END IF
      END IF
C
      NEEDFD = NCDIFF .GT. 0 .OR. NFDIFF .GT. 0
      IF (NEEDFD) THEN
C        ===============================================================
C        Compute the missing gradient elements.
C        ===============================================================
         CALL E04XAY(MODE,MSGNP,LVLDER,N,NCNLN,LDCJ,LDCJU,BIGBND,EPSRF,
     *               FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,C1,CJDX,
     *               CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,DX,IUSER,
     *               USER)
C
         IF (MODE.LT.0) GO TO 80
C
         IF (LFDSET.GT.0) THEN
            CENTRL = LVLDIF .EQ. 2
            CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,NEEDC,BL,BU,C,
     *                  C1,CJDX,CJAC,CJACU,GRAD,GRADU,HFORWD,HCNTRL,X,W,
     *                  LENW,IUSER,USER)
C
            IF (MODE.LT.0) GO TO 80
         END IF
      END IF
C
      INFORM = INFOJ + INFOG
      RETURN
C
C     The user requested termination.
C
   80 INFORM = MODE
      RETURN
C
C
C     End of  E04UCY. (NPCHKD)
C
99999 FORMAT (/' All Jacobian elements have been set.',/' Derivative l',
     *       'evel increased to ',I4)
99998 FORMAT (/' The user sets ',I6,'   out of',I6,'   Jacobian elemen',
     *       'ts.',/' Each iteration, ',I6,'   Jacobian elements will ',
     *       'be estimated numerically.')
99997 FORMAT (/' All objective gradient elements have been set.',/' De',
     *       'rivative level increased to ',I4)
99996 FORMAT (/' The user sets ',I6,'   out of',I6,'   objective gradi',
     *       'ent elements.',/' Each iteration, ',I6,'   gradient elem',
     *       'ents will be estimated numerically.')
      END
      SUBROUTINE E04MFK(MSGLVL,NFREE,LDA,N,NCLIN,NCNLN,NCTOTL,BIGBND,
     *                  NAMED,NAMES,LENNAM,NACTIV,ISTATE,KACTIV,KX,A,BL,
     *                  BU,C,CLAMDA,RLAMDA,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFK   creates the expanded Lagrange multiplier vector CLAMDA.
C     If  MSGLVL .eq. 1  or  MSGLVL .ge. 10,  E04MFK  prints  x,  A*x,
C     c(x),  their bounds, the multipliers, and the residuals (distance
C     to the nearer bound).
C
C     This version of  E04MFK  dated  14-Dec-90.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LCMDBG
      PARAMETER         (LCMDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           LDA, LENNAM, MSGLVL, N, NACTIV, NCLIN, NCNLN,
     *                  NCTOTL, NFREE
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(NCTOTL), BU(NCTOTL), C(*),
     *                  CLAMDA(NCTOTL), RLAMDA(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      INTEGER           ICMDBG(LCMDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RES, RES2, RLAM, V, WLAM
      INTEGER           IP, IS, J, K, NFIXED, NPLIN, NZ
      CHARACTER*2       LS
      CHARACTER*5       ID3
      CHARACTER*8       ID4
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(7)
      CHARACTER*5       ID(3)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Data statements ..
      DATA              ID(1)/'VARBL'/
      DATA              ID(2)/'LNCON'/
      DATA              ID(3)/'NLCON'/
      DATA              LSTATE(1)/'--'/, LSTATE(2)/'++'/
      DATA              LSTATE(3)/'FR'/, LSTATE(4)/'LL'/
      DATA              LSTATE(5)/'UL'/, LSTATE(6)/'EQ'/
      DATA              LSTATE(7)/'TF'/
C     .. Executable Statements ..
C
      NPLIN = N + NCLIN
      NZ = NFREE - NACTIV
C
C     Expand multipliers for bounds, linear and nonlinear constraints
C     into the  CLAMDA  array.
C
      CALL F06FBF(NCTOTL,ZERO,CLAMDA,1)
      NFIXED = N - NFREE
      DO 20 K = 1, NACTIV + NFIXED
         IF (K.LE.NACTIV) THEN
            J = KACTIV(K) + N
            RLAM = RLAMDA(NACTIV-K+1)
         ELSE
            J = KX(NZ+K)
            RLAM = RLAMDA(K)
         END IF
         CLAMDA(J) = RLAM
   20 CONTINUE
C
      IF (MSGLVL.NE.1 .AND. MSGLVL.LT.10) RETURN
C
      WRITE (REC,FMT=99999)
      CALL X04BAY(IPRINT,4,REC)
      ID3 = ID(1)
C
      DO 40 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
         WLAM = CLAMDA(J)
         IS = ISTATE(J)
         LS = LSTATE(IS+3)
         IF (J.LE.N) THEN
C
C           Section 1 -- the variables  x.
C           ------------------------------
            K = J
            V = X(J)
C
         ELSE IF (J.LE.NPLIN) THEN
C
C           Section 2 -- the linear constraints  A*x.
C           -----------------------------------------
            IF (J.EQ.N+1) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(2)
            END IF
C
            K = J - N
            V = DDOT(N,A(K,1),LDA,X,1)
         ELSE
C
C           Section 3 -- the nonlinear constraints  c(x).
C           ---------------------------------------------
C
            IF (J.EQ.NPLIN+1) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(IPRINT,4,REC)
               ID3 = ID(3)
            END IF
C
            K = J - NPLIN
            V = C(K)
         END IF
C
C        Print a line for the j-th variable or constraint.
C        -------------------------------------------------
         RES = V - B1
         RES2 = B2 - V
         IF (ABS(RES).GT.ABS(RES2)) RES = RES2
         IP = 1
         IF (B1.LE.(-BIGBND)) IP = 2
         IF (B2.GE.BIGBND) IP = IP + 2
         IF (NAMED) THEN
C
            ID4 = NAMES(J)
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99996) ID4, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99995) ID4, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99994) ID4, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99993) ID4, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
C
         ELSE
C
            IF (IP.EQ.1) THEN
               WRITE (REC,FMT=99992) ID3, K, LS, V, B1, B2, WLAM, RES
            ELSE IF (IP.EQ.2) THEN
               WRITE (REC,FMT=99991) ID3, K, LS, V, B2, WLAM, RES
            ELSE IF (IP.EQ.3) THEN
               WRITE (REC,FMT=99990) ID3, K, LS, V, B1, WLAM, RES
            ELSE
               WRITE (REC,FMT=99989) ID3, K, LS, V, WLAM, RES
            END IF
            CALL X04BAF(IPRINT,REC(1))
         END IF
   40 CONTINUE
      RETURN
C
C     End of  E04MFK.  (CMPRNT)
C
99999 FORMAT (//1X,'Varbl',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99998 FORMAT (//1X,'L Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99997 FORMAT (//1X,'N Con',1X,'State',5X,'Value',5X,'Lower Bound',3X,
     *       'Upper Bound',4X,'Lagr Mult',3X,'Residual',/)
99996 FORMAT (1X,A4,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99995 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,2G12.4)
99994 FORMAT (1X,A4,4X,A2,1X,1P,2G14.6,5X,'None',5X,1P,2G12.4)
99993 FORMAT (1X,A4,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,2G12.4)
99992 FORMAT (1X,A1,I3,4X,A2,1X,1P,3G14.6,1P,2G12.4)
99991 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',5X,1P,G14.6,1P,
     *       2G12.4)
99990 FORMAT (1X,A1,I3,4X,A2,1X,1P,2G14.6,5X,'None',5X,2G12.4)
99989 FORMAT (1X,A1,I3,4X,A2,1X,1P,G14.6,5X,'None',10X,'None',5X,1P,
     *       2G12.4)
      END
      SUBROUTINE E04NCX(UNITQ,INFORM,NZ,NFREE,NRANK,NRES,NGQ,N,LDZY,LDA,
     *                  LDR,LDT,ISTATE,KX,CONDMX,A,R,T,RES,GQ,ZY,W,C,S,
     *                  MSGLVL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1074 (JUL 1993).
C
C     ******************************************************************
C     E04NCX updates the factor R as KX is reordered to reflect the
C     status of the bound constraints given by ISTATE.  KX is reordered
C     so that the fixed variables come last.  One of two alternative
C     are used to reorder KX. One method needs fewer accesses to KX, the
C     other gives a matrix Rz with more rows and columns.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written  30-December-1985.
C     Level-2 matrix routines added 18-May-1988.
C     This version dated 27-April-1993.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDMX
      INTEGER           INFORM, LDA, LDR, LDT, LDZY, MSGLVL, N, NFREE,
     *                  NGQ, NRANK, NRES, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(N), GQ(N,*), R(LDR,*), RES(N,*),
     *                  S(N), T(LDT,*), W(N), ZY(LDZY,*)
      INTEGER           ISTATE(*), KX(N)
C     .. Local Scalars ..
      INTEGER           IADD, IFIX, J, J2, JADD, K, L, LSTART, NACTV,
     *                  NFIXED
C     .. External Subroutines ..
      EXTERNAL          E04NBU, E04NCV
C     .. Executable Statements ..
      NFIXED = N - NFREE
C
      IF (NRANK.LT.N .AND. NRANK.GT.0) THEN
C        ---------------------------------------------------------------
C        R is specified but singular.  Try and keep the dimension of Rz
C        as large as possible.
C        ---------------------------------------------------------------
         NACTV = 0
         NFREE = N
         NZ = N
C
         J = N
C        +       WHILE (J .GT. 0  .AND.  N-NFREE .LT. NFIXED) DO
   20    IF (J.GT.0 .AND. N-NFREE.LT.NFIXED) THEN
            IF (ISTATE(J).GT.0) THEN
               JADD = J
               DO 40 IFIX = NFREE, 1, -1
                  IF (KX(IFIX).EQ.JADD) GO TO 60
   40          CONTINUE
C
C              Add bound JADD.
C
   60          CALL E04NCV(UNITQ,INFORM,IFIX,IADD,JADD,NACTV,NZ,NFREE,
     *                     NRANK,NRES,NGQ,N,LDA,LDZY,LDR,LDT,KX,CONDMX,
     *                     A,R,T,RES,GQ,ZY,W,C,S,MSGLVL)
C
               NFREE = NFREE - 1
               NZ = NZ - 1
            END IF
            J = J - 1
            GO TO 20
C           +       END WHILE
         END IF
      ELSE
C        ---------------------------------------------------------------
C        R is of full rank,  or is not specified.
C        ---------------------------------------------------------------
         IF (NFIXED.GT.0) THEN
C
C           Order KX so that the free variables come first.
C
            LSTART = NFREE + 1
            DO 120 K = 1, NFREE
               J = KX(K)
               IF (ISTATE(J).GT.0) THEN
                  DO 80 L = LSTART, N
                     J2 = KX(L)
                     IF (ISTATE(J2).EQ.0) GO TO 100
   80             CONTINUE
C
  100             KX(K) = J2
                  KX(L) = J
                  LSTART = L + 1
C
                  IF (NRANK.GT.0) CALL E04NBU(N,NRES,NRANK,LDR,K,L,R,
     *                                 RES,C,S)
               END IF
  120       CONTINUE
C
         END IF
         NZ = NFREE
      END IF
C
      RETURN
C
C     End of  E04NCX. (LSBNDS)
C
      END
      SUBROUTINE E04MFR(JOB,MSGLVL,N,NCLIN,NMOVED,ITER,NUMINF,ISTATE,
     *                  BIGBND,AX,BL,BU,FEATOL,FEATLU,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFR performs most of the manoeuvres associated with degeneracy.
C     the degeneracy-resolving strategy operates in the following way.
C
C     Over a cycle of iterations, the feasibility tolerance FEATOL
C     increases slightly (from TOLX0 to TOLX1 in steps of TOLINC).
C     This ensures that all steps taken will be positive.
C
C     After KDEGEN consecutive iterations, variables within
C     FEATOL of their bounds are set exactly on their bounds and x is
C     recomputed to satisfy the general constraints in the working set.
C     FEATOL is then reduced to TOLX0 for the next cycle of iterations.
C
C     FEATLU  is the array of user-supplied feasibility tolerances.
C     FEATOL  is the array of current feasibility tolerances.
C
C     If JOB = 'I', E04MFR initializes the parameters in
C     common block CE04MF:
C
C     TOLX0   is the minimum (scaled) feasibility tolerance.
C     TOLX1   is the maximum (scaled) feasibility tolerance.
C     TOLINC  is the scaled increment to the current FEATOL.
C     KDEGEN  is the expand frequency (specified by the user).
C             it is the frequency of resetting FEATOL to (scaled) TOLX0.
C     NDEGEN  counts the number of degenerate steps (incremented
C             by E04MFS).
C     ITNFIX  is the last iteration at which a JOB = 'E' or 'O' entry
C             caused an x to be put on a constraint.
C     NFIX(j) counts the number of times a JOB = 'O' entry has
C             caused the variables to be placed on the working set,
C             where j=1 if infeasible, j=2 if feasible.
C
C     TOLX0*FEATLU and TOLX1*FEATLU are both close to the feasibility
C     Tolerance FEATLU specified by the user.  (They must both be less
C     than FEATLU.)
C
C
C     If JOB = 'E',  E04MFR has been called after a cycle of KDEGEN
C     iterations.  Constraints in the working set are examined to see if
C     any are off their bounds by an amount approaching FEATOL.  NMOVED
C     returns how many.  If NMOVED is positive,  x  is moved onto the
C     constraints in the working set.  It is assumed that the calling
C     routine will then continue iterations.
C
C
C     If JOB = 'O',  E04MFR is being called after a subproblem has been
C     judged optimal, infeasible or unbounded.  Constraint violations
C     are examined as above.
C
C     E04MFR is based on MINOS 5.2 routine M5DGEN,
C
C     Original version dated  19-April 1988.
C     This version of  E04MFR  dated 13-Dec-1990.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT6
      PARAMETER         (ZERO=0.0D+0,POINT6=0.6D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           ITER, MSGLVL, N, NCLIN, NMOVED, NUMINF
      CHARACTER         JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATLU(N+NCLIN), FEATOL(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      DOUBLE PRECISION  TOLINC, TOLX0
      INTEGER           IPRINT, ISUMM, ITNFIX, KDEGEN, LINES1, LINES2,
     *                  NDEGEN, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, EPSMCH, TOLX1, TOLZ
      INTEGER           IS, J, MAXFIX
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04MF/TOLX0, TOLINC, KDEGEN, NDEGEN, ITNFIX,
     *                  NFIX
C     .. Save statement ..
      SAVE              /AX02ZA/, TOLZ
C     .. Executable Statements ..
C
      NMOVED = 0
      IF (JOB.EQ.'i' .OR. JOB.EQ.'I') THEN
C        ---------------------------------------------------------------
C        Job = 'Initialize'.
C        Initialize at the start of each linear problem.
C        KDEGEN  is the expand frequency      and
C        FEATLU  are the user-supplied feasibility tolerances.
C        They are not changed.
C        ---------------------------------------------------------------
         EPSMCH = WMACH(3)
C
         NDEGEN = 0
         ITNFIX = 0
         NFIX(1) = 0
         NFIX(2) = 0
         TOLX0 = 0.5D+0
         TOLX1 = 0.99D+0
         TOLZ = EPSMCH**POINT6
C
         IF (KDEGEN.LT.9999999) THEN
            TOLINC = (TOLX1-TOLX0)/KDEGEN
         ELSE
            TOLINC = ZERO
         END IF
C
         DO 20 J = 1, N + NCLIN
            FEATOL(J) = TOLX0*FEATLU(J)
   20    CONTINUE
      ELSE
C        ---------------------------------------------------------------
C        JOB = 'End of cycle' or 'Optimal'.
C        initialize local variables MAXFIX and TOLZ.
C        ---------------------------------------------------------------
         MAXFIX = 2
C
         IF (JOB.EQ.'o' .OR. JOB.EQ.'O') THEN
C           ------------------------------------------------------------
C           JOB = 'Optimal'.
C           return with NMOVED = 0 if the last call was at the same
C           iteration,  or if there have already been MAXFIX calls with
C           the same state of feasibility.
C           ------------------------------------------------------------
            IF (ITNFIX.EQ.ITER) RETURN
            IF (NUMINF.GT.0) THEN
               J = 1
            ELSE
               J = 2
            END IF
C
            IF (NFIX(J).GE.MAXFIX) RETURN
            NFIX(J) = NFIX(J) + 1
         END IF
C
C        Reset FEATOL to its minimum value.
C
         DO 40 J = 1, N + NCLIN
            FEATOL(J) = TOLX0*FEATLU(J)
   40    CONTINUE
C
C        Count the number of times a variable is moved a nontrivial
C        distance onto its bound.
C
         ITNFIX = ITER
C
         DO 60 J = 1, N
            IS = ISTATE(J)
            IF (IS.GT.0 .AND. IS.LT.4) THEN
               IF (IS.EQ.1) THEN
                  D = ABS(X(J)-BL(J))
               ELSE
                  D = ABS(X(J)-BU(J))
               END IF
C
               IF (D.GT.TOLZ) NMOVED = NMOVED + 1
            END IF
   60    CONTINUE
C
         IF (NMOVED.GT.0) THEN
C
C           Some variables were moved onto their bounds.
C
            IF (MSGLVL.GT.0) THEN
               WRITE (REC,FMT=99999) ITER, NMOVED
               CALL X04BAF(IPRINT,REC)
            END IF
         END IF
      END IF
C
C     End of E04MFR.  (CMDGEN)
C
99999 FORMAT (' Itn',I6,' --',I7,'  variables moved to their bounds.')
      END
      SUBROUTINE E04NBZ(NERROR,MSGLVL,LCRASH,LIWORK,LWORK,LITOTL,LWTOTL,
     *                  N,NCLIN,NCNLN,ISTATE,KX,NAMED,NAMES,BIGBND,BL,
     *                  BU,X,M,LDA,LDR,LDCJ,LDFJ,NERR,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-539 (FEB 1987).
C     MARK 13 REVISED. IER-592 (MAR 1988).
C     MARK 16 REVISED. IER-1059 (JUL 1993).
C
C     ******************************************************************
C     E04NBZ   checks the input data for E04UCF and E04UPF.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written 10-May-1980.
C     Fortran 77 version written  5-October-1984.
C     This version of E04NBZ dated  21-Mar-93.
C     ******************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND
      INTEGER           IFAIL, LCRASH, LDA, LDCJ, LDFJ, LDR, LITOTL,
     *                  LIWORK, LWORK, LWTOTL, M, MSGLVL, N, NCLIN,
     *                  NCNLN, NERR, NERROR
      LOGICAL           NAMED
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN), X(N)
      INTEGER           ISTATE(N+NCLIN+NCNLN), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2
      INTEGER           IS, J, K, L
      LOGICAL           OK
C     .. Local Arrays ..
      CHARACTER*5       ID(3)
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Data statements ..
      DATA              ID(1), ID(2), ID(3)/'Varbl', 'L Con', 'N Con'/
C     .. Executable Statements ..
C
      NERROR = 0
C
C     ------------------------------------------------------------------
C     Check  M.
C     ------------------------------------------------------------------
      IF (M.LE.0) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99997) M
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  N.
C     ------------------------------------------------------------------
      IF (N.LE.0) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99996) N
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  NCLIN and NCNLN.
C     ------------------------------------------------------------------
      IF (NCLIN.LT.0 .OR. NCNLN.LT.0) THEN
         IF (NCLIN.LT.0) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99995) NCLIN
               CALL X04BAY(NERR,3,REC)
            END IF
         END IF
C
         IF (NCNLN.LT.0) THEN
            NERROR = NERROR + 1
            IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
               WRITE (REC,FMT=99994) NCNLN
               CALL X04BAY(NERR,3,REC)
            END IF
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDA.
C     ------------------------------------------------------------------
      IF (LDA.LT.MAX(1,NCLIN)) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99993) LDA, NCLIN
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDCJ.
C     ------------------------------------------------------------------
      IF (LDCJ.LT.MAX(1,NCNLN)) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99992) LDCJ, NCNLN
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDFJ.
C     ------------------------------------------------------------------
      IF (LDFJ.LT.M) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99991) LDFJ, M
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check  LDR.
C     ------------------------------------------------------------------
      IF (LDR.LT.N) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99990) LDR, N
            CALL X04BAY(NERR,3,REC)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Check if there is enough workspace to solve the problem.
C     ------------------------------------------------------------------
      OK = LITOTL .LE. LIWORK .AND. LWTOTL .LE. LWORK
      IF ( .NOT. OK) THEN
         NERROR = NERROR + 1
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            WRITE (REC,FMT=99999) LIWORK, LWORK, LITOTL, LWTOTL
            CALL X04BAY(NERR,3,REC)
            WRITE (REC,FMT=99998)
            CALL X04BAY(NERR,2,REC)
         END IF
      ELSE IF (MSGLVL.GT.0) THEN
         WRITE (REC,FMT=99999) LIWORK, LWORK, LITOTL, LWTOTL
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      IF (NERROR.EQ.0) THEN
C
C        ---------------------------------------------------------------
C        Check the bounds on all variables and constraints.
C        ---------------------------------------------------------------
         DO 20 J = 1, N + NCLIN + NCNLN
            B1 = BL(J)
            B2 = BU(J)
            OK = B1 .LT. B2 .OR. (B1.EQ.B2 .AND. ABS(B1).LT.BIGBND)
            IF ( .NOT. OK) THEN
               NERROR = NERROR + 1
               IF (J.GT.N+NCLIN) THEN
                  K = J - N - NCLIN
                  L = 3
               ELSE IF (J.GT.N) THEN
                  K = J - N
                  L = 2
               ELSE
                  K = J
                  L = 1
               END IF
               IF (NAMED) THEN
                  IF (B1.EQ.B2) THEN
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99989) NAMES(J), J, J, B1, BIGBND
                        CALL X04BAY(NERR,4,REC)
                     END IF
                  ELSE
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99988) NAMES(J), J, B1, J, B2
                        CALL X04BAY(NERR,3,REC)
                     END IF
                  END IF
               ELSE
                  IF (B1.EQ.B2) THEN
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99987) ID(L), K, J, J, B1, BIGBND
                        CALL X04BAY(NERR,4,REC)
                     END IF
                  ELSE
                     IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                        WRITE (REC,FMT=99986) ID(L), K, J, B1, J, B2
                        CALL X04BAY(NERR,3,REC)
                     END IF
                  END IF
               END IF
            END IF
   20    CONTINUE
C
C        ---------------------------------------------------------------
C        If warm start, check  ISTATE.
C        ---------------------------------------------------------------
         IF (LCRASH.EQ.1) THEN
            DO 40 J = 1, N + NCLIN + NCNLN
               IS = ISTATE(J)
               OK = IS .GE. (-2) .AND. IS .LE. 4
               IF ( .NOT. OK) THEN
                  NERROR = NERROR + 1
                  IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
                     WRITE (REC,FMT=99985) J, J, IS
                     CALL X04BAY(NERR,3,REC)
                  END IF
               END IF
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NBZ. (CMCHK)
C
99999 FORMAT (/' Workspace provided is     IWORK(',I8,'),  WORK(',I8,
     *       ').',/' To solve problem we need  IWORK(',I8,'),  WORK(',
     *       I8,').')
99998 FORMAT (/' ** Not enough workspace to solve problem.')
99997 FORMAT (/' ** On entry, M.le.0:',/'    M = ',I16)
99996 FORMAT (/' ** On entry, N.le.0:',/'    N = ',I16)
99995 FORMAT (/' ** On entry, NCLIN.lt.0:',/'    NCLIN = ',I16)
99994 FORMAT (/' ** On entry, NCNLN.lt.0:',/'    NCNLN = ',I16)
99993 FORMAT (/' ** On entry, LDA.lt.max(1,NCLIN):',/'    LDA = ',I16,
     *       '   NCLIN = ',I16)
99992 FORMAT (/' ** On entry, LDCJ.lt.max(1,NCNLN):',/'    LDCJ = ',I16,
     *       '   NCNLN = ',I16)
99991 FORMAT (/' ** On entry, LDFJ.lt.M:',/'    LDFJ = ',I16,'   M = ',
     *       I16)
99990 FORMAT (/' ** On entry, LDR.lt.N:',/'    LDR = ',I16,'   N = ',
     *       I16)
99989 FORMAT (/' ** On entry, the equal bounds on  ',A8,'  are infinit',
     *       'e (because',/'    BL(',I4,').eq.beta and BU(',I4,').eq.b',
     *       'eta, but abs(beta).ge.bigbnd):',/'    beta =',G16.7,' bi',
     *       'gbnd =',G16.7)
99988 FORMAT (/' ** On entry, the bounds on  ',A8,'  are inconsistent:',
     *       /'    BL(',I4,') =',G16.7,'   BU(',I4,') =',G16.7)
99987 FORMAT (/' ** On entry, the equal bounds on  ',A5,I3,'  are infi',
     *       'nite (because',/'    BL(',I4,').eq.beta and BU(',I4,').e',
     *       'q.beta, but abs(beta).ge.bigbnd):',/'    beta =',G16.7,
     *       ' bigbnd =',G16.7)
99986 FORMAT (/' ** On entry, the bounds on  ',A5,I3,'  are inconsiste',
     *       'nt:',/'    BL(',I4,') =',G16.7,'   BU(',I4,') =',G16.7)
99985 FORMAT (/' ** On entry with a Warm Start, ISTATE(',I4,') is out ',
     *       'of range:',/'    ISTATE(',I4,') = ',I16)
      END
      SUBROUTINE E04UCS(COLD,N,NCLIN,NCNLN,NCTOTL,NACTIV,NFREE,NZ,
     *                  ISTATE,KACTIV,BIGBND,TOLACT,BL,BU,C)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1087 (JUL 1993).
C
C     ******************************************************************
C     E04UCS  adds indices of nonlinear constraints to the initial
C     working set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   14-February 1985.
C     This version of  E04UCS  dated 14-November-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLACT
      INTEGER           N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE, NZ
      LOGICAL           COLD
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(NCTOTL), BU(NCTOTL), C(*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, BIGLOW, BIGUPP, CMIN, RES, RESL, RESU,
     *                  TOOBIG
      INTEGER           I, IMIN, IS, J, LINACT, NFIXED, NLNACT, NPLIN
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Subroutines ..
      EXTERNAL          X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      LINACT = NACTIV
      NPLIN = N + NCLIN
C
C     If a cold start is being made, initialize the status of the QP
C     working set.  First,  if  BL(j) = BU(j),  set ISTATE(j)=3.
C
      IF (COLD) THEN
         DO 20 J = NPLIN + 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J).EQ.BU(J)) ISTATE(J) = 3
   20    CONTINUE
      END IF
C
C     Increment NACTIV and KACTIV.
C     Ensure that the number of bounds and general constraints in the
C     QP  working set does not exceed N.
C
      DO 40 J = NPLIN + 1, NCTOTL
         IF (NFIXED+NACTIV.EQ.N) ISTATE(J) = 0
         IF (ISTATE(J).GT.0) THEN
            NACTIV = NACTIV + 1
            KACTIV(NACTIV) = J - N
         END IF
   40 CONTINUE
C
      IF (COLD) THEN
C
C        ---------------------------------------------------------------
C        If a cold start is required, an attempt is made to add as many
C        nonlinear constraints as possible to the working set.
C        ---------------------------------------------------------------
C        The following loop finds the most violated constraint.  If
C        there is room in KACTIV, it will be added to the working set
C        and the process will be repeated.
C
C
         IS = 1
         BIGLOW = -BIGBND
         BIGUPP = BIGBND
         TOOBIG = TOLACT + TOLACT
C
C        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
   60    IF (IS.GT.0 .AND. NFIXED+NACTIV.LT.N) THEN
            IS = 0
            CMIN = TOLACT
C
            DO 80 I = 1, NCNLN
               J = NPLIN + I
               IF (ISTATE(J).EQ.0) THEN
                  B1 = BL(J)
                  B2 = BU(J)
                  RESL = TOOBIG
                  RESU = TOOBIG
                  IF (B1.GT.BIGLOW) RESL = ABS(C(I)-B1)/(ONE+ABS(B1))
                  IF (B2.LT.BIGUPP) RESU = ABS(C(I)-B2)/(ONE+ABS(B2))
                  RES = MIN(RESL,RESU)
                  IF (RES.LT.CMIN) THEN
                     CMIN = RES
                     IMIN = I
                     IS = 1
                     IF (RESL.GT.RESU) IS = 2
                  END IF
               END IF
   80       CONTINUE
C
            IF (IS.GT.0) THEN
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = NCLIN + IMIN
               J = NPLIN + IMIN
               ISTATE(J) = IS
            END IF
            GO TO 60
C           end while
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     An initial working set has now been selected.
C     ------------------------------------------------------------------
      NLNACT = NACTIV - LINACT
      NZ = NFREE - NACTIV
      IF (NPDBG .AND. INPDBG(1).GT.0) THEN
         WRITE (REC,FMT=99999) NFIXED, LINACT, NLNACT
         CALL X04BAY(IPRINT,4,REC)
      END IF
C
      RETURN
C
C
C     End of  E04UCS. (NPCRSH)
C
99999 FORMAT (/' //E04UCS//  Working set selected....',/' //E04UCS// N',
     *       'FIXED LINACT NLNACT     ',/' //E04UCS//',3I7)
      END
      SUBROUTINE E04MFU(PRBTYP,HEADER,RSET,MSGLVL,ITER,ISDEL,JDEL,JADD,
     *                  N,NCLIN,NACTIV,NFREE,NZ,NRZ,LDR,LDT,ISTATE,ALFA,
     *                  CONDRZ,CONDT,DRZZ,GFNORM,GRZNRM,NUMINF,SUMINF,
     *                  NOTOPT,OBJLP,TRUSML,AX,R,T,X,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ==================================================================
C     E04MFU  prints various levels of output for E04MFZ.
C
C           msg        cumulative result
C           ---        -----------------
C
C       .le.  0        no output.
C
C       .eq.  1        nothing now (but full output later).
C
C       .eq.  5        one terse line of output.
C
C       .ge. 10        same as 5 (but full output later).
C
C       .ge. 20        constraint status,  x  and  Ax.
C
C       .ge. 30        diagonals of  T  and  R.
C
C
C     Debug printing is controlled by the logical variable  LCDBG.
C     LCDBG  is set true after  IDBG  iterations.  Then, the amount of
C     output is determined by a string of binary digits of the form
C     SVT  (stored in the integer array  ILCDBG).
C
C     S  set 'on'  gives information from the max step routine E04MFS.
C     V  set 'on'  gives various vectors in E04MFZ  and its auxiliaries.
C     T  set 'on'  gives a trace of which routine was called and an
C                  indication of the progress of the run.
C
C     Based on a version of lcprt, written by PEG, 16-February-1987.
C     This version of  E04MFU  dated  11-Nov-92.
C     ==================================================================
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      INTEGER           MLINE1, MLINE2
      PARAMETER         (MLINE1=50000,MLINE2=50000)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, CONDRZ, CONDT, DRZZ, GFNORM, GRZNRM,
     *                  OBJLP, SUMINF, TRUSML
      INTEGER           ISDEL, ITER, JADD, JDEL, LDR, LDT, MSGLVL, N,
     *                  NACTIV, NCLIN, NFREE, NOTOPT, NRZ, NUMINF, NZ
      LOGICAL           HEADER, RSET
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), R(LDR,*), T(LDT,*), WORK(N), X(N)
      INTEGER           ISTATE(*)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           LCDBG
C     .. Arrays in Common ..
      INTEGER           ILCDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  OBJ
      INTEGER           ITN, J, JJ, K, KADD, KDEL, KK, NDF
      LOGICAL           NEWSET, PRTHDR
      CHARACTER*2       LADD, LDEL
      CHARACTER*15      LMCHAR
C     .. Local Arrays ..
      CHARACTER*2       LSTATE(0:5)
      CHARACTER*120     REC(5)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, MOD
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /EE04MF/ILCDBG, LCDBG
C     .. Data statements ..
      DATA              LSTATE(0), LSTATE(1), LSTATE(2)/'  ', 'L ',
     *                  'U '/
      DATA              LSTATE(3), LSTATE(4), LSTATE(5)/'E ', 'F ',
     *                  'A '/
C     .. Executable Statements ..
C
      IF (MSGLVL.GE.15) THEN
         IF (ISUMM.GE.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
            CALL X04BAY(ISUMM,5,REC)
         END IF
      END IF
C
      IF (MSGLVL.GE.5) THEN
C        ---------------------------------------------------------------
C        Some printing required.  Set up information for the terse line.
C        ---------------------------------------------------------------
         ITN = MOD(ITER,10000)
         NDF = MOD(NRZ,10000)
C
         IF (JDEL.NE.0) THEN
            IF (NOTOPT.GT.0) THEN
               WRITE (LMCHAR,FMT='( I5, 1P,D10.2 )') NOTOPT, TRUSML
            ELSE
               WRITE (LMCHAR,FMT='( 5X, 1P,D10.2 )') TRUSML
            END IF
C
            IF (JDEL.GT.0) THEN
               KDEL = ISDEL
C
            ELSE IF (JDEL.LT.0) THEN
               JDEL = NZ - NRZ + 1
               KDEL = 5
            END IF
         ELSE
            JDEL = 0
            KDEL = 0
            LMCHAR = '               '
         END IF
C
         IF (JADD.GT.0) THEN
            KADD = ISTATE(JADD)
         ELSE
            KADD = 0
         END IF
C
         LDEL = LSTATE(KDEL)
         LADD = LSTATE(KADD)
C
         IF (NUMINF.GT.0) THEN
            OBJ = SUMINF
         ELSE
            OBJ = OBJLP
         END IF
C
C        ---------------------------------------------------------------
C        If necessary, print a header.
C        Print a single line of information.
C        ---------------------------------------------------------------
         IF (ISUMM.GE.0) THEN
C           -----------------------------------
C           Terse line for the Monitoring file.
C           -----------------------------------
            NEWSET = LINES1 .GE. MLINE1
            PRTHDR = MSGLVL .GE. 15 .OR. HEADER .OR. NEWSET
C
            IF (PRTHDR) THEN
               WRITE (REC,FMT=99997)
               CALL X04BAY(ISUMM,3,REC)
               LINES1 = 0
            END IF
C
            WRITE (REC,FMT=99995) ITN, JDEL, LDEL, JADD, LADD, ALFA,
     *        NUMINF, OBJ, N - NFREE, NACTIV, NZ - NRZ, NDF, GRZNRM,
     *        LMCHAR, CONDT
            CALL X04BAF(ISUMM,REC(1))
            LINES1 = LINES1 + 1
         END IF
C
         IF (IPRINT.GE.0 .AND. ISUMM.NE.IPRINT) THEN
C           ------------------------------
C           Terse line for the Print file.
C           ------------------------------
            NEWSET = LINES2 .GE. MLINE2
            PRTHDR = HEADER .OR. NEWSET
C
            IF (PRTHDR) THEN
               WRITE (REC,FMT=99998)
               CALL X04BAY(IPRINT,3,REC)
               LINES2 = 0
            END IF
C
            WRITE (REC,FMT=99996) ITN, ALFA, NUMINF, OBJ, GRZNRM
            CALL X04BAF(IPRINT,REC(1))
            LINES2 = LINES2 + 1
         END IF
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99994) PRBTYP
               CALL X04BAY(ISUMM,3,REC)
               WRITE (REC,FMT=99993)
               CALL X04BAY(ISUMM,2,REC)
               DO 20 J = 1, N, 5
                  WRITE (REC,FMT=99992) (X(JJ),ISTATE(JJ),JJ=J,
     *              MIN(J+4,N))
                  CALL X04BAF(ISUMM,REC(1))
   20          CONTINUE
               IF (NCLIN.GT.0) THEN
                  WRITE (REC,FMT=99991)
                  CALL X04BAY(ISUMM,2,REC)
                  DO 40 K = 1, NCLIN, 5
                     WRITE (REC,FMT=99992) (AX(KK),ISTATE(N+KK),KK=K,
     *                 MIN(K+4,NCLIN))
                     CALL X04BAF(ISUMM,REC(1))
   40             CONTINUE
               END IF
C
               IF (MSGLVL.GE.30) THEN
C                 ------------------------------------------------------
C                 Print the diagonals of  T  and  R.
C                 ------------------------------------------------------
                  IF (NACTIV.GT.0) THEN
                     CALL DCOPY(NACTIV,T(1,NZ+1),LDT+1,WORK,1)
                     WRITE (REC,FMT=99990) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 60 J = 1, NACTIV, 5
                        WRITE (REC,FMT=99989) (WORK(JJ),JJ=J,
     *                    MIN(J+4,NACTIV))
                        CALL X04BAF(ISUMM,REC(1))
   60                CONTINUE
                  END IF
                  IF (RSET .AND. NRZ.GT.0) THEN
                     WRITE (REC,FMT=99988) PRBTYP
                     CALL X04BAY(ISUMM,2,REC)
                     DO 80 J = 1, NRZ, 5
                        WRITE (REC,FMT=99989) (R(JJ,JJ),JJ=J,
     *                    MIN(J+4,NRZ))
                        CALL X04BAF(ISUMM,REC(1))
   80                CONTINUE
                  END IF
               END IF
               WRITE (REC,FMT=99987)
               CALL X04BAY(ISUMM,4,REC)
            END IF
         END IF
      END IF
C
      HEADER = .FALSE.
      JDEL = 0
      JADD = 0
      ALFA = ZERO
C
C
C     End of E04MFU.  (LPPRNT)
C
      RETURN
C
99999 FORMAT (///' ',A2,' iteration',I5,/' =================')
99998 FORMAT (//'  Itn     Step Ninf Sinf/Objective  Norm Gz')
99997 FORMAT (//'  Itn Jdel  Jadd      Step Ninf  Sinf/Objective  Bnd ',
     *       ' Lin  Art   Zr  Norm Gz NOpt    Min Lm   Cond T')
99996 FORMAT (I5,1P,D9.1,I5,D15.6,D9.1)
99995 FORMAT (I5,I5,A1,I5,A1,1P,D9.1,I5,D16.8,4I5,D9.1,A15,D9.1)
99994 FORMAT (/' Values and status of the ',A2,' constraints',/' -----',
     *       '----------------------------------')
99993 FORMAT (/' Variables...')
99992 FORMAT (1X,5(1P,D15.6,I5))
99991 FORMAT (/' General linear constraints...')
99990 FORMAT (/' Diagonals of ',A2,' working set factor T')
99989 FORMAT (1P,5D15.6)
99988 FORMAT (/' Diagonals of ',A2,' triangle Rz        ')
99987 FORMAT (///' ---------------------------------------------------',
     *       '--------------------------------------------')
      END
      SUBROUTINE E04UDR(UNITQ,N,NFREE,NZ,NQ,NROWR,IPERM,KX,GQ,R,ZY,WORK,
     *                  QRWORK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-793 (DEC 1989).
C     MARK 16 REVISED. IER-1094 (JUL 1993).
C
C     ******************************************************************
C     E04UDR  bounds the condition estimator of the transformed Hessian.
C     On exit, R is of the form
C                  ( DRz   0     )
C                  (  0  sigma*I )
C     where D is a diagonal matrix such that DRz has a bounded condition
C     number,  I is the identity matrix and sigma  is the geometric mean
C     of the largest and smallest elements of DRz. The QR factorization
C     with interchanges is used to give diagonals of DRz that are
C     decreasing in modulus.
C
C     Systems Optimization Laboratory, Stanford University.
C
C     Original version of E04UDR dated  4-August-1986.
C     This version dated  14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           N, NFREE, NQ, NROWR, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  GQ(N), QRWORK(2*N), R(NROWR,*), WORK(N),
     *                  ZY(NQ,*)
      INTEGER           IPERM(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DRMAX, DRMIN, RCNDBD, RFROBN
      LOGICAL           NPDBG
C     .. Arrays in Common ..
      INTEGER           INPDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  DRGM, DRGS, GJMAX, SCLE, SUMSQ
      INTEGER           INFO, J, JMAX, JSAVE, NRANK
C     .. External Functions ..
      DOUBLE PRECISION  F06BMF
      INTEGER           F06KLF
      EXTERNAL          F06BMF, F06KLF
C     .. External Subroutines ..
      EXTERNAL          DSWAP, F01QFF, F06FBF, F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Common blocks ..
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /FE04UC/INPDBG, NPDBG
C     .. Executable Statements ..
C
C     ==================================================================
C     Bound the condition estimator of Q'HQ.
C     ==================================================================
      IF (NZ.GT.1) THEN
C        ---------------------------------------------------------------
C        Refactorize Rz.  Interchanges are used to give diagonals
C        of decreasing magnitude.
C        ---------------------------------------------------------------
         DO 20 J = 1, NZ - 1
            CALL F06FBF(NZ-J,ZERO,R(J+1,J),1)
   20    CONTINUE
C
         CALL F01QFF('Column iterchanges',NZ,NZ,R,NROWR,WORK,IPERM,
     *               QRWORK,INFO)
C
         DO 40 J = 1, NZ
            JMAX = IPERM(J)
            IF (JMAX.GT.J) THEN
               IF (UNITQ) THEN
                  JSAVE = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J) = JSAVE
               ELSE
                  CALL DSWAP(NFREE,ZY(1,JMAX),1,ZY(1,J),1)
               END IF
C
               GJMAX = GQ(JMAX)
               GQ(JMAX) = GQ(J)
               GQ(J) = GJMAX
            END IF
   40    CONTINUE
      END IF
C
      DRGM = ONE
C
      IF (NZ.GT.0) THEN
         NRANK = F06KLF(NZ,R,NROWR+1,ONE/RCNDBD)
         DRGM = HALF*SQRT(ABS(R(1,1)*R(NRANK,NRANK)))
         DRGS = ABS(R(1,1))/RCNDBD
C
         IF (NZ.GT.NRANK) THEN
            DO 60 J = NRANK + 1, NZ
               CALL F06FBF(J-1,ZERO,R(1,J),1)
   60       CONTINUE
            CALL F06FBF(NZ-NRANK,DRGS,R(NRANK+1,NRANK+1),NROWR+1)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Reset the range-space partition of the Hessian.
C     ------------------------------------------------------------------
      IF (NZ.LT.N) THEN
         DO 80 J = NZ + 1, N
            CALL F06FBF(J,ZERO,R(1,J),1)
   80    CONTINUE
         CALL F06FBF(N-NZ,DRGM,R(NZ+1,NZ+1),NROWR+1)
      END IF
C
C     Recompute the Frobenius norm of R.
C
      SCLE = SQRT(DBLE(N-NZ))*DRGM
      SUMSQ = ONE
      DO 100 J = 1, NZ
         CALL F06FJF(J,R(1,J),1,SCLE,SUMSQ)
  100 CONTINUE
      RFROBN = F06BMF(SCLE,SUMSQ)
C
      RETURN
C
C     End of  E04UDR. (NPRSET)
C
      END
      SUBROUTINE E04MFT(START,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,
     *                  LDA,ISTATE,KACTIV,KX,BIGBND,TOLACT,A,AX,BL,BU,
     *                  FEATOL,X,WX,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFT  computes the quantities  ISTATE (optionally),  KACTIV,
C     NACTIV,  nz  and  NFREE  associated with the working set at x.
C
C     The computation depends upon the value of the input parameter
C     START,  as follows...
C
C     START = 'COLD'  An initial working set will be selected. First,
C                     nearly-satisfied or violated bounds are added.
C                     Next,  general linear constraints are added that
C                     have small residuals.
C
C     START = 'WARM'  The quantities KACTIV, NACTIV and NFREE are
C                     initialized from ISTATE,  specified by the user.
C
C     If VERTEX is true, an artificial vertex is defined by fixing some
C     variables on their bounds.  Infeasible variables selected for the
C     artificial vertex are fixed at their nearest bound.  Otherwise,
C     the variables are unchanged.
C
C     Values of ISTATE(j)....
C
C        - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     Original version written by  PEG, 31-October-1984.
C     This version of  E04MFT  dated 14-May-1992.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLACT
      INTEGER           LDA, N, NACTIV, NARTIF, NCLIN, NCTOTL, NFREE
      LOGICAL           VERTEX
      CHARACTER*4       START
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), WORK(N), WX(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           CMDBG
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           ICMDBG(LDBG)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, BIGLOW, BIGUPP, COLMIN, COLSIZ, FLMAX,
     *                  RESIDL, RESL, RESMIN, RESU, TOL, TOOBIG
      INTEGER           I, IMIN, IS, J, JFIX, JFREE, JMIN, K
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04MFN, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /FE04NB/ICMDBG, CMDBG
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      FLMAX = WMACH(7)
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ------------------------------------------------------------------
C     Move the variables inside their bounds.
C     ------------------------------------------------------------------
C
      DO 20 J = 1, N
         B1 = BL(J)
         B2 = BU(J)
         TOL = FEATOL(J)
C
         IF (B1.GT.BIGLOW) THEN
            IF (X(J).LT.B1-TOL) X(J) = B1
         END IF
C
         IF (B2.LT.BIGUPP) THEN
            IF (X(J).GT.B2+TOL) X(J) = B2
         END IF
   20 CONTINUE
C
      CALL DCOPY(N,X,1,WX,1)
C
      IF (CMDBG) THEN
         IF (ICMDBG(1).GT.0) THEN
            WRITE (REC,FMT=99999) START, NCLIN, NCTOTL
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (ICMDBG(2).GT.0) CALL E04MFN('E04MFT',
     *                            'variables before crash',WX,N)
      END IF
C
      NFREE = N
      NACTIV = 0
      NARTIF = 0
C
      IF (START.EQ.'COLD') THEN
         DO 40 J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J).EQ.BU(J)) ISTATE(J) = 3
   40    CONTINUE
C
      ELSE IF (START.EQ.'WARM') THEN
         DO 60 J = 1, NCTOTL
            IF (ISTATE(J).GT.3 .OR. ISTATE(J).LT.0) ISTATE(J) = 0
            IF (BL(J).NE.BU(J) .AND. ISTATE(J).EQ.3) ISTATE(J) = 0
   60    CONTINUE
      END IF
C
C     Define NFREE and KACTIV.
C     Ensure that the number of bounds and general constraints in the
C     working set does not exceed N.
C
      DO 80 J = 1, NCTOTL
         IF (NACTIV.EQ.NFREE) ISTATE(J) = 0
C
         IF (ISTATE(J).GT.0) THEN
            IF (J.LE.N) THEN
               NFREE = NFREE - 1
C
               IF (ISTATE(J).EQ.1) THEN
                  WX(J) = BL(J)
               ELSE IF (ISTATE(J).GE.2) THEN
                  WX(J) = BU(J)
               END IF
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
   80 CONTINUE
C
C     ------------------------------------------------------------------
C     If a cold start is required,  attempt to add as many
C     constraints as possible to the working set.
C     ------------------------------------------------------------------
      IF (START.EQ.'COLD') THEN
C
C        See if any bounds are violated or nearly satisfied.
C        If so,  add these bounds to the working set and set the
C        variables exactly on their bounds.
C
         J = N
C        +       while (j .ge. 1  .and.  NACTIV .lt. NFREE) do
  100    IF (J.GE.1 .AND. NACTIV.LT.NFREE) THEN
            IF (ISTATE(J).EQ.0) THEN
               B1 = BL(J)
               B2 = BU(J)
               IS = 0
               IF (B1.GT.BIGLOW) THEN
                  IF (WX(J)-B1.LE.(ONE+ABS(B1))*TOLACT) IS = 1
               END IF
               IF (B2.LT.BIGUPP) THEN
                  IF (B2-WX(J).LE.(ONE+ABS(B2))*TOLACT) IS = 2
               END IF
               IF (IS.GT.0) THEN
                  ISTATE(J) = IS
                  IF (IS.EQ.1) WX(J) = B1
                  IF (IS.EQ.2) WX(J) = B2
                  NFREE = NFREE - 1
               END IF
            END IF
            J = J - 1
            GO TO 100
C           +       end while
         END IF
C
C        ---------------------------------------------------------------
C        The following loop finds the linear constraint (if any) with
C        smallest residual less than or equal to TOLACT  and adds it
C        to the working set.  This is repeated until the working set
C        is complete or all the remaining residuals are too large.
C        ---------------------------------------------------------------
C        First, compute the residuals for all the constraints not in the
C        working set.
C
         IF (NCLIN.GT.0 .AND. NACTIV.LT.NFREE) THEN
            DO 120 I = 1, NCLIN
               IF (ISTATE(N+I).LE.0) AX(I) = DDOT(N,A(I,1),LDA,WX,1)
  120       CONTINUE
C
            IS = 1
            TOOBIG = TOLACT + TOLACT
C
C           +          while (is .gt. 0  .and.  NACTIV .lt. NFREE) do
  140       IF (IS.GT.0 .AND. NACTIV.LT.NFREE) THEN
               IS = 0
               RESMIN = TOLACT
C
               DO 160 I = 1, NCLIN
                  J = N + I
                  IF (ISTATE(J).EQ.0) THEN
                     B1 = BL(J)
                     B2 = BU(J)
                     RESL = TOOBIG
                     RESU = TOOBIG
                     IF (B1.GT.BIGLOW) RESL = ABS(AX(I)-B1)/(ONE+ABS(B1)
     *                                        )
                     IF (B2.LT.BIGUPP) RESU = ABS(AX(I)-B2)/(ONE+ABS(B2)
     *                                        )
                     RESIDL = MIN(RESL,RESU)
                     IF (RESIDL.LT.RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN = I
                        IS = 1
                        IF (RESL.GT.RESU) IS = 2
                     END IF
                  END IF
  160          CONTINUE
C
               IF (IS.GT.0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 140
C              +          end while
            END IF
         END IF
      END IF
C
      IF (VERTEX .AND. NACTIV.LT.NFREE) THEN
C        ---------------------------------------------------------------
C        Find an initial vertex by temporarily fixing some variables.
C        ---------------------------------------------------------------
C        Compute lengths of columns of selected linear constraints
C        (just the ones corresponding to variables eligible to be
C        temporarily fixed).
C
         DO 200 J = 1, N
            IF (ISTATE(J).EQ.0) THEN
               COLSIZ = ZERO
               DO 180 K = 1, NCLIN
                  IF (ISTATE(N+K).GT.0) COLSIZ = COLSIZ + ABS(A(K,J))
  180          CONTINUE
               WORK(J) = COLSIZ
            END IF
  200    CONTINUE
C
C        Find the  NARTIF  smallest such columns.
C        This is an expensive loop.  Later we can replace it by a
C        4-pass process (say), accepting the first col that is within
C        t  of  COLMIN, where  t = 0.0, 0.001, 0.01, 0.1 (say).
C        (This comment written in 1980).
C
C        +       while (NACTIV .lt. NFREE) do
  220    IF (NACTIV.LT.NFREE) THEN
            COLMIN = FLMAX
            DO 240 J = 1, N
               IF (ISTATE(J).EQ.0) THEN
                  IF (NCLIN.EQ.0) GO TO 260
                  COLSIZ = WORK(J)
                  IF (COLMIN.GT.COLSIZ) THEN
                     COLMIN = COLSIZ
                     JMIN = J
                  END IF
               END IF
  240       CONTINUE
            J = JMIN
C
C           Fix x(j) at its current value.
C
  260       ISTATE(J) = 4
            NARTIF = NARTIF + 1
            NFREE = NFREE - 1
            GO TO 220
C           +       end while
         END IF
      END IF
C
      JFREE = 1
      JFIX = NFREE + 1
      DO 280 J = 1, N
         IF (ISTATE(J).LE.0) THEN
            KX(JFREE) = J
            JFREE = JFREE + 1
         ELSE
            KX(JFIX) = J
            JFIX = JFIX + 1
         END IF
  280 CONTINUE
C
      IF (CMDBG) THEN
         IF (ICMDBG(1).GT.0) THEN
            WRITE (REC,FMT=99998) N - NFREE, NACTIV, NARTIF
            CALL X04BAY(IPRINT,4,REC)
         END IF
         IF (ICMDBG(2).GT.0) CALL E04MFN('E04MFT',
     *                            'variables after  crash',WX,N)
      END IF
      RETURN
C
C     End of  E04MFT.  (CMCRSH)
C
99999 FORMAT (/' //E04MFT// START NCLIN NCTOTL',/' //E04MFT// ',A4,I6,
     *       I7)
99998 FORMAT (/' //E04MFT// Working set selected ...  ',/' //E04MFT// ',
     *       'NFIXED NACTIV NARTIF      ',/' //E04MFT// ',I6,2I7)
      END

      SUBROUTINE F06PFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := A'*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PFF/DTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := A*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PFF (DTRMV ).
C
      END

      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END
      SUBROUTINE F06ECF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06ECF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06ECF. ( DAXPY )
C
      END



      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END

      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP, TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY, M4, N4
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JX = KX
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 60, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 50, I = 1, M
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( I, J ) )
     $                        + TEMP2*A( I, J + 1 ) )
     $                        + TEMP3*A( I, J + 2 ) )
     $                        + TEMP4*A( I, J + 3 ) )
   50             CONTINUE
               END IF
               JX = JX + 4*INCX
   60       CONTINUE
C**** Clean-up loop ****************************************************
            DO 80, J = N4 + 1, N, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 70, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 100, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 90, I = 1, M
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( I, J ) )
     $                         + TEMP2*A( I, J + 1 ) )
     $                         + TEMP3*A( I, J + 2 ) )
     $                         + TEMP4*A( I, J + 3 ) )
                     IY = IY + INCY
   90             CONTINUE
               END IF
               JX = JX + 4*INCX
  100       CONTINUE
C**** Clean-up loop ****************************************************
            DO 120, J = N4 + 1, N, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 110, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY = IY + INCY
  110             CONTINUE
               END IF
               JX = JX + INCX
  120       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 140, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 130, I = 1, N
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( J, I ) )
     $                        + TEMP2*A( J + 1, I ) )
     $                        + TEMP3*A( J + 2, I ) )
     $                        + TEMP4*A( J + 3, I ) )
  130             CONTINUE
               END IF
               JX = JX + 4*INCX
  140       CONTINUE
C**** Clean-up loop ****************************************************
            DO 160, J = M4 + 1, M, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 150, I = 1, N
                     Y( I ) = Y( I ) + TEMP*A( J, I )
  150             CONTINUE
               END IF
               JX = JX + INCX
  160       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 180, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 170, I = 1, N
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( J, I ) )
     $                         + TEMP2*A( J + 1, I ) )
     $                         + TEMP3*A( J + 2, I ) )
     $                         + TEMP4*A( J + 3, I ) )
                     IY = IY + INCY
  170             CONTINUE
               END IF
               JX = JX + 4*INCX
  180       CONTINUE
C**** Clean-up loop ****************************************************
            DO 200, J = M4 + 1, M, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 190, I = 1, N
                     Y( IY ) = Y( IY ) + TEMP*A( J, I )
                     IY = IY + INCY
  190             CONTINUE
               END IF
               JX = JX + INCX
  200       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END

      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
      SUBROUTINE F06EGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EGF. ( DSWAP )
C
      END
      DOUBLE PRECISION FUNCTION F06EAF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DDOT
      ENTRY                     DDOT  ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               I, IX, IY
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06EAF = SUM
      RETURN
C
C     End of F06EAF. ( DDOT )
C
      END

      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END



      SUBROUTINE F06EFF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DCOPY ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EFF. ( DCOPY )
C
      END


      DOUBLE PRECISION FUNCTION F06EJF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DNRM2
      ENTRY                     DNRM2 ( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
C     ..
C
C  F06EJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     F06EJF := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      NORM, SCALE, SSQ
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. External Subroutines ..
      EXTERNAL              F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
         CALL F06FJF( N, X, INCX, SCALE, SSQ )
         NORM  = F06BMF( SCALE, SSQ )
      END IF
C
      F06EJF = NORM
      RETURN
C
C     End of F06EJF. ( DNRM2 )
C
      END

      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END

      SUBROUTINE F01QFF(PIVOT,M,N,A,LDA,ZETA,PERM,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QFF  finds  a  QR factorization  of  the  real  m by n  matrix  A,
C  incorporating  column interchanges,  so that  A  is reduced to  upper
C  triangular form  by means of  orthogonal transformations  and  column
C  permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where  Q  is an  m by m  orthogonal matrix, R  is  a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used  to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
C  triangular part of  A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  The elements of  R
C  are returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',   p = min( m, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'   ( Column interchanges )
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'   ( Scaled column interchanges )
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of DIMENSION at least ( n ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then  ZETA( k ) = 0.0,
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
C           When n.gt.m the elements  ZETA( m + 1 ), ZETA( m + 2 ), ...,
C           ZETA( n )  are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least  min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if columns  k  and  j  were  interchanged  at the  kth step.
C           Note that there are  min( m, n ) permutations.
C
C  WORK   - REAL array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length of the jth column of the permuted matrix A*P'.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit, IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to   -1  indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C        B := Q*B   and   B := Q'*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary linear algebra routine  F01QDF. The operation  B := Q*B can
C  be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'No transpose', 'Separate', M, MIN( M, N ), A, LDA,
C    $             ZETA, K, B, LDB, WORK, IFAIL )
C
C  and  B := Q'*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'Transpose', 'Separate', M, MIN( M, N ), A, LDA,
C    $             ZETA, K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine F01QDF for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  F06QJF.  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL F06QJF( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL F06QJF( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be  replaced  by  N  in the above  two calls.  The operation
C  B := B*P  can be obtained by the call:
C
C     CALL F06QJF( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL F06QJF( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine F06QJF for further details.
C
C  Operations involving  the matrix  R  can readily be performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine DTRSV.
C  If  R is nearly singular then the  NAG library routine  F02WUF can be
C  used to determine the singular value decomposition of  R.  Operations
C  involving  the matrix  X  can also be  performed by the  Level 2 BLAS
C  routines.  Matrices  of  the  form  ( R  X )  can  be  factorized  as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T  is  upper triangular and  S  is orthogonal,  using the  NAG
C  Library routine  F01QGF.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-March-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  LAMDA, ONE, ZERO
      PARAMETER         (LAMDA=1.0D-2,ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER*1       PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*), ZETA(*)
      INTEGER           PERM(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, MAXNRM, NORM, TEMP, TOL
      INTEGER           IERR, J, JMAX, K, LA
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSWAP, F06FRF, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF ((PIVOT.NE.'C') .AND. (PIVOT.NE.'c') .AND. (PIVOT.NE.'S')
     *    .AND. (PIVOT.NE.'s')) CALL P01ABW(PIVOT,'PIVOT',IFAIL,IERR,
     *                               SRNAME)
      IF (M.LT.0) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Compute eps and the initial column norms.
C
      IF (MIN(M,N).EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      EPS = X02AJF()
      DO 20 J = 1, N
         WORK(J) = DNRM2(M,A(1,J),1)
         WORK(J+N) = WORK(J)
   20 CONTINUE
C
C     Perform the factorization. TOL is the tolerance for F06FRF.
C
      LA = LDA
      DO 120 K = 1, MIN(M,N)
C
C        Find the pivot column.
C
         MAXNRM = ZERO
         JMAX = K
         IF ((PIVOT.EQ.'C') .OR. (PIVOT.EQ.'c')) THEN
            DO 40 J = K, N
               IF (WORK(J+N).GT.MAXNRM) THEN
                  MAXNRM = WORK(J+N)
                  JMAX = J
               END IF
   40       CONTINUE
         ELSE
            DO 60 J = K, N
               IF (WORK(J).GT.ZERO) THEN
                  IF (K.LE.1) THEN
                     JMAX = J
                     GO TO 80
                  ELSE IF ((WORK(J+N)/WORK(J)).GT.MAXNRM) THEN
                     MAXNRM = WORK(J+N)/WORK(J)
                     JMAX = J
                  END IF
               END IF
   60       CONTINUE
   80       CONTINUE
         END IF
         PERM(K) = JMAX
         IF (JMAX.GT.K) THEN
            CALL DSWAP(M,A(1,K),1,A(1,JMAX),1)
            TEMP = WORK(K)
            WORK(K) = WORK(JMAX)
            WORK(JMAX) = TEMP
            WORK(JMAX+N) = WORK(K+N)
         END IF
         TOL = EPS*WORK(K)
         IF (K.LT.M) THEN
C
C           Use a Householder reflection to zero the kth column of A.
C           First set up the reflection.
C
            CALL F06FRF(M-K,A(K,K),A(K+1,K),1,TOL,ZETA(K))
            IF (K.LT.N) THEN
               IF (ZETA(K).GT.ZERO) THEN
                  IF ((K+1).EQ.N) LA = M - K + 1
C
C                 Temporarily store beta and put zeta( k ) in a( k, k ).
C
                  TEMP = A(K,K)
                  A(K,K) = ZETA(K)
C
C                 We now perform the operation  A := Q( k )*A.
C
C                 Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )
C                 part of A.
C
C                 First  form   work = B'*u.  ( work  is  stored  in the
C                 elements ZETA( k + 1 ), ..., ZETA( n ). )
C
                  CALL DGEMV('Transpose',M-K+1,N-K,ONE,A(K,K+1),LA,
     *                       A(K,K),1,ZERO,ZETA(K+1),1)
C
C                 Now form  B := B - u*work'.
C
                  CALL DGER(M-K+1,N-K,-ONE,A(K,K),1,ZETA(K+1),1,A(K,K+1)
     *                      ,LA)
C
C                 Restore beta.
C
                  A(K,K) = TEMP
               END IF
C
C              Update  the  unreduced  column  norms.  Use  the  Linpack
C              criterion for when to recompute the norms, except that we
C              retain  the original column lengths throughout  and use a
C              smaller lamda.
C
               DO 100 J = K + 1, N
                  IF (WORK(J+N).GT.ZERO) THEN
                     TEMP = ABS(A(K,J))/WORK(J+N)
                     TEMP = MAX((ONE+TEMP)*(ONE-TEMP),ZERO)
                     NORM = TEMP
                     TEMP = ONE + LAMDA*TEMP*(WORK(J+N)/WORK(J))**2
                     IF (TEMP.GT.ONE) THEN
                        WORK(J+N) = WORK(J+N)*SQRT(NORM)
                     ELSE
                        WORK(J+N) = DNRM2(M-K,A(K+1,J),1)
                     END IF
                  END IF
  100          CONTINUE
            END IF
         END IF
  120 CONTINUE
C
C     Set the final  ZETA  when  m.le.n.
C
      IF (M.LE.N) ZETA(M) = ZERO
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QFF. ( SGEQRP )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
      SUBROUTINE F06FDF( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06FDF performs the operation
C
C     y := alpha*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     .. External Subroutines ..
      EXTERNAL           F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( ALPHA.EQ.ZERO ).AND.( INCY.NE.0 ) )THEN
            CALL F06FBF( N, ZERO, Y, ABS( INCY ) )
         ELSE
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FDF. ( SSCMV )
C
      END

      SUBROUTINE F01QCF(M,N,A,LDA,ZETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QCF  finds  the  QR factorization  of the real  m by n,  m .ge. n,
C  matrix A,  so that  A is reduced to upper triangular form by means of
C  orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m orthogonal matrix and  R  is an  n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used  to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
C  triangular part of  A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  The elements of  R
C  are returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( n )*Q( n - 1 )*...*Q( 1 ) )'.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the upper triangular matrix R and the  M by N strictly lower
C           triangular  part   of   A   will  contain  details   of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of DIMENSION at least ( n ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then  ZETA( k ) = 0.0,
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On successful  exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  If  on  entry,  IFAIL  was  either  -1 or 0  then  further diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C        B := Q*B   and   B := Q'*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  NAG Library routine  F01QDF. The  operation  B := Q*B can be obtained
C  by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'No transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  and  B := Q'*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'Transpose', 'Separate', M, N, A, LDA, ZETA,
C    $             K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine F01QDF for further
C  details.
C
C  The first k columns of the orthogonal matrix Q can either be obtained
C  by setting  B to the first k columns of the unit matrix and using the
C  first of the above two calls,  or by calling the  NAG Library routine
C  F01QEF, which overwrites the k columns of Q on the first k columns of
C  the array A.  Q is obtained by the call:
C
C     CALL F01QEF( 'Separate', M, N, K, A, LDA, ZETA, WORK, IFAIL )
C
C  As above WORK must be a k element array.  If K is larger than N, then
C  A must have been declared to have at least K columns.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV  (see Chapter F06), but note
C  that no test for  near singularity  of  R  is incorporated in DTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02WUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-December-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           IERR, K, LA
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, F06FRF, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the factorization.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      LA = LDA
      DO 20 K = 1, MIN(M-1,N)
C
C        Use a  Householder reflection  to  zero the  kth column  of  A.
C        First set up the reflection.
C
         CALL F06FRF(M-K,A(K,K),A(K+1,K),1,ZERO,ZETA(K))
         IF ((ZETA(K).GT.ZERO) .AND. (K.LT.N)) THEN
            IF ((K+1).EQ.N) LA = M - K + 1
C
C           Temporarily  store  beta and  put  zeta( k )  in  a( k, k ).
C
            TEMP = A(K,K)
            A(K,K) = ZETA(K)
C
C           We now perform the operation  A := Q( k )*A.
C
C           Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )  part
C           of  A.
C
C           First form   work = B'*u.  ( work  is stored in the elements
C           ZETA( k + 1 ), ..., ZETA( n ). )
C
            CALL DGEMV('Transpose',M-K+1,N-K,ONE,A(K,K+1),LA,A(K,K),1,
     *                 ZERO,ZETA(K+1),1)
C
C           Now form  B := B - u*work'.
C
            CALL DGER(M-K+1,N-K,-ONE,A(K,K),1,ZETA(K+1),1,A(K,K+1),LA)
C
C           Restore beta.
C
            A(K,K) = TEMP
         END IF
   20 CONTINUE
C
C     Set the final  ZETA  when  m.eq.n.
C
      IF (M.EQ.N) ZETA(N) = ZERO
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QCF. ( SGEQR )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
      DOUBLE PRECISION FUNCTION F06RJF(NORM,UPLO,DIAG,M,N,A,LDA,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANTR
C     ENTRY                DLANTR(NORM,UPLO,DIAG,M,N,A,LDA,WORK)
C
C  Purpose
C  =======
C
C  DLANTR  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  trapezoidal or triangular matrix A.
C
C  Description
C  ===========
C
C  DLANTR returns the value
C
C     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
C              (
C              ( norm1(A),         NORM = '1', 'O' or 'o'
C              (
C              ( normI(A),         NORM = 'I' or 'i'
C              (
C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
C
C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
C  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
C
C  Arguments
C  =========
C
C  NORM    (input) CHARACTER*1
C          Specifies the value to be returned in DLANTR as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower trapezoidal.
C          = 'U':  Upper trapezoidal
C          = 'L':  Lower trapezoidal
C          Note that A is triangular instead of trapezoidal if M = N.
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A has unit diagonal.
C          = 'N':  Non-unit diagonal
C          = 'U':  Unit diagonal
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0, and if
C          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0, and if
C          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The trapezoidal matrix A (A is triangular if M = N).
C          If UPLO = 'U', the leading m by n upper trapezoidal part of
C          the array A contains the upper trapezoidal matrix, and the
C          strictly lower triangular part of A is not referenced.
C          If UPLO = 'L', the leading m by n lower trapezoidal part of
C          the array A contains the lower trapezoidal matrix, and the
C          strictly upper triangular part of A is not referenced.  Note
C          that when DIAG = 'U', the diagonal elements of A are not
C          referenced and are assumed to be one.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(M,1).
C
C  WORK    (workspace) REAL array, dimension (LWORK),
C          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
C          referenced.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, ZERO
      PARAMETER                        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER                          LDA, M, N
      CHARACTER                        DIAG, NORM, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(LDA,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J
      LOGICAL                          UDIAG
C     .. External Subroutines ..
      EXTERNAL                         F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      IF (MIN(M,N).EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M' .OR. NORM.EQ.'m')) THEN
C
C        Find max(abs(A(i,j))).
C
         IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
            VALUE = ONE
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 40 J = 1, N
                  DO 20 I = 1, MIN(M,J-1)
                     VALUE = MAX(VALUE,ABS(A(I,J)))
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = J + 1, M
                     VALUE = MAX(VALUE,ABS(A(I,J)))
   60             CONTINUE
   80          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 120 J = 1, N
                  DO 100 I = 1, MIN(M,J)
                     VALUE = MAX(VALUE,ABS(A(I,J)))
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 I = J, M
                     VALUE = MAX(VALUE,ABS(A(I,J)))
  140             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE IF (((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1')) THEN
C
C        Find norm1(A).
C
         VALUE = ZERO
         UDIAG = (DIAG.EQ.'U' .OR. DIAG.EQ.'u')
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 220 J = 1, N
               IF ((UDIAG) .AND. (J.LE.M)) THEN
                  SUM = ONE
                  DO 180 I = 1, J - 1
                     SUM = SUM + ABS(A(I,J))
  180             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 200 I = 1, MIN(M,J)
                     SUM = SUM + ABS(A(I,J))
  200             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  220       CONTINUE
         ELSE
            DO 280 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 240 I = J + 1, M
                     SUM = SUM + ABS(A(I,J))
  240             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 260 I = J, M
                     SUM = SUM + ABS(A(I,J))
  260             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  280       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
C
C        Find normI(A).
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 300 I = 1, M
                  WORK(I) = ONE
  300          CONTINUE
               DO 340 J = 1, N
                  DO 320 I = 1, MIN(M,J-1)
                     WORK(I) = WORK(I) + ABS(A(I,J))
  320             CONTINUE
  340          CONTINUE
            ELSE
               DO 360 I = 1, M
                  WORK(I) = ZERO
  360          CONTINUE
               DO 400 J = 1, N
                  DO 380 I = 1, MIN(M,J)
                     WORK(I) = WORK(I) + ABS(A(I,J))
  380             CONTINUE
  400          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 420 I = 1, N
                  WORK(I) = ONE
  420          CONTINUE
               DO 440 I = N + 1, M
                  WORK(I) = ZERO
  440          CONTINUE
               DO 480 J = 1, N
                  DO 460 I = J + 1, M
                     WORK(I) = WORK(I) + ABS(A(I,J))
  460             CONTINUE
  480          CONTINUE
            ELSE
               DO 500 I = 1, M
                  WORK(I) = ZERO
  500          CONTINUE
               DO 540 J = 1, N
                  DO 520 I = J, M
                     WORK(I) = WORK(I) + ABS(A(I,J))
  520             CONTINUE
  540          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 560 I = 1, M
            VALUE = MAX(VALUE,WORK(I))
  560    CONTINUE
      ELSE IF (((NORM.EQ.'F' .OR. NORM.EQ.'f'))
     *         .OR. ((NORM.EQ.'E' .OR. NORM.EQ.'e'))) THEN
C
C        Find normF(A).
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = MIN(M,N)
               DO 580 J = 2, N
                  CALL F06FJF(MIN(M,J-1),A(1,J),1,SCALE,SUM)
  580          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 600 J = 1, N
                  CALL F06FJF(MIN(M,J),A(1,J),1,SCALE,SUM)
  600          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = MIN(M,N)
               DO 620 J = 1, N
                  CALL F06FJF(M-J,A(MIN(M,J+1),J),1,SCALE,SUM)
  620          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 640 J = 1, N
                  CALL F06FJF(M-J+1,A(J,J),1,SCALE,SUM)
  640          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RJF = VALUE
      RETURN
C
C     End of F06RJF (DLANTR)
C
      END
      SUBROUTINE F06BCF( T, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-602 (MAR 1988).
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S, T
C     ..
C
C  F06BCF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  c is always non-negative and s has the same sign as t, so that
C
C     c = 1.0/sqrt( 1.0 + t**2 ),   s = t/sqrt( 1.0 + t**2 ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 28-February-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified 19-August-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C        No longer sets s to zero when t is less than eps.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST, EPS, REPS, RRTEPS, RTEPS
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF
      EXTERNAL           X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE               FIRST, EPS, REPS, RTEPS, RRTEPS
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         EPS    =  X02AJF( )
         REPS   =  1/EPS
         RTEPS  =  SQRT( EPS )
         RRTEPS =  1/RTEPS
      END IF
C
      ABST = ABS( T )
      IF( ABST.LT.RTEPS )THEN
         C = ONE
         S = T
      ELSE IF( ABST.GT.RRTEPS )THEN
         C = 1/ABST
         S = SIGN( ONE, T )
      ELSE
         C = 1/SQRT( 1 + ABST**2 )
         S = C*T
      END IF
C
      RETURN
C
C     End of F06BCF. ( SCSG )
C
      END
      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FJF returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
C
C     End of F06FJF. ( SSSQ )
C
      END
      SUBROUTINE F06FLF( N, X, INCX, XMAX, XMIN )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   XMAX, XMIN
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FLF returns the values xmax and xmin given by
C
C     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
C             i                              i
C
C  If n is less than unity then xmax and xmin are returned as zero.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         XMAX = ZERO
         XMIN = ZERO
      ELSE
         XMAX = ABS( X( 1 ) )
         XMIN = XMAX
         DO 10 IX = 1 + INCX, 1 + ( N - 1 )*INCX, INCX
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            XMIN = MIN( XMIN, ABS( X( IX ) ) )
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06FLF. ( SCOND )
C
      END
      SUBROUTINE F06DBF( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            CONST, INCX, N
C     .. Array Arguments ..
      INTEGER            X( * )
C     ..
C
C  F06DBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 18-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.0 )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = 0
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06DBF. ( ILOAD )
C
      END

      INTEGER FUNCTION F06KLF( N, X, INCX, TOL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION         TOL
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06KLF finds the first element of the n element vector x for which
C
C     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
C
C  and returns the value ( k - 1 ) in the function name F06KLF. If no
C  such k exists then F06KLF is returned as n.
C
C  If tol is supplied as less than zero then the value epsmch, where
C  epsmch is the relative machine precision, is used in place of tol.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION         ZERO
      PARAMETER              ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION         TL, XMAX
      INTEGER                  IX, K
C     .. External Functions ..
      DOUBLE PRECISION         X02AJF
      EXTERNAL                 X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, MAX
C     ..
C     .. Executable Statements ..
      K = 0
      IF( N.GE.1 )THEN
         IX = 1
         IF( TOL.LT.ZERO )THEN
            TL = X02AJF( )
         ELSE
            TL = TOL
         END IF
         XMAX = ABS( X( IX ) )
C
C+       WHILE( K.LT.N )LOOP
   10    IF   ( K.LT.N )THEN
            IF( ABS( X( IX ) ).LE.TL*XMAX )
     $         GO TO 20
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            K    = K  + 1
            IX   = IX + INCX
            GO TO 10
         END IF
C+       END WHILE
C
      END IF
C
   20 F06KLF = K
      RETURN
C
C     End of F06KLF. ( ISRANK )
C
      END
      SUBROUTINE F06QTF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QTF performs the transformation
C
C     R := P*U*Q'  when  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     R := Q*U*P'  when  SIDE = 'R' or 'r'  ( Right-hand side ),
C
C  where  U and R  are  n by n  upper  triangular  matrices,   P  is  an
C  orthogonal matrix,  consisting of a given sequence of plane rotations
C  to be  applied  in  planes  k1 to k2,  and  Q  is  a  unitary  matrix
C  consisting of a sequence of plane rotations, applied in planes  k1 to
C  k2,  chosen to make  R  upper triangular.
C
C  When  SIDE = 'L' or 'l'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k2 - 1 )*...*Q( k1 + 1 )*Q( k1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  When  SIDE = 'R' or 'r'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k1 )*Q( k1 + 1 )*...*Q( k2 - 1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  The  upper  triangular  matrix  U  must  be  supplied  in the  n by n
C  leading upper triangular part of  A,  and this  is overwritten by the
C  upper triangular matrix  R.  The cosine  and  sine  that  define  the
C  plane rotation matrix  P( k )  must be supplied in  c( k ) and s( k )
C  respectively,  and  the two by two rotation part of  P( k ),  T( k ),
C  is assumed to be of the form
C
C     T( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The cosine  and  sine that define  Q( k )  are overwritten on  c( k )
C  and  s( k )  respectively and the two by two rotation part of  Q( k )
C  will have the form of  T( k )  above.
C
C  If  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
C  k2  is greater than  n  then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 26-November-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, FILL, STEMP, TEMP
      INTEGER            I, I1, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the left-hand transformations,  column by column,  to the
C        triangular part of  U,  but not to  anywhere  that would  cause
C        fill.
C
         DO 20 J = K1 + 1, N
C
C           Apply  P( k1 ) ... P( j - 1 )  to column j.
C
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J - 1, K2 - 1 )
               A( I, J ) = S( I )*A( I + 1, J ) + C( I )*AIJ
               AIJ = C( I )*A( I + 1, J ) - S( I )*AIJ
   10       CONTINUE
            A( I, J ) = AIJ
   20    CONTINUE
C
C           Now apply each  left-hand tranformation  to form the fill-in
C           elements and apply a  right-hand transformation to eliminate
C           the fill-in element.
C
         DO 40 J = K1, K2 - 1
C
C           Apply  P( j )  to the jth diagonal element  and the  fill-in
C           position.
C
            FILL = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
C
C           Now  set up  the rotation  Q( j )  to eliminate the  fill-in
C           element,  and  apply  Q( j )  to  the  jth  and  ( j + 1 )th
C           columns.
C
            CALL F06BAF( A( J + 1, J + 1 ), FILL, CTEMP, STEMP )
            C( J ) = CTEMP
            S( J ) = -STEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMP = -STEMP
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        We intermingle the  left and right hand transformations so that
C        at the kth step we form
C
C           A := Q( k )*A*P( k )'.
C
C        First  apply  the  transformations  in  columns  k2 back to k1.
C
         DO 60 J = K2 - 1, K1, -1
C
C           First apply  P( j ).
C
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               CTEMP = C( J )
               STEMP = S( J )
               DO 50 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   50          CONTINUE
C
C              Next form the fill-in element  a( j + 1, j )  by applying
C              P( j ).
C
               FILL = S( J )*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = C( J )*A( J + 1, J + 1 )
C
C              Now set up the rotation  Q( j )  to eliminate the fill-in
C              element.
C
               CALL F06BAF( A( J, J ), FILL, C( J ), S( J ) )
            END IF
   60    CONTINUE
C
C        Finally  apply  Q( k2 - 1 ) ... Q( k1 )  to columns  n  back to
C        ( k1 + 1 ).
C
         DO 80 J = N, K1 + 1, -1
            I1 = MIN( K2, J )
            AIJ = A( I1, J )
            DO 70 I = I1 - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   70       CONTINUE
            A( K1, J ) = AIJ
   80    CONTINUE
      END IF
      RETURN
C
C     End of F06QTF. ( SUTSQR )
C
      END

      SUBROUTINE F06FCF( N, D, INCD, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCD, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), X( * )
C     ..
C
C  F06FCF performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, ID, IX
C     .. External Subroutines ..
      EXTERNAL           DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCD.EQ.0 ).AND.( INCX.NE.0 ) )THEN
            CALL DSCAL( N, D( 1 ), X, ABS( INCX ) )
         ELSE IF( ( INCD.EQ.INCX ).AND.( INCD.GT.0 ) )THEN
            DO 10, ID = 1, 1 + ( N - 1 )*INCD, INCD
               X( ID ) = D( ID )*X( ID )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCD.GT.0 )THEN
               DO 20, ID = 1, 1 + ( N - 1 )*INCD, INCD
                  X( IX ) = D( ID )*X( IX )
                  IX      = IX              + INCX
   20          CONTINUE
            ELSE
               ID = 1 - ( N - 1 )*INCD
               DO 30, I = 1, N
                  X( IX ) = D( ID )*X( IX )
                  ID      = ID              + INCD
                  IX      = IX              + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FCF. ( SDSCL )
C
      END
      SUBROUTINE F06FQF( PIVOT, DIRECT, N, ALPHA, X, INCX, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
      CHARACTER*1        DIRECT, PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * ), X( * )
C     ..
C
C  F06FQF generates the parameters of an orthogonal matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P*( alpha ) = ( beta ),
C          (   x   )   (   0  )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P*(   x   ) = (   0  ),
C          ( alpha ) = ( beta )
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c( k ), and sine, s( k ) that define
C  the matrix P( k ), such that the two by two rotation part of P( k ),
C  R( k ), has the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  On entry, ALPHA must contain  the scalar alpha and on exit, ALPHA is
C  overwritten by beta. The cosines and sines are returned in the arrays
C  C and S and the vector x is overwritten by the tangents of the plane
C  rotations ( t( k ) = s( k )/c( k ) ).
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 19-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
            IX = 1 + ( N - 1 )*INCX
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
               DO 10, I = N, 2, -1
                  CALL F06BAF( X( IX - INCX ), X( IX ), C( I ), S( I ) )
                  IX = IX - INCX
   10          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( 1 ), S( 1 ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( alpha ) := (  c  s )*( alpha  )
C                 (   0   )    ( -s  c ) ( x( i ) )
C
C              which is equivalent to
C
C                 (   0   ) := ( c  -s )*( x( i ) )
C                 ( alpha )    ( s   c ) ( alpha  )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 20, I = N, 1, -1
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      - INCX
   20          CONTINUE
            END IF
         ELSE IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
            IX = 1
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( x( i + 1 ) ) := (  c  s )*( x( i + 1 ) )
C                 (    0       )    ( -s  c ) ( x( i )     )
C
C              which is equivalent to
C
C                 (    0       ) := ( c  -s )*( x( i )     )
C                 ( x( i + 1 ) )    ( s   c ) ( x( i + 1 ) )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 30, I = 1, N - 1
                  CALL F06BAF( X( IX + INCX ), X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      + INCX
   30          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( N ), S( N ) )
               S( N )  = -S( N )
               X( IX ) = -X( IX )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
               DO 40, I = 1, N
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  IX = IX + INCX
   40          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FQF. ( SSROTG )
C
      END
      DOUBLE PRECISION FUNCTION F06BLF( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
      LOGICAL                           FAIL
C     ..
C
C  F06BLF returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
C
C  where  flmax  is a large value, via the function name. In addition if
C  a/b would overflow then  fail is returned as true, otherwise  fail is
C  returned as false.
C
C  Note that when  a and b  are both zero, fail is returned as true, but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that  abs( div ) = flmax.
C
C  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSB, DIV, FLMAX, FLMIN
      LOGICAL               FIRST
C     .. External Functions ..
      DOUBLE PRECISION      X02AMF
      EXTERNAL              X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SIGN
C     .. Save statement ..
      SAVE                  FIRST, FLMIN, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         DIV = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            FLMIN =  X02AMF( )
            FLMAX =  1/FLMIN
         END IF
C
         IF( B.EQ.ZERO )THEN
            DIV  =  SIGN( FLMAX, A )
            FAIL = .TRUE.
         ELSE
            ABSB = ABS( B )
            IF( ABSB.GE.ONE )THEN
               FAIL = .FALSE.
               IF( ABS( A ).GE.ABSB*FLMIN )THEN
                  DIV = A/B
               ELSE
                  DIV = ZERO
               END IF
            ELSE
               IF( ABS( A ).LE.ABSB*FLMAX )THEN
                  FAIL = .FALSE.
                  DIV  =  A/B
               ELSE
                  FAIL = .TRUE.
                  DIV  = FLMAX
                  IF( ( ( A.LT.ZERO ).AND.( B.GT.ZERO ) ).OR.
     $                ( ( A.GT.ZERO ).AND.( B.LT.ZERO ) )     )
     $               DIV = -DIV
               END IF
            END IF
         END IF
      END IF
C
      F06BLF = DIV
      RETURN
C
C     End of F06BLF. ( SDIV )
C
      END
      SUBROUTINE F06QSF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QSF restores an upper spiked matrix  H to upper triangular form by
C  applying a sequence of plane rotations, in planes  k1 up to k2,  from
C  either the left, or the right.
C
C  The matrix  H is assumed to have non-zero elements only in the spiked
C  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
C  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
C  elements s( k ).
C
C  When  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     H  is  assumed  to have a  row spike  and is restored to the upper
C     triangular matrix  R as
C
C        R = P*H,
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C     P( k )  being a  plane rotation  matrix for the  ( k, k2 )  plane.
C
C  When  SIDE = 'R' or 'r'  ( Right-hand side )
C
C     H  is assumed to have a  column spike and is restored to the upper
C     triangular matrix R as
C
C        R = H*P',
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C     P( k ) being a plane rotation matrix for the  ( k1, k + 1 ) plane.
C
C  The  two by two  rotation  part of  P( k ),  Q( k ),  is of  the form
C
C     Q( k ) = (  c( k )  s( k ) )
C              ( -s( k )  c( k ) )
C
C  and  c( k ) and s( k ) are returned in the kth elements of the arrays
C  C and S respectively.
C
C  The upper triangular part of the matrix  H must be supplied in the  n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore H to upper triangular form by annihilating the elements
C        in  the  spike  of  H.  The  jth rotation  is  chosen  so  that
C
C        ( h( j, j ) ) := (  c  s )*( h( j , j ) ).
C        (     0     )    ( -s  c ) ( h( k2, j ) )
C
C        Apply the rotations in columns k1 up to ( k2 - 1 ).
C
         DO 20 J = K1, K2 - 1
            SPIKE = S( J )
            DO 10 I = K1, J - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   10       CONTINUE
C
C           Set up the rotation.
C
            CALL F06BAF( A( J, J ), SPIKE, C( J ), S( J ) )
   20    CONTINUE
C
C        Apply the rotations to columns k2 up to n.
C
         DO 40 J = K2, N
            TEMP = A( K2, J )
            DO 30 I = K1, K2 - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   30       CONTINUE
            A( K2, J ) = TEMP
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore H to upper triangular form by annihilating the spike of
C        H. The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )  ),
C           (     0     )    ( -s  c ) ( h( j, k1 ) )
C
C        which can be expressed as
C
C           ( 0  h( j, j ) ) := ( h( j, k1 )  h( j, j ) )*(  c  s ).
C                                                         ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 70 J = K2, K1 + 1, -1
            CALL F06BAF( A( J, J ), S( J - 1 ), CTEMP, STEMP )
            STEMP = -STEMP
            S( J - 1 ) = STEMP
            C( J - 1 ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = J - 1, K1 + 1, -1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   50          CONTINUE
               DO 60 I = K1, 1, -1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   60          CONTINUE
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QSF. ( SUSQR )
C
      END
      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END
      SUBROUTINE F06DFF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      INTEGER            X( * ), Y( * )
C     ..
C
C  F06DFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 10-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06DFF. ( ICOPY )
C
      END
      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION      FLMAX, FLMIN, NORM, SQT
      LOGICAL               FIRST
C     .. External Functions ..
      DOUBLE PRECISION      X02AMF
      EXTERNAL              X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      SAVE                  FIRST, FLMAX
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         FLMIN =  X02AMF( )
         FLMAX =  1/FLMIN
      END IF
C
      SQT = SQRT( SSQ )
      IF( SCALE.LT.FLMAX/SQT )THEN
         NORM = SCALE*SQT
      ELSE
         NORM = FLMAX
      END IF
C
      F06BMF = NORM
      RETURN
C
C     End of F06BMF. ( SNORM )
C
      END
      SUBROUTINE F06QVF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QVF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The rotations are applied
C  in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = C( J )*A( J, J )
               S( J ) = -S( J )*A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the   additional  sub-diagonal  elements,   storing
C        h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP = S( J )
               CTEMP = C( J )
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
               S( J ) = STEMP*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = CTEMP*A( J + 1, J + 1 )
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QVF. ( SUTSRH )
C
      END
      SUBROUTINE F06QKF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   PERM( * ), B( LDB, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QKF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and  index( i )  interchanged. Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This  routine is intended  for use in conjunction with  Nag auxiliary
C  routines  that  perform  interchange  operations,  such  as  sorting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the value of n.  N must be at least
C           zero.  When  N = 0  then an  immediate  return  is effected.
C
C           Unchanged on exit.
C
C  PERM   - REAL             array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain  the  n indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine. It is assumed that the statement
C           INDEX = PERM( I )  returns the correct integer in  INDEX, so
C           that,  if necessary,  PERM( I )  should contain a real value
C           slightly larger than  INDEX.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the number of rows of  B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array  of  DIMENSION ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and before  entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On exit,   B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 11-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
      INTEGER            I, J, L
      DOUBLE PRECISION   TEMP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QKF. ( SGEAPR )
C
      END
      SUBROUTINE F06QNZ(SIDE,N,K1,K2,S,A,LDA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  F06QNZ applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K1, K2, LDA, N
      CHARACTER*1       SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AIJ, TEMP
      INTEGER           I, J
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF ((MIN(N,K1).LT.1) .OR. (K2.LE.K1) .OR. (K2.GT.N)) RETURN
      IF ((SIDE.EQ.'L') .OR. (SIDE.EQ.'l')) THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 40 J = N, K1, -1
            IF (J.GE.K2) THEN
               AIJ = A(K2,J)
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = ZERO
               S(J) = A(J,J)
            END IF
            DO 20 I = MIN(K2,J) - 1, K1, -1
               TEMP = A(I,J)
               A(I+1,J) = TEMP
               AIJ = AIJ
   20       CONTINUE
            A(K1,J) = AIJ
   40    CONTINUE
      ELSE IF ((SIDE.EQ.'R') .OR. (SIDE.EQ.'r')) THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 80 J = K1, K2 - 1
            DO 60 I = 1, J
               TEMP = A(I,J+1)
               A(I,J+1) = A(I,J)
               A(I,J) = TEMP
   60       CONTINUE
            S(J) = A(J+1,J+1)
            A(J+1,J+1) = ZERO
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QNZ. ( SUTSRH )
C
      END
      SUBROUTINE F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, ZETA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FRF generates details of a generalized Householder reflection such
C  that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA, EPS, SCALE, SSQ
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF
      EXTERNAL           X02AJF
C     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
C     .. Save statement ..
      SAVE               EPS, FIRST
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         ZETA = ZERO
      ELSE IF( ( N.EQ.1 ).AND.( X( 1 ).EQ.ZERO ) )THEN
         ZETA = ZERO
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            EPS   =  X02AJF( )
         END IF
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         IF( N.EQ.1 )THEN
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            IF( ALPHA.EQ.ZERO )THEN
               ZETA   =  ONE
               ALPHA  =  ABS ( X( 1 ) )
               X( 1 ) = -SIGN( ONE, X( 1 ) )
            ELSE IF( ABS( X( 1 ) ).LE.MAX( EPS*ABS( ALPHA ), TOL ) )THEN
               ZETA   =  ZERO
            ELSE
               IF( ABS( ALPHA ).GE.ABS( X( 1 ) ) )THEN
                  BETA = ABS( ALPHA ) *SQRT( 1 + ( X( 1 )/ALPHA )**2 )
               ELSE
                  BETA = ABS( X( 1 ) )*SQRT( 1 + ( ALPHA/X( 1 ) )**2 )
               END IF
               ZETA = SQRT( ( ABS( ALPHA ) + BETA )/BETA )
               IF( ALPHA.GE.ZERO )
     $            BETA = -BETA
               X( 1 ) = -X( 1 )/( ZETA*BETA )
               ALPHA  = BETA
            END IF
         ELSE
C
C           Now P is larger than 2 by 2.
C
            SSQ   = ONE
            SCALE = ZERO
            CALL F06FJF( N, X, INCX, SCALE, SSQ )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            IF( ( SCALE.EQ.ZERO ).OR.
     $          ( SCALE.LE.MAX( EPS*ABS( ALPHA ), TOL ) ) )THEN
               ZETA  = ZERO
            ELSE IF( ALPHA.EQ.ZERO )THEN
               ZETA  = ONE
               ALPHA = SCALE*SQRT( SSQ )
               CALL DSCAL( N, -1/ALPHA, X, INCX )
            ELSE
               IF( SCALE.LT.ABS( ALPHA ) )THEN
                  BETA = ABS( ALPHA )*SQRT( 1 + SSQ*( SCALE/ALPHA )**2 )
               ELSE
                  BETA = SCALE       *SQRT( SSQ +   ( ALPHA/SCALE )**2 )
               END IF
               ZETA = SQRT( ( BETA + ABS( ALPHA ) )/BETA )
               IF( ALPHA.GT.ZERO )
     $            BETA = -BETA
               CALL DSCAL( N, -1/( ZETA*BETA ), X, INCX )
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FRF. ( SGRFG )
C
      END
      SUBROUTINE F06QZZ(HESS,N,K1,K2,C,S,A,LDA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  F06QZZ  either applies a  given sequence  of  plane rotations  to the
C  right of the n by n reverse lower triangular matrix T, to transform T
C  to a  reverse lower Hessenberg matrix  H, or restores a reverse lower
C  Hessenberg matrix H to reverse lower triangular form T, by applying a
C  sequence of plane rotations from the right.
C
C  The rotations are applied  in planes k1 up to k2.
C
C  When   HESS = 'C' or 'c',   ( Create ),  then   the   reverse   lower
C  Hessenberg matrix, H, is formed as
C
C     H = T*P',
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  T must be supplied in the n by n reverse lower triangular
C  part  of the array  A,  and this is overwritten by the  reverse lower
C  triangular part of  H.
C
C  The super-diagonal elements of  H, h( n - k, k ), are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C  When   HESS = 'R' or 'r',   ( Remove ),  then   the   reverse   lower
C  Hessenberg matrix  H  is  assumed  to  have  non-zero  super-diagonal
C  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only and  h( n - k, k ) must be supplied in  s( k ). H is restored to
C  the reverse lower triangular matrix T as
C
C     T = H*P',
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation for the  ( k, k + 1 ) plane. The cosine
C  and  sine  that  define  P( k )  are  returned  in  c( k ) and s( k )
C  respectively.  The  two by two  rotation part of  P( k ),  R( k ), is
C  of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The reverse lower triangular part of the matrix H must be supplied in
C  the  n by n  reverse  lower  triangular  part  of  A,   and  this  is
C  overwritten by the reverse triangular matrix T.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C  When   n = 7, k1 = 2 and k2 = 5   then  T  and  H  are  of  the  form
C
C     T = ( 0  0  0  0  0  0  X ),   H = ( 0  0  0  0  0  0  X ).
C         ( 0  0  0  0  0  X  X )        ( 0  0  0  0  X  X  X )
C         ( 0  0  0  0  X  X  X )        ( 0  0  0  X  X  X  X )
C         ( 0  0  0  X  X  X  X )        ( 0  0  X  X  X  X  X )
C         ( 0  0  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
C         ( 0  X  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
C         ( X  X  X  X  X  X  X )        ( X  X  X  X  X  X  X )
C
C
C  This routine  is  principally intended  for use  with the  non-linear
C  optimization routines such as E04UCF, in order to help vectorization.
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 10-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K1, K2, LDA, N
      CHARACTER*1       HESS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CTEMP, STEMP, SUPH, TEMP
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF ((MIN(N,K1).LT.1) .OR. (K2.LE.K1) .OR. (K2.GT.N)) RETURN
      IF ((HESS.EQ.'C') .OR. (HESS.EQ.'c')) THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the  additional  super-diagonal  elements,  storing
C        h( n - j, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF ((C(J).NE.ONE) .OR. (S(J).NE.ZERO)) THEN
               STEMP = S(J)
               CTEMP = C(J)
               S(J) = STEMP*A(N-J,J+1)
               A(N-J,J+1) = CTEMP*A(N-J,J+1)
               DO 20 I = N - J + 1, N
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J) = STEMP*TEMP + CTEMP*A(I,J)
   20          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF ((HESS.EQ.'R') .OR. (HESS.EQ.'r')) THEN
C
C        Restore  H to reverse lower triangular form by annihilating the
C        super-diagonal elements of  H.  The  jth rotation  is chosen so
C        that
C
C          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
C          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )
C
C        which can be expressed as
C
C           ( 0  h( n - j, n - j ) ) :=
C
C               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
C                                                            ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           R( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 80 J = K2 - 1, K1, -1
            SUPH = S(J)
            CALL F06BAF(A(N-J,J+1),SUPH,CTEMP,STEMP)
            STEMP = -STEMP
            S(J) = STEMP
            C(J) = CTEMP
            IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
               DO 60 I = N - J + 1, N
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J) = STEMP*TEMP + CTEMP*A(I,J)
   60          CONTINUE
            END IF
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QZZ.
C
      END
      SUBROUTINE F06QWF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QWF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular  matrix  U  to
C  transform  U  to an upper spiked matrix. The rotations are applied in
C  planes k1 up to k2.
C
C  The upper spiked matrix, H, is formed as
C
C     H = P*U,   when   SIDE = 'L' or 'l',  ( Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
C  formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.
C
C  The cosine and sine that define  P( k ), k = k1, k1 + 1, ..., k2 - 1,
C  must be  supplied  in  c( k ) and s( k ) respectively. The two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of H.
C
C  When  SIDE = 'L' or 'l'  then a  row spike  is  generated  in  H  and
C  when  SIDE = 'R' or 'r'  then a  column spike is generated. For a row
C  spike the elements  h( k2, k )  and for a  column spike  the elements
C  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k2.
C
         DO 20 J = N, K2, -1
            TEMP = A( K2, J )
            DO 10 I = K2 - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            A( K2, J ) = TEMP
   20    CONTINUE
C
C        Form  the spike  and apply the rotations in columns  ( k2 - 1 )
C        back to k1.
C
         DO 40 J = K2 - 1, K1, -1
            SPIKE = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
            DO 30 I = J - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   30       CONTINUE
            S( J ) = SPIKE
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply the  plane rotations to columns  ( k1 + 1 ) up to k2  and
C        form the spike.
C
         DO 70 J = K1 + 1, K2
            CTEMP = C( J - 1 )
            STEMP = S( J - 1 )
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = 1, K1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   50          CONTINUE
               DO 60 I = K1 + 1, J - 1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   60          CONTINUE
               S( J - 1 ) = STEMP*A( J, J )
               A( J, J ) = CTEMP*A( J, J )
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QWF. ( SUTSRS )
C
      END
      SUBROUTINE F06QHF( MATRIX, M, N, CONST, DIAG, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      DOUBLE PRECISION   CONST, DIAG
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  F06QHF forms the m by n matrix A given by
C
C     a( i, j ) = (  diag  i.eq.j,
C                 (
C                 ( const  i.ne.j.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  if   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C                             and only  elements  for which  i.le.j  are
C                             referenced,
C  if   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C                             and only  elements  for which  i.ge.j  are
C                             referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               A( I, J ) = CONST
   10       CONTINUE
   20    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 30 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   30       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 50 J = 1, N
            DO 40 I = 1, MIN( M, J )
               A( I, J ) = CONST
   40       CONTINUE
   50    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 60 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   60       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 80 J = 1, MIN( M, N )
            DO 70 I = J, M
               A( I, J ) = CONST
   70       CONTINUE
   80    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 90 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   90       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QHF. ( SMLOAD )
C
      END
      SUBROUTINE F06QRF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QRF restores an upper Hessenberg matrix H to upper triangular form
C  by  applying a sequence of  plane rotations  from either the left, or
C  the right.  The matrix  H  is assumed to have  non-zero  sub-diagonal
C  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).
C
C  H is restored to the upper triangular matrix R either as
C
C     R = P*H,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  or as
C
C     R = H*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
C  plane.  The cosine and sine that define P( k ) are returned in c( k )
C  and  s( k )  respectively.  The two by two  rotation part of  P( k ),
C  Q( k ), is of the form
C
C     Q( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The upper triangular part of the matrix  H  must be supplied in the n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, SUBH, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
C           (     0     )    ( -s  c ) ( h( j + 1, j ) )
C
C        Apply the rotations in columns k1 up to n.
C
         DO 20 J = K1, N
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J, K2 ) - 1
               TEMP = A( I + 1, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               AIJ = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            IF( J.LT.K2 )THEN
C
C              Set up the rotation.
C
               SUBH = S( J )
               CALL F06BAF( AIJ, SUBH, C( J ), S( J ) )
               A( J, J ) = AIJ
            ELSE
               A( K2, J ) = AIJ
            END IF
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
C           (         0         )    ( -s  c ) ( h( j + 1, j )     )
C
C        which can be expressed as
C
C           ( 0  h( j + 1, j + 1 ) ) :=
C
C               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
C                                                    ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 40 J = K2 - 1, K1, -1
            SUBH = S( J )
            CALL F06BAF( A( J + 1, J + 1 ), SUBH, CTEMP, STEMP )
            STEMP = -STEMP
            S( J ) = STEMP
            C( J ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30 I = J, 1, -1
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QRF. ( SUHQR )
C
      END
      SUBROUTINE F06QFF( MATRIX, M, N, A, LDA, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            M, N, LDA, LDB
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  F06QFF  copies  the  m by n  matrix  A  into  the  m by n  matrix  B.
C
C  If   MATRIX = 'G' or 'g'   then  A  and  B  are  regarded as  general
C                             matrices,
C  if   MATRIX = 'U' or 'u'   then  A  and  B  are  regarded  as   upper
C                             triangular,  and only  elements  for which
C                             i.le.j  are referenced,
C  if   MATRIX = 'L' or 'l'   then  A  and  B  are  regarded  as   lower
C                             triangular,  and only  elements  for which
C                             i.ge.j  are referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 40 J = 1, N
            DO 30 I = 1, MIN( M, J )
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 60 J = 1, MIN( M, N )
            DO 50 I = J, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QFF. ( SMCOPY )
C
      END
      SUBROUTINE F06BAF( A, B, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, S
C     ..
C
C  Nag  Fortran 77  version of the  SROTG BLAS,  except that c is always
C  returned as non-negative and  b  is overwritten by the tangent of the
C  angle that defines the plane rotation.
C
C  c and s are given as
C
C     c = 1.0/sqrt( 1.0 + t**2 ),   s = c*t   where   t = b/a.
C
C  When  abs( b ) .le. eps*abs( a ),  where  eps is the relative machine
C  precision as  returned by routine  X02AJF,  then  c and s  are always
C  returned as
C
C     c = 1.0  and  s = 0.0
C
C  and when  abs( a ) .le. eps*abs( b ) then c and s are always returned
C  as
C
C     c = 0.0  and  s = sign( t ).
C
C  Note that t is always returned as  b/a, unless this would overflow in
C  which  case the value  sign( t )*flmax  is returned,  where  flmax is
C  the value given by  1/X02AMF( ).
C
C  c and s  can be reconstructed from the tangent,  t,  by a call to the
C  Nag basic linear algebra routine F06BCF.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 3-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   T
      LOGICAL            FAIL
C     .. External Functions ..
      DOUBLE PRECISION   F06BLF
      EXTERNAL           F06BLF
C     .. External Subroutines ..
      EXTERNAL           F06BCF
C     ..
C     .. Executable Statements ..
      IF( B.EQ.ZERO )THEN
         C  = ONE
         S  = ZERO
      ELSE
         T  = F06BLF( B, A, FAIL )
         CALL F06BCF( T, C, S )
         A  = C*A + S*B
         B  = T
      END IF
C
      RETURN
C
C     End of F06BAF. ( SROTGC )
C
      END
      SUBROUTINE F06FBF( N, CONST, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FBF performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06FBF. ( SLOAD )
C
      END

      SUBROUTINE P01ABY(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABY increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.LT.1 )CALL P01ABY( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 23-February-1984.  Sven.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM, N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(2)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
C
C     End of P01ABY.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is ',I6)
      END


      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
c     STOP
      END

      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      INTEGER FUNCTION P01ACF(IFAIL,IERROR,SRNAME,VARBNM,NREC,REC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     P01ACF is the error-handling routine for the F06 AND F07
C     Chapters of the NAG Fortran Library. It is a slightly modified
C     version of P01ABF.
C
C     P01ACF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ACF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ACF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME, VARBNM
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR, VARLEN
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, LEN, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
         VARLEN = 0
         DO 20 I = LEN(VARBNM), 1, -1
            IF (VARBNM(I:I).NE.' ') THEN
               VARLEN = I
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 60 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   60       CONTINUE
            IF (IFAIL.NE.-13) THEN
               IF (VARLEN.NE.0) THEN
                  WRITE (MESS,FMT=99999) SRNAME, VARBNM(1:VARLEN),
     *              IERROR
               ELSE
                  WRITE (MESS,FMT=99998) SRNAME
               END IF
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ACF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': ',A,
     *       ' =',I6)
99998 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A)
      END
      SUBROUTINE P01ABW(N,NAME,INFORM,IERR,SRNAME)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     P01ABW increases the value of IERR by 1 and, if
C
C        ( mod( INFORM, 10 ).ne.1 ).or.( mod( INFORM/10, 10 ).ne.0 )
C
C     writes a message on the current error message channel giving the
C     value of N, a message to say that N is invalid and the strings
C     NAME and SRNAME.
C
C     NAME must be the name of the actual argument for N and SRNAME must
C     be the name of the calling routine.
C
C     This routine is intended for use when N is an invalid input
C     parameter to routine SRNAME. For example
C
C        IERR = 0
C        IF( N.NE.'Valid value' )
C     $     CALL P01ABW( N, 'N', IDIAG, IERR, SRNAME )
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
C     .. Scalar Arguments ..
      INTEGER           IERR, INFORM
      CHARACTER*(*)     N
      CHARACTER*(*)     NAME, SRNAME
C     .. Local Scalars ..
      INTEGER           NERR
C     .. Local Arrays ..
      CHARACTER*65      REC(3)
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IERR = IERR + 1
      IF ((MOD(INFORM,10).NE.1) .OR. (MOD(INFORM/10,10).NE.0)) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) NAME, SRNAME, N
         CALL X04BAF(NERR,' ')
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
         CALL X04BAF(NERR,REC(3))
      END IF
      RETURN
C
C
C     End of P01ABW.
C
99999 FORMAT (' *****  Parameter  ',A,'  is invalid in routine  ',A,
     *  '  ***** ',/8X,'Value supplied is',/8X,A)
      END


      INTEGER FUNCTION X02BHF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     2
      RETURN
      END
      INTEGER FUNCTION X02BJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, p.
C
C     .. Executable Statements ..
      X02BJF =    53
      RETURN
      END
      LOGICAL FUNCTION X02DAF(X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS .FALSE. IF THE SYSTEM SETS UNDERFLOWING QUANTITIES
C     TO ZERO, WITHOUT ANY ERROR INDICATION OR UNDESIRABLE WARNING
C     OR SYSTEM OVERHEAD.
C     RETURNS .TRUE. OTHERWISE, IN WHICH CASE CERTAIN LIBRARY
C     ROUTINES WILL TAKE SPECIAL PRECAUTIONS TO AVOID UNDERFLOW
C     (USUALLY AT SOME COST IN EFFICIENCY).
C
C     X IS A DUMMY ARGUMENT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        X
C     .. Executable Statements ..
      X02DAF = .FALSE.
      RETURN
      END
      SUBROUTINE X04ABF(I,NADV)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-830 (DEC 1989).
C      IF I = 0, SETS NADV TO CURRENT ADVISORY MESSAGE UNIT NUMBER
C     (STORED IN NADV1).
C     IF I = 1, CHANGES CURRENT ADVISORY MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NADV.
C
C     .. Scalar Arguments ..
      INTEGER           I, NADV
C     .. Local Scalars ..
      INTEGER           NADV1
C     .. Save statement ..
      SAVE              NADV1
C     .. Data statements ..
      DATA              NADV1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NADV = NADV1
      IF (I.EQ.1) NADV1 = NADV
      RETURN
      END
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END

      SUBROUTINE X04BAY(NOUT,NREC,REC)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAY outputs NREC records on device NOUT, by calling X04BAF.
C     If NREC is 0 then no records are output.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT, NREC
C     .. Array Arguments ..
      CHARACTER*(*)     REC(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Executable Statements ..
      DO 20 I = 1, NREC
         CALL X04BAF(NOUT,REC(I))
   20 CONTINUE
      RETURN
C
C     End of X04BAY.
C
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END

      SUBROUTINE X02ZAZ
C     MARK 16 REVISED. IER-1046 (JUN 1993).
C
C***********************************************************************
C
C     NAG version of the Stanford routine MCHPAR.
C     Sven Hammarling, NAG Central Office.
C
C     X02ZAZ sets machine parameters as follows:
C
C     WMACH(  1 ) = nbase  = base of floating-point arithmetic.
C     WMACH(  2 ) = ndigit = no. of base ( nbase ) digits in the
C                            mantissa
C     WMACH(  3 ) = eps    = relative machine accuracy. (X02AJF.)
C     WMACH(  4 ) = rteps  = sqrt( eps ).
C     WMACH(  5 ) = rmin   = small positive floating-point number whose
C                             reciprocal does not overflow.
C     WMACH(  6 ) = rtrmin = sqrt( rmin ).
C     WMACH(  7 ) = rmax   = 1/rmin
C     WMACH(  8 ) = rtrmax = sqrt( rmax ).
C     WMACH(  9 ) = undflw = 0 if underflow is not fatal, +ve otherwise.
C     WMACH( 10 ) = nin    = input  stream unit number. ( 5.)
C     WMACH( 11 ) = nout   = output stream unit number.
C                          = advisory message unit number. ( X04ABF.)
C     WMACH( 12 ) = nerr   = error    message unit number. ( X04AAF.)
C     WMACH( 13 )
C     WMACH( 14 )   Not currently used.
C     WMACH( 15 )
C
C     Note that constants that represent integers may hold a number just
C     less than the integer, so that the integer should be recovered by
C     adding, say, 0.25. e.g.
C
C     IBASE = WMACH( 1 ) + 0.25
C
C***********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D+0)
C     .. Arrays in Common ..
      DOUBLE PRECISION WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION EPS, RMAX, RMIN, UNDFLW
      INTEGER          NBASE, NDIGIT, NERR, NOUT
      LOGICAL          FIRST
C     .. External Functions ..
      DOUBLE PRECISION X02AJF, X02AMF
      INTEGER          X02BHF, X02BJF
      LOGICAL          X02DAF
      EXTERNAL         X02AJF, X02AMF, X02BHF, X02BJF, X02DAF
C     .. External Subroutines ..
      EXTERNAL         X04AAF, X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC        SQRT
C     .. Common blocks ..
      COMMON           /AX02ZA/WMACH
C     .. Save statement ..
      SAVE             /AX02ZA/, FIRST
C     .. Data statements ..
      DATA             FIRST/.TRUE./
C     .. Executable Statements ..
C
      IF (FIRST) THEN
         FIRST = .FALSE.
C
         IF (X02DAF(ZERO)) THEN
            UNDFLW = 1
         ELSE
            UNDFLW = 0
         END IF
         NBASE = X02BHF()
         NDIGIT = X02BJF()
         EPS = X02AJF()
         RMIN = X02AMF()
         RMAX = 1/RMIN
C
         WMACH(1) = NBASE
         WMACH(2) = NDIGIT
         WMACH(3) = EPS
         WMACH(4) = SQRT(EPS)
         WMACH(5) = RMIN
         WMACH(6) = SQRT(RMIN)
         WMACH(7) = RMAX
         WMACH(8) = SQRT(RMAX)
         WMACH(9) = UNDFLW
      END IF
      CALL X04ABF(0,NOUT)
      WMACH(10) = 5
      WMACH(11) = NOUT
      CALL X04AAF(0,NERR)
      WMACH(12) = NERR
      RETURN
C
C     End of  X02ZAZ. (MCHPAR)
C
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /1.11022302462516D-16 /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      subroutine e04uef(string)
c----------------------------------------------------------------------
c     e04uef  loads the option supplied in string into the relevant
c     element of iprmls, rprmls, iprmnp or rprmnp.

      character*(*)     string
c     .. scalars in common ..
      logical           newopt

      double precision wmach
      common/ cstmch /wmach(9)
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
      common/ cstmch /wmach(9)

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