      call junk3
      end 

      subroutine junk 

c C:\Program Files (x86)\VNI\imsl\fnl700\winin111e64\include\dll
      INCLUDE 'link_fnl_shared.h'

      USE LCONF_INT

      IMPLICIT NONE
! Declaration of variables
      INTEGER NCON, NEQ, NVAR, LDA, NACT, IPRINT, INFO, MCON, MVAR,DCON,
     *        DVAR,DVDC
      PARAMETER (NCON=2, NEQ=0, NVAR=3,MCON=10,MVAR=10,DVAR=MVAR-NVAR,
     *           DCON = MCON - NCON, DVDC=DVAR*DCON)
!
      INTEGER MAXFCN, NOUT, IACT(NCON+2*NVAR)
      double precision A(MCON,MVAR), ACC, B(MCON), OBJ, ALAMDA(MVAR), 
     *SOL(MVAR), XGUESS(MVAR), XLB(MVAR), XUB(MVAR), gfinal,
     * WK(MVAR**2+11*MVAR+MCON)
      EXTERNAL FCN
!
! Set values for the following problem.
!
! Min -X(1)*X(2)*X(3)
!
! -X(1) - 2*X(2) - 2*X(3) .LE. 0
! X(1) + 2*X(2) + 2*X(3) .LE. 72
!
! 0 .LE. X(1) .LE. 20
! 0 .LE. X(2) .LE. 11
! 0 .LE. X(3) .LE. 42
!
      DATA B/0.0, 72.0,DCON*0./
      DATA XLB/MVAR*0.0/, XUB/20.0, 11.0, 42.0,DVAR*0./, 
     *     XGUESS/MVAR*10.0/
      DATA ACC/0.0/, MAXFCN/400/

      nout = 6

      CALL UMACH (2, NOUT)


      LDA = NVAR
      ACC = 1e-8
      iprint = -1

      A(1,1) = -1
      A(1,2) = -2
      A(1,3) = -2

      A(2,1) = 1
      A(2,2) = 2
      A(2,3) = 2

c     CALL DLCONF (FCN, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
c    *XGUESS, ACC, MAXFCN, SOL, OBJ, NACT, IACT, ALAMDA)

c     CALL dL2ONF (FCN, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
c    *XGUESS, ACC, MAXFCN, SOL, gfinal, NACT, IACT, ALAMDA, IPRINT, 
c    *  INFO, WK)

      CALL lconf (FCN, NEQ, A, B, XLB, XUB, XGUESS, XGUESS=XGUESS, 
     *            MAXFCN=MAXFCN, ACC=1d-8, OBJ= gfinal, NACT= NACT,
     *            IACT = IACT, ALAMDA = ALAMDA)

      WRITE (NOUT,99998) 'Solution:'
      WRITE (NOUT,99999) XGUESS
      WRITE (NOUT,99998) 'Function value at solution:'
      WRITE (NOUT,99999) gfinal
      WRITE (NOUT,99998) 'Number of function evaluations:', MAXFCN

99998 FORMAT (//, ' ', A, I4)
99999 FORMAT (1X, 5F16.6)
      END
!
      SUBROUTINE FCN (N, X, F)
      INTEGER N
      double precision X(*), F
!
      F = -(X(1))*(X(2))*(X(3))

      write (*,*) x(1:3)
      RETURN
      END


      subroutine minime (ids)
c-----------------------------------------------------------------------
c     INCLUDE 'link_fnl_shared.h'

c     USE LCONF_INT

      implicit none

      integer h4,h9,m0,m1,m2,m3,m4,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15
      integer m16,m17,m18,h0
      integer msp,mst,mdim,ms1
      parameter (h4=5,h9=30,h0=h9+1)
      parameter (m0=12,m1=60,m2=8,m3=3,m4=96,m6=6,m7=15,m8=9,m9=10,
     *           m10=6,m11=11,m12=4,m13=8,m14=2,m15=85,m16=6,m17=5,
     *           m18=6)
      parameter (mst=4,mdim=8,msp=mdim+6,ms1=msp-1)

      integer ids, i, j, k, ncon, ocon, nact, iact(2*m4+2*m11),
     *        maxfcn, neq, gfinal

      external gsol2

      double precision az(2*m11,m4), bz(2*m11), lambda(m4), plb(2*m11),
     *                 pub(2*m11)

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)
c in perp.h
      integer msite, zsp
      double precision zcoef, zmult
      common/ cxt1n /zcoef(0:m0,m11,m10,h0),zmult(h0,m10),
     *               msite(h0),zsp(h0,m10)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision dcoef, scoef
      common/ cxt1r /dcoef(0:m0,m11,m10,h9),scoef(m4,h9)
c-----------------------------------------------------------------------
c                                 create constraint matrix z
c                                 --------------------------
c                                 initialize bounds
      do i = 1, nstot(ids)
         pub(i) = 1d1
         plb(i) = -1d1
      end do 
c                                 1 linear constraint, closure:
      neq = 1

      do k = 1, nstot(ids)

         az(neq,k) = 1d0

      end do

      bz(ncon) = 1d0
c                                 ncon - neq linear inequalities,
c                                 to be counted:
      ncon = neq
c                                 for each site
      do i = 1, msite(ids)
c                                 for each species
         do j = 1, zsp(ids,i)
c                                 initial az, bz
            ncon = ncon + 1
c                                 both Temkin and non-Temkin have
c                                 -Az*p <= 0 constraints:
            bz(ncon) = dcoef(0,j,i,ids)

            do k = 1, nstot(ids)

               az(ncon,k) = 0d0

            end do

            do k = 1, lterm(j,i,ids)
 
               az(ncon,ksub(k,j,i,ids)) = -dcoef(k,j,i,ids)

            end do
c                                 non-Temkin have the Az*p <= 1 constraint
            if (zmult(ids,i).ne.0d0) then

               ocon = ncon

               ncon = ncon + 1

               bz(ncon) = -bz(ocon) + 1d0

               do k = 1, nstot(ids)

                  az(ncon,k) = az(ocon,k)

               end do

            end if

         end do 

      end do

      maxfcn = 400

c     CALL lconf (gsol2,neq,az,bz,plb,pub,pp,xguess=pp,
c    *            maxfcn = maxfcn, acc = 1d-8, obj = gfinal,
c    *            nact = nact, iact = iact, alamda = lambda)

      end 

      subroutine gsol2 (nstot,pp,gval) 
      integer nstot
      double precision pp(*), gval


      end

      subroutine junk2 

      INCLUDE 'link_fnl_shared.h'

      USE LCONF_INT

      IMPLICIT NONE

      include 'perplex_parameters.h'

      integer ids, i, j, k, ncon, ocon, nact, iact(2*m4+2*m11),
     *        maxfcn, neq


! Declaration of variables
      INTEGER NVAR, LDA, IPRINT, INFO, MCON, MVAR

      PARAMETER (MCON=10,MVAR=10)
!
      INTEGER NOUT
      double precision A(MCON,MVAR), ACC, B(MCON), OBJ, ALAMDA(MVAR), 
     *SOL(MVAR), XGUESS(MVAR), XLB(MVAR), XUB(MVAR), gfinal,
     * WK(MVAR**2+11*MVAR+MCON)
      EXTERNAL fcn2
!
! Set values for the following problem.
!
      DATA ACC/0.0/, MAXFCN/400/

      nout = 6


      ACC = 1e-8
      iprint = -1
c                                 closure
      neq = 1
      ncon = neq + 4
      nvar = 2
      LDA = NVAR

      xlb(1) = 1d-16
      xlb(2) = 1d-16
      xub(1) = 1d0 - xlb(1)
      xub(2) = 1d0 - xlb(2)

      A(1,1) = 1
      A(1,2) = 1
      B(1)   = 1d0
c                                  x1<=1
      A(2,1) = 1
      A(2,2) = 0
      B(2)   = xub(1)
c                                  -x1<=0
      A(3,1) = -1
      A(3,2) = 0
      B(2)   = -xlb(1)
c                                  x2<=1
      A(4,1) = 0
      A(4,2) = 1
      B(4)   = xub(2)
c                                  -x2<=0
      A(5,1) = 0
      A(5,2) = -1
      B(5)   = -xlb(2)

 
      xguess(1) = 0.1
      xguess(2) = 1d0 - xguess(1)

      write (*,*) 'peep'

c     CALL DLCONF (FCN, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
c    *XGUESS, ACC, MAXFCN, SOL, OBJ, NACT, IACT, ALAMDA)

      CALL DL2ONF (FCN2, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
     *XGUESS, ACC, MAXFCN, SOL, OBJ, NACT, IACT, ALAMDA, IPRINT, INFO, 
     *WK)

      CALL lconf (fcn2, NEQ, A, B, XLB, XUB, XGUESS, XGUESS=XGUESS, 
     *            MAXFCN=MAXFCN, ACC=1d-8, OBJ= gfinal, NACT= NACT,
     *            IACT = IACT, ALAMDA = ALAMDA)

      WRITE (NOUT,99998) 'Solution:'
      WRITE (NOUT,99999) XGUESS
      WRITE (NOUT,99998) 'Function value at solution:'
      WRITE (NOUT,99999) gfinal
      WRITE (NOUT,99998) 'Number of function evaluations:', MAXFCN
      STOP
99998 FORMAT (//, ' ', A, I4)
99999 FORMAT (1X, 5F16.6)
      END

      subroutine junk3 

      INCLUDE 'link_fnl_shared.h'

      USE LCONF_INT

      IMPLICIT NONE

      include 'perplex_parameters.h'

      integer ids, i, j, k, ncon, ocon, nact, iact(2*m4+2*m11),
     *        maxfcn, neq


! Declaration of variables
      INTEGER NVAR, LDA, IPRINT, INFO, MCON, MVAR

      PARAMETER (MCON=10,MVAR=10)
!
      INTEGER NOUT
      double precision A(MCON,MVAR), ACC, B(MCON), OBJ, ALAMDA(MVAR), 
     *SOL(MVAR), XGUESS(MVAR), XLB(MVAR), XUB(MVAR), gfinal,
     * WK(MVAR**2+11*MVAR+MCON)
      EXTERNAL fcn3
!
! Set values for the following problem.
!
      DATA ACC/0.0/, MAXFCN/400/

      nout = 6


      ACC = 1e-8
      iprint = -1
c                                 closure
      neq = 0
      ncon = 2
      nvar = 1
      LDA = NVAR

      xlb(1) = 1d-7
      xub(1) = 1d0 - 1d-7
c                                  x1<=1
      A(1,1) = 1
      A(1,2) = 0
      B(1)   = 1
c                                  -x1<=0
      A(2,1) = -1
      A(2,2) = 0
      B(2)   = 0

      xguess(1) = 0.1

      write (*,*) 'peep'

c     CALL DLCONF (FCN, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
c    *XGUESS, ACC, MAXFCN, SOL, OBJ, NACT, IACT, ALAMDA)

c     CALL DL2ONF (FCN2, NVAR, NCON, NEQ, A, LDA, B, XLB, XUB,
c    *XGUESS, ACC, MAXFCN, SOL, OBJ, NACT, IACT, ALAMDA, IPRINT, INFO, 
c    *WK)

      CALL lconf (fcn3, NEQ, A, B, XLB, XUB, XGUESS, XGUESS=XGUESS, 
     *            MAXFCN=MAXFCN, ACC=1d-8, OBJ= gfinal, NACT= NACT,
     *            IACT = IACT, ALAMDA = ALAMDA)

      WRITE (NOUT,99998) 'Solution:'
      WRITE (NOUT,99999) XGUESS
      WRITE (NOUT,99998) 'Function value at solution:'
      WRITE (NOUT,99999) gfinal
      WRITE (NOUT,99998) 'Number of function evaluations:', MAXFCN

99998 FORMAT (//, ' ', A, I4)
99999 FORMAT (1X, 5F16.6)
      END

      SUBROUTINE FCN2 (N, X, F)
      INTEGER N
      double precision X(*), F
!
      F = x(1)*dlog(x(1)) + x(2)*dlog(x(2))

      write (*,*) x(1:2), x(1) + x(2)
      write (*,*) f
      RETURN
      END

      SUBROUTINE FCN3 (N, X, F)
      INTEGER N
      double precision X(*), F

      x(2) = 1 - x(1)
!
      F = x(1)*dlog(x(1)) + x(2)*dlog(x(2))

      write (*,*) x(1:2), x(1) + x(2)
      write (*,*) f
      RETURN
      END