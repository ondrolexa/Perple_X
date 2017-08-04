
c integrate an arbitrary function f between limits a..b

      implicit none
      
      INTEGER          LW, LIW, ifail, nbaf
      PARAMETER        (LW=40000,LIW=LW/4)

      double precision ABSERR, EPSABS, epsr, vol, W(LW)
      double precision xmin,xmax,f,ans

      INTEGER          IW(LIW)

      EXTERNAL         f

      common/ parm2 /nbaf
      common/ parm3 /epsr
c----------------------------------------------------------------------
      nbaf = 32
      epsr = 1d-3
      EPSABS = 0d0
      IFAIL = -1

      do 

         write (*,*) 'enter xmin..xmax'
         read (*,*) xmin, xmax

         CALL D01AJF(f,xmin,xmax,EPSABS,epsr,ans,ABSERR,
     *               W,LW,IW,LIW,IFAIL)


         IF (IFAIL.NE.0) WRITE (*,*) 'vol IFAIL = ', IFAIL

         write (*,*) ans

      end do 

      end 



      double precision FUNCTION f (x)
c----------------------------------------------------------------------
      implicit none 
      double precision x,c0,tmax,c,a,cp
c----------------------------------------------------------------------

      c0 = 80d0
      tmax = 4d3
      c = 3d5
      a = 0.6d0

      cp = c0*( 1d0 + a*exp( -(x-tmax)**2/2d0/c))
	f = cp/x

      end 
