c----------------------------------------------------------------------
c   ctransf is a program to read a vertex thermo-data file and
c   rewrite the data in a new file with transformed components.  
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i 

      character*8 name

      logical eof

      character*8 names
      common/ cst8 /names(k1)

      integer eos
      common/ cst303 /eos(k10)

      integer ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),ikind,icmpn,ieos

      integer ic
      common/ cst42 /ic(k0)

      integer ilam,jlam,idiso,lamin,idsin
      double precision tm,td
      common/ cst202 /tm(m7,m6),td(m8),ilam,jlam,idiso,lamin,idsin

      integer ltyp,lct,lmda,idis
      common/ cst204 /ltyp(k10),lct(k10),lmda(k10),idis(k10)

      integer iam
      common/ cst4 /iam
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 6
c                                 version info
      call vrsion

      write (*,1000)
c                                 assign data files
      call sopen 
c                                 Read THERMODYNAMIC DATA file (N2):
c                                 read the data base header
      call topn2 (5)
c                                 mock pointers
      do i = 1, icmpn
         ic(i) = i
      end do
c                                 read and echo data cards with
c                                 component conversion
      do 

         call getphi (name,eof)

         names(k10) = name

         eos(k10) = ieos

         lct(k10) = ilam

         ltyp(k10) = jlam

         idis(k10) = idiso

         if (eof) exit

         if (ieos.eq.12.or.ieos.eq.14) then 
            write (*,1010) name
            cycle
         end if 
c                                 output new data
         call outdat (n8,k10,0)

      end do 
      
1000  format (//,'NO is the default answer to all Y/N prompts',/)
1010  format (//,'**warning ver000** ctransf cannot reformat CALPHAD ',
     *          'format data',/,'the data for ',a,' will not be ',
     *          'written to ctransf.dat',//)

      end

      subroutine grxn (g)
c--------------------------------------------------------------------
c a dummy routine to allow rk to be linked with rlib.f
c--------------------------------------------------------------------
      implicit none
      double precision g
      g = g
      end