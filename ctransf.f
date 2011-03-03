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

      integer idh2o,idco2,ikind,icmpn,icout,ieos
      double precision comp,tot
      common/ cst43 /comp(k0),tot,icout(k0),idh2o,idco2,
     *               ikind,icmpn,ieos

      integer ic
      common/ cst42 /ic(k0)
c----------------------------------------------------------------------
c                                 version info
      call vrsion

      write (*,1000)
c                                 assign data files
      call sopen (1)
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

         if (eof) exit
c                                 output new data
         call outdat (n8,k10,0)

      end do 
      
1000  format (//,'NO is the default answer to all Y/N prompts',/)

      end

      subroutine grxn (g)
c--------------------------------------------------------------------
c a dummy routine to allow rk to be linked with rlib.f
c--------------------------------------------------------------------
      implicit none
      double precision g
      g = g
      end