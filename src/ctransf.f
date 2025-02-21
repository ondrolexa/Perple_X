c----------------------------------------------------------------------
c   ctransf is a program to read a vertex thermo-data file and
c   rewrite the data in a new file with transformed components.  
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i 

      character*8 name

      logical eof

      integer eos
      common/ cst303 /eos(k10)

      character specie*4
      integer ins, isp
      common/ cxt33 /isp,ins(nsp),specie(nsp)

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

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp
c----------------------------------------------------------------------- 
c                                 iam is a flag indicating the Perple_X program
      iam = 6
c                                 perplexwrap.f flags
      getInput = .true.
      sWarn = .false.
c                                 version info
      call vrsion (6)

      write (*,1000)
c                                 assign data files
      call sopen 

      write (*,1030)
c                                 Read THERMODYNAMIC DATA file (N2):
c                                 read the data base header
      call topn2 (5)
c                                 disable HSC conversion
      hscon = .false. 
c                                 mock pointers
      do i = 1, icmpn
         ic(i) = i
      end do
c                                 this is necessary for getphi to transform GH = G0
      icomp = icmpn
c                                 read and echo data cards with
c                                 component conversion
      do 

         call getphi (name,.true.,eof)

         names(k10) = name
c                                  7.1.8 archaic trap
         if (ieos.gt.200.and.ieos.lt.203) call errdbg (name//
     *       ' has an invalid EoS specification, execution terminated.')

         eos(k10) = ieos

         lct(k10) = ilam

         ltyp(k10) = jlam

         idis(k10) = idiso

         if (eof) exit

         if (ieos.eq.12.or.ieos.eq.14.or.ieos.eq.17) then 
            write (*,1010) name
            cycle
         end if 
c                                 output new data
         call outdat (n8,k10,0)

      end do 

      write (*,1020)

1000  format (//,'NO is the default answer to all Y/N prompts',/)
1010  format (//,'**warning ver000** ctransf cannot reformat CALPHAD ',
     *          'format data',/,'the data for ',a,' will not be ',
     *          'written to ctransf.dat',//)
1020  format (/,'The transformed dataset has been written to file: ',
     *          'ctransf.dat',/)
1030  format (/,'The transformed dataset will be written to file: ',
     *          'ctransf.dat',/)

      end
