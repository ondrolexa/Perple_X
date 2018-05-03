      subroutine kill01 (site,im)
c---------------------------------------------------------------------
c reform - counts the number of species that can be respresented for a 
c solution given the present endmembers.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      character*10 sname

      logical first

      integer kill,ikill,jkill,kill1,i,j,kosp(mst,msp),kill2,
     *        k,im,idsp,ksp(mst),site

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer ostot
      common/ junk /ostot

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf
c----------------------------------------------------------------------


c                                 count the number of species
c                                 missing on site
      ksp = 0

      do
         do i = 1, isp(site)
            if (kdsol(istot+i).eq.0) then 
               ksp = ksp + 1
               call killsp (site,i)
               exit 
            end if 
         end do
         if (i.gt.isp(site)) exit
      end do 

      end 
