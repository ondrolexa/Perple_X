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
               call nkillsp (site,i)
               exit 
            end if 
         end do
         if (i.gt.isp(site)) exit
      end do 

      end 


      subroutine killsp (ikill,jkill)
c---------------------------------------------------------------------
c killsp - eliminates species jkill from site ikill in a solution model
c and reformulates the model accordingly
c---------------------------------------------------------------------
      implicit none
  
      include 'perplex_parameters.h'

      logical skip, bad, dead

      integer jsp,jtic,morder,
     *        i,j,ikill,jkill,kill,kdep,jdqf,ktic,jold,
     *        iwas(m4),i2ni(m4),kwas(m4),
     *        k,l,itic,ijkill(m4),
     *        j2oj(msp),j2nj(msp),i2oi(m4),maxord,mord
c                                 dqf variables
      integer indq,idqf
      double precision dqf
      common/ cst222 /dqf(m3,m4),indq(m4),idqf

      logical depend,laar,order,fluid,macro,recip
      common/ cst160 /depend,laar,order,fluid,macro,recip
c                                 local input variables
      integer iddeps,norder,nr
      double precision depnu,denth
      common/ cst141 /depnu(j4,j3),denth(j3,3),iddeps(j4,j3),norder,
     *                nr(j3)

      integer jsmod
      double precision vlaar
      common/ cst221 /vlaar(m3,m4),jsmod

      integer jmsol,kdsol
      common/ cst142 /jmsol(m4,mst),kdsol(m4)

      integer nsub,nttyp,nterm,nspm1,nsite
      double precision acoef,smult,a0
      common/ cst107 /a0(m10,m11),acoef(m10,m11,m0),smult(m10),
     *      nsite,nspm1(m10),nterm(m10,m11),nsub(m10,m11,m0,m12),
     *      nttyp(m10,m11,m0)

      logical stck, norf
      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),reach,iend(m4),
     *      isub(m1,m2,2),imd(msp,mst),insp(m4),ist(mst),isp(mst),
     *      rkord(m18),isite,iterm,iord,istot,jstot,kstot,xtyp,stck,norf

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      integer iorig,jnsp,iy2p
      common / cst159 /iorig(m4),jnsp(m4),iy2p(m4)

      integer mdep,idep,jdep,ndph
      double precision nu,y2p
      common/ cst146 /nu(m15,j4),y2p(m4,m15),mdep,jdep(m15),
     *                idep(m15,j4),ndph(m15)

      double precision yin
      common/ cst50 /yin(ms1,mst)

      integer limn,limt,limid,jimid,jimt
      double precision limc,jimc
      common/ cxt30 /limc(j6+2,j5,j3),limid(m0,j5,j3),jimid(j3,j5,j3),
     *               limn(j3),limt(j5,j3),jimc(j3,j5,j3),jimt(j5,j3)
c----------------------------------------------------------------------
      do i = 1, isite

         if (i.ne.ikill) then
c                                 nothing happens  
            do j = 1, isp(i)  
              j2oj(j) = j 
            end do
         else 
c                                 on a site where we kill a species
            jsp = isp(i) - 1
c                                 should also check and store subdivsion
c                                 model here (i.e., some ternary model
c                                 choices are not valid as binary):

c                                 make a pointer from the new species
c                                 index to the old index
            jtic = 0 

            do j = 1, isp(i)
               if (j.ne.jkill) then 
                  jtic = jtic + 1 
c                              pointer from new j to old j
                  j2oj(jtic) = j
c                              pointer from old j to new j
                  j2nj(j) = jtic
               end if 
            end do
c                              now reload
            isp(i) = jsp

            if (jsp.gt.1) then   
c                              now shift subdivision ranges
               do j = 1, jsp - 1
                  xmn(i,j) = xmn(i,j2oj(j))
                  xmx(i,j) = xmx(i,j2oj(j))
                  xnc(i,j) = xnc(i,j2oj(j))
                  imd(j,i) = imd(j2oj(j),i)

                  if (imd(j,i).gt.0) yin(j,i) = yin(j2oj(j),i)

               end do
            else 
               xmn(i,j) = 1d0
               xmx(i,j) = 1d0
               xnc(i,j) = 1d0                 
            end if
         end if 
      end do 

      kdep = 0 

      do i = 1, istot

         if (depend.and.kdsol(i).eq.-2) then 
c                                create an array which gives the
c                                original locations of the dependent
c                                endmembers, need this to be able to
c                                reorder the y2p array:
            kdep = kdep + 1
            iwas(kdep) = i 
         end if 
c                                 kill endmembers with the species
c                                 to be deleted:
         if (jmsol(i,ikill).eq.jkill) kdsol(i) = -3

      end do 
c                                 check for dependent endmembers
      call redep (-3)
c                                 at this point all ordered endmembers to 
c                                 be killed are flagged by kdsol = -3.

c                                 now check the ordered species
      morder = 0 

      if (order) then 
c                                 first check if the ordered endmember
c                                 may be stable 
         do k = 1, norder      
c                                 check if a missing constituent 
            bad = .false.

            do j = 1, nr(k)
               if (kdsol(iddeps(j,k)).eq.-3) then 
                  bad = .true.
                  exit 
               end if 
            end do

            if (bad) then 
c                                 add species to the kill list
               kdsol(istot+k) = -3
            
            else 

               morder = morder + 1
               kdsol(istot+k) = -1
               kwas(morder) = k

            end if 

         end do  

      end if 
c                                figure out which dependent endmembers have
c                                been killed:
      if (depend) then 

         if (kdep.ne.mdep) call error (54,dqf(1,1),mdep,'KILLSP')

         mdep = 0 

         do i = 1, kdep
            if (kdsol(iwas(i)).ne.-3) then 
               mdep = mdep + 1
c                                 iwas is now the original index of the 
c                                 dependent endmember, and mdep is the reset
c                                 value of the dependent endmember counter
               iwas(mdep) = i                
            end if
         end do 

      end if 
 
      itic = 0 
      jtic = 0
      kill = 0 
      
      do i = 1, istot + norder 

         if (kdsol(i).ge.-2) then 
c                                 replacement for istot (itic)
            itic = itic + 1
c                                 pointers from new to old endmember index (i2oi)
            i2oi(itic) = i
c                                 pointers from new to old endmember index (i2ni)
            i2ni(i) = itic                                
c                                 pointer to original species index
            iorig(itic) = iorig(i)
c                                 number of missing endmembers (jtic)
            if (kdsol(i).eq.0) jtic = jtic + 1
c                                reset the kdsol array
            kdsol(itic) = kdsol(i)

            if (i.gt.istot) cycle
c                                 reset the species pointers (jmsol) 
            do j = 1, isite
               if (j.eq.ikill) then
                  jmsol(itic,j) = j2nj(jmsol(i,j))
               else
                  jmsol(itic,j) = jmsol(i,j)
               end if 
            end do
         else 
c                                 kill records the killed endmembers
            kill = kill + 1
            ijkill(kill) = i

         end if 
      end do  
c                                reset total and present counters
      istot = itic - morder 
c            
      jstot = itic - jtic - morder
c                                --------------------------------------
c                                excess terms:
      itic = 0 
      maxord = 0 
    
      do i = 1, iterm 
c                                check for forbidden terms (i.e., terms
c                                with a missing endmember
         skip = .false. 
c                                 macroscopic formulation
         do j = 1, kill
c                                 check if subscript points to a killed 
c                                 endmember
            do k = 1, iord
               if (isub(i,k,1).eq.0) then
                  cycle 
               else if (isub(i,k,1).eq.ijkill(j)) then
                  skip = .true.
                  exit 
               end if 
            end do 

            if (skip) exit

         end do 

         if (skip) cycle 
c                               the term is acceptable
         itic = itic + 1

         mord = iord

         do j = 1, iord
            if (isub(i,j,1).eq.0) then
               isub(itic,j,1) = 0 
            else 
               isub(itic,j,1) = i2ni(isub(i,j,1))
            end if 
         end do 

         if (xtyp.eq.0) then 
c                                save the coefficient
            do j = 1, m3
               wg(itic,j) = wg(i,j)
            end do 
c                                find highest order term
            if (mord.gt.maxord) maxord = mord

         else
c                                 redlich kistler
            rkord(itic) = rkord(i)

            do j = 1, rkord(itic)
               do k = 1, m16
                  wk(k,j,itic) = wk(k,j,i)
               end do
            end do 

            maxord = 2

         end if 

      end do     
c                                reset counters, iord is not reset
      iterm = itic
      iord = maxord                         
c                                --------------------------------------
c                                van laar volume functions
      if (laar) then
         do i = 1, istot + morder 
            do j = 1, m3
               vlaar(j,i) = vlaar(j,i2oi(i))
            end do 
         end do  
      end if 
c                                 --------------------------------------
c                                 dqf corrections, this is sloppy since
c                                 uses istot instead of kstot
      if (idqf.gt.0) then 

         jdqf = 0 
c                                 check if a retained species has a dqf
c                                 correction
         do j = 1, idqf
c                                 the itoi index must be in the inner loop
c                                 in case the values of indq are not sequential
            do i = 1, istot
               if (indq(j).eq.i2oi(i)) then 
c                                 found a dqf'd endmember
                  jdqf = jdqf + 1
                  indq(jdqf) = i
                  do k = 1, m3
                     dqf(k,jdqf) = dqf(k,j)
                  end do 
                  exit 
               end if
            end do 

            if (jdqf.eq.idqf) exit
 
         end do 

         idqf = jdqf 

      end if 
c                                 --------------------------------------
c                                 configurational entropy model

c                                 site fractions as a function of bulk
c                                 y's and dependent species y:
      do i = 1, nsite 
c                                 for each species, read function to define 
c                                 the site fraction of the species and eliminate
c                                 killed species

c                                 species counter is incremented in advance
c                                 and must be decremented before saving the 
c                                 final value:
         jtic = 1 

         do j = 1, nspm1(i)

            ktic = 0 
c                                 for each term:
            do k = 1, nterm(i,j)
c                                 macroscopic formulation:
c                                 note: 4th index (nttyp) is only used
c                                 for bragg-williams models.
               dead = .false.
               do l = 1, kill 
                  if (nsub(i,j,k,1).eq.ijkill(l)) then 
                     dead = .true.
                     exit 
                  end if 
               end do 

               if (.not.dead) then
c                                 the term has survived (and therefore 
c                                 also the species):
c                                 don't save nttyp since this is always
c                                 1 for the macroscopic formulation
                  ktic = ktic + 1
c                                 but my dear peanut brained friend, do 
c                                 not forget to move the pointer:
                  nsub(i,jtic,ktic,1) = i2ni(nsub(i,j,k,1)) 
                  acoef(i,jtic,ktic) = acoef(i,j,k)
               end if  
            end do 
c                                 ktic is the number of terms representing
c                                 the jth species, we won't count species 
c                                 with no terms because the endmember configurational
c                                 entropy is assumed to be implicit. 
         if (ktic.gt.0) then 
c                                 increment the species counter
            nterm(i,jtic) = ktic 
            a0(i,jtic) = a0(i,j)
            jtic = jtic + 1
         end if 

      end do

      nspm1(i) = jtic - 1

      end do 
c                                 ---------------------------------------
c                                 ordered species:
      if (order) then 

         norder = morder

         if (morder.eq.0) then 
c                                 there are no ordered species left
            order = .false.

            if (depend) then

               jsmod = 7

            else if (jsmod.eq.27) then 
c                                 special case, green et al 2016 melt model
c                                 converts to normal HP melt model
               jsmod = 24

            else 
c                                 why jsmod = 2?
               jsmod = 2

            end if 

         else 
c                                 shift the ordered species pointers 
c                                 and data to eliminate kill ordered
c                                 species.    
            do j = 1, morder 

               jold = kwas(j)

               do i = 1, 3
                  denth(j,i) = denth(jold,i)
               end do 

               nr(j) = nr(jold)

               do i = 1, nr(j)
                  iddeps(i,j) = i2ni(iddeps(i,jold))
                  depnu(i,j) = depnu(i,jold)
               end do

               itic = 1 

               do i = 1, limn(jold)
c                                 eliminate absent species from 
c                                 stoichiometric p0 limits
                  ktic = 0 

                  if (limt(i,jold).gt.0) then 
   
                     do k = 1, limt(i,jold)

                        skip = .false.

                        do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                           if (limid(k,i,jold).eq.ijkill(l)) then
                              skip = .true.
                              exit 
                           end if
  
                        end do  

                        if (skip) cycle

                        ktic = ktic + 1

                        limid(ktic,itic,j) = i2ni(limid(k,i,jold))
                        limc(ktic,itic,j) = limc(k,i,jold)

                     end do 

                     if (ktic.eq.0) cycle

                     limt(itic,j) = ktic

                  else 
c                                 constant bounds
                     limt(itic,j) = -1
                     k = 1

                  end if   

                  limc(ktic+1,itic,j) = limc(k,i,jold)
                  limc(ktic+2,itic,j) = limc(k+1,i,jold)
 
c                                 now check the p terms, this assumes
c                                 there are no p terms if there are no
c                                 p0 terms, which maynot be true?
                  ktic = 0 

                  do k = 1, jimt(i,jold)

                     skip = .false.

                     do l = 1, kill
c                                 check if limid points to a killed 
c                                 endmember
                        if (jimid(k,i,jold).eq.ijkill(l)) then
                           skip = .true.
                           exit 
                        end if
  
                     end do  

                     if (skip) cycle

                     ktic = ktic + 1

                     jimid(ktic,itic,j) = i2ni(jimid(k,i,jold))
                     jimc(ktic,itic,j) = jimc(k,i,jold)

                  end do 

                  jimt(itic,j) = ktic

                  itic = itic + 1

               end do 

               limn(j) = itic - 1

            end do 

         end if  

      end if 
c                                 --------------------------------------
c                                 dependent endmember properties, the
      if (depend) then 
c                                 dependent endmembers have been reordered
c                                 in redep, but are still expressed in
c                                 terms of the old indices, so reset the 
c                                 indices:
         do i = 1, mdep
            jdep(i) = i2ni(jdep(i))
            do j = 1, ndph(i)
               idep(i,j) = i2ni(idep(i,j))
            end do 
         end do 

      end if 

      end