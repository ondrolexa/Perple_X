      subroutine yclos1 (clamda,x,is,jphct,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement. this routine is only called as preparation
c for iterative refinement.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, i, j, k, is(*), idsol(k5), kdv(k19+1), nsol, 
     *        mpt, iam, id, is1, left, right, inc, jdsol(k5,k5), 
     *        kdsol(k5), max, is0, mcpd

      external ffirst

      logical solvus, quit

      double precision clamda(*), x(*),  slam(k19+1)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isoct
      common/ cst79 /isoct

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer ikp
      common/ cst61 /ikp(k1)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k5),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      npt = 0 
      mcpd = 0 
      nsol = 0
      inc = istct - 1
      is1 = isoct + 2 
      is0 = is1 - 1
      quit = .true.
c                                 solvus_tolerance_II, 0.25
      soltol = nopt(25)

      do i = 1, jphct

         if (is(i).ne.1) then 
c                                 make a list of found phases:
            id = i + inc
c                                 currently endmember compositions are not 
c                                 refined (this is probably a mistake, but 
c                                 seems to work fine), so use id > ipoint
c                                 to identify a refineable point

c                                 modified 3/4/2015 to refine endmember compositions. JADC
c           if (ikp(id).ne.0) then  #refine endmember change
            if (ikp(id).ne.0.and.id.gt.ipoint) then

               quit = .false.

               do j = 1, nsol
                  if (ikp(id).eq.idsol(j)) then 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = id
                     goto 10
                  end if 
               end do
c                                 new phase, add to list
               nsol = nsol + 1
               idsol(nsol) = ikp(id)
               jdsol(1,nsol) = id
               kdsol(nsol) = 1

            end if 
c                                 new point, add to list
10          npt = npt + 1
            jdv(npt) = i

         end if 

      end do

      do i = 1, is1
         slam(i) = 1d99
         kdv(i) = 0 
      end do 
c                                 perp 6.6.3, make a list of metastable
c                                 phases, this list includes two compounds
c                                 and the least metastable composition of
c                                 each solution.      
      do 20 i = 1, jphct

         if (is(i).ne.1) cycle 

         id = i + inc 
         iam = ikp(id)
c                                modified to allow endmembers, Mar 4, 2015. JADC
c        if (iam.ne.0.and.id.gt.0*ipoint) then   #refine endmember change
         if (iam.ne.0.and.id.gt.ipoint) then  

            if (clamda(i).lt.slam(iam)) then
c                                the composition is more stable
c                                than the previous composition 
c                                of the solution. check if it's 
c                                one of the stable solutions                               
               do j = 1, nsol
                  if (iam.eq.idsol(j)) then
c                                it's already stable, only accept
c                                it if its further than the solvus
c                                tolerance from any of the stable
c                                compositions.
                     do k = 1, kdsol(j)
                        if (.not.solvus(jdsol(k,j),id,iam)) goto 20
                     end do 
                  end if
               end do
c                                the composition is acceptable
               slam(iam) = clamda(i)
               kdv(iam) = i

            end if

         else 
c                                a compound, save lowest 2, changed from 1
c                                march 4, 2015, JADC
            if (clamda(i).lt.slam(is0)) then 
c                                put the old min into is1
               kdv(is1) = kdv(is0)
               slam(is1) = slam(is0)
c                                and save the current in is1
               kdv(is0) = i
               slam(is0) = clamda(i)

               mcpd = mcpd + 1

            end if 
 
         end if 

20    continue

      if (mcpd.gt.2) mcpd = 2
c                                 load the metastable points into
c                                 kdv, the mcpd metastable points
c                                 will be last in this list
      mpt = 0 

      do i = 1, is1

         if (kdv(i).eq.0) cycle
         mpt = mpt + 1
         kdv(mpt) = kdv(i)
         slam(mpt) = slam(i)

      end do 

      if (mpt-mcpd.le.iopt(12)) then 
c                                 less metastable refinement points than
c                                 iopt(12)
            max = mpt

      else 
c                                 sort the metastable points to
c                                 find the most stable iopt(12) points
         left = 1
         right = mpt-mcpd
         max = iopt(12)

         call ffirst (slam,kdv,left,right,max,k19+1,ffirst)

      end if 
 
      do i = 1, max

         jdv(npt+i) = kdv(i)
c                                 a metastable solution to be refined
c        if (kdv(i)+inc.gt.ipoint) quit = .false.
         if (ikp(kdv(i)+inc).ne.0) quit = .false.
      end do
c                                 and load the compounds
      do i = 1, mcpd
         jdv(npt+max+i) = kdv(mpt-mcpd+i)
      end do

      if (quit) then 
c                                 zero mode filter and 
c                                 save amounts for final processing
         mpt = npt
         npt = 0 

         do i = 1, mpt
            if (x(jdv(i)).lt.nopt(9)) cycle 
            npt = npt + 1
            jdv(npt) = jdv(i)
            amt(npt) = x(jdv(i)) 
         end do 

      else 

         npt = npt + max + mcpd
c                                 sort the phases, why? don't know, but it's 
c                                 necessary
         call sortin 

      end if 

      end 



      subroutine yclos3 (clamda,x,is,iter,opt)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. yclos3 serves the same 
c function as yclos2, but adopts the strategy of yclos1: a maximum of
c iopt(12) + 2 refinement points will be selected, these include the 
c 2 least metastable compounds and 1 least stable composition for 
c each solution phase.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, is(*), id, jmin(k19), opt, mpt, iter, tic

      double precision clamda(*), clam(k19), x(*)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k5),ctotal,jdv(k19),npt,fulrnk

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer kkp, np, ncpd, ntot
      double precision cp3, amt
      common/ cxt15 /cp3(k0,k5),amt(k5),kkp(k5),np,ncpd,ntot

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      save tic
      data tic/0/
c----------------------------------------------------------------------
c                                 npt is the number of points refined
c                                 from the previous cycle, opt is the
c                                 number of points in the original 
c                                 solution.
      do i = 1, npt
         jmin(i) = 0 
         clam(i) = 1d99
      end do 

      npt = 0

      do i = 1, jphct
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 check the stability of all points 
         if (is(i).ne.1.and.x(i).gt.0d0) then 
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i

         else if (clamda(i).lt.clam(id)) then 
c                                 find the nearest phase           
            jmin(id) = i
            clam(id) = clamda(i)

         end if 

      end do 

      if (iter.le.iopt(10)) then
c                                 at this point there is one metastable refinement
c                                 point for each original point.



c                                 if not done iterating, add the metastable
c                                 phases
         do i = 1, opt
            if (jmin(i).eq.0) cycle 
            npt = npt + 1
            jdv(npt) = jmin(i)
         end do
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin
c                                 make a pointer to the original refinement 
c                                 point
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do 

      else  
  
         mpt = npt 
         npt = 0  
c                                 check zero modes the amounts
         do i = 1, mpt

            if (x(jdv(i)).ge.nopt(9)) then 
               npt = npt + 1
               amt(npt) = x(jdv(i))
               jdv(npt) = jdv(i)
            else if (lopt(13).and.x(jdv(i)).lt.-nopt(9)
     *                             .and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

         if (npt.ne.icp) then
 
            fulrnk = .false.
         
         else 

            fulrnk = .true.

         end if 


      end if 

      end