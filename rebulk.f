      subroutine rebulk (static)
c----------------------------------------------------------------------
c upon successful completion of an optimization with either static or
c dynamic pseudocompounds rebulk:
c     1) loads the generic arrays cp3, cptot, ctot3, kkp and x3
c     2) computes the amounts of saturated component phases.
c     3) computes the dependent potentials.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k, tictoc

      logical static

      double precision c(k5)

      integer iff,idss,ifug,ifyn,isyn
      common/ cst10  /iff(2),idss(h5),ifug,ifyn,isyn
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,ctot3
      common/ cxt15 /cp3(k5,k5),ctot3(k5),kkp(k5),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k5),ctotal,jdv(k19),npt,fulrnk

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop

      save tictoc
      data tictoc/0/
c----------------------------------------------------------------------
      do i = 1, npt

         if (static) then 

            id = iam(jdv(i))
c                                 set identifier flag
            if (i.le.ipoint) then
               kkp(i) = -id
            else  
               kkp(i) = ikp(id)
            end if 

            cptot(i) = ctot(id)
           
            do j = 1, jbulk
               if (j.gt.icp.and.usv) exit
               cp3(j,i) = cp(j,id)
            end do
c                                 set the x3 array
            if (ikp(id).ne.0) call setx3 (i,id,ikp(id))

            if (usv) then 
c                                 usv calculations are currently
c                                 not set up for solutions and use
c                                 only static optimization
               call pt4sv (jam(jdv(i)))

               call getusv (u,cp3(icp1,i),cp3(jbulk,i),id,bad)

            end if 

         else 
c                                 getcmp assigns cp3, cptot, x3, and kkp
            call getcmp (i,jdv(i),jkp(jdv(i)))

         end if          
c                                 convert normalized amounts to molar 
c                                 amounts
         amt(i) = ctotal*amt(i)/cptot(i)

      end do 

      if (jbulk.gt.icp) then  
c                                 get the amounts of the saturated phases:
         do i = jbulk-icp+1, jbulk
c                                 initialize bulk                                 
            c(i) = cblk(i)
c                                 subtract the amount of the component in the 
c                                 phases in the thermodynamic c-space
            do j = 1, npt 
               c(i) = c(i)- amt(j)*cp3(i,j)
            end do 

         end do 

         do i = jbulk, jbulk-icp+1, -1
c                                  cycle through the saturated phases
            npt = npt + 1
            id = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(id).eq.0) then 
               kkp(npt) = -id
            else 
               write (*,1000) 
               stop 
            end if 
c                                  amount of the staurated phase
            amt(npt) = c(i)/cp(i,id)
c                                  warn on undersaturation
            if (amt(npt).lt.nopt(9)) then 
               if (amt(npt).lt.-nopt(9).and.tictoc.lt.5) 
     *                             call warn (2,amt(npt),i,'REBULK')
               npt = npt - 1
               tictoc = tictoc + 1
               exit
            end if 
c                                  remove the saturated phase from 
c                                  the bulk composition.
            do j = icp+1, i - 1
               c(j) = c(j) - amt(npt)*cp(j,id)
            end do 
c                                  load the saturated phase composition 
            do j = 1, jbulk
               cp3(j,i) = cp(j,id)
            end do           

         end do

      end if 

      ntot = npt

      if (usv.or.jpot.eq.0) then
c                                 compute chemical potentials
         if (npt.ne.jbulk) then 
c                                 not full rank
            do i = 1, hcp
               mu(i) = nopt(7)
            end do
          
         else 

         end if 

      end if 


1000  format (/,'**error ver901** solutions not allowed in saturated ',
     *   'component composition space',/,'in adaptive optimization ',
     *   'calculations, this limitation will be removed upon request.')

      end
