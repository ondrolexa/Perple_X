        subroutine chopit (ycum,jst,jsp,ksite,ids)
c---------------------------------------------------------------------
c subroutine to do cartesian or transform subdivision of species
c jst through jsp on site k of solution ids. ycum is the smallest
c fraction possible (i.e., if the minimum bound for some species 
c is > 0). the npair coordinate sets are loaded into xy(mdim,k1).
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer mres
 
      parameter (mres=12000)

      integer mode, i, ind(ms1), iy(ms1), jsp, ksite, indx, iexit, 
     *        ieyit, j, ids, jst

      double precision y(ms1,mres), ycum, ymax, dy, ync, res, ylmn,
     *                 ylmx, yloc, x, unstch, strtch, yreal

      logical odd

      integer ntot,npairs
      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24),ntot,npairs
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
      do i = 1, jsp
c                                 generate coordinates for i'th component
         iy(i) = 1
         y(i,1) = xmn(ksite,i)
         ync = xnc(ksite,i)

         if (ync.eq.0d0) cycle

         mode = imdg(i,ksite,ids)
c                                 avoid impossible compositions 'cause a min > 0
         if (i.gt.1) then 

            ycum = ycum + xmn(ksite,i-1)
c                                 1-ycum is the smallest fraction possible
            if (1d0-ycum.lt.0d0) then 
c                                 inconsistent limits
               call error (999,ycum,jsp,'cartes')

            else
c                                 the smallest fraction possible is lt
c                                 than xmax
               ymax = xmx(ksite,i)

            end if 
         else 
            ymax = xmx(ksite,i)
         end if 
c                                 two means of extracting y-range, cartesian
c                                 imod = 0 and transformation imod = 1
         if (mode.eq.0) then 
c                                 cartesian
            do while (y(i,iy(i)).lt.ymax)
               iy(i) = iy(i) + 1
               if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))
               y(i,iy(i)) = y(i,iy(i)-1) + ync
               if (dabs(y(i,iy(i))-ymax).lt.nopt(5)) then
                  y(i,iy(i)) = ymax
               else if (y(i,iy(i)).gt.ymax) then
                  y(i,iy(i)) = ymax
               end if 
            end do

         else 
c                                 conformal x is the cartesian coordinate
c                                 y is the real coordinate.
            if (mode.lt.4) then 
               odd = .false.
            else
               odd = .true.
            end if 

            res = 0d0
c                                 there are as many as intv(mode)
c                                 intervals to cycle through
            do j = 1, intv(mode)
c                                 odd or even interval?
               odd = .not.odd
c                                 interval limits              
               ylmn = yint(j,i,ksite,ids)
               ylmx = yint(j+1,i,ksite,ids)
c                                 which interval are we starting from?
               if (y(i,iy(i)).gt.ylmx-nopt(5)) cycle
c
               dy = ylmx - ylmn
c                                 pathological case y = ylmn
               if (dabs(y(i,iy(i))-ylmn).lt.nopt(5))  
     *                                y(i,iy(i)) = ylmn + nopt(5)

               if (res.eq.0d0) then 
c                                 the current value is in interval j
c                                 convert to raw y (varies from 0 ->1 
c                                 over the local interval)
                  yloc = (y(i,iy(i))-ylmn) / dy
c                                 convert to conformal x
                  if (odd) then 
                     x = unstch(yloc)
                  else 
                     x = 1d0 - unstch(1d0-yloc)
                  end if 

               else
c                                 have jumped from an earlier interval
                  x = res - ync / yfrc(j-1,i,ksite,ids)
c                 if (x.lt.0d0) x = 0d0

               end if                 
c                                 now generate all compositions in
c                                 local interval
               do while (x.le.1d0) 
c                                 increment conformal x
                  x = x + ync / yfrc(j,i,ksite,ids)
c                                 compute yreal
                  if (x.le.1d0) then 
                     if (odd) then 
                        yreal = ylmn + strtch(x) * dy
                     else
                        yreal = ylmx - strtch(1d0-x) * dy
                     end if 
 
                     iy(i) = iy(i) + 1
                     if (iy(i).gt.mres) 
     *                  call error (50,ync,mres,fname(ids))
c                                 check if in bounds
                     if (dabs(yreal-ymax).lt.nopt(5).or.
     *                                     yreal.gt.ymax) then
                        res = 0d0
                        y(i,iy(i)) = ymax
                        exit
                     else 
                        y(i,iy(i)) = yreal
                     end if 

                  else if (x.gt.1d0.and.j.eq.intv(mode)) then
c                                 at the last interval
                     iy(i) = iy(i) + 1
                     y(i,iy(i)) = ymax                
                     exit

                  else 
c                                the point is in the next interval
                     res = x - 1d0 
                     exit 

                  end if 
c                                 coordinate generating end do 
               end do 
c                                 if y is at ymax exit
               if (y(i,iy(i)).ge.ymax) exit
c                                 interval loop end do 
            end do 

         end if 
c                                 add last point if necessary, can it be? 
c                                 certainly not for conformal. 
         if (y(i,iy(i)).lt.ymax) then
            iy(i) = iy(i) + 1
            if (iy(i).gt.mres) call error (50,ync,mres,fname(ids))
            y(i,iy(i)) = ymax
         end if          
 
      end do
c                                 the first coordinate
      npairs = 1

      do i = 1, jsp
         ind(i) = 1
         simp(i) = y(i,1)
      end do
c                                 now make the array index run over all
c                                 values increasing the last index fastest
      iexit = 0 
      ieyit = 0 
      dy = 0d0

      do while (iexit.eq.0)
c                                 figure out which index to increment
         do i = jsp, 1, -1

            if (ind(i).lt.iy(i).and.ieyit.eq.0) then
c                                 this is the one to increment
               ind(i) = ind(i) + 1
               indx = i 
               exit 

            else if (i.gt.1) then 
c                                 saturated the index
               ind(i) = 1
               ieyit = 0 
               
            else
c                                 saturated first index, done.
                return 

            end if
 
         end do 
c                                 ok now we have the indices, check
c                                 the composition
         ycum = 0d0

         do i = 1, jsp 
            ycum = ycum + y(i,ind(i))  
         end do 

         if (ycum.gt.1d0) then

            ieyit = 1
c                                 here is where it gets messy:
            if (indx.eq.1) then
c                                 we're at the first point, and already
c                                 over the top
               iexit = 1
               cycle

            else if ( y(indx,ind(indx)) - y(indx,ind(indx)-1)
     *               - ycum + 1d0    .gt. nopt(5) ) then          
c                                 reset the current variable
c                                 and max loop index
               dy =  1d0 - ycum 

            else
c                                 must have just hit on the last increment
               cycle 

            end if 

         else if (ind(indx).eq.iy(indx)) then 

            ieyit = 1

         end if 

         npairs = npairs + 1
         j = (npairs-1)*jsp

         if (j+jsp.gt.k13) call error (180,ycum,k13,
     *                      'CARTES increase parameter k13')

         do i = 1, jsp
            simp(j+i) = y(i,ind(i))
         end do

         simp(j+indx) = simp(j+indx) + dy

         dy = 0d0 

      end do

      end 

      subroutine carteq (ids)
c---------------------------------------------------------------------
c subroutine to cartesian or transform subdivision on a single site
c solution with charge balance. called by subdiv. 
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, jsp, j, k, l, ids, qpairs

      double precision ycum, sum, q, ratio

      integer ntot,npairs
      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24),ntot,npairs

      integer iend,isub,imd,insp,ist,isp,isite,iterm,iord,istot,jstot,
     *        kstot,rkord,xtyp
      double precision wg,wk,reach
      common/ cst108 /wg(m1,m3),wk(m16,m17,m18),
     *      reach,iend(m4),isub(m1,m2,2),
     *      imd(msp,mst),insp(m4),ist(mst),isp(mst),rkord(m18),isite,
     *      iterm,iord,istot,jstot,kstot,xtyp
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)

      integer nq,nn,ns
      common/ cst337 /nq,nn,ns

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision thermo, uf, us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)
c----------------------------------------------------------------------

      if (nq.ne.0) then 
c                                 do the first nq-1 species independently
         ycum = 0d0
         jsp = nq - 1
 
         call chopit (ycum,1,jsp,1,ids)
c                                 at this point xy(mdim,k1) contains all 
c                                 possible compositions of the nq-1 species,
c                                 use charge balance to get the nqth species
         qpairs = 1 

         do i = 1, npairs

            q = 0d0
            sum = 0d0
            k = (i-1)*jsp
            l = (qpairs - 1)*(nq+nn+ns - 1)

            do j = 1, jsp
               q = q + thermo(6,knsp(j,ids))*simp(k+j)
               sum = sum + simp(k+j)
               prism(l+j) = simp(k+j)
            end do
c                                 charge ratio
            ratio = q/thermo(6,knsp(nq,ids))
c                                 the net charge has the same sign as the nqth
c                                 species or its amount violates closure, reject:
            if (ratio.gt.0d0.or.sum-ratio.ge.1e0) cycle
c                                 the amount of the species determined by charge balance
            prism(l+nq) = -ratio
c                                 the total moles of charged species
            prism(l+nq+1) = sum - ratio

            qpairs = qpairs + 1

         end do

         qpairs = qpairs - 1

      end if 

      end
