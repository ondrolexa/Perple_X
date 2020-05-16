
      subroutine minime (ids)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical bad

      integer ids, i, j, k, l, nvar, iter, iwork(m22),
     *        istuff(1),istate(m21), ifail, nclin, ntot

      double precision gsol1, gval, sum, objgrd(m19),
     *                 bl(m21), bu(m21), gfinal, dbz, ppp(m19), 
     *                 clamda(m21),r(m19,m19),work(m23),stuff(1),
c                                 dummies for NCNLN > 0
     *                 c(1),cjac(1,1)

      integer jds, tphct
      double precision az, bz
      common/ cstexp /az(m20,m19),bz(m20),jds,tphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)

      integer lterm, ksub
      common/ cxt1i /lterm(m11,m10,h9),ksub(m0,m11,m10,h9)

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      external gsol2, dummy

      external gsol1
c-----------------------------------------------------------------------
c                                 create constraint matrix z
c                                 --------------------------
      ntot = nstot(ids)
      nvar = ntot - 1
c                                 initialize bounds
      do i = 1, nvar
         bu(i) = 1d20
         bl(i) = -1d20
      end do 

      if (.not.mus) then 
         write (*,*) 'no mus'
         call errpau
      end if 

      jds = ids
      tphct = jphct
c----------------------------------------------------------
      write (*,1000) 0, (pp(i),i=1,nstot(jds))
      write (*,1000) 0, (p0a(i),i=1,nstot(jds))
      write (*,1000) 0, (pa(i),i=1,nstot(jds))
      gval = gsol1 (jds)
      write (*,1000) gval, (pp(i),i=1,nstot(jds))
      write (*,1000) gval, (p0a(i),i=1,nstot(jds))
      write (*,1000) gval, (pa(i),i=1,nstot(jds))

      sum = 0d0 

      do i = 1, nstot(jds)
         sum = sum + pa(i)
         p0a(i) = pa(i)
         pp(i) = pa(i)
      end do

      call makepp (jds)

      write (*,*) sum, g2(tphct)
      
      g2(tphct) = gsol1 (jds)

      call csol (tphct,jds,bad)

      gval = g2(tphct)

      do i = 1, icp
         gval = gval - cp2(i,tphct)*mu(i)
      end do

      write (*,1000) g2(tphct), (cp2(i,tphct),i=1,icp)

      write (*,1000) gval, (pp(i),i=1,nstot(jds))
      write (*,1000) gval, (p0a(i),i=1,nstot(jds))
      write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g14.6,1x,12(f8.5,1x))
c----------------------------------------------------------
      call makp2y (ids)

      call maky2x (ids)
c                                 inequalities, to be counted:
      nclin = 0
c                                 for each site
      do i = 1, msite(ids)
c                                 for each species
         do j = 1, zsp(ids,i) + 1
c                                 initial az, bz
            nclin = nclin + 1
            dbz = 0d0
c                                 both Temkin and non-Temkin have
c                                 Az*p >= 0 constraints:
            bl(nvar+nclin) = -dcoef(0,j,i,ids)

            do k = 1, nvar

               az(nclin,k) = 0d0

            end do

            do k = 1, lterm(j,i,ids)

               if (ksub(k,j,i,ids).ne.ntot) then 

                  az(nclin,ksub(k,j,i,ids)) = az(nclin,k) 
     *                                     + dcoef(k,j,i,ids)

               else 
c                                 decompose p(ntot) into 1 - p(1) - ...- p(nvar)
                  dbz = dcoef(k,j,i,ids)

                  do l = 1, nvar

                     az(nclin,l) = az(nclin,l) - dbz

                  end do

               end if 

            end do

            bl(nvar+nclin) = bz(nvar+nclin) - dbz
c                                 non-Temkin have the Az*p <= 1 constraint
            if (zmult(ids,i).ne.0d0) then

               bu(nvar+nclin) =  1d0 - dcoef(0,j,i,ids) + dbz

            else 

               bu(nvar+nclin) = 1d20

            end if

         end do 

      end do

      ifail = -1
c auto
c      CALL E04UEF ('difference interval = 1d-3')
c eps = 1e-5
c sqrt(eps)
c      CALL E04UEF ('linear feasibility tolerance = 0.0032')
c sqrt(eps)
c      CALL E04UEF ('feasibility tolerance = 0.0032')
c eps**(0.8)
c      CALL E04UEF ('optimality tolerance = 1d-4')
c eps**(0.9)
c      CALL E04UEF ('function precision = 0.0000316')
c eps = 1e-4
c sqrt(eps)
c      CALL E04UEF ('linear feasibility tolerance = 0.01')
c sqrt(eps)
c      CALL E04UEF ('feasibility tolerance = 0.01')
c eps**(0.8)
c      CALL E04UEF ('optimality tolerance = 0.00063')
c eps**(0.9)
c      CALL E04UEF ('function precision = 0.00025')
       CALL E04UEF ('derivative level = 0')

      call e04ucf (nvar,nclin,0,m20,1,m19,az,bl,bu,dummy,gsol2,
     *             iter,istate,c,cjac,clamda,gfinal,objgrd,r,ppp,iwork,
     *             m22,work,m23,istuff,stuff,ifail)

      end

      subroutine gsol2 (mode,nvar,ppp,gval,objgrd,istart,istuff,stuff)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, nvar, istart, mode, istuff(*)

      logical bad

      double precision ppp(*), gval, gsol1, psum, objgrd(*), stuff(*)

      external gsol1

      integer jds, tphct
      double precision az, bz
      common/ cstexp /az(m20,m19),bz(m20),jds,tphct

      logical mus
      double precision mu
      common/ cst330 /mu(k8),mus

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer jphct
      double precision g2, cp2, c2tot
      common/ cxt12 /g2(k21),cp2(k5,k21),c2tot(k21),jphct

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c-----------------------------------------------------------------------
      psum = 0d0

      do i = 1, nvar
         psum = psum + ppp(i)
         p0a(i) = ppp(i)
      end do

      p0a(nstot(jds)) = 1d0 - psum
      
      do i = 1, nstot(jds)
         pp(i)  = p0a(i)
         pa(i)  = p0a(i)
      end do

      call makepp (jds)

      do i = 1, 36
         psum = 0d0
         do j = 1, nstot(jds)-1
            psum = psum + pp(j)*az(i,j)
         end do 
c        write (*,1010) psum, bz(i)
c        if (psum.gt.bz(i)) write (*,*) 'oink'
      end do 

      write (*,1000) 0, (pa(i),i=1,nstot(jds))

      g2(tphct) = gsol1(jds)
      gval = g2(tphct)

      call csol (tphct,jds,bad)

      write (*,1000) 0, (cp2(i,tphct),i=1,icp)

      do i = 1, icp
         gval = gval - cp2(i,tphct)*mu(i)
      end do

      write (*,1000) gval, (pa(i),i=1,nstot(jds))

1000  format (g12.6,1x,12(f8.5,1x))
1010  format (2(g14.7,2x))

      end

      subroutine dummy
      end 

      subroutine makp2y (id)
c----------------------------------------------------------------------
c subroutine to construct the Ay matrix for the p to y conversion
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, k, l, m, id

      double precision ay(h9,m14,m4), xt

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
c                                 the Ay matrix in Ay*y = p
c                                 y is dimension mstot < m4
c                                 p is dimension nstot < m14
c                                 Ay(nstot,mstot)
c                                 assemble Ay
      do i = 1, nstot(id)
         do j = 1, mstot(id)
            ay(id,i,j) = 0d0 
         end do 
      end do
c                                 set independent disordered
c                                 endmember coefficients
      do i = 1, lstot(id)
c                                 independent endmember i has
c                                 index k in the y array
         ay(id,i,knsp(i,id)) = 1d0

      end do
c                                 fill in contributions from the
c                                 dependent endmembers:
      do k = 1, ndep(id)

         j = knsp(lstot(id)+k,id)
c                                 y(j) decomposes
c                                 to p(i) as y2pg(k,k)*y()
         do i = 1, nstot(id)

             ay(id,i,j) = ay(id,i,j) + y2pg(k,i,id)

         end do

      end do
c                                  test
      do i = 1, nstot(id)
         xt = 0d0
         do j = 1, mstot(id)
            xt = xt + ay(id,i,j)*y(j)
         end do 

         write (*,1000) i, xt, p0a(i), pa(i)

      end do

1000  format (i,3(3x,g14.6))

      end


      subroutine maky2x (id)
c----------------------------------------------------------------------
c subroutine to construct the Ax matrices for the y to x conversion
c of each subcomposition.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ii, i, j, k, l, m, id, xcomps

      double precision ax(h9,h4,mst*msp,m4), xt

      double precision z, pa, p0a, x, w, y, wl, pp
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(h4,mst,msp),w(m1),
     *              wl(m17,m18),pp(m4)
c----------------------------------------------------------------------
      do ii = 1, pop1(id)
c                                 for subcomposition ii 
c                                 the Ax matrix in Ax*y = x
c                                 m = sum( ispg(1:istg) )
c                                 n = pvert(2) - pvert(1) + 1

         xcomps = 0

         do i = 1, istg(id,ii)
            xcomps = xcomps + ispg(id,ii,i)
         end do 

         do i = 1, xcomps
            do j = 1, pvert(id,ii,2) - pvert(id,ii,2) + 1
               ax(id,ii,i,j) = 0d0 
            end do 
         end do

         do k = pvert(id,ii,1), pvert(id,ii,2)
c                                 for each endmember in the subcomposition
c                                 load the column of ax:
            j = k - pvert(id,ii,1) + 1

            i = 0

            do l = 1, istg(id,ii)

               do m = 1, ispg(id,ii,l)
c                                 kmsol indicates the species on 
c                                 the site m of endmember k
                  if (kmsol(id,k,l).eq.m) then

                     ax(id,ii,i+m,j) = 1d0

                     exit 

                  end if

               end do

               i = i + ispg(id,ii,l)

            end do 

         end do 

      end do
c                                  test
      do ii = 1, poly(id)
c                                  get the polytope weight
         if (pop1(id).eq.1) then 

            pwt(ii) = 1d0

         else

            pwt(ii) = 0d0 

            do k = pvert(id,ii,1), pvert(id,ii,2)
               pwt(ii) = pwt(ii) + y(k)
            end do

         end if

         xcomps = 0

         do i = 1, istg(id,ii)
            xcomps = xcomps + ispg(ii,id,i)
         end do

         l = 1
         m = 1

         do i = 1, xcomps

            xt = 0d0
c                                 the algebra
            do k = pvert(id,ii,1), pvert(id,ii,2)

               j = k - pvert(id,ii,1) + 1

               xt = xt + ax(id,ii,i,j)*y(k)

            end do

            write (*,1000) ii,l,m,xt,x(ii,l,m)

            m = m + 1

            if (m.gt.ispg(id,ii,l)) then 
               l = l + 1
               m = 1
            end if

         end do

      end do

1000  format (3(i2,1x),3(3x,g14.6))

      end
